### Baseline hazards:
# 1:5 = c("exponential", "weibull", "gompertz", "weibull-weibull (1)", "weibull-weibull (2)")

### Frailty distributions:
# 1:2 = c("Gamma", "Normal")

### Simulate data
make_data <- function(n_individuals, n_clusters, fv, fv_dist = 1:2, treatment_effect, distribution = 1:5, pars, maxt = 5, scenario) {
  # Get RNG seed
  seed <- .Random.seed

  # Make covariates data.frame
  X <- data.frame(
    id = 1:(n_individuals * n_clusters),
    grpid = rep(1:n_clusters, each = n_individuals),
    trt = stats::rbinom(n_individuals * n_clusters, 1, 0.5)
  )

  # Make frailty
  if (fv_dist == 1) {
    X[["frvec"]] <- rep(log(stats::rgamma(n_clusters, shape = 1 / fv, scale = fv)), each = n_individuals)
  } else {
    X[["frvec"]] <- rep(stats::rnorm(n_clusters, mean = 0, sd = sqrt(fv)), each = n_individuals)
  }

  # Make survival indicator and status variable
  if (distribution == 1) {
    S <- simsurv::simsurv(dist = "exponential", lambdas = pars[["lambda"]], x = X, betas = c(trt = treatment_effect, frvec = 1), maxt = maxt)
  } else if (distribution == 2) {
    S <- simsurv::simsurv(dist = "weibull", lambdas = pars[["lambda"]], gammas = pars[["gamma"]], x = X, betas = c(trt = treatment_effect, frvec = 1), maxt = maxt)
  } else if (distribution == 3) {
    S <- simsurv::simsurv(dist = "gompertz", lambdas = pars[["lambda"]], gammas = pars[["gamma"]], x = X, betas = c(trt = treatment_effect, frvec = 1), maxt = maxt)
  } else {
    S <- simsurv::simsurv(dist = "weibull", mixture = TRUE, lambdas = pars[["lambda"]], gammas = pars[["gamma"]], pmix = pars[["pmix"]], x = X, betas = c(trt = treatment_effect, frvec = 1), maxt = maxt, interval = c(0, 500))
  }

  # Make data frame
  out <- merge(X, S)
  out <- dplyr::rename(out, t = eventtime, d = status)
  out[["scenario"]] <- scenario
  out <- dplyr::arrange(out, grpid)

  # Winsorising tiny values of t (smaller than one day, e.g. 1 / 365.242)
  # We also add a tiny amount of white noise, and we make sure the resulting values are positive
  out[["t"]][out[["t"]] < 1 / 365.242] <- 1 / 365.242 + rnorm(length(out[["t"]][out[["t"]] < 1 / 365.242]), mean = 0, sd = 1e-4)
  out[["t"]] <- abs(out[["t"]])

  # Attach .Random.seed as attribute
  attr(out, "seed") <- seed

  # Return simulated data frame
  return(out)
}

### Packages
library("dplyr")
library("tidyr")
library("simsurv")

### Define DGMs
# Simulation factors (not fully factorial)
dgms <- crossing(
  data.frame(
    n_individuals = c(2, 10, 50, 250),
    n_clusters = c(750, 100, 50, 15)
  ),
  treatment_effect = -0.50,
  fv = c(0.25, 0.75, 1.25),
  fv_dist = 1:2,
  baseline = 1:5
) %>%
  full_join(
    # Parameters for the different baseline hazard functions
    data_frame(
      baseline = 1:5,
      pars = list(
        list(lambda = 0.5),
        list(lambda = 0.5, gamma = 0.8),
        list(lambda = 0.5, gamma = 0.2),
        list(lambda = c(0.3, 0.5), gamma = c(1.5, 2.5), pmix = 0.7),
        list(lambda = c(0.5, 0.5), gamma = c(1.3, 0.7), pmix = 0.5)
      )
    ),
    "baseline"
  )

### Generate B datasets for a given scenario:
d <- as.numeric(Sys.getenv("PBS_ARRAYID"))
B <- 1000
out <- vector(mode = "list", length = B)
for (b in 1:B) {
  out[[b]] <- make_data(
    n_individuals = dgms[["n_individuals"]][d],
    n_clusters = dgms[["n_clusters"]][d],
    fv = dgms[["fv"]][d],
    fv_dist = dgms[["fv_dist"]][d],
    treatment_effect = dgms[["treatment_effect"]][d],
    distribution = dgms[["baseline"]][d],
    pars = dgms[["pars"]][[d]],
    maxt = 5,
    scenario = d
  )
}

# Save dataset
saveRDS(object = out, file = paste0("Data/simdata_", sprintf("%03d", d), ".RDS"))
