### Required packages
library(tidyverse)
library(pracma)
library(rsimsum)

### All jobs
looktab <- tidyr::crossing(
  scenario = seq(90),
  model = seq(22)
)

### Load results
# make a single file for better storage, faster loading, etc.
# done <- lapply(X = list.files(path = "Results", pattern = "^jout_|^multiout_", full.names = TRUE), FUN = readRDS)
# saveRDS(object = done, file = "Results.RDS")
done <- readRDS(file = "Results.RDS")

# Data cleaning
# Some results were duplicated because of quirks with the HPC server,
# but point estimates and standard errors are the same.
# Taking the average of them.
results <- dplyr::bind_rows(done) %>%
  dplyr::group_by(scenario, model, i, par) %>%
  dplyr::summarise(
    coef = mean(coef, na.rm = TRUE),
    se = mean(se, na.rm = TRUE),
    time = mean(time, na.rm = TRUE)
  ) %>%
  dplyr::ungroup()

# Some estimates SEs are Inf, set to maximum integer that this machine can represent
results <- results %>%
  mutate(se = ifelse(se == Inf, .Machine$integer.max, se))

# Set to NA values too large (|value| > 10 times the (robust)-standardised value), both coef and se
results <- results %>%
  dplyr::group_by(scenario, model, par) %>%
  dplyr::arrange(scenario, model, par) %>%
  dplyr::mutate(
    coef.me = median(coef, na.rm = TRUE),
    se.me = median(se, na.rm = TRUE),
    coef.iqr = fivenum(coef, na.rm = TRUE)[4] - fivenum(coef, na.rm = TRUE)[2],
    se.iqr = fivenum(se, na.rm = TRUE)[4] - fivenum(se, na.rm = TRUE)[2],
    coef.std = (coef - coef.me) / coef.iqr,
    se.std = (se - se.me) / se.iqr
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    coef = ifelse(abs(coef.std) > 10, NA, coef),
    se = ifelse(abs(se.std) > 10, NA, se)
  )

# Set to NA when only one of coef or se is NA
results <- results %>%
  dplyr::mutate(isna = is.na(coef) | is.nan(coef) | is.na(se) | is.nan(se)) %>%
  dplyr::group_by(scenario, model, i) %>%
  dplyr::mutate(anyna = max(isna)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    coef = ifelse(anyna == 1, NA, coef),
    se = ifelse(anyna == 1, NA, se)
  )

# Remove variables not needed anymore
results <- results %>%
  dplyr::select(-coef.me, -se.me, -coef.iqr, -se.iqr, -coef.std, -se.std, -isna, -anyna)

### Calculate summary statistics
# Split results
lofresults <- split(results, f = lapply(c("scenario", "par"), function(f) results[[f]]))

# All DGMs
dgms <- tidyr::crossing(
  tibble::tibble(
    n_individuals = c(2, 150),
    n_clusters = c(750, 20)
  ),
  treatment_effect = -0.50,
  fv = c(0.25, 0.75, 1.25),
  fv_dist = 1:3, # c("Gamma", "Log-Normal", "Mixture Normal")
  baseline = 1:5 # c("Exponential", "Weibull", "Gompertz", "Weibull-Weibull (1)", "Weibull-Weibull (2)")
) %>%
  dplyr::full_join(
    # Parameters for the different baseline hazard functions
    tibble::tibble(
      baseline = 1:5, # c("Exponential", "Weibull", "Gompertz", "Weibull-Weibull (1)", "Weibull-Weibull (2)")
      pars = list(
        list(lambda = 0.5),
        list(lambda = 0.5, gamma = 0.8),
        list(lambda = 0.5, gamma = 0.2),
        list(lambda = c(0.3, 0.5), gamma = c(1.5, 2.5), pmix = 0.7),
        list(lambda = c(0.5, 0.5), gamma = c(1.3, 0.7), pmix = 0.5)
      )
    ),
    "baseline"
  ) %>%
  dplyr::mutate(scenario = dplyr::row_number())

# Make true LLE for each scenario
true_ms <- Vectorize(function(t, X, distribution, fv_dist, fv, beta = -0.50) {
  if (distribution == 1) {
    lambda <- dgms$pars[dgms$baseline == distribution][[1]]$lambda
    if (fv_dist == 1) {
      (1 - fv * (-lambda * exp(beta * X) * t))^(-1 / fv)
    } else {
      S <- exp(-lambda * exp(beta * X) * t)
      if (fv_dist == 2) {
        fn <- function(eta) S^exp(eta) * dnorm(eta, mean = 0, sd = sqrt(fv))
        pracma::quadinf(f = fn, xa = -Inf, xb = Inf)$Q
      } else if (fv_dist == 3) {
        fn <- function(eta) S^exp(eta) * LaplacesDemon::dnormm(x = eta, p = c(0.5, 0.5), mu = c(-3 * sqrt(fv), +3 * sqrt(fv)), sigma = c(sqrt(fv), sqrt(fv)))
        pracma::quadinf(f = fn, xa = -Inf, xb = Inf)$Q
      }
    }
  } else if (distribution == 2) {
    lambda <- dgms$pars[dgms$baseline == distribution][[1]]$lambda
    gamma <- dgms$pars[dgms$baseline == distribution][[1]]$gamma
    if (fv_dist == 1) {
      (1 - fv * (-lambda * exp(beta * X) * t^gamma))^(-1 / fv)
    } else {
      S <- exp(-lambda * exp(beta * X) * t^gamma)
      if (fv_dist == 2) {
        fn <- function(eta) S^exp(eta) * dnorm(eta, mean = 0, sd = sqrt(fv))
        pracma::quadinf(f = fn, xa = -Inf, xb = Inf)$Q
      } else if (fv_dist == 3) {
        fn <- function(eta) S^exp(eta) * LaplacesDemon::dnormm(x = eta, p = c(0.5, 0.5), mu = c(-3 * sqrt(fv), +3 * sqrt(fv)), sigma = c(sqrt(fv), sqrt(fv)))
        pracma::quadinf(f = fn, xa = -Inf, xb = Inf)$Q
      }
    }
  } else if (distribution == 3) {
    lambda <- dgms$pars[dgms$baseline == distribution][[1]]$lambda
    gamma <- dgms$pars[dgms$baseline == distribution][[1]]$gamma
    if (fv_dist == 1) {
      (1 - fv * ((-lambda * exp(beta * X) / gamma) * (exp(gamma * t) - 1)))^(-1 / fv)
    } else {
      S <- exp((-lambda * exp(beta * X) / gamma) * (exp(gamma * t) - 1))
      if (fv_dist == 2) {
        fn <- function(eta) S^exp(eta) * dnorm(eta, mean = 0, sd = sqrt(fv))
        pracma::quadinf(f = fn, xa = -Inf, xb = Inf)$Q
      } else if (fv_dist == 3) {
        fn <- function(eta) S^exp(eta) * LaplacesDemon::dnormm(x = eta, p = c(0.5, 0.5), mu = c(-3 * sqrt(fv), +3 * sqrt(fv)), sigma = c(sqrt(fv), sqrt(fv)))
        pracma::quadinf(f = fn, xa = -Inf, xb = Inf)$Q
      }
    }
  } else if (distribution == 4) {
    lambda1 <- dgms$pars[dgms$baseline == distribution][[1]]$lambda[1]
    lambda2 <- dgms$pars[dgms$baseline == distribution][[1]]$lambda[2]
    gamma1 <- dgms$pars[dgms$baseline == distribution][[1]]$gamma[1]
    gamma2 <- dgms$pars[dgms$baseline == distribution][[1]]$gamma[2]
    pmix <- dgms$pars[dgms$baseline == distribution][[1]]$pmix
    if (fv_dist == 1) {
      (1 - fv * log((pmix * exp(-lambda1 * t^gamma1) + (1 - pmix) * exp(-lambda2 * t^gamma2))^exp(beta * X)))^(-1 / fv)
    } else {
      S <- ((pmix * exp(-lambda1 * t^gamma1) + (1 - pmix) * exp(-lambda2 * t^gamma2))^exp(beta * X))
      if (fv_dist == 2) {
        fn <- function(eta) S^exp(eta) * dnorm(eta, mean = 0, sd = sqrt(fv))
        pracma::quadinf(f = fn, xa = -Inf, xb = Inf)$Q
      } else if (fv_dist == 3) {
        fn <- function(eta) S^exp(eta) * LaplacesDemon::dnormm(x = eta, p = c(0.5, 0.5), mu = c(-3 * sqrt(fv), +3 * sqrt(fv)), sigma = c(sqrt(fv), sqrt(fv)))
        pracma::quadinf(f = fn, xa = -Inf, xb = Inf)$Q
      }
    }
  } else if (distribution == 5) {
    lambda1 <- dgms$pars[dgms$baseline == distribution][[1]]$lambda[1]
    lambda2 <- dgms$pars[dgms$baseline == distribution][[1]]$lambda[2]
    gamma1 <- dgms$pars[dgms$baseline == distribution][[1]]$gamma[1]
    gamma2 <- dgms$pars[dgms$baseline == distribution][[1]]$gamma[2]
    pmix <- dgms$pars[dgms$baseline == distribution][[1]]$pmix
    if (fv_dist == 1) {
      (1 - fv * log((pmix * exp(-lambda1 * t^gamma1) + (1 - pmix) * exp(-lambda2 * t^gamma2))^exp(beta * X)))^(-1 / fv)
    } else {
      S <- ((pmix * exp(-lambda1 * t^gamma1) + (1 - pmix) * exp(-lambda2 * t^gamma2))^exp(beta * X))
      if (fv_dist == 2) {
        fn <- function(eta) S^exp(eta) * dnorm(eta, mean = 0, sd = sqrt(fv))
        pracma::quadinf(f = fn, xa = -Inf, xb = Inf)$Q
      } else if (fv_dist == 3) {
        fn <- function(eta) S^exp(eta) * LaplacesDemon::dnormm(x = eta, p = c(0.5, 0.5), mu = c(-3 * sqrt(fv), +3 * sqrt(fv)), sigma = c(sqrt(fv), sqrt(fv)))
        pracma::quadinf(f = fn, xa = -Inf, xb = Inf)$Q
      }
    }
  }
})

# Make LLE
lle <- dplyr::distinct(dgms, fv, fv_dist, baseline)
lle$lle <- sapply(1:nrow(lle), function(i) {
  cat(i, "\n")
  t <- seq(0, 5, length.out = 1000)
  ms1 <- with(lle, true_ms(t = t, X = 1, distribution = baseline[i], fv_dist = fv_dist[i], fv = fv[i]))
  ms0 <- with(lle, true_ms(t = t, X = 0, distribution = baseline[i], fv_dist = fv_dist[i], fv = fv[i]))
  ms1f <- stats::splinefun(x = t, y = ms1)
  ms0f <- stats::splinefun(x = t, y = ms0)
  is1 <- pracma::quadinf(f = ms1f, xa = 0, xb = 5)$Q
  is0 <- pracma::quadinf(f = ms0f, xa = 0, xb = 5)$Q
  is1 - is0
})

# Merge back
dgms <- dgms %>%
  dplyr::left_join(lle, c("fv", "fv_dist", "baseline"))

# Make data.frame with true values
true <- tidyr::crossing(
  scenario = seq(90),
  par = c("fv", "lle", "trt")
) %>%
  dplyr::left_join(dplyr::distinct(dgms, scenario, fv) %>%
    dplyr::mutate(par = "fv"),
  by = c("scenario", "par")
  ) %>%
  dplyr::left_join(dplyr::distinct(dgms, scenario, lle) %>%
    dplyr::mutate(par = "lle"),
  by = c("scenario", "par")
  ) %>%
  dplyr::left_join(dplyr::distinct(dgms, scenario, treatment_effect) %>%
    dplyr::mutate(par = "trt"),
  by = c("scenario", "par")
  ) %>%
  dplyr::mutate(true = dplyr::case_when(
    par == "fv" ~ fv,
    par == "lle" ~ lle,
    par == "trt" ~ treatment_effect
  )) %>%
  dplyr::select(scenario, par, true)

# Call simsum on each object of the list to summarise the results
lofsummary <- lapply(seq_along(lofresults), function(i) {
  tv <- dplyr::left_join(lofresults[[i]],
    true,
    by = c("scenario", "par")
  ) %>%
    dplyr::distinct(true) %>%
    dplyr::pull()
  s <- rsimsum::simsum(lofresults[[i]], true = tv, estvarname = "coef", se = "se", methodvar = "model", by = "scenario", ref = "1")
  smr <- summary(s)
  out <- rsimsum::get_data(smr) %>%
    dplyr::filter(stat %in% c("nsim", "thetamean", "thetamedian", "se2mean", "se2median", "bias", "empse", "mse", "cover"))
  out[["par"]] <- unique(lofresults[[i]][["par"]])
  out
})

# Bind each scenario in a data.frame
summary <- dplyr::bind_rows(lofsummary) %>%
  dplyr::arrange(par, stat, scenario, model)

# Factorise
summary <- summary %>%
  dplyr::left_join(
    dplyr::mutate(dgms, scenario = as.character(scenario)),
    by = "scenario"
  ) %>%
  dplyr::select(-pars) %>%
  dplyr::mutate(
    model = factor(model, levels = 1:22, labels = c("Cox, Gamma", "Cox, log-Normal", "Exp, Gamma", "Exp, log-Normal", "Weibull, Gamma", "Weibull, log-Normal", "Gompertz, Gamma", "Gompertz, log-Normal", "RP (3), Gamma", "RP (3), log-Normal", "RP (5), Gamma", "RP (5), log-Normal", "RP (9), Gamma", "RP (9), log-Normal", "RP (P), Gamma", "RP (P), log-Normal", "FP (W), Gamma", "FP (W), log-Normal", "FP (k=10), Gamma", "FP (k=10), log-Normal", "FP (k=10000), Gamma", "FP (k=10000), log-Normal")),
    model2 = factor(model, levels = levels(model)[c(1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22)]),
    fv = factor(fv, levels = c(0.25, 0.75, 1.25)),
    fv_dist = factor(fv_dist, levels = 1:3, labels = c("Gamma", "log-Normal", "Mixture Normal")),
    baseline = factor(baseline, levels = 1:5, labels = c("Exponential", "Weibull", "Gompertz", "Weibull-Weibull (1)", "Weibull-Weibull (2)"))
  ) %>%
  tidyr::separate(model, into = c("mbaseline", "mfrailty"), sep = ", ", remove = FALSE) %>%
  dplyr::mutate(
    mbaseline = factor(mbaseline, levels = c("Cox", "Exp", "Weibull", "Gompertz", "RP (3)", "RP (5)", "RP (9)", "RP (P)", "FP (W)", "FP (k=10)", "FP (k=10000)")),
    mfrailty = factor(mfrailty, levels = c("Gamma", "log-Normal"))
  ) %>%
  dplyr::mutate(par = factor(par, levels = c("trt", "lle", "fv"), labels = c("Treatment effect", "LLE", "Frailty variance"))) %>%
  dplyr::mutate(sig = dplyr::case_when(
    stat == "bias" & (lower > 0 | upper < 0) ~ "*",
    stat == "cover" & (lower > 0.95 | upper < 0.95) ~ "*",
    TRUE ~ ""
  )) %>%
  dplyr::arrange(scenario, par, stat, model)

# Export results
saveRDS(object = summary, file = "summary.RDS")

# Export results for exploring convergence issues
convergence <- results %>%
  dplyr::filter(par == "trt") %>%
  dplyr::left_join(dgms, by = "scenario") %>%
  dplyr::mutate(ss = paste0(n_clusters, "c x ", n_individuals, "i")) %>%
  dplyr::select(-time, -treatment_effect, -pars, -lle, -par, -n_clusters, -n_individuals) %>%
  dplyr::mutate(
    model = factor(model, levels = 1:22, labels = c("Cox, Gamma", "Cox, log-Normal", "Exp, Gamma", "Exp, log-Normal", "Weibull, Gamma", "Weibull, log-Normal", "Gompertz, Gamma", "Gompertz, log-Normal", "RP (3), Gamma", "RP (3), log-Normal", "RP (5), Gamma", "RP (5), log-Normal", "RP (9), Gamma", "RP (9), log-Normal", "RP (P), Gamma", "RP (P), log-Normal", "FP (W), Gamma", "FP (W), log-Normal", "FP (k=10), Gamma", "FP (k=10), log-Normal", "FP (k=10000), Gamma", "FP (k=10000), log-Normal")),
    fv = factor(fv, levels = c(0.25, 0.75, 1.25)),
    fv_dist = factor(fv_dist, levels = 1:3, labels = c("Gamma", "log-Normal", "Mixture Normal")),
    baseline = factor(baseline, levels = 1:5, labels = c("Exponential", "Weibull", "Gompertz", "Weibull-Weibull (1)", "Weibull-Weibull (2)")),
    ss = factor(ss, labels = c("20c x 150i", "750c x 2i"))
  ) %>%
  tidyr::separate(model, into = c("mbaseline", "mfrailty"), sep = ", ", remove = FALSE) %>%
  dplyr::mutate(
    mbaseline = factor(mbaseline, levels = c("Cox", "Exp", "Weibull", "Gompertz", "RP (3)", "RP (5)", "RP (9)", "RP (P)", "FP (W)", "FP (k=10)", "FP (k=10000)")),
    mfrailty = factor(mfrailty, levels = c("Gamma", "log-Normal"))
  ) %>%
  dplyr::mutate(converged = as.numeric(!is.na(coef) & !is.na(se))) %>%
  dplyr::mutate(nonconverged = 1 - converged) %>%
  dplyr::mutate(si = as.numeric(factor(paste(scenario, i))))
saveRDS(object = convergence, file = "convergence.RDS")
