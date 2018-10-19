### Required packages
library(tidyverse)
library(pracma)
library(rsimsum)

### All jobs
looktab <- expand.grid(
  scenario = seq(120),
  model = seq(22)
)

### Load results
# make a single file for better storage, faster loading, etc.
# done <- lapply(list.files("Results", pattern = "^out_|^multiout_", full.names = TRUE), readRDS)
# saveRDS(done, "Results.RDS")
done <- readRDS("Results.RDS")

# Data cleaning
# Some results were duplicated because of quirks with the HPC server,
# but point estimates and standard errors are the same.
# Taking the average of them.
results <- bind_rows(done) %>%
  group_by(scenario, model, i, par) %>%
  summarise(
    coef = mean(coef, na.rm = TRUE),
    se = mean(se, na.rm = TRUE),
    time = mean(time, na.rm = TRUE)
  ) %>%
  ungroup()

# Set to NA values too large (|value| > 10 times the standardised value), both coef and se
results <- group_by(results, scenario, model, par) %>%
  arrange(scenario, model, par) %>%
  mutate(
    coef.std = (coef - median(coef, na.rm = TRUE)) / (fivenum(coef, na.rm = TRUE)[4] - fivenum(coef, na.rm = TRUE)[2]),
    se.std = (se - median(se, na.rm = TRUE)) / (fivenum(se, na.rm = TRUE)[4] - fivenum(se, na.rm = TRUE)[2])
  ) %>%
  ungroup() %>%
  mutate(
    coef = ifelse(abs(coef.std) > 10, NA, coef),
    se = ifelse(abs(se.std) > 10, NA, se)
  )

# Set to NA when only one of coef or se is NA
results <- results %>%
  mutate(isna = is.na(coef) | is.nan(coef) | is.na(se) | is.nan(se)) %>%
  group_by(scenario, model, i) %>%
  mutate(anyna = max(isna)) %>%
  ungroup() %>%
  mutate(
    coef = ifelse(anyna == 1, NA, coef),
    se = ifelse(anyna == 1, NA, se)
  )

# Remove variables not needed anymore
results <- results %>%
  select(-coef.std, -se.std, -isna, -anyna)

### Calculate summary statistics
# Split results
lofresults <- split(results, f = lapply(c("scenario", "par"), function(f) results[[f]]))

# All DGMs
dgms <- crossing(
  data.frame(
    n_individuals = c(2, 10, 50, 250),
    n_clusters = c(750, 100, 50, 15)
  ),
  treatment_effect = -0.50,
  fv = c(0.25, 0.75, 1.25),
  fv_dist = 1:2, # c("Gamma", "Log-Normal")
  baseline = 1:5 # c("Exponential", "Weibull", "Gompertz", "Weibull-Weibull (1)", "Weibull-Weibull (2)")
) %>%
  full_join(
    # Parameters for the different baseline hazard functions
    data_frame(
      baseline = 1:5, # c("Exponential", "Weibull", "Gompertz", "Weibull-Weibull (1)", "Weibull-Weibull (2)"),
      pars = list(
        list(lambda = 0.5),
        list(lambda = 0.5, gamma = 0.8),
        list(lambda = 0.5, gamma = 0.2),
        list(lambda1 = 0.3, lambda2 = 0.5, p1 = 1.5, p2 = 2.5, pmix = 0.7),
        list(lambda1 = 0.5, lambda2 = 0.5, p1 = 1.3, p2 = 0.7, pmix = 0.5)
      )
    ),
    "baseline"
  ) %>%
  mutate(scenario = row_number())

# Make true LLE for each scenario
true_ms <- Vectorize(function(t, X, distribution, fv_dist, fv, beta = -0.50) {
  if (distribution == 1) {
    lambda <- 0.5
    if (fv_dist == 1) {
      (1 - fv * (-lambda * exp(beta * X) * t))^(-1 / fv)
    } else {
      S <- exp(-lambda * exp(beta * X) * t)
      fn <- function(eta) S^exp(eta) * dnorm(eta, mean = 0, sd = sqrt(fv))
      quadinf(f = fn, xa = -Inf, xb = Inf)$Q
    }
  } else if (distribution == 2) {
    lambda <- 0.5
    gamma <- 0.8
    if (fv_dist == 1) {
      (1 - fv * (-lambda * exp(beta * X) * t^gamma))^(-1 / fv)
    } else {
      S <- exp(-lambda * exp(beta * X) * t^gamma)
      fn <- function(eta) S^exp(eta) * dnorm(eta, mean = 0, sd = sqrt(fv))
      quadinf(f = fn, xa = -Inf, xb = Inf)$Q
    }
  } else if (distribution == 3) {
    lambda <- 0.5
    gamma <- 0.2
    if (fv_dist == 1) {
      (1 - fv * ((-lambda * exp(beta * X) / gamma) * (exp(gamma * t) - 1)))^(-1 / fv)
    } else {
      S <- exp((-lambda * exp(beta * X) / gamma) * (exp(gamma * t) - 1))
      fn <- function(eta) S^exp(eta) * dnorm(eta, mean = 0, sd = sqrt(fv))
      quadinf(f = fn, xa = -Inf, xb = Inf)$Q
    }
  } else if (distribution == 4) {
    lambda1 <- 0.3
    lambda2 <- 0.5
    gamma1 <- 1.5
    gamma2 <- 2.5
    pmix <- 0.7
    if (fv_dist == 1) {
      (1 - fv * log((pmix * exp(-lambda1 * t^gamma1) + (1 - pmix) * exp(-lambda2 * t^gamma2))^exp(beta * X)))^(-1 / fv)
    } else {
      S <- ((pmix * exp(-lambda1 * t^gamma1) + (1 - pmix) * exp(-lambda2 * t^gamma2))^exp(beta * X))
      fn <- function(eta) S^exp(eta) * dnorm(eta, mean = 0, sd = sqrt(fv))
      quadinf(f = fn, xa = -Inf, xb = Inf)$Q
    }
  } else if (distribution == 5) {
    lambda1 <- 0.5
    lambda2 <- 0.5
    gamma1 <- 1.3
    gamma2 <- 0.7
    pmix <- 0.5
    if (fv_dist == 1) {
      (1 - fv * log((pmix * exp(-lambda1 * t^gamma1) + (1 - pmix) * exp(-lambda2 * t^gamma2))^exp(beta * X)))^(-1 / fv)
    } else {
      S <- ((pmix * exp(-lambda1 * t^gamma1) + (1 - pmix) * exp(-lambda2 * t^gamma2))^exp(beta * X))
      fn <- function(eta) S^exp(eta) * dnorm(eta, mean = 0, sd = sqrt(fv))
      quadinf(f = fn, xa = -Inf, xb = Inf)$Q
    }
  }
})

# Make LLE
lle <- distinct(dgms, fv, fv_dist, baseline)
lle$lle <- sapply(1:nrow(lle), function(i) {
  cat(i, "\n")
  t <- seq(0, 5, length.out = 1000)
  ms1 <- with(lle, true_ms(t = t, X = 1, distribution = baseline[i], fv_dist = fv_dist[i], fv = fv[i]))
  ms0 <- with(lle, true_ms(t = t, X = 0, distribution = baseline[i], fv_dist = fv_dist[i], fv = fv[i]))
  ms1f <- splinefun(x = t, y = ms1)
  ms0f <- splinefun(x = t, y = ms0)
  is1 <- quadinf(f = ms1f, xa = 0, xb = 5)$Q
  is0 <- quadinf(f = ms0f, xa = 0, xb = 5)$Q
  is1 - is0
})

# Merge back
dgms <- dgms %>%
  left_join(lle, c("fv", "fv_dist", "baseline"))

# Make data.frame with true values
true <- crossing(
  scenario = 1:120,
  par = c("fv", "lle", "trt")
) %>%
  left_join(distinct(dgms, scenario, fv) %>%
    mutate(par = "fv"),
  by = c("scenario", "par")
  ) %>%
  left_join(distinct(dgms, scenario, lle) %>%
    mutate(par = "lle"),
  by = c("scenario", "par")
  ) %>%
  left_join(distinct(dgms, scenario, treatment_effect) %>%
    mutate(par = "trt"),
  by = c("scenario", "par")
  ) %>%
  mutate(true = case_when(
    par == "fv" ~ fv,
    par == "lle" ~ lle,
    par == "trt" ~ treatment_effect
  )) %>%
  select(scenario, par, true)

# Call simsum on each object of the list to summarise the results
lofsummary <- lapply(seq_along(lofresults), function(i) {
  tv <- left_join(lofresults[[i]],
    true,
    by = c("scenario", "par")
  ) %>%
    distinct(true) %>%
    pull()
  s <- simsum(lofresults[[i]], true = tv, estvarname = "coef", se = "se", methodvar = "model", by = "scenario", ref = "1")
  smr <- summary(s)
  out <- get_data(smr) %>%
    filter(stat %in% c("nsim", "thetamean", "thetamedian", "se2mean", "se2median", "bias", "empse", "mse", "cover"))
  out[["par"]] <- unique(lofresults[[i]][["par"]])
  out
})

# Bind each scenario in a data.frame
summary <- bind_rows(lofsummary) %>%
  mutate(
    model = as.numeric(model),
    scenario = as.numeric(scenario)
  ) %>%
  arrange(stat, scenario, model)

# Factorise
summary <- summary %>%
  left_join(dgms, by = "scenario") %>%
  select(-pars) %>%
  mutate(
    model = factor(model, levels = 1:22, labels = c("Cox, Gamma", "Cox, Normal", "Exp, Gamma", "Exp, Normal", "Weibull, Gamma", "Weibull, Normal", "Gompertz, Gamma", "Gompertz, Normal", "RP (3), Gamma", "RP (3), Normal", "RP (5), Gamma", "RP (5), Normal", "RP (9), Gamma", "RP (9), Normal", "RP (P), Gamma", "RP (P), Normal", "FP (W), Gamma", "FP (W), Normal", "FP (k=10), Gamma", "FP (k=10), Normal", "FP (k=10000), Gamma", "FP (k=10000), Normal")),
    fv = factor(fv, levels = c(0.25, 0.75, 1.25)),
    fv_dist = factor(fv_dist, levels = 1:2, labels = c("Gamma", "Normal")),
    baseline = factor(baseline, levels = 1:5, labels = c("Exponential", "Weibull", "Gompertz", "Weibull-Weibull (1)", "Weibull-Weibull (2)"))
  ) %>%
  separate(model, into = c("mbaseline", "mfrailty"), sep = ", ", remove = FALSE) %>%
  mutate(
    mbaseline = factor(mbaseline, levels = c("Cox", "Exp", "Weibull", "Gompertz", "RP (3)", "RP (5)", "RP (9)", "RP (P)", "FP (W)", "FP (k=10)", "FP (k=10000)")),
    mfrailty = factor(mfrailty, levels = c("Gamma", "Normal"))
  ) %>%
  mutate(par = factor(par, levels = c("trt", "lle", "fv"), labels = c("Treatment effect", "LLE", "Frailty variance"))) %>%
  mutate(sig = case_when(
    stat == "bias" & (lower > 0 | upper < 0) ~ "*",
    stat == "cover" & (lower > 0.95 | upper < 0.95) ~ "*",
    TRUE ~ ""
  )) %>%
  arrange(scenario, par, stat, model)

# Export results
saveRDS(summary, "summary.RDS")
