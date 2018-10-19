### Function to estimate all models of interest for a given scenario and produce absolute risk estimates

fcoxme <- function(data, B = 1000, ptime, null_obj, progress = FALSE) {
  # fit model
  fit <- tryCatch(
    coxme::coxme(Surv(t, d) ~ trt + (1 | grpid), data = data),
    error = function(e)
      return(null_obj)
  )

  # break out of the function if the fit did not converge
  if ("null_obj" %in% class(fit)) {
    return(null_obj)
  }

  # boostrap SE for the frailty variance
  # NB: re-sample the clusters!!! (silly you)
  # Set up progress bar if required
  if (progress) {
    cat("\nBootstrap:\n")
    pb <- txtProgressBar(min = 0, max = B, style = 3)
  }
  b <- vapply(
    X = 1:B,
    FUN = function(bbb) {
      boot_fit <- NULL
      # putting each bootstrap fit in a while loop:
      # the fit function was failing at times, should that happen the fit is repeated with a different bootstrap sample
      while (is.null(boot_fit)) {
        i <- sample(
          x = unique(data[["grpid"]]),
          size = length(unique(data[["grpid"]])),
          replace = TRUE
        )
        out <- lapply(seq_along(i), function(u) {
          d <- data[data[["grpid"]] == i[u], ]
          d[["grpid"]] <- u
          d
        })
        boot_data <- do.call("rbind.data.frame", out) %>%
          arrange(grpid)
        boot_fit <- tryCatch(
          coxme::coxme(Surv(t, d) ~ trt + (1 | grpid), data = boot_data),
          error = function(e)
            return(NULL)
        )
      }
      boot_theta <- VarCorr(boot_fit)[["grpid"]][["Intercept"]]
      # Update progress bar
      if (progress) setTxtProgressBar(pb, bbb)
      return(boot_theta)
    },
    numeric(1)
  )

  # Close progress bar
  if (progress) close(pb)

  llef <- function(object, newdata) {
    # extract frailty variance and regression coefficient
    beta <- coef(object)[["trt"]]
    theta <- VarCorr(object)[["grpid"]][["Intercept"]]

    # make breslow estimates of H0 for all possible times
    # this is less general, but waaaay quicker for my application as I know all the possible t in advance
    ti <- sort(unique(c(data[["t"]], ptime)))
    h0i <- vapply(
      ti, function(x) {
        di <- sum(data[["d"]][data[["t"]] == x])
        deni <- sum(exp(object[["linear.predictor"]][data[["t"]] >= x]))
        return(di / deni)
      },
      numeric(1)
    )
    H0i <- cumsum(h0i)

    ms <- function(t, X = 0) {
      S <- (exp(-H0i[ti == t])^exp(beta * X))
      fn <-
        function(eta)
          S^exp(eta) * stats::dnorm(eta, mean = 0, sd = sqrt(theta))
      pracma::quadinf(f = fn, xa = -Inf, xb = Inf)[["Q"]]
    }

    # obtain the integrated difference in survival curves
    pr <- newdata
    pr[["s"]] <-
      vapply(
        X = 1:nrow(pr),
        FUN = function(i)
          ms(t = pr[["t"]][i], X = pr[["trt"]][i]),
        FUN.VALUE = numeric(1)
      )

    nsf_0 <-
      stats::splinefun(ptime, pr[["s"]][pr[["trt"]] == 0], method = "natural")
    nsf_1 <-
      stats::splinefun(ptime, pr[["s"]][pr[["trt"]] == 1], method = "natural")
    ims_0 <- pracma::quadinf(f = nsf_0, xa = 0, xb = 5)[["Q"]]
    ims_1 <- pracma::quadinf(f = nsf_1, xa = 0, xb = 5)[["Q"]]
    return(ims_1 - ims_0)
  }

  lle <- rstpm2::predictnl(object = fit, fun = llef, newdata = expand.grid(t = ptime, trt = 0:1))

  # make object to return
  out <- data.frame(
    par = c("trt", "fv", "lle"),
    coef = c(coef(fit)[["trt"]], VarCorr(fit)[["grpid"]][["Intercept"]], lle[["fit"]]),
    se = c(sqrt(vcov(fit))[1, 1], sqrt(var(b)), lle[["se.fit"]])
  )
  return(out)
}

fstpm2 <- function(data, df, RandDist, ptime, null_obj) {
  # fit model
  fit <- tryCatch(
    rstpm2::stpm2(
      Surv(t, d) ~ trt,
      data = data,
      cluster = data[["grpid"]],
      RandDist = RandDist,
      df = df
    ),
    error = function(e)
      return(null_obj)
  )

  # break out of the function if the fit did not converge
  if ("null_obj" %in% class(fit)) {
    return(null_obj)
  }

  # make summary object (has frailty variance and not log-theta)
  smr <- tryCatch(
    summary(fit),
    error = function(e)
      return(null_obj)
  )

  # break out of the function if the fit converged to bad values (e.g. hessian not invertible)
  if ("null_obj" %in% class(smr)) {
    return(null_obj)
  }

  # if the object has NaN or NA values, break out and return the null object
  if (any(is.nan(smr@coef[, "Estimate"]) | is.na(smr@coef[, "Estimate"]))) {
    return(null_obj)
  }
  if (any(is.nan(smr@coef[, "Std. Error"]) | is.na(smr@coef[, "Std. Error"]))) {
    return(null_obj)
  }
  if (is.nan(smr@theta[["theta"]]) | is.na(smr@theta[["theta"]])) {
    return(null_obj)
  }
  if (is.nan(smr@theta[["se.theta"]]) | is.na(smr@theta[["se.theta"]])) {
    return(null_obj)
  }

  # make integrated difference of survival curves
  llef <- function(object, newdata) {
    s <-
      predict(object, newdata, type = "margsurv")
    pr <- cbind(newdata, s)
    nsf_0 <-
      stats::splinefun(ptime, pr[["s"]][pr[["trt"]] == 0], method = "natural")
    nsf_1 <-
      stats::splinefun(ptime, pr[["s"]][pr[["trt"]] == 1], method = "natural")
    ims_0 <- pracma::quadinf(f = nsf_0, xa = 0, xb = 5)[["Q"]]
    ims_1 <- pracma::quadinf(f = nsf_1, xa = 0, xb = 5)[["Q"]]
    return(ims_1 - ims_0)
  }
  lle <- rstpm2::predictnl(object = fit, fun = llef, newdata = expand.grid(t = ptime, trt = 0:1))

  # make object to return
  out <- data.frame(
    par = c("trt", "fv", "lle"),
    coef = c(smr@coef["trt", "Estimate"], smr@theta[["theta"]], lle[["fit"]]),
    se = c(smr@coef["trt", "Std. Error"], smr@theta[["se.theta"]], lle[["se.fit"]])
  )
  return(out)
}

fpstpm2 <- function(data, RandDist, ptime, null_obj) {
  # fit model
  fit <- tryCatch(
    rstpm2::pstpm2(
      Surv(t, d) ~ trt,
      data = data,
      cluster = data[["grpid"]],
      RandDist = RandDist
    ),
    error = function(e)
      return(null_obj)
  )

  # break out of the function if the fit did not converge
  if ("null_obj" %in% class(fit)) {
    return(null_obj)
  }

  # make summary object (has frailty variance and not log-theta)
  smr <- tryCatch(
    summary(fit),
    error = function(e)
      return(null_obj)
  )

  # break out of the function if the fit converged to bad values (e.g. hessian not invertible)
  if ("null_obj" %in% class(smr)) {
    return(null_obj)
  }

  # if the object has NaN or NA values, break out and return the null object
  if (any(is.nan(smr@coef[, "Estimate"]) | is.na(smr@coef[, "Estimate"]))) {
    return(null_obj)
  }
  if (any(is.nan(smr@coef[, "Std. Error"]) | is.na(smr@coef[, "Std. Error"]))) {
    return(null_obj)
  }
  if (is.nan(smr@theta[["theta"]]) | is.na(smr@theta[["theta"]])) {
    return(null_obj)
  }
  if (is.nan(smr@theta[["se.theta"]]) | is.na(smr@theta[["se.theta"]])) {
    return(null_obj)
  }

  # make integrated difference of survival curves
  llef <- function(object, newdata) {
    s <- predict(object, newdata = newdata, type = "margsurv")
    pr <- cbind(newdata, s)
    nsf_0 <- stats::splinefun(ptime, pr[["s"]][pr[["trt"]] == 0], method = "natural")
    nsf_1 <- stats::splinefun(ptime, pr[["s"]][pr[["trt"]] == 1], method = "natural")
    ims_0 <- pracma::quadinf(f = nsf_0, xa = 0, xb = 5)[["Q"]]
    ims_1 <- pracma::quadinf(f = nsf_1, xa = 0, xb = 5)[["Q"]]
    return(ims_1 - ims_0)
  }
  lle <- rstpm2::predictnl(object = fit, fun = llef, newdata = expand.grid(t = ptime, trt = 0:1))

  # make object to return
  out <- data.frame(
    par = c("trt", "fv", "lle"),
    coef = c(smr@coef["trt", "Estimate"], smr@theta[["theta"]], lle[["fit"]]),
    se = c(smr@coef["trt", "Std. Error"], smr@theta[["se.theta"]], lle[["se.fit"]])
  )
  return(out)
}

fparfm <- function(data, dist, frailty, ptime, null_obj) {
  # fit model
  fit <- tryCatch(
    parfm::parfm(
      Surv(t, d) ~ trt,
      cluster = "grpid",
      data = data,
      dist = dist,
      frailty = frailty
    ),
    error = function(e)
      return(null_obj),
    warning = function(w)
      return(null_obj)
  )

  # there is some weird bug in parfm...
  sink()

  # break out of the function if the fit did not converge
  if ("null_obj" %in% class(fit)) {
    return(null_obj)
  }

  llef <- function(object, newdata) {
    # extract frailty variance and regression coefficient
    theta <-
      ifelse(attr(object, "frailty") == "gamma", coef(object)[["theta"]], coef(object)[["sigma2"]])
    beta <- coef(object)[["trt"]]

    # define function for computing marginal survival. if frailty = gamma using analytical formulae, numerical integration otherwise
    ms <- switch(
      attr(object, "dist"),
      "exponential" = if (attr(object, "frailty") == "gamma") {
        function(t, X = 0)
          (1 - theta * (-coef(object)[["lambda"]] * exp(beta * X) * t))^(-1 / theta)
      } else {
        function(t, X = 0) {
          S <- (exp(-coef(object)[["lambda"]] * exp(beta * X) * t))
          fn <-
            function(eta)
              S^exp(eta) * stats::dnorm(eta, mean = 0, sd = sqrt(theta))
          pracma::quadinf(
            f = fn,
            xa = -Inf,
            xb = Inf
          )[["Q"]]
        }
      },
      "weibull" = if (attr(object, "frailty") == "gamma") {
        function(t, X = 0)
          (1 - theta * (-coef(object)[["lambda"]] * exp(beta * X) * t^coef(object)[["rho"]]))^(-1 / theta)
      } else {
        function(t, X = 0) {
          S <-
            (exp(-coef(object)[["lambda"]] * exp(beta * X) * t^coef(object)[["rho"]]))
          fn <-
            function(eta)
              S^exp(eta) * stats::dnorm(eta, mean = 0, sd = sqrt(theta))
          pracma::quadinf(
            f = fn,
            xa = -Inf,
            xb = Inf
          )[["Q"]]
        }
      },
      "gompertz" = if (attr(object, "frailty") == "gamma") {
        function(t, X = 0)
          (1 - theta * ((-coef(object)[["lambda"]] * exp(beta * X) / coef(object)[["gamma"]]) * (exp(coef(object)[["gamma"]] * t) - 1)))^(-1 / theta)
      } else {
        function(t, X = 0) {
          S <-
            (exp((-coef(object)[["lambda"]] * exp(beta * X) / coef(object)[["gamma"]]) * (exp(coef(object)[["gamma"]] * t) - 1)))
          fn <-
            function(eta)
              S^exp(eta) * stats::dnorm(eta, mean = 0, sd = sqrt(theta))
          pracma::quadinf(
            f = fn,
            xa = -Inf,
            xb = Inf
          )[["Q"]]
        }
      }
    )

    # make predictions for integrated survival difference
    pr <- newdata
    pr[["s"]] <-
      vapply(
        X = 1:nrow(pr),
        FUN = function(i)
          ms(t = pr[["t"]][i], X = pr[["trt"]][i]),
        FUN.VALUE = numeric(1)
      )
    nsf_0 <-
      stats::splinefun(ptime, pr[["s"]][pr[["trt"]] == 0], method = "natural")
    nsf_1 <-
      stats::splinefun(ptime, pr[["s"]][pr[["trt"]] == 1], method = "natural")
    ims_0 <- pracma::quadinf(f = nsf_0, xa = 0, xb = 5)[["Q"]]
    ims_1 <- pracma::quadinf(f = nsf_1, xa = 0, xb = 5)[["Q"]]
    return(ims_1 - ims_0)
  }

  pnl <- function(object = fit, fun = llef, newdata = expand.grid(t = ptime, trt = 0:1)) {
    local1 <- function(coef, newdata, ...) {
      coef(object) <- coef
      fun(object, newdata = newdata, ...)
    }
    coef <- coef(object)
    Sigma <- vcov(object)
    est <- local1(coef, newdata = newdata)
    gd <- rstpm2:::grad(func = local1, x = coef, newdata = newdata)
    se.fit <- as.vector(sqrt(colSums(gd * (Sigma %*% gd))))
    data.frame(fit = est, se.fit = se.fit)
  }
  lle <- pnl()

  # make object to return
  out <- data.frame(
    par = c("trt", "fv", "lle"),
    coef = c(coef(fit)[["trt"]], ifelse(attr(fit, "frailty") == "gamma", coef(fit)[["theta"]], coef(fit)[["sigma2"]]), lle[["fit"]]),
    se = c(fit["trt", "SE"], ifelse(
      attr(fit, "frailty") == "gamma", fit["theta", "SE"], fit["sigma2", "SE"]
    ), lle[["se.fit"]])
  )
  return(out)
}

ffrailtyEM <- function(data, ptime, null_obj) {
  # fit model
  cc <- frailtyEM::emfrail_control(se_adj = FALSE, ca_test = FALSE)
  fit <- tryCatch(
    frailtyEM::emfrail(
      Surv(t, d) ~ trt + cluster(grpid),
      distribution = emfrail_dist(dist = "gamma"),
      data = data,
      control = cc
    ),
    error = function(e)
      return(null_obj)
  )

  # break out of the function if the fit did not converge
  if ("null_obj" %in% class(fit)) {
    return(null_obj)
  }

  # make summary object
  smr <- tryCatch(
    summary(fit),
    error = function(e)
      return(null_obj)
  )

  # break out of the function if the fit converged to bad values (e.g. hessian not invertible)
  if ("null_obj" %in% class(smr)) {
    return(null_obj)
  }

  # if the object has NaN values, break out and return the null object
  if (any(is.nan(smr[["coefmat"]][, "coef"]))) {
    return(null_obj)
  }
  if (any(is.nan(smr[["coefmat"]][, "se(coef)"]))) {
    return(null_obj)
  }
  if (any(is.nan(smr[["fr_var"]]["fr_var"]))) {
    return(null_obj)
  }
  if (any(is.nan(smr[["fr_var"]]["se_fr_var"]))) {
    return(null_obj)
  }

  llef <- function(object, newdata) {
    # make predictions
    preds <- predict(
      object,
      newdata = expand.grid(trt = 0:1),
      type = "marginal",
      quantity = "survival",
      conf_int = NULL
    )
    preds <- dplyr::bind_rows(
      dplyr::mutate(preds[[1]], trt = 0),
      dplyr::mutate(preds[[2]], trt = 1),
      newdata
    ) %>%
      dplyr::distinct(time, trt, .keep_all = TRUE) %>%
      dplyr::arrange(trt, time) %>%
      dplyr::group_by(trt) %>%
      dplyr::mutate(survival_m = ifelse(is.na(survival_m) &
        time == dplyr::first(time), 1, survival_m)) %>%
      dplyr::mutate(survival_m = zoo::na.locf(survival_m)) %>%
      dplyr::ungroup()
    preds[["lp"]] <- NULL

    # obtain the integrated difference in survival curves
    pr <- dplyr::filter(preds, time %in% ptime)
    nsf_0 <-
      stats::splinefun(ptime, pr[["survival_m"]][pr[["trt"]] == 0], method = "natural")
    nsf_1 <-
      stats::splinefun(ptime, pr[["survival_m"]][pr[["trt"]] == 1], method = "natural")
    ims_0 <- pracma::quadinf(f = nsf_0, xa = 0, xb = 5)[["Q"]]
    ims_1 <- pracma::quadinf(f = nsf_1, xa = 0, xb = 5)[["Q"]]
    return(ims_1 - ims_0)
  }

  pnl <- function(object = fit, fun = llef, newdata = expand.grid(time = ptime, trt = 0:1)) {
    local1 <- function(coef, newdata, ...) {
      coef(object) <- coef
      fun(object, newdata = newdata, ...)
    }
    coef <- coef(object)
    Sigma <- vcov(object)
    est <- local1(coef, newdata = newdata)
    gd <- rstpm2:::grad(func = local1, x = coef, newdata = newdata)
    se.fit <- as.vector(sqrt(colSums(gd * (Sigma %*% gd))))
    data.frame(fit = est, se.fit = se.fit)
  }
  lle <- pnl()

  # make object to return
  out <- data.frame(
    par = c("trt", "fv", "lle"),
    coef = c(coef(fit)["trt"][[1]], smr[["fr_var"]]["fr_var"][[1]], lle[["fit"]]),
    se = c(smr[["coefmat"]]["trt", "se(coef)"], smr[["fr_var"]]["se_fr_var"][[1]], lle[["se.fit"]])
  )
  return(out)
}

ffrailtypack <- function(data, hazard, kappa = NULL, RandDist, ptime, null_obj) {
  # fit model
  if (hazard == "Weibull") {
    fit <- tryCatch(
      frailtypack::frailtyPenal(
        Surv(t, d) ~ trt + cluster(grpid),
        data = data,
        RandDist = RandDist,
        hazard = "Weibull",
        print.times = FALSE
      ),
      error = function(e)
        return(null_obj),
      warning = function(w)
        return(null_obj)
    )
  } else {
    fit <- tryCatch(
      frailtypack::frailtyPenal(
        Surv(t, d) ~ trt + cluster(grpid),
        data = data,
        RandDist = RandDist,
        hazard = "Splines",
        kappa = kappa,
        n.knots = 12,
        print.times = FALSE
      ),
      error = function(e)
        return(null_obj),
      warning = function(w)
        return(null_obj)
    )
  }

  # break out of the function if the fit did not converge
  if ("null_obj" %in% class(fit)) {
    return(null_obj)
  }

  llef <- function(object) {
    predata <- data.frame(trt = 0:1, grpid = 1)
    # make integrated difference in survival curves
    probj <-
      frailtypack::prediction(
        fit = object,
        data = predata,
        t = 1e-10,
        window = (ptime - 2 * 1e-10)
      )
    pr <- data.frame(
      t = rep(ptime, 2),
      s = 1 - c(probj[["pred"]][1, ], probj[["pred"]][2, ]),
      trt = c(rep(predata[["trt"]][1], length(ptime)), rep(predata[["trt"]][2], length(ptime)))
    )
    nsf_0 <-
      stats::splinefun(ptime, pr[["s"]][pr[["trt"]] == 0], method = "natural")
    nsf_1 <-
      stats::splinefun(ptime, pr[["s"]][pr[["trt"]] == 1], method = "natural")
    ims_0 <- pracma::quadinf(f = nsf_0, xa = 0, xb = 5)[["Q"]]
    ims_1 <- pracma::quadinf(f = nsf_1, xa = 0, xb = 5)[["Q"]]
    return(ims_1 - ims_0)
  }

  sink(tempfile())
  lle <- rstpm2::predictnl.default(object = fit, fun = llef)
  sink()

  # make object to return
  out <- data.frame(
    par = c("trt", "fv", "lle"),
    coef = c(fit[["coef"]]["trt"], ifelse(RandDist == "Gamma", fit[["theta"]], fit[["sigma2"]]), lle[["fit"]]),
    se = c(sqrt(fit[["varH"]]), sqrt(fit[["varTheta"]][1]), lle[["se.fit"]])
  )
  return(out)
}

### Custom methods for coef, coef<-, and vcov
### Required by predictnl
coef.parfm <- function(object) {
  object[, "ESTIMATE"]
}

`coef<-.parfm` <- function(x, value) {
  x[, "ESTIMATE"] <- value
  x
}

vcov.parfm <- function(object, ...) {
  x <- inv(attr(object, "FisherI"))
  rownames(x)[grepl("trt", rownames(x))] <- "trt"
  colnames(x)[grepl("trt", colnames(x))] <- "trt"
  x
}

`coef<-.emfrail` <- function(x, value) {
  x$coef <- value
  x
}

vcov.emfrail <- function(object, ...) {
  frailtyEM:::vcov.emfrail(object, type = "regular")[seq_along(object$coef), seq_along(object$coef)]
}

coef.frailtyPenal <- function(object) {
  object$coef
}

`coef<-.frailtyPenal` <- function(x, value) {
  x$coef <- value
  x
}

vcov.frailtyPenal <- function(object, ...) {
  object$varH
}

### Function to fit a given model
make_models <- function(data, i, ptime, model) {
  # Null object to return in case of errors:
  null_obj <- data.frame(
    par = c("trt", "fv", "lle"),
    coef = rep(NA, 3),
    se = rep(NA, 3)
  )
  class(null_obj) <- c(class(null_obj), "null_obj")

  # Required packages:
  if (!requireNamespace("pacman")) {
    install.packages("pacman")
  }
  pacman::p_load(
    "dplyr",
    "survival",
    "rstpm2",
    "parfm",
    "tidyr",
    "zoo",
    "frailtyEM",
    "coxme",
    "pracma",
    "frailtypack"
  )

  tstart <- Sys.time()

  # Models:
  out <- if (model == 1) {
    ffrailtyEM(
      data = data,
      ptime = ptime,
      null_obj = null_obj
    )
  } else if (model == 2) {
    fcoxme(
      data = data,
      ptime = ptime,
      null_obj = null_obj
    )
  } else if (model == 3) {
    fparfm(
      data = data,
      dist = "exponential",
      frailty = "gamma",
      ptime = ptime,
      null_obj = null_obj
    )
  } else if (model == 4) {
    fparfm(
      data = data,
      dist = "exponential",
      frailty = "lognormal",
      ptime = ptime,
      null_obj = null_obj
    )
  } else if (model == 5) {
    fparfm(
      data = data,
      dist = "weibull",
      frailty = "gamma",
      ptime = ptime,
      null_obj = null_obj
    )
  } else if (model == 6) {
    fparfm(
      data = data,
      dist = "weibull",
      frailty = "lognormal",
      ptime = ptime,
      null_obj = null_obj
    )
  } else if (model == 7) {
    fparfm(
      data = data,
      dist = "gompertz",
      frailty = "gamma",
      ptime = ptime,
      null_obj = null_obj
    )
  } else if (model == 8) {
    fparfm(
      data = data,
      dist = "gompertz",
      frailty = "lognormal",
      ptime = ptime,
      null_obj = null_obj
    )
  } else if (model == 9) {
    fstpm2(
      data = data,
      df = 3,
      RandDist = "Gamma",
      ptime = ptime,
      null_obj = null_obj
    )
  } else if (model == 10) {
    fstpm2(
      data = data,
      df = 3,
      RandDist = "LogN",
      ptime = ptime,
      null_obj = null_obj
    )
  } else if (model == 11) {
    fstpm2(
      data = data,
      df = 5,
      RandDist = "Gamma",
      ptime = ptime,
      null_obj = null_obj
    )
  } else if (model == 12) {
    fstpm2(
      data = data,
      df = 5,
      RandDist = "LogN",
      ptime = ptime,
      null_obj = null_obj
    )
  } else if (model == 13) {
    fstpm2(
      data = data,
      df = 9,
      RandDist = "Gamma",
      ptime = ptime,
      null_obj = null_obj
    )
  } else if (model == 14) {
    fstpm2(
      data = data,
      df = 9,
      RandDist = "LogN",
      ptime = ptime,
      null_obj = null_obj
    )
  } else if (model == 15) {
    fpstpm2(
      data = data,
      RandDist = "Gamma",
      ptime = ptime,
      null_obj = null_obj
    )
  } else if (model == 16) {
    fpstpm2(
      data = data,
      RandDist = "LogN",
      ptime = ptime,
      null_obj = null_obj
    )
  } else if (model == 17) {
    ffrailtypack(
      data = data,
      hazard = "Weibull",
      RandDist = "Gamma",
      ptime = ptime,
      null_obj = null_obj
    )
  } else if (model == 18) {
    ffrailtypack(
      data = data,
      hazard = "Weibull",
      RandDist = "LogN",
      ptime = ptime,
      null_obj = null_obj
    )
  } else if (model == 19) {
    ffrailtypack(
      data = data,
      hazard = "Splines",
      kappa = 10,
      RandDist = "Gamma",
      ptime = ptime,
      null_obj = null_obj
    )
  } else if (model == 20) {
    ffrailtypack(
      data = data,
      hazard = "Splines",
      kappa = 10,
      RandDist = "LogN",
      ptime = ptime,
      null_obj = null_obj
    )
  } else if (model == 21) {
    ffrailtypack(
      data = data,
      hazard = "Splines",
      kappa = 10000,
      RandDist = "Gamma",
      ptime = ptime,
      null_obj = null_obj
    )
  } else if (model == 22) {
    ffrailtypack(
      data = data,
      hazard = "Splines",
      kappa = 10000,
      RandDist = "LogN",
      ptime = ptime,
      null_obj = null_obj
    )
  }

  tstop <- Sys.time()

  # Add some info to return
  out[["time"]] <- difftime(tstop, tstart, units = "mins")
  out[["model"]] <- model
  out[["i"]] <- i
  out[["scenario"]] <- unique(data[["scenario"]])

  # Return results
  return(out)
}
