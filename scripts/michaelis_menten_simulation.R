## =========================================================
## michaelis_menten_simulation.R
##
## Streamlined manuscript simulation script for the
## heteroscedastic Michaelis-Menten paper.
##
## This script keeps only the functions needed to:
##   1. run the manuscript simulation,
##   2. save raw/summary/reference outputs,
##   3. export 12 scenario-specific manuscript tables, and
##   4. export 4 manuscript-ready tables grouped by
##      error distribution (Normal / Skewed).
##
## Default output file names remain compatible with the
## current manuscript build pipeline.
## =========================================================

V_TRUE <- 100
K_TRUE <- 20
X_DEFAULT <- seq(1, 100, length.out = 100)

MANUSCRIPT_ERROR_CASES <- c("Normal", "Skewed")
MANUSCRIPT_MODEL_NAMES <- c(
  "NLS",
  "Semi_log",
  "Semi_x_half",
  "Semi_x_third",
  "Bayes_logvar"
)
MANUSCRIPT_FILE_PREFIX <- "simulation_real_manuscript_final"
MANUSCRIPT_SCENARIO_ORDER <- c("MM", "Exp", "Hill")

manuscript_model_names <- function(use_bayes = TRUE) {
  keep_models <- MANUSCRIPT_MODEL_NAMES

  if (!isTRUE(use_bayes)) {
    keep_models <- setdiff(keep_models, "Bayes_logvar")
  }

  keep_models
}

## =========================================================
## 1. Shared helpers
## =========================================================

mm_mean <- function(x, vmax, km) {
  vmax * x / (km + x)
}

mm_grad <- function(x, vmax, km) {
  cbind(
    Vmax = x / (km + x),
    Km = -(vmax * x) / (km + x)^2
  )
}

safe_pos <- function(x, lower = 1e-10, upper = 1e10) {
  x[!is.finite(x)] <- upper
  pmin(pmax(x, lower), upper)
}

safe_inverse <- function(mat, ridge = 1e-8, max_tries = 6L) {
  p <- nrow(mat)
  eye <- diag(p)

  for (j in 0:max_tries) {
    out <- tryCatch(
      solve(mat + ridge * 10^j * eye),
      error = function(e) NULL
    )
    if (!is.null(out) && all(is.finite(out))) {
      return(out)
    }
  }

  stop("Unable to invert matrix.")
}

get_wald_interval <- function(beta, vcov_beta, level = 0.95) {
  alpha <- 1 - level
  z_val <- stats::qnorm(1 - alpha / 2)
  se <- sqrt(pmax(diag(vcov_beta), 0))

  list(
    lower = beta - z_val * se,
    upper = beta + z_val * se
  )
}

get_posterior_interval <- function(draws, level = 0.95) {
  alpha <- 1 - level
  q <- stats::quantile(
    draws,
    probs = c(alpha / 2, 1 - alpha / 2),
    names = FALSE,
    na.rm = TRUE
  )

  c(lower = unname(q[1]), upper = unname(q[2]))
}

round_numeric_df <- function(x, digits = 5) {
  numeric_cols <- vapply(x, is.numeric, logical(1))
  x[numeric_cols] <- lapply(x[numeric_cols], round, digits = digits)
  x
}

interval_score <- function(lower, upper, truth, level = 0.95) {
  alpha <- 1 - level
  width <- upper - lower
  below <- as.numeric(truth < lower)
  above <- as.numeric(truth > upper)

  width +
    (2 / alpha) * (lower - truth) * below +
    (2 / alpha) * (truth - upper) * above
}

safe_mean <- function(x) {
  if (!length(x)) {
    return(NA_real_)
  }

  x <- x[is.finite(x)]
  if (!length(x)) {
    return(NA_real_)
  }

  mean(x)
}

se_calibration_ratio <- function(estimates, standard_errors) {
  ok <- is.finite(estimates) & is.finite(standard_errors) & standard_errors > 0
  if (sum(ok) < 2L) {
    return(NA_real_)
  }

  mc_sd <- stats::sd(estimates[ok])
  mean_se <- mean(standard_errors[ok])

  if (!is.finite(mc_sd) || !is.finite(mean_se) || mean_se <= 0) {
    return(NA_real_)
  }

  mc_sd / mean_se
}

scenario_from_dgp <- function(dgp_name) {
  if (grepl("^DGP_MM_", dgp_name)) {
    return("MM")
  }
  if (grepl("^DGP_Exp_", dgp_name)) {
    return("Exp")
  }
  if (grepl("^DGP_Hill_", dgp_name)) {
    return("Hill")
  }
  NA_character_
}

error_from_dgp <- function(dgp_name) {
  if (grepl("_Normal$", dgp_name)) {
    return("Normal")
  }
  if (grepl("_Skewed$", dgp_name)) {
    return("Skewed")
  }
  NA_character_
}

sort_manuscript_rows <- function(dat) {
  dat[
    order(
      match(dat$Scenario, MANUSCRIPT_SCENARIO_ORDER),
      match(dat$Model, MANUSCRIPT_MODEL_NAMES)
    ),
    ,
    drop = FALSE
  ]
}

## =========================================================
## 2. Variance structures and working weights
## =========================================================

var_mm <- function(x, s0 = 1, t2 = 9, c = 20) {
  s0 + t2 * x / (c + x)
}

var_exp <- function(x, s0 = 1, t2 = 9, a = 0.05) {
  s0 + t2 * (1 - exp(-a * x))
}

var_hill <- function(x, s0 = 1, t2 = 9, c = 20, alpha = 2) {
  s0 + t2 * x^alpha / (c^alpha + x^alpha)
}

h_log <- function(x) {
  log(x + 1)
}

h_pow <- function(a) {
  force(a)
  function(x) x^a
}

## =========================================================
## 3. Error generation and DGP library
## =========================================================

draw_standardized_error <- function(n,
                                    error_case = MANUSCRIPT_ERROR_CASES,
                                    skew_shape = 4) {
  error_case <- match.arg(error_case, MANUSCRIPT_ERROR_CASES)

  if (error_case == "Normal") {
    return(stats::rnorm(n))
  }

  z <- stats::rgamma(n, shape = skew_shape, scale = 1)
  (z - skew_shape) / sqrt(skew_shape)
}

make_dgp_library <- function(error_cases = MANUSCRIPT_ERROR_CASES,
                             s0 = 1,
                             t2 = 9,
                             c = 20,
                             a = 0.05,
                             alpha = 2,
                             skew_shape = 4) {
  error_cases <- match.arg(
    error_cases,
    MANUSCRIPT_ERROR_CASES,
    several.ok = TRUE
  )

  variance_specs <- list(
    DGP_MM = list(
      label = "MM bounded variance",
      variance = function(x) var_mm(x, s0 = s0, t2 = t2, c = c)
    ),
    DGP_Exp = list(
      label = "Exponential bounded variance",
      variance = function(x) var_exp(x, s0 = s0, t2 = t2, a = a)
    ),
    DGP_Hill = list(
      label = "Hill bounded variance",
      variance = function(x) var_hill(x, s0 = s0, t2 = t2, c = c, alpha = alpha)
    )
  )

  out <- list()
  for (v_name in names(variance_specs)) {
    for (err in error_cases) {
      dgp_name <- paste(v_name, err, sep = "_")
      spec <- variance_specs[[v_name]]
      out[[dgp_name]] <- list(
        label = paste(spec$label, err, "errors"),
        error_case = err,
        skew_shape = skew_shape,
        variance = spec$variance
      )
    }
  }

  out
}

simulate_response <- function(x, dgp_spec,
                              vmax_true = V_TRUE,
                              km_true = K_TRUE) {
  v_true <- dgp_spec$variance(x)
  z <- draw_standardized_error(
    n = length(x),
    error_case = dgp_spec$error_case,
    skew_shape = dgp_spec$skew_shape
  )

  list(
    y = mm_mean(x, vmax_true, km_true) + sqrt(v_true) * z,
    v_true = v_true,
    z = z
  )
}

## =========================================================
## 4. Weighted profile estimators
## =========================================================

weighted_score_function <- function(km, x, y, h_vals) {
  w <- 1 / safe_pos(h_vals)
  left_1 <- sum(w * x * y / (km + x))
  right_1 <- sum(w * x^2 / (km + x)^3)
  left_2 <- sum(w * x * y / (km + x)^2)
  right_2 <- sum(w * x^2 / (km + x)^2)
  left_1 * right_1 - left_2 * right_2
}

weighted_vmax_hat <- function(km, x, y, h_vals) {
  w <- 1 / safe_pos(h_vals)
  num <- sum(w * x * y / (km + x))
  den <- sum(w * x^2 / (km + x)^2)
  num / den
}

find_km <- function(score_fun, interval = c(1e-4, 500), grid_n = 400L) {
  grid <- exp(seq(log(interval[1]), log(interval[2]), length.out = grid_n))
  vals <- vapply(
    grid,
    function(k) {
      out <- tryCatch(score_fun(k), error = function(e) NA_real_)
      if (length(out) != 1L || !is.finite(out)) {
        return(NA_real_)
      }
      out
    },
    numeric(1)
  )

  ok_pair <- is.finite(vals[-length(vals)]) & is.finite(vals[-1])
  sign_change <- vals[-length(vals)] * vals[-1] <= 0
  idx <- which(ok_pair & sign_change)

  if (length(idx) > 0L) {
    bracket <- idx[1]
    return(
      stats::uniroot(
        score_fun,
        interval = c(grid[bracket], grid[bracket + 1])
      )$root
    )
  }

  obj <- function(log_k) {
    val <- tryCatch(score_fun(exp(log_k)), error = function(e) Inf)
    if (!is.finite(val)) {
      return(Inf)
    }
    val^2
  }

  exp(stats::optimize(obj, interval = log(interval))$minimum)
}

fit_mm_weighted <- function(x, y, h_fun, interval = c(1e-4, 500)) {
  h_vals <- safe_pos(h_fun(x))
  km_hat <- find_km(
    function(k) weighted_score_function(k, x, y, h_vals),
    interval = interval
  )
  vmax_hat <- weighted_vmax_hat(km_hat, x, y, h_vals)

  eps <- y - mm_mean(x, vmax_hat, km_hat)
  gamma_hat <- mean(eps^2 / h_vals)

  grad <- mm_grad(x, vmax_hat, km_hat)
  info <- crossprod(grad, grad / h_vals) / length(x)
  vcov_beta <- gamma_hat * safe_inverse(info) / length(x)

  beta <- c(Vmax = unname(vmax_hat), Km = unname(km_hat))
  list(
    beta = beta,
    vcov_beta = vcov_beta,
    scale = gamma_hat,
    variance = gamma_hat * h_vals
  )
}

fit_mm_nls <- function(x, y) {
  dat <- data.frame(x = x, y = y)
  fit <- stats::nls(
    y ~ (V * x) / (K + x),
    data = dat,
    start = list(V = 1.1 * max(y), K = stats::median(x)),
    control = stats::nls.control(maxiter = 200, warnOnly = TRUE)
  )

  cf <- stats::coef(fit)
  beta <- c(Vmax = unname(cf[["V"]]), Km = unname(cf[["K"]]))
  eps <- y - mm_mean(x, beta["Vmax"], beta["Km"])
  sigma2_hat <- mean(eps^2)
  vcov_beta <- as.matrix(stats::vcov(fit))
  dimnames(vcov_beta) <- list(c("Vmax", "Km"), c("Vmax", "Km"))

  list(
    beta = beta,
    vcov_beta = vcov_beta,
    scale = sigma2_hat,
    variance = rep(sigma2_hat, length(x))
  )
}

## =========================================================
## 5. Bayesian estimator via brms
## =========================================================

require_bayes_packages <- function() {
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("Package 'brms' is required for the Bayesian method.")
  }
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("Package 'posterior' is required for the Bayesian method.")
  }

  invisible(TRUE)
}

get_draw_column <- function(draws_df, candidates, pattern) {
  for (nm in candidates) {
    if (nm %in% names(draws_df)) {
      return(nm)
    }
  }

  hit <- grep(pattern, names(draws_df), value = TRUE)
  if (length(hit) > 0L) {
    return(hit[1])
  }

  stop("Could not find posterior draw column matching: ", pattern)
}

fit_bayes_brms_logvar <- function(x, y,
                                  iter = 1000,
                                  warmup = floor(iter / 2),
                                  chains = 2,
                                  cores = chains,
                                  seed = 2026,
                                  backend = NULL,
                                  adapt_delta = 0.95,
                                  max_treedepth = 12,
                                  refresh = 0,
                                  interval_level = 0.95) {
  require_bayes_packages()

  log_x1 <- safe_pos(log(x + 1))
  dat <- data.frame(
    y = y,
    x = x,
    log_sigma_offset = 0.5 * log(log_x1)
  )

  brms_formula <- brms::bf(
    y ~ exp(logVmax) * x / (exp(logKm) + x),
    logVmax ~ 1,
    logKm ~ 1,
    sigma ~ 1 + offset(log_sigma_offset),
    nl = TRUE
  )

  priors <- c(
    brms::set_prior(
      sprintf("normal(%0.8f, 0.5)", log(V_TRUE)),
      class = "b",
      nlpar = "logVmax"
    ),
    brms::set_prior(
      sprintf("normal(%0.8f, 0.5)", log(K_TRUE)),
      class = "b",
      nlpar = "logKm"
    ),
    brms::set_prior("normal(0, 2)", class = "Intercept", dpar = "sigma")
  )

  fit_args <- list(
    formula = brms_formula,
    data = dat,
    family = stats::gaussian(),
    prior = priors,
    iter = iter,
    warmup = warmup,
    chains = chains,
    cores = cores,
    seed = seed,
    refresh = refresh,
    control = list(
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth
    )
  )

  if (!is.null(backend)) {
    fit_args$backend <- backend
  }

  fit <- do.call(brms::brm, fit_args)
  draws <- posterior::as_draws_df(fit)

  v_col <- get_draw_column(
    draws,
    candidates = c("b_logVmax_Intercept"),
    pattern = "^b_logVmax_Intercept$"
  )
  k_col <- get_draw_column(
    draws,
    candidates = c("b_logKm_Intercept"),
    pattern = "^b_logKm_Intercept$"
  )
  sigma_col <- get_draw_column(
    draws,
    candidates = c("b_sigma_Intercept", "sigma_Intercept"),
    pattern = "sigma.*Intercept"
  )

  vmax_draw <- exp(draws[[v_col]])
  km_draw <- exp(draws[[k_col]])
  gamma_draw <- exp(2 * draws[[sigma_col]])

  v_int <- get_posterior_interval(vmax_draw, level = interval_level)
  k_int <- get_posterior_interval(km_draw, level = interval_level)

  beta <- c(
    Vmax = mean(vmax_draw, na.rm = TRUE),
    Km = mean(km_draw, na.rm = TRUE)
  )
  vcov_beta <- stats::cov(cbind(Vmax = vmax_draw, Km = km_draw), use = "complete.obs")
  gamma_hat <- mean(gamma_draw, na.rm = TRUE)

  list(
    beta = beta,
    vcov_beta = vcov_beta,
    scale = gamma_hat,
    variance = gamma_hat * log_x1,
    lower = c(Vmax = unname(v_int["lower"]), Km = unname(k_int["lower"])),
    upper = c(Vmax = unname(v_int["upper"]), Km = unname(k_int["upper"])),
    brmsfit = fit
  )
}

## =========================================================
## 6. Model library and result collection
## =========================================================

make_model_library <- function(use_bayes = TRUE,
                               bayes_iter = 1000,
                               bayes_warmup = floor(bayes_iter / 2),
                               bayes_chains = 2,
                               bayes_cores = bayes_chains,
                               bayes_backend = NULL,
                               bayes_adapt_delta = 0.95,
                               bayes_max_treedepth = 12,
                               interval_level = 0.95) {
  models <- list(
    NLS = list(
      label = "NLS",
      is_bayes = FALSE,
      fit = function(x, y, seed) fit_mm_nls(x, y)
    ),
    Semi_log = list(
      label = "Semi: log(x+1)",
      is_bayes = FALSE,
      fit = function(x, y, seed) fit_mm_weighted(x, y, h_log)
    ),
    Semi_x_half = list(
      label = "Semi: x^(1/2)",
      is_bayes = FALSE,
      fit = function(x, y, seed) fit_mm_weighted(x, y, h_pow(1 / 2))
    ),
    Semi_x_third = list(
      label = "Semi: x^(1/3)",
      is_bayes = FALSE,
      fit = function(x, y, seed) fit_mm_weighted(x, y, h_pow(1 / 3))
    )
  )

  if (isTRUE(use_bayes)) {
    models$Bayes_logvar <- list(
      label = "Bayesian: Var = gamma log(x+1)",
      is_bayes = TRUE,
      fit = function(x, y, seed) {
        fit_bayes_brms_logvar(
          x = x,
          y = y,
          iter = bayes_iter,
          warmup = bayes_warmup,
          chains = bayes_chains,
          cores = bayes_cores,
          seed = seed,
          backend = bayes_backend,
          adapt_delta = bayes_adapt_delta,
          max_treedepth = bayes_max_treedepth,
          refresh = 0,
          interval_level = interval_level
        )
      }
    )
  }

  models
}

collect_one_fit <- function(fit_obj, true_v, true_beta,
                            interval_level = 0.95) {
  beta_hat <- fit_obj$beta[c("Vmax", "Km")]
  vcov_beta <- as.matrix(fit_obj$vcov_beta)
  se_beta <- sqrt(pmax(diag(vcov_beta)[c("Vmax", "Km")], 0))

  if (!is.null(fit_obj$lower) && !is.null(fit_obj$upper)) {
    lower <- fit_obj$lower[c("Vmax", "Km")]
    upper <- fit_obj$upper[c("Vmax", "Km")]
  } else {
    interval <- get_wald_interval(beta_hat, vcov_beta, level = interval_level)
    lower <- interval$lower[c("Vmax", "Km")]
    upper <- interval$upper[c("Vmax", "Km")]
  }

  c(
    Vmax_hat = unname(beta_hat["Vmax"]),
    Km_hat = unname(beta_hat["Km"]),
    scale_hat = as.numeric(fit_obj$scale)[1],
    se_Vmax = unname(se_beta["Vmax"]),
    se_Km = unname(se_beta["Km"]),
    lower_Vmax = unname(lower["Vmax"]),
    upper_Vmax = unname(upper["Vmax"]),
    lower_Km = unname(lower["Km"]),
    upper_Km = unname(upper["Km"]),
    cover_Vmax = as.numeric(lower["Vmax"] <= true_beta["Vmax"] &&
      true_beta["Vmax"] <= upper["Vmax"]),
    cover_Km = as.numeric(lower["Km"] <= true_beta["Km"] &&
      true_beta["Km"] <= upper["Km"]),
    var_mse = mean((fit_obj$variance - true_v)^2)
  )
}

is_valid_fit_metrics <- function(metrics) {
  needed <- c("Vmax_hat", "Km_hat", "se_Vmax", "se_Km", "var_mse")
  all(is.finite(as.numeric(metrics[needed])))
}

make_result_row <- function(replicate,
                            dgp,
                            model,
                            model_label,
                            error_case,
                            success,
                            runtime = NA_real_,
                            Vmax_hat = NA_real_,
                            Km_hat = NA_real_,
                            scale_hat = NA_real_,
                            se_Vmax = NA_real_,
                            se_Km = NA_real_,
                            lower_Vmax = NA_real_,
                            upper_Vmax = NA_real_,
                            lower_Km = NA_real_,
                            upper_Km = NA_real_,
                            cover_Vmax = NA_real_,
                            cover_Km = NA_real_,
                            var_mse = NA_real_) {
  data.frame(
    replicate = replicate,
    dgp = dgp,
    model = model,
    model_label = model_label,
    error_case = error_case,
    success = success,
    runtime = runtime,
    Vmax_hat = Vmax_hat,
    Km_hat = Km_hat,
    scale_hat = scale_hat,
    se_Vmax = se_Vmax,
    se_Km = se_Km,
    lower_Vmax = lower_Vmax,
    upper_Vmax = upper_Vmax,
    lower_Km = lower_Km,
    upper_Km = upper_Km,
    cover_Vmax = cover_Vmax,
    cover_Km = cover_Km,
    var_mse = var_mse,
    stringsAsFactors = FALSE
  )
}

## =========================================================
## 7. Summaries
## =========================================================

summarize_raw_results <- function(raw_table,
                                  true_beta,
                                  interval_level = 0.95) {
  split_key <- interaction(raw_table$dgp, raw_table$model, drop = TRUE)

  out <- do.call(
    rbind,
    lapply(split(raw_table, split_key), function(dat) {
      ok <- dat$success == 1
      fit_dat <- dat[ok, , drop = FALSE]

      if (nrow(fit_dat) == 0L) {
        return(data.frame(
          DGP = dat$dgp[1],
          Model = dat$model[1],
          Model_Label = dat$model_label[1],
          Error = dat$error_case[1],
          Success_Rate = 0,
          N = 0L,
          Vmax_Bias = NA_real_,
          Vmax_Abs_Bias = NA_real_,
          Vmax_RMSE = NA_real_,
          Km_Bias = NA_real_,
          Km_Abs_Bias = NA_real_,
          Km_RMSE = NA_real_,
          Vmax_Mean_CI_Length = NA_real_,
          Km_Mean_CI_Length = NA_real_,
          Vmax_Interval_Score = NA_real_,
          Km_Interval_Score = NA_real_,
          Vmax_SE_Calibration_Ratio = NA_real_,
          Km_SE_Calibration_Ratio = NA_real_,
          Var_MSE = NA_real_,
          Vmax_CP = NA_real_,
          Km_CP = NA_real_,
          Mean_Runtime = NA_real_
        ))
      }

      vmax_error <- fit_dat$Vmax_hat - true_beta["Vmax"]
      km_error <- fit_dat$Km_hat - true_beta["Km"]
      vmax_ci_length <- fit_dat$upper_Vmax - fit_dat$lower_Vmax
      km_ci_length <- fit_dat$upper_Km - fit_dat$lower_Km
      vmax_interval_score <- interval_score(
        fit_dat$lower_Vmax,
        fit_dat$upper_Vmax,
        truth = true_beta["Vmax"],
        level = interval_level
      )
      km_interval_score <- interval_score(
        fit_dat$lower_Km,
        fit_dat$upper_Km,
        truth = true_beta["Km"],
        level = interval_level
      )

      data.frame(
        DGP = fit_dat$dgp[1],
        Model = fit_dat$model[1],
        Model_Label = fit_dat$model_label[1],
        Error = fit_dat$error_case[1],
        Success_Rate = mean(dat$success),
        N = nrow(fit_dat),
        Vmax_Bias = mean(vmax_error, na.rm = TRUE),
        Vmax_Abs_Bias = mean(abs(vmax_error), na.rm = TRUE),
        Vmax_RMSE = sqrt(mean(vmax_error^2, na.rm = TRUE)),
        Km_Bias = mean(km_error, na.rm = TRUE),
        Km_Abs_Bias = mean(abs(km_error), na.rm = TRUE),
        Km_RMSE = sqrt(mean(km_error^2, na.rm = TRUE)),
        Vmax_Mean_CI_Length = mean(vmax_ci_length, na.rm = TRUE),
        Km_Mean_CI_Length = mean(km_ci_length, na.rm = TRUE),
        Vmax_Interval_Score = mean(vmax_interval_score, na.rm = TRUE),
        Km_Interval_Score = mean(km_interval_score, na.rm = TRUE),
        Vmax_SE_Calibration_Ratio = se_calibration_ratio(
          fit_dat$Vmax_hat,
          fit_dat$se_Vmax
        ),
        Km_SE_Calibration_Ratio = se_calibration_ratio(
          fit_dat$Km_hat,
          fit_dat$se_Km
        ),
        Var_MSE = mean(fit_dat$var_mse, na.rm = TRUE),
        Vmax_CP = mean(fit_dat$cover_Vmax, na.rm = TRUE),
        Km_CP = mean(fit_dat$cover_Km, na.rm = TRUE),
        Mean_Runtime = mean(fit_dat$runtime, na.rm = TRUE)
      )
    })
  )

  row.names(out) <- NULL
  out
}

make_dgp_tables <- function(summary_table,
                            bias_type = c("signed", "absolute"),
                            digits = 5) {
  bias_type <- match.arg(bias_type)
  bias_cols <- if (identical(bias_type, "absolute")) {
    c("Vmax_Abs_Bias", "Km_Abs_Bias")
  } else {
    c("Vmax_Bias", "Km_Bias")
  }
  bias_label <- if (identical(bias_type, "absolute")) "Abs_Bias" else "Bias"

  dgp_split <- split(summary_table, summary_table$DGP)
  dgp_order <- names(make_dgp_library())
  dgp_order <- intersect(dgp_order, names(dgp_split))

  out <- setNames(vector("list", length(dgp_order)), dgp_order)

  for (dgp_name in dgp_order) {
    x <- dgp_split[[dgp_name]]
    x <- x[order(match(x$Model, MANUSCRIPT_MODEL_NAMES)), , drop = FALSE]

    point_table <- data.frame(
      Model = x$Model,
      Model_Label = x$Model_Label,
      Bias_Vmax = x[[bias_cols[1]]],
      Vmax_RMSE = x$Vmax_RMSE,
      Bias_Km = x[[bias_cols[2]]],
      Km_RMSE = x$Km_RMSE,
      stringsAsFactors = FALSE
    )
    names(point_table)[names(point_table) == "Bias_Vmax"] <- paste0(bias_label, "_Vmax")
    names(point_table)[names(point_table) == "Bias_Km"] <- paste0(bias_label, "_Km")

    interval_table <- data.frame(
      Model = x$Model,
      Model_Label = x$Model_Label,
      Coverage_Vmax = x$Vmax_CP,
      Mean_CI_Length_Vmax = x$Vmax_Mean_CI_Length,
      Interval_Score_Vmax = x$Vmax_Interval_Score,
      SE_Calibration_Ratio_Vmax = x$Vmax_SE_Calibration_Ratio,
      Coverage_Km = x$Km_CP,
      Mean_CI_Length_Km = x$Km_Mean_CI_Length,
      Interval_Score_Km = x$Km_Interval_Score,
      SE_Calibration_Ratio_Km = x$Km_SE_Calibration_Ratio,
      Var_MSE = x$Var_MSE,
      Mean_Runtime = x$Mean_Runtime,
      stringsAsFactors = FALSE
    )

    out[[dgp_name]] <- list(
      point = round_numeric_df(point_table, digits = digits),
      interval_variance = round_numeric_df(interval_table, digits = digits),
      bias_type = bias_type
    )
  }

  out
}

make_error_group_tables <- function(summary_table,
                                    bias_type = c("signed", "absolute"),
                                    digits = 5) {
  bias_type <- match.arg(bias_type)
  bias_col_vmax <- if (identical(bias_type, "absolute")) "Vmax_Abs_Bias" else "Vmax_Bias"
  bias_col_km <- if (identical(bias_type, "absolute")) "Km_Abs_Bias" else "Km_Bias"
  bias_label <- if (identical(bias_type, "absolute")) "Abs Bias" else "Bias"

  dat <- summary_table
  dat$Scenario <- vapply(dat$DGP, scenario_from_dgp, character(1))
  dat$Error_Group <- vapply(dat$DGP, error_from_dgp, character(1))
  dat <- dat[is.finite(match(dat$Model, MANUSCRIPT_MODEL_NAMES)), , drop = FALSE]

  build_point <- function(error_case) {
    x <- dat[dat$Error_Group == error_case, , drop = FALSE]
    x <- sort_manuscript_rows(x)
    out <- data.frame(
      Scenario = x$Scenario,
      Method = x$Model_Label,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    out[[sprintf("%s(Vmax)", bias_label)]] <- x[[bias_col_vmax]]
    out[["RMSE(Vmax)"]] <- x$Vmax_RMSE
    out[[sprintf("%s(Km)", bias_label)]] <- x[[bias_col_km]]
    out[["RMSE(Km)"]] <- x$Km_RMSE
    round_numeric_df(out, digits = digits)
  }

  build_diag <- function(error_case) {
    x <- dat[dat$Error_Group == error_case, , drop = FALSE]
    x <- sort_manuscript_rows(x)
    out <- data.frame(
      Scenario = x$Scenario,
      Method = x$Model_Label,
      "Coverage(Vmax)" = x$Vmax_CP,
      "Mean CI Length(Vmax)" = x$Vmax_Mean_CI_Length,
      "Interval Score(Vmax)" = x$Vmax_Interval_Score,
      "SE Calibration Ratio(Vmax)" = x$Vmax_SE_Calibration_Ratio,
      "Coverage(Km)" = x$Km_CP,
      "Mean CI Length(Km)" = x$Km_Mean_CI_Length,
      "Interval Score(Km)" = x$Km_Interval_Score,
      "SE Calibration Ratio(Km)" = x$Km_SE_Calibration_Ratio,
      "Var_MSE" = x$Var_MSE,
      "Mean Runtime" = x$Mean_Runtime,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    round_numeric_df(out, digits = digits)
  }

  list(
    point_normal = build_point("Normal"),
    point_skewed = build_point("Skewed"),
    diagnostics_normal = build_diag("Normal"),
    diagnostics_skewed = build_diag("Skewed"),
    bias_type = bias_type
  )
}

## =========================================================
## 8. Output writers
## =========================================================

write_simulation_outputs <- function(sim_obj,
                                     summary_file = paste0(MANUSCRIPT_FILE_PREFIX, "_summary.csv"),
                                     raw_file = paste0(MANUSCRIPT_FILE_PREFIX, "_raw.csv"),
                                     dgp_file = paste0(MANUSCRIPT_FILE_PREFIX, "_dgp_reference.csv"),
                                     model_file = paste0(MANUSCRIPT_FILE_PREFIX, "_model_reference.csv"),
                                     dir = getwd()) {
  utils::write.csv(sim_obj$summary, file.path(dir, summary_file), row.names = FALSE)
  utils::write.csv(sim_obj$raw, file.path(dir, raw_file), row.names = FALSE)
  utils::write.csv(sim_obj$dgp_reference, file.path(dir, dgp_file), row.names = FALSE)
  utils::write.csv(sim_obj$model_reference, file.path(dir, model_file), row.names = FALSE)

  invisible(sim_obj)
}

write_dgp_table_outputs <- function(dgp_tables,
                                    save_dir = getwd(),
                                    dir_name = "table_final",
                                    file_prefix = MANUSCRIPT_FILE_PREFIX) {
  out_dir <- file.path(save_dir, dir_name)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  bias_suffix <- if (length(dgp_tables) > 0L &&
    identical(dgp_tables[[1]]$bias_type, "absolute")) {
    "abs_bias"
  } else {
    "signed_bias"
  }

  for (dgp_name in names(dgp_tables)) {
    tables <- dgp_tables[[dgp_name]]
    safe_name <- gsub("[^A-Za-z0-9]+", "_", dgp_name)

    point_file <- paste0(
      file_prefix, "_", safe_name,
      "_table1_point_estimation_", bias_suffix, ".csv"
    )
    interval_file <- paste0(
      file_prefix, "_", safe_name,
      "_table2_interval_variance_diagnostics.csv"
    )

    utils::write.csv(
      tables$point,
      file.path(out_dir, point_file),
      row.names = FALSE
    )
    utils::write.csv(
      tables$interval_variance,
      file.path(out_dir, interval_file),
      row.names = FALSE
    )
  }

  invisible(dgp_tables)
}

write_error_group_tables <- function(group_tables,
                                     save_dir = getwd(),
                                     dir_name = file.path("table_final", "notion_csv_4tables"),
                                     file_prefix = MANUSCRIPT_FILE_PREFIX) {
  out_dir <- file.path(save_dir, dir_name)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  if (identical(group_tables$bias_type, "absolute")) {
    point_normal_file <- paste0(file_prefix, "_table1_point_normal_abs_bias.csv")
    point_skewed_file <- paste0(file_prefix, "_table1_point_skewed_abs_bias.csv")
  } else {
    point_normal_file <- paste0(file_prefix, "_table1_point_normal.csv")
    point_skewed_file <- paste0(file_prefix, "_table1_point_skewed.csv")
  }

  utils::write.csv(
    group_tables$point_normal,
    file.path(out_dir, point_normal_file),
    row.names = FALSE
  )
  utils::write.csv(
    group_tables$point_skewed,
    file.path(out_dir, point_skewed_file),
    row.names = FALSE
  )
  utils::write.csv(
    group_tables$diagnostics_normal,
    file.path(out_dir, paste0(file_prefix, "_table2_diagnostics_normal.csv")),
    row.names = FALSE
  )
  utils::write.csv(
    group_tables$diagnostics_skewed,
    file.path(out_dir, paste0(file_prefix, "_table2_diagnostics_skewed.csv")),
    row.names = FALSE
  )

  invisible(group_tables)
}

write_all_manuscript_outputs <- function(summary_table,
                                         save_dir = getwd(),
                                         file_prefix = MANUSCRIPT_FILE_PREFIX,
                                         bias_type = c("signed", "absolute"),
                                         digits = 5) {
  bias_type <- match.arg(bias_type)

  dgp_tables <- make_dgp_tables(
    summary_table = summary_table,
    bias_type = bias_type,
    digits = digits
  )
  grouped_tables <- make_error_group_tables(
    summary_table = summary_table,
    bias_type = bias_type,
    digits = digits
  )

  write_dgp_table_outputs(
    dgp_tables,
    save_dir = save_dir,
    file_prefix = file_prefix
  )
  write_error_group_tables(
    grouped_tables,
    save_dir = save_dir,
    file_prefix = file_prefix
  )

  invisible(list(
    dgp_tables = dgp_tables,
    grouped_tables = grouped_tables
  ))
}

rebuild_manuscript_tables_from_csv <- function(summary_csv,
                                               save_dir = dirname(summary_csv),
                                               file_prefix = MANUSCRIPT_FILE_PREFIX,
                                               bias_type = c("signed", "absolute"),
                                               digits = 5) {
  bias_type <- match.arg(bias_type)
  summary_table <- utils::read.csv(summary_csv, stringsAsFactors = FALSE)

  write_all_manuscript_outputs(
    summary_table = summary_table,
    save_dir = save_dir,
    file_prefix = file_prefix,
    bias_type = bias_type,
    digits = digits
  )
}

## =========================================================
## 9. Main simulation runner
## =========================================================

run_michaelis_menten_simulation <- function(n_sim = 1000L,
                                            n_obs = 100L,
                                            x = NULL,
                                            seed = 2026,
                                            vmax_true = V_TRUE,
                                            km_true = K_TRUE,
                                            error_cases = MANUSCRIPT_ERROR_CASES,
                                            dgp_names = NULL,
                                            use_bayes = TRUE,
                                            bayes_every = 1L,
                                            bayes_iter = 2000,
                                            bayes_warmup = floor(bayes_iter / 2),
                                            bayes_chains = 4,
                                            bayes_cores = bayes_chains,
                                            bayes_backend = "cmdstanr",
                                            bayes_adapt_delta = 0.95,
                                            bayes_max_treedepth = 12,
                                            interval_level = 0.95,
                                            save_dir = getwd(),
                                            file_prefix = MANUSCRIPT_FILE_PREFIX,
                                            verbose = TRUE,
                                            show_errors = FALSE,
                                            digits = 5,
                                            bias_type = c("signed", "absolute")) {
  bias_type <- match.arg(bias_type)

  if (is.null(x)) {
    x <- seq(1, 100, length.out = n_obs)
  } else {
    n_obs <- length(x)
  }

  if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  }

  error_cases <- match.arg(
    error_cases,
    MANUSCRIPT_ERROR_CASES,
    several.ok = TRUE
  )

  if (isTRUE(use_bayes)) {
    require_bayes_packages()
  }

  set.seed(seed)
  true_beta <- c(Vmax = vmax_true, Km = km_true)
  dgp_library <- make_dgp_library(error_cases = error_cases)

  if (is.null(dgp_names)) {
    dgp_names <- names(dgp_library)
  }
  dgp_names <- match.arg(dgp_names, names(dgp_library), several.ok = TRUE)
  dgp_reference <- dgp_library[dgp_names]

  model_library <- make_model_library(
    use_bayes = use_bayes,
    bayes_iter = bayes_iter,
    bayes_warmup = bayes_warmup,
    bayes_chains = bayes_chains,
    bayes_cores = bayes_cores,
    bayes_backend = bayes_backend,
    bayes_adapt_delta = bayes_adapt_delta,
    bayes_max_treedepth = bayes_max_treedepth,
    interval_level = interval_level
  )

  rows <- list()
  row_id <- 1L
  progress_every <- max(1L, floor(n_sim / 10))

  for (dgp_name in dgp_names) {
    dgp <- dgp_library[[dgp_name]]

    if (isTRUE(verbose)) {
      message("Running ", dgp_name, " (", dgp$label, ")")
    }

    for (sim_id in seq_len(n_sim)) {
      if (isTRUE(verbose) &&
          (sim_id == 1L || sim_id == n_sim || sim_id %% progress_every == 0L)) {
        message("  simulation ", sim_id, "/", n_sim)
      }

      sim_dat <- simulate_response(
        x = x,
        dgp_spec = dgp,
        vmax_true = vmax_true,
        km_true = km_true
      )

      for (model_name in names(model_library)) {
        model <- model_library[[model_name]]

        if (isTRUE(model$is_bayes) && sim_id %% bayes_every != 0L) {
          next
        }

        start_time <- proc.time()[["elapsed"]]
        fit <- tryCatch(
          model$fit(x, sim_dat$y, seed = seed + 10000L * sim_id + row_id),
          error = function(e) {
            if (isTRUE(show_errors)) {
              message(
                "[", dgp_name, "] ",
                model_name,
                " failed at simulation ",
                sim_id,
                ": ",
                conditionMessage(e)
              )
            }
            NULL
          }
        )
        runtime <- proc.time()[["elapsed"]] - start_time

        if (is.null(fit)) {
          rows[[row_id]] <- make_result_row(
            replicate = sim_id,
            dgp = dgp_name,
            model = model_name,
            model_label = model$label,
            error_case = dgp$error_case,
            success = 0,
            runtime = runtime
          )
        } else {
          metrics <- collect_one_fit(
            fit,
            true_v = sim_dat$v_true,
            true_beta = true_beta,
            interval_level = interval_level
          )

          if (!is_valid_fit_metrics(metrics)) {
            rows[[row_id]] <- make_result_row(
              replicate = sim_id,
              dgp = dgp_name,
              model = model_name,
              model_label = model$label,
              error_case = dgp$error_case,
              success = 0,
              runtime = runtime
            )
          } else {
            rows[[row_id]] <- make_result_row(
              replicate = sim_id,
              dgp = dgp_name,
              model = model_name,
              model_label = model$label,
              error_case = dgp$error_case,
              success = 1,
              runtime = runtime,
              Vmax_hat = unname(metrics["Vmax_hat"]),
              Km_hat = unname(metrics["Km_hat"]),
              scale_hat = unname(metrics["scale_hat"]),
              se_Vmax = unname(metrics["se_Vmax"]),
              se_Km = unname(metrics["se_Km"]),
              lower_Vmax = unname(metrics["lower_Vmax"]),
              upper_Vmax = unname(metrics["upper_Vmax"]),
              lower_Km = unname(metrics["lower_Km"]),
              upper_Km = unname(metrics["upper_Km"]),
              cover_Vmax = unname(metrics["cover_Vmax"]),
              cover_Km = unname(metrics["cover_Km"]),
              var_mse = unname(metrics["var_mse"])
            )
          }
        }

        row_id <- row_id + 1L
      }
    }
  }

  raw_table <- do.call(rbind, rows)
  summary_unrounded <- summarize_raw_results(
    raw_table,
    true_beta = true_beta,
    interval_level = interval_level
  )
  summary_table <- round_numeric_df(summary_unrounded, digits = digits)

  out <- list(
    raw = raw_table,
    summary_unrounded = summary_unrounded,
    summary = summary_table,
    dgp_reference = data.frame(
      dgp = names(dgp_reference),
      label = vapply(dgp_reference, `[[`, character(1), "label"),
      error_case = vapply(dgp_reference, `[[`, character(1), "error_case"),
      stringsAsFactors = FALSE
    ),
    model_reference = data.frame(
      model = names(model_library),
      label = vapply(model_library, `[[`, character(1), "label"),
      stringsAsFactors = FALSE
    ),
    truth = true_beta,
    x = x
  )

  write_simulation_outputs(
    out,
    summary_file = paste0(file_prefix, "_summary.csv"),
    raw_file = paste0(file_prefix, "_raw.csv"),
    dgp_file = paste0(file_prefix, "_dgp_reference.csv"),
    model_file = paste0(file_prefix, "_model_reference.csv"),
    dir = save_dir
  )

  out$manuscript_tables <- write_all_manuscript_outputs(
    summary_table = summary_unrounded,
    save_dir = save_dir,
    file_prefix = file_prefix,
    bias_type = bias_type,
    digits = digits
  )

  invisible(out)
}

## =========================================================
## Example usage
## =========================================================

## source("michaelis_menten_simulation.R")
##
## Rebuild manuscript tables from an existing summary CSV:
## rebuild_manuscript_tables_from_csv(
##   "results/chapter3_simulation_n50/simulation_real_manuscript_n50_summary.csv"
## )
##
## Full manuscript run:
## res <- run_michaelis_menten_simulation(
##   n_sim = 1000,
##   n_obs = 50,
##   use_bayes = FALSE,
##   save_dir = "results/chapter3_simulation_n50",
##   file_prefix = "simulation_real_manuscript_n50"
## )
