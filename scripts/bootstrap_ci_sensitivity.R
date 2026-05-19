find_repo_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)

  repeat {
    if (dir.exists(file.path(current, "scripts")) &&
        dir.exists(file.path(current, "data")) &&
        dir.exists(file.path(current, "results"))) {
      return(current)
    }

    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not locate the inferMM-cils-repro project root.")
    }
    current <- parent
  }
}

repo_path <- function(...) {
  file.path(find_repo_root(), ...)
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

se_calibration_ratio <- function(estimates, standard_errors) {
  ok <- is.finite(estimates) & is.finite(standard_errors) & standard_errors > 0
  if (sum(ok) < 2L) {
    return(NA_real_)
  }

  stats::sd(estimates[ok]) / mean(standard_errors[ok])
}

collect_ci_metrics <- function(fit,
                               truth,
                               method = c("wald", "bootstrap"),
                               level = 0.95,
                               B = 99L,
                               seed = NULL) {
  method <- match.arg(method)

  timing <- system.time({
    ci <- confint(
      fit,
      level = level,
      method = method,
      B = B,
      seed = seed
    )
  })[["elapsed"]]

  if (identical(method, "wald")) {
    se <- sqrt(pmax(diag(vcov(fit)), 0))
  } else {
    draws <- attr(ci, "bootstrap_draws")
    se <- apply(draws, 2, stats::sd)
  }

  c(
    lower_Vmax = ci["Vmax", "lower"],
    upper_Vmax = ci["Vmax", "upper"],
    lower_Km = ci["Km", "lower"],
    upper_Km = ci["Km", "upper"],
    se_Vmax = se["Vmax"],
    se_Km = se["Km"],
    cover_Vmax = as.numeric(ci["Vmax", "lower"] <= truth["Vmax"] &&
      truth["Vmax"] <= ci["Vmax", "upper"]),
    cover_Km = as.numeric(ci["Km", "lower"] <= truth["Km"] &&
      truth["Km"] <= ci["Km", "upper"]),
    runtime = timing
  )
}

run_bootstrap_ci_sensitivity <- function(n = 100L,
                                         n_sim = 1000L,
                                         B = 99L,
                                         level = 0.95,
                                         seed = 20260504L,
                                         save_csv = NULL) {
  if (!requireNamespace("inferMM", quietly = TRUE)) {
    stop("Package 'inferMM' must be installed before running this script.")
  }
  library(inferMM)

  truths <- c(Vmax = 100, Km = 20)
  variance_shapes <- c("mm", "exp", "hill")
  error_shapes <- c("normal", "skewed")
  rows <- vector("list", length = n_sim * length(variance_shapes) * length(error_shapes) * 2L)
  row_id <- 1L

  set.seed(seed)

  for (variance_shape in variance_shapes) {
    for (error_shape in error_shapes) {
      for (replicate in seq_len(n_sim)) {
        sim_seed <- seed + 10000L * match(variance_shape, variance_shapes) +
          1000L * match(error_shape, error_shapes) + replicate

        dat <- simulate_mm_data(
          n = n,
          seed = sim_seed,
          variance_shape = variance_shape,
          error = error_shape
        )

        fit <- fit_mm(dat$x, dat$y, variance = "sqrt")
        beta_hat <- coef(fit)

        wald_metrics <- collect_ci_metrics(
          fit = fit,
          truth = truths,
          method = "wald",
          level = level
        )
        boot_metrics <- collect_ci_metrics(
          fit = fit,
          truth = truths,
          method = "bootstrap",
          level = level,
          B = B,
          seed = sim_seed + 500000L
        )

        rows[[row_id]] <- data.frame(
          Scenario = toupper(variance_shape),
          Error = tools::toTitleCase(error_shape),
          Method = "Semi: S^(1/2), Wald CI",
          replicate = replicate,
          Vmax_hat = beta_hat["Vmax"],
          Km_hat = beta_hat["Km"],
          lower_Vmax = wald_metrics["lower_Vmax"],
          upper_Vmax = wald_metrics["upper_Vmax"],
          lower_Km = wald_metrics["lower_Km"],
          upper_Km = wald_metrics["upper_Km"],
          se_Vmax = wald_metrics["se_Vmax"],
          se_Km = wald_metrics["se_Km"],
          cover_Vmax = wald_metrics["cover_Vmax"],
          cover_Km = wald_metrics["cover_Km"],
          ci_runtime = wald_metrics["runtime"],
          stringsAsFactors = FALSE
        )
        row_id <- row_id + 1L

        rows[[row_id]] <- data.frame(
          Scenario = toupper(variance_shape),
          Error = tools::toTitleCase(error_shape),
          Method = "Semi: S^(1/2), Wild bootstrap CI",
          replicate = replicate,
          Vmax_hat = beta_hat["Vmax"],
          Km_hat = beta_hat["Km"],
          lower_Vmax = boot_metrics["lower_Vmax"],
          upper_Vmax = boot_metrics["upper_Vmax"],
          lower_Km = boot_metrics["lower_Km"],
          upper_Km = boot_metrics["upper_Km"],
          se_Vmax = boot_metrics["se_Vmax"],
          se_Km = boot_metrics["se_Km"],
          cover_Vmax = boot_metrics["cover_Vmax"],
          cover_Km = boot_metrics["cover_Km"],
          ci_runtime = boot_metrics["runtime"],
          stringsAsFactors = FALSE
        )
        row_id <- row_id + 1L
      }
    }
  }

  raw <- do.call(rbind, rows)
  summary <- do.call(
    rbind,
    lapply(split(raw, interaction(raw$Scenario, raw$Error, raw$Method, drop = TRUE)), function(dat) {
      data.frame(
        Scenario = dat$Scenario[1],
        Error = dat$Error[1],
        Method = dat$Method[1],
        Coverage_Vmax = mean(dat$cover_Vmax),
        Mean_CI_Length_Vmax = mean(dat$upper_Vmax - dat$lower_Vmax),
        Interval_Score_Vmax = mean(interval_score(dat$lower_Vmax, dat$upper_Vmax, truths["Vmax"], level = level)),
        SECR_Vmax = se_calibration_ratio(dat$Vmax_hat, dat$se_Vmax),
        Coverage_Km = mean(dat$cover_Km),
        Mean_CI_Length_Km = mean(dat$upper_Km - dat$lower_Km),
        Interval_Score_Km = mean(interval_score(dat$lower_Km, dat$upper_Km, truths["Km"], level = level)),
        SECR_Km = se_calibration_ratio(dat$Km_hat, dat$se_Km),
        Mean_CI_Runtime = mean(dat$ci_runtime),
        stringsAsFactors = FALSE
      )
    })
  )

  row.names(summary) <- NULL

  if (!is.null(save_csv)) {
    utils::write.csv(summary, save_csv, row.names = FALSE)
  }

  list(raw = raw, summary = summary)
}

if (sys.nframe() == 0L) {
  out <- run_bootstrap_ci_sensitivity(
    n_sim = 1000L,
    B = 49L,
    save_csv = repo_path(
      "results",
      "chapter3_simulation",
      "table",
      "simulation_real_manuscript_bootstrap_ci_sensitivity.csv"
    )
  )
  print(out$summary)
}
