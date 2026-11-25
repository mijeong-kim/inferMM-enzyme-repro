# Run from BiomJ_rev: Rscript analysis/revision_benchmark.R [replications]
# All inputs are frozen copies; this script never writes outside BiomJ_rev.
args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args)) as.integer(args[1]) else 1000L
stopifnot(is.finite(R), R >= 10)
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
root <- dirname(dirname(normalizePath(sub("^--file=", "", file_arg))))
out <- file.path(root, "analysis", "results")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
loader <- new.env(parent = globalenv())
sys.source(file.path(root, "examples/load_example_package.R"), loader)
loader$load_example_package(root)
pkg <- asNamespace("inferMM")
stopifnot(requireNamespace("nlme", quietly = TRUE))
# Prespecified seven-concentration subset; the full observed SDL panel has eight.
sdl_levels <- c(7.8125, 15.625, 31.25, 62.5, 125, 250, 500)
designs <- list(Grid20 = seq(1, 100, length.out = 20),
                SDL28 = rep(sdl_levels / max(sdl_levels) * 100, each = 4))
vfun <- list(Constant = function(s) rep(6, length(s)),
             MM = function(s) 1 + 9 * s / (20 + s),
             Exp = function(s) 1 + 9 * (1 - exp(-0.05 * s)),
             Hill = function(s) 1 + 9 * s^2 / (20^2 + s^2))
truth <- c(Vmax = 100, Km = 20)
model_names <- c("Constant", "Log", "Sqrt", "Cuberoot", "AIC screen",
                 "Constant HC3", "GNLS power")

fit_record <- function(fit, x, y, v, method) {
  if (is.null(fit)) return(data.frame(method = method, ok = FALSE,
    Vmax = NA_real_, Km = NA_real_, seV = NA_real_, seK = NA_real_,
    coverV = NA_real_, coverK = NA_real_, widthV = NA_real_, widthK = NA_real_,
    var_mse = NA_real_, selected = NA_character_, score_error = NA_real_))
  beta <- fit$coefficients
  se <- sqrt(diag(fit$vcov))
  z <- qnorm(.975)
  h <- fit$variance_fun(x)
  # The normalized residual measures stationarity independently of score scale.
  g <- pkg$mm_grad(x, beta[1], beta[2])
  residual <- y - pkg$mm_mean(x, beta[1], beta[2])
  sc <- colSums(g * (residual / h))
  sc_scale <- colSums(abs(g) * (abs(y) / h))
  score_error <- max(abs(sc) / pmax(sc_scale, 1))
  good <- all(is.finite(c(beta, se))) && all(beta > 0) && all(se > 0)
  data.frame(method = method, ok = good, Vmax = beta[1], Km = beta[2],
    seV = se[1], seK = se[2],
    coverV = abs(beta[1] - truth[1]) <= z * se[1],
    coverK = abs(beta[2] - truth[2]) <= z * se[2],
    widthV = 2 * z * se[1], widthK = 2 * z * se[2],
    var_mse = mean((fit$variance_hat - v)^2),
    selected = fit$variance_model, score_error = score_error, row.names = NULL)
}

run_case <- function(case) {
  di <- case$design; sc <- case$scenario; ei <- case$error
  x <- designs[[di]]; v <- vfun[[sc]](x)
  mu <- pkg$mm_mean(x, truth[1], truth[2])
  set.seed(9078092L + case$id)
  rows <- vector("list", R)
  for (r in seq_len(R)) {
    innovation <- if (ei == "Normal") rnorm(length(x)) else
      (rgamma(length(x), shape = 4) - 4) / 2
    y <- mu + sqrt(v) * innovation
    fits <- lapply(c("constant", "log", "sqrt", "cuberoot"), function(h)
      tryCatch(pkg$fit_mm(x, y, variance = h), error = function(e) NULL))
    records <- lapply(seq_along(fits), function(j)
      fit_record(fits[[j]], x, y, v, model_names[j]))
    criterion <- vapply(fits, function(f) if (is.null(f)) Inf else f$quasi_aic,
                        numeric(1))
    selected <- if (any(is.finite(criterion))) fits[[which.min(criterion)]] else NULL
    records[[5]] <- fit_record(selected, x, y, v, model_names[5])
    hc <- fits[[1]]
    if (!is.null(hc)) {
      hc <- tryCatch({
        J <- pkg$mm_grad(x, hc$coefficients[1], hc$coefficients[2])
        bread <- solve(crossprod(J))
        lev <- rowSums((J %*% bread) * J)
        adjusted <- hc$residuals / pmax(1 - lev, 1e-8)
        hc$vcov <- bread %*% crossprod(J, J * adjusted^2) %*% bread
        hc
      }, error = function(e) NULL)
    }
    records[[6]] <- fit_record(hc, x, y, v, model_names[6])
    gn <- tryCatch(suppressWarnings({
      start <- if (is.null(fits[[1]])) c(100, 20) else fits[[1]]$coefficients
      f <- nlme::gnls(y ~ exp(lV) * x / (exp(lK) + x),
        data = data.frame(x, y),
        start = c(lV = log(unname(start[1])), lK = log(unname(start[2]))),
        weights = nlme::varPower(form = ~ x),
        control = nlme::gnlsControl(maxIter = 100, nlsMaxIter = 50))
      be <- exp(coef(f)); names(be) <- names(truth)
      delta <- diag(be, 2)
      covariance <- delta %*% vcov(f) %*% delta
      exponent <- as.numeric(coef(f$modelStruct$varStruct, unconstrained = FALSE))
      hf <- local({p <- 2 * exponent; function(s) s^p})
      list(coefficients = be, vcov = covariance, variance_fun = hf,
           variance_hat = f$sigma^2 * hf(x), variance_model = "fitted power")
    }), error = function(e) NULL)
    records[[7]] <- fit_record(gn, x, y, v, model_names[7])
    rows[[r]] <- cbind(design = di, scenario = sc, error = ei,
                      replicate = r, do.call(rbind, records))
  }
  raw <- do.call(rbind, rows)
  write.csv(raw, file.path(out, paste0("benchmark_", di, "_", sc, "_", ei, "_raw.csv")),
            row.names = FALSE)
  cat(di, sc, ei, "completed", R, "replications\n")
  raw
}

cases <- expand.grid(design = names(designs), scenario = names(vfun),
                     error = c("Normal", "Skewed"), stringsAsFactors = FALSE)
cases$id <- seq_len(nrow(cases))
case_list <- split(cases, cases$id)
# Separate deterministic streams for each case; scheduling does not change results.
cores <- if (.Platform$OS.type == "windows") 1L else 4L
raw_list <- parallel::mclapply(case_list, run_case, mc.cores = cores,
                              mc.set.seed = FALSE, mc.preschedule = FALSE)
stopifnot(all(vapply(raw_list, is.data.frame, logical(1))))
raw <- do.call(rbind, raw_list)
write.csv(raw, file.path(out, "revision_benchmark_raw.csv"), row.names = FALSE)
groups <- split(raw, interaction(raw$design, raw$scenario, raw$error,
                                raw$method, drop = TRUE))
summ <- do.call(rbind, lapply(groups, function(d) {
  z <- d[d$ok, ]; n <- nrow(z)
  m <- function(a) if (n) mean(a) else NA_real_
  cpv <- m(z$coverV); cpk <- m(z$coverK)
  data.frame(design = d$design[1], scenario = d$scenario[1], error = d$error[1],
    method = d$method[1], successful = n, attempted = nrow(d),
    biasV = m(z$Vmax - 100), rmseV = sqrt(m((z$Vmax - 100)^2)),
    biasK = m(z$Km - 20), rmseK = sqrt(m((z$Km - 20)^2)),
    cpV = cpv, mcseV = sqrt(cpv * (1-cpv) / n),
    cpK = cpk, mcseK = sqrt(cpk * (1-cpk) / n),
    widthV = m(z$widthV), widthK = m(z$widthK), var_mse = m(z$var_mse),
    max_score_error = if(n) max(z$score_error) else NA_real_)
}))
row.names(summ) <- NULL
write.csv(summ, file.path(out, "revision_benchmark_summary.csv"), row.names = FALSE)
selection <- as.data.frame(xtabs(~design + scenario + error + selected,
                                raw[raw$method == "AIC screen" & raw$ok, ]))
write.csv(selection, file.path(out, "selection_counts.csv"), row.names = FALSE)
paired <- merge(raw[raw$method == "Constant", ],
                raw[raw$method %in% c("Log", "Sqrt", "Cuberoot"), ],
                by = c("design", "scenario", "error", "replicate"),
                suffixes = c("_baseline", "_weighted"))
paired <- paired[paired$ok_baseline & paired$ok_weighted, ]
pg <- split(paired, interaction(paired$design, paired$scenario, paired$error,
                                paired$method_weighted, drop = TRUE))
ps <- do.call(rbind, lapply(pg, function(d) {
  dv <- d$coverV_weighted - d$coverV_baseline
  dk <- d$coverK_weighted - d$coverK_baseline
  data.frame(design=d$design[1], scenario=d$scenario[1], error=d$error[1],
    method=d$method_weighted[1], paired_n=nrow(d), difference_cpV=mean(dv),
    mcse_difference_cpV=sd(dv)/sqrt(nrow(d)), difference_cpK=mean(dk),
    mcse_difference_cpK=sd(dk)/sqrt(nrow(d)))
}))
write.csv(ps, file.path(out, "paired_coverage_differences.csv"), row.names=FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
write.csv(data.frame(S = designs$SDL28), file.path(out, "SDL28_design.csv"), row.names=FALSE)

# Figure requested by Reviewer 1; all three generating functions on common axes.
png(file.path(root, "figures", "variance_functions_rev1.png"),
    width = 1800, height = 1050, res = 200)
par(mar = c(4.4, 4.5, 1.5, .7), las = 1, cex = 1.05)
s <- seq(1, 100, length.out = 400)
matplot(s, cbind(vfun$MM(s), vfun$Exp(s), vfun$Hill(s)), type="l", lwd=2.5,
        lty=c(1,2,3), col=c("#155c8a", "#b14d19", "#30784c"),
        xlab="Substrate concentration", ylab="Conditional variance", ylim=c(1,10))
legend("bottomright", c("MM: Michaelis-Menten saturation", "Exp: exponential saturation",
                       "Hill: Hill-type saturation"), lty=c(1,2,3), lwd=2.5,
       col=c("#155c8a", "#b14d19", "#30784c"), bty="n", cex=.85)
dev.off()
cat("All revision benchmarks completed.\n")
