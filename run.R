# Launch each analysis in an isolated subprocess. Never write to expected/.
args <- commandArgs(TRUE)
mode <- if (length(args)) args[1] else "help"
file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
home <- dirname(normalizePath(file))
Sys.setenv(INFERMM_REPRO_ROOT = home)
rscript <- file.path(R.home("bin"), "Rscript")
package_loader <- new.env(parent = globalenv())
sys.source(file.path(home, "examples/load_example_package.R"), package_loader)
for (part in c("outputs", "results", "sections", "figures"))
  dir.create(file.path(home, part), showWarnings = FALSE)
writeLines(capture.output(sessionInfo()), file.path(home, "outputs", paste0(mode, "_session.txt")))
call <- function(file, args = character()) {
  status <- system2(rscript, c(shQuote(file), shQuote(args)))
  if (status != 0L) stop("Analysis failed: ", basename(file), "; exit ", status)
}
dependencies <- function(packages) {
  target <- read.csv(file.path(home, "dependency_versions.csv"))
  if (getRversion() != "4.6.1") warning("Reference results used R 4.6.1.")
  for (package in packages) {
    if (!requireNamespace(package, quietly = TRUE)) stop("Install dependency: ", package)
    expected <- target$version[match(package, target$package)]
    if (!is.na(expected) && as.character(packageVersion(package)) != expected)
      stop("Reference version required for ", package, ": ", expected,
           "; installed ", as.character(packageVersion(package)))
  }
}
if (mode == "all") {
  package_loader$check_bundled_package(home, "0.1.2", "library012")
  for (part in c("single50", "single100", "cluster", "real", "benchmark", "matched", "higher_moments", "bootstrap20", "bootstrap30", "projected"))
    call(file.path(home, "run.R"), part)
} else if (mode %in% c("single50", "single100")) {
  env <- new.env(parent = globalenv())
  sys.source(file.path(home, "R/michaelis_menten_simulation.R"), env)
  n <- if (mode == "single50") 50L else 100L
  env$run_michaelis_menten_simulation(n_sim = 1000L, n_obs = n,
    seed = 2026L, use_bayes = FALSE, save_dir = file.path(home, "outputs", mode),
    file_prefix = paste0("rerun_n", n), verbose = TRUE)
} else if (mode == "cluster") {
  dependencies("nlme")
  env <- new.env(parent = globalenv())
  sys.source(file.path(home, "R/realistic_mm_empirical_design_simulation.R"), env)
  env$PROJECT_ROOT <- home
  env$SDL_TEMPLATE_PATH <- file.path(home, "data/sdl_demo.csv")
  env$ALVES_TEMPLATE_PATH <- file.path(home, "data/raw/Alves_data.csv")
  templates <- env$build_clustered_template_suite()
  out <- file.path(home, "outputs/cluster")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  methods <- c("pooled_constant", "pooled_sqrt", "clustered_sqrt",
               "nlme_homoscedastic", "nlme_varConstPower")
  results <- lapply(seq_along(templates), function(i) {
    log <- character()
    result <- withCallingHandlers(env$run_alves_empirical_clustered_simulation(
      template = templates[[i]], n_sim = 1000L,
      seed = 20260509L + (i - 1L) * 1000L, method_names = methods),
      warning = function(w) { log <<- c(log, conditionMessage(w)); invokeRestart("muffleWarning") })
    result$raw$template_name <- names(templates)[i]
    write.csv(as.data.frame(table(log)), file.path(out, paste0("cluster_", i, "_warnings.csv")), row.names = FALSE)
    result$raw
  })
  rawpath <- file.path(out, "cluster_raw.csv")
  write.csv(do.call(rbind, results), rawpath, row.names = FALSE)
  call(file.path(home, "R/summarize_cluster.R"), c(rawpath, out))
} else if (mode == "real") {
  dependencies(c("dplyr", "tidyr", "stringr", "purrr", "ggplot2", "knitr"))
  work <- file.path(home, "work/real")
  for (part in c("data/raw", "scripts", "results")) dir.create(file.path(work, part), recursive = TRUE, showWarnings = FALSE)
  stopifnot(all(file.copy(list.files(file.path(home, "data/raw"), full.names = TRUE), file.path(work, "data/raw"), overwrite = TRUE)))
  code <- file.path(work, "scripts/real_data_extracted.R")
  knitr::purl(file.path(home, "scripts/real_data.Rmd"), output = code, quiet = TRUE)
  setwd(work)
  sys.source(code, new.env(parent = globalenv()))
} else if (mode %in% c("benchmark", "matched")) {
  dependencies(if (mode == "matched") c("nlme", "drc") else "nlme")
  script <- if (mode == "benchmark") "revision_benchmark.R" else "matched_comparisons.R"
  call(file.path(home, "analysis", script), if (mode == "benchmark") "1000" else c("1000", "final"))
  render <- if (mode == "benchmark") "render_benchmark_tables.R" else "render_matched_tables.R"
  call(file.path(home, "analysis", render))
} else if (mode == "higher_moments") {
  env <- new.env(parent=globalenv())
  sys.source(file.path(home,"examples/load_example_package.R"),env)
  lib <- env$load_example_package(home)
  call(file.path(home,"analysis/higher_moments/benchmark.R"),
       c(file.path(home,"outputs/higher_moments"),"1000",lib,"3"))
  call(file.path(home,"analysis/render_higher_moment_tables.R"))
} else if (mode %in% c("bootstrap20", "bootstrap30", "bootstrap_smoke")) {
  package_loader$check_bundled_package(home, "0.1.2", "library012")
  bootstrap_mode <- switch(mode, bootstrap20 = "full", bootstrap30 = "full30", bootstrap_smoke = "smoke")
  call(file.path(home, "run_bootstrap.R"), bootstrap_mode)
} else if (mode == "projected") {
  package_loader$check_bundled_package(home, "0.1.2", "library012")
  call(file.path(home, "examples/validate_ma_genton.R"))
} else {
  cat("Rscript run.R all\nModules: single50 single100 cluster real benchmark matched higher_moments bootstrap20 bootstrap30 projected\n")
  if (mode != "help") stop("Unknown module: ", mode)
}
