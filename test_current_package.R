# Optional API tests for the supplied current archive, separate from table reruns.
file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
home <- dirname(normalizePath(file))
stopifnot(requireNamespace("testthat", quietly = TRUE))
out <- file.path(home, "outputs/package_tests")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
untar(file.path(home, "packages/inferMM_0.1.2.tar.gz"), exdir = out)
result <- testthat::test_local(file.path(out, "inferMM"), reporter = "summary", stop_on_failure = FALSE)
summary <- as.data.frame(result)
write.csv(summary[setdiff(names(summary), "result")], file.path(out, "test_summary.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
stopifnot(sum(summary$failed) == 0L, sum(summary$error) == 0L,
          sum(summary$warning) == 0L)
cat("Passed", sum(summary$passed), "current-package expectations.\n")
