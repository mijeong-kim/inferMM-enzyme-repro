# Regenerate manuscript tables from the executed benchmark, without hand entry.
arg <- grep("^--file=", commandArgs(FALSE), value=TRUE)[1]
root <- dirname(dirname(normalizePath(sub("^--file=", "", arg))))
dat <- read.csv(file.path(root,"analysis/results/revision_benchmark_summary.csv"))
methods <- c("Constant","Log","Sqrt","Cuberoot","AIC screen","Constant HC3","GNLS power")
labels <- c("Constant","Log","Square root","Cube root","AIC screen","Constant + HC3","GNLS power")
fmt <- function(x) sprintf("%.3f",x)
blue <- "\\ifshowrevisions\\color{blue}\\fi"

main <- c("\\begin{table}[H]", "\\centering", blue,
 "\\caption{\\rev{Additional constant-variance experiment under Normal errors. Each design uses 1,000 datasets with true variance 6. CP denotes coverage of a nominal 95\\% interval. The constant and HC3 rows share point estimates. All fits in these two Normal scenarios returned finite estimates.}}",
 "\\label{tab:revision-homoscedastic}", "\\small",
 "\\setlength{\\tabcolsep}{3pt}",
 "\\begin{tabular}{lrrrr}", "\\toprule",
 "Method & RMSE($V_{\\max}$) & RMSE($K_m$) & CP($V_{\\max}$) & CP($K_m$) \\\\", "\\midrule")
for(di in c("Grid20","SDL28")) {
  main <- c(main, paste0("\\multicolumn{5}{l}{",if(di=="Grid20")
      "$n=20$: 20 equally spaced concentrations" else "$n=28$: seven levels, four replicates", "} \\\\"))
  for(k in seq_along(methods)) {
    z <- subset(dat, design==di & scenario=="Constant" & error=="Normal" & method==methods[k])
    main <- c(main,paste0(paste(c(labels[k],fmt(unlist(z[c("rmseV","rmseK","cpV","cpK")]))),
                        collapse=" & ")," \\\\"))
  }
  if(di=="Grid20") main <- c(main,"\\midrule")
}
main <- c(main,"\\bottomrule","\\end{tabular}","\\end{table}")
writeLines(main,file.path(root,"sections/benchmark_main_table_rev1.tex"))

supp <- character()
for(di in c("Grid20","SDL28")) for(er in c("Normal","Skewed")) {
  supp <- c(supp,"\\clearpage","\\begin{landscape}","\\begin{table}[H]",
   "\\centering",blue,
   paste0("\\caption{\\rev{Additional benchmark: ",di,", ",er,
    " errors. Each scenario has 1,000 attempted replications. N is the number of successful fits; summaries condition on success. CP denotes 95\\% coverage; MIL denotes mean interval length. Failed fits are retained in the replication-level output.}}"),
   paste0("\\label{tab:revision-",di,"-",er,"}"),"\\scriptsize",
   "\\setlength{\\tabcolsep}{3pt}","\\begin{tabular}{llrrrrrrrrr}","\\toprule",
   "Scenario & Method & N & RMSE($V_{\\max}$) & RMSE($K_m$) & CP($V_{\\max}$) & CP($K_m$) & MIL($V_{\\max}$) & MIL($K_m$) & Var\\_MSE & MCSE(CP$_{K_m}$) \\\\","\\midrule")
  for(sc in c("Constant","MM","Exp","Hill")) {
    for(k in seq_along(methods)) {
      z <- subset(dat,design==di & scenario==sc & error==er & method==methods[k])
      values <- c(if(k==1) sc else "",labels[k],z$successful,
        fmt(unlist(z[c("rmseV","rmseK","cpV","cpK","widthV","widthK","var_mse","mcseK")])))
      supp <- c(supp,paste0(paste(values,collapse=" & ")," \\\\"))
    }
    if(sc!="Hill") supp <- c(supp,"\\midrule")
  }
  supp <- c(supp,"\\bottomrule","\\end{tabular}","\\end{table}","\\end{landscape}")
}
writeLines(supp,file.path(root,"sections/benchmark_supp_tables_rev1.tex"))
source(file.path(root,"analysis/highlight_result_tables.R"))
highlight_result_tables(file.path(root,"sections/benchmark_main_table_rev1.tex"),
                        highlight_coverage = TRUE)
highlight_result_tables(file.path(root,"sections/benchmark_supp_tables_rev1.tex"),
                        highlight_coverage = TRUE, highlight_variance_mse = TRUE)
cat("Generated main table and four supplementary tables.\n")
