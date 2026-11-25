# Generate only new comparison tables from executed replication-level summaries.
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
root <- dirname(dirname(normalizePath(sub("^--file=", "", file_arg))))
args <- commandArgs(TRUE)
input <- if (length(args)) normalizePath(args[1]) else
  file.path(root, "analysis/results/matched_final")
setwd(root)
s <- read.csv(file.path(input, "summary.csv"))
eq <- read.csv(file.path(input, "equivalence_summary.csv"))
pc <- read.csv(file.path(input, "paired_coverage.csv"))
stopifnot(all(s$attempted==1000L),nrow(s)==46L)
labels <- c(Profile="Profile",WNLS_low="WNLS: low start",
  WNLS_central="WNLS: central start",WNLS_high="WNLS: high start",
  Profile_sqrt="Profile sqrt",Profile_screen="Profile screen",
  Profile_screen_mu2="Profile screen + $\\widehat{\\mu}^2$",
  GNLS_power="GNLS power",TBS_log="Log-TBS (mean corrected)",
  Pooled_sqrt="Pooled sqrt",Cluster_profile="Cluster profile",NLME_matched="NLME matched")
fmt <- function(x,d=3) ifelse(is.finite(x),formatC(x,format="f",digits=d),"--")
fmt_best <- function(x, d = 3, target = 0) {
  shown <- fmt(x, d)
  valid <- is.finite(x)
  distance <- rep(Inf, length(x))
  scale <- 10^d
  # Compare in integer units of the last displayed digit to keep ties exact.
  distance[valid] <- abs(round(as.numeric(shown[valid]) * scale) - round(target * scale))
  best <- valid & distance == min(distance)
  shown[best] <- paste0("\\textbf{", shown[best], "}")
  shown
}
row <- function(x) paste0(paste(x,collapse=" & ")," \\\\")
lines <- character()
table_start <- function(caption,label,columns,header) {
  lines <<- c(lines,"\\clearpage","\\begin{landscape}","\\begin{table}[H]",
    "\\centering","\\ifshowrevisions\\color{blue}\\fi",
    paste0("\\caption{\\rev{",caption,"}}"),paste0("\\label{",label,"}"),
    "\\footnotesize","\\setlength{\\tabcolsep}{4pt}",
    paste0("\\begin{tabular}{",columns,"}"),"\\toprule",row(header),"\\midrule")
}
table_end <- function(note="") {
  lines <<- c(lines,"\\bottomrule","\\end{tabular}",
    if(nzchar(note)) c("\\par\\medskip",paste0("\\begin{minipage}{0.96\\linewidth}\\footnotesize ",note,"\\end{minipage}")),
    "\\end{table}","\\end{landscape}")
}
table_start("Same-weight numerical equivalence and computation, SDL28. Each variance shape has 1,000 datasets shared by all four implementations. Both methods use the residual scale RSS/$n$ for the reported covariance comparison. Bold values mark the shortest median CPU time for each variance shape, including ties at the displayed precision.",
  "tab:matched-weighting","llrrrrr",c("Shape","Implementation","Usable","Converged","CPU (ms)","Max rel. $\\Delta\\beta$","Max rel. $\\Delta$SE"))
for(id in 1:2) {
  group <- s[s$id == id, ]
  cpu <- fmt_best(group$cpu_median_ms, d = 1)
  for(m in c("Profile","WNLS_low","WNLS_central","WNLS_high")) {
    i <- match(m, group$method)
    z <- group[i, ]; e<-eq[eq$id==id&eq$method==m,]
    lines<-c(lines,row(c(if(m=="Profile") z$scenario else "",labels[m],z$usable,z$converged,
      cpu[i],if(nrow(e)) formatC(e$max_rel_beta,format="e",digits=2) else "--",
      if(nrow(e)) formatC(e$max_rel_se,format="e",digits=2) else "--")))
  }
  if(id==1) lines<-c(lines,"\\midrule")
}
table_end("Relative differences are maxima over datasets where both the profile fit and the indicated WNLS fit converged. WNLS starts use 0.1, 1, or 10 times the median distinct concentration for $K_m$, with the corresponding data-based $V_{\\max}$ start. CPU is median process CPU time per attempted fit, not an algorithmic complexity bound. Timings are aggregated from the attempt-level records of the archived reference execution. The profile implementation evaluates a 400-point root-search grid.")

table_start("Point estimation and variance recovery in the transformation comparison. Both generating mechanisms have the same original-scale conditional mean parameters $(V_{\\max},K_m)=(100,20)$. Log-TBS includes estimated lognormal mean correction. Bold entries indicate the smallest absolute bias and RMSE within each condition, including ties at the displayed precision.",
  "tab:matched-transform-point","lllrrrrr",c("Design","Mechanism","Method","Bias($V_{\\max}$)","RMSE($V_{\\max}$)","Bias($K_m$)","RMSE($K_m$)","Var RMSE"))
transform_methods <- c("Profile_sqrt", "Profile_screen", "Profile_screen_mu2",
                       "GNLS_power", "TBS_log")
for(id in 3:6) {
  group <- s[s$id == id, ]
  metrics <- c("biasV", "rmseV", "biasK", "rmseK", "var_rmse")
  point <- vapply(group[metrics], fmt_best, character(nrow(group)))
  for(m in transform_methods) {
    i <- match(m, group$method)
    z <- group[i, ]
    lines<-c(lines,row(c(if(m=="Profile_sqrt") z$design else "",if(m=="Profile_sqrt") z$scenario else "",
      labels[m],point[i, ])))
  }
  if(id<6) lines<-c(lines,"\\midrule")
}
table_end("Var RMSE is the square root of the average squared error of the fitted original-scale variance over concentrations and usable replications. The square-root profile is correctly specified for Additive data and Log-TBS for Lognormal data. The expanded screen adds a variance shape proportional to the squared pilot mean curve to the original four candidates.")

table_start("Interval calibration and computation in the transformation comparison. CP is empirical coverage of nominal 95\\% Wald intervals. MIL is mean interval length. Every cell has 1,000 attempted datasets. Bold values mark the coverage closest to 95\\% for each parameter and the shortest median CPU time within each condition, including ties at the displayed precision.",
  "tab:matched-transform-inference","lllrrrrrrr",c("Design","Mechanism","Method","Usable","Conv.","CP($V_{\\max}$)","CP($K_m$)","MIL($V_{\\max}$)","MIL($K_m$)","CPU (ms)"))
for(id in 3:6) {
  group <- s[s$id == id, ]
  coverage <- vapply(group[c("cpV", "cpK")], fmt_best, character(nrow(group)), target = 0.95)
  cpu <- fmt_best(group$cpu_median_ms, d = 1)
  for(m in transform_methods) {
    i <- match(m, group$method)
    z <- group[i, ]
    lines<-c(lines,row(c(if(m=="Profile_sqrt") z$design else "",if(m=="Profile_sqrt") z$scenario else "",
      labels[m],z$usable,z$converged,coverage[i, ],fmt(z$widthV),fmt(z$widthK),cpu[i])))
  }
  if(id<6) lines<-c(lines,"\\midrule")
}
table_end("The expanded screen uses HC3 standard errors and $t_{n-2}$ critical values. The other methods use normal critical values. Both screens treat the selected shape as fixed, including the pilot curve when used. Log-TBS includes first-order log-variance uncertainty in the retransformed mean parameter. CPU includes four fits for the original screen and five for the expanded screen. Coverage MCSE is about 0.007 near 0.95. Paired GNLS and Log-TBS contrasts are given below.")


table_start("Matched clustered comparison: point estimation and variance recovery. All methods use $h(S)=\\sqrt S$; the clustered methods share the single additive random effect on $V_{\\max}$ and have true $\\tau^2=100$, $\\gamma=1.5$. Each cluster has the SDL28 layout. Bold entries indicate the smallest absolute bias and RMSE within each condition, including ties at the displayed precision.",
  "tab:matched-cluster-point","rllrrrrrr",c("Clusters","Innovations","Method","Bias($V_{\\max}$)","RMSE($V_{\\max}$)","Bias($K_m$)","RMSE($K_m$)","RMSE($\\tau^2$)","Var RMSE"))
for(id in 7:12) {
  group <- s[s$id == id, ]
  metrics <- c("biasV", "rmseV", "biasK", "rmseK", "rmseTau", "var_rmse")
  point <- vapply(group[metrics], fmt_best, character(nrow(group)))
  for(m in c("Pooled_sqrt","Cluster_profile","NLME_matched")) {
    i <- match(m, group$method)
    z <- group[i, ]
    lines<-c(lines,row(c(if(m=="Pooled_sqrt") z$clusters else "",if(m=="Pooled_sqrt") z$error else "",labels[m],
      point[i, ])))
  }
  if(id<12) lines<-c(lines,"\\midrule")
}
table_end("Skewed random effects and residual innovations are independent standardized Gamma(4,1) draws; no within-sample centering or standardization is used. Var RMSE compares marginal variances, including the random-effect contribution for both clustered fits. The pooled fit ignores dependence and does not estimate $\\tau^2$.")

table_start("Matched clustered comparison: interval calibration and computation. CP refers to nominal 95\\% intervals; each cell has 1,000 attempted datasets. Conv. records the implemented convergence criterion, not proof of a global optimum. Bold values mark the coverage closest to 95\\% for each parameter and the shortest median CPU time within each condition, including ties at the displayed precision.",
  "tab:matched-cluster-inference","rllrrrrrrr",c("Clusters","Innovations","Method","Usable","Conv.","CP($V_{\\max}$)","CP($K_m$)","MIL($V_{\\max}$)","MIL($K_m$)","CPU (ms)"))
for(id in 7:12) {
  group <- s[s$id == id, ]
  coverage <- vapply(group[c("cpV", "cpK")], fmt_best, character(nrow(group)), target = 0.95)
  cpu <- fmt_best(group$cpu_median_ms, d = 1)
  for(m in c("Pooled_sqrt","Cluster_profile","NLME_matched")) {
    i <- match(m, group$method)
    z <- group[i, ]
    lines<-c(lines,row(c(if(m=="Pooled_sqrt") z$clusters else "",if(m=="Pooled_sqrt") z$error else "",labels[m],
      z$usable,z$converged,coverage[i, ],fmt(z$widthV),fmt(z$widthK),cpu[i])))
  }
  if(id<12) lines<-c(lines,"\\midrule")
}
table_end("All intervals use normal cutoffs, including NLME. The profile uses the mean-information covariance approximation; NLME uses its package-reported fixed-effect covariance. CPU includes pooled-profile initialization for each clustered method. Timings are aggregated from the attempt-level records of the archived reference execution. Coverage and CPU summaries are reported for every method and cluster-count setting.")


table_start("Paired coverage differences on datasets with usable fits from both compared methods. Each difference is comparator minus reference. The Monte Carlo standard error (MCSE) uses the paired indicator differences, not two independent-binomial approximations.",
  "tab:matched-paired","llllrrrrr",c("Experiment","Design/$m$","Mechanism","Comparator","Pairs","$\\Delta$CP($V_{\\max}$)","MCSE","$\\Delta$CP($K_m$)","MCSE"))
for(id in 3:12) {
  methods<-if(id<=6) c("TBS_log","GNLS_power") else "NLME_matched"
  for(m in methods) {
    a<-pc[pc$id==id&pc$comparator==m&pc$parameter=="V",]
    b<-pc[pc$id==id&pc$comparator==m&pc$parameter=="K",]
    z<-s[s$id==id,][1,]
    lines<-c(lines,row(c(if(id<=6) "Single curve" else "Clustered",
      if(id<=6) z$design else z$clusters,if(id<=6) z$scenario else z$error,labels[m],
      a$paired_n,fmt(a$cp_difference),fmt(a$mcse),fmt(b$cp_difference),fmt(b$mcse))))
  }
  if(id==6) lines<-c(lines,"\\midrule")
}
table_end("Reference: Profile screen for independent curves, Cluster profile for clustered data. A positive difference means more coverage, not necessarily better calibration if coverage exceeds 0.95. Raw output also retains pooled-cluster contrasts, model-selection choices, all warnings, and all failures.")
writeLines(lines,"sections/matched_comparison_tables_rev1.tex")
cat("Wrote six matched-comparison supplementary tables from",nrow(s),"summary rows.\n")
