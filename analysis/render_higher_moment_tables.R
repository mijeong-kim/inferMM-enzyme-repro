args <- commandArgs(TRUE)
file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])
home <- dirname(dirname(normalizePath(file)))
input <- if (length(args)) args[1] else file.path(home,"outputs/higher_moments")
target <- if (length(args)>1L) args[2] else file.path(home,"sections/higher_moment_tables_rev1.tex")
s <- read.csv(file.path(input,"summary.csv"))
stopifnot(nrow(s)==12L,all(s$attempted==1000L))
labels <- c(WNLS="WNLS",Supplied="Supplied moments",Estimated="Estimated moments")
fmt <- function(x) formatC(x,format="f",digits=3)
fmt_best <- function(x, score = identity) {
  shown <- fmt(x)
  scores <- score(as.numeric(shown))
  best <- scores == min(scores)
  shown[best] <- paste0("\\textbf{", shown[best], "}")
  shown
}
fmt_closest <- function(x, target) {
  shown <- fmt(x)
  distance <- abs(round(1000 * as.numeric(shown)) - round(1000 * target))
  best <- distance == min(distance)
  shown[best] <- paste0("\\rev{\\textbf{", shown[best], "}}")
  shown
}
row <- function(x) paste0(paste(x,collapse=" & ")," \\\\")
lines <- character()
start <- function(caption,label,columns,header) {
  lines <<- c(lines,"\\clearpage","\\begin{landscape}","\\begin{table}[H]",
    "\\centering","\\ifshowrevisions\\color{blue}\\fi",
    paste0("\\caption{\\rev{",caption,"}}"),paste0("\\label{",label,"}"),
    "\\footnotesize","\\setlength{\\tabcolsep}{3pt}",
    paste0("\\begin{tabular}{",columns,"}"),"\\toprule",row(header),"\\midrule")
}
finish <- function(note) {
  lines <<- c(lines,"\\bottomrule","\\end{tabular}","\\par\\medskip",
    paste0("\\begin{minipage}{0.96\\linewidth}\\footnotesize ",note,"\\end{minipage}"),
    "\\end{table}","\\end{landscape}")
}
start("Direct higher-moment comparison: point estimation. Every method uses the same Michaelis--Menten mean and correct fixed variance shape, with 1,000 shared datasets per cell. SD is the empirical standard deviation; VR is the paired WNLS variance divided by the method variance. Bold entries indicate the smallest absolute bias, RMSE and SD, and the largest VR for each condition and parameter, including ties at the displayed precision.",
  "tab:higher-moment-point","rllrrrrrrrr",
  c("$n$","Errors","Method","Bias($V_{\\max}$)","RMSE($V_{\\max}$)","SD($V_{\\max}$)",
    "Bias($K_m$)","RMSE($K_m$)","SD($K_m$)","VR($V_{\\max}$)","VR($K_m$)"))
for (id in c(1L,3L,2L,4L)) {
  condition <- s[s$case_id==id,]
  point <- cbind(
    fmt_best(condition$V_bias, abs), fmt_best(condition$V_rmse),
    fmt_best(condition$V_empirical_sd), fmt_best(condition$K_bias, abs),
    fmt_best(condition$K_rmse), fmt_best(condition$K_empirical_sd),
    fmt_best(condition$V_paired_variance_ratio, function(x) -x),
    fmt_best(condition$K_paired_variance_ratio, function(x) -x)
  )
  for (method in c("WNLS","Supplied","Estimated")) {
    index <- match(method, condition$method)
    a <- condition[index,]
    lines <- c(lines,row(c(if(method=="WNLS") a$n else "",
      if(method=="WNLS") if(a$error_law=="Normal") "Normal" else "Gamma(4)" else "",
      labels[method],point[index,])))
  }
  if(id!=4L) lines <- c(lines,"\\midrule")
}
finish("Supplied and Estimated moments use \\texttt{fit\\_mm\\_moments()} with default \\texttt{covariance=\"model\"}. Supplied moments are the true population values and provide an additional-information reference; Estimated moments are obtained from pilot residuals. Gamma errors are centered and standardized. A VR greater than 1 indicates smaller sampling variance than WNLS on paired successful datasets. All conditions and methods are reported.")
start("Direct higher-moment comparison: model-based standard errors and 95\\% Wald intervals. CP is empirical coverage; SE/SD is mean reported standard error divided by empirical standard deviation; MIL is mean interval length. Bold values mark CP closest to 0.95 and SE/SD closest to 1 for each condition and parameter, including ties at the displayed precision. At 1, the mean reported standard error matches the empirical standard deviation.",
  "tab:higher-moment-inference","rllrrrrrrr",
  c("$n$","Errors","Method","Usable","CP($V_{\\max}$)","CP($K_m$)",
    "SE/SD($V_{\\max}$)","SE/SD($K_m$)","MIL($V_{\\max}$)","MIL($K_m$)"))
for (id in c(1L,3L,2L,4L)) {
  condition <- s[s$case_id==id,]
  calibration <- cbind(
    fmt_closest(condition$V_coverage, 0.95), fmt_closest(condition$K_coverage, 0.95),
    fmt_closest(condition$V_se_sd_ratio, 1), fmt_closest(condition$K_se_sd_ratio, 1)
  )
  for (method in c("WNLS","Supplied","Estimated")) {
    index <- match(method, condition$method)
    a <- condition[index,]
    lines <- c(lines,row(c(if(method=="WNLS") a$n else "",
      if(method=="WNLS") if(a$error_law=="Normal") "Normal" else "Gamma(4)" else "",
      labels[method],a$successful,calibration[index,],fmt(a$V_mean_width),fmt(a$K_mean_width))))
  }
  if(id!=4L) lines <- c(lines,"\\midrule")
}
finish("All 1,000 attempted fits per method and cell were usable and converged, with no warnings. Intervals use the same standard-normal critical value. WNLS uses RSS/$n$; higher-moment fits use the mean block of the joint inverse-information covariance. Coverage MCSE is about 0.007 near 0.95 (at most 0.016). Raw results include every attempt, seed, warning, root diagnostic and interval. Paired variance-ratio Monte Carlo uncertainty and true-model asymptotic variances are retained with the code and data.")
dir.create(dirname(target),recursive=TRUE,showWarnings=FALSE)
writeLines(lines,target)
cat("Wrote two higher-moment tables from",nrow(s),"rows.\n")
