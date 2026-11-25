# Run from BiomJ_rev: Rscript analysis/matched_comparisons.R [1000] [run-name]
args <- commandArgs(trailingOnly = TRUE)
B <- if (length(args)) as.integer(args[1]) else 1000L
tag <- if (length(args) > 1L) args[2] else "final"
stopifnot(is.finite(B), B >= 10L, grepl("^[A-Za-z0-9_-]+$", tag))
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
root <- dirname(dirname(normalizePath(sub("^--file=", "", file_arg))))
.libPaths(c(file.path(root, "analysis", "Rlib"), .libPaths()))
stopifnot(requireNamespace("nlme", quietly = TRUE),
          requireNamespace("drc", quietly = TRUE))
out <- file.path(root, "analysis", "results", paste0("matched_", tag))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
loader <- new.env(parent = globalenv())
sys.source(file.path(root, "examples/load_example_package.R"), loader)
loader$load_example_package(root)
pkg <- asNamespace("inferMM")
# Retain the original SDL28/SDL56 protocol when package data include eight levels.
levels <- c(7.8125, 15.625, 31.25, 62.5, 125, 250, 500)
levels <- 100 * levels / max(levels)
designs <- list(SDL28 = rep(levels, each = 4), SDL56 = rep(levels, each = 8))
truth <- c(Vmax = 100, Km = 20)

standard <- function(f, h = NULL) {
  list(beta = f$coefficients, covariance = f$vcov, variance = f$variance_hat,
       tau2 = if (is.null(f$tau2)) NA_real_ else unname(f$tau2),
       gamma = unname(f$scale), selected = f$variance_model,
       converged = if (is.null(f$converged)) TRUE else isTRUE(f$converged),
       iterations = if (is.null(f$iterations)) NA_integer_ else f$iterations,
       score = if (is.null(h)) NA_real_ else {
         J <- pkg$mm_grad(f$x, f$coefficients[1], f$coefficients[2])
         max(abs(colSums(J * (f$residuals / h))) /
               pmax(colSums(abs(J) * (abs(f$y) / h)), 1))
       })
}

profile <- function(x, y, shape) {
  f <- pkg$fit_mm(x, y, variance = shape)
  ans <- standard(f, f$variance_fun(x))
  ans$converged <- ans$score < 1e-5
  ans
}

screen <- function(x, y, include_mean_square = FALSE) {
  shapes <- c("constant", "log", "sqrt", "cuberoot")
  fit_shape <- function(h) {
    tryCatch(pkg$fit_mm(x, y, variance = h), error = function(e) NULL)
  }
  fs <- setNames(lapply(shapes, fit_shape), shapes)
  if (include_mean_square && !is.null(fs$sqrt)) {
    pilot_k <- unname(fs$sqrt$coefficients[2])
    shape <- pkg$variance_function(function(s) (s / (pilot_k + s))^2, label = "mean_square")
    fs["mean_square"] <- list(fit_shape(shape))
  }
  aic <- vapply(fs, function(f) if (is.null(f)) Inf else f$quasi_aic, numeric(1))
  if (!any(is.finite(aic))) stop("All screening candidates failed")
  f <- fs[[which.min(aic)]]
  h <- f$variance_fun(x)
  ans <- standard(f, h)
  ans$converged <- ans$score < 1e-5
  if (include_mean_square) {
    J <- pkg$mm_grad(x, ans$beta[1], ans$beta[2])
    weighted_J <- J / sqrt(h)
    inverse_information <- solve(crossprod(weighted_J))
    leverage <- rowSums((weighted_J %*% inverse_information) * weighted_J)
    if (any(!is.finite(leverage) | leverage >= 1)) stop("Invalid leverage for HC3 covariance")
    score <- J * (f$residuals / (h * (1 - leverage)))
    ans$covariance <- inverse_information %*% crossprod(score) %*% inverse_information
    ans$critical <- qt(.975, df = length(y) - 2L)
  }
  ans
}

wnls <- function(x, y, h, factor) {
  k0 <- median(unique(x)) * factor
  v0 <- max(y) * (k0 + max(x)) / max(x)
  f <- stats::nls(y ~ exp(a) * x / (exp(b) + x),
    data = data.frame(x, y, h), weights = 1/h,
    start = list(a = log(max(v0, 1e-6)), b = log(k0)),
    control = stats::nls.control(maxiter = 200, tol = 1e-6, warnOnly = TRUE))
  be <- exp(coef(f)); names(be) <- names(truth)
  r <- y - pkg$mm_mean(x, be[1], be[2]); ga <- mean(r^2/h)
  J <- pkg$mm_grad(x, be[1], be[2])
  # Match the package's RSS/n scale rather than the nls RSS/(n-2) default.
  cv <- ga * solve(crossprod(J, J/h))
  list(beta = be, covariance = cv, variance = ga*h, gamma = ga,
       converged = isTRUE(f$convInfo$isConv), iterations = f$convInfo$finIter,
       score = max(abs(colSums(J*(r/h))) / pmax(colSums(abs(J)*(abs(y)/h)), 1)))
}

tbs_log <- function(x, y) {
  if (any(y <= 0)) stop("Log-TBS requires strictly positive responses; no shift applied")
  k0 <- median(unique(x)); v0 <- max(y)*(k0+max(x))/max(x)
  f <- drc::drm(y ~ x, data = data.frame(x, y), fct = drc::MM.2(),
    bcVal = 0, start = c(v0, k0),
    control = drc::drmc(noMessage = TRUE, maxIt = 500, relTol = 1e-9))
  be <- unname(coef(f)); cv <- unname(vcov(f)); n <- length(y)
  logres <- log(y) - log(pkg$mm_mean(x, be[1], be[2]))
  s2 <- mean(logres^2)
  # A lognormal mean, not the inverse-transformed location, is the common target.
  locV <- be[1]; be[1] <- be[1]*exp(s2/2)
  delta <- diag(c(exp(s2/2), 1)); cv <- delta %*% cv %*% delta
  cv[1, 1] <- cv[1, 1] + (be[1]/2)^2 * 2*s2^2/n
  mu <- pkg$mm_mean(x, be[1], be[2])
  list(beta = be, covariance = cv, variance = mu^2*expm1(s2),
       gamma = s2, locationV = locV, converged = isTRUE(f$fit$convergence),
       selected = "log (lambda=0)")
}

gnls_power <- function(x, y) {
  k0 <- median(unique(x)); v0 <- max(y)*(k0+max(x))/max(x)
  f <- nlme::gnls(y ~ exp(a)*x/(exp(b)+x), data = data.frame(x, y),
    start = c(a=log(max(v0,1e-6)), b=log(k0)),
    weights = nlme::varPower(form = ~x),
    control = nlme::gnlsControl(maxIter=100, nlsMaxIter=50))
  be <- exp(coef(f)); p <- as.numeric(coef(f$modelStruct$varStruct,
                                          unconstrained=FALSE))*2
  D <- diag(be,2)
  list(beta=be, covariance=D%*%vcov(f)%*%D,
       variance=f$sigma^2*x^p, gamma=f$sigma^2, converged=TRUE,
       selected=paste0("power=",format(p,digits=8)))
}

cluster_profile <- function(d) {
  standard(pkg$cluster_mm(d, s="x", v="y", cluster="id", variance="sqrt"))
}

cluster_nlme <- function(d) {
  start <- pkg$fit_mm(d$x, d$y, variance="sqrt")$coefficients
  d$h <- sqrt(d$x)
  f <- nlme::nlme(y ~ Vmax*x/(exp(lK)+x), data=d,
    fixed=Vmax+lK~1, random=Vmax~1|id,
    start=c(Vmax=unname(start[1]),lK=log(unname(start[2]))),
    weights=nlme::varFixed(~h), method="ML",
    control=nlme::nlmeControl(pnlsTol=1e-3,msMaxIter=200,
                             pnlsMaxIter=30,returnObject=FALSE))
  be <- nlme::fixef(f); be[2] <- exp(be[2]); D <- diag(c(1,be[2]))
  tau2 <- as.numeric(nlme::VarCorr(f)["Vmax","Variance"])
  z <- d$x/(be[2]+d$x)
  list(beta=be, covariance=D%*%vcov(f)%*%D,
       variance=tau2*z^2+f$sigma^2*sqrt(d$x), tau2=tau2,
       gamma=f$sigma^2, converged=TRUE)
}

capture_fit <- function(fn, x, y, true_var, true_tau2=NA_real_) {
  warnings <- character(); error <- ""
  before <- proc.time()
  ans <- tryCatch(withCallingHandlers(fn(), warning=function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  }), error=function(e) {error <<- conditionMessage(e); NULL})
  timing <- proc.time()-before
  ok <- !is.null(ans) && length(ans$beta)==2 &&
    all(is.finite(ans$beta)) && all(ans$beta>0) &&
    all(is.finite(ans$covariance)) && all(diag(ans$covariance)>0)
  pick <- function(name, default=NA_real_) if(is.null(ans[[name]])) default else ans[[name]]
  be <- if(ok) as.numeric(ans$beta) else rep(NA_real_,2)
  se <- if(ok) sqrt(diag(ans$covariance)) else rep(NA_real_,2)
  converged <- ok && isTRUE(pick("converged",FALSE))
  critical <- pick("critical", qnorm(.975))
  data.frame(ok=ok,converged=converged,Vmax=be[1],Km=be[2],
    seV=se[1],seK=se[2],coverV=abs(be[1]-truth[1])<=critical*se[1],
    coverK=abs(be[2]-truth[2])<=critical*se[2],
    widthV=2*critical*se[1],widthK=2*critical*se[2],
    var_mse=if(ok) mean((ans$variance-true_var)^2) else NA_real_,
    tau2=pick("tau2"),gamma=pick("gamma"),true_tau2=true_tau2,
    locationV=pick("locationV"),score=pick("score"),
    objective=if(ok) sum((y-pkg$mm_mean(x,be[1],be[2]))^2/sqrt(x)) else NA_real_,
    iterations=pick("iterations"),selected=as.character(pick("selected","")),
    cpu_seconds=unname(sum(timing[c("user.self","sys.self")])),
    elapsed_seconds=unname(timing["elapsed"]),
    warning=paste(unique(warnings),collapse=" | "),fit_error=error,
    stringsAsFactors=FALSE)
}

if (Sys.getenv("MATCHED_DEFINITIONS_ONLY") != "1") {
cases <- rbind(
  data.frame(experiment="Weighting",design="SDL28",scenario=c("Constant","Sqrt"),error="Normal",clusters=0L),
  expand.grid(experiment="Transformation",design=names(designs),
              scenario=c("Additive","Lognormal"),error="Normal",clusters=0L,
              stringsAsFactors=FALSE),
  expand.grid(experiment="Cluster",design="SDL28",scenario="Sqrt",
              error=c("Normal","Skewed"),clusters=c(3L,6L,30L),
              stringsAsFactors=FALSE))
cases$id <- seq_len(nrow(cases)); cases$seed <- 907809200L+cases$id
write.csv(cases,file.path(out,"case_manifest.csv"),row.names=FALSE)
write.csv(data.frame(S=designs$SDL28),file.path(out,"SDL28_design.csv"),row.names=FALSE)

run_case <- function(i) {
  ca <- cases[i,]; x0 <- designs[[ca$design]]
  # Warm all code paths before timing, then reset this case's simulation stream.
  set.seed(1L)
  invisible(try(profile(x0,100*x0/(20+x0)+rnorm(length(x0)),"sqrt"),silent=TRUE))
  set.seed(ca$seed)
  rows <- vector("list",B)
  for(r in seq_len(B)) {
    x <- x0; n <- length(x); mu <- pkg$mm_mean(x,truth[1],truth[2])
    tau2 <- NA_real_
    if(ca$experiment=="Weighting") {
      h <- if(ca$scenario=="Constant") rep(1,n) else sqrt(x)
      v <- 1.5*h; y <- mu+sqrt(v)*rnorm(n)
      shape <- if(ca$scenario=="Constant") "constant" else "sqrt"
      jobs <- list(Profile=function() profile(x,y,shape),
        WNLS_low=function() wnls(x,y,h,.1),
        WNLS_central=function() wnls(x,y,h,1),
        WNLS_high=function() wnls(x,y,h,10))
    } else if(ca$experiment=="Transformation") {
      if(ca$scenario=="Additive") {
        v <- 1.5*sqrt(x); y <- mu+sqrt(v)*rnorm(n)
      } else {
        sigma <- .25; y <- mu*exp(sigma*rnorm(n)-sigma^2/2)
        v <- mu^2*expm1(sigma^2)
      }
      jobs <- list(Profile_sqrt=function() profile(x,y,"sqrt"),
        Profile_screen=function() screen(x,y),
        GNLS_power=function() gnls_power(x,y),
        TBS_log=function() tbs_log(x,y),
        Profile_screen_mu2=function() screen(x,y,include_mean_square=TRUE))
    } else {
      m <- ca$clusters; x <- rep(x0,m); id <- factor(rep(seq_len(m),each=n))
      draw <- function(N) if(ca$error=="Normal") rnorm(N) else (rgamma(N,shape=4)-4)/2
      b <- 10*draw(m); z <- x/(20+x)
      y <- (100+rep(b,each=n))*z+sqrt(1.5*sqrt(x))*draw(length(x))
      tau2 <- 100; v <- tau2*z^2+1.5*sqrt(x); d <- data.frame(x,y,id)
      jobs <- list(Pooled_sqrt=function() profile(x,y,"sqrt"),
                   Cluster_profile=function() cluster_profile(d),
                   NLME_matched=function() cluster_nlme(d))
    }
    # Rotate the order without consuming random numbers or favoring one method.
    order <- ((seq_along(jobs)+r-2L)%%length(jobs))+1L
    if(r==1L) invisible(lapply(jobs,function(fn) try(suppressWarnings(fn()),silent=TRUE)))
    rr <- lapply(order,function(j) cbind(method=names(jobs)[j],
      capture_fit(jobs[[j]],x,y,v,tau2)))
    rows[[r]] <- cbind(ca,replicate=r,n=length(x),min_y=min(y),do.call(rbind,rr))
  }
  raw <- do.call(rbind,rows); rownames(raw) <- NULL
  write.csv(raw,file.path(out,sprintf("case_%02d_raw.csv",i)),row.names=FALSE)
  cat("Finished",i,ca$experiment,ca$design,ca$scenario,ca$error,ca$clusters,"B=",B,"\n")
  raw
}

cores <- if(.Platform$OS.type=="windows") 1L else 3L
results <- parallel::mclapply(seq_len(nrow(cases)),run_case,mc.cores=cores,
                             mc.set.seed=FALSE,mc.preschedule=FALSE)
stopifnot(all(vapply(results,is.data.frame,logical(1))))
raw <- do.call(rbind,results)
write.csv(raw,file.path(out,"raw.csv"),row.names=FALSE)
summarize <- function(d) {
  z <- d[d$ok,]; nn <- nrow(z)
  avg <- function(a) if(length(a)&&any(is.finite(a))) mean(a,na.rm=TRUE) else NA_real_
  cpv <- avg(z$coverV); cpk <- avg(z$coverK)
  data.frame(d[1,c("id","experiment","design","scenario","error","clusters","method")],
    attempted=nrow(d),usable=nn,converged=sum(d$converged),warnings=sum(nzchar(d$warning)),
    biasV=avg(z$Vmax-100),rmseV=sqrt(avg((z$Vmax-100)^2)),
    biasK=avg(z$Km-20),rmseK=sqrt(avg((z$Km-20)^2)),
    cpV=cpv,cpK=cpk,mcseV=sqrt(cpv*(1-cpv)/nn),mcseK=sqrt(cpk*(1-cpk)/nn),
    widthV=avg(z$widthV),widthK=avg(z$widthK),var_rmse=sqrt(avg(z$var_mse)),
    biasTau=avg(z$tau2-z$true_tau2),rmseTau=sqrt(avg((z$tau2-z$true_tau2)^2)),
    mean_seV=avg(z$seV),emp_sdV=sd(z$Vmax),mean_seK=avg(z$seK),emp_sdK=sd(z$Km),
    cpu_median_ms=1000*median(d$cpu_seconds),elapsed_median_ms=1000*median(d$elapsed_seconds),
    boundary_tau=sum(z$tau2<=1.01e-8,na.rm=TRUE),row.names=NULL)
}
groups <- split(raw,interaction(raw$id,raw$method,drop=TRUE))
summary <- do.call(rbind,lapply(groups,summarize))
write.csv(summary,file.path(out,"summary.csv"),row.names=FALSE)
conv <- raw[raw$converged,]
write.csv(do.call(rbind,lapply(split(conv,interaction(conv$id,conv$method,drop=TRUE)),summarize)),
          file.path(out,"converged_only_summary.csv"),row.names=FALSE)
paired <- merge(raw[raw$experiment=="Weighting"&raw$method=="Profile",],
                raw[raw$experiment=="Weighting"&raw$method!="Profile",],
                by=c("id","replicate"),suffixes=c("_profile","_wnls"))
paired$relative_beta_difference <- pmax(abs(paired$Vmax_profile-paired$Vmax_wnls)/
  pmax(abs(paired$Vmax_profile),1),abs(paired$Km_profile-paired$Km_wnls)/pmax(abs(paired$Km_profile),1))
paired$relative_se_difference <- pmax(abs(paired$seV_profile-paired$seV_wnls)/
  pmax(abs(paired$seV_profile),1e-8),abs(paired$seK_profile-paired$seK_wnls)/pmax(abs(paired$seK_profile),1e-8))
write.csv(paired,file.path(out,"weighting_pairs.csv"),row.names=FALSE)
eq <- do.call(rbind,lapply(split(paired,interaction(paired$id,paired$method_wnls,drop=TRUE)),function(d) {
  z<-d[d$converged_profile&d$converged_wnls,]
  data.frame(id=d$id[1],shape=d$scenario_profile[1],method=d$method_wnls[1],paired_n=nrow(z),
    max_rel_beta=max(z$relative_beta_difference),max_rel_se=max(z$relative_se_difference),
    median_rel_beta=median(z$relative_beta_difference))
}))
write.csv(eq,file.path(out,"equivalence_summary.csv"),row.names=FALSE)
# Paired RMSE-error and coverage contrasts: comparator minus profile on shared fits.
contrasts <- list(); k <- 0L
for(i in cases$id[cases$experiment!="Weighting"]) {
  d<-raw[raw$id==i,]; ref<-if(d$experiment[1]=="Cluster") "Cluster_profile" else "Profile_screen"
  a<-d[d$method==ref,]
  for(m in setdiff(unique(d$method),ref)) {
    p<-merge(a,d[d$method==m,],by="replicate",suffixes=c("_ref","_alt"))
    p<-p[p$ok_ref&p$ok_alt,]; nn<-nrow(p)
    for(parameter in c("V","K")) {
      delta<-as.numeric(p[[paste0("cover",parameter,"_alt")]])-as.numeric(p[[paste0("cover",parameter,"_ref")]])
      k<-k+1L; contrasts[[k]]<-data.frame(id=i,reference=ref,comparator=m,
        parameter=parameter,paired_n=nn,cp_difference=mean(delta),mcse=sd(delta)/sqrt(nn))
    }
  }
}
write.csv(do.call(rbind,contrasts),file.path(out,"paired_coverage.csv"),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
inputs<-c(file.path(root,"analysis","matched_comparisons.R"),
          file.path(root,"packages","inferMM_0.1.2.tar.gz"))
write.csv(data.frame(file=sub(paste0(root,"/"),"",inputs,fixed=TRUE),
                    md5=unname(tools::md5sum(inputs))),file.path(out,"input_md5.csv"),row.names=FALSE)
cat("Completed",nrow(cases)*B,"datasets and",nrow(raw),"method attempts.\n")
}
