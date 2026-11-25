args <- commandArgs(TRUE)
stopifnot(length(args) >= 3L)
output <- normalizePath(args[1], mustWork = FALSE)
replications <- as.integer(args[2])
library_path <- normalizePath(args[3])
cores <- if (length(args) >= 4L) as.integer(args[4]) else 3L
stopifnot(replications >= 20L, replications <= 1000L)
dir.create(output, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(library_path, .libPaths()))
library(inferMM)
stopifnot(as.character(packageVersion("inferMM")) == "0.1.2")
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
cases <- expand.grid(n = c(50L, 200L), error = c("Normal", "Gamma4"),
                     stringsAsFactors = FALSE)
cases$case_id <- seq_len(nrow(cases))
cases$replications <- replications
write.csv(cases, file.path(output, "case_manifest.csv"), row.names = FALSE)
source_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
writeLines(capture.output(sessionInfo()), file.path(output, "sessionInfo.txt"))
write.csv(data.frame(file = basename(source_file), md5 = unname(tools::md5sum(source_file))),
          file.path(output, "source_fingerprint.csv"), row.names = FALSE)

mean_mm <- function(x, V, K) V*x/(K+x)
gradient_mm <- function(x, V, K) cbind(Vmax=x/(K+x), Km=-V*x/(K+x)^2)
mean_log <- function(x, logV, logK) {
  value <- exp(logV)*x/(exp(logK)+x)
  attr(value,"gradient") <- cbind(logV=value, logK=-value*exp(logK)/(exp(logK)+x))
  value
}

independent_information <- function(x, theta, skewness, kurtosis) {
  v <- theta[3]*sqrt(x)
  g <- gradient_mm(x, theta[1], theta[2])
  information <- matrix(0, 3, 3)
  for (i in seq_along(x)) {
    omega <- matrix(c(v[i], skewness[i]*v[i]^1.5,
                      skewness[i]*v[i]^1.5, (kurtosis[i]-1)*v[i]^2), 2, 2)
    derivative <- rbind(c(g[i, ], 0), c(0, 0, sqrt(x[i]))) %*% diag(theta)
    information <- information + crossprod(derivative, solve(omega, derivative))
  }
  diag(theta) %*% solve(information) %*% diag(theta)
}

fit_wnls <- function(x, y) {
  h <- sqrt(x)
  k0 <- median(unique(x))
  v0 <- max(y)*(k0+max(x))/max(x)
  fit <- nls(y ~ mean_log(x, logV, logK), weights=1/h,
             start=list(logV=log(v0), logK=log(k0)),
             control=nls.control(maxiter=200L, tol=1e-6,
                                 minFactor=2^-20, warnOnly=TRUE))
  if (!isTRUE(fit$convInfo$isConv)) stop("WNLS did not converge")
  beta <- setNames(exp(coef(fit)), c("Vmax", "Km"))
  residual <- y-mean_mm(x, beta[1], beta[2])
  gamma <- mean(residual^2/h)
  g <- gradient_mm(x, beta[1], beta[2])
  covariance <- gamma*solve(crossprod(g/sqrt(h)))
  list(coefficients=beta, vcov=covariance, scale=gamma,
       convergence=list(converged=TRUE,
                        score_norm=max(abs(colSums(g*residual/h))/
                          sqrt(colSums(g^2/h)*sum(residual^2/h))),
                        iterations=fit$convInfo$finIter))
}

attempt <- function(method, x, y, error, validate) {
  warnings <- character()
  failure <- NULL
  start <- proc.time()[[3]]
  fit <- tryCatch(withCallingHandlers({
    if (method == "WNLS") {
      fit_wnls(x, y)
    } else if (method == "Supplied") {
      fit_mm_moments(x, y, variance="sqrt",
                     moments=if (error == "Normal") "normal" else "gamma",
                     shape=if (error == "Normal") NULL else 4,
                     covariance="model")
    } else {
      fit_mm_moments(x, y, variance="sqrt", moments="estimate", covariance="model")
    }
  }, warning=function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error=function(e) {
    failure <<- paste(class(e)[1], conditionMessage(e), sep=": ")
    NULL
  })
  out <- data.frame(method=method, success=FALSE, Vmax=NA_real_, Km=NA_real_,
                    gamma=NA_real_, seV=NA_real_, seK=NA_real_,
                    lowerV=NA_real_, upperV=NA_real_, lowerK=NA_real_, upperK=NA_real_,
                    skewness=NA_real_, kurtosis=NA_real_, score_norm=NA_real_,
                    iterations=NA_integer_, independent_cov_rel_error=NA_real_,
                    elapsed_seconds=proc.time()[[3]]-start,
                    warning_count=length(warnings), warnings=paste(warnings,collapse=" | "),
                    error=if (is.null(failure)) "" else failure)
  if (is.null(fit)) return(out)
  beta <- unname(fit$coefficients)
  se <- sqrt(diag(fit$vcov))
  if (any(!is.finite(c(beta,se,fit$scale))) || any(c(beta,se,fit$scale)<=0)) {
    out$error <- "Invalid parameter or covariance after fitting"
    return(out)
  }
  if (method != "WNLS") {
    out$skewness <- mean(fit$standardized_moments$skewness)
    out$kurtosis <- mean(fit$standardized_moments$kurtosis)
    out$score_norm <- fit$convergence$score_norm
    if (validate) {
      independent <- independent_information(x, c(beta,fit$scale),
          fit$standardized_moments$skewness, fit$standardized_moments$kurtosis)
      out$independent_cov_rel_error <- max(abs(independent-fit$vcov_joint_model)/
                                                pmax(1e-12,abs(independent)))
      stopifnot(out$independent_cov_rel_error < 1e-7)
    }
  }
  out$success <- TRUE
  out$Vmax <- beta[1]; out$Km <- beta[2]; out$gamma <- fit$scale
  out$seV <- se[1]; out$seK <- se[2]
  out$score_norm <- fit$convergence$score_norm
  q <- qnorm(0.975)
  out$lowerV <- beta[1]-q*se[1]; out$upperV <- beta[1]+q*se[1]
  out$lowerK <- beta[2]-q*se[2]; out$upperK <- beta[2]+q*se[2]
  out$iterations <- fit$convergence$iterations
  out
}

run_case <- function(index) {
  case <- cases[index, ]
  x <- seq(1, 100, length.out=case$n)
  mu <- mean_mm(x, 100, 20)
  sigma <- sqrt(1.5*sqrt(x))
  rows <- vector("list", replications)
  for (replicate in seq_len(replications)) {
    seed <- 20261024L+10000L*case$case_id+replicate
    set.seed(seed)
    z <- if (case$error == "Normal") rnorm(case$n) else (rgamma(case$n,shape=4)-4)/2
    y <- mu+sigma*z
    fitted <- lapply(c("WNLS","Supplied","Estimated"), attempt,
                     x=x, y=y, error=case$error, validate=replicate<=10L)
    result <- do.call(rbind, fitted)
    result$case_id <- case$case_id; result$n <- case$n; result$error_law <- case$error
    result$replicate <- replicate; result$seed <- seed
    result$sum_y <- sum(y); result$sum_y_squared <- sum(y^2)
    rows[[replicate]] <- result
    if (replicate %% 100L == 0L || replicate == replications) {
      write.csv(do.call(rbind,rows[seq_len(replicate)]),
                file.path(output,sprintf("case_%02d_raw.csv",index)),row.names=FALSE)
      cat(sprintf("case %d (%s, n=%d): %d/%d datasets\n",
                  index,case$error,case$n,replicate,replications))
      flush.console()
    }
  }
  do.call(rbind, rows)
}

started <- Sys.time()
if (.Platform$OS.type != "windows" && cores > 1L) {
  results <- parallel::mclapply(seq_len(nrow(cases)),run_case,
                              mc.cores=min(cores,nrow(cases)),mc.set.seed=FALSE)
} else {
  results <- lapply(seq_len(nrow(cases)),run_case)
}
if (any(vapply(results,inherits,logical(1),"try-error"))) stop("A case worker failed")
raw <- do.call(rbind,results)
write.csv(raw,file.path(output,"raw.csv"),row.names=FALSE)

summaries <- list()
for (id in cases$case_id) for (method in c("WNLS","Supplied","Estimated")) {
  a <- raw[raw$case_id==id & raw$method==method,]
  usable <- a[a$success,]
  reference <- raw[raw$case_id==id & raw$method=="WNLS" & raw$success,]
  paired <- merge(usable,reference,by="replicate",suffixes=c("","_ref"))
  row <- data.frame(case_id=id,n=a$n[1],error_law=a$error_law[1],method=method,
                     attempted=nrow(a),successful=nrow(usable),paired_n=nrow(paired),
                     warnings=sum(a$warning_count),failures=sum(!a$success))
  for (par in c("Vmax","Km")) {
    truth <- if (par=="Vmax") 100 else 20
    letter <- if (par=="Vmax") "V" else "K"
    value <- usable[[par]]
    sd_value <- sd(value)
    se <- usable[[paste0("se",letter)]]
    cover <- usable[[paste0("lower",letter)]]<=truth &
      usable[[paste0("upper",letter)]]>=truth
    width <- usable[[paste0("upper",letter)]]-usable[[paste0("lower",letter)]]
    u <- (paired[[par]]-mean(paired[[par]]))^2
    v <- (paired[[paste0(par,"_ref")]]-mean(paired[[paste0(par,"_ref")]]))^2
    ratio <- mean(v)/mean(u)
    ratio_se <- sd(v/mean(v)-u/mean(u))/sqrt(nrow(paired))
    values <- c(bias=mean(value)-truth,rmse=sqrt(mean((value-truth)^2)),
                empirical_sd=sd_value,mean_se=mean(se),se_sd_ratio=mean(se)/sd_value,
                coverage=mean(cover),coverage_mcse=sqrt(mean(cover)*(1-mean(cover))/length(cover)),
                mean_width=mean(width),paired_variance_ratio=ratio,
                ratio_mc_lower=ratio*exp(-qnorm(.975)*ratio_se),
                ratio_mc_upper=ratio*exp(qnorm(.975)*ratio_se))
    for (name in names(values)) row[[paste0(letter,"_",name)]] <- values[name]
  }
  summaries[[length(summaries)+1L]] <- row
}
summary <- do.call(rbind,summaries)
write.csv(summary,file.path(output,"summary.csv"),row.names=FALSE)
theory <- lapply(seq_len(nrow(cases)),function(i) {
  x <- seq(1,100,length.out=cases$n[i]); h <- sqrt(x)
  g <- gradient_mm(x,100,20)
  wnls <- solve(crossprod(g/sqrt(1.5*h)))
  s <- if (cases$error[i]=="Normal") 0 else 1
  k <- if (cases$error[i]=="Normal") 3 else 4.5
  joint <- independent_information(x,c(100,20,1.5),rep(s,length(x)),rep(k,length(x)))
  data.frame(case_id=i,n=length(x),error_law=cases$error[i],
             V_wnls_variance=wnls[1,1],V_efficient_variance=joint[1,1],
             V_ratio=wnls[1,1]/joint[1,1],K_wnls_variance=wnls[2,2],
             K_efficient_variance=joint[2,2],K_ratio=wnls[2,2]/joint[2,2])
})
write.csv(do.call(rbind,theory),file.path(output,"theoretical_variances.csv"),row.names=FALSE)
writeLines(c(paste("started:",started),paste("finished:",Sys.time()),
             paste("datasets:",nrow(raw)/3),paste("fits:",nrow(raw))),
           file.path(output,"execution_record.txt"))
print(summary[,c("n","error_law","method","successful","V_rmse","K_rmse",
                  "V_coverage","K_coverage","V_paired_variance_ratio","K_paired_variance_ratio")],
      row.names=FALSE)
