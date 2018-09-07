#'
#' function to calculate all the internal consitency estimates
#' takes as input both datasets, either matrix or frame, and covariance matrices
#' so far input of a covariance matrix is not supported
#'
#' @export
brel <- function(x, boot.n = 200, interval = .95, n.iter = 2e3, n.burnin = 50,
                estimates = c("alpha", "lambda2", "lambda6", "glb", "omega"), supr.warnings = TRUE,
                omega.freq.method = "pa", omega.conf.int.type = "boot", omega.bay.cov.samp = FALSE,
                prior.samp = FALSE, item.dropped = FALSE, alpha.int.analytic = FALSE,
                bayes = TRUE, freq = TRUE, para.boot = TRUE, boot.interval.type = "basic", jags = FALSE) {
  if (supr.warnings) {
    options(warn = - 1)
  }
  estimates <- match.arg(estimates, several.ok = T)
  default <- c("alpha", "lambda2", "lambda6", "glb", "omega")
  mat <- match(default, estimates)
  estimates <- estimates[mat]
  estimates <- estimates[!is.na(estimates)]

  sum.res <- list()
  sum.res$call <- match.call()

  if (sum(is.na(x)) > 0) {
    return("missing values in data detected, please remove and run again")
  }
  data <- NULL
  sigma <- NULL
  if (ncol(x) == nrow(x)){
    return("so far input of a covariance matrix is not supported")
    if (sum(v[lower.tri(v)] != t(v)[lower.tri(v)]) > 0) {return("input matrix is not symmetric")}
    if (sum(eigen(x)$values < 0) > 0) {return("input matrix is not positive definite")}
    if (freq) {return("bootstrap confidence interval estimation requires a dataset")}
    if ("omega" %in% estimates) {return("omega can only be calculated with a dataset as input")}
    sigma <- x
  } else{
    data <- scale(x, scale = F)
    sigma <- cov(data)
  }

  if("glb" %in% estimates){
    control <- Rcsdp:::csdp.control(printlevel = 0)
    Rcsdp:::write.control.file(control)
  }
  if (bayes){
    if (jags){
      sum.res$bay <- jagsFun(data, n.iter, n.burnin, estimates, interval, omega.bay.cov.samp)
      sum.res$omega.pa <- omega.bay.cov.samp
    }
    else{
      sum.res$bay <- gibbsFun(data, n.iter, n.burnin, estimates, interval, omega.bay.cov.samp, item.dropped)
      sum.res$omega.pa <- omega.bay.cov.samp
    }
  }
  sum.res$freq.true <- FALSE
  if(freq){
    if (para.boot){
      sum.res$freq <- freqFun_para(data, boot.n, estimates, interval, omega.freq.method, omega.conf.int.type, item.dropped,
                                   alpha.int.analytic)
    } else{
    sum.res$freq <- freqFun_nonpara(data, boot.n, estimates, interval, omega.freq.method, omega.conf.int.type, item.dropped,
                                    alpha.int.analytic)
    }
    sum.res$freq.true <- TRUE
    sum.res$omega.freq.method <- omega.freq.method
    sum.res$omega.conf.int.type <- omega.conf.int.type
    sum.res$alpha.int.analytic <- alpha.int.analytic
    if (omega.freq.method == "pa" && omega.conf.int.type == "alg"){
      sum.res$omega.conf.type <- "boot"
        print("algebraic confidence interval for omega not available with method PA")
    }
  }
  sum.res$item.dropped <- item.dropped
  if("glb" %in% estimates)
    unlink("param.csdp")

  if (prior.samp) {
    sum.res$priors <- priorSamp(ncol(data), estimates)
  }


  sum.res$estimates <- estimates
  sum.res$n.iter <- n.iter
  sum.res$n.burnin <- n.burnin
  sum.res$boot.interval.type <- boot.interval.type
  sum.res$interval <- interval

  class(sum.res) = 'bayesrel'
  options(warn = 0)
  return(sum.res)
}
