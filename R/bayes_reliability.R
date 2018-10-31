#'
#' calculate all the internal consistency estimates
#' takes as input both datasets, either matrix or frame, and covariance matrices
#' so far input of a covariance matrix is not supported
#'
#' @export
rel <- function(x, estimates = c("alpha", "lambda2", "lambda6", "glb", "omega"), interval = .95,
                omega.freq.method = "pa", omega.fit = FALSE,
                alpha.int.analytic = FALSE, n.obs = NULL,
                bayes = TRUE, freq = TRUE, cor.mat.out = TRUE,
                para.boot = FALSE, prior.samp = FALSE, item.dropped = FALSE,
                n.iter = 2e3, n.burnin = 50, boot.n = 1000, supr.warnings = TRUE) {
  if (supr.warnings) {
    options(warn = - 1)
  }
  estimates <- match.arg(estimates, several.ok = T)
  default <- c("alpha", "lambda2", "lambda4", "lambda6", "glb", "omega")
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
    if (is.null(n.obs)) {return("number of observations needs to be specified when entering a covariance matrix")}
    # return("so far input of a covariance matrix is not supported")
    if (sum(x[lower.tri(x)] != t(x)[lower.tri(x)]) > 0) {return("input matrix is not symmetric")}
    if (sum(eigen(x)$values < 0) > 0) {return("input matrix is not positive definite")}
    sigma <- x
    data <- MASS::mvrnorm(n.obs, rep(0, p), sigma, empirical = TRUE)
  } else{
    data <- scale(x, scale = F)
    sigma <- cov(data)
  }

  if("glb" %in% estimates){
    control <- Rcsdp:::csdp.control(printlevel = 0)
    Rcsdp:::write.control.file(control)
  }
  if (bayes){
    sum.res$bay <- gibbsFun(data, n.iter, n.burnin, estimates, interval, item.dropped, omega.fit)
  }
  if (omega.fit) {omega.freq.method <- "cfa"}
  if(freq){
    if (para.boot){
      sum.res$freq <- freqFun_para(data, boot.n, estimates, interval, omega.freq.method, item.dropped,
                                   alpha.int.analytic, omega.fit)
    } else{
      sum.res$freq <- freqFun_nonpara(data, boot.n, estimates, interval, omega.freq.method, item.dropped,
                                    alpha.int.analytic, omega.fit)
    }
    sum.res$omega.freq.method <- omega.freq.method
  }

  if("glb" %in% estimates)
    unlink("param.csdp")

  if (prior.samp) {
    sum.res$priors <- priorSamp(ncol(data), estimates)
  }

  sum.res$estimates <- estimates
  sum.res$n.iter <- n.iter
  sum.res$n.burnin <- n.burnin
  sum.res$interval <- interval
  if (cor.mat.out) {
    sum.res$cor.mat.out <- cor(data)
  }

  class(sum.res) = 'bayesrel'
  options(warn = 0)
  return(sum.res)
}
