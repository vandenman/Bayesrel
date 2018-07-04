#'
#' function to calculate all the internal consitency estimates
#'
#'
#' @export
brel <- function(raw.data, boot.n = 200, interval = .95, boot.interval.type = "basic",
                jags = FALSE, n.iter = 2e3, n.burnin = 50, freq = TRUE,
                estimates = c("alpha", "l2", "l4", "l6", "glb", "omega"), supr.warnings = TRUE,
                omega.freq.method = "pa", omega.conf.type = "boot", omega.cov.samp = TRUE,
                returnSamples = FALSE, prior.samp = FALSE) {
  if (supr.warnings) {
    options(warn = - 1)
  }
  estimates <- match.arg(estimates, several.ok = T)
  sum.res <- list()
  sum.res$call <- match.call()

  data <- scale(raw.data, scale = F)

  if("glb" %in% estimates){
    control <- Rcsdp:::csdp.control(printlevel = 0)
    Rcsdp:::write.control.file(control)
  }

  if (jags){
    sum.res$bay <- jagsFun(data, n.iter, n.burnin, estimates, interval, omega.cov.samp)
    sum.res$bayes.method <- "jags"
    sum.res$omega.pa <- omega.cov.samp
  }
  else{
    sum.res$bay <- gibbsFun(data, n.iter, n.burnin, estimates, interval, omega.cov.samp, returnSamples)
    sum.res$bayes.method <- "gibbs"
    sum.res$omega.pa <- omega.cov.samp
  }

  sum.res$freq.true <- FALSE
  if(freq){
    # sum.res$freq <- freqFun(data, boot.n, boot.interval.type, estimates, interval, omega.freq.method, omega.conf.type)
    # #this was the command with the boot package
    sum.res$freq <- freqFun2(data, boot.n, estimates, interval, omega.freq.method, omega.conf.type)
    sum.res$freq.true <- TRUE
    sum.res$omega.freq.method <- omega.freq.method
    sum.res$omega.conf.type <- omega.conf.type
    if (omega.freq.method == "pa" && omega.conf.type == "alg"){
      sum.res$omega.conf.type <- "boot"
        print("algebraic confidence interval for omega not available with method PA")
    }
  }
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
