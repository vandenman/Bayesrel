#' functions to return the jags sampled estimates of internal concsistency
#' together with credible intervals as well as the sampled posterior distribution

bayesFun <- function(data, n.iter, n.burnin, estimates, interval, omega.cov.samp, return.cov.samples){

  bay.cov <- bayesianRel(data, inits = NULL, n.iter, n.burnin,
                       n.thin = 1, n.chains = 1, DIC = FALSE)
# large object that contains the posterior cov matrices:
  jC <- bay.cov$bayes.res$BUGSoutput$sims.list$C

  res <- list()
  if (return.cov.samples){
    res$samp$C <- jC
  }
  if ("alpha" %in% estimates){
    res$samp$bayes.alpha <- coda::as.mcmc(apply(jC, MARGIN = 1, applyAlpha))
    int <- coda::HPDinterval(res$samp$bayes.alpha, prob = interval)
    res$cred$low$bayes.alpha <- int[1]
    res$cred$up$bayes.alpha <- int[2]
    res$est$bayes.alpha<- median(res$samp$bayes.alpha)
  }

  if ("l2" %in% estimates){
    res$samp$bayes.l2 <- coda::as.mcmc(apply(jC, MARGIN = 1, applyL2))
    int <- coda::HPDinterval(res$samp$bayes.l2, prob = interval)
    res$cred$low$bayes.l2 <- int[1]
    res$cred$up$bayes.l2 <- int[2]
    res$est$bayes.l2 <- median(res$samp$bayes.l2)
  }

  if ("l4" %in% estimates){
    res$samp$bayes.l4 <- coda::as.mcmc(apply(jC, MARGIN = 1, applyL4))
    int <- coda::HPDinterval(res$samp$bayes.l4, prob = interval)
    res$cred$low$bayes.l4 <- int[1]
    res$cred$up$bayes.l4 <- int[2]
    res$est$bayes.l4 <- median(res$samp$bayes.l4)
  }

  if ("l6" %in% estimates){
    res$samp$bayes.l6 <- coda::as.mcmc(apply(jC, MARGIN = 1, applyL6))
    int <- coda::HPDinterval(res$samp$bayes.l6, prob = interval)
    res$cred$low$bayes.l6 <- int[1]
    res$cred$up$bayes.l6 <- int[2]
    res$est$bayes.l6 <- median(res$samp$bayes.l6)
  }

  if ("glb" %in% estimates){
    res$samp$bayes.glb <- coda::as.mcmc(apply(jC, MARGIN = 1, applyGlb))
    int <- coda::HPDinterval(res$samp$bayes.glb, prob = interval)
    res$cred$low$bayes.glb <- int[1]
    res$cred$up$bayes.glb <- int[2]
    res$est$bayes.glb<- median(res$samp$bayes.glb)
  }

  # special case omega ----------------------------------------------------------------
  if ("omega" %in% estimates){
    if (omega.cov.samp){
      res$samp$bayes.omega <- coda::as.mcmc(apply(jC, MARGIN = 1, applyOmega_boot_pa))
      int <- coda::HPDinterval(res$samp$bayes.omega, prob = interval)
      res$cred$low$bayes.omega <- int[1]
      res$cred$up$bayes.omega <- int[2]
      res$est$bayes.omega <- median(res$samp$bayes.omega)
    }
    else{
    bay.omega <- bayesianOmega_h(data, inits = NULL, n.iter, n.burnin,
                                 n.thin = 1, n.chains = 1, DIC = FALSE)
    res$samp$bayes.omega <- coda::as.mcmc(bay.omega$bayes.res$BUGSoutput$sims.list$omega_h)
    int <- coda::HPDinterval(res$samp$bayes.omega, prob = interval)
    res$cred$low$bayes.omega <- int[1]
    res$cred$up$bayes.omega <- int[2]
    res$est$bayes.omega <- median(res$samp$bayes.omega)
    }
  }
  return(res)
}
