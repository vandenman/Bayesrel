#' functions to return the jags sampled estimates of internal concsistency
#' together with credible intervals as well as the sampled posterior distribution

jagsFun <- function(data, n.iter, n.burnin, estimates, interval){

  bay.cov <- bayesianRel(data, inits = NULL, n.iter, n.burnin,
                         n.thin = 1, n.chains = 1, DIC = FALSE)
  # large object that contains all the cov matrices:
  jC <- bay.cov$jags.res$BUGSoutput$sims.list$C
  res <- list()

  if ("alpha" %in% estimates){
    res$samp$jags.alpha <- coda::as.mcmc(apply(jC, MARGIN = 1, applyAlpha))
    int <- coda::HPDinterval(res$samp$jags.alpha, prob = interval)
    res$cred$low$jags.alpha <- int[1]
    res$cred$up$jags.alpha <- int[2]
    res$est$jags.alpha<- median(res$samp$jags.alpha)
  }

  if ("l2" %in% estimates){
    res$samp$jags.l2 <- coda::as.mcmc(apply(jC, MARGIN = 1, applyL2))
    int <- coda::HPDinterval(res$samp$jags.l2, prob = interval)
    res$cred$low$jags.l2 <- int[1]
    res$cred$up$jags.l2 <- int[2]
    res$est$jags.l2 <- median(res$samp$jags.l2)
  }

  if ("l6" %in% estimates){
    res$samp$jags.l6 <- coda::as.mcmc(apply(jC, MARGIN = 1, applyL6))
    int <- coda::HPDinterval(res$samp$jags.l6, prob = interval)
    res$cred$low$jags.l6 <- int[1]
    res$cred$up$jags.l6 <- int[2]
    res$est$jags.l6 <- median(res$samp$jags.l6)
  }

  if ("glb" %in% estimates){
    res$samp$jags.glb <- coda::as.mcmc(apply(jC, MARGIN = 1, applyGlb))
    int <- coda::HPDinterval(res$samp$jags.glb, prob = interval)
    res$cred$low$jags.glb <- int[1]
    res$cred$up$jags.glb <- int[2]
    res$est$jags.glb<- median(res$samp$jags.glb)
  }

  # special case omega ----------------------------------------------------------------
  if ("omega" %in% estimates){
    bay.omega <- bayesianOmega_h(data, inits = NULL, n.iter, n.burnin,
                                 n.thin = 1, n.chains = 1, DIC = FALSE)
    res$samp$jags.omega <- coda::as.mcmc(bay.omega$jags.res$BUGSoutput$sims.list$omega_h)
    int <- coda::HPDinterval(res$samp$jags.omega, prob = interval)
    res$cred$low$jags.omega <- int[1]
    res$cred$up$jags.omega <- int[2]
    res$est$jags.omega <- median(res$samp$jags.omega)
  }
  return(res)
}
