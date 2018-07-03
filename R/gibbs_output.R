#' this function calls on other functions in order to return the gibbs sampled estimates
#' and the credible intervals together with the posterior distribution objects
#' to be passed on for forther analysis

gibbsFun <- function(data, n.iter, n.burnin, estimates, interval, omega.cov.samp, returnSamples = FALSE){
  if ("alpha" %in% estimates || "l2" %in% estimates || "l6" %in% estimates || "glb" %in% estimates || omega.cov.samp){
    C <- covSamp2(data, n.iter, n.burnin)
  } else {
  	C <- NULL
  }
  if (returnSamples) {
  	res <- list(samp = list(C = C))
  } else {
  	res <- list()
  }

  if ("alpha" %in% estimates){
    res$samp$gibbs.alpha <- coda::as.mcmc(apply(C, MARGIN = 1, applyAlpha))
    int <- coda::HPDinterval(res$samp$gibbs.alpha, prob = interval)
    res$cred$low$gibbs.alpha <- int[1]
    res$cred$up$gibbs.alpha <- int[2]
    res$est$gibbs.alpha<- median(res$samp$gibbs.alpha)
  }

  if ("l2" %in% estimates){
    res$samp$gibbs.l2 <- coda::as.mcmc(apply(C, MARGIN = 1, applyL2))
    int <- coda::HPDinterval(res$samp$gibbs.l2, prob = interval)
    res$cred$low$gibbs.l2 <- int[1]
    res$cred$up$gibbs.l2 <- int[2]
    res$est$gibbs.l2<- median(res$samp$gibbs.l2)
  }

  if ("l6" %in% estimates){
    res$samp$gibbs.l6 <- coda::as.mcmc(apply(C, MARGIN = 1, applyL6))
    int <- coda::HPDinterval(res$samp$gibbs.l6, prob = interval)
    res$cred$low$gibbs.l6 <- int[1]
    res$cred$up$gibbs.l6 <- int[2]
    res$est$gibbs.l6<- median(res$samp$gibbs.l6)
  }

  if ("glb" %in% estimates){
    res$samp$gibbs.glb <- coda::as.mcmc(apply(C, MARGIN = 1, applyGlb))
    int <- coda::HPDinterval(res$samp$gibbs.glb, prob = interval)
    res$cred$low$gibbs.glb <- int[1]
    res$cred$up$gibbs.glb <- int[2]
    res$est$gibbs.glb<- median(res$samp$gibbs.glb)
  }

  # special case omega -----------------------------------------------------------------
  if ("omega" %in% estimates){
    if (omega.cov.samp){
      res$samp$gibbs.omega <- coda::as.mcmc(apply(C, MARGIN = 1, applyOmega_boot_pa))
      int <- coda::HPDinterval(res$samp$gibbs.omega, prob = interval)
      res$cred$low$gibbs.omega <- int[1]
      res$cred$up$gibbs.omega <- int[2]
      res$est$gibbs.omega <- median(res$samp$gibbs.omega)
    }
    else{
    om.samp <- omegaSampler(data, n.iter, n.burnin)
    res$samp$gibbs.omega <- coda::as.mcmc(om.samp)
    int <- coda::HPDinterval(res$samp$gibbs.omega, prob = interval)
    res$cred$low$gibb.omega <- int[1]
    res$cred$up$gibbs.omega<- int[2]
    res$est$gibbs.omega <- median(res$samp$gibbs.omega)
    }
  }

  return(res)

}
