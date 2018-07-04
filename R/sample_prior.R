

priorSamp <- function(p, estimates, n.samp = 2e3){

  v0 <- p
  k0 <- 1e-10
  t <- diag(p)
  T0 <- solve(t/k0)
  m <- array(0, c(n.samp, p, p))
  out <- list()
  for (i in 1:n.samp){
    m[i, , ] <- LaplacesDemon::rinvwishart(v0, T0)
  }

  if ("alpha" %in% estimates){
    priora <- apply(m, MARGIN = 1, bayesrel:::applyAlpha)
    out$priorAlpha <- quantiles(priora[priora >= 0])
  }
  if ("l2" %in% estimates){
    priorl2 <- apply(m, MARGIN = 1, bayesrel:::applyL2)
    out$priorLambda2 <- quantiles(priorl2[priorl2 >= 0])
  }
  if ("l4" %in% estimates){
    priorl4 <- apply(m, MARGIN = 1, bayesrel:::applyL4)
    out$priorLambda4 <- quantiles(priorl4[priorl4 >= 0])
  }
  if ("l6" %in% estimates){
    priorl6 <- apply(m, MARGIN = 1, bayesrel:::applyL6)
    out$priorLambda6 <- quantiles(priorl6[priorl6 >= 0])
  }
  if ("glb" %in% estimates){
    control <- Rcsdp:::csdp.control(printlevel = 0)
    Rcsdp:::write.control.file(control)
    priorglb <- apply(m, MARGIN = 1, bayesrel:::applyGlb)
    out$priorGlb <- quantiles(priorglb[priorglb >= 0])
    unlink("param.csdp")
  }
  if ("omega" %in% estimates){
    prioromega <- apply(m, MARGIN = 1, bayesrel:::applyOmega_boot_pa)
    out$priorOmega <- quantiles(priorOmega[priorOmega >= 0])
  }

  return(out)

}


