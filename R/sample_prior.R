#' this function samples priors for the estimates and the number of indicators

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
    priora <- apply(m, MARGIN = 1, bayesrel:::applyalpha)
    out$prioralpha <- quantiles(priora[priora >= 0])
  }
  if ("lambda2" %in% estimates){
    priorl2 <- apply(m, MARGIN = 1, bayesrel:::applyl2)
    out$priorlambda2 <- quantiles(priorl2[priorl2 >= 0])
  }
  if ("lambda4" %in% estimates){
    priorl4 <- apply(m, MARGIN = 1, bayesrel:::applyl4)
    out$priorlambda4 <- quantiles(priorl4[priorl4 >= 0])
  }
  if ("lambda6" %in% estimates){
    priorl6 <- apply(m, MARGIN = 1, bayesrel:::applyl6)
    out$priorlambda6 <- quantiles(priorl6[priorl6 >= 0])
  }
  if ("glb" %in% estimates){
    control <- Rcsdp:::csdp.control(printlevel = 0)
    Rcsdp:::write.control.file(control)
    priorglb <- apply(m, MARGIN = 1, bayesrel:::applyglb)
    out$priorglb <- quantiles(priorglb[priorglb >= 0])
    unlink("param.csdp")
  }
  if ("omega" %in% estimates){
    prioromega <- apply(m, MARGIN = 1, bayesrel:::applyomega_pa)
    out$prioromega <- quantiles(prioromega[prioromega >= 0])
  }

  return(out)

}


