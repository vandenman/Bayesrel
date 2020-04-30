# the basic functions for calculating and bootstrapping the internal consistency estimates

#######  measures functions ##########

applyalpha <- function(M){
  p <- ncol(M)
  a <- (p/(p-1))*(1-(sum(diag((M)))/sum(M)))
  return(a)
}

applylambda2 <- function(M){
  p <- ncol(M)
  M0 <- M
  diag(M0) <- 0
  lambda2 <- (sum(M0) + sqrt(p/(p-1) * sum(M0^2))) / sum(M)
  return(lambda2)
}

applylambda4 <- function(M){
  if (ncol(M) < 15) {l4 <- MaxSplitExhaustive(M)}
  else {l4 <- quant.lambda4(M)}
  return(l4)
}


applylambda6 <- function(M){
  M <- cov2cor(M)
  smc <- 1 - (1 / diag(solve(M)))
  lambda6 <- 1 - (sum(1 - (smc)) / sum(M))
  return(lambda6)
}

# applyglb <- function(M){
#   gl <- glbOnArray(M)
#   return(gl)
# }

applyglb <- function(M){
  gl <- glb.algebraic2(M)
  return(gl)
}

applyomega_cfa_data <- function(data, interval, pairwise){
  out <- omegaFreqData(data, interval=interval, omega.int.analytic=T, pairwise=pairwise)
  om <- out$omega
  return(om)
}

applyomega_cfa_cov <- function(cv, interval, omega.int.analytic, pairwise, n.boot){
  data <- MASS::mvrnorm(500, numeric(ncol(cv)), cv)
  out <- omegaFreqData(data, interval, omega.int.analytic, pairwise, n.boot)
  om <- out$omega
  return(om)
}

applyomega_pfa <- function(m){
  f <- princFac(m)
  l_fa <- f$loadings
  er_fa <- f$err_var
  om <- sum(l_fa)^2 / (sum(l_fa)^2 + sum(er_fa))
  if (om < 0 || om > 1 || is.na(om)) om <- NA
  return(om)
}


