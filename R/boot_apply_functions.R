# the basic functions for calculating and bootstrapping the internal consistency estimates

#######  measures functions ##########

applyalpha <- function(M, callback = function(){}){
  p <- ncol(M)
  a <- (p/(p-1))*(1-(sum(diag((M)))/sum(M)))
  callback()
  return(a)
}

applylambda2 <- function(M, callback = function(){}){
  p <- ncol(M)
  M0 <- M
  diag(M0) <- 0
  lambda2 <- (sum(M0) + sqrt(p/(p-1) * sum(M0^2))) / sum(M)
  callback()
  return(lambda2)
}

applylambda4 <- function(M, callback = function(){}){
  if (ncol(M) < 15) {l4 <- MaxSplitExhaustive(M)}
  else {l4 <- quant.lambda4(M)}
  callback()
  return(l4)
}


applylambda6 <- function(M, callback = function(){}){
  M <- cov2cor(M)
  smc <- try(1 - (1 / diag(solve(M))), silent = T)
  if (class(smc) == "try-error") {
    lambda6 <- NA; warning("singular bootstrapped covariance matrices encountered")
  } else {
    lambda6 <- 1 - (sum(1 - (smc)) / sum(M))
  }
  callback()
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

applyomega_cfa_data <- function(data, interval, pairwise, callback = function(){}){
  out <- omegaFreqData(data, interval=interval, omega.int.analytic=T, pairwise=pairwise)
  om <- out$omega
  callback()
  return(om)
}

applyomega_cfa_cov <- function(cv, interval, omega.int.analytic, pairwise, n.boot){
  data <- MASS::mvrnorm(500, numeric(ncol(cv)), cv)
  out <- omegaFreqData(data, interval, omega.int.analytic, pairwise, n.boot)
  om <- out$omega
  return(om)
}

applyomega_pfa <- function(m, callback = function(){}){
  f <- princFac(m)
  l_fa <- f$loadings
  er_fa <- f$err_var
  om <- sum(l_fa)^2 / (sum(l_fa)^2 + sum(er_fa))
  if (om < 0 || om > 1 || is.na(om)) om <- NA
  callback()
  return(om)
}




