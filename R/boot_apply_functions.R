#' the basic functions for calculating and bootstrapping the internal consistency estimates


#######  measures functions ##########

applyalpha <- function(M){
  p <- ncol(M)
  a <- (p/(p-1))*(1-(sum(diag((M)))/sum(M)))
  return(a)
}

applyl2 <- function(M){
  p <- ncol(M)
  M0 <- M
  diag(M0) <- 0
  l2 <- (sum(M0) + sqrt(p/(p-1) * sum(M0^2))) / sum(M)
  return(l2)
}

applyl4 <- function(M){
  if (ncolM < 15) l4 <- MaxSplitExhaustive(M)
  else l4 <- quant.lambda4(M)
  return(l4)
}


applyl6 <- function(M){
  l6 <- 1 - (sum(1 - (1 - (1 / diag(solve(M))))) / sum(M))
  return(l6)
}

applyglb <- function(M){
  gl <- glb.algebraic2(M)$glb
  return(gl)
}

applyomega_cfa <- function(data){
  lav.mod.file <- lavOneFile(data)
  colnames(data) <- lav.mod.file$names
  p <- ncol(data)
  fitdata <- lavaan::cfa(model = lav.mod.file$model, data = data, std.lv = T)
  params <- lavaan::parameterestimates(fitdata)
  om <- sum(params$est[1:p])^2 / (sum(params$est[1:p])^2 + sum(params$est[(p+1):(p+p)]))
  if (om < 0 || om > 1 || is.na(om)) om <- NA
  return(om)
}

applyomega_pa <- function(m){
  f <- princFac(m)
  l.fa <- f$loadings
  er.fa <- f$err.var
  om <- sum(l.fa)^2 / (sum(l.fa)^2 + sum(er.fa))
  if (om < 0 || om > 1 || is.na(om)) om <- NA
  return(om)
}

applyomega_alg <- function(data, interval){
  lav.mod.file <- lavOneFile(data)
  colnames(data) <- lav.mod.file$names
  p <- ncol(data)
  fitdata <- lavaan::cfa(model = lav.mod.file$model, data = data, std.lv = T)
  params <- lavaan::parameterestimates(fitdata, ci = TRUE, level = interval)
  om <- sum(params$est[1:p])^2 / (sum(params$est[1:p])^2 + sum(params$est[(p+1):(p+p)]))
  if (om < 0 || om > 1 || is.na(om)) om <- NA
  om.low <- sum(params$ci.lower[1:p])^2 / (sum(params$ci.lower[1:p])^2 + sum(params$ci.lower[(p+1):(p+p)]))
  if (om.low < 0 || om.low > 1 || is.na(om.low)) om.low <- NA
  om.up <- sum(params$ci.upper[1:p])^2 / (sum(params$ci.upper[1:p])^2 + sum(params$ci.upper[(p+1):(p+p)]))
  if (om.up < 0 || om.up > 1 || is.na(om.up)) om.up <- NA
  oms <- c(om, om.low, om.up)
  return(oms)
}

omegaBasic <- function(l, e){
  o <- sum(l)^2 / (sum(l)^2 + sum(e))
  return(o)
}


