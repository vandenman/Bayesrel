#' the basic functions for calculating and bootstrapping the internal consistency estimates

#########   boot functions ##########

bootAlpha <- function(data, indices){
  a <- applyAlpha(cov(data[indices, ]))
  return(a)
}

#freq conf intervall with bootstrapping
bootL2 <- function(data, indices){
  l2 <- applyL2(cov(data[indices, ]))
  return(l2)
}

#freq conf intervall with bootstrapping
bootL6 <- function(data, indices){
  l6 <- applyL6(cov(data[indices, ]))
  return(l6)
}

#freq conf intervall with bootstrapping
bootL4 <- function(data, indices){
  l4 <- applyL4(cov(data[indices, ]))
  return(l4)
}

# define function for bootstrapping
bootGlb <- function(data, indices){
  gl <- glb.algebraic2(cov(data[indices, ]))$glb
  return(gl)
}

bootOmega_cfa <- function(data, indices){
  om <- applyOmega_boot_cfa(data[indices, ])
  return(om)
}
bootOmega_pa <- function(data, indices){
  om <- applyOmega_boot_pa(cov(data[indices, ]))
  return(om)
}

#######  measures functions ##########

applyAlpha <- function(M){
  p <- ncol(M)
  a <- (p/(p-1))*(1-(sum(diag((M)))/sum(M)))
  return(a)
}

applyL2 <- function(M){
  p <- ncol(M)
  M0 <- M
  diag(M0) <- 0
  l2 <- (sum(M0) + sqrt(p/(p-1) * sum(M0^2))) / sum(M)
  return(l2)
}

applyL4 <- function(M){
  l4 <- MaxSplitHalfHad12(M)
  return(l4)
}

applyL6 <- function(M){
  l6 <- 1 - (sum(1 - (1 - (1 / diag(solve(M))))) / sum(M))
  return(l6)
}

applyGlb <- function(M){
  gl <- glb.algebraic2(M)$glb
  return(gl)
}

applyOmega_boot_cfa <- function(data){
  lav.mod.file <- lavOneFile(data)
  colnames(data) <- lav.mod.file$names
  p <- ncol(data)
  fitdata <- lavaan::cfa(model = lav.mod.file$model, data = data, std.lv = T)
  params <- lavaan::parameterestimates(fitdata)
  om <- sum(params$est[1:p])^2 / (sum(params$est[1:p])^2 + sum(params$est[(p+1):(p+p)]))
  if (om < 0 || om > 1 || is.na(om)) om <- NA
  return(om)
}

applyOmega_boot_pa <- function(m){
  f <- princFac(m)
  l.fa <- f$loadings
  er.fa <- f$err.var
  om <- sum(l.fa)^2 / (sum(l.fa)^2 + sum(er.fa))
  if (om < 0 || om > 1 || is.na(om)) om <- NA
  return(om)
}

applyOmega_alg <- function(data, interval){
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

omegaBase <- function(l, e){
  o <- sum(l)^2 / (sum(l)^2 + sum(e))
  return(o)
}

quantiles <- function(samp){
  q <- quantile(samp, probs = seq(0, 1, length.out = 200))
  return(q)
}
