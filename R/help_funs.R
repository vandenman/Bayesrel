# forces a quadratic matrix to be symmetrical
# source:
# R package fifer from Dustin fifer
# https://www.rdocumentation.org/packages/fifer
make.symmetric <- function(a, lower.tri=TRUE){
  if (lower.tri){
    ind <- upper.tri(a)
    a[ind] <- t(a)[ind]
  } else {
    ind <- lower.tri(a)
    a[ind] <- t(a)[ind]
  }
  a
}

# computes alpha analytical interval with given bounds
# source:
# Bonett, D. G., & Wright, T. A. (2015). Cronbach’s alpha reliability:
# Interval estimation, hypothesis testing, and sample size planning. Journal of Organizational Behavior,36(1), 3–15.
#
ciAlpha <- function(palpha, n, V){
  p <- ncol(V)
  z <- qnorm(1 - palpha/2)
  b <- log(n/ (n - 1))
  j <- matrix(rep(1, p), nrow = p, ncol = 1)
  a0 <- t(j) %*% V %*% j
  t1 <- sum(diag(V))
  t2 <- sum(diag(V %*% V))
  a1 <- a0 ^ 3
  a2 <- a0 * (t2 + t1^2)
  a3 <- 2 * t1 * t(j) %*% (V %*% V) %*% j
  a4 <- 2 * p^2 / (a1 * (p - 1)^2)
  r <- (p/ (p-1)) * (1 - t1 / a0)
  var <- a4 * (a2 - a3) / (n - 3)
  ll <- 1 - exp(log(1 - r) - b + z * sqrt(var / (1 - r)^2))
  ul <- 1 - exp(log(1 - r) - b - z * sqrt(var / (1 - r)^2))
  out <- c(ll, ul)
  return(out)
}

# does quantile averaging and returns 2000 datapoints
quantiles <- function(samp){
  q <- quantile(samp, probs = seq(0, 1, length.out = 2e3))
  return(q)
}


se <- function(x) {
  b <- length(x)
  se <- sqrt(1/(b-1) * sum((x - mean(x))^2))
  se
}


# create lavaan cfa one factor model file from data

lavOneFile <- function(data){
  p <- ncol(data)
  v <- 0
  for(i in 1:p){
    v[i] <- paste0("x", i)
  }
  v <- paste0(v, collapse = "+")
  mod <- paste0("g=~", v) # dynamic lavaan model file

  # column names specify
  names <- 0
  for(i in 1:p){
    names[i] <- paste0("x",i)
  }
  return(list(names = names, model = mod))
}


# create multivariate sample with mu and sigma given
# source: MASS package
# Venables, W. N. & Ripley, B. D. (2002) Modern Applied Statistics with S. Fourth Edition. Springer, New York. ISBN 0-387-95457-0
#
mvrnorm2 <- function(n = 1, mu, Sigma, tol = 1e-06, empirical = TRUE){
  p <- length(mu)
  if (!all(dim(Sigma) == c(p, p)))
    stop("incompatible arguments")
  eS <- eigen(Sigma, symmetric = TRUE)
  ev <- eS$values
  if (!all(ev >= -tol * abs(ev[1L])))
    stop("'Sigma' is not positive definite")
  X <- matrix(rnorm(p * n), n)
  if (empirical) {
    X <- scale(X, TRUE, FALSE)
    X <- X %*% svd(X, nu = 0)$v
    X <- scale(X, FALSE, TRUE)
  }
  X <- drop(mu) + eS$vectors %*% diag(sqrt(pmax(ev, 0)), p) %*%
    t(X)
  nm <- names(mu)
  if (is.null(nm) && !is.null(dn <- dimnames(Sigma)))
    nm <- dn[[1L]]
  dimnames(X) <- list(nm, NULL)
  if (n == 1)
    drop(X)
  else t(X)
}

