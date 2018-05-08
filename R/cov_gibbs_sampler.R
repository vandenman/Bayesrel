#' this function uses gibbs sampling to estimate the posterior distribution
#' of a sample's covariance matrix

covSamp <- function(data, n.iter, n.burnin){
  n <- nrow(data)
  p <- ncol(data)
  # posterior covariance matrix ---------------------------------------------------
  k0 <- 1e-10
  v0 <- p
  t <- diag(p)
  T0 <- solve(t/k0) # inverse scale matrix, prior
  mu0 <- rep(0, p) # prior means
  ym <- apply(data, 2, mean)
  # https://en.wikipedia.org/wiki/Normal-inverse-Wishart_distribution, murphy 2007
  mun <- (k0 * mu0 + n * ym) / (k0 + n)
  kn <- k0 + n
  vn <- v0 + n
  S <- 0
  for (i in 1:n){
    M <- (data[i, ] - ym) %*% t(data[i, ] - ym)
    S <- S + M
  }
  Tn <- T0 + S + (k0 * n / (k0 + n)) * (ym - mu0) %*% t(ym - mu0)
  # drawing samples from posterior:
  c.post <- array(0, c(n.iter, p, p))
  for ( i in 1:n.iter){
    c.post[i, , ] <- LaplacesDemon::rinvwishart(vn, Tn) # sample from inverse wishart
  }
  c.post <- c.post[(n.burnin + 1):n.iter, , ]

  return(c.post)
}
