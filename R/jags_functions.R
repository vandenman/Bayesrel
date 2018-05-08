#' functions for sampling using JAGS
#' separate sampling for the covariance matrix and the cfa model

modelfile <- function() {

  # priors
  V ~ dwish(R, p)
  for (i in 1:p) {
    mu[i] ~ dnorm(0, 1E-4)
  }

  # likelihood
  for (i in 1:n) {
    x[i, ] ~ dmnorm(mu, V)
  }

  # raw cronbachs alpha -- adapted from psych::alpha
  C <- inverse(V)
}

bayesianRel = function(x, ...) {

  data <- list(x = x, R = diag(ncol(x)), n = nrow(x), p = ncol(x))
  params <- c("C")
  capture.output(
  res <- R2jags::jags(data, parameters.to.save = params,
                model.file = modelfile, ...)
  )
  return(list(jags.res = res))

}


##### special case Omega #################### okamoto 2013
modelOmega_h <- function() {
  # priors
  for (i in 1:p) {
    lambda[i] ~ dunif(0, 10)
    psi[i] ~ dgamma(.001, .001)
  }
  for(i in 1:n) {
    f[i] ~ dnorm(0, 1)
  }

  # likelihood
  for(i in 1:n) {
    for(j in 1:p) {
      x[i,j] ~ dnorm(lambda[j]*f[i], psi[j])
    }
  }

  # Omega
  tmp <- pow(sum(lambda), 2)
  omega_h <- tmp  / (tmp + (sum(1/psi)))
}

bayesianOmega_h <- function(x, ...) {

  data <- list(x = x,
              n = nrow(x), p = ncol(x))
  params <- c("lambda", "omega_h",  "psi", "f")
  capture.output(
  res <- R2jags::jags(data, parameters.to.save = params,
              model.file = modelOmega_h, ...)
  )
  return(list(jags.res = res,
              omega_h = mean(res$BUGSoutput$sims.list$omega_h)))
}




