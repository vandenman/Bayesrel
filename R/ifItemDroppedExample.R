p <- 40
nobs <- 10

# some covariance matrix
aa <- matrix(runif(p^2), p, p)
sigma <- t(aa) %*% aa

x <- mvtnorm::rmvnorm(n=500, sigma=sigma)

# covariance are equal
cov0 <- var(x)[1:(p-1), 1:(p-1)]
cov1 <- var(x[, 1:(p-1)])
all.equal(cov0, cov1)

# I made a small adjustment to store
# C in res$samp$C in gibbsFun. Perhaps it is nice
#  to allow this with a boolean that is FALSE by default
# (e.g., storeCovSamples = FALSE).

n.iter <- 1e4

# the full thing
e0 <- bayesrel::brel(x, freq = FALSE, estimates = "alpha", n.iter = n.iter, returnSamples = TRUE)
# extract the covariance samples but drop one variable
C0 <- e0$bay$samp$C[, 1:(p-1), 1:(p-1)] # estimates
# calculate alpha on reduced covariance samples
e0a <- coda::as.mcmc(apply(C0, MARGIN = 1, bayesrel:::applyAlpha))
# directly calculate alpha on reduced dataset
e1 <- bayesrel::brel(x[, 1:(p-1)], freq = FALSE, estimates = "alpha", n.iter = n.iter, returnSamples = TRUE)
e1a <- e1$bay$samp$gibbs.alpha

# density plot
par(bty = "n", las = 1)
# plot(density(e0a))
# lines(density(e1a), col = 2)

# very convincing qq-plot
step <- 1e-3
probs <- seq(step, 1 - step, step)
q0 <- quantile(e0a, probs = probs)
q1 <- quantile(e1a, probs = probs)
plot(q0, q1)
abline(0, 1)
