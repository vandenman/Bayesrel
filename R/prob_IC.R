#' takes a mcmc posterior sample of any of the internal consistency estimates
#' and calculates any given probability of the estimate being bigger
#' or smaller than an arbitrary value
#' @export

probic <- function(x, low.bound){
  obj <- ecdf(x)
  est.prob <- 1 - obj(low.bound)
  return(est.prob)
}

# probrel <- function(x, estimate, low.bound){
#   posi <- grep(estimate, x$estimates, ignore.case = T)
#   samp <- coda::as.mcmc(unlist(x$bay$samp[posi]))
#   obj <- ecdf(samp)
#   est.prob <- 1 - obj(low.bound)
#   return(est.prob)
# }
