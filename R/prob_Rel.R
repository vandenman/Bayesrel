#' takes a mcmc posterior sample of any of the internal consistency estimates
#' and calculates any given probability of the estimate being bigger
#' or smaller than an arbitrary value
#' @export
probrel <- function(x, low.bound){
  obj <- ecdf(x)
  est.prob <- 1 - obj(low.bound)
  return(est.prob)
}

