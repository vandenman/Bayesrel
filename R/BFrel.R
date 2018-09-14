
#' gives the odds ratio of an estimate being bigger than low.bound between the prior and the posterior
#' needs a posterior and a prior sample of the estimate as input, and the critical value of the hypothesis
#'
#' @export
bfrel <- function(xprior, xpost, low.bound){
  probrel(xprior, low.bound)/probrel(xpost, low.bound)
}
