#' probability of estimate being bigger than threshold
#' @description
#' takes a mcmc posterior sample of any of the single test reliability estimates
#' and calculates any given probability of the estimate being bigger
#' or smaller than an arbitrary value
#'
#' @param x A strel output object (list)
#' @param estimate A character string indicating what estimate to plot from the strel output object
#' @param low.bound A number for the threshold to be tested against
#'
#' @examples
#' p_strel(strel(asrm, "lambda2"), "lambda2", .80)
#' @export


p_strel <- function(x, estimate, low.bound){
  posi <- grep(estimate, x$estimates, ignore.case = T)
  samp <- unlist(x$Bayes$samp[posi])
  obj <- ecdf(samp)
  est_prob <- 1 - obj(low.bound)
  return(est_prob)
}

# probic <- function(x, low.bound){
#   obj <- ecdf(x)
#   est_prob <- 1 - obj(low.bound)
#   return(est_prob)
# }
