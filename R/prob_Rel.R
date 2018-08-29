#' this function takes a mcmc posterior sample of any of the internal consistency estimates
#' and calculates any given probability of the estimate being bigger
#' or smaller than an arbitrary value
#' @export
probRel <- function(obj, low.bound){

  if (0 < low.bound && low.bound < 1){
    obj <- ecdf(obj)
    est.prob <- 1 - obj(low.bound)
    return(est.prob)
  }
  else{
    return("please supply a value between 0 and 1")
  }

}

