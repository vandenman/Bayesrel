#' create lavaan cfa one factor model file from data

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
