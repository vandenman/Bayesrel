
#'@export
print.bayesrel <- function(x, ...){
  if (x$freq.true){
    est <- cbind(as.data.frame(as.matrix(x$bay$est)),
                 as.data.frame(as.matrix(x$freq$est)))
    colnames(est) <- c("bayes", "frequentist")
  }
  else{
    est <- as.data.frame(as.matrix(x$bay$est))
    colnames(est) <- "bayes"
  }
  row.names(est) <- x$estimates
  cat("Call: \n")
  print.default(x$call)
  cat("\n")
  cat("Point Estimates of Internal Consistency Measures: \n")
  cat("\n")
  print(est)

}

#'@export
summary.bayesrel <- function(x, ...){

  out.matrix <- list()
  if (x$freq.true){
    out.matrix$est <- rbind(as.data.frame(as.matrix(x$bay$est)),
                            as.data.frame(as.matrix(x$freq$est)))
    out.matrix$int$low <- rbind(as.data.frame(as.matrix(x$bay$cred$low)),
                                as.data.frame(as.matrix(x$freq$conf$low)))
    out.matrix$int$up <- rbind(as.data.frame(as.matrix(x$bay$cred$up)),
                               as.data.frame(as.matrix(x$freq$conf$up)))
    out.matrix$omega.freq.method <- x$omega.freq.method
    out.matrix$omega.conf.int.type <- x$omega.conf.int.type
    out.matrix$alpha.int.analytic <- x$alpha.int.analytic
  }
  else{
    out.matrix$est <- as.data.frame(as.matrix(x$bay$est))
    out.matrix$int$low <- as.data.frame(as.matrix(x$bay$cred$low))
    out.matrix$int$up <- as.data.frame(as.matrix(x$bay$cred$up))
  }
  out.matrix$call <- x$call
  out.matrix$n.iter <- x$n.iter
  out.matrix$n.burnin <- x$n.burnin
  out.matrix$boot.interval.type <- x$boot.interval.type
  out.matrix$interval <- x$interval
  out.matrix$estimates <- x$estimates
  out.matrix$omega.pa <- x$omega.pa
  out.matrix$freq.true <- x$freq.true
  out.matrix$ifitem$bool <- x$item.dropped
  if (x$item.dropped){
    out.matrix$ifitem$bay.tab <- x$bay$ifitem$est
    if (x$freq.true) {
      out.matrix$ifitem$freq.tab <- x$freq$ifitem
      }
  }
  class(out.matrix) <- "summary.bayesrel"
  out.matrix
}

#'@export
print.summary.bayesrel <- function(x, ...){
  n.row <- length(x$est$V1)
  mat <- data.frame(matrix(0, n.row, 0))
  mat[, 1] <- x$est
  mat[, 2] <- '   '
  mat[, 3] <- x$int$low
  mat[, 4] <- x$int$up
  row.names(mat) <- row.names(x$est)
  colnames(mat) <- c("estimate", '', "interval.low", "interval.up")

  cat("Call: \n")
  print.default(x$call)
  cat("Results: \n")
  print(mat, right = F)
  cat("\n")
  cat("iterations: ")
  print.default(x$n.iter)
  cat("burnin: ")
  print.default(x$n.burnin)
  cat("uncertainty interval:")
  print.default(x$interval)
  if (x$omega.pa){
    cat("bayesian omega is calculated with posterior sampling from covariance matrices\n")
  }
  if (x$freq.true){
    cat("confidence intervals are estimated with bootstrapping interval type:")
    print.default(x$boot.interval.type)
    if ("omega" %in% x$estimates){
      cat("frequentist omega method is:")
      print.default(x$omega.freq.method)
      cat("omega confidence interval is estimated with:")
      print.default(x$omega.conf.int.type)
    }
    if (x$alpha.int.analytic){
      cat("alpha confidence interval is analytically computed")
    }
  }

  if (x$ifitem$bool){
    n.row <- length(unlist(x$ifitem$bay.tab[1])) + 1
    n.col <- length(x$ifitem$bay.tab)
    mat.ifitem.bay <- data.frame(matrix(0, n.row, n.col))
    mat.ifitem.bay[1, ] <- as.vector(unlist(x$est))[1:n.col]
    for (i in 1:n.col){
      mat.ifitem.bay[2:n.row, i] <- as.vector(unlist(x$ifitem$bay.tab[i]))
    }
    colnames(mat.ifitem.bay) <- x$estimates
    names <- NULL
    for(i in 1:(n.row-1)){
      names[i] <- paste0("x", i)
    }
    row.names(mat.ifitem.bay) <- c("original", names)

    if (x$freq.true){
      mat.ifitem.freq <- data.frame(matrix(0, n.row, n.col))
      mat.ifitem.freq[1, ] <- as.vector(unlist(x$est)[(n.col+1):(n.col*2)])
      for (i in 1:n.col){
        mat.ifitem.freq[2:n.row, i] <- as.vector(unlist(x$ifitem$freq.tab[i]))
      }
      colnames(mat.ifitem.freq) <- x$estimates
      row.names(mat.ifitem.freq) <- c("original", names)
    }

    cat("\n")
    cat("Bayesian coefficient if item dropped: \n")
    print(mat.ifitem.bay)
    if (x$freq.true){
      cat("\n")
      cat("Frequentist coefficient if item dropped: \n")
      print(mat.ifitem.freq)
    }
  }

}

