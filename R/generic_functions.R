
#'@export
print.strel <- function(x, ...){
  if (!is.null(x$freq)){
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
summary.strel <- function(object, ...){

  out_matrix <- list()
  if (!is.null(object$freq)){
    out_matrix$est <- rbind(as.data.frame(as.matrix(object$bay$est)),
                            as.data.frame(as.matrix(object$freq$est)))
    out_matrix$int$low <- rbind(as.data.frame(as.matrix(object$bay$cred$low)),
                                as.data.frame(as.matrix(object$freq$conf$low)))
    out_matrix$int$up <- rbind(as.data.frame(as.matrix(object$bay$cred$up)),
                               as.data.frame(as.matrix(object$freq$conf$up)))
    out_matrix$omega_freq_method <- object$omega_freq_method
  } else{
    out_matrix$est <- as.data.frame(as.matrix(object$bay$est))
    out_matrix$int$low <- as.data.frame(as.matrix(object$bay$cred$low))
    out_matrix$int$up <- as.data.frame(as.matrix(object$bay$cred$up))
  }
  out_matrix$call <- object$call
  out_matrix$n.iter <- object$n.iter
  out_matrix$n.burnin <- object$n.burnin
  out_matrix$interval <- object$interval
  out_matrix$estimates <- object$estimates
  out_matrix$ifitem$bay_tab <- object$bay$ifitem$est
  out_matrix$ifitem$freq_tab <- object$freq$ifitem

  class(out_matrix) <- "summary.strel"
  out_matrix
}

#'@export
print.summary.strel <- function(x, ...){
  n_row <- length(x$est$V1)
  mat <- data.frame(matrix(0, n_row, 0))
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

  if (length(grep("freq", x$est)) > 0){
    if ("omega" %in% x$estimates){
      cat("frequentist omega method is:")
      print.default(x$omega_freq_method)
      cat("omega confidence interval is estimated with:")
      if (x$omega_freq_method == "pfa") {print.default("bootstrap")}
      if (x$omega_freq_method == "cfa") {print.default("maximum likelihood z-value")}
    }
    # if (!is.null(x$fit.indices)){
    #   options(scipen = 999)
    #   cat("\nFrequentist fit of 1-factor-model for omega is:\n")
    #   print.default(as.matrix(x$fit.indices))
    # }
  }


  if (!is.null(x$ifitem$bay_tab)){
    n_row <- length(unlist(x$ifitem$bay_tab[1])) + 1
    n_col <- length(x$ifitem$bay_tab)
    mat_ifitem_bay <- data.frame(matrix(0, n_row, n_col))
    mat_ifitem_bay[1, ] <- as.vector(unlist(x$est))[1:n_col]
    for (i in 1:n_col){
      mat_ifitem_bay[2:n_row, i] <- as.vector(unlist(x$ifitem$bay_tab[i]))
    }
    colnames(mat_ifitem_bay) <- x$estimates
    names <- NULL
    for(i in 1:(n_row-1)){
      names[i] <- paste0("x", i)
    }
    row.names(mat_ifitem_bay) <- c("original", names)

    if (length(grep("freq", x$est)) > 0){
      mat_ifitem_freq <- data.frame(matrix(0, n_row, n_col))
      mat_ifitem_freq[1, ] <- as.vector(unlist(x$est)[(n_col+1):(n_col*2)])
      for (i in 1:n_col){
        mat_ifitem_freq[2:n_row, i] <- as.vector(unlist(x$ifitem$freq_tab[i]))
      }
      colnames(mat_ifitem_freq) <- x$estimates
      row.names(mat_ifitem_freq) <- c("original", names)
    }

    cat("\n")
    cat("Bayesian coefficient if item dropped: \n")
    print(mat_ifitem_bay)
    if (length(grep("freq", x$est)) > 0){
      cat("\n")
      cat("Frequentist coefficient if item dropped: \n")
      print(mat_ifitem_freq)
    }
  }

}

