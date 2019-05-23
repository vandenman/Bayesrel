
#' graphical posterior predictive check for the 1-factor omega model, based on covariance matrix eigenvalues
#'
#' @description
#' gives posterior predictive check for the 1-factor model:
#' comparison between model implied covariance matrix and sample covariance matrix
#' also displays frequentist fit indices
#'
#' @param x A strel output object (list)
#'
#' @examples fit.omega(strel(cavalini[1:200, ], "omega"))
#'
#' @export
fit.omega <- function(x){
  if (!("omega" %in% x$estimates)) {return("please run the analysis with omega as an estimate")}

  if (!is.null(x$freq$fit.omega)){
    print(x$freq$fit.omega)}

  sigma <- cov(x$data)
  lambda <- x$bay$loadings
  psi <- x$bay$resid.var
  cimpl <- lambda %*% t(lambda) + diag(psi)
  ymax <- max(eigen(cimpl)$values, eigen(sigma)$values) * 1.3
  plot(eigen(sigma)$values, axes = F, ylim = c(0, ymax), ylab = "Eigenvalue - Size", xlab = "Eigenvalue - No.")
  axis(side = 1, at = seq(1:ncol(sigma)))
  axis(side = 2)
  title(main = "Posterior Predictive Check for Omega 1-Factor-Model")

  for (i in 1:1e3) {
    dtmp <- MASS::mvrnorm(nrow(x$data), rep(0, ncol(sigma)), cimpl)
    lines(eigen(cov(dtmp))$values, type = "l", col = "gray")

  }
  lines(eigen(sigma)$values, col = "black", type = "p")
  lines(eigen(sigma)$values, col = "black", lwd = 2)
  legend(ncol(sigma)/3, ymax*(2/3), legend = c("Dataset Covariance Matrix", "Simulated Data from Model Implied Covariance Matrix"),
         col=c("black", "gray"), lwd = c(2, 2), box.lwd = 0)
}
