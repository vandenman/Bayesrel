# gives freq omega, and loadings and errors
#

omegaFreqData <- function(data){
  p <- ncol(data)
  file <- lavOneFile(data)
  colnames(data) <- file$names
  mod <- file$model
  fit <- try(lavaan::cfa(mod, data), silent = TRUE)
  params <- try(lavaan::standardizedsolution(fit), silent = TRUE)
  if ("try-error" %in% class(params)) {
    load <- resid <- omega <- om.low <- om.up <- fit.tmp <- indic <- NA
  } else {
    load <- params$est.std[1:p]
    resid <- params$est.std[(1+p) : (p*2)]
    omega <- omegaBasic(load, resid)

    load.low <- params$ci.lower[1:p]
    resid.low <- params$ci.lower[(1+p) : (p*2)]
    om.low <- omegaBasic(load.low, resid.low)
    load.up <- params$ci.upper[1:p]
    resid.up <- params$ci.upper[(1+p) : (p*2)]
    om.up <- omegaBasic(load.up, resid.up)

    fit.tmp <- lavaan::fitMeasures(fit)
    indic <- c(fit.tmp["chisq"], fit.tmp["df"], fit.tmp["pvalue"],
               fit.tmp["rmsea"], fit.tmp["rmsea.ci.lower"], fit.tmp["rmsea.ci.upper"],
               fit.tmp["srmr"])
  }

  return(list(omega = omega, loadings = load, errors = resid,
              omega.lower = om.low, omega.upper = om.up, indices = indic))
}

