# gives freq omega, and loadings and errors
#

omegaFreqData <- function(data){
  p <- ncol(data)
  file <- lavOneFile(data)
  colnames(data) <- file$names
  mod <- file$model
  fit <- try(lavaan::cfa(mod, data, std.lv = T), silent = TRUE)
  params <- try(lavaan::parameterestimates(fit), silent = TRUE)
  if ("try-error" %in% class(params)) {
    load <- resid <- omega <- om_low <- om_up <- fit_tmp <- indic <- NA
  } else {
    load <- params$est[1:p]
    resid <- params$est[(p+2) : (p*2+1)]
    omega <- omegaBasic(load, resid)

    load_low <- params$ci.lower[1:p]
    resid_low <- params$ci.lower[(p+2) : (p*2+1)]
    om_low <- omegaBasic(load_low, resid_low)
    load_up <- params$ci.upper[1:p]
    resid_up <- params$ci.upper[(p+2) : (p*2+1)]
    om_up <- omegaBasic(load_up, resid_up)

    fit_tmp <- lavaan::fitMeasures(fit)
    indic <- c(fit_tmp["chisq"], fit_tmp["df"], fit_tmp["pvalue"],
               fit_tmp["rmsea"], fit_tmp["rmsea.ci.lower"], fit_tmp["rmsea.ci.upper"],
               fit_tmp["srmr"])
  }

  return(list(omega = omega, loadings = load, errors = resid,
              omega_lower = om_low, omega_upper = om_up, indices = indic))
}

