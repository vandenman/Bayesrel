#' this function calls on other functions in order to return the fuentist estimates
#' and bootstrapped cinfidence intervals, which can be of different types

freqFun<- function(data, boot.n, boot.interval.type, estimates, interval, omega.freq.method, omega.conf.type){
  p <- ncol(data)
  n <- nrow(data)
  res <- list()

  if ("alpha" %in% estimates){
    res$est$freq.alpha <- applyAlpha(cov(data))
    f.alpha.boot.obj <- boot::boot(data, statistic = bootAlpha, R = boot.n)
    tmp <- boot::boot.ci(f.alpha.boot.obj, conf = interval, type = boot.interval.type)[[4]]
    res$ci$low$freq.alpha <- tmp[1, 4]
    res$ci$up$freq.alpha <- tmp[1, 5]
  }
  if ("l2" %in% estimates){
    res$est$freq.l2 <- applyL2(cov(data))
    f.l2.boot.obj <- boot::boot(data, statistic = bootL2, R = boot.n)
    tmp <- boot::boot.ci(f.l2.boot.obj, conf = interval, type = boot.interval.type)[[4]]
    res$ci$low$freq.l2 <- tmp[1, 4]
    res$ci$up$freq.l2 <- tmp[1, 5]
  }

  if ("l6" %in% estimates){
    res$est$freq.l6 <- applyL6(cov(data))
    f.l6.boot.obj <- boot::boot(data, statistic = bootL6, R = boot.n)
    tmp <- boot::boot.ci(f.l6.boot.obj, conf = interval, type = boot.interval.type)[[4]]
    res$ci$low$freq.l6 <- tmp[1, 4]
    res$ci$up$freq.l6 <- tmp[1, 5]
  }

  if ("glb" %in% estimates){
    res$est$freq.glb <- applyGlb(cov(data))
    f.glb.boot.obj <- boot::boot(data, statistic = bootGlb, R = boot.n)
    tmp <- boot::boot.ci(f.glb.boot.obj, conf = interval, type = boot.interval.type)[[4]]
    res$ci$low$freq.glb <- tmp[1, 4]
    res$ci$up$freq.glb <- tmp[1, 5]
  }

  #omega --------------------------------------------------------------------------
  if ("omega" %in% estimates){
    if (omega.freq.method == "cfa"){
      res$est$freq.omega <- applyOmega_boot(data)
      if (omega.conf.type == "boot"){
        f.omega.boot.obj <- boot::boot(data, statistic = bootOmega_cfa, R = boot.n)
        tmp <- boot::boot.ci(f.omega.boot.obj, conf = interval, type = boot.interval.type)[[4]]
        if (tmp[1, 4] < 0 || tmp[1, 4] > 1 || is.na(tmp[1, 4])) {tmp[1, 4] <- NA}
        if (tmp[1, 5] < 0 || tmp[1, 5] > 1 || is.na(tmp[1, 5])) {tmp[1, 5] <- NA}
        res$ci$low$freq.omega <- tmp[1, 4]
        res$ci$up$freq.omega <- tmp[1, 5]
      }
      else{
        omega.alg <- applyOmega_alg(data, interval)
        # res$est$freq.omega.alg <- omega.alg[1]
        res$ci$low$freq.omega <- omega.alg[2]
        res$ci$up$freq.omega <- omega.alg[3]
      }
    }
    if (omega.freq.method == "pa"){
      res$est$freq.omega <- applyOmega_boot_pa(data, omega.freq.method)
      f.omega.boot.obj <- boot::boot(data, statistic = bootOmega_pa, R = boot.n)
      tmp <- boot::boot.ci(f.omega.boot.obj, conf = interval, type = boot.interval.type)[[4]]
      if (tmp[1, 4] < 0 || tmp[1, 4] > 1 || is.na(tmp[1, 4])) {tmp[1, 4] <- NA}
      if (tmp[1, 5] < 0 || tmp[1, 5] > 1 || is.na(tmp[1, 5])) {tmp[1, 5] <- NA}
      res$ci$low$freq.omega <- tmp[1, 4]
      res$ci$up$freq.omega <- tmp[1, 5]
    }
  }
  return(res)
}
