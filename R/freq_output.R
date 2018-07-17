#' this function calls on other functions in order to return the fuentist estimates
#' and bootstrapped cinfidence intervals, which can be of different types

freqFun<- function(data, boot.n, boot.interval.type, estimates, interval, omega.freq.method,
                   omega.conf.int.type){
  p <- ncol(data)
  n <- nrow(data)
  res <- list()

  if ("alpha" %in% estimates){
    res$est$freq.alpha <- applyalpha(cov(data))
    f.alpha.boot.obj <- boot::boot(data, statistic = bootAlpha, R = boot.n)
    if (length(unique(round(f.alpha.boot.obj$t, 4))) == 1){
      res$ci$low$freq.alpha <- 1
      res$ci$up$freq.alpha <- 1
    }
    else{
      tmp <- boot::boot.ci(f.alpha.boot.obj, conf = interval, type = boot.interval.type)[[4]]
      res$ci$low$freq.alpha <- tmp[1, 4]
      res$ci$up$freq.alpha <- tmp[1, 5]
    }
    res$boot$alpha <- f.alpha.boot.obj$t
  }
  if ("lambda2" %in% estimates){
    res$est$freq.l2 <- applyl2(cov(data))
    f.l2.boot.obj <- boot::boot(data, statistic = bootL2, R = boot.n)
    if (length(unique(round(f.l2.boot.obj$t, 4))) == 1){
      res$ci$low$freq.l2 <- 1
      res$ci$up$freq.l2 <- 1
    }
    else{
      tmp <- boot::boot.ci(f.l2.boot.obj, conf = interval, type = boot.interval.type)[[4]]
      res$ci$low$freq.l2 <- tmp[1, 4]
      res$ci$up$freq.l2 <- tmp[1, 5]
    }
    res$boot$l2 <- f.l2.boot.obj$t
  }

  if ("lambda6" %in% estimates){
    res$est$freq.l6 <- applyl6(cov(data))
    f.l6.boot.obj <- boot::boot(data, statistic = bootL6, R = boot.n)
    if (length(unique(round(f.l6.boot.obj$t, 4))) == 1){
      res$ci$low$freq.l6 <- 1
      res$ci$up$freq.l6 <- 1
    }
    else{
      tmp <- boot::boot.ci(f.l6.boot.obj, conf = interval, type = boot.interval.type)[[4]]
      res$ci$low$freq.l6 <- tmp[1, 4]
      res$ci$up$freq.l6 <- tmp[1, 5]
    }
    res$boot$l6 <- f.l6.boot.obj$t
  }

  if ("glb" %in% estimates){
    res$est$freq.glb <- applyglb(cov(data))
    f.glb.boot.obj <- boot::boot(data, statistic = bootGlb, R = boot.n)
    if (length(unique(round(f.glb.boot.obj$t, 4))) == 1){
      res$ci$low$freq.glb <- 1
      res$ci$up$freq.glb <- 1
    }
    else{
      tmp <- boot::boot.ci(f.glb.boot.obj, conf = interval, type = boot.interval.type)[[4]]
      res$ci$low$freq.glb <- tmp[1, 4]
      res$ci$up$freq.glb <- tmp[1, 5]
    }
    res$boot$glb <- f.glb.boot.obj$t
  }

  #omega --------------------------------------------------------------------------
  if ("omega" %in% estimates){
    if (omega.freq.method == "cfa"){
      res$est$freq.omega <- applyomega_boot_cfa(data)
      if (omega.conf.int.type == "alg"){
        omega.alg <- applyomega_alg(data, interval)
        # res$est$freq.omega.alg <- omega.alg[1]
        res$ci$low$freq.omega <- omega.alg[2]
        res$ci$up$freq.omega <- omega.alg[3]
      }
      else{
        f.omega.boot.obj <- boot::boot(data, statistic = bootOmega_cfa, R = boot.n)
        tmp <- boot::boot.ci(f.omega.boot.obj, conf = interval, type = boot.interval.type)[[4]]
        if (is.na(tmp[1, 4])) {tmp[1, 4] <- NA}
        if (is.na(tmp[1, 5])) {tmp[1, 5] <- NA}
        res$ci$low$freq.omega <- tmp[1, 4]
        res$ci$up$freq.omega <- tmp[1, 5]
        res$boot$omega <- f.omega.boot.obj$t
      }
    }
    if (omega.freq.method == "pa"){
      res$est$freq.omega <- applyomega_boot_pa(cov(data))
      f.omega.boot.obj <- boot::boot(data, statistic = bootOmega_pa, R = boot.n)
      tmp <- boot::boot.ci(f.omega.boot.obj, conf = interval, type = boot.interval.type)[[4]]
      if (tmp[1, 4] < 0 || tmp[1, 4] > 1 || is.na(tmp[1, 4])) {tmp[1, 4] <- NA}
      if (tmp[1, 5] < 0 || tmp[1, 5] > 1 || is.na(tmp[1, 5])) {tmp[1, 5] <- NA}
      res$ci$low$freq.omega <- tmp[1, 4]
      res$ci$up$freq.omega <- tmp[1, 5]
      res$boot$omega <- f.omega.boot.obj$t
    }
    res$omega.freq.method <- omega.freq.method
  }
  return(res)
}
