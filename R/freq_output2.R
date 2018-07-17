

#' this function calls on other functions in order to return the fuentist estimates
#' and bootstrapped confidence intervals, which can be of different types

freqFun2<- function(data, boot.n, estimates, interval, omega.freq.method,
                   omega.conf.int.type, if.item.dropped){
  p <- ncol(data)
  n <- nrow(data)
  res <- list()
  boot.data <- array(0, c(boot.n, n, p))
  boot.cov <- array(0, c(boot.n, p, p))
  for (i in 1:boot.n){
    boot.data[i, , ] <- as.matrix(data[sample.int(nrow(data), size = n, replace = TRUE), ])
    boot.cov[i, , ] <- cov(boot.data[i, , ])
  }
  res$boot$C <- boot.cov
  if (if.item.dropped){
    Ctmp <- array(0, c(p, p - 1, p - 1))
    Dtmp <- array(0, c(p, n, p - 1))
    for (i in 1:p){
      Ctmp[i, , ] <- cov(data)[-i, -i]
      Dtmp[i, , ] <- data[, -i]
    }
  }
  if ("alpha" %in% estimates){
    res$est$freq.alpha <- applyalpha(cov(data))
    alpha.obj <- apply(boot.cov, 1, applyalpha)
    if (length(unique(round(alpha.obj, 4))) == 1){
      res$ci$low$freq.alpha <- 1
      res$ci$up$freq.alpha <- 1
    }
    else{
      res$ci$low$freq.alpha <- quantile(alpha.obj, probs = (1 - interval)/2, na.rm = T)
      res$ci$up$freq.alpha <- quantile(alpha.obj, probs = interval + (1 - interval)/2, na.rm = T)
    }
    res$boot$alpha <- alpha.obj
    if (if.item.dropped){
      res$ifitem$alpha <- apply(Ctmp, 1, applyalpha)
    }
  }
  if ("lambda2" %in% estimates){
    res$est$freq.l2 <- applyl2(cov(data))
    l2.obj <- apply(boot.cov, 1, applyl2)
    if (length(unique(round(l2.obj, 4))) == 1){
      res$ci$low$freq.l2 <- 1
      res$ci$up$freq.l2 <- 1
    }
    else{
      res$ci$low$freq.l2 <- quantile(l2.obj, probs = (1 - interval)/2, na.rm = T)
      res$ci$up$freq.l2 <- quantile(l2.obj, probs = interval + (1 - interval)/2, na.rm = T)
    }
    res$boot$l2 <- l2.obj
    if (if.item.dropped){
      res$ifitem$l2 <- apply(Ctmp, 1, applyl2)
    }
  }

  if ("lambda4" %in% estimates){
    res$est$freq.l4 <- applyl4(cov(data))
    l4.obj <- apply(boot.cov, 1, applyl4)
    if (length(unique(round(l4.obj, 4))) == 1){
      res$ci$low$freq.l4 <- 1
      res$ci$up$freq.l4 <- 1
    }
    else{
      res$ci$low$freq.l4 <- quantile(l4.obj, probs = (1 - interval)/2, na.rm = T)
      res$ci$up$freq.l4 <- quantile(l4.obj, probs = interval + (1 - interval)/2, na.rm = T)
    }
    res$boot$l4 <- l4.obj
    if (if.item.dropped){
      res$ifitem$l4 <- apply(Ctmp, 1, applyl4)
    }
  }

  if ("lambda6" %in% estimates){
    res$est$freq.l6 <- applyl6(cov(data))
    l6.obj <- apply(boot.cov, 1, applyl6)
    if (length(unique(round(l6.obj, 4))) == 1){
      res$ci$low$freq.l6 <- 1
      res$ci$up$freq.l6 <- 1
    }
    else{
      res$ci$low$freq.l6 <- quantile(l6.obj, probs = (1 - interval)/2, na.rm = T)
      res$ci$up$freq.l6 <- quantile(l6.obj, probs = interval + (1 - interval)/2, na.rm = T)
    }
    res$boot$l6 <- l6.obj
    if (if.item.dropped){
      res$ifitem$l6 <- apply(Ctmp, 1, applyl6)
    }
  }
  if ("glb" %in% estimates){
    res$est$freq.glb <- applyglb(cov(data))
    glb.obj <- apply(boot.cov, 1, applyglb)
    if (length(unique(round(glb.obj, 4))) == 1){
      res$ci$low$freq.glb <- 1
      res$ci$up$freq.glb <- 1
    }
    else{
      res$ci$low$freq.glb <- quantile(glb.obj, probs = (1 - interval)/2, na.rm = T)
      res$ci$up$freq.glb <- quantile(glb.obj, probs = interval + (1 - interval)/2, na.rm = T)
    }
    res$boot$glb <- glb.obj
    if (if.item.dropped){
      res$ifitem$glb <- apply(Ctmp, 1, applyglb)
    }
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
        res$boot$omega <- NULL
      }
      else{
        omega.obj <- apply(boot.data, 1, applyomega_boot_cfa)
        if (length(unique(round(omega.obj, 4))) == 1){
          res$ci$low$freq.omega <- 1
          res$ci$up$freq.omega <- 1
        }
        else{
          res$ci$low$freq.omega <- quantile(omega.obj, probs = (1 - interval)/2, na.rm = T)
          res$ci$up$freq.omega <- quantile(omega.obj, probs = interval + (1 - interval)/2, na.rm = T)
        }
        res$boot$omega <- omega.obj
      }
      if (if.item.dropped){
        res$ifitem$omega <- apply(Dtmp, 1, applyomega_boot_cfa)
      }
    }
    if (omega.freq.method == "pa"){
      omega.obj <- apply(boot.cov, 1, applyomega_boot_pa)
      if (length(unique(round(omega.obj, 4))) == 1){
        res$ci$low$freq.omega <- 1
        res$ci$up$freq.omega <- 1
      }
      else{
        res$ci$low$freq.omega <- quantile(omega.obj, probs = (1 - interval)/2, na.rm = T)
        res$ci$up$freq.omega <- quantile(omega.obj, probs = interval + (1 - interval)/2, na.rm = T)
      }
      res$boot$omega <- omega.obj
      res$est$freq.omega <- applyomega_boot_pa(cov(data))
      if (if.item.dropped){
        res$ifitem$omega <- apply(Ctmp, 1, applyomega_boot_pa)
      }
    }
    res$omega.freq.method <- omega.freq.method
  }
  return(res)
}
