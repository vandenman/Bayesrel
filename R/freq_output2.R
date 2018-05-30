

#' this function calls on other functions in order to return the fuentist estimates
#' and bootstrapped confidence intervals, which can be of different types

freqFun2<- function(data, boot.n, estimates, interval, omega.freq.method,
                   omega.conf.type){
  p <- ncol(data)
  n <- nrow(data)
  res <- list()
  boot.data <- array(0, c(boot.n, 1000, p))
  boot.cov <- array(0, c(boot.n, p, p))
  for (i in 1:boot.n){
    boot.data[i, , ] <- as.matrix(dplyr::sample_n(as.data.frame(data), size = 1000, replace=TRUE))
    boot.cov[i, , ] <- cov(boot.data[i, , ])
  }
  if ("alpha" %in% estimates){
    res$est$freq.alpha <- applyAlpha(cov(data))
    alpha.obj <- apply(boot.cov, 1, applyAlpha)
    if (length(unique(round(alpha.obj, 4))) == 1){
      res$ci$low$freq.alpha <- 1
      res$ci$up$freq.alpha <- 1
    }
    else{
      res$ci$low$freq.alpha <- quantile(alpha.obj, probs = (1 - interval)/2)
      res$ci$up$freq.alpha <- quantile(alpha.obj, probs = interval + (1 - interval)/2)
    }
    res$boot$alpha <- alpha.obj
  }
  if ("l2" %in% estimates){
    res$est$freq.l2 <- applyL2(cov(data))
    l2.obj <- apply(boot.cov, 1, applyL2)
    if (length(unique(round(l2.obj, 4))) == 1){
      res$ci$low$freq.l2 <- 1
      res$ci$up$freq.l2 <- 1
    }
    else{
      res$ci$low$freq.l2 <- quantile(l2.obj, probs = (1 - interval)/2)
      res$ci$up$freq.l2 <- quantile(l2.obj, probs = interval + (1 - interval)/2)
    }
    res$boot$l2 <- l2.obj
  }

  if ("l6" %in% estimates){
    res$est$freq.l6 <- applyL6(cov(data))
    l6.obj <- apply(boot.cov, 1, applyL6)
    if (length(unique(round(l6.obj, 4))) == 1){
      res$ci$low$freq.l6 <- 1
      res$ci$up$freq.l6 <- 1
    }
    else{
      res$ci$low$freq.l6 <- quantile(l6.obj, probs = (1 - interval)/2)
      res$ci$up$freq.l6 <- quantile(l6.obj, probs = interval + (1 - interval)/2)
    }
    res$boot$l6 <- l6.obj
  }
  if ("glb" %in% estimates){
    res$est$freq.glb <- applyGlb(cov(data))
    glb.obj <- apply(boot.cov, 1, applyGlb)
    if (length(unique(round(glb.obj, 4))) == 1){
      res$ci$low$freq.glb <- 1
      res$ci$up$freq.glb <- 1
    }
    else{
      res$ci$low$freq.glb <- quantile(glb.obj, probs = (1 - interval)/2)
      res$ci$up$freq.glb <- quantile(glb.obj, probs = interval + (1 - interval)/2)
    }
    res$boot$glb <- glb.obj
  }

  #omega --------------------------------------------------------------------------
  if ("omega" %in% estimates){
    if (omega.freq.method == "cfa"){
      res$est$freq.omega <- applyOmega_boot_cfa(data)
      if (omega.conf.type == "alg"){
        omega.alg <- applyOmega_alg(data, interval)
        # res$est$freq.omega.alg <- omega.alg[1]
        res$ci$low$freq.omega <- omega.alg[2]
        res$ci$up$freq.omega <- omega.alg[3]
        res$boot$omega <- NULL
      }
      else{
        omega.obj <- apply(boot.data, 1, applyOmega_boot_cfa)
        if (length(unique(round(omega.obj, 4))) == 1){
          res$ci$low$freq.omega <- 1
          res$ci$up$freq.omega <- 1
        }
        else{
          res$ci$low$freq.omega <- quantile(omega.obj, probs = (1 - interval)/2)
          res$ci$up$freq.omega <- quantile(omega.obj, probs = interval + (1 - interval)/2)
        }
        res$boot$omega <- omega.obj
      }
    }
    if (omega.freq.method == "pa"){
      omega.obj <- apply(boot.cov, 1, applyOmega_boot_pa)
      if (length(unique(round(omega.obj, 4))) == 1){
        res$ci$low$freq.omega <- 1
        res$ci$up$freq.omega <- 1
      }
      else{
        res$ci$low$freq.omega <- quantile(omega.obj, probs = (1 - interval)/2)
        res$ci$up$freq.omega <- quantile(omega.obj, probs = interval + (1 - interval)/2)
      }
      res$boot$omega <- omega.obj
      res$est$freq.omega <- applyOmega_boot_pa(cov(data))
    }
    res$omega.freq.method <- omega.freq.method
  }
  return(res)
}
