# this function calls on other functions in order to return the sampled estimates
# and the credible intervals together with the posterior distribution objects
# to be passed on for forther analysis

gibbsFun <- function(data, n.iter, n.burnin, estimates, interval, item.dropped){
  p <- ncol(data)
  res <- list()
  if ("alpha" %in% estimates || "lambda2" %in% estimates || "lambda4" %in% estimates || "lambda6" %in% estimates ||
      "glb" %in% estimates){
    C <- covSamp(data, n.iter, n.burnin)
    if (item.dropped) {
      Ctmp <- array(0, c(p, n.iter - n.burnin, p - 1, p - 1))
      for (i in 1:p){
        Ctmp[i, , , ] <- C[, -i, -i]
      }
    }
  } else {
    C = NULL
  }
  res$covsamp <- C

  if ("alpha" %in% estimates){
    res$samp$bayes_alpha <- coda::mcmc(apply(C, MARGIN = 1, applyalpha))
    int <- coda::HPDinterval(res$samp$bayes_alpha, prob = interval)
    res$cred$low$bayes_alpha <- int[1]
    res$cred$up$bayes_alpha <- int[2]
    res$est$bayes_alpha<- median(res$samp$bayes_alpha)
    if (item.dropped){
      res$ifitem$samp$alpha <- coda::mcmc(apply(Ctmp, c(2, 1), applyalpha))
      res$ifitem$est$alpha <- apply(res$ifitem$samp$alpha, 2, median)
    }
  }

  if ("lambda2" %in% estimates){
    res$samp$bayes_l2 <- coda::mcmc(apply(C, MARGIN = 1, applyl2))
    int <- coda::HPDinterval(res$samp$bayes_l2, prob = interval)
    res$cred$low$bayes_l2 <- int[1]
    res$cred$up$bayes_l2 <- int[2]
    res$est$bayes_l2<- median(res$samp$bayes_l2)
    if (item.dropped){
      res$ifitem$samp$l2 <- coda::mcmc(apply(Ctmp, c(2, 1), applyl2))
      res$ifitem$est$l2 <- apply(res$ifitem$samp$l2, 2, median)
    }
  }

  if ("lambda4" %in% estimates){
    res$samp$bayes_l4 <- coda::mcmc(apply(C, MARGIN = 1, applyl4))
    int <- coda::HPDinterval(res$samp$bayes_l4, prob = interval)
    res$cred$low$bayes_l4 <- int[1]
    res$cred$up$bayes_l4 <- int[2]
    res$est$bayes_l4<- median(res$samp$bayes_l4)
    if (item.dropped){
      res$ifitem$samp$l4 <- coda::mcmc(apply(Ctmp, c(2, 1), applyl4))
      res$ifitem$est$l4 <- apply(res$ifitem$samp$l4, 2, median)
    }
  }

  if ("lambda6" %in% estimates){
    res$samp$bayes_l6 <- coda::mcmc(apply(C, MARGIN = 1, applyl6))
    int <- coda::HPDinterval(res$samp$bayes_l6, prob = interval)
    res$cred$low$bayes_l6 <- int[1]
    res$cred$up$bayes_l6 <- int[2]
    res$est$bayes_l6<- median(res$samp$bayes_l6)
    if (item.dropped){
      res$ifitem$samp$l6 <- coda::mcmc(apply(Ctmp, c(2, 1), applyl6))
      res$ifitem$est$l6 <- apply(res$ifitem$samp$l6, 2, median)
    }
  }

  if ("glb" %in% estimates){
    res$samp$bayes_glb <- coda::mcmc(glbOnArray(C))
    if (sum(is.na(res$samp$bayes_glb) > 0)) {
      int <- c(NA, NA)
    } else {
      int <- coda::HPDinterval(res$samp$bayes_glb, prob = interval)
    }
    res$cred$low$bayes_glb <- int[1]
    res$cred$up$bayes_glb <- int[2]
    res$est$bayes_glb<- median(res$samp$bayes_glb)
    if (item.dropped){
      res$ifitem$samp$glb <- coda::mcmc(glbOnArray(Ctmp))
      res$ifitem$est$glb <- apply(res$ifitem$samp$glb, 2, median)
    }
  }

  # special case omega -----------------------------------------------------------------
  if ("omega" %in% estimates){
    om_samp <- omegaSampler(data, n.iter, n.burnin)
    res$samp$bayes_omega <- coda::mcmc(om_samp$omega)
    res$loadings <- apply(om_samp$lambda, 2, mean)
    res$resid_var <- apply(om_samp$psi, 2, mean)
    # res$loadings <- coda::mcmc(om_samp$lambda)
    # res$resid_var <- coda::mcmc(om_samp$psi)
    int <- coda::HPDinterval(res$samp$bayes_omega, prob = interval)
    res$cred$low$bayes_omega <- int[1]
    res$cred$up$bayes_omega<- int[2]
    res$est$bayes_omega <- mean(res$samp$bayes_omega)

    if (item.dropped){
      om_samp_ifitem <- matrix(0, n.iter - n.burnin, p)
      for (i in 1:p){
        tmp <- data[-i, -i]
        om_samp_ifitem[, i] <- omegaSampler(tmp, n.iter, n.burnin)$omega
      }
      res$ifitem$samp$omega <- coda::mcmc(om_samp_ifitem)
      res$ifitem$est$omega <- apply(om_samp_ifitem, 2, mean)
    }
  }

  return(res)

}
