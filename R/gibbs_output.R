# this function calls on other functions in order to return the sampled estimates
# and the credible intervals together with the posterior distribution objects
# to be passed on for forther analysis

gibbsFun <- function(data, estimates, n.iter, n.burnin, thin, n.chains, interval, item.dropped, pairwise, callback){
  p <- ncol(data)
  res <- list()
  if ("alpha" %in% estimates || "lambda2" %in% estimates || "lambda4" %in% estimates || "lambda6" %in% estimates ||
      "glb" %in% estimates){
    tmp_out <- covSamp(data, n.iter, n.burnin, thin, n.chains, pairwise)
    callback()
    C <- tmp_out$cov_mat
    res$covsamp <- C
    res$data_mis_samp_cov <- tmp_out$dat_mis_samp_cov
    if (item.dropped) {
      Ctmp <- array(0, c(n.chains, (n.iter - n.burnin)/thin, p, p - 1, p - 1))
      for (i in 1:p){
        Ctmp[, , i, , ] <- C[ , , -i, -i]
      }
    }
  }

  if ("alpha" %in% estimates){
    res$samp$Bayes_alpha <- coda::mcmc(apply(C, MARGIN = c(1, 2), applyalpha))
    int <- coda::HPDinterval(chainSmoker(res$samp$Bayes_alpha), prob = interval)
    res$cred$low$Bayes_alpha <- int[1]
    res$cred$up$Bayes_alpha <- int[2]
    res$est$Bayes_alpha <- mean(res$samp$Bayes_alpha)
    if (item.dropped){
      res$ifitem$samp$alpha <- (apply(Ctmp, c(1, 2, 3), applyalpha))
      res$ifitem$est$alpha <- apply(res$ifitem$samp$alpha, 3, mean)
      res$ifitem$cred$alpha <- coda::HPDinterval(chainSmoker(res$ifitem$samp$alpha), prob = interval)
    }
    callback()
  }

  if ("lambda2" %in% estimates){
    res$samp$Bayes_lambda2 <- coda::mcmc(apply(C, MARGIN = c(1, 2), applylambda2))
    int <- coda::HPDinterval(chainSmoker(res$samp$Bayes_lambda2), prob = interval)
    res$cred$low$Bayes_lambda2 <- int[1]
    res$cred$up$Bayes_lambda2 <- int[2]
    res$est$Bayes_lambda2<- mean(res$samp$Bayes_lambda2)
    if (item.dropped){
      res$ifitem$samp$lambda2 <- apply(Ctmp, c(1, 2, 3), applylambda2)
      res$ifitem$est$lambda2 <- apply(res$ifitem$samp$lambda2, 3, mean)
      res$ifitem$cred$lambda2 <- coda::HPDinterval(chainSmoker(res$ifitem$samp$lambda2), prob = interval)
    }
    callback()
  }

  if ("lambda4" %in% estimates){
    res$samp$Bayes_lambda4 <- coda::mcmc(apply(C, MARGIN = c(1, 2), applylambda4))
    int <- coda::HPDinterval(chainSmoker(res$samp$Bayes_lambda4), prob = interval)
    res$cred$low$Bayes_lambda4 <- int[1]
    res$cred$up$Bayes_lambda4 <- int[2]
    res$est$Bayes_lambda4<- mean(res$samp$Bayes_lambda4)
    if (item.dropped){
      res$ifitem$samp$lambda4 <- (apply(Ctmp, c(1, 2, 3), applylambda4))
      res$ifitem$est$lambda4 <- apply(res$ifitem$samp$lambda4, 3, mean)
      res$ifitem$cred$lambda4 <- coda::HPDinterval(chainSmoker(res$ifitem$samp$lambda4), prob = interval)
    }
    callback()
  }

  if ("lambda6" %in% estimates){
    res$samp$Bayes_lambda6 <- coda::mcmc(apply(C, MARGIN = c(1, 2), applylambda6))
    int <- coda::HPDinterval(chainSmoker(res$samp$Bayes_lambda6), prob = interval)
    res$cred$low$Bayes_lambda6 <- int[1]
    res$cred$up$Bayes_lambda6 <- int[2]
    res$est$Bayes_lambda6<- mean(res$samp$Bayes_lambda6)
    if (item.dropped){
      res$ifitem$samp$lambda6 <- (apply(Ctmp, c(1, 2, 3), applylambda6))
      res$ifitem$est$lambda6 <- apply(res$ifitem$samp$lambda6, 3, mean)
      res$ifitem$cred$lambda6 <- coda::HPDinterval(chainSmoker(res$ifitem$samp$lambda6), prob = interval)
    }
    callback()
  }

  if ("glb" %in% estimates){
    res$samp$Bayes_glb <- coda::mcmc(apply(C, c(1, 2), glbOnArray))
    callback()
    if (sum(is.na(res$samp$Bayes_glb) > 0)) {
      int <- c(NA, NA)
    } else {
      int <- coda::HPDinterval(chainSmoker(res$samp$Bayes_glb), prob = interval)
    }
    res$cred$low$Bayes_glb <- int[1]
    res$cred$up$Bayes_glb <- int[2]
    res$est$Bayes_glb <- mean(res$samp$Bayes_glb)
    if (item.dropped){
      res$ifitem$samp$glb <- (apply(Ctmp, c(1, 2, 3), glbOnArray))
      res$ifitem$est$glb <- apply(res$ifitem$samp$glb, 3, mean)
      res$ifitem$cred$glb <- coda::HPDinterval(chainSmoker(res$ifitem$samp$glb), prob = interval)
    }
    callback()
  }

  # special case omega -----------------------------------------------------------------
  if ("omega" %in% estimates){
    om_samp <- omegaSampler(data, n.iter, n.burnin, thin, n.chains, pairwise)
    res$samp$Bayes_omega <- om_samp$omega
    res$data_mis_samp_fm <- om_samp$dat_mis_samp_fm
    res$loadings <- apply((om_samp$lambda), 3, mean)
    res$resid_var <- apply((om_samp$psi), 3, mean)

    int <- coda::HPDinterval(chainSmoker(res$samp$Bayes_omega), prob = interval)
    res$cred$low$Bayes_omega <- int[1]
    res$cred$up$Bayes_omega<- int[2]
    res$est$Bayes_omega <- mean(res$samp$Bayes_omega)

    if (item.dropped){
      om_samp_ifitem <- array(0, c(n.chains, (n.iter - n.burnin)/thin, p))
      for (i in 1:p){
        tmp <- data[-i, -i]
        om_samp_ifitem[, , i] <- omegaSampler(tmp, n.iter, n.burnin, thin, n.chains, pairwise)$omega
      }
      res$ifitem$samp$omega <- (om_samp_ifitem)
      res$ifitem$est$omega <- apply(om_samp_ifitem, 3, mean)
      res$ifitem$cred$omega <- coda::HPDinterval(chainSmoker(res$ifitem$samp$omega), prob = interval)
    }
    callback()
  }

  return(res)

}
