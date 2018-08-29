


# data <- MASS::mvrnorm(100, rep(0, 9), psych::sim.hierarchical())
#
# a <- Sys.time()
# x <- bayesrel::brel(data, item.dropped = T, prior.samp = T)
# Sys.time() - a
#
# summary(x)
#
# plotBrel(x, "lambda2", criteria = F, blackwhite = F, greek = T)
# plotIfItem_one(x, "lambda6", 1, criteria = F, blackwhite = F, greek = T)
# plotIfItem_all(x, "lambda2")
