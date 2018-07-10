
data <- MASS::mvrnorm(200, rep(0, 4), psych::sim.congeneric())
a <- Sys.time()
bb <- brel(data, if.item.dropped = T, freq = F, estimates = "alpha")
Sys.time() - a

summary(bb)


