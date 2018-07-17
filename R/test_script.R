
data <- MASS::mvrnorm(200, rep(0, 4), psych::sim.congeneric())
a <- Sys.time()
bb <- brel(data)
Sys.time() - a

summary(bb)

