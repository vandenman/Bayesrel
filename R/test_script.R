
dat <- MASS::mvrnorm(100, rep(0, 4), psych::sim.congeneric())

bb <- brel(dat)
bb
summary(bb)
