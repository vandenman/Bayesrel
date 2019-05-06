


set.seed(1234)

test_that("Reliability estimates are correct", {

  data(cavalini, package = "Bayesrel")
  ee <- Bayesrel::strel(cavalini[1:200, ], estimates = c("lambda2", "omega"), n.iter=200, n.boot=200)

  expect_equal(ee$bay$est$bayes.l2, 0.8024656, tolerance = 1e-2)
  expect_equal(ee$freq$est$freq.l2, 0.7998767, tolerance = 1e-3)
  expect_equal(ee$bay$est$bayes.omega, 0.7994719, tolerance = 1e-2)
  expect_equal(ee$freq$est$freq.omega, 0.7976158, tolerance = 1e-3)

})

set.seed(1234)

test_that("Bayes glb is correct", {

  data(cavalini, package = "Bayesrel")
  ee <- Bayesrel::strel(cavalini[1:200, ], estimates = "glb", n.iter = 200, freq = F)

  expect_equal(ee$bay$est$bayes.glb, 0.8711535, tolerance = 1e-2)


})

set.seed(1234)

test_that("Bayes Alpha is correct", {

  data(cavalini, package = "Bayesrel")
  ee <- Bayesrel::strel(cavalini, estimates = "alpha", n.iter = 200, freq = F, item.dropped = T)

  expect_equal(ee$bay$ifitem$est$alpha, c(0.7687316, 0.7325940, 0.7622475, 0.7567935, 0.7558746, 0.7517239, 0.7526119, 0.7619934),
               tolerance = 1e-2)

})

