


set.seed(1234)

test_that("Reliability estimates are correct", {

  data(cavalini, package = "Bayesrel")
  ee <- Bayesrel::strel(asrm, estimates = c("lambda2", "omega"), n.iter=500, n.boot=200)

  expect_equal(ee$bay$est$bayes_l2, 0.7962999, tolerance = 1e-2)
  expect_equal(ee$freq$est$freq_l2, 0.7960336, tolerance = 1e-2)
  expect_equal(ee$bay$est$bayes_omega, 0.7915491, tolerance = 1e-2)
  expect_equal(ee$freq$est$freq_omega, 0.7928258, tolerance = 1e-2)

  expect_equal(ee$bay$cred$low$bayes_l2, 0.7112308, tolerance = 1e-2)
  expect_equal(ee$bay$cred$up$bayes_omega, 0.8356543, tolerance = 1e-2)


})

set.seed(1234)

test_that("Bayes glb is correct", {

  data(cavalini, package = "Bayesrel")
  ee <- Bayesrel::strel(asrm, estimates = "glb", n.iter = 500, freq = F)

  expect_equal(ee$bay$est$bayes_glb, 0.8541144, tolerance = 1e-2)


})

set.seed(1234)

test_that("Bayes Alpha if item deleted is correct", {

  data(cavalini, package = "Bayesrel")
  ee <- Bayesrel::strel(asrm, estimates = "alpha", n.iter = 250, freq = F, item.dropped = T)

  expect_equal(as.numeric(ee$bay$ifitem$est$alpha), c(0.7229561, 0.7291643, 0.7894970, 0.7597513, 0.7320029),
               tolerance = 1e-2)

})

