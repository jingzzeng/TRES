context("Covariance matrix")

test_that("Covariance matrix", {
  data("bat")
  x <- bat$x
  y <- bat$y
  fit_std <- TRR.fit(x, y, method="standard")
  expect_equal(fit_std$Sigma[[1]][1:5], c(0.35754933,0.00745360,-0.10871563,-0.09793919,0.03908298))
})
