context("Test TRR and TPR with .sim function, bat and square.")
library("TRES")

testthat::skip('skip')
# skip_if_not(as.numeric(strsplit(packageDescription("TRES")$Version, "\\.")[[1]][3]) == 1)

test_that("TRR and TPR works with .sim function", {
  ## TRR
  set.seed(1)
  r <- c(10, 10, 10)
  u <- c(2, 2, 2)
  p <- 5
  n <- 100
  dat <- TRR_sim(r = r, p = p, u = u, n = n)
  x <- dat$x
  y <- dat$y
  B <- dat$coefficients
  fit_std <- TRR.fit(x, y, method="standard")
  fit_1D <- TRR.fit(x, y, u, method="1D")
  expect_equal(fit_std$fitted.values@data[1:5,2,2,2], c(0.60646833, -0.25093892, -0.19145049, -0.22467693, -0.05699922))
  expect_equal(fit_1D$fitted.values@data[1:5,2,2,2], c(-0.022359406, -0.026888949, -0.055316725, -0.064956123, -0.036918492))
  expect_equal(rTensor::fnorm(B-stats::coef(fit_std)), 9.9016367)
  expect_equal(rTensor::fnorm(B-stats::coef(fit_1D)), 0.22109907)

  ## TPR
  rm(list = ls())
  p <- c(10, 10, 10)
  u <- c(1, 1, 1)
  r <- 5
  n <- 200
  dat <- TPR_sim(p = p, r = r, u = u, n = n)
  x <- dat$x
  y <- dat$y
  B <- dat$coefficients
  fit_std <- TPR.fit(x, y, method="standard")
  fit_pls <- TPR.fit(x, y, u, method="PLS")
  expect_equal(fit_std$fitted.values[3, 51:55], c(8.6413847, -2.2877337, 10.5641821, 8.4774664, 5.4915451))
  expect_equal(fit_pls$fitted.values[3, 51:55], c(3.7562096, -3.0230371, 10.7840372, 6.2212578, 2.7932273))
  expect_equal(rTensor::fnorm(B-stats::coef(fit_std)), 2516.9442)
  expect_equal(rTensor::fnorm(B-stats::coef(fit_pls)), 9.34154)
})

test_that("TRR and TPR works with bat and square datasets", {
  ## bat
  set.seed(1)
  data("bat")
  x <- bat$x
  y <- bat$y
  fit_std <- TRR.fit(x, y, method="standard")
  expect_equal(fit_std$fitted.values@data[1:5,3,2], c(-0.00362305969, 0.00162395994, 0.00073736371, 0.00332445728, -0.00217930245))

  ## square
  data("square")
  x <- square$x
  y <- square$y
  fit_std <- TPR.fit(x, y, method="standard")
  expect_equal(fit_std$fitted.values[1, 51:55], c(15.090976910, -0.046982797, 12.441433679, -51.758687874, -25.743955453))
})
