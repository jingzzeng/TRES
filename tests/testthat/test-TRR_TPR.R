context("Test TRR and TPR with .sim function, bat and square.")
# library("TRES")

# testthat::skip('skip')
# skip_if_not(as.numeric(strsplit(packageDescription("TRES")$Version, "\\.")[[1]][3]) == 1)

test_that("TRR and TPR works with .sim function", {
  ## TRR
  set.seed(1)
  r <- c(10, 10, 10)
  u <- c(2, 2, 2)
  p <- 5
  n <- 100
  dat <- TRRsim(r = r, p = p, u = u, n = n)
  x <- dat$x
  y <- dat$y
  B <- dat$coefficients
  fit_std <- TRR.fit(x, y, method="standard")
  fit_1D <- TRR.fit(x, y, u, method="1D")
  expect_equal(fit_std$fitted.values@data[1:5,2,2,2], c(0.60646833, -0.25093892, -0.19145049, -0.22467693, -0.05699922))
  expect_equal(fit_1D$fitted.values@data[1:5,2,2,2], c(-0.022359405833, -0.026888949671, -0.055316716989, -0.064956126162, -0.036918476815))
  expect_equal(rTensor::fnorm(B-stats::coef(fit_std)), 9.9016367)
  expect_equal(rTensor::fnorm(B-stats::coef(fit_1D)), 0.22109903471)

  ## TPR
  set.seed(1)
  p <- c(10, 10, 10)
  u <- c(1, 1, 1)
  r <- 5
  n <- 200
  dat <- TPRsim(p = p, r = r, u = u, n = n)
  x <- dat$x
  y <- dat$y
  B <- dat$coefficients
  fit_std <- TPR.fit(x, y, method="standard")
  fit_pls <- TPR.fit(x, y, u, method="PLS")
  expect_equal(fit_std$fitted.values[3, 51:55], c(-8.7722310529, 5.4961530259, 9.2623739155, 1.8764567351, -3.1972341381))
  expect_equal(fit_pls$fitted.values[3, 51:55], c(-11.7888643451, 16.9736184146, 1.1881597301, 8.2860773355, -1.9265788809))
  expect_equal(rTensor::fnorm(B-stats::coef(fit_std)), 2424.2593189)
  expect_equal(rTensor::fnorm(B-stats::coef(fit_pls)), 8.1809065472)
})

test_that("TRR and TPR works with bat and square datasets", {
  ## bat
  set.seed(1)
  data("bat")
  x <- bat$x
  y <- bat$y
  fit_std <- TRR.fit(x, y, method="standard")
  expect_equal(fit_std$fitted.values@data[1:5,3,2], c(0.005919749, -0.004850022, -0.001173755, -0.002218994, -0.003877271))

  ## square
  data("square")
  x <- square$x
  y <- square$y
  fit_std <- TPR.fit(x, y, method="standard")
  expect_equal(fit_std$fitted.values[1, 51:55], c(15.090983362073, -0.046981801747, 12.441415239454, -51.758700615955, -25.743953633943))
})
