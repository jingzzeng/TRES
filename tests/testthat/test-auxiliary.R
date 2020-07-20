context("Test auxiliary functions")
# library("TRES")

# testthat::skip('skip')
# testthat::skip_if_not(as.numeric(strsplit(packageDescription("TRES")$Version, "\\.")[[1]][3]) == 1)

test_that("Tenv_Pval function works", {
  set.seed(1)
  data("bat")
  x <- bat$x
  y <- bat$y
  fit_std <- TRR.fit(x, y, method="standard")
  fit <- Tenv_Pval(x, y, fit_std$coefficients)
  expect_equal(fit$se@data[51:55,3,1], c(0.059671626039, 0.067443096107, 0.063656446666, 0.068024423057, 0.057238323458))
  expect_equal(fit$p_val@data[6:10,10,1], c(0.56343181958, 0.51710260368, 0.30713200748, 0.66371425932 , 0.99487291620))
})

test_that("TRRdim works", {
  set.seed(1)
  r <- c(10, 10, 10)
  u <- c(2, 2, 2)
  p <- 5
  n <- 100
  dat <- TRRsim(r = r, p = p, u = u, n = n)
  x <- dat$x
  y <- dat$y
  result <- TRRdim(x, y)
  expect_equal(result$mse, 2465.4516196)
  expect_equal(result$bicval, c(-1061.34245959, -990.29055815, -1126.47314475))
  expect_equal(result$u, c(2,2,2))
})

test_that("PMSE works", {
  set.seed(1)
  p <- c(10, 10, 10)
  u <- c(1, 1, 1)
  r <- 5
  n <- 200
  dat <- TPRsim(p = p, r = r, u = u, n = n)
  x <- dat$x
  y <- dat$y
  fit_std <- TPR.fit(x, y, u, method="standard")
  result <- PMSE(x, y, fit_std$coefficients)
  expect_equal(result$mse, 384.69608034)
  expect_equal(result$pred[3, 51:55], c(-8.7722310529, 5.4961530259, 9.2623739155, 1.8764567351, -3.1972341381))
})
