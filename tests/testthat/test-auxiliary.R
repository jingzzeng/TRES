context("Test auxiliary functions")
library("TRES")

testthat::skip('skip')
# testthat::skip_if_not(as.numeric(strsplit(packageDescription("TRES")$Version, "\\.")[[1]][3]) == 1)

test_that("Tenv_Pval function works", {
  set.seed(1)
  data("bat")
  x <- bat$x
  y <- bat$y
  fit_std <- TRR.fit(x, y, method="standard")
  fit <- Tenv_Pval(x, y, fit_std$coefficients)
  expect_equal(fit$se@data[51:55,3,1], c(0.68712685, 0.74488826, 0.66162475, 0.77511230, 0.66944955))
  expect_equal(fit$p_val@data[6:10,10,1], c(0.075987612, 0.710086051, 0.959610679, 0.276910535, 0.644239021))
})

test_that("TensEnv_dim works", {
  set.seed(1)
  r <- c(10, 10, 10)
  u <- c(2, 2, 2)
  p <- 5
  n <- 100
  dat <- TRR_sim(r = r, p = p, u = u, n = n)
  x <- dat$x
  y <- dat$y
  result <- TensEnv_dim(x, y)
  expect_equal(result, c(2,2,2))
})

test_that("PMSE works", {
  set.seed(1)
  p <- c(10, 10, 10)
  u <- c(1, 1, 1)
  r <- 5
  n <- 200
  dat <- TPR_sim(p = p, r = r, u = u, n = n)
  x <- dat$x
  y <- dat$y
  fit_std <- TPR.fit(x, y, u, method="standard")
  result <- PMSE(x, y, fit_std$coefficients)
  expect_equal(result$mse, 384.69601)
  expect_equal(result$pred[3, 51:55], c(-8.7722278, 5.4961524, 9.2623723, 1.8764552, -3.1972383))
})
