context("Covariance matrix")

test_that("Covariance matrix", {
  data("bat")
  x <- bat$x
  y <- bat$y
  fit_1d1 <- TRR.fit(bat$x, bat$y, u = c(14,14), method = "1D")
  expect_equal(fit_1d1$Gamma[[1]][1:5], c(-0.00014302,0.00068032,0.00039628,0.00031333,0.00006079))
})
