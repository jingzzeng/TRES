context("Covariance matrix")

test_that("Covariance matrix", {
  data("bat")
  x <- bat$x
  y <- bat$y
  fit_1d1 <- TRR.fit(bat$x, bat$y, u = c(14,14), method = "1D")
  # expect_equal(fit_1d1$Gamma[[1]][1:5], c(-0.00014302,0.00068032,0.00039628,0.00031333,0.00006079))
  expect_equal(fit_1d1$M_list[[1]][1:5], c(0.32888006,0.00685595,-0.09999852,-0.09008617,0.03594920))
  expect_equal(fit_1d1$U_list[[2]][1:5], c(0.00963257,0.00692476,-0.00068983,-0.00242814,0.00285420))
})
