context("Covariance matrix")

test_that("Covariance matrix", {
  data("bat")
  x <- bat$x
  y <- bat$y
  fit_1d1 <- TRR.fit(bat$x, bat$y, u = c(14,14), method = "1D")
  print(fit_1d1$call)
  # expect_equal(fit_1d1$Gamma[[1]][1:5], c(-0.00014302,0.00068032,0.00039628,0.00031333,0.00006079))
  expect_equal(fit_1d1$M_list[[1]][1:5], c(0.32888006,0.00685595,-0.09999852,-0.09008617,0.03594920))
  expect_equal(fit_1d1$U_list[[1]][1:5], c(0.02336009,0.00036831,-0.00771248,-0.00700672,0.00223131))
  expect_equal(fit_1d1$W0_list[[1]][[1]][1:5], c(0.00199160,-0.00068083,-0.00099282,-0.00080295,0.00002528))
  expect_equal(fit_1d1$Gamma[[1]][1:5,1], c(0.00014302,-0.00068032,-0.00039628,-0.00031333,-0.00006079))
})

