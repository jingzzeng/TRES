context("Test reproducible examples in the paper.")
# library("TRES")

## Skip if the version is not 1.1.1
# testthat::skip('skip')
# testthat::skip_if_not(as.numeric(strsplit(packageDescription("TRES")$Version, "\\.")[[1]][3]) == 1)


test_that("Section 3.1", {
  set.seed(1)
  data("bat")
  fit_ols1 <- TRR.fit(bat$x, bat$y, method = 'standard')
  fit_pls1 <- TRR.fit(bat$x, bat$y, u = c(14, 14), method = 'PLS')
  dist_ols1 <- rTensor::fnorm(coef(fit_ols1) - bat$coefficients)
  dist_pls1 <- rTensor::fnorm(coef(fit_pls1) - bat$coefficient)
  Pdist_pls1 <- rep(NA_real_, 2)
  for (i in 1:2){
    Pdist_pls1[i] <- subspace(bat$Gamma[[i]],fit_pls1$Gamma[[i]])
  }
  Pdist_pls1 <- sum(Pdist_pls1)
  expect_equal(c(dist_ols1, dist_pls1, Pdist_pls1), c(0.97204002561, 1.19191180720, 1.46626921297))
})

test_that("Section 3.2", {
  set.seed(1)
  data("square")
  fit_ecd2 <- TPR.fit(square$x, square$y, u = c(2, 2), method = "ECD")
  fit_pls2 <- TPR.fit(square$x, square$y, u = c(2, 2), method = "PLS")
  fit_ols2 <- TPR.fit(square$x, square$y, method = "standard")
  dist_ecd2 <- rTensor::fnorm(coef(fit_ecd2) - square$coefficients)
  dist_pls2 <- rTensor::fnorm(coef(fit_pls2) - square$coefficients)
  dist_ols2 <- rTensor::fnorm(coef(fit_ols2) - square$coefficients)
  Pdist_ecd2 <- rep(NA_real_, 2)
  Pdist_pls2 <- rep(NA_real_, 2)
  Pdist_ols2 <- rep(NA_real_, 2)
  for (i in 1:2){
    Pdist_ecd2[i] <- subspace(square$Gamma[[i]],fit_ecd2$Gamma[[i]])
    Pdist_pls2[i] <- subspace(square$Gamma[[i]],fit_pls2$Gamma[[i]])
  }
  Pdist_ecd2 <- sum(Pdist_ecd2)
  Pdist_pls2 <- sum(Pdist_pls2)
  expect_equal(c(dist_ecd2, dist_pls2, dist_ols2, Pdist_ecd2, Pdist_pls2), c(5.80500045743, 5.59056824394, 2225.63779502817, 1.41448470979, 0.15631639801))
})

test_that("Section 3.5", {
  set.seed(1)
  data("EEG")
  u_eeg <- c(1,1)
  fit_eeg_ols <- TRR.fit(EEG$x, EEG$y, method = "standard")
  fit_eeg_1D <- TRR.fit(EEG$x, EEG$y, u_eeg, method = "1D")
  testthat::expect_equal(fit_eeg_ols$coefficients@data[8:12,8,1], c(-0.366630830, -0.261239935, -0.219068713, -0.446891195, -0.052921259))
  testthat::expect_equal(fit_eeg_1D$coefficients@data[8:12,25,1], c(0.89687379902, 0.90835127523, 0.38472666044, 0.41908589959, 0.47720497460))
})

test_that("Section 4.3", {
  set.seed(1)
  p <- 20
  u <- 5
  # Model M1
  data <- MenvU_sim(p, u, jitter = 1e-5)
  Gamma <- data$Gamma
  M <- data$M
  U <- data$U

  G1 <- simplsMU(M, U, u)
  G2 <- ECD(M, U, u)
  G3 <- manifold1D(M, U, u)
  G4 <- OptM1D(M, U, u)
  G5 <- manifoldFG(M, U, u)
  G6 <- OptMFG(M, U, u)
  d1 <- subspace(G1, Gamma)
  d2 <- subspace(G2, Gamma)
  d3 <- subspace(G3, Gamma)
  d4 <- subspace(G4, Gamma)
  d5 <- subspace(G5, Gamma)
  d6 <- subspace(G6, Gamma)

  # The randomly generated matrix
  A <- matrix(runif(p*u), p, u)
  G7 <- manifoldFG(M, U, u, Gamma_init = A)
  G8 <- OptMFG(M, U, u, Gamma_init = A)
  d7 <- subspace(G7, Gamma)
  d8 <- subspace(G8, Gamma)
  expect_equal(c(d1,d2,d3,d4,d5,d6,d7,d8), c(1.4656091082e-09, 5.7760373252e-13, 5.1506432311e-10, 1.9246746873e-09, 4.5243985165e-10, 1.7448436515e-09, 4.4721359550e-01, 6.0979971922e-01))
})


test_that("Section 4.4", {
  set.seed(1)
  p <- 50
  u <- 5
  n0 <- c(50,200)
  u_est3 <- rep(NA_integer_, length(n0))
  for (i in seq_along(n0)){
    n <- n0[i]
    data <- MenvU_sim(p, u, jitter = 1e-5, wishart = TRUE, n = n)
    M <- data$M
    U <- data$U
    output <- oneD_bic(M, U, n, maxdim = p/2)
    u_est3[i] <- output$u
  }
  expect_equal(u_est3, c(9,5))
})


