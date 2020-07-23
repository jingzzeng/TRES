context("Test reproducible examples in the paper.")

testthat::skip('skip')
## Set up RNG kind
RNGkind("L'Ecuyer-CMRG")

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
  fit_ols2 <- TPR.fit(square$x, square$y, method = "standard")
  fit_fg2 <- TPR.fit(square$x, square$y, u = c(2, 2), method = "FG")
  fit_1d2 <- TPR.fit(square$x, square$y, u = c(2, 2), method = "1D")
  fit_ecd2 <- TPR.fit(square$x, square$y, u = c(2, 2), method = "ECD")
  fit_pls2 <- TPR.fit(square$x, square$y, u = c(2, 2), method = "PLS")
  dist_ols2 <- rTensor::fnorm(coef(fit_ols2) - square$coefficients)
  dist_fg2 <- rTensor::fnorm(coef(fit_fg2) - square$coefficients)
  dist_1d2 <- rTensor::fnorm(coef(fit_1d2) - square$coefficients)
  dist_ecd2 <- rTensor::fnorm(coef(fit_ecd2) - square$coefficients)
  dist_pls2 <- rTensor::fnorm(coef(fit_pls2) - square$coefficients)
  Pdist_ols2 <- rep(NA_real_, 2)
  Pdist_fg2 <- rep(NA_real_, 2)
  Pdist_1d2 <- rep(NA_real_, 2)
  Pdist_ecd2 <- rep(NA_real_, 2)
  Pdist_pls2 <- rep(NA_real_, 2)
  for (i in 1:2){
    Pdist_fg2[i] <- subspace(square$Gamma[[i]],fit_fg2$Gamma[[i]])
    Pdist_1d2[i] <- subspace(square$Gamma[[i]],fit_1d2$Gamma[[i]])
    Pdist_ecd2[i] <- subspace(square$Gamma[[i]],fit_ecd2$Gamma[[i]])
    Pdist_pls2[i] <- subspace(square$Gamma[[i]],fit_pls2$Gamma[[i]])
  }
  Pdist_fg2 <- sum(Pdist_fg2)
  Pdist_1d2 <- sum(Pdist_1d2)
  Pdist_ecd2 <- sum(Pdist_ecd2)
  Pdist_pls2 <- sum(Pdist_pls2)
  expect_equal(c(dist_ols2, dist_1d2, dist_ecd2, dist_pls2, dist_ols2, Pdist_fg2, Pdist_1d2, Pdist_ecd2, Pdist_pls2), c(2225.63779502817, 5.79782569895, 5.79782589677, 5.59056824394, 2225.63779502817, 1.41425438027, 1.41425441003, 1.41425440173, 0.15631639801))
})

test_that("Section 3.5", {
  set.seed(1)
  data("EEG")
  u_eeg <- c(1,1)
  fit_eeg_ols <- TRR.fit(EEG$x, EEG$y, method = "standard")
  fit_eeg_1D <- TRR.fit(EEG$x, EEG$y, u_eeg, method = "1D")
  testthat::expect_equal(fit_eeg_ols$coefficients@data[8:12,8,1], c(-0.08669501, -0.10257079, -0.03880391, 0.02420581,  0.08765343))
  testthat::expect_equal(fit_eeg_1D$coefficients@data[8:12,25,1], c(1.1223038, 1.0996184, 0.4879741, 0.5078613, 0.6050595), tolerance = 1e-4)
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
  expect_equal(c(d1,d2,d3,d4,d5,d6,d7,d8), c(1.2998375801e-08, 4.1871487339e-13, 1.3689449324e-07, 2.5737304665e-08, 1.3151709806e-07, 2.4338580287e-08, 6.3245553203e-01, 7.2155626963e-01))
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
  expect_equal(u_est3, c(11,5))
})


