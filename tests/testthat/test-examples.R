context("Test reproducible examples in the paper.")
library("TRES")

## Skip if the version is not 1.1.1
testthat::skip('skip')
# testthat::skip_if_not(as.numeric(strsplit(packageDescription("TRES")$Version, "\\.")[[1]][3]) == 1)


test_that("Section 3.1", {
  set.seed(1)
  data("bat")
  u1 <- c(14, 14)
  fit_ols1 <- TRR.fit(bat$x, bat$y, method = 'standard')
  fit_pls1 <- TRR.fit(bat$x, bat$y, u1, method = 'PLS')
  dist_ols1 <- rTensor::fnorm(coef(fit_ols1) - bat$coefficients)
  dist_pls1 <- rTensor::fnorm(coef(fit_pls1) - bat$coefficient)
  Pdist_pls1 <- NULL
  for (i in 1:2){
    Pdist_pls1 <- c(Pdist_pls1, subspace(bat$Gamma[[i]],fit_pls1$Gamma[[i]]))
  }
  Pdist_pls1 <- sum(Pdist_pls1)
  expect_equal(c(dist_ols1, dist_pls1, Pdist_pls1), c(7.1396255, 9.4505919, 1.4238539))
})

test_that("Section 3.2", {
  set.seed(1)
  data("square")
  u2 <- c(2, 2)
  fit_ecd2 <- TPR.fit(square$x, square$y, u2, method = "ECD")
  fit_pls2 <- TPR.fit(square$x, square$y, u2, method = "PLS")
  dist_ecd2 <- rTensor::fnorm(coef(fit_ecd2) - square$coefficients)
  dist_pls2 <- rTensor::fnorm(coef(fit_pls2) - square$coefficients)
  Pdist_ecd2 <- NULL
  Pdist_pls2 <- NULL
  for (i in 1:2){
    Pdist_ecd2 <- c(Pdist_ecd2, subspace(square$Gamma[[i]],fit_ecd2$Gamma[[i]]))
    Pdist_pls2 <- c(Pdist_pls2, subspace(square$Gamma[[i]],fit_pls2$Gamma[[i]]))
  }
  Pdist_ecd2 <- sum(Pdist_ecd2)
  Pdist_pls2 <- sum(Pdist_pls2)
  expect_equal(c(dist_ecd2, dist_pls2, Pdist_ecd2, Pdist_pls2), c(5.80500046, 5.59056824, 1.41448471, 0.15631639))
})

test_that("Section 3.5", {
  set.seed(1)
  data("EEG")
  u_eeg <- c(1,1)
  fit_eeg_ols <- TRR.fit(EEG$x, EEG$y, method = "standard")
  fit_eeg_1D <- TRR.fit(EEG$x, EEG$y, u_eeg, method = "1D")
  testthat::expect_equal(fit_eeg_ols$coefficients@data[8:12,8,1], c(-0.366630830, -0.261239935, -0.219068713, -0.446891195, -0.052921259))
  testthat::expect_equal(fit_eeg_1D$coefficients@data[8:12,25,1], c(0.89855351, 0.91008623, 0.38449156, 0.41894440, 0.47826560))
})

test_that("Section 4.3", {
  set.seed(1)
  p <- 20
  u <- 5
  tmp <- matrix(runif(p*u), p, u)
  Gamma <- qr.Q(qr(tmp))
  Gamma0 <- qr.Q(qr(Gamma),complete=T)[,(u+1):p]
  A <- matrix(runif(u^2), u, u)
  Omega <- A%*%t(A)
  A <- matrix(runif((p-u)^2), p-u, p-u)
  Omega0 <- A%*%t(A)
  A <- matrix(runif(u^2), u, u)
  Phi <- A%*%t(A)

  # Model M1
  U <- Gamma%*%Phi%*%t(Gamma)
  M <- Gamma%*%Omega%*%t(Gamma) + Gamma0%*%Omega0%*%t(Gamma0)
  M <- M + 0.00001*diag(1,p,p)
  G1 <- EnvMU(M, U, u)
  G2 <- ECD(M, U, u)
  G3 <- manifold1D(M, U, u)
  G4 <- OptimballGBB1D(M, U, u)
  G5 <- manifoldFG(M, U, u, Gamma_init = G3)
  G6 <- OptStiefelGBB(Gamma_init = G4, opts=NULL, fun = FGfun, M, U)$Gamma
  d1 <- subspace(G5, Gamma)
  d2 <- subspace(G6, Gamma)

  # The randomly generated matrix
  A <- matrix(runif(p*u), p, u)
  G7 <- manifoldFG(M, U, u, A)
  G8 <- OptStiefelGBB(Gamma_init = A, opts=NULL, FGfun, M, U)$Gamma
  d3 <- subspace(G7, Gamma)
  d4 <- subspace(G8, Gamma)
  expect_equal(c(d1,d2,d3,d4), c(4.5243985e-10, 1.7391552e-09, 4.4721360e-01, 6.0979188e-01))
})


test_that("Section 4.4", {
  set.seed(1)
  p <- 50
  u <- 5
  # Generate M and U from Model (M1)
  tmp <- matrix(runif(p*u), p, u)
  Gamma <- qr.Q(qr(tmp))
  Gamma0 <- qr.Q(qr(Gamma),complete=T)[,(u+1):p]
  A <- matrix(runif(u^2), u, u)
  Omega <- A%*%t(A)
  A <- matrix(runif((p-u)^2), p-u, p-u)
  Omega0 <- A%*%t(A)
  A <- matrix(runif(u^2), u, u)
  Phi <- A%*%t(A)

  # Use model M1
  U <- Gamma%*%Phi%*%t(Gamma)
  M <- Gamma%*%Omega%*%t(Gamma) + Gamma0%*%Omega0%*%t(Gamma0)
  M <- M + 0.00001*diag(1,p,p)

  n0 <- c(50,200)
  u_est3 <- NULL
  for (n in n0){
    X <- MASS::mvrnorm(n, rep(0,p), M)
    Mhat <- (t(X)%*%X)/n
    X <- MASS::mvrnorm(n, rep(0,p), U)
    Uhat <- (t(X)%*%X)/n
    output <- ballGBB1D_bic(Mhat, Uhat, n, maxdim = p/2)
    u_est3 <- c(u_est3, output$u)
  }
  expect_equal(u_est3, c(9,5))
})


