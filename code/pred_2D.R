# Section 3.2 + Section 3.3
library(TRES)
# Convert a matlab file to R file
library(R.matlab)
source('~/Documents/GitHub/TRES_code/R/TPR.R')
source('~/Documents/GitHub/TRES_code/R/PMSE.R')
source('~/Documents/GitHub/TRES_code/R-more/print.Tenv.R')
source('~/Documents/GitHub/TRES_code/R-more/fitted.Tenv.R')
source('~/Documents/GitHub/TRES_code/R-more/residuals.Tenv.R')
source('~/Documents/GitHub/TRES_code/R-more/predict.Tenv.R')
source('~/Documents/GitHub/TRES_code/R-more/plot.Tenv.R')


set.seed(1)
B2 <- readMat("~/Documents/GitHub/TRES_jss/code/prB.mat")$B
r <- 1
m <- 2
p2 <- dim(B2)
u2 <- c(2, 2)
n <- 200

# Generate Gamma and Sigma
Omega <- Omega0 <- Gamma <- Gamma0 <- Sig <- Sigsqrtm <- U <- NULL
for (i in 1:m) {
  if(i == 1){
    Gamma[[i]] <- svd(B2)$u[,1:u2[i]]
  }
  else if(i == 2){
    Gamma[[i]] <- svd(B2)$v[,1:u2[i]]
  }
  Gamma0[[i]] <- qr.Q(qr(Gamma[[i]]), complete=T)[, (u2[i]+1):p2[i]]
  Omega[[i]] <- diag(u2[i])
  Omega0[[i]] <- 0.01*diag(p2[i]-u2[i])
  Sig[[i]] <- Gamma[[i]] %*% Omega[[i]] %*% t(Gamma[[i]]) +
    Gamma0[[i]] %*% Omega0[[i]] %*% t(Gamma0[[i]])
  Sig[[i]] <- Sig[[i]]/norm(Sig[[i]], type = "F")
  Sigsqrtm[[i]] <- pracma::sqrtm(Sig[[i]])$B
}

# Convert matrix B to three-way tensor
B2 <- array(B2, c(p2, 1))
B2 <- as.tensor(B2)

# Generate data
Epsilon2 <- matrix(rnorm(r*n),r,n)
Xn2 <- array(rnorm(prod(p2)*n), c(p2, n))
Xn2 <- as.tensor(Xn2)
Xn2 <-  ttl(Xn2, Sigsqrtm, ms = c(1:m))

Y_tmp <- matrix(NA, 1, 200)
for (j in 1:n) {
  tmp <-  B2@data[, , 1]*Xn2@data[, , j]
  Y_tmp[, j] <- sum(tmp)
}
Yn2 <-  Y_tmp + Epsilon2

# dimension selection
# This function is time-consuming. The execution time is around 20 mins.
# u_est2 <- TensPLS_cv2d3d(Xn2, Yn2, maxdim = 32, nfolds = 5)

# Fitting the TPR model with different methods
fit_ols2 <- TPR(Yn2, Xn2, method = 'standard')
fit_fg2 <- TPR(Yn2, Xn2, u2, method = 'FG')
fit_1D2 <- TPR(Yn2, Xn2, u2, method = '1D')
fit_pls2 <- TPR(Yn2, Xn2, u2, method = 'PLS')
fit_ecd2 <- TPR(Yn2, Xn2, u2, method = 'ECD')

# The estimated coefficient B
B_ols2 <- coefficients(fit_ols2)
B_fg2 <- coefficients(fit_fg2)
B_1D2 <- coefficients(fit_1D2)
B_pls2 <- coefficients(fit_pls2)
B_ecd2 <- coefficients(fit_ecd2)

# Coefficient plot
par(mfrow=c(2,2),  mar = c(2.5,3.1,2.1,2.1), mgp = c(1.5, 0.5, 0))
image(-B2@data[, , 1], axes=TRUE, col = grey(seq(0, 1, length = 256)))
title('True coefficient matrix')

plot(fit_ols2, main = 'OLS')
plot(fit_fg2, main = 'FG')
plot(fit_1D2, main = '1D')
plot(fit_pls2, main = 'PLS')
plot(fit_ecd2, main = 'ECD')

# Distances between estimated coefficient and the true one.
dist_ols2 <- fnorm(B_ols2 - B2)
dist_1D2 <- fnorm(B_1D2 - B2)
dist_pls2 <- fnorm(B_pls2 - B2)

# dist_1D_ecd = 0.000376, which means the 1D estimator is close to the ECD estimator.
dist_1D_ecd <- fnorm(B_1D2 - B_ecd2)

# Subspace distance
Pdist_1D2 <- NULL
Pdist_pls2 <- NULL
for (i in 1:m){
  Pdist_1D2 <- c(Pdist_1D2, subspace(Gamma[[i]],fit_1D2$Gamma_hat[[i]]))
  Pdist_pls2 <- c(Pdist_pls2, subspace(Gamma[[i]],fit_pls2$Gamma_hat[[i]]))
}
Pdist_1D2 <- sum(Pdist_1D2)
Pdist_pls2 <- sum(Pdist_pls2)

c(dist_ols = dist_ols2, dist_1D = dist_1D2, dist_pls = dist_pls2)
c(Pdist_1D = Pdist_1D2, Pdist_pls = Pdist_pls2)
