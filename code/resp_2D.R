# The TRR model: estimation, p-value plot, dimension selection
# Section 3.1 + Section 3.3 + Section 3.4
library(TRES)
# Convert a matlab file to R file
library(R.matlab)
source('~/Documents/GitHub/TRES_code/R/TRR.R')
source('~/Documents/GitHub/TRES_code/R/PMSE.R')
source('~/Documents/GitHub/TRES_code/R-more/print.Tenv.R')
source('~/Documents/GitHub/TRES_code/R-more/fitted.Tenv.R')
source('~/Documents/GitHub/TRES_code/R-more/residuals.Tenv.R')
source('~/Documents/GitHub/TRES_code/R-more/predict.Tenv.R')
source('~/Documents/GitHub/TRES_code/R-more/plot.Tenv.R')
source('~/Documents/GitHub/TRES_code/R-more/summary.Tenv.R')
source('~/Documents/GitHub/TRES_code/R-more/print.summary.Tenv.R')

set.seed(1)
B1 <- readMat("~/Documents/GitHub/TRES_jss/code/B.mat")$B
r1 <- dim(B1)
m <- length(r1)
p <- 1
u1 <- c(14, 14)
n <- 20

Omega <- Omega0 <- Gamma <- Gamma0 <- Sig <- Sigsqrtm <- U <- NULL

# Generate Gamma and Sigma
for (i in 1:m) {
  if(i == 1){
    Gamma[[i]] <- svd(B1)$u[,1:u1[i]]
  }
  else if(i == 2){
    Gamma[[i]] <- svd(B1)$v[,1:u1[i]]
  }
  Gamma0[[i]] <- qr.Q(qr(Gamma[[i]]), complete=T)[, (u1[i]+1):r1[i]]
  A <- matrix(runif(u1[i]^2), u1[i], u1[i])
  Omega[[i]] <- A %*% t(A)
  A <- matrix(runif((r1[i]-u1[i])^2), (r1[i]-u1[i]), (r1[i]-u1[i]))
  Omega0[[i]] <- A %*% t(A)
  Sig[[i]] <- Gamma[[i]] %*% Omega[[i]] %*% t(Gamma[[i]]) +
    Gamma0[[i]] %*% Omega0[[i]] %*% t(Gamma0[[i]])
  Sig[[i]] <- Sig[[i]]/norm(Sig[[i]], type="F")
  Sigsqrtm[[i]] <- pracma::sqrtm(Sig[[i]])$B
}

# Convert matrix B to three-way tensor
B1 <- array(B1, c(r1, 1))
B1 <- as.tensor(B1)

# Generate data
Xn1 <- matrix(rep(c(1, 0), each=10), 1, n)/fnorm(B1)
Epsilon1 <- array(rnorm(prod(r1)*n), c(r1, n))
Epsilon1 <- as.tensor(Epsilon1)
Epsilon1 <- ttl(Epsilon1, Sigsqrtm, ms=1:m)
Yn1 <- Epsilon1 + ttm(B1, t(Xn1), m+1)

# Dimension selection
# u_est1 <-TensEnv_dim(Yn1, Xn1, bic_max = 32)

# Fitting the TRR model with different methods
fit_ols1 <- TRR(Yn1, Xn1, method = 'standard')
fit_fg1 <- TRR(Yn1, Xn1, u1, method = 'FG')
fit_1D1 <- TRR(Yn1, Xn1, u1, method = '1D')
fit_pls1 <- TRR(Yn1, Xn1, u1, method = 'PLS')
fit_ecd1 <- TRR(Yn1, Xn1, u1, method = 'ECD')

# The estimated coefficient B
B_ols1 <- coefficients(fit_ols1)
B_fg1 <- coefficients(fit_fg1)
B_1D1 <- coefficients(fit_1D1)
B_pls1 <- coefficients(fit_pls1)
B_ecd1 <- coefficients(fit_ecd1)

# Coefficient plots
par(mfrow=c(3,2),  mar = c(2.5,3.1,2.1,2.1), mgp = c(1.5, 0.5, 0))
image(-B1@data[, , 1], axes=TRUE, col = grey(seq(0, 1, length = 256)))
title('True coefficient matrix')

plot(fit_ols1, main = 'OLS')
plot(fit_fg1, main = 'FG')
plot(fit_1D1, main = '1D')
plot(fit_pls1, main = 'PLS')
plot(fit_ecd1, main = 'ECD')

# P-values
pvalue_1D1 <-  Tenv_Pval(Yn1, Xn1, B_est = B_1D1)
pvalue_pls1 <-  Tenv_Pval(Yn1, Xn1, B_est = B_pls1)

# P-value plots
par(mfrow=c(1,3),  mar = c(2.5,3.1,2.1,2.1), mgp = c(1.5, 0.5, 0))

image((pvalue_1D1$P_OLS[, , 1]>0.05), axes=TRUE, col = grey(seq(0, 1, length = 256)))
title('OLS')

image((pvalue_1D1$P_val[, , 1]>0.05), axes=TRUE, col = grey(seq(0, 1, length = 256)))
title('1D')

image((pvalue_pls1$P_val[, , 1]>0.05), axes=TRUE, col = grey(seq(0, 1, length = 256)))
title('PLS')

# Distances between estimated coefficient and the true one.
dist_ols1 <- fnorm(B_ols1 - B1)
dist_1D1 <- fnorm(B_1D1 - B1)
dist_pls1 <- fnorm(B_pls1 - B1)

# dist_1D_ecd = 0.156, which means the 1D estimator is close to the ECD estimator.
dist_1D_ecd <- fnorm(B_ecd1 - B_1D1)

# Subspace distance
Pdist_1D1 <- NULL
Pdist_pls1 <- NULL
for (i in 1:m){
  Pdist_1D1 <- c(Pdist_1D1, subspace(Gamma[[i]],fit_1D1$Gamma_hat[[i]]))
  Pdist_pls1 <- c(Pdist_pls1, subspace(Gamma[[i]],fit_pls1$Gamma_hat[[i]]))
}
Pdist_1D1 <- sum(Pdist_1D1)
Pdist_pls1 <- sum(Pdist_pls1)

c(dist_ols = dist_ols1, dist_1D = dist_1D1, dist_pls = dist_pls1)
c(Pdist_1D = Pdist_1D1, Pdist_pls = Pdist_pls1)
