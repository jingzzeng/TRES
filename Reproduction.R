# Set digits option
options(digits = 4)

# # Install the package if it is not intalled yet
# install.packages("TRES")
library(TRES)

## ----------------------------- Section 3.1 ----------------------------- ##

# The TRR model: estimation, coefficient plots, coefficients distance and subspace distance

set.seed(1)

# laod dataset "bat"
data("bat")
dim(bat$x)
bat$y
bat$coefficients
names(bat$Gamma)

# u1 <- c(14, 14)

# Fitting the TRR model with different methods
fit_ols1 <- TRR.fit(bat$x, bat$y, method = 'standard')
fit_fg1 <- TRR.fit(bat$x, bat$y, u = c(14,14), method = 'FG')
fit_1d1 <- TRR.fit(bat$x, bat$y, u = c(14,14), method = '1D')
fit_ecd1 <- TRR.fit(bat$x, bat$y, u = c(14,14), method = 'ECD')
fit_pls1 <- TRR.fit(bat$x, bat$y, u = c(14,14), method = 'PLS')

fit_1d1
coef(fit_1d1)
fitted(fit_1d1)
residuals(fit_1d1)
summary(fit_1d1)
predict(fit_1d1, bat$x)

# True coefficient plots (p-value plots are also generated)
true_plt <- t(-bat$coefficients@data[, , 1]) # switch the sign so that pattern area is highlighted
image(x = 1:nrow(true_plt), y = 1:ncol(true_plt), z=true_plt, ylim = c(ncol(true_plt), 1), col = grey(seq(0, 1, length = 256)), xlab = "", ylab="", main='True coefficient matrix', cex.main=2, cex.axis = 2)

# Coefficient plots for each estimators (p-value plots are also generated)
plot(fit_ols1, cex.main = 2, cex.axis= 2)
plot(fit_fg1, cex.main = 2, cex.axis= 2)
plot(fit_1d1, cex.main = 2, cex.axis= 2)
plot(fit_ecd1, cex.main = 2, cex.axis= 2)
plot(fit_pls1, cex.main = 2, cex.axis= 2)

# Distances between estimated coefficient and the true one.
dist_ols1 <- rTensor::fnorm(coef(fit_ols1) - bat$coefficients)
dist_fg1 <- rTensor::fnorm(coef(fit_fg1) - bat$coefficient)
dist_1d1 <- rTensor::fnorm(coef(fit_1d1) - bat$coefficient)
dist_ecd1 <- rTensor::fnorm(coef(fit_ecd1) - bat$coefficient)
dist_pls1 <- rTensor::fnorm(coef(fit_pls1) - bat$coefficient)

# The 1D estimator is close to the ECD estimator.
dist_1D_ecd <- rTensor::fnorm(coef(fit_ecd1) - coef(fit_1d1))

c(dist_ols = dist_ols1, dist_fg = dist_fg1, dist_1d = dist_1d1, dist_ecd = dist_ecd1, dist_pls = dist_pls1)

# Subspace distance
Pdist_fg1 <- NULL
Pdist_1d1 <- NULL
Pdist_ecd1 <- NULL
Pdist_pls1 <- NULL
for (i in 1:2){
  Pdist_fg1 <- c(Pdist_fg1, subspace(bat$Gamma[[i]],fit_fg1$Gamma[[i]]))
  Pdist_1d1 <- c(Pdist_1d1, subspace(bat$Gamma[[i]],fit_1d1$Gamma[[i]]))
  Pdist_ecd1 <- c(Pdist_ecd1, subspace(bat$Gamma[[i]],fit_ecd1$Gamma[[i]]))
  Pdist_pls1 <- c(Pdist_pls1, subspace(bat$Gamma[[i]],fit_pls1$Gamma[[i]]))
}
Pdist_fg1 <- sum(Pdist_fg1)
Pdist_1d1 <- sum(Pdist_1d1)
Pdist_ecd1 <- sum(Pdist_ecd1)
Pdist_pls1 <- sum(Pdist_pls1)

c(Pdist_fg = Pdist_fg1, Pdist_1d = Pdist_1d1, Pdist_ecd = Pdist_ecd1, Pdist_pls = Pdist_pls1)


## ----------------------------- Section 3.2 ----------------------------- ##

# The TPR model: estimation, coefficient plots, coefficients distance and subspace distance

set.seed(1)

# Load dataset "square"
data("square")
square$x
dim(square$y)
square$coefficients
names(square$Gamma)
# u2 <- c(2, 2)

# Fitting the TPR model with different methods
fit_ols2 <- TPR.fit(square$x, square$y, method = "standard")
fit_fg2 <- TPR.fit(square$x, square$y, u = c(2, 2), method = "FG")
fit_1d2 <- TPR.fit(square$x, square$y, u = c(2, 2), method = "1D")
fit_ecd2 <- TPR.fit(square$x, square$y, u = c(2, 2), method = "ECD")
fit_pls2 <- TPR.fit(square$x, square$y, u = c(2, 2), method = "PLS")

fit_1d2
coef(fit_1d2)
dim(fitted(fit_1d2))
dim(residuals(fit_1d2))

# True coefficient plot
true_plt <- t(-square$coefficients@data[, , 1]) # switch the sign so that pattern area is highlighted.
image(x = 1:nrow(true_plt), y = 1:ncol(true_plt), z=true_plt, ylim = c(ncol(true_plt), 1), col = grey(seq(0, 1, length = 256)), xlab = "", ylab="", main='True coefficient matrix', cex.main=2, cex.axis = 2)

# Coefficient plots for each estimators
plot(fit_ols2, cex.main = 2, cex.axis= 2)
plot(fit_fg2, cex.main = 2, cex.axis= 2)
plot(fit_1d2, cex.main = 2, cex.axis= 2)
plot(fit_ecd2, cex.main = 2, cex.axis= 2)
plot(fit_pls2, cex.main = 2, cex.axis= 2)

# Distances between estimated coefficient and the true one.
dist_ols2 <- rTensor::fnorm(coef(fit_ols2) - square$coefficients)
dist_fg2 <- rTensor::fnorm(coef(fit_fg2) - square$coefficients)
dist_1d2 <- rTensor::fnorm(coef(fit_1d2) - square$coefficients)
dist_ecd2 <- rTensor::fnorm(coef(fit_ecd2) - square$coefficients)
dist_pls2 <- rTensor::fnorm(coef(fit_pls2) - square$coefficients)

# The 1D estimator is close to the ECD estimator.
dist_1d_ecd <- rTensor::fnorm(coef(fit_1d2) - coef(fit_ecd2))

c(dist_ols = dist_ols2, dist_fg = dist_fg2, dist_1d = dist_1d2, dist_ecd = dist_ecd2, dist_pls = dist_pls2)

# Subspace distance
Pdist_fg2 <- NULL
Pdist_1d2 <- NULL
Pdist_ecd2 <- NULL
Pdist_pls2 <- NULL
for (i in 1:2){
  Pdist_fg2 <- c(Pdist_fg2, subspace(square$Gamma[[i]],fit_fg2$Gamma[[i]]))
  Pdist_1d2 <- c(Pdist_1d2, subspace(square$Gamma[[i]],fit_1d2$Gamma[[i]]))
  Pdist_ecd2 <- c(Pdist_ecd2, subspace(square$Gamma[[i]],fit_ecd2$Gamma[[i]]))
  Pdist_pls2 <- c(Pdist_pls2, subspace(square$Gamma[[i]],fit_pls2$Gamma[[i]]))
}
Pdist_fg2 <- sum(Pdist_fg2)
Pdist_1d2 <- sum(Pdist_1d2)
Pdist_ecd2 <- sum(Pdist_ecd2)
Pdist_pls2 <- sum(Pdist_pls2)

c(Pdist_fg = Pdist_fg2, Pdist_1d = Pdist_1d2, Pdist_ecd = Pdist_ecd2, Pdist_pls = Pdist_pls2)


## ----------------------------- Section 3.3 ----------------------------- ##

set.seed(1)

# Dimension selection for both TRR and TPR models
u_est1 <-TensEnv_dim(bat$x, bat$y, maxdim = 32)
u_est2 <- TensPLS_cv2d3d(square$x, square$y, maxdim = 16, nfolds = 5)

u_est1
u_est2


## ----------------------------- Section 3.4 ----------------------------- ##

# P-values for TRR estimators
summary(fit_ols1)$p_val
summary(fit_fg1)$p_val
summary(fit_1d1)$p_val
summary(fit_ecd1)$p_val
summary(fit_pls1)$p_val

# P-value plots for TRR estimators
plot(fit_ols1, cex.main = 2, cex.axis= 2, level = 0.05)
plot(fit_fg1, cex.main = 2, cex.axis= 2, level = 0.05)
plot(fit_1d1, cex.main = 2, cex.axis= 2, level = 0.05)
plot(fit_ecd1, cex.main = 2, cex.axis= 2, level = 0.05)
plot(fit_pls1, cex.main = 2, cex.axis= 2, level = 0.05)

## ----------------------------- Section 3.5 ----------------------------- ##

set.seed(1)
data("EEG")
u_eeg <- TensEnv_dim(EEG$x, EEG$y)
fit_eeg_ols <- TRR.fit(EEG$x, EEG$y, method = "standard")
fit_eeg_1D <- TRR.fit(EEG$x, EEG$y, u_eeg, method = "1D")
fit_eeg_pls <- TRR.fit(EEG$x, EEG$y, u_eeg, method = "PLS")

plot(fit_eeg_ols, xlab = "Time", ylab = "Channels", cex.main = 2, cex.axis= 2, cex.lab=1.5)
plot(fit_eeg_1D, xlab = "Time", ylab = "Channels", cex.main = 2, cex.axis= 2, cex.lab = 1.5)
plot(fit_eeg_pls, xlab = "Time", ylab = "Channels", cex.main = 2, cex.axis= 2, cex.lab = 1.5)


## ----------------------------- Section 4.3, Table 3 ----------------------------- ##

# Compare the execution time and estimation accuracy for each functions in each model

set.seed(1)
times <- 50

for (m in c("M1", "M2", "M3")){
  # Simulation for each model
  output <- lapply(seq_len(times), function(i){
    p <- 20
    u <- 5
    # Construct U and M
    tmp <- matrix(runif(p*u), p, u)
    Gamma <- qr.Q(qr(tmp))
    Gamma0 <- qr.Q(qr(Gamma),complete=T)[,(u+1):p]

    A <- matrix(runif(u^2), u, u)
    Omega <- A%*%t(A)

    A <- matrix(runif((p-u)^2), p-u, p-u)
    Omega0 <- A%*%t(A)

    A <- matrix(runif(u^2), u, u)
    Phi <- A%*%t(A)

    if(m == "M1"){
      ## Model (M1)
      U <- Gamma%*%Phi%*%t(Gamma)
      M <- Gamma%*%Omega%*%t(Gamma) + Gamma0%*%Omega0%*%t(Gamma0)
      M <- M + 0.00001*diag(1,p,p)
    }else if(m == "M2"){
      ## Model (M2)
      U <- Gamma%*%Phi%*%t(Gamma)
      M <- Gamma%*%t(Gamma) + 0.01*Gamma0%*%t(Gamma0)
    }else if(m == "M3"){
      ## Model (M3)
      U <- Gamma%*%Phi%*%t(Gamma)
      M <- 0.01*Gamma%*%t(Gamma) + Gamma0%*%t(Gamma0)
    }

    start_time <- Sys.time()
    Ghat_pls <- EnvMU(M, U, u)
    end_time <- Sys.time()
    exe_time_1 <- difftime(end_time, start_time, units = 'secs')
    dist_1 <- subspace(Ghat_pls, Gamma)

    start_time <- Sys.time()
    Ghat_ecd <- ECD(M, U, u)
    end_time <- Sys.time()
    exe_time_2 <- difftime(end_time, start_time, units = 'secs')
    dist_2 <- subspace(Ghat_ecd, Gamma)

    start_time <- Sys.time()
    Ghat_mani1D <- manifold1D(M, U, u)
    end_time <- Sys.time()
    exe_time_3 <- difftime(end_time, start_time, units = 'secs')
    dist_3 <- subspace(Ghat_mani1D, Gamma)

    start_time <- Sys.time()
    Ghat_feasi1D <- OptimballGBB1D(M, U, u)
    end_time <- Sys.time()
    exe_time_4 <- difftime(end_time, start_time, units = 'secs')
    dist_4 <- subspace(Ghat_feasi1D, Gamma)

    start_time <- Sys.time()
    Ghat_maniFG <- manifoldFG(M, U, u, Ghat_mani1D)
    end_time <- Sys.time()
    exe_time_5 <- difftime(end_time, start_time, units = 'secs')
    dist_5 <- subspace(Ghat_maniFG, Gamma)

    start_time <- Sys.time()
    Ghat_feasiFG <- OptStiefelGBB(Ghat_feasi1D, opts=NULL, FGfun, M, U)$Gamma
    end_time <- Sys.time()
    exe_time_6 <- difftime(end_time, start_time, units = 'secs')
    dist_6 <- subspace(Ghat_feasiFG, Gamma)

    list(c(exe_time_1, exe_time_2, exe_time_3, exe_time_4, exe_time_5, exe_time_6),
        c(dist_1, dist_2, dist_3, dist_4, dist_5, dist_6))

  })

  exe_time <- do.call(rbind, lapply(output, "[[", 1))
  dist <-  do.call(rbind, lapply(output, "[[", 2))

  # Average execution time and standard error for each method
  mean_time <- apply(exe_time, 2, mean)
  se_time <- apply(exe_time, 2, sd)/sqrt(times)

  # Average subspace distance and standard error for each method
  mean_dist <- apply(dist, 2, mean)
  se_dist <- apply(dist, 2, sd)/sqrt(times)

  cat("--------------------------------------------------------------------------------\n")
  cat("Model:", m, "\n")
  cat("Averaged execution time (standard error)\n")
  tmp <- paste0(format(mean_time, digits = 2), '(', format(se_time, scientific = TRUE), ')')
  names(tmp) <-  c('PLS', 'ECD', '1D_Mani', '1D_Feasi', 'FG_Mani', 'FG_Feasi')
  print(tmp, quote=FALSE, print.gap=2L)
  cat("\n--------------------------------------------------------------------------------\n")

  cat("--------------------------------------------------------------------------------\n")
  cat("Model:", m, "\n")
  cat("Estimation accuracy (standard error)\n")
  tmp <- paste0(format(mean_dist, scientific = 2), '(', format(se_dist, scientific = TRUE), ')')
  names(tmp) <-  c('PLS', 'ECD', '1D_Mani', '1D_Feasi', 'FG_Mani', 'FG_Feasi')
  print(tmp, quote=FALSE, print.gap=2L)
  cat("\n--------------------------------------------------------------------------------\n")

}


## ----------------------------- Section 4.3 ----------------------------- ##

set.seed(1)

# Compare the performance of the FG algorithm with different initial values: 1D estimator and randomly generated matrix.
p <- 20
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

c(d1,d2,d3,d4)


## ----------------------------- Section 4.4 ----------------------------- ##

set.seed(1)

# Dimension selection with different sample size
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

n0 <- c(50,70,100, 200, 400, 800)
u_est3 <- NULL
for (n in n0){
  # Generate Wishart sample Mhat and Uhat
  X <- MASS::mvrnorm(n, rep(0,p), M)
  Mhat <- (t(X)%*%X)/n
  X <- MASS::mvrnorm(n, rep(0,p), U)
  Uhat <- (t(X)%*%X)/n
  output <- ballGBB1D_bic(Mhat, Uhat, n, maxdim = p/2)
  u_est3 <- c(u_est3, output$u)
}

u_est3
