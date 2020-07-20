# Set digits option
options(digits = 3)

# # Install the package if it is not intalled yet
# install.packages("TRES")
library(TRES)
## ----------------------------- Section 3.1 ----------------------------- ##

######### REMOVE
# set.seed(1)
# data("bat")
# str(bat)
# # FG and weighted-FG
# fit_fg <- TRR.fit(bat$x, bat$y, u = c(14,14), method = 'FG')
# fit_wtfg <- TRR.fit(bat$x, bat$y, u = c(14,14), method = 'wtFG')
#
# # 1D and weighted-1D
# fit_wt1D <- TRR.fit(bat$x, bat$y, u = c(14,14), method = 'wt1D')
# fit_1D <- TRR.fit(bat$x, bat$y, u = c(14,14), method = '1D')
#
# # plots
# plot(fit_fg)
# plot(fit_wtfg)
# plot(fit_1D)
# plot(fit_wt1D)
#
# # estimation error for coefficients
# dist_fg <- rTensor::fnorm(coef(fit_fg) - bat$coefficient)
# dist_wtfg <- rTensor::fnorm(coef(fit_wtfg) - bat$coefficient)
# dist_1D <- rTensor::fnorm(coef(fit_1D) - bat$coefficient)
# dist_wt1D <- rTensor::fnorm(coef(fit_wt1D) - bat$coefficient)
#
# set.seed(1)
# data("square")
# str(square)

# Fitting the TPR model with different methods
# fit_ols2 <- TPR.fit(square$x, square$y, method = "standard")
# fit_fg2 <- TPR.fit(square$x, square$y, u = c(2, 2), method = "FG")
# fit_1d2 <- TPR.fit(square$x, square$y, u = c(2, 2), method = "1D")
############

# The TRR model: estimation, coefficient plots, coefficients distance and subspace distance
set.seed(1)

# laod dataset "bat"
data("bat")
str(bat)

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
true_B <- bat$coefficients@data[, , 1] # switch the sign so that pattern area is highlighted
image(x = 1:nrow(true_B), y = 1:ncol(true_B), z=-t(true_B), ylim = c(ncol(true_B), 1), col = grey(seq(0, 1, length = 256)), xlab = "", ylab="", main='True coefficient matrix', cex.main=2, cex.axis = 2)
box()

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

# c(dist_ols = dist_ols1, dist_fg = dist_fg1, dist_1d = dist_1d1, dist_ecd = dist_ecd1, dist_pls = dist_pls1)
c(dist_ols1, dist_1d1, dist_pls1)

# Subspace distance
Pdist_fg1 <- rep(NA_real_, 2)
Pdist_1d1 <- rep(NA_real_, 2)
Pdist_ecd1 <- rep(NA_real_, 2)
Pdist_pls1 <- rep(NA_real_, 2)
for (i in 1:2){
  Pdist_fg1[i] <- subspace(bat$Gamma[[i]],fit_fg1$Gamma[[i]])
  Pdist_1d1[i] <- subspace(bat$Gamma[[i]],fit_1d1$Gamma[[i]])
  Pdist_ecd1[i] <- subspace(bat$Gamma[[i]],fit_ecd1$Gamma[[i]])
  Pdist_pls1[i] <- subspace(bat$Gamma[[i]],fit_pls1$Gamma[[i]])
}
Pdist_fg1 <- sum(Pdist_fg1)
Pdist_1d1 <- sum(Pdist_1d1)
Pdist_ecd1 <- sum(Pdist_ecd1)
Pdist_pls1 <- sum(Pdist_pls1)

# c(Pdist_fg = Pdist_fg1, Pdist_1d = Pdist_1d1, Pdist_ecd = Pdist_ecd1, Pdist_pls = Pdist_pls1)
c(Pdist_1d1, Pdist_pls1)


## ----------------------------- Section 3.2 ----------------------------- ##

# The TPR model: estimation, coefficient plots, coefficients distance and subspace distance
set.seed(1)
# Load dataset "square"
data("square")
str(square)

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
true_B <- square$coefficients@data[, , 1] # switch the sign so that pattern area is highlighted
image(x = 1:nrow(true_B), y = 1:ncol(true_B), z=-t(true_B), ylim = c(ncol(true_B), 1), col = grey(seq(0, 1, length = 256)), xlab = "", ylab="", main='True coefficient matrix', cex.main=2, cex.axis = 2)
box()

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

# c(dist_ols = dist_ols2, dist_fg = dist_fg2, dist_1d = dist_1d2, dist_ecd = dist_ecd2, dist_pls = dist_pls2)
c(dist_ols2, dist_1d2, dist_pls2)

# Subspace distance
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

# c(Pdist_fg = Pdist_fg2, Pdist_1d = Pdist_1d2, Pdist_ecd = Pdist_ecd2, Pdist_pls = Pdist_pls2)
c(Pdist_1d2, Pdist_pls2)


## ----------------------------- Section 3.3 ----------------------------- ##
set.seed(1)

# Dimension selection for both TRR and TPR models
uhat1 <- TRRdim(bat$x, bat$y, maxdim = 32)
uhat1

uhat2 <- TPRdim(square$x, square$y, maxdim = 16)
uhat2
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
str(EEG)
u_eeg <- TRRdim(EEG$x, EEG$y)
u_eeg

fit_eeg_ols <- TRR.fit(EEG$x, EEG$y, method = "standard")
fit_eeg_1D <- TRR.fit(EEG$x, EEG$y, u_eeg$u, method = "1D")
fit_eeg_pls <- TRR.fit(EEG$x, EEG$y, u_eeg$u, method = "PLS")

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
    if(m == "M1"){
      data <- MenvU_sim(p, u, jitter = 1e-5)
      Gamma <- data$Gamma
      M <- data$M
      U <- data$U
    }else if(m == "M2"){
      Omega <- diag(1, u, u)
      Omega0 <- diag(0.01, p-u, p-u)
      data <- MenvU_sim(p, u, Omega = Omega, Omega0 =  Omega0)
      Gamma <- data$Gamma
      M <- data$M
      U <- data$U
    }else if(m == "M3"){
      Omega <- diag(0.01, u, u)
      Omega0 <- diag(1, p-u, p-u)
      data <- MenvU_sim(p, u, Omega = Omega, Omega0 =  Omega0)
      Gamma <- data$Gamma
      M <- data$M
      U <- data$U
    }

    start_time <- Sys.time()
    Ghat_pls <- simplsMU(M, U, u)
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
    Ghat_feasi1D <- OptM1D(M, U, u)
    end_time <- Sys.time()
    exe_time_4 <- difftime(end_time, start_time, units = 'secs')
    dist_4 <- subspace(Ghat_feasi1D, Gamma)

    start_time <- Sys.time()
    Ghat_maniFG <- manifoldFG(M, U, u)
    end_time <- Sys.time()
    exe_time_5 <- difftime(end_time, start_time, units = 'secs')
    dist_5 <- subspace(Ghat_maniFG, Gamma)

    start_time <- Sys.time()
    Ghat_feasiFG <- OptMFG(M, U, u)
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

# Generate Gamma, M and U from Model (M1)
data <- MenvU_sim(p, u, jitter = 1e-5)
Gamma <- data$Gamma
M <- data$M
U <- data$U

G <- vector("list", 8)
G[[1]] <- simplsMU(M, U, u)
G[[2]] <- ECD(M, U, u)
G[[3]] <- manifold1D(M, U, u)
G[[4]] <- OptM1D(M, U, u)
G[[5]] <- manifoldFG(M, U, u)
G[[6]] <- OptMFG(M, U, u)

d <- rep(NA_real_, 8)
for (i in 1:6){
  d[i] <- subspace(G[[i]], Gamma)
}
d[1:6]

# The randomly generated matrix
A <- matrix(runif(p*u), p, u)
G[[7]] <- manifoldFG(M, U, u, Gamma_init = A)
G[[8]] <- OptMFG(M, U, Gamma_init = A)
for (i in 7:8){
  d[i] <- subspace(G[[i]], Gamma)
}
d[5:8]


## ----------------------------- Section 4.4 ----------------------------- ##
set.seed(1)
p <- 50
u <- 5
n0 <- c(50, 70, 100, 200, 400, 800)
uhat3 <- rep(NA_integer_, length(n0))
for (i in seq_along(n0)){
  n <- n0[i]
  data <- MenvU_sim(p, u, jitter = 1e-5, wishart = TRUE, n = n)
  M <- data$M
  U <- data$U
  output <- oneD_bic(M, U, n, maxdim = p/2)
  uhat3[i] <- output$u
}
uhat3
