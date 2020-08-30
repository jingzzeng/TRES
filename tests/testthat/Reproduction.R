## -------------------------------------------- Note ------------------------------------------- ##
# The results in paper "TRES: An R Package for Tensor Regression and Envelope Algorithms"
# can be reproduced by this file on macOS High Sierra 10.13.6 and Windows 10. The results 
# on Linux using R 3.6.x might be slightly different due to the inconsistency of random seed.
## --------------------------------------------------------------------------------------------- ##

# Set digits option
options(digits = 3)

# # Install the package if it is not intalled yet
# install.packages("TRES")
library(TRES)

## Set up RNG kind
RNGkind("L'Ecuyer-CMRG")
## ----------------------------- Section 3.1 ----------------------------- ##

# The TRR model: estimation, coefficient plots, coefficients distance and subspace distance
set.seed(1)
# load dataset "bat"
data("bat")
str(bat)

# Fitting the TRR model with different methods
fit_ols1 <- TRR.fit(bat$x, bat$y, method = "standard")
fit_1d1 <- TRR.fit(bat$x, bat$y, u = c(14,14), method = "1D")
fit_pls1 <- TRR.fit(bat$x, bat$y, u = c(14,14), method = "PLS")

fit_1d1
coef(fit_1d1)
fitted(fit_1d1)
residuals(fit_1d1)
summary(fit_1d1)
predict(fit_1d1, bat$x)

# True coefficient plots (p-value plots are also generated)
true_B <- bat$coefficients@data[, , 1] # switch the sign so that pattern area is highlighted
image(x = 1:nrow(true_B), y = 1:ncol(true_B), z=-t(true_B), ylim = c(ncol(true_B), 1), col = grey(seq(0, 1, length = 256)), xlab = "", ylab="", main='True coefficient matrix', cex.main=2, cex.axis = 2)
graphics::box()

# Coefficient plots for each estimators (p-value plots are also generated)
plot(fit_ols1, cex.main = 2, cex.axis= 2)
plot(fit_1d1, cex.main = 2, cex.axis= 2)
plot(fit_pls1, cex.main = 2, cex.axis= 2)

# Distances between estimated coefficient and the true one.
dist_ols1 <- rTensor::fnorm(coef(fit_ols1) - bat$coefficients)
dist_1d1 <- rTensor::fnorm(coef(fit_1d1) - bat$coefficient)
dist_pls1 <- rTensor::fnorm(coef(fit_pls1) - bat$coefficient)

c(dist_ols1, dist_1d1, dist_pls1)

# Subspace distance
Pdist_1d1 <- rep(NA_real_, 2)
Pdist_pls1 <- rep(NA_real_, 2)
for (i in 1:2){
  Pdist_1d1[i] <- subspace(bat$Gamma[[i]],fit_1d1$Gamma[[i]])
  Pdist_pls1[i] <- subspace(bat$Gamma[[i]],fit_pls1$Gamma[[i]])
}
Pdist_1d1 <- sum(Pdist_1d1)
Pdist_pls1 <- sum(Pdist_pls1)

c(Pdist_1d1, Pdist_pls1)


## ----------------------------- Section 3.2 ----------------------------- ##

# The TPR model: estimation, coefficient plots, coefficients distance and subspace distance
set.seed(1)
# Load dataset "square"
data("square")
str(square)

# Fitting the TPR model with different methods
fit_ols2 <- TPR.fit(square$x, square$y, method = "standard")
fit_1d2 <- TPR.fit(square$x, square$y, u = c(2, 2), method = "1D")
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
plot(fit_1d2, cex.main = 2, cex.axis= 2)
plot(fit_pls2, cex.main = 2, cex.axis= 2)

# Distances between estimated coefficient and the true one.
dist_ols2 <- rTensor::fnorm(coef(fit_ols2) - square$coefficients)
dist_1d2 <- rTensor::fnorm(coef(fit_1d2) - square$coefficients)
dist_pls2 <- rTensor::fnorm(coef(fit_pls2) - square$coefficients)

c(dist_ols2, dist_1d2, dist_pls2)

# Subspace distance
Pdist_1d2 <- rep(NA_real_, 2)
Pdist_pls2 <- rep(NA_real_, 2)
for (i in 1:2){
  Pdist_1d2[i] <- subspace(square$Gamma[[i]],fit_1d2$Gamma[[i]])
  Pdist_pls2[i] <- subspace(square$Gamma[[i]],fit_pls2$Gamma[[i]])
}
Pdist_1d2 <- sum(Pdist_1d2)
Pdist_pls2 <- sum(Pdist_pls2)

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
summary(fit_1d1)$p_val
summary(fit_pls1)$p_val

# P-value plots for TRR estimators
plot(fit_ols1, cex.main = 2, cex.axis= 2, level = 0.05)
plot(fit_1d1, cex.main = 2, cex.axis= 2, level = 0.05)
plot(fit_pls1, cex.main = 2, cex.axis= 2, level = 0.05)

## ----------------------------- Section 3.5 ----------------------------- ##
set.seed(1)
library(ggplot2) # For density plots

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

# The material part of y
Gamma <- lapply(fit_eeg_1D$Gamma, t)
material <- rTensor::ttl(EEG$y, Gamma, 1:2)
material <- drop(material@data)
material <- material/sd(material)
data_m <- data.frame(data = material, class = as.factor(EEG$x))
g1 <- ggplot(data_m, aes(x = data))+
  geom_density(aes(linetype = class))+
  labs(title = "Material information in response")+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        title = element_text(size = 18),
        legend.position = "none",
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))
print(g1)

# The immaterial part of y.
Gamma0 <- lapply(fit_eeg_1D$Gamma, function(x){
  i <- sample(2:64, size = 1)
  t(qr.Q(qr(x),complete=TRUE)[, i, drop = FALSE])
})
immaterial <- rTensor::ttl(EEG$y, Gamma0, 1:2)
immaterial <- drop(immaterial@data)
immaterial <- immaterial/sd(immaterial)
data_im <- data.frame(data = immaterial, class = as.factor(EEG$x))
g2 <- ggplot(data_im, aes(x = data))+
  geom_density(aes(linetype = class))+
  labs(title = "Immaterial information in response") +
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        title = element_text(size = 18),
        legend.position = "none",
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))
print(g2)
## ----------------------------- Section 4.3, Table 3 ----------------------------- ##

# Compare the execution time and estimation accuracy for each functions in each model

set.seed(1)
times <- 50

## Calculate the standard error of median based on 1000 boostrap samples.
bootse <- function(x){
  result <- sapply(seq_len(1000), function(i){median(sample(x, replace = TRUE))})
  sd(result)
}

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
    Ghat_OptM1D <- OptM1D(M, U, u)
    end_time <- Sys.time()
    exe_time_4 <- difftime(end_time, start_time, units = 'secs')
    dist_4 <- subspace(Ghat_OptM1D, Gamma)

    start_time <- Sys.time()
    Ghat_maniFG <- manifoldFG(M, U, u)
    end_time <- Sys.time()
    exe_time_5 <- difftime(end_time, start_time, units = 'secs')
    dist_5 <- subspace(Ghat_maniFG, Gamma)

    start_time <- Sys.time()
    Ghat_OptMFG <- OptMFG(M, U, u)
    end_time <- Sys.time()
    exe_time_6 <- difftime(end_time, start_time, units = 'secs')
    dist_6 <- subspace(Ghat_OptMFG, Gamma)

    list(c(exe_time_1, exe_time_2, exe_time_3, exe_time_4, exe_time_5, exe_time_6),
        c(dist_1, dist_2, dist_3, dist_4, dist_5, dist_6))

  })

  exe_time <- do.call(rbind, lapply(output, "[[", 1))
  dist <-  do.call(rbind, lapply(output, "[[", 2))

  # Average execution time
  median_time <- apply(exe_time, 2, median)
  ## Standard error based on 1000 bootstrap samples
  se_time <- apply(exe_time, 2, bootse)

  # Average subspace distance
  median_dist <- apply(dist, 2, median)
  ## Standard error based on 1000 bootstrap samples
  se_dist <- apply(dist, 2, bootse)

  cat("--------------------------------------------------------------------------------\n")
  cat("Model:", m, "\n")
  cat("Median execution time (standard error)\n")
  tmp <- paste0(format(median_time, digits = 2), '(', format(se_time, scientific = TRUE), ')')
  names(tmp) <-  c('PLS', 'ECD', '1D_Mani', '1D_OptM', 'FG_Mani', 'FG_OptM')
  print(tmp, quote=FALSE, print.gap=2L)
  cat("\n--------------------------------------------------------------------------------\n")

  cat("--------------------------------------------------------------------------------\n")
  cat("Model:", m, "\n")
  cat("Estimation accuracy (standard error)\n")
  tmp <- paste0(format(median_dist, scientific = 2), '(', format(se_dist, scientific = TRUE), ')')
  names(tmp) <-  c('PLS', 'ECD', '1D_Mani', '1D_OptM', 'FG_Mani', 'FG_OptM')
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
