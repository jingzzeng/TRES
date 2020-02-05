## TODO:
## 1. Change the date in DESCRIPTION
## 2. We use Cholesky decomposition to calculate inverse, sqrt of matrix. Use crossprod and tcrossprod to replace the matrix multiplicatio.

## Test for TPR.fit
rm(list = ls())

# The dimension of predictor
p <- c(10, 10, 10)
# The envelope dimensions u.
u <- c(1, 1, 1)
# The dimension of response
r <- 5
# The sample size
n <- 200

set.seed(1)
# Simulate the data with \code{\link{TPR_sim}}.
dat <- TPR_sim(p = p, r = r, u = u, n = n)
x <- dat$x
y <- dat$y

# Xn <- new.env()
# Xn$Xn <- x
# Xn$Yn <- y

fit_1D <- TPR.fit(x, y, u = u, method="1D")

x <- list(x = x, y = y)

cp_results_TPR <- microbenchmark(
  fit_FG <- TPR.fit(x, u = u, method="FG"),
  fit_1D <- TPR.fit(x, u = u, method="1D"),
  fit_ECD <- TPR.fit(x, u = u, method="ECD"),
  fit_PLS <- TPR.fit(x, u = u, method="PLS"),
  fit_FG_chol <- TPR2.fit(x, u = u, method="FG"),
  fit_1D_chol <- TPR2.fit(x, u = u, method="1D"),
  fit_ECD_chol <- TPR2.fit(x, u = u, method="ECD"),
  fit_PLS_chol <- TPR2.fit(x, u = u, method="PLS"),
  times = 50
)


## Test for TRR.fit

rm(list=ls())
source("../TRR2.fit.R")
# The dimension of response
r <- c(10, 10, 10)
# The envelope dimensions u.
u <- c(2, 2, 2)
# The dimension of predictor
p <- 5
# The sample size
n <- 100

set.seed(1)
#
# Simulate the data with \code{\link{TRR_sim}}.
dat <- TRR_sim(r = r, p = p, u = u, n = n)
x <- dat$x
y <- dat$y

# Xn <- new.env()
# Xn$Xn <- x
# Xn$Yn <- y

x <- list(x = x, y = y)

cp_results <- microbenchmark(
  fit_FG <- TRR.fit(x, u = u, method="FG"),
  fit_1D <- TRR.fit(x, u = u, method="1D"),
  fit_ECD <- TRR.fit(x, u = u, method="ECD"),
  fit_PLS <- TRR.fit(x, u = u, method="PLS"),
  fit_FG_chol <- TRR2.fit(x, u = u, method="FG"),
  fit_1D_chol <- TRR2.fit(x, u = u, method="1D"),
  fit_ECD_chol <- TRR2.fit(x, u = u, method="ECD"),
  fit_PLS_chol <- TRR2.fit(x, u = u, method="PLS"),
  times = 50
)
# summary(res_1D)

# ## lm
# library(alr4)
# attach(oldfaith)
# fit <- lm(Duration~Interval)

## EEG
rm(list = ls())
set.seed(1)
library(R.matlab)
eeg <- readMat("data_EEG_LiBing.mat")

x <- eeg$y.vec
y <- eeg$X.mat
ind0 <- sample(which(x == 0), size = length(which(x==0))*0.7)
ind1 <- sample(which(x == 1), size = length(which(x==1))*0.7)
ind <- sort(c(ind0, ind1))
x <- x[ind]
y <- y[,,ind]

n <- dim(y)[3]
muy <- apply(y, c(1,2), mean)
tmp1 <- lapply(1:n, function(x) muy)
tmp1 <- array(unlist(tmp1), dim(y))
y <- y - tmp1

Cmat <- matrix(0, dim(y)[1], dim(y)[2])
for (j in seq_len(dim(y)[2])){
  idx <- ((j-1)*4 + 1):(j*4)
  Cmat[idx, j] <- 0.25
}

y <- ttm(rTensor::as.tensor(y), t(Cmat), 1)
y <- aperm(y@data, perm = c(2,1,3))
y <- as.tensor(y)

u_eeg <- TensEnv_dim(x, y)
u_eeg <- c(1,1)
fit_eeg_ols <- TRR.fit(x, y, method = "standard")
fit_eeg_1D <- TRR.fit(x, y, u_eeg, method = "1D")
fit_eeg_pls <- TRR.fit(x, y, u_eeg, method = "PLS")

plot(fit_eeg_ols, xlab = "Time", ylab = "Channels", yticks = seq(64,0, length.out=5))
plot(fit_eeg_1D, xlab = "Time", ylab = "Channels", yticks = seq(64,0, length.out=5))
plot(fit_eeg_pls, xlab = "Time", ylab = "Channels", yticks = seq(64,0, length.out=5))

# lattice::wireframe(drop(fit_ols$coefficients@data), col='gray', xlab = 'Time', ylab = 'Channel', zlab = list('Voltage (mV)', rot = 90), scales=list(arrows=FALSE), zlim = c(-5,5), aspect = c(1,0.5), panel.aspect = 0.5)

# lattice::wireframe(drop(fit_FG$coefficients@data), col='gray', xlab = 'Time', ylab = 'Channel', zlab = list('Voltage (mV)', rot = 90), scales=list(arrows=FALSE), zlim = c(-5,5), aspect = c(1,0.5), panel.aspect = 0.5)

## --------------- Compare cholesky -------------------##

