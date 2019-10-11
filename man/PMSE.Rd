\name{PMSE}
\alias{PMSE}

\title{
Predictions and Mean squared error for tensor predictor regression
}
\description{
Evaluate tensor predictor regression through prediction mean squared error.
}
\usage{
PMSE(Xn, Yn, Bhat)
}

\arguments{
  \item{Xn}{A predictor tensor of dimension \eqn{p_1\times p_2\times\cdots\times p_m \times n}, where \eqn{n} is the sample size. Matrix and array are also allowed.}
  \item{Yn}{A response matrix of dimension \eqn{r \times n}, where \eqn{n} is the sample size. Vector is also allowed.}
  \item{Bhat}{An estimation of coefficient tensor of dimension \eqn{p_1\times p_2\times\cdots\times p_m \times r}. Vector, matrix and array are also allowed as long as the dimensions match the ones of \code{Xn} and \code{Yn}.}
}

\value{
\item{mse}{Mean squared error. Defined as \eqn{1/n\sum_{i=1}^n(\|\mathbf{Y}_i-\hat{\mathbf{Y}}_i)\|_F}, where \eqn{\hat{\mathbf{Y}}_i} is the predictions and \eqn{\|\cdot\|_F} is the Frobenius norm of tensor.}
\item{Yhat}{The predictions of tensor predictor regression.}
}


\examples{

rm(list = ls())

p <- c(10, 10, 10)
u <- c(1, 1, 1)
m <- 3; r <- 5; n <- 200
eta <- array(runif(prod(u,r)), c(u,r))
eta <- rTensor::as.tensor(eta)

Gamma <- Gamma0 <- Omega <- Omega0 <- Sig <- Sigsqrtm <- NULL
for(i in 1:m) {
  tmp <- matrix(runif(p[i]*u[i]), p[i], u[i])
  Gamma[[i]] <- qr.Q(qr(tmp))
  Gamma0[[i]] <- qr.Q(qr(tmp), complete=TRUE)[, (u[i]+1):p[i]]
  Omega[[i]] <- diag(u[i])
  Omega0[[i]] <- 0.01*diag(p[i]-u[i])
  Sig[[i]] <- Gamma[[i]] \%*\% Omega[[i]] \%*\% t(Gamma[[i]])+
    Gamma0[[i]] \%*\% Omega0[[i]] \%*\% t(Gamma0[[i]])
  Sig[[i]] <- 2*Sig[[i]]/norm(Sig[[i]], type="F")
  Sigsqrtm[[i]] <- pracma::sqrtm(Sig[[i]])$B
}

B <- rTensor::ttl(eta,Gamma, ms = c(1:m))
A <- matrix(runif(r^2), r, r)
SigY <- A \%*\% t(A)
SigY <- SigY/norm(SigY, type="F")


##generate data
Epsilon <- MASS::mvrnorm(n, mu=rep(0, r), Sigma=SigY)
tmp2 <- array(rnorm(prod(p, n)), c(p, n))
Xn <- rTensor::as.tensor(tmp2)
Xn <- rTensor::ttl(Xn, Sigsqrtm, ms = c(1:m))
vecXn <- matrix(Xn@data, prod(p), n)
Y_tmp <- matrix(NA, r, 200)
tmp <- array(NA, c(p, r))
for (j in 1:n) {
  for (s in 1:r) {
    tmp[, , , s] <-  B@data[, , , s]*Xn@data[, , , j]
  }
  Y_tmp[, j] <- apply(tmp, 4, sum)
}
Yn <- Y_tmp + t(Epsilon)

res_std <- TPR(Xn, Yn, u, method="standard")
PMSE(Xn, Yn, res_std$coefficients)

}
