\name{Tenv_Pval}
\alias{Tenv_Pval}

\title{
The \eqn{p}-values of each elements in the tensor response regression coefficient estimator
}
\description{
Obtain \eqn{p}-values of each elements in the tensor regression coefficient estimator. Two-sided t-tests are applied on the coefficent estimator, where the estimator use the asymptotic covariance of the OLS estimator.
}
\usage{
Tenv_Pval(Xn, Yn, B_est)
}

\arguments{
  \item{Xn}{The response tensor instance \eqn{ r_1\times r_2\times \cdots \times r_m}.}
  \item{Yn}{A vector predictor of dimension \eqn{p}.}
  \item{B_est}{
The tensor regression coefficient estimator for tensor response regression.
}
}

\value{
\item{P_OLS}{The p-value tensor of OLS estimator.}
\item{P_val}{The p-value tensor of \code{B_est}.}
\item{se}{The standard error tensor of \code{B_est}.}
}


\examples{
rm(list=ls())
r <- c(10, 10, 10)
m <- length(r)
u <- c(2, 2, 2)
p <- 5
n <- 100
set.seed(1)
eta <- array(runif(prod(u)*p), c(u, p))
eta <- rTensor::as.tensor(eta)
Gamma <- Gamma0 <- Omega <- Omega0 <- Sig <- Sigsqrtm <- NULL
for (i in 1:m){
  tmp <- matrix(runif(r[i]*u[i]), r[i], u[i])
  Gamma[[i]] <- qr.Q(qr(tmp))
  Gamma0[[i]] <- qr.Q(qr(Gamma[[i]]),complete=TRUE)[,(u[i]+1):r[i]]
  A <- matrix(runif(u[i]^2), u[i], u[i])
  Omega[[i]] <- A \%*\% t(A)
  A <- matrix(runif((r[i]-u[i])^2), (r[i]-u[i]), (r[i]-u[i]))
  Omega0[[i]] <- A \%*\% t(A)
  Sig[[i]] <- Gamma[[i]] \%*\% Omega[[i]] \%*\% t(Gamma[[i]])+
    Gamma0[[i]] \%*\% Omega0[[i]] \%*\% t(Gamma0[[i]])
  Sig[[i]] <- 10*Sig[[i]]/norm(Sig[[i]], type="F")+0.01*diag(r[i])
  ttt
  Sigsqrtm[[i]] <- pracma::sqrtm(Sig[[i]])$B
}
B <- rTensor::ttl(eta, Gamma, ms=1:m)
Xn <- matrix(rnorm(p*n), p, n)
Xn_inv <- MASS::ginv(Xn \%*\% t(Xn)) \%*\% Xn
Epsilon <- array(rnorm(prod(r)*n), c(r, n))
Epsilon <- rTensor::as.tensor(Epsilon)
Epsilon <- rTensor::ttl(Epsilon, Sigsqrtm, ms=1:m)
Yn <- Epsilon + rTensor::ttm(B, t(Xn), m+1)

res_1D <- TRR(Xn, Yn, u, method="1D")
Tenv_Pval(Xn, Yn, res_1D$coefficients)
}
