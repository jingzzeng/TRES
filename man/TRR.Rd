\name{TRR}
\alias{TRR}

\title{
Tensor response regression
}
\description{
This function is used for estimation of tensor response regression. The avaiable method including standard OLS type estimation, PLS type of estimation as well as envelope estimation with 1D and ECD approaches.
}
\usage{
TRR(Yn, Xn, u, method)
}

\arguments{
  \item{Yn}{The response tensor instance or dimension \eqn{r_1\times r_2\times\cdots\times r_m \times n}, where \eqn{n} is the sample size.}
  \item{Xn}{The predictor matrix of dimension \eqn{p \times n}.}
  \item{u}{The dimension of envelope subspace. \eqn{u=(u_1,\cdots, u_m)}.}
  \item{method}{The method used for estimation of tensor response regression. There are four possible choices.}
  \itemize{
\item{\code{"standard"}}: The standard OLS type estimation.
\item{\code{"1D"}}: Envelope estimation with one dimenaional optimization approaches by 1D algorithm.
\item{\code{"ECD"}}: Envelope estimation with one dimenaional optimization approaches by ECD algorithm.
\item{\code{"PLS"}}: The SIMPLS-type estimation without mannifold optimization.
}

}

\value{
\item{Bhat}{The estimation of regression coefficient tensor.}
\item{Gamma_hat}{The estimation of envelope subspace basis.}
\item{Sig}{A matrix lists of \eqn{\boldsymbol{\Sigma}_k}, where \eqn{\boldsymbol{\Sigma}=\boldsymbol{\Sigma}_m \otimes \cdots \otimes \boldsymbol{\Sigma}_1}.}
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
eta <- as.tensor(eta)

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
  Sigsqrtm[[i]] <- pracma::sqrtm(Sig[[i]])$B
}
B <- ttl(eta, Gamma, ms=1:m)
Xn <- matrix(rnorm(p*n), p, n)
Xn_inv <- MASS::ginv(Xn \%*\% t(Xn)) \%*\% Xn
Epsilon <- array(rnorm(prod(r)*n), c(r, n))
Epsilon <- as.tensor(Epsilon)
Epsilon <- ttl(Epsilon, Sigsqrtm, ms=1:m)
Yn <- Epsilon + ttm(B, t(Xn), m+1)

res_ECD = TRR(Yn, Xn, u, method="ECD")
res_1D = TRR(Yn, Xn, u, method="1D")
res_pls = TRR(Yn, Xn, u, method="PLS")
res_std = TRR(Yn, Xn, u, method="standard")

rTensor::fnorm(B-res_ECD$Bhat)
rTensor::fnorm(B-res_1D$Bhat)
rTensor::fnorm(B-res_pls$Bhat)
rTensor::fnorm(B-res_std$Bhat)
}
