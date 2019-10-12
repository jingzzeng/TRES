\name{TRR}
\alias{TRR}

\title{
Tensor response regression
}
\description{
This function is used for estimation of tensor response regression. The avaiable method including standard OLS type estimation, PLS type of estimation as well as envelope estimation with FG, 1D and ECD approaches.
}
\usage{
TRR(Xn, Yn, method=c('standard', 'FG', '1D', 'ECD', 'PLS'),
    u=NULL, Gamma_init=NULL)
}

\arguments{
  \item{Xn}{The predictor matrix of dimension \eqn{p \times n}. Vector of length \eqn{n} is acceptable.}
  \item{Yn}{The response tensor instance with dimension \eqn{r_1\times r_2\times\cdots\times r_m \times n}, where \eqn{n} is the sample size. Array with the same dimensions and matrix with dimension \eqn{r\times n} are acceptable.}
  \item{u}{The dimension of envelope subspace. \eqn{u=(u_1,\cdots, u_m)}. Used for methods "FG", "1D", "ECD" and "PLS". User can use \code{\link{TensEnv_dim}} to select dimension.}
  \item{method}{The method used for estimation of tensor response regression. There are four possible choices.}
  \itemize{
\item{\code{"standard"}}: The standard OLS type estimation.
\item{\code{"FG"}}: Envelope estimation with full Grassmannian (FG) algorithm.
\item{\code{"1D"}}: Envelope estimation with one dimensional optimization approaches by 1D algorithm.
\item{\code{"ECD"}}: Envelope estimation with one dimensional optimization approaches by ECD algorithm.
\item{\code{"PLS"}}: The SIMPLS-type estimation without manifold optimization.
}
  \item{Gamma_init}{A list specifying the initial estimation of envelope subspace basis for "FG" method. If unspecified, use the estimation from "1D" algorithm.}

}

\value{
\item{Xn}{The input predictor matrix.}
\item{Yn}{The input response tensor.}
\item{method}{The method used.}
\item{coefficients}{The estimation of regression coefficient tensor.}
\item{Gamma_hat}{The estimation of envelope subspace basis.}
\item{Sig}{A matrix lists of \eqn{\boldsymbol{\Sigma}_k}, where \eqn{\boldsymbol{\Sigma}=\boldsymbol{\Sigma}_m \otimes \cdots \otimes \boldsymbol{\Sigma}_1}.}
\item{fitted.values}{The fitted response tensor.}
\item{residuals}{The residuals tensor.}
}

\examples{
rm(list=ls())

r <- c(10, 10, 10)
m <- length(r)
# The envelope dimensions u.
u <- c(2, 2, 2)
p <- 5
n <- 100

set.seed(1)
eta <- array(runif(prod(u)*p), c(u, p))
eta <- rTensor::as.tensor(eta)

# Gamma is the list of envelopes.
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
B <- rTensor::ttl(eta, Gamma, ms=1:m)
Xn <- matrix(rnorm(p*n), p, n)
Xn_inv <- MASS::ginv(Xn \%*\% t(Xn)) \%*\% Xn
Epsilon <- array(rnorm(prod(r)*n), c(r, n))
Epsilon <- rTensor::as.tensor(Epsilon)
Epsilon <- rTensor::ttl(Epsilon, Sigsqrtm, ms=1:m)
Yn <- Epsilon + rTensor::ttm(B, t(Xn), m+1)

res_std <- TRR(Xn, Yn, method="standard")
res_fg <- TRR(Xn, Yn, u, method="FG")
res_1D <- TRR(Xn, Yn, u, method="1D")
res_pls <- TRR(Xn, Yn, u, method="PLS")
res_ECD <- TRR(Xn, Yn, u, method="ECD")

rTensor::fnorm(B-stats::coef(res_std))
rTensor::fnorm(B-stats::coef(res_fg))
rTensor::fnorm(B-stats::coef(res_1D))
rTensor::fnorm(B-stats::coef(res_pls))
rTensor::fnorm(B-stats::coef(res_ECD))

## Use dataset bat
data("bat")
Xn <- bat$Xn
Yn <- bat$Yn
res_std <- TRR(Xn, Yn, method="standard")
}

\seealso{
\code{\link{plot.Tenv}, \link{predict.Tenv}, \link{vcov.Tenv}, \link{summary.Tenv}, \link{TensEnv_dim}}
}
