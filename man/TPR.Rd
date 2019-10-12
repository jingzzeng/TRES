\name{TPR}
\alias{TPR}

\title{
Tensor Predictor Regression
}
\description{
This function is used for estimation of tensor predictor regression. The avaiable method including standard OLS type estimation, PLS type of estimation as well as envelope estimation with 1D and ECD approaches.
}
\usage{
TPR(Xn, Yn, method=c('standard', 'FG', '1D', 'ECD', 'PLS'),
    u=NULL, Gamma_init=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Xn}{The predictor tensor instance of dimension \eqn{p_1\times p_2\times\cdots\times p_m \times n}, where \eqn{n} is the sample size. Array with the same dimensions and matrix with dimension \eqn{p\times n} are acceptable.}
  \item{Yn}{The response matrix of dimension \eqn{r \times n}, where \eqn{n} is the sample size. Vector of length \eqn{n} is acceptable}
  \item{u}{The dimension of envelope subspace. \eqn{u=(u_1,\cdots, u_m)}. Used for methods "FG", "1D", "ECD" and "PLS". User can use \code{\link{TensPLS_cv2d3d}} to select dimension.}
  \item{method}{The method used for estimation of tensor response regression. There are four possible choices.}
  \itemize{
  \item{\code{"standard"}}: The standard OLS type estimation.
  \item{\code{"FG"}}: Envelope estimation with full Grassmannian (FG) algorithm.
  \item{\code{"1D"}}: Envelope estimation with one dimenaional optimization approaches by 1D algorithm.
  \item{\code{"ECD"}}: Envelope estimation with one dimenaional optimization approaches by ECD algorithm.
  \item{\code{"PLS"}}: The SIMPLS-type estimation without mannifold optimization.
}
  \item{Gamma_init}{A list specifying the initial estimation of envelope subspace basis for "FG" method. If unspecified, use the estimation from "1D" algorithm.}
}

\value{
\item{Xn}{The input predictor tensor.}
\item{Yn}{The input response matrix.}
\item{method}{The method used.}
\item{coefficients}{The estimation of regression coefficient tensor.}
\item{Gamma_hat}{The estimation of envelope subspace basis.}
\item{Sig}{A matrix lists of \eqn{\boldsymbol{\Sigma}_\mathrm{xk}}, where \eqn{\boldsymbol{\Sigma}_\mathrm{x}=\boldsymbol{\Sigma}_\mathrm{xm} \otimes \cdots \otimes \boldsymbol{\Sigma}_\mathrm{x1}}.}
\item{fitted.values}{The fitted response matrix.}
\item{residuals}{The residuals matrix.}
}

\note{
This function supports tensor predictor regression for 1-dimensional, 2-dimensional or 3-dimensional tensor predictor in current version.
}


\examples{
rm(list = ls())

p <- c(10, 10, 10)
# The envelope dimensions u.
u <- c(1, 1, 1)
m <- 3; r <- 5; n <- 200
eta <- array(runif(prod(u,r)), c(u,r))
eta <- rTensor::as.tensor(eta)

# Gamma is the list of envelopes.
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
Yn <-  Y_tmp + t(Epsilon)


res_std = TPR(Xn, Yn, method="standard")
res_FG = TPR(Xn, Yn, u, method="FG")
res_1D = TPR(Xn, Yn, u, method="1D")
res_pls = TPR(Xn, Yn, u, method="PLS")
res_ECD = TPR(Xn, Yn, u, method="ECD")

rTensor::fnorm(B-stats::coef(res_std))
rTensor::fnorm(B-stats::coef(res_FG))
rTensor::fnorm(B-stats::coef(res_1D))
rTensor::fnorm(B-stats::coef(res_pls))
rTensor::fnorm(B-stats::coef(res_ECD))

## Use dataset square
data("square")
Xn <- square$Xn
Yn <- square$Yn
res_std <- TPR(Xn, Yn, method="standard")
}

\seealso{
\code{\link{plot.Tenv}, \link{predict.Tenv}, \link{vcov.Tenv}, \link{summary.Tenv}, \link{TensPLS_cv2d3d}}
}
