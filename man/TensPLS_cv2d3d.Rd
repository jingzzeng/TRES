\name{TensPLS_cv2d3d}
\alias{TensPLS_cv2d3d}
\title{
Envelope dimension by cross-validation for tensor predictor regression
}
\description{
This function obtain the envelope dimension by cross-validation for tensor predictor regression.
}
\usage{
TensPLS_cv2d3d(Xn, Yn, maxdim=10, nfolds=5)
}

\arguments{
  \item{Xn}{A predictor tensor instance.}
  \item{Yn}{The response vector.}
  \item{maxdim}{The largest dimension to be considered for selection.}
  \item{nfolds}{Number of folds for cross-validation.}
}

\value{
\item{mincv}{The minimum sum of squared error.}
\item{u}{The envelope subspace dimension selected.}
}

\references{
Zhang, X., & Li, L. (2017). Tensor Envelope Partial Least-Squares Regression. Technometrics, 59(4), 426-436.
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

TensPLS_cv2d3d(Xn, Yn, maxdim = 10) # The estimated envelope dimensions are the same as u.

## Use dataset square, but it is time-consuming
\dontrun{
  data("square")
  Xn <- square$Xn
  Yn <- square$Yn
  # check the dimension of Xn
  dim(Xn)
  # use 32 as the maximal envelope dimension
  TensPLS_cv2d3d(Xn, Yn, maxdim=32)
}
}
