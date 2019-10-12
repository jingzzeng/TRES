\name{TensEnv_dim}
\alias{TensEnv_dim}

\title{
Envelope dimension selection for tensor response regression
}
\description{
This function use the 1D-BIC criterion by by Zhang, X., & Mai, Q. (2018) to select envelope dimensions in tensor response regression.
}
\usage{
TensEnv_dim(Xn, Yn, multiD=1, bic_max=10, opts=NULL)
}

\arguments{
  \item{Xn}{A vector predictor of dimension p.}
  \item{Yn}{The response tensor instance \eqn{r_1\times \cdots \times r_m}.}
  \item{multiD}{The parameter in \code{ballGBB1D_bic}.}
  \item{bic_max}{The maximum envelope dimension to be considered, the parameter in \code{ballGBB1D_bic}.}
  \item{opts}{The parameter in \code{ballGBB1D_bic}.}

}

\value{
\item{u}{The envelope dimension of \eqn{(u_1, u_2,\cdots,u_m)}.}
}

\examples{
rm(list=ls())

r <- c(10, 10, 10)
m <- length(r)
# The envelope dimensions u
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

TensEnv_dim(Xn, Yn) # The estimated envelope dimensions are the same as u.

## Use dataset bat, but it is time-consuming
\dontrun{
  data("bat")
  Xn <- bat$Xn
  Yn <- bat$Yn
  # check the dimension of Yn
  dim(Yn)
  # use 32 as the maximal envelope dimension
  TensEnv_dim(Xn, Yn, bic_max=32)
}
}


