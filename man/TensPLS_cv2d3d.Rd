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
# The dimension of predictor
p <- c(10, 10, 10)
# The envelope dimensions u.
u <- c(1, 1, 1)
# The dimension of response
r <- 5
# The sample size
n <- 200

dat <- TPR_sim(p = p, r = r, u = u, n = n)
Xn <- dat$Xn
Yn <- dat$Yn

## It is time-consuming
\dontrun{
  TensPLS_cv2d3d(Xn, Yn, maxdim = 5) # The estimated envelope dimensions are the same as u.
}

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

\seealso{
\code{\link{TPR_sim}}.
}
