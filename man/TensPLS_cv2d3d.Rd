\name{TensPLS_cv2d3d}
\alias{TensPLS_cv2d3d}
\title{
Envelope dimension by cross-validation for tensor predictor regression (TPR).
}
\description{
Select the envelope dimension by cross-validation for tensor predictor regression.
}
\usage{
TensPLS_cv2d3d(x, y, maxdim=10, nfolds=5)
}

\arguments{
  \item{x}{The predictor tensor instance of dimension \eqn{p_1\times p_2\times\cdots\times p_m \times n}, where \eqn{n} is the sample size. Array with the same dimensions and matrix with dimension \eqn{p\times n} are acceptable.}
  \item{y}{The response matrix of dimension \eqn{r \times n}, where \eqn{n} is the sample size. Vector of length \eqn{n} is acceptable.}
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
x <- dat$x
y <- dat$y

## It is time-consuming
\dontrun{
  TensPLS_cv2d3d(x, y, maxdim = 5) # The estimated envelope dimensions are the same as u.
}

## Use dataset square, but it is time-consuming
\dontrun{
  data("square")
  x <- square$x
  y <- square$y
  # check the dimension of x
  dim(x)
  # use 32 as the maximal envelope dimension
  TensPLS_cv2d3d(x, y, maxdim=32)
}
}

\seealso{
\code{\link{TPR_sim}}.
}
