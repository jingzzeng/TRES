\name{TensPLS_cv2d3d}
\alias{TensPLS_cv2d3d}
\title{
Envelope dimension by cross-validation for tensor predictor regression
}
\description{
This function obtain the envelope dimension by cross-validation for tensor predictor regression.
}
\usage{
TensPLS_cv2d3d(X0, Y0, maxdim, nfolds)
}

\arguments{
   \item{X0}{A predictor tensor instance.}
  \item{Y0}{The response vector.}
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
