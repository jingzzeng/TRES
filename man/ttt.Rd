\name{ttt}
\alias{ttt}

\title{
Inner product of two tensors
}
\description{
Computes the inner product of two tensors on the specified modes.
}
\usage{
ttt(X, Y, ms)
}

\arguments{
  \item{X}{A tensor instance.}
  \item{Y}{A tensor instance.}
  \item{ms}{The indices of the modes to compute on.}
}

\value{
Return a matrix, the inner product of tensors \code{X} and \code{Y} on the modes specified by the vector \code{ms}.
}

\details{
Two tensors \code{X} and \code{Y} are unfolded into two matrix with the modes in \code{ms} onto the columns. Then implement matrix multiplication on the two matrices. The sizes of the dimensions specified by \code{ms} must match, e.g., the number of columns of the two unfolded matrices must match.
}

\examples{
 X = rTensor::as.tensor(array(runif(24), c(3, 4, 2)));
 Y = rTensor::as.tensor(array(runif(24), c(3, 4, 2)));
 Z = ttt(X, Y, 1:2)
}
