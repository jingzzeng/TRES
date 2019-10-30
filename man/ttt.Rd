\name{ttt}
\alias{ttt}

\title{
Inner product of tensor X and Y
}
\description{
\code{ttt(X,Y,ms)} computes the inner product of tensors \code{X} and \code{Y} in the
dimensions specified by the vector \code{ms}. The sizes of the dimensions specified by \code{ms} must match, which is \code{size(X, ms)} must equal \code{size(Y, ms)}.
}
\usage{
ttt(X, Y, ms)
}

\arguments{
  \item{X}{A tensor instance.}
  \item{Y}{A tensor instance.}
  \item{ms}{The indices of the modes to map onto the column space.}
}

\value{
\item{mat}{The inner product of tensors \code{X} and \code{Y} in the dimensions specified by the vector \code{ms}, which is a matrix.}
}

\examples{
 X = rTensor::as.tensor(array(runif(24), c(3, 4, 2)));
 Y = rTensor::as.tensor(array(runif(24), c(3, 4, 2)));
 Z = ttt(X, Y, 1:2)
}
