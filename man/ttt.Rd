\name{ttt}
\alias{ttt}

\title{
Inner product of tensor X and Y
}
\description{
\code{ttt(X,Y,dims)} computes the inner product of tensors \code{X} and \code{Y} in the
dimensions specified by the vector \code{dims}. The sizes of the dimensions specified by \code{dims} must match, which is \code{size(X, dims)} must equal \code{size(Y, dims)}.
}
\usage{
ttt(X, Y, dims)
}

\arguments{
  \item{X}{A tensor instance.}
  \item{Y}{A tensor instance.}
  \item{dims}{The indices of the modes to map onto the column space.}
}

\value{
\item{mat}{The inner product of tensors \eqn{X} and \eqn{Y} in the dimensions specified by the vector dims, which is a matrix.}
}

\examples{
 X = as.tensor(array(runif(24), c(3, 4, 2)));
 Y = as.tensor(array(runif(24), c(3, 4, 2)));
 Z = ttt(X, Y, 1:2)
}
