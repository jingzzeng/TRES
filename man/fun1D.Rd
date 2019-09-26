\name{fun1D}
\alias{fun1D}
\title{
The 1D objective function and its gradient
}
\description{
The objective function and its gradient value that defined in equation (4.1) of Cook, R. D., & Zhang, X. (2016). A speicial case of \code{FGfun} where \code{W} is a one-dimensional vector.\cr
}
\usage{
fun1D(W, M, U)
}
\arguments{
 \item{M}{M matrix in the envelope objective function. A \eqn{p}-by-\eqn{p} positive semi-definite matrix.}
 \item{U}{U matrix in the envelope objective function. A \eqn{p}-by-\eqn{p} positive semi-definite matrix.}
  \item{W}{A vector of \eqn{p} by 1.}
}
\details{
This is the objective function and its gradient for the constrained optimization in the 1D algorithm.
}
\value{
\item{F}{The value of objective function given \code{W}.}
\item{G}{The value of the gradient function given \code{W}.}
}
\references{
Cook, R. D., & Zhang, X. (2016). Algorithms for envelope estimation. Journal of Computational and Graphical Statistics, 25(1), 284-300.
}
