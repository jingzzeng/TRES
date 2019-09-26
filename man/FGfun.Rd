\name{FGfun}
\alias{FGfun}
\title{
The Objective function and its gradient}
\description{
The generic objective function listed below for estimating \code{M} envelope of \code{span(U)}
and its gradient value. For the detailed description, see Cook, R. D., & Zhang, X. (2016).

\deqn{F(\Gamma)=\log|\Gamma^T M \Gamma|+\log|\Gamma^T(M+U)^{-1}\Gamma|}
}
\usage{
FGfun(W, M, U)
}
\arguments{
\item{M}{M matrix in the envelope objective function. A \eqn{p}-by-\eqn{p} positive semi-definite matrix.}
 \item{U}{U matrix in the envelope objective function. A \eqn{p}-by-\eqn{p} positive semi-definite matrix.}
  \item{W}{A vector of \eqn{p} by \eqn{u}.}
}
\details{
This is the objective function and its gradient for estimating \code{M}-envelope contains \code{span(U)}, where \code{M > 0} and is symmetric, the
dimension of the envelope is \eqn{u}.
}
\value{
\item{F}{The value of objective function given \code{W}.}
\item{G}{The value of the gradient function given \code{W}.}
}
\references{
Cook, R. D., & Zhang, X. (2016). Algorithms for envelope estimation. Journal of Computational and Graphical Statistics, 25(1), 284-300.
}
