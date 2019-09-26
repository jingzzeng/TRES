\name{OptimballGBB1D}
\alias{OptimballGBB1D}
\title{
Estimate the envelope subspace (Feasi 1D)
}
\description{
The 1D algorithm to estimate the envelope subspace with specified dimension based on Wen and Yin (2013).
}
\usage{
OptimballGBB1D(M, U, u, opts=NULL)
}

\arguments{
 \item{M}{M matrix in the envelope objective function. A \eqn{p}-by-\eqn{p} positive semi-definite matrix.}
 \item{U}{U matrix in the envelope objective function. A \eqn{p}-by-\eqn{p} positive semi-definite matrix.}
  \item{u}{Dimension of the envelope. An integer between 0 and \eqn{p}.}
  \item{opts}{Option structure with fields:\cr
            \code{"record = 0"} -- no print out.\cr
            \code{"mxitr"} -- max number of iterations.\cr
            \code{"xtol"} -- stop control for \eqn{||X_k - X_{k-1}||}.\cr
            \code{"gtol"} -- stop control for the projected gradient.\cr
            \code{"ftol"} -- stop control for \eqn{\frac{|F_k - F_{k-1}|}{(1+|F_{k-1}|)}} usually with \code{max{xtol, gtol} > ftol}.\cr

The default values are: \code{"xtol"=1e-08; "gtol"=1e-08;}\cr
\code{"ftol"=1e-12; "mxitr"=500.}
}
}
\details{
Estimate \code{M}-envelope contains \code{span(U)}
where \code{M > 0} and is symmetric. The
dimension of the envelope is \code{u}.
}
\value{
\item{Ghat}{The orthogonal basis of the envelope subspace with each column represent the sequential direction. For example, the 1st column is the most informative direction.}
}
\references{
Wen, Z., & Yin, W. (2013). A feasible method for optimization with orthogonality constraints. Mathematical Programming, 142(1-2), 397-434.
}
\examples{
##simulate two matrices M and U with an envelope structure#
data <- MenvU_sim(n=200, p=20, u=5, Omega=NULL, Omega0=NULL, Phi=NULL)
Mhat <- data$Mhat
Uhat <- data$Uhat
G <- data$Gamma

Ghat_1D <- OptimballGBB1D(Mhat, Uhat, u=5)
subspace(Ghat_1D, G)
}
