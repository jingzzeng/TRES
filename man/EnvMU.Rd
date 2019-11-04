\name{EnvMU}
\alias{EnvMU}

\title{
Estimate envelope subspace basis
}
\description{
SIMPLS-type algorithm for estimating the envelope subspace without manifold optimization.
}
\usage{
EnvMU(M, U, u)
}
\arguments{
 \item{M}{M matrix in the envelope objective function. An \eqn{p}-by-\eqn{p} positive semi-definite matrix.}
 \item{U}{U matrix in the envelope objective function. An \eqn{p}-by-\eqn{p} positive semi-definite matrix.}
 \item{u}{The envelope dimension.}
}
\details{
Estimate \code{M}-envelope of \code{span(U)} where \code{M > 0}. The dimension of the envelope is \code{u}. This algorithm is a generalization of De Jong, S. (1993) and Cook, R. D., Helland, I. S., & Su, Z. (2013). It generalizes from predictor envelopes to an arbitrary \code{M}-envelope of \code{span(U)}.
}

\value{
Return the orthogonal basis of the envelope subspace.
}
\references{
Cook, R. D., Helland, I. S., Su, Z. (2013). Envelopes and partial least squares regression. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 75(5), 851-877.

De Jong, S. (1993). SIMPLS: an alternative approach to partial least squares regression. Chemometrics and intelligent laboratory systems, 18(3), 251-263.
}
