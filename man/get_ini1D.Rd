\name{get_ini1D}
\alias{get_ini1D}
\title{
Initial value for 1D algorithm
}
\description{
This function gives the initial value for 1D algorithm, which is searched from the eigenvectors of matrices \code{M} and \code{M+U} corresponding to their objective function value.
}
\usage{
get_ini1D(M, U)
}

\arguments{
 \item{M}{M matrix in the envelope objective function. A \eqn{p}-by-\eqn{p} positive semi-definite matrix.}
 \item{U}{U matrix in the envelope objective function. A \eqn{p}-by-\eqn{p} positive semi-definite matrix.}
}

\value{
\item{W0}{Initial value for the 1D algorithm}
}

