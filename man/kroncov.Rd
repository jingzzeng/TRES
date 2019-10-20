\name{kroncov}
\alias{kroncov}

\title{
The covariance estimation of a tensor random variable
}
\description{
This function estimates the covariance of a tensor random variable. We assume the covariance of the tensor random variable has a separable Kronecker covariance structure, i.e. \eqn{\boldsymbol{\Sigma}=\boldsymbol{\Sigma}_{m}\otimes\cdots\otimes\boldsymbol{\Sigma}_{1}}. This algorithm is described in Manceur, A. M., & Dutilleul, P. (2013).
}
\usage{
kroncov(Tn)
}

\arguments{
  \item{Tn}{A \eqn{p_1\times\cdots p_m\times n} data array, where \eqn{n} is the sample size.}
}

\value{
\item{lambda}{The normalizing constant.}
\item{S}{A matrix lists with each element being the individual estimation of the separable Kronecker covariance element \eqn{\boldsymbol{\Sigma}_m,\ldots,\boldsymbol{\Sigma}_1}.
}
}
\references{
Manceur, A. M., & Dutilleul, P. (2013). Maximum likelihood estimation for the tensor normal distribution: Algorithm, minimum sample size, and empirical bias and dispersion. Journal of Computational and Applied Mathematics, 239, 37-49.
}
