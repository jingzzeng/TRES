\name{Tenv}
\alias{Tenv}
\title{
Tensor response envelope estimator
}
\description{
This function gives the tensor envelope estimator for tensor response regression.
}
\usage{
Tenv(Yn, Xn, u, opts=NULL)
}
\arguments{
  \item{Yn}{The response tensor instance \eqn{ r_1\times r_2\times \cdots \times r_m \times n}, where \eqn{n} is the sample size.}
   \item{Xn}{The predictor matrix of dimension \eqn{p \times n}.}
  \item{u}{The dimension of envelope subspace. \eqn{u=(u_1,\cdots,u_m)}.}
  \item{opts}{The option structure for Feasi. See function \code{OptimballGBB1D}.}
}

\value{
\item{Btil}{The ordinary least square estimator (OLS).}
\item{Bhat}{The envelope based estimator.}
\item{PGamma}{The projection matrix onto envelope subspace.}
}
\references{
Li, L., & Zhang, X. (2017). Parsimonious tensor response regression. Journal of the American Statistical Association, 112(519), 1131-1146.
}

