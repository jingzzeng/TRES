\name{TensPLS_fit}
\alias{TensPLS_fit}

\title{
Tensor envelope partial least squares (PLS) regression
}
\description{
This function estimates the factor matrix \eqn{W_k, k=1,\cdots,m} in tensor PLS algorithm for tensor predictor regression, see Zhang, X., & Li, L. (2017).
}
\usage{
TensPLS_fit(Xn, Yn, SigX, u)
}
\arguments{
  \item{Xn}{A predictor tensor of dimension \eqn{p_1\times \cdots \times p_m}.}
  \item{Yn}{The response vector of dimension \eqn{r}.}
  \item{SigX}{A matrix lists \eqn{\boldsymbol{\Sigma}_k, k=1,\cdots, m}, which determins the estimation of covariance matrix \eqn{\boldsymbol{\Sigma}=\boldsymbol{\Sigma}_m \otimes \cdots \otimes \boldsymbol{\Sigma}_1}.}
  \item{u}{The dimension of envelope subspace, \eqn{u=(u_1,\cdots,u_m)}.}
}

\value{
\item{Gamma}{The estimation of factor matrix \eqn{W_k, k=1,\cdots,m}.}
\item{PGamma}{The projection matrix \eqn{W_k(W_k'\boldsymbol{\Sigma}_k W_k)^{-1}W_k'\boldsymbol{\Sigma}_k, k=1,\cdots,m}.}
}
\references{
Zhang, X., & Li, L. (2017). Tensor Envelope Partial Least-Squares Regression. Technometrics, 59(4), 426-436.
}

