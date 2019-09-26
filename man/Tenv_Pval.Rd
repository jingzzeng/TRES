\name{Tenv_Pval}
\alias{Tenv_Pval}

\title{
The \eqn{p}-values of each elements in the tensor response regression coefficient estimator
}
\description{
Obtain \eqn{p}-values of each elements in the tensor regression coefficient estimator. Two-sided t-tests are applied on the coefficent estimator, where the estimator use the asymptotic covariance of the OLS estimator.
}
\usage{
Tenv_Pval(Yn, Xn, B_est)
}

\arguments{
  \item{Yn}{The response tensor instance \eqn{ r_1\times r_2\times \cdots \times r_m}.}
  \item{Xn}{A vector predictor of dimension \eqn{p}.}
  \item{B_est}{
The tensor regression coefficient estimator for tensor response regression.
}
}

\value{
\item{P_OLS}{The p-value of OLS estimator.}
\item{P_val}{The p-value of \code{B_est}.}
\item{se_mat}{The standard error matrix of \eqn{\mathrm{vec}(\hat{B})}.}
}
