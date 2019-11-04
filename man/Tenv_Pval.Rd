\name{Tenv_Pval}
\alias{Tenv_Pval}

\title{
The \eqn{p}-value and standard error of each element in the estimator of the tensor response regression (TRR) coefficient
}
\description{
Obtain \eqn{p}-value of each element in the tensor regression coefficient estimator. Two-sided t-tests are implemented on the coefficient estimator, where asymptotic covariance of the OLS estimator is used.
}
\usage{
Tenv_Pval(Xn, Yn, Bhat)
}

\arguments{
  \item{Xn}{The response tensor instance \eqn{ r_1\times r_2\times \cdots \times r_m}.}
  \item{Yn}{A vector predictor of dimension \eqn{p}.}
  \item{Bhat}{The estimator of tensor regression coefficient.}
}

\note{
The \eqn{p}-value and the standard error of estimated coefficient are not provided for tensor predictor regression since they depend on \eqn{\widehat{\mathrm{cov}}^{-1}\{\mathrm{vec}(\mathbf{X})\}} which is unavailable due to the ultra-high dimension of \eqn{\mathrm{vec}(\mathbf{X})}.
}

\value{
\item{p_ols}{The p-value tensor of OLS estimator.}
\item{p_val}{The p-value tensor of \code{Bhat}.}
\item{se}{The standard error tensor of \code{Bhat}.}
}


\examples{
## Use dataset bat
data("bat")
Xn <- bat$Xn
Yn <- bat$Yn
res_std <- TRR.fit(Xn, Yn, method="standard")
Tenv_Pval(Xn, Yn, res_std$coefficients)
}
