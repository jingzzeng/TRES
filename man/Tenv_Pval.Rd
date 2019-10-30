\name{Tenv_Pval}
\alias{Tenv_Pval}

\title{
The \eqn{p}-values of each elements in the tensor response regression coefficient estimator
}
\description{
Obtain \eqn{p}-values of each elements in the tensor regression coefficient estimator. Two-sided t-tests are applied on the coefficient estimator, where the estimator use the asymptotic covariance of the OLS estimator.
}
\usage{
Tenv_Pval(Xn, Yn, Bhat)
}

\arguments{
  \item{Xn}{The response tensor instance \eqn{ r_1\times r_2\times \cdots \times r_m}.}
  \item{Yn}{A vector predictor of dimension \eqn{p}.}
  \item{Bhat}{
The tensor regression coefficient estimator for tensor response regression.
}
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
res_std <- TRR(Xn, Yn, method="standard")
Tenv_Pval(Xn, Yn, res_std$coefficients)
}
