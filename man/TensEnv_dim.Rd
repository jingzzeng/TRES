\name{TensEnv_dim}
\alias{TensEnv_dim}

\title{
Envelope dimension selection for tensor response regression
}
\description{
This function use the 1D-BIC criterion by by Zhang, X., & Mai, Q. (2018) to select envelope dimensions in tensor response regression.
}
\usage{
TensEnv_dim(Yn, Xn, multiD=1, bic_max=10, opts=NULL)
}

\arguments{
  \item{Xn}{A vector predictor of dimension p.}
  \item{Yn}{The response tensor instance \eqn{r_1\times \cdots \times r_m}.}
  \item{multiD}{The parameter in \code{ballGBB1D_bic}.}
  \item{bic_max}{The maximum envelope dimension to be considered, the parameter in \code{ballGBB1D_bic}.}
  \item{opts}{The parameter in \code{ballGBB1D_bic}.}

}

\value{
\item{u}{The envelope dimension of \eqn{(u_1, u_2,\cdots,u_m)}.}
}

