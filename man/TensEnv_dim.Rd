\name{TensEnv_dim}
\alias{TensEnv_dim}

\title{
Envelope dimension selection for tensor response regression (TRR)
}
\description{
This function uses the 1D-BIC criterion proposed by Zhang, X., & Mai, Q. (2018) to select envelope dimensions in tensor response regression.
}
\usage{
TensEnv_dim(Xn, Yn, multiD=1, maxdim=10, opts=NULL)
}

\arguments{
  \item{Xn}{A vector predictor of dimension p.}
  \item{Yn}{The response tensor instance \eqn{r_1\times \cdots \times r_m}.}
  \item{multiD}{The parameter in \code{ballGBB1D_bic}.}
  \item{maxdim}{The maximum envelope dimension to be considered, the parameter in \code{ballGBB1D_bic}.}
  \item{opts}{The parameter in \code{ballGBB1D_bic}.}

}

\value{
Return the envelope dimension of \eqn{(u_1, u_2,\cdots,u_m)}.
}

\examples{
rm(list=ls())

# The dimension of response
r <- c(10, 10, 10)
# The envelope dimensions u.
u <- c(2, 2, 2)
# The dimension of predictor
p <- 5
# The sample size
n <- 100

# Simulate the data with \code{\link{TRR_sim}}.
dat <- TRR_sim(r = r, p = p, u = u, n = n)
Xn <- dat$Xn
Yn <- dat$Yn

TensEnv_dim(Xn, Yn) # The estimated envelope dimensions are the same as u.

## Use dataset bat, but it is time-consuming
\dontrun{
  data("bat")
  Xn <- bat$Xn
  Yn <- bat$Yn
  # check the dimension of Yn
  dim(Yn)
  # use 32 as the maximal envelope dimension
  TensEnv_dim(Xn, Yn, maxdim=32)
}
}

\seealso{
\code{\link{TRR_sim}}.
}

