\name{TensEnv_dim}
\alias{TensEnv_dim}

\title{
Envelope dimension selection for tensor response regression (TRR)
}
\description{
This function uses the 1D-BIC criterion proposed by Zhang, X., & Mai, Q. (2018) to select envelope dimensions in tensor response regression.
}
\usage{
TensEnv_dim(x, y, multiD=1, maxdim=10, opts=NULL)
}

\arguments{
  \item{x}{The predictor matrix of dimension \eqn{p \times n}. Vector of length \eqn{n} is acceptable.}
  \item{y}{The response tensor instance with dimension \eqn{r_1\times r_2\times\cdots\times r_m \times n}, where \eqn{n} is the sample size. Array with the same dimensions and matrix with dimension \eqn{r\times n} are acceptable.}
  \item{multiD}{The parameter in \code{ballGBB1D_bic}. Default is 1.}
  \item{maxdim}{The maximum envelope dimension to be considered. Default is 10.}
  \item{opts}{The parameter in \code{ballGBB1D_bic}.}

}

\value{
Return the envelope dimension of \eqn{(u_1, u_2,\cdots,u_m)}.
}

\seealso{Arguments \code{multiD} and \code{opts} are the same as the ones in \code{\link{ballGBB1D_bic}}.

\code{\link{TRR_sim}}.
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

# Simulate the data with TRR_sim.
dat <- TRR_sim(r = r, p = p, u = u, n = n)
x <- dat$x
y <- dat$y

TensEnv_dim(x, y) # The estimated envelope dimensions are the same as u.

## Use dataset bat, but it is time-consuming
\dontrun{
  data("bat")
  x <- bat$x
  y <- bat$y
  # check the dimension of y
  dim(y)
  # use 32 as the maximal envelope dimension
  TensEnv_dim(x, y, maxdim=32)
}
}

