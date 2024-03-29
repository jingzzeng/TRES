% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TRRdim.R
\name{TRRdim}
\alias{TRRdim}
\title{Envelope dimension selection for tensor response regression (TRR)}
\usage{
TRRdim(x, y, C = NULL, maxdim = 10, ...)
}
\arguments{
\item{x}{The predictor matrix of dimension \eqn{p \times n}. Vector of length \eqn{n} is acceptable.}

\item{y}{The response tensor instance with dimension \eqn{r_1\times r_2\times\cdots\times r_m \times n}, where \eqn{n} is the sample size. Array with the same dimensions and matrix with dimension \eqn{r\times n} are acceptable.}

\item{C}{The parameter passed to \code{\link{oneD_bic}}. Default is \code{nrow(x) = p}.}

\item{maxdim}{The maximum envelope dimension to be considered. Default is 10.}

\item{...}{Additional arguments passed to \code{\link{oneD_bic}}.}
}
\value{
\item{bicval}{The minimal BIC values for each mode.}
\item{u}{The optimal envelope subspace dimension \eqn{(u_1, u_2,\cdots,u_m).}}
\item{mse}{The prediction mean squared error using the selected envelope basis.}
}
\description{
This function uses the 1D-BIC criterion proposed by Zhang, X., & Mai, Q. (2018) to select envelope dimensions in tensor response regression. Refer to \code{\link{oneD_bic}} for more details.
}
\details{
See \code{\link{oneD_bic}} for more details on the definition of 1D-BIC criterion and on the arguments \eqn{C} and the additional arguments.

Let \eqn{B} denote the estimated envelope with the selected dimension \code{u}, then the prediction is \eqn{\hat{Y}_i = B \bar{\times}_{(m+1)} X_i} for each observation. And the mean squared error is defined as \eqn{1/n\sum_{i=1}^n||Y_i-\hat{Y}_i||_F^2}, where \eqn{||\cdot||_F} denotes the Frobenius norm.
}
\examples{
# The dimension of response
r <- c(10, 10, 10)
# The envelope dimensions u.
u <- c(2, 2, 2)
# The dimension of predictor
p <- 5
# The sample size
n <- 100

# Simulate the data with TRRsim.
dat <- TRRsim(r = r, p = p, u = u, n = n)
x <- dat$x
y <- dat$y

TRRdim(x, y) # The estimated envelope dimensions are the same as u.

## Use dataset bat. (time-consuming)
\donttest{
data("bat")
x <- bat$x
y <- bat$y
# check the dimension of y
dim(y)
# use 32 as the maximal envelope dimension
TRRdim(x, y, maxdim=32)
}

}
\references{
Li, L. and Zhang, X., 2017. Parsimonious tensor response regression. Journal of the American Statistical Association, 112(519), pp.1131-1146.

Zhang, X. and Mai, Q., 2018. Model-free envelope dimension selection. Electronic Journal of Statistics, 12(2), pp.2193-2216.
}
\seealso{
\code{\link{oneD_bic}}, \code{\link{TRRsim}}.
}
