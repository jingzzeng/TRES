\name{PMSE}
\alias{PMSE}

\title{
Predictions and Mean squared error for tensor predictor regression (TPR).
}
\description{
Evaluate tensor predictor regression through prediction mean squared error.
}
\usage{
PMSE(x, y, Bhat)
}

\arguments{
  \item{x}{A predictor tensor of dimension \eqn{p_1\times p_2\times\cdots\times p_m \times n}, where \eqn{n} is the sample size. Matrix and array are also allowed.}
  \item{y}{A response matrix of dimension \eqn{r \times n}, where \eqn{n} is the sample size. Vector is also allowed.}
  \item{Bhat}{An estimation of coefficient tensor of dimension \eqn{p_1\times p_2\times\cdots\times p_m \times r}. Vector, matrix and array are also allowed as long as the dimensions match the ones of \code{x} and \code{y}.}
}

\value{
\item{mse}{Mean squared error. Defined as \eqn{1/n\sum_{i=1}^n\|\mathbf{Y}_i-\hat{\mathbf{Y}}_i\|_F^2}, where \eqn{\hat{\mathbf{Y}}_i} is the prediction and \eqn{\|\cdot\|_F} is the Frobenius norm of tensor.}
\item{pred}{The predictions of tensor predictor regression.}
}


\examples{

rm(list = ls())
# The dimension of predictor
p <- c(10, 10, 10)
# The envelope dimensions u.
u <- c(1, 1, 1)
# The dimension of response
r <- 5
# The sample size
n <- 200

# Simulate the data with \code{\link{TPR_sim}}.
dat <- TPR_sim(p = p, r = r, u = u, n = n)
x <- dat$x
y <- dat$y

fit_std <- TPR.fit(x, y, u, method="standard")
result <- PMSE(x, y, fit_std$coefficients)
}

\seealso{
\code{\link{TPR.fit}, \link{TPR_sim}}.
}
