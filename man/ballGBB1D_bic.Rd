\name{ballGBB1D_bic}
\alias{ballGBB1D_bic}
\title{
Envelope dimension selection based on 1D-BIC
}
\description{
This function selects envelope subspace dimension using 1D-BIC proposed by Zhang, X., & Mai, Q. (2018). The constrained optimization in the 1D algorithm is based on the method proposed by Wen and Yin (2013).
}
\usage{
ballGBB1D_bic(M, U, n, multiD=1, maxdim=10, opts=NULL)
}
\arguments{
\item{M}{M matrix in the envelope objective function. A \eqn{p}-by-\eqn{p} positive semi-definite matrix.}
 \item{U}{U matrix in the envelope objective function. A \eqn{p}-by-\eqn{p} positive semi-definite matrix.}
 \item{n}{Sample size.}
 \item{multiD}{A constant, the default value is 1.}
 \item{maxdim}{The maximum dimension to consider, \code{maxdim} is smaller than \eqn{p}, the default value is 10.}
 \item{opts}{Option structure for GBB algorithm. See function \code{\link{OptStiefelGBB}}.}
}

\value{
\item{bicval}{The BIC values for different envelope dimensions.}
\item{u}{The dimension selected which corresponds to the smallest BIC values.}
}
\references{
Zhang, X., & Mai, Q. (2018). Model-free envelope dimension selection. Electronic Journal of Statistics, 12(2), 2193-2216.\cr
Wen, Z., & Yin, W. (2013). A feasible method for optimization with orthogonality constraints. Mathematical Programming, 142(1-2), 397-434.
}

\examples{

##simulate two matrices M and U with an envelope structure#
data <- MenvU_sim(n=200, p=20, u=5)
Mhat <- data$Mhat
Uhat <- data$Uhat

fit <- ballGBB1D_bic(Mhat, Uhat, n=200)

## visualization
plot(1:10, fit$bicval, type="o", xlab="Envelope Dimension", ylab="BIC values",
main="Envelope Dimension Selection")
}
