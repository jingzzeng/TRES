\name{ECD}
\alias{ECD}
\title{
ECD algorithm for estimating the envelope subspace
}
\description{
Estimate the envelope subspace with specified dimension based on ECD algorithm that described in Cook, R. D., & Zhang, X. (2018).
}
\usage{
ECD(M, U, u, maxiter=500, epsilon=1e-08)
}

\arguments{
 \item{M}{M matrix in the envelope objective function. A \eqn{p}-by-\eqn{p} positive semi-definite matrix.}
 \item{U}{U matrix in the envelope objective function. A \eqn{p}-by-\eqn{p} positive semi-definite matrix.}
 \item{u}{Envelope dimension. An integer between 0 and \eqn{p}.}
 \item{maxiter}{Maximum number of iterations.}
 \item{epsilon}{Convergence criterion. \eqn{|F_k - F_{k-1}|< \epsilon}, where \eqn{F_k}   is the objective function.}
}
\details{
Estimate \code{M}-envelope contains \code{span(U)}
where \code{M > 0} and is symmetric. The
dimension of the envelope is \code{u}.
}
\value{
\item{Ghat}{The orthogonal basis of the envelope subspace with each column represent the sequential direction. For example, the 1st column is the most informative direction.}
}
\references{
Cook, R. D., & Zhang, X. (2018). Fast envelope algorithms. Statistica Sinica, 28(3), 1179-1197.
}

\examples{

##simulate two matrices M and U with an envelope structure#
data <- MenvU_sim(n=200, p=20, u=5)
Mhat <- data$Mhat
Uhat <- data$Uhat

Ghat_ECD <- ECD(Mhat, Uhat, u=5)
}
