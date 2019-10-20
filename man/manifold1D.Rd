\name{manifold1D}
\alias{manifold1D}

\title{
Estimate the envelope subspace (ManifoldOptim 1D)
}
\description{
The 1D algorithm to estimate the envelope subspace with specified dimension based on R package "ManifoldOptim".
}
\usage{
manifold1D(M, U, u, params=NULL)
}

\arguments{
  \item{M}{M matrix in the envelope objective function. A \eqn{p}-by-\eqn{p} positive semi-definite matrix.}
 \item{U}{U matrix in the envelope objective function. A \eqn{p}-by-\eqn{p} positive semi-definite matrix.}
  \item{u}{Dimension of the envelope. An integer between 0 and \eqn{p}.}
  \item{params}{Option structure with fields:\cr
            \code{"max_iter"} -- max number of iterations.\cr
            \code{"tol"} -- Tolerance used to assess convergence. See Huang et al. (2018) for details on
how this is used.\cr
            \code{"method"} -- Name of optimization method supported by R package \pkg{ManifoldOptim}.
\itemize{
\item{\code{"LRBFGS"}}: Limited-memory RBFGS
\item{\code{"LRTRSR1"}}: Limited-memory RTRSR1
\item{\code{"RBFGS"}}: Riemannian BFGS
\item{\code{"RBroydenFamily"}}: Riemannian Broyden family
\item{\code{"RCG"}}: Riemannian conjugate gradients
\item{\code{"RNewton"}}: Riemannian line-search Newton
\item{\code{"RSD"}}: Riemannian steepest descent
\item{\code{"RTRNewton"}}: Riemannian trust-region Newton
\item{\code{"RTRSD"}}: Riemannian trust-region steepest descent
\item{\code{"RTRSR1"}}: Riemannian trust-region symmetric rank-one update
\item{\code{"RWRBFGS"}}: Riemannian BFGS
}

\code{"check"} -- Should internal manifold object check inputs and print summary message before optimization (TRUE or FALSE).\cr

The default values are: \code{"max_iter"=500; "tol"=1e-08;}\cr
\code{"method"="RCG"; "check"="False"}.
}
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
Huang, W., Absil, P. A., Gallivan, K. A., & Hand, P. (2018). ROPTLIB: an object-oriented C++ library for optimization on Riemannian manifolds. ACM Transactions on Mathematical Software (TOMS), 44(4), 43.
}

\examples{
##simulate two matrices M and U with an envelope structure#
data <- MenvU_sim(n=200, p=20, u=5)
Mhat <- data$Mhat
Uhat <- data$Uhat
G <- data$Gamma

Ghat_1D <- manifold1D(Mhat, Uhat, u=5, params=NULL)
subspace(Ghat_1D, G)
}
