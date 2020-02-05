\name{OptStiefelGBB}
\alias{OptStiefelGBB}

\title{
Optimization on Stiefel manifold (GBB)
}
\description{
Curvilinear search algorithm for optimization on Stiefel manifold based on Wen and Yin (2013). Used for estimating the envelope subspace.
}
\usage{
OptStiefelGBB(Gamma_init, opts=NULL, fun, ...)
}

\arguments{
  \item{Gamma_init}{initial envelope subspace basis.}
   \item{opts}{Option structure with fields:\cr
            \code{"record = 0"} -- no print out.\cr
            \code{"maxiter"} -- max number of iterations.\cr
            \code{"xtol"} -- stop control for \eqn{||X_k - X_{k-1}||}.\cr
            \code{"gtol"} -- stop control for the projected gradient.\cr
            \code{"ftol"} -- stop control for \eqn{\frac{|F_k - F_{k-1}|}{(1+|F_{k-1}|)}} usually with \code{max{xtol, gtol} > ftol}.\cr

The default values are: \code{"xtol"=1e-08; "gtol"=1e-08; "ftol"=1e-12; "maxiter"=500.}
}
\item{fun}{Objective function and its gradient:\cr
         \code{fun(Gamma_init,  data1, data2)}\cr
         \code{data1, data2} are additional data.
         }
\item{...}{Additional input for \code{fun}, Calling syntax:\cr
         \code{OptStiefelGBB(Gamma_init, opts, fun, data1, data2)}.}
}

\value{
\item{Gamma}{The estimation of the envelope subspace basis.}
\item{Out}{Output information, including estimation error, function value, iteration times etc.}
}
\references{
Wen, Z., & Yin, W. (2013). A feasible method for optimization with orthogonality constraints. Mathematical Programming, 142(1-2), 397-434.
}

\examples{
fun <- function(X,  W) {
      G <- -2*(W \%*\% X)
      F <- -sum(diag(t(X) \%*\% W \%*\% X))
      return(list(F = F, G = G))
}
n <- 1000
k <- 6
W <- matrix(rnorm(n^2), n, n)
W <- t(W) \%*\% W

opts<-c()
opts$record <- 0;
opts$maxiter  <- 1000;
opts$xtol <- 1e-5;
opts$gtol <- 1e-5;
opts$ftol <- 1e-8;

Gamma0 <- matrix(rnorm(n*k), n, k);
Gamma0 <- qr.Q(qr(Gamma0));

eva <- OptStiefelGBB(Gamma0, opts, fun, W)
Gamma <- eva$Gamma
out <- eva$out
out$fval <- -2*out$fval;
}
