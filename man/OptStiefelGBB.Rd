\name{OptStiefelGBB}
\alias{OptStiefelGBB}

\title{
Optimization on Stiefel manifold (GBB)
}
\description{
Curvilinear search algorithm for optimization on Stiefel manifold based on Wen and Yin (2013). Used for estimating the envelope subspace.
}
\usage{
OptStiefelGBB(X, opts=NULL, fun, ...)
}

\arguments{
  \item{X}{n by k matrix such that \eqn{X'X = I}}
   \item{opts}{Option structure with fields:\cr
            \code{"record = 0"} -- no print out.\cr
            \code{"mxitr"} -- max number of iterations.\cr
            \code{"xtol"} -- stop control for \eqn{||X_k - X_{k-1}||}.\cr
            \code{"gtol"} -- stop control for the projected gradient.\cr
            \code{"ftol"} -- stop control for \eqn{\frac{|F_k - F_{k-1}|}{(1+|F_{k-1}|)}} usually with \code{max{xtol, gtol} > ftol}.\cr

The default values are: \code{"xtol"=1e-08; "gtol"=1e-08; "ftol"=1e-12; "mxitr"=500.}
}
\item{fun}{Objective function and its gradient:\cr
         \code{fun(X,  data1, data2)}\cr
         \code{data1, data2} are additional data.
         }
\item{...}{Additional input for \code{fun}, Calling syntax:\cr
         \code{OptStiefelGBB(X0, fun, opts, data1, data2)}.}
}

\value{
\item{X}{Solution}
\item{Out}{Output information, including estimation error, function value, iteration times etc.}
}
\references{
Wen, Z., & Yin, W. (2013). A feasible method for optimization with orthogonality constraints. Mathematical Programming, 142(1-2), 397-434.
}

\examples{
fun <- function(X,  W) {
      G = -2*(W \%*\% X)
      F = -sum(diag(t(X) \%*\% W \%*\% X))
      return(list(F = F, G = G))
}
n = 1000; k = 6;
W = matrix(rnorm(n^2), n, n); W = t(W) \%*\% W

opts=c()
opts$record = 0;
opts$mxitr  = 1000;
opts$xtol = 1e-5;
opts$gtol = 1e-5;
opts$ftol = 1e-8;

X0 = matrix(rnorm(n*k), n, k);
X0 = qr.Q(qr(X0));

eva <- OptStiefelGBB(X0, opts, fun, W)
X <- eva$X
out <- eva$out
out$fval = -2*out$fval;

}
