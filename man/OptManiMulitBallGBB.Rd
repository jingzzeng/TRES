\name{OptManiMulitBallGBB}
\alias{OptManiMulitBallGBB}

\title{
Line search algorithm for optimization on manifold
}
\description{
Line search algorithm for optimization on Stiefel manifold based on Wen and Yin (2013). Used for the 1D-algorithm to estimate the envelope subspace.
}
\usage{
OptManiMulitBallGBB(X, opts=NULL, fun, ...)
}
\arguments{
  \item{X}{n by k matrix such that \eqn{X'X = I}}
   \item{opts}{Option structure with fields:\cr
            \code{"record = 0"} -- no print out.\cr
            \code{"mxitr"} -- max number of iterations.\cr
            \code{"xtol"} -- stop control for \eqn{||X_k - X_{k-1}||}.\cr
            \code{"gtol"} -- stop control for the projected gradient.\cr
            \code{"ftol"} -- stop control for \eqn{\frac{|F_k - F_{k-1}|}{(1+|F_{k-1}|)}} usually with \code{max{xtol, gtol} > ftol}.
}
\item{fun}{Objective function and its gradient:\cr
         \code{fun(X,  data1, data2)}\cr
         \code{data1, data2} are additional data.
         }
\item{...}{Additional input for \code{fun}, Calling syntax:\cr
        \code{OptManiMulitBallGBB(X0, fun, opts, data1, data2)}.}
}

\value{
\item{X}{Solution.}
\item{g}{Gradient of X.}
\item{Out}{Output information, include estimation error, function value, iteration times etc.}
}
\references{
Wen, Z., & Yin, W. (2013). A feasible method for optimization with orthogonality constraints. Mathematical Programming, 142(1-2), 397-434.
}

