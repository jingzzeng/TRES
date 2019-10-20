\name{TPR}
\alias{TPR}

\title{
Tensor Predictor Regression
}
\description{
This function is used for estimation of tensor predictor regression. The available method including standard OLS type estimation, PLS type of estimation as well as envelope estimation with 1D and ECD approaches.
}
\usage{
TPR(Xn, Yn, method=c('standard', 'FG', '1D', 'ECD', 'PLS'),
    u=NULL, Gamma_init=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Xn}{The predictor tensor instance of dimension \eqn{p_1\times p_2\times\cdots\times p_m \times n}, where \eqn{n} is the sample size. Array with the same dimensions and matrix with dimension \eqn{p\times n} are acceptable.}
  \item{Yn}{The response matrix of dimension \eqn{r \times n}, where \eqn{n} is the sample size. Vector of length \eqn{n} is acceptable}
  \item{u}{The dimension of envelope subspace. \eqn{u=(u_1,\cdots, u_m)}. Used for methods "FG", "1D", "ECD" and "PLS". User can use \code{\link{TensPLS_cv2d3d}} to select dimension.}
  \item{method}{The method used for estimation of tensor response regression. There are four possible choices.}
  \itemize{
  \item{\code{"standard"}}: The standard OLS type estimation.
  \item{\code{"FG"}}: Envelope estimation with full Grassmannian (FG) algorithm.
  \item{\code{"1D"}}: Envelope estimation with one dimenaional optimization approaches by 1D algorithm.
  \item{\code{"ECD"}}: Envelope estimation with one dimenaional optimization approaches by ECD algorithm.
  \item{\code{"PLS"}}: The SIMPLS-type estimation without mannifold optimization.
}
  \item{Gamma_init}{A list specifying the initial estimation of envelope subspace basis for "FG" method. If unspecified, use the estimation from "1D" algorithm.}
}

\value{
\item{Xn}{The input predictor tensor.}
\item{Yn}{The input response matrix.}
\item{method}{The method used.}
\item{coefficients}{The estimation of regression coefficient tensor.}
\item{Gamma_hat}{The estimation of envelope subspace basis.}
\item{Sig}{A matrix lists of \eqn{\boldsymbol{\Sigma}_\mathrm{xk}}, where \eqn{\boldsymbol{\Sigma}_\mathrm{x}=\boldsymbol{\Sigma}_\mathrm{xm} \otimes \cdots \otimes \boldsymbol{\Sigma}_\mathrm{x1}}.}
\item{fitted.values}{The fitted response matrix.}
\item{residuals}{The residuals matrix.}
}

\note{
This function supports tensor predictor regression for 1-dimensional, 2-dimensional or 3-dimensional tensor predictor in current version.
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
Xn <- dat$Xn
Yn <- dat$Yn
B <- dat$Bhat

res_std = TPR(Xn, Yn, method="standard")
res_FG = TPR(Xn, Yn, u, method="FG")
res_pls = TPR(Xn, Yn, u, method="PLS")

rTensor::fnorm(B-stats::coef(res_std))
rTensor::fnorm(B-stats::coef(res_FG))
rTensor::fnorm(B-stats::coef(res_pls))

## Use dataset square, but it is time-consuming
\dontrun{
  data("square")
  Xn <- square$Xn
  Yn <- square$Yn
  res_std <- TPR(Xn, Yn, method="standard")
  }
}

\seealso{
\code{\link{plot.Tenv}, \link{predict.Tenv}, \link{vcov.Tenv}, \link{summary.Tenv}, \link{TensPLS_cv2d3d}, \link{TPR_sim}}.
}
