#' Generate simulation data for tensor response regression (TRR)
#'
#' This function is used to generate simulation data used in tensor response regression.
#'
#' The tensor response regression model is of the form,
#' \deqn{\mathbf{Y} = \mathbf{B}\bar{\times}_{(m+1)}\mathbf{X} +\boldsymbol{\varepsilon}}
#' where predictor \eqn{\mathbf{X} \in R^{p}}, response \eqn{\mathbf{Y} \in R^{r_1\times \cdots\times r_m}}, and
#' the error term is tensor normal distributed
#' \deqn{\boldsymbol{\varepsilon}\sim TN(0;\boldsymbol{\Sigma}_1,\dots,\boldsymbol{\Sigma}_m)}
#' According to the tensor envelope structure, we have
#' \deqn{\mathbf{B} = [\Theta;\boldsymbol{\Gamma}_1,\ldots,\boldsymbol{\Gamma}_m,\mathbf{I}_p], \quad \mbox{for some } \boldsymbol{\Theta} \in R^{u_1\times\cdots\times u_m \times p}}
#' \deqn{\boldsymbol{\Sigma}_k = \boldsymbol{\Gamma}_k\boldsymbol{\Omega}_k\boldsymbol{\Gamma}_k^{T}+\boldsymbol{\Gamma}_{0k}\boldsymbol{\Omega}_{0k}\boldsymbol{\Gamma}_{0k}^\top,
#' \quad \mbox{for some } \boldsymbol{\Omega}_k, \boldsymbol{\Omega}_{0k},\ k=1,\ldots,m.}
#'
#' @note The length of r must match the length of u, and each element of u must be less than the corresponding element in r.
#'
#' @param r The dimension of response, a vector with length larger than 2.
#' @param p The dimension of predictor, a scale.
#' @param u The structural dimension of envelopes at each mode, a vector with the same length as r.
#' @param n The sample size.
#' @return
#' \item{Xn}{The predictor of dimension \eqn{p\times n}}
#' \item{Yn}{The response of dimension \eqn{r_1\times \cdots\times r_m \times n}}
#' \item{Gamma}{The envelope subspace basis of dimension \eqn{r_k \times u_k, \ k=1,\ldots,m}}
#' \item{coefficients}{The tensor coefficients of dimension \eqn{r_1\times \cdots\times r_m \times p}}
#' \item{Sigma}{The covariance matrix of error}
#' \item{p, r, u}{The input \code{p,r,u}}
#'
#' @examples
#' r <- c(10, 10, 10)
#' u <- c(2, 2, 2)
#' p <- 5
#' n <- 100
#' dat <- TRR_sim(r = r, p = p, u = u, n = n)
#' Xn <- dat$Xn
#' Yn <- dat$Yn
#' @seealso \code{\link{TPR_sim}}.
#' @references Li L, Zhang X (2017). “Parsimonious Tensor Response Regression.” Journal of the American Statistical Association, 112(519), 1131–1146.

#' @export
#' @importFrom stats rnorm
TRR_sim <- function(r, p, u, n){
  if(length(r) != length(u)){stop("The length of r must match the length of u.")}
  if(any(r < u)){stop("u must be less than r element-wise.")}
  m <- length(r)

  # Gamma is the list of envelopes.
  Gamma <- Gamma0 <- Omega <- Omega0 <- Sig <- Sigsqrtm <- NULL
  for (i in 1:m){
    tmp <- matrix(runif(r[i]*u[i]), r[i], u[i])
    Gamma[[i]] <- qr.Q(qr(tmp))
    Gamma0[[i]] <- qr.Q(qr(Gamma[[i]]),complete=TRUE)[,(u[i]+1):r[i]]
    A <- matrix(runif(u[i]^2), u[i], u[i])
    Omega[[i]] <- A %*% t(A)
    A <- matrix(runif((r[i]-u[i])^2), (r[i]-u[i]), (r[i]-u[i]))
    Omega0[[i]] <- A %*% t(A)
    Sig[[i]] <- Gamma[[i]] %*% Omega[[i]] %*% t(Gamma[[i]])+
      Gamma0[[i]] %*% Omega0[[i]] %*% t(Gamma0[[i]])
    Sig[[i]] <- 10*Sig[[i]]/norm(Sig[[i]], type="F")+0.01*diag(r[i])
    Sigsqrtm[[i]] <- pracma::sqrtm(Sig[[i]])$B
  }
  eta <- array(runif(prod(u)*p), c(u, p))
  eta <- rTensor::as.tensor(eta)
  B <- rTensor::ttl(eta, Gamma, ms=1:m)
  Xn <- matrix(rnorm(p*n), p, n)
  Epsilon <- array(rnorm(prod(r)*n), c(r, n))
  Epsilon <- rTensor::as.tensor(Epsilon)
  Epsilon <- rTensor::ttl(Epsilon, Sigsqrtm, ms=1:m)
  Yn <- Epsilon + rTensor::ttm(B, t(Xn), m+1)

  output <- list(Xn = Xn, Yn = Yn, Gamma = Gamma, coefficients = B, Sigma = Sig, p = p, r = r, u = u)
  output
}
