#' Generate simulation data for tensor predictor regression (TPR)
#'
#' This function is used to generate simulation data used in tensor prediction regression.
#'
#' The tensor predictor regression model is of the form,
#' \deqn{\mathbf{Y} = \mathbf{B}_{(m+1)}\mathrm{vec}(\mathbf{X}) +\boldsymbol{\varepsilon}}
#' where response \eqn{\mathbf{Y} \in R^{r}}, predictor \eqn{\mathbf{X} \in R^{p_1\times \cdots\times p_m}}, and
#' the error term is multivariate normal distributed. The predictor is tensor normal distributed,
#' \deqn{\mathbf{X}\sim TN(0;\boldsymbol{\Sigma}_1,\dots,\boldsymbol{\Sigma}_m)}
#' According to the tensor envelope structure, we have
#' \deqn{\mathbf{B} = [\Theta;\boldsymbol{\Gamma}_1,\ldots,\boldsymbol{\Gamma}_m,\mathbf{I}_p], \quad \mbox{for some } \boldsymbol{\Theta} \in R^{u_1\times\cdots\times u_m \times p}}
#' \deqn{\boldsymbol{\Sigma}_k = \boldsymbol{\Gamma}_k\boldsymbol{\Omega}_k\boldsymbol{\Gamma}_k^{T}+\boldsymbol{\Gamma}_{0k}\boldsymbol{\Omega}_{0k}\boldsymbol{\Gamma}_{0k}^\top,
#' \quad \mbox{for some } \boldsymbol{\Omega}_k, \boldsymbol{\Omega}_{0k},\ k=1,\ldots,m.}
#'
#' @note The length of p must match the length of u, and each element of u must be less than the corresponding element in p.
#'
#' @param p The dimension of predictor, a vector in the form of \eqn{(p_1,\cdots, p_m)}.
#' @param r The dimension of response, a scale.
#' @param u The structural dimension of envelopes at each mode, a vector with the same length as p.
#' @param n The sample size.
#' @return
#' \item{x}{The predictor of dimension \eqn{p_1\times \cdots\times p_m \times n}}
#' \item{y}{The response of dimension \eqn{r\times n}}
#' \item{Gamma}{The envelope subspace basis of dimension \eqn{p_k \times u_k, \ k=1,\ldots,m}}
#' \item{coefficients}{The tensor coefficients of dimension \eqn{p_1\times \cdots\times p_m \times r}}
#' \item{Sigma}{A lists of estimated covariance matrices at each mode for the tensor predictors, i.e., \eqn{\boldsymbol{\Sigma}_1,\dots,\boldsymbol{\Sigma}_m}}
#' \item{p, r, u}{The input \code{p,r,u}}
#'
#' @examples
#' p <- c(10, 10, 10)
#' u <- c(1, 1, 1)
#' r <- 5
#' n <- 200
#' dat <- TPR_sim(p = p, r = r, u = u, n = n)
#' x <- dat$x
#' y <- dat$y
#' fit_std <- TPR.fit(x, y, method="standard")
#'
#' @seealso \code{\link{TPR.fit}, \link{TRR_sim}}.
#' @references Zhang, X., Li, L. (2017). Tensor Envelope Partial Least-Squares Regression. Technometrics, 59(4), 426-436.

#' @export
#' @importFrom stats rnorm
TPR_sim <- function(p, r, u, n){
  if(length(p) != length(u)){stop("The length of p must match the length of u.")}
  if(any(p < u)){stop("u must less than p element-wise.")}
  m <- length(p)
  Gamma <- Gamma0 <- Omega <- Omega0 <- Sig <- Sigsqrtm <- NULL
  for(i in 1:m) {
    tmp <- matrix(runif(p[i]*u[i]), p[i], u[i])
    Gamma[[i]] <- qr.Q(qr(tmp))
    Gamma0[[i]] <- qr.Q(qr(tmp), complete=TRUE)[, (u[i]+1):p[i]]
    Omega[[i]] <- diag(u[i])
    Omega0[[i]] <- 0.01*diag(p[i]-u[i])
    Sig[[i]] <- Gamma[[i]] %*% Omega[[i]] %*% t(Gamma[[i]])+
      Gamma0[[i]] %*% Omega0[[i]] %*% t(Gamma0[[i]])
    Sig[[i]] <- 2*Sig[[i]]/norm(Sig[[i]], type="F")
    Sigsqrtm[[i]] <- pracma::sqrtm(Sig[[i]])$B
  }
  A <- matrix(runif(r^2), r, r)
  SigY <- tcrossprod(A)
  SigY <- SigY/norm(SigY, type="F")

  ##generate data
  Epsilon <- MASS::mvrnorm(n, mu=rep(0, r), Sigma=SigY)
  tmp2 <- array(rnorm(prod(p, n)), c(p, n))
  x <- rTensor::as.tensor(tmp2)
  x <- rTensor::ttl(x, Sigsqrtm, ms = c(1:m))
  vecx <- matrix(x@data, prod(p), n)
  eta <- array(runif(prod(u,r)), c(u,r))
  eta <- rTensor::as.tensor(eta)
  B <- rTensor::ttl(eta,Gamma, ms = c(1:m))
  vecB <- vecB <- matrix(B@data, prod(p), r)
  Y_tmp <- crossprod(vecB, vecx)
  y <- Y_tmp + t(Epsilon)

  output <- list(x = x, y = y, Gamma = Gamma, coefficients = B, Sigma = Sig, p = p, r = r, u = u)
  output
}
