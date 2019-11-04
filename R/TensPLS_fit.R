#' @title Tensor envelope partial least squares (PLS) regression
#' @description Estimates the factor matrix \eqn{W_k, k=1,\cdots,m} in tensor PLS algorithm for tensor predictor regression (TPR), see Zhang, X., & Li, L. (2017).
#'
#' @param  Xn A predictor tensor of dimension \eqn{p_1\times \cdots \times p_m}.
#' @param Yn The response vector of dimension \eqn{r}.
#' @param u The dimension of envelope subspace, \eqn{u=(u_1,\cdots,u_m)}.
#' @param SigX A matrix lists \eqn{\boldsymbol{\Sigma}_k, k=1,\cdots, m}, which determines the estimation of covariance matrix \eqn{\boldsymbol{\Sigma}=\boldsymbol{\Sigma}_m \otimes \cdots \otimes \boldsymbol{\Sigma}_1}.
#'
#' @return
#' \describe{
#'   \item{Gamma}{The estimation of factor matrix \eqn{W_k, k=1,\cdots,m}.}
#'   \item{PGamma}{The projection matrix \eqn{W_k(W_k'\boldsymbol{\Sigma}_k W_k)^{-1}W_k'\boldsymbol{\Sigma}_k, k=1,\cdots,m}.}
#' }
#'
#' @references Zhang, X., & Li, L. (2017). Tensor Envelope Partial Least-Squares Regression. Technometrics, 59(4), 426-436.

#' @name TensPLS_fit-deprecated
#' @usage TensPLS_fit(Xn, Yn, SigX, u)
#' @seealso \code{\link{TRES-deprecated}}
#' @keywords internal
NULL

#' @rdname TRES-deprecated
#' @section \code{TensPLS_fit}:
#' For \code{TensPLS_fit}, use \code{\link{TPR.fit}} with \code{method = "PLS"}.


#' @import rTensor
#' @importFrom pracma kron sqrtm
#' @importFrom stats cov
TensPLS_fit <- function(Xn, Yn, SigX, u) {
  ss <- dim(Xn)
  len <- length(ss)
  n <- ss[len]
  p <- ss[1:(len-1)]
  m <- length(p)


  ##center the data
  muy <- as.matrix(apply(Yn, 1, mean))
  Yn <- Yn - muy[, rep(1, n)]
  mux <- apply(Xn@data, c(1:m), mean)
  ttmp <- lapply(1:n, function(x) mux)
  ttmp <- array(unlist(ttmp), c(p, n))
  ttmp2 <- Xn@data - ttmp

  Xn <- rTensor::as.tensor(ttmp2)


  SigY <- (n-1)*cov(t(Yn))/n

  Sinvhalf <- NULL
  for (i in 1:m) {
    Sinvhalf[[i]] <- pracma::sqrtm(SigX[[i]])$Binv
  }

  Sinvhalf[[m+1]] <- pracma::sqrtm(SigY)$Binv

  C <- rTensor::ttm(Xn, Yn, m+1)/n

  Gamma <- PGamma <- NULL
  for (i in 1:m) {
    M <- SigX[[i]]
    idx <- c(1:(m+1))[-i]

    Ck <- rTensor::ttl(C, Sinvhalf[idx], ms = idx)

    U <- rTensor::unfold(Ck, row_idx = i, col_idx = idx)@data
    Uk <- U %*% t(U)

    Gamma[[i]] <- EnvMU(M, Uk, u[i])
    tmp3 <- t(Gamma[[i]]) %*% SigX[[i]] %*% Gamma[[i]]
    PGamma[[i]] <- Gamma[[i]] %*% solve(tmp3) %*% t(Gamma[[i]]) %*% SigX[[i]]
  }

  return(list(Gamma = Gamma, PGamma = PGamma))
}
