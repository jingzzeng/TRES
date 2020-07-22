#' @title Tensor envelope partial least squares (PLS) regression
#' @description Estimates the factor matrix \eqn{W_k, k=1,\cdots,m} in tensor PLS algorithm for tensor predictor regression (TPR), see Zhang, X., & Li, L. (2017).
#'
#' @param  x A predictor tensor of dimension \eqn{p_1\times \cdots \times p_m}.
#' @param y The response vector of dimension \eqn{r}.
#' @param u The dimension of envelope subspace, \eqn{u=(u_1,\cdots,u_m)}.
#' @param SigX A matrix lists \eqn{\boldsymbol{\Sigma}_k, k=1,\cdots, m}, which determines the estimation of covariance matrix \eqn{\boldsymbol{\Sigma}=\boldsymbol{\Sigma}_m \otimes \cdots \otimes \boldsymbol{\Sigma}_1}.
#'
#' @return
#'   \item{Gamma}{The estimation of factor matrix \eqn{W_k, k=1,\cdots,m}.}
#'   \item{PGamma}{The projection matrix \eqn{W_k(W_k'\boldsymbol{\Sigma}_k W_k)^{-1}W_k'\boldsymbol{\Sigma}_k, k=1,\cdots,m}.}
#'
#' @references Zhang, X., & Li, L. (2017). Tensor Envelope Partial Least-Squares Regression. Technometrics, 59(4), 426-436.

#' @name TensPLS_fit-deprecated
#' @usage TensPLS_fit(x, y, SigX, u)
#' @seealso \code{\link{TRES-deprecated}}
#' @keywords internal
NULL

#' @rdname TRES-deprecated
#' @section \code{TensPLS_fit}:
#' For \code{TensPLS_fit}, use \code{\link{TPR.fit}} with \code{method = "PLS"}.

#' @export
#' @import rTensor
#' @importFrom pracma kron sqrtm
#' @importFrom stats cov
TensPLS_fit <- function(x, y, SigX, u) {
  ss <- dim(x)
  len <- length(ss)
  n <- ss[len]
  p <- ss[1:(len-1)]
  m <- length(p)

  ##center the data
  muy <- as.matrix(apply(y, 1, mean))
  y <- y - muy[, rep(1, n)]
  mux <- apply(x@data, c(1:m), mean)
  ttmp <- lapply(1:n, function(x) mux)
  ttmp <- array(unlist(ttmp), c(p, n))
  ttmp2 <- x@data - ttmp

  x <- as.tensor(ttmp2)
  SigY <- (n-1)*cov(t(y))/n

  Sinvhalf <- NULL
  for (i in 1:m) {
    Sinvhalf[[i]] <- sqrtm(SigX[[i]])$Binv
  }

  Sinvhalf[[m+1]] <- sqrtm(SigY)$Binv

  C <- ttm(x, y, m+1)/n

  Gamma <- PGamma <- NULL
  for (i in 1:m) {
    M <- SigX[[i]]
    idx <- c(1:(m+1))[-i]

    Ck <- ttl(C, Sinvhalf[idx], ms = idx)

    U <- unfold(Ck, row_idx = i, col_idx = idx)@data
    Uk <- tcrossprod(U)

    Gamma[[i]] <- simplsMU(M, Uk, u[i])
    tmp3 <- t(Gamma[[i]]) %*% SigX[[i]] %*% Gamma[[i]]
    PGamma[[i]] <- Gamma[[i]] %*% chol2inv(chol(tmp3)) %*% t(Gamma[[i]]) %*% SigX[[i]]
  }

  return(list(Gamma = Gamma, PGamma = PGamma))
}
