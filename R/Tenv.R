#' @title Tensor response envelope estimator
#' @description This function gives the tensor envelope estimator for tensor response regression (TRR).
#'
#' @param  Xn The predictor matrix of dimension \eqn{p \times n}.
#' @param Yn The response tensor instance \eqn{ r_1\times r_2\times \cdots \times r_m \times n}, where \eqn{n} is the sample size.
#' @param u The dimension of envelope subspace. \eqn{u=(u_1,\cdots,u_m)}.
#' @param opts The option structure for Feasi. See function \code{OptimballGBB1D}.
#'
#' @return
#'   \item{Btil}{The ordinary least square estimator (OLS).}
#'   \item{Bhat}{The envelope based estimator.}
#'   \item{PGamma}{The projection matrix onto envelope subspace.}
#'
#' @references Li, L., & Zhang, X. (2017). Parsimonious tensor response regression. Journal of the American Statistical Association, 112(519), 1131-1146.

#' @name Tenv-deprecated
#' @usage Tenv(Xn, Yn, u, opts=NULL)
#' @seealso \code{\link{TRES-deprecated}}
#' @keywords internal
NULL

#' @rdname TRES-deprecated
#' @section \code{Tenv}:
#' For \code{Tenv}, use \code{\link{TRR.fit}} with \code{method = "1D"}.
#'
#' @export
#' @import rTensor
#' @import MASS

Tenv <- function(Xn, Yn, u, opts=NULL){
  .Deprecated("TRR.fit", package = "TRES")
  ss <- dim(Yn)
  len <- length(ss)
  n <- ss[len]
  r <- ss[1:(len-1)]
  m <- length(r)
  prodr <- prod(r)
  p <- dim(Xn)[1]
  mux <- as.matrix(apply(Xn, 1, mean))
  Xn <- Xn -mux[, rep(1, n)]
  muy <- apply(Yn@data, c(1:m), mean)

  tmp1 <- lapply(1:n, function(x) muy)
  tmp1 <- array(unlist(tmp1), c(r, n))
  tmp2 <- Yn@data-tmp1

  Yn <- rTensor::as.tensor(tmp2)
  Xn_inv <- MASS::ginv(tcrossprod(Xn)) %*% Xn
  Btil <- rTensor::ttm(Yn, Xn_inv, m+1)
  En <- Yn - rTensor::ttm(Btil, t(Xn), m+1)

  res <- kroncov(En)
  lambda <- res$lambda
  Sig <- res$S

  Sinvhalf <- NULL
  for (i in 1:m) {
    Sinvhalf[[i]] <- sqrtm(Sig[[i]])$Binv
  }
  Gamma1 <- PGamma <- NULL
  for (i in 1:m) {
    M <- lambda*Sig[[i]]
    idx <-  c(1:(m+1))[-i]
    len <- length(idx)
    if (len > 1) {
      Ysn <- rTensor::ttl(Yn, Sinvhalf[c(idx[1:(len-1)])], ms=idx[1:(len-1)])
    }else {
      Ysn <- rTensor::ttl(Yn, Sinvhalf, ms=1)
    }
    idxprod <- (r[i]/n)/prodr
    YsnYsn <- ttt(Ysn, Ysn, ms=idx)@data*idxprod
    U <- YsnYsn - M
    Gamma1[[i]] <- OptimballGBB1D(M, U, u[i], opts)
    PGamma[[i]] <- tcrossprod(Gamma1[[i]])
  }
  tp <- ttl(Yn, PGamma, ms=1:m)
  Bhat <- ttm(tp, Xn_inv, m+1)
  return(list(Btil = Btil, Bhat = Bhat, PGamma = PGamma))
}
