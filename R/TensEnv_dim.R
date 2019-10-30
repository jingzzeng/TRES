#' @export
#' @import rTensor
#' @import MASS
#' @importFrom pracma kron sqrtm

TensEnv_dim <- function(Xn, Yn, multiD=1, maxdim=10, opts=NULL){
  ss <- dim(Yn)
  len <- length(ss)
  n <- ss[len]
  r <- ss[1:(len-1)]
  m <- length(r)
  prodr <- prod(r)
  p <- dim(Xn)[1]
  multiD <- p
  u <- r

  Xn_inv <- MASS::ginv(Xn %*% t(Xn)) %*% Xn
  Btil <- rTensor::ttm(Yn, Xn_inv, m+1)
  En <- Yn - rTensor::ttm(Btil, t(Xn), m+1)

  res <- kroncov(En)
  lambda <- res$lambda
  Sig <- res$S

  Sinvhalf <- NULL
  for (i in 1:m) {
    Sinvhalf[[i]] <- pracma::sqrtm(Sig[[i]])$Binv
  }


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
    res <- ballGBB1D_bic(M, U, n, multiD, maxdim, opts)
    u[i] <- res$u
  }
  return(u)
}
