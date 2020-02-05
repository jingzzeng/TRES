#' @export
#' @importFrom rTensor ttm ttl
#' @importFrom pracma kron sqrtm

TensEnv_dim <- function(x, y, multiD=1, maxdim=10, opts=NULL){
  if(!is.matrix(x)){
    if(is.vector(x)){
      x <- t(as.matrix(x))
    }
    else stop("x should be vector or matrix.")
  }
  if(!inherits(y, "Tensor")){
    if(is.matrix(y) || inherits(y, "array")){
      y <- as.tensor(y)
    }
    else stop("y should be matrix, array or Tensor.")
  }
  ss <- dim(y)
  len <- length(ss)
  n <- ss[len]
  r <- ss[1:(len-1)]
  m <- length(r)
  prodr <- prod(r)
  p <- dim(x)[1]
  multiD <- p
  u <- r

  x_inv <- chol2inv(chol(tcrossprod(x))) %*% x
  Btil <- ttm(y, x_inv, m+1)
  En <- y - ttm(Btil, t(x), m+1)

  res <- kroncov(En)
  lambda <- res$lambda
  Sig <- res$S

  Sinvhalf <- NULL
  for (i in 1:m) {
    Sinvhalf[[i]] <- sqrtm(Sig[[i]])$Binv
  }


  for (i in 1:m) {
    M <- lambda*Sig[[i]]
    idx <-  c(1:(m+1))[-i]
    len <- length(idx)
    if (len > 1) {
      Ysn <- ttl(y, Sinvhalf[c(idx[1:(len-1)])], ms=idx[1:(len-1)])
    }else {
      Ysn <- ttl(y, Sinvhalf, ms=1)
    }

    idxprod <- (r[i]/n)/prodr
    YsnYsn <- ttt(Ysn, Ysn, ms=idx)@data*idxprod
    U <- YsnYsn - M
    res <- ballGBB1D_bic(M, U, n, multiD, maxdim, opts)
    u[i] <- res$u
  }
  return(u)
}
