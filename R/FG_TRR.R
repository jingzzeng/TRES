#' @export
#'

# This function gives FG estimation of tensor response regression
FG_TRR <- function(Yn, Xn, Gamma_init) {

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
    YsnYsn <- ttt(Ysn, Ysn, dims=idx)@data*idxprod
    U <- YsnYsn - M
    Gamma1[[i]] <- OptStiefelGBB(Gamma_init[[i]], opts=NULL, FGfun, M, U)$X
    PGamma[[i]] <- Gamma1[[i]] %*% t(Gamma1[[i]])
  }
  tp <- ttl(Yn, PGamma, ms=1:m)
  Bhat <- ttm(tp, Xn_inv, m+1)

  return(list(Bhat=Bhat, Gamma_hat=Gamma1))
}
