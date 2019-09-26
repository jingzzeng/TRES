#' @export
#'

# This function gives FG estimation of tensor predictor regression
FG_TPR <- function(Yn, Xn, Gamma_init) {
  ss <- dim(Xn)
  len <- length(ss)
  n <- ss[len]
  p <- ss[1:(len-1)]
  m <- length(p)
  r <- dim(Yn)[1]

  ##center the data
  muy <- as.matrix(apply(Yn, 1, mean))
  Yn <- Yn - muy[, rep(1, n)]
  mux <- apply(Xn@data, c(1:m), mean)
  ttmp <- lapply(1:n, function(x) mux)
  ttmp <- array(unlist(ttmp), c(p, n))
  ttmp2 <- Xn@data - ttmp

  Xn <- rTensor::as.tensor(ttmp2)
  vecXn <- matrix(Xn@data, prod(p), n)

  res <- kroncov(Xn)
  Sigx <- res$S; lambda <- res$lambda
  Sigx[[1]] <- lambda*Sigx[[1]]

  Sinvhalf <- NULL
  for (i in 1:m) {
    Sinvhalf[[i]] <- pracma::sqrtm(Sigx[[i]])$Binv
  }
  SigY <- (n-1)*cov(t(Yn))/n
  Sinvhalf[[m+1]] <- pracma::sqrtm(SigY)$Binv

  C <- ttm(Xn, Yn, m+1)/n
  Gamma1 <- PGamma <- NULL
  for (i in 1:m) {

    idx <- c(1:(m+1))[-i]

    Ck <- ttl(C, Sinvhalf[idx], ms = idx)

    U <- unfold(Ck, row_idx = i, col_idx = idx)@data
    idxprod <- (p[i]/r)/prod(p)

    Uk <- idxprod*U %*% t(U)

    Gamma1[[i]] <- OptStiefelGBB(Gamma_init[[i]], opts=NULL, FGfun, Sigx[[i]], Uk)$X

    tmp8 <- t(Gamma1[[i]]) %*% Sigx[[i]] %*% Gamma1[[i]]
    PGamma[[i]] <- Gamma1[[i]] %*% solve(tmp8) %*% t(Gamma1[[i]]) %*% Sigx[[i]]
  }

  if(length(dim(Xn))==4) {
    tmp9 <- pracma::kron(PGamma[[2]], PGamma[[1]])
    Bhat_env <- pracma::kron(PGamma[[3]], tmp9) %*% vecXn %*% t(Yn)/n
  }else if(length(dim(Xn))==3) {
    tmp9 <- pracma::kron(PGamma[[2]], PGamma[[1]])
    Bhat_env <- tmp9 %*% vecXn %*% t(Yn)/n
  }else if(length(dim(Xn))==2) {
    Bhat_env <- PGamma[[1]] %*% vecXn %*% t(Yn)/n
  }
  Bhat <- array(Bhat_env, c(p, r))

  return(list(Bhat=Bhat, Gamma_hat=Gamma1))

}
