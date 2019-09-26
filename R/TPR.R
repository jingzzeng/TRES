#' @export
#'

# This function gives all the estimation of tensor predictor regression
# The tensor predictor should be 2-dimensional or 3-dimensional

TPR <- function(Yn, Xn, u, method) {
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

  ## fit ls
  res <- kroncov(Xn)
  Sigx <- res$S; lambda <- res$lambda
  Sigx[[1]] <- lambda*Sigx[[1]]

  if(length(dim(Xn))==4){
    tmp6 <- pracma::kron(ginv(Sigx[[3]]), ginv(Sigx[[2]]))
    Btil <- pracma::kron(tmp6, ginv(Sigx[[1]])) %*% vecXn %*% t(Yn)/n
  }else if(length(dim(Xn))==3) {
    tmp6 <- pracma::kron(ginv(Sigx[[2]]), ginv(Sigx[[1]]))
    Btil <- tmp6 %*% vecXn %*% t(Yn)/n
  }else if(length(dim(Xn))==2) {
    Btil <- ginv(Sigx[[1]])%*%vecXn%*%t(Yn)/n
  }
  Btil <- array(Btil, c(p, r))

  if(method == "standard") {
    Bhat = Btil
    Gamma1 = NULL
  }

  if(method == "PLS") {
    res_PLS <- TensPLS_fit(Xn, Yn, Sigx, u)
    Gamma1 <- res_PLS$Gamma; pghat <- res_PLS$PGamma

    if(length(dim(Xn))==4) {
      tmp7 <- pracma::kron(pghat[[2]], pghat[[1]])
      Bhat_pls <- pracma::kron(pghat[[3]], tmp7) %*% vecXn %*% t(Yn)/n
    }else if(length(dim(Xn))==3) {
      tmp7 <- pracma::kron(pghat[[2]], pghat[[1]])
      Bhat_pls <- tmp7 %*% vecXn %*% t(Yn)/n
    }else if(length(dim(Xn))==2) {
      Bhat_pls <- pghat[[1]]%*%vecXn%*%t(Yn)/n
    }
    Bhat <- array(Bhat_pls, c(p, r))
  }

  if(method == "1D") {
    Sinvhalf <- NULL
    for (i in 1:m) {
      Sinvhalf[[i]] <- pracma::sqrtm(Sigx[[i]])$Binv
    }
    SigY <- (n-1)*cov(t(Yn))/n
    Sinvhalf[[m+1]] <- pracma::sqrtm(SigY)$Binv

    C <- ttm(Xn, Yn, m+1)/n
    Gamma1 <- PGamma <- NULL
    for (i in 1:m) {
      M <- Sigx[[i]]
      idx <- c(1:(m+1))[-i]

      Ck <- ttl(C, Sinvhalf[idx], ms = idx)

      U <- unfold(Ck, row_idx = i, col_idx = idx)@data


      Uk <- U %*% t(U)

      Gamma1[[i]] <- OptimballGBB1D(M, Uk, u[i])

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
      Bhat_env <- PGamma[[1]] %*%vecXn%*%t(Yn)/n
    }
    Bhat <- array(Bhat_env, c(p, r))
  }

  if(method == "ECD") {
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


      Uk <- U %*% t(U)

      Gamma1[[i]] <- ECD(Sigx[[i]], Uk, u[i])

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
      Bhat_env <- PGamma[[1]]%*%vecXn%*%t(Yn)/n
    }
    Bhat <- array(Bhat_env, c(p, r))
  }
  return(list(Bhat=Bhat, Gamma_hat=Gamma1, Sigx=Sigx))
}
