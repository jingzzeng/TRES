#' @export
#'

# This function gives all the estimation of tensor predictor regression
# The tensor predictor should be 2-dimensional or 3-dimensional

TPR.fit <- function(Xn, Yn, method=c('standard', 'FG', '1D', 'ECD', 'PLS'), u=NULL, Gamma_init=NULL) {
  cl <- match.call()
  method <- match.arg(method)
  if(!is.matrix(Yn)){
    if(is.vector(Yn)){
      Yn <- t(as.matrix(Xn))
    }
    else stop("Yn should be vector or matrix.")
  }
  if(!inherits(Xn, "Tensor")){
    if(is.matrix(Xn) || inherits(Xn, "array")){
      Xn <- as.tensor(Xn)
    }
    else stop("Xn should be matrix, array or Tensor.")
  }
  Xn_old <- Xn
  Yn_old <- Yn
  method <- match.arg(method)
  ss <- dim(Xn)
  len <- length(ss)
  n <- ss[len]
  if(n != dim(Yn)[2]){stop("Unmatched dimension.")}
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
  lambda <- res$lambda
  Sigx <- res$S
  Sigx[[1]] <- lambda*Sigx[[1]]


  if(method == "standard") {
    if(length(dim(Xn))==4){
      tmp6 <- pracma::kron(MASS::ginv(Sigx[[3]]), MASS::ginv(Sigx[[2]]))
      Btil <- pracma::kron(tmp6, MASS::ginv(Sigx[[1]])) %*% vecXn %*% t(Yn)/n
    }else if(length(dim(Xn))==3) {
      tmp6 <- pracma::kron(MASS::ginv(Sigx[[2]]), MASS::ginv(Sigx[[1]]))
      Btil <- tmp6 %*% vecXn %*% t(Yn)/n
    }else if(length(dim(Xn))==2) {
      Btil <- MASS::ginv(Sigx[[1]])%*%vecXn%*%t(Yn)/n
    }
    Btil <- array(Btil, c(p, r))
    Bhat <- Btil
    Gamma1 <- NULL
  }else{
    if(missing(u)){stop("A user-defined u is required.")}
    if(method == "PLS") {
      res_PLS <- TensPLS_fit(Xn, Yn, Sigx, u)
      Gamma1 <- res_PLS$Gamma; PGamma <- res_PLS$PGamma

      if(length(dim(Xn))==4) {
        tmp7 <- pracma::kron(PGamma[[2]], PGamma[[1]])
        Bhat_pls <- pracma::kron(PGamma[[3]], tmp7) %*% vecXn %*% t(Yn)/n
      }else if(length(dim(Xn))==3) {
        tmp7 <- pracma::kron(PGamma[[2]], PGamma[[1]])
        Bhat_pls <- tmp7 %*% vecXn %*% t(Yn)/n
      }else if(length(dim(Xn))==2) {
        Bhat_pls <- PGamma[[1]]%*%vecXn%*%t(Yn)/n
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
        idx <- c(1:(m+1))[-i]
        Ck <- ttl(C, Sinvhalf[idx], ms = idx)
        U <- unfold(Ck, row_idx = i, col_idx = idx)@data
        Uk <- U %*% t(U)
        Gamma1[[i]] <- OptimballGBB1D(Sigx[[i]], Uk, u[i])
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

    if(method=='FG'){
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
        if(missing(Gamma_init)){
          init <-  OptimballGBB1D(Sigx[[i]], Uk, u[i], opts=NULL)
        }else{
          init <- Gamma_init[[i]]
        }
        Gamma1[[i]] <- OptStiefelGBB(init, opts=NULL, FGfun, Sigx[[i]], Uk)$Gamma
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
    }
  }

  Bhat <- as.tensor(Bhat)
  tp1 <- matrix(Bhat@data, nrow = c(prod(p)))
  tp2 <- matrix(Xn_old@data, c(prod(p), n))
  fitted.values <- t(tp1) %*% tp2
  residuals <- Yn_old - fitted.values
  output <- list(Xn=Xn_old, Yn=Yn_old, method = method, coefficients=Bhat, Gamma=Gamma1, Sigx=Sigx, fitted.values = fitted.values, residuals=residuals)
  class(output) <- "Tenv"
  output$call <- cl
  output
}
