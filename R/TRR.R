#' @export

# This function gives all the estimation of tensor response regression
TRR <- function(Xn, Yn, method=c('standard', 'FG', '1D', 'ECD', 'PLS'), u=NULL, Gamma_init=NULL) {
  cl <- match.call()
  method <- match.arg(method)
  if(!is.matrix(Xn)){
    if(is.vector(Xn)){
      Xn <- t(as.matrix(Xn))
    }
    else stop("Xn should be vector or matrix.")
  }
  if(!inherits(Yn, "Tensor")){
    if(is.matrix(Yn) || inherits(Yn, "array")){
      Yn <- as.tensor(Yn)
    }
    else stop("Yn should be matrix, array or Tensor.")
  }
  Xn_old <- Xn
  Yn_old <- Yn
  ss <- dim(Yn)
  len <- length(ss)
  n <- ss[len]
  if(n != dim(Xn)[2]){stop("Unmatched dimension.")}
  r <- ss[1:(len-1)]
  m <- length(r)
  prodr <- prod(r)
  p <- dim(Xn)[1]
  ##center the data
  mux <- as.matrix(apply(Xn, 1, mean))
  Xn <- Xn-mux[, rep(1, n)]
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

  if(method == "standard") {
    Bhat = Btil
    Gamma1 = NULL
  }else {
    if(missing(u)){stop("A user-defined u is required.")}
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

      if (method == "1D") {
        Gamma1[[i]] <- OptimballGBB1D(M, U, u[i], opts=NULL)
      }else if(method == "ECD") {
        Gamma1[[i]] <- ECD(M, U, u[i])
      }else if(method == "PLS") {
        Gamma1[[i]] <- EnvMU(M, U, u[i])
      }else if(method == "FG"){
        if(missing(Gamma_init)){
          init <-  OptimballGBB1D(M, U, u[i], opts=NULL)
        }else{
          init <- Gamma_init[[i]]
        }
        Gamma1[[i]] <- OptStiefelGBB(init, opts=NULL, FGfun, M, U)$X
      }
      PGamma[[i]] <- Gamma1[[i]] %*% t(Gamma1[[i]])
    }
    tp <- rTensor::ttl(Yn, PGamma, ms=1:m)
    Bhat <- rTensor::ttm(tp, Xn_inv, m+1)

  }
  m <- Bhat@num_modes
  fitted.values <- rTensor::ttm(Bhat, t(Xn_old), m)
  residuals <- Yn_old - fitted.values
  output <- list(Xn=Xn_old, Yn=Yn_old, method = method, coefficients=Bhat, Gamma_hat=Gamma1, Sig=Sig, fitted.values = fitted.values, residuals = residuals)
  class(output) <- "Tenv"
  output$call <- cl
  output
}
