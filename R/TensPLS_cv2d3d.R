#' @export
#' @import rTensor
#' @importFrom pracma kron sqrtm
TensPLS_cv2d3d <- function(x, y, maxdim=10, nfolds=5) {
  if(!is.matrix(y)){
    if(is.vector(y)){
      y <- t(as.matrix(y))
    }
    else stop("y should be vector or matrix.")
  }
  if(!inherits(x, "Tensor")){
    if(is.matrix(x) || inherits(x, "array")){
      x <- as.tensor(x)
    }
    else stop("x should be matrix, array or Tensor.")
  }
  ss <- dim(x)
  len <- length(ss)
  n <- ss[len]
  p <- ss[1:(len-1)]
  r <- dim(y)[1]
  m <- length(p)
  vecX0 <- matrix(x@data, c(prod(p), r))
  idx <- sample(1:n, n, replace = FALSE)
  Ntest <- floor(n/nfolds)
  Ntrain <- n - Ntest
  cv_sse <- matrix(0, c(maxdim, 1))

  for (i in 1:nfolds) {
    testid <- c(1:Ntest) + (i-1)*Ntest
    testid <- idx[testid]
    Ytrain <- y
    vecXtrain <- vecX0
    vecXtrain_cv <- vecXtrain[, -testid]
    Ytrain_cv <- matrix(Ytrain[, -testid], 1, Ntrain)
    mu_vecX <- as.matrix(apply(vecXtrain_cv, 1, mean))
    mu_Y <- as.matrix(apply(Ytrain_cv, 1, mean))
    Ytrain_cv <- Ytrain_cv - mu_Y[, rep(1, Ntrain)]
    vecXtrain_cv <- vecXtrain_cv - mu_vecX[, rep(1, Ntrain)]

    tp <- array(vecXtrain_cv, c(p, Ntrain))
    Xtrain <- rTensor::as.tensor(tp)
    Ytest <- matrix(y[, testid], 1, Ntest)
    vecXtest <- vecX0[, testid]
    Ytest <- Ytest - mu_Y[, rep(1, Ntest)]
    vecXtest <- vecXtest - mu_vecX[, rep(1, Ntest)]
    tp2 <- array(vecXtest, c(p, Ntest))
    Xtest <- rTensor::as.tensor(tp2)
    ##Fit TPLS ##
    res <- kroncov(Xtrain)
    lambda <- res$lambda
    SigX <- res$S
    SigX[[1]] <- lambda*SigX[[1]]
    res <- TensPLS_fit(Xtrain, Ytrain_cv, SigX, (maxdim*matrix(1, m, 1)))
    Gamma <- res$Gamma
    PGamma <- res$PGamma

    Ghat <- NULL
    for (k in 1:maxdim) {
      for (j in 1:m) {
         if (k==p[j]) {
           Ghat[[j]] <- diag(p[j])
          }else {
            Gtmp <- Gamma[[j]]
            Ghat[[j]] <- Gtmp[, 1:k]
            tmp <- t(Ghat[[j]]) %*% SigX[[j]] %*% Ghat[[j]]
            PGamma[[j]] <- Ghat[[j]] %*% chol2inv(chol(tmp)) %*% t(Ghat[[j]])
          }
          if (m == 2) {
            tmp2 <- pracma::kron(PGamma[[2]], PGamma[[1]])
            Bhat_pls <- tmp2 %*% vecXtrain_cv %*% t(Ytrain_cv)/Ntrain
          }else if(m == 3){
            tmp2 <- pracma::kron(PGamma[[2]], PGamma[[1]])
            Bhat_pls <- pracma::kron(PGamma[[3]], tmp2) %*% vecXtrain_cv %*% t(Ytrain_cv)/Ntrain
          }else if(m == 1){
            Bhat_pls <- PGamma[[1]] %*% vecXtrain_cv %*% t(Ytrain_cv)/Ntrain
          }
      }
      ehat <- crossprod(Bhat_pls, vecXtest) - Ytest
      cv_sse[k] <- cv_sse[k] + sum(diag(tcrossprod(ehat)))
    }
  }
  mincv <- min(Re(cv_sse)); u <- which.min(Re(cv_sse))
  return(list(mincv=mincv, u=u))
}
