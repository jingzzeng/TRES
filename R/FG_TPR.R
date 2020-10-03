#' @title Envelope estimation of tensor predictor regression (TPR) with the full Grassmannian optimization
#' @description This function is used for envelope estimation of tensor predictor regression with the full Grassmannian (FG) optimization.
#'
#' @param Xn The predictor tensor instance or dimension \eqn{p_1\times p_2\times\cdots\times p_m \times n}, where \eqn{n} is the sample size.
#' @param Yn The predictor matrix of dimension \eqn{r \times n}.
#' @param Gamma_init The initial estimation of envelope subspace basis, can be derived from \code{TPR.fit}.
#'
#' @return
#' \describe{
#'  \item{Bhat}{The estimation of regression coefficient tensor.}
#'  \item{Gamma_hat}{The FG estimation of envelope subspace basis.}
#' }
#' @examples
#' rm(list = ls())
#' p <- c(10, 10, 10)
#' u <- c(1, 1, 1)
#' m <- 3; r <- 5; n <- 200
#' eta <- array(runif(prod(u,r)), c(u,r))
#' eta <- rTensor::as.tensor(eta)
#'
#' Gamma <- Gamma0 <- Omega <- Omega0 <- Sig <- Sigsqrtm <- NULL
#' for(i in 1:m) {
#'   tmp <- matrix(runif(p[i]*u[i]), p[i], u[i])
#'   Gamma[[i]] <- qr.Q(qr(tmp))
#'   Gamma0[[i]] <- qr.Q(qr(tmp), complete=TRUE)[, (u[i]+1):p[i]]
#'   Omega[[i]] <- diag(u[i])
#'   Omega0[[i]] <- 0.01*diag(p[i]-u[i])
#'   Sig[[i]] <- Gamma[[i]] %*% Omega[[i]] %*% t(Gamma[[i]])+
#'     Gamma0[[i]] %*% Omega0[[i]] %*% t(Gamma0[[i]])
#'   Sig[[i]] <- 2*Sig[[i]]/norm(Sig[[i]], type="F")
#'   Sigsqrtm[[i]] <- pracma::sqrtm(Sig[[i]])$B
#' }
#'
#' B <- rTensor::ttl(eta,Gamma, ms = c(1:m))
#' A <- matrix(runif(r^2), r, r)
#' SigY <- A %*% t(A)
#' SigY <- SigY/norm(SigY, type="F")

#' ##generate data
#' Epsilon <- MASS::mvrnorm(n, mu=rep(0, r), Sigma=SigY)
#' tmp2 <- array(rnorm(prod(p, n)), c(p, n))
#' Xn <- rTensor::as.tensor(tmp2)
#' Xn <- rTensor::ttl(Xn, Sigsqrtm, ms = c(1:m))
#' vecXn <- matrix(Xn@data, prod(p), n)
#' Y_tmp <- matrix(NA, r, n)
#' tmp <- array(NA, c(p, r))
#' for (j in 1:n) {
#'   for (s in 1:r) {
#'     tmp[, , , s] <-  B@data[, , , s]*Xn@data[, , , j]
#'   }
#'   Y_tmp[, j] <- apply(tmp, 4, sum)
#' }
#' Yn <-  Y_tmp + t(Epsilon)
#'
#' # use the result of 1D method as the initial value
#' res_1D = TPR.fit(Xn, Yn, u, method="1D")
#' \dontrun{
#'   res_FG = FG_TPR(Xn, Yn, Gamma_init=res_1D$Gamma_hat)
#' }
#'
#' @name FG_TPR-deprecated
#' @usage FG_TPR(Xn, Yn, Gamma_init)
#' @seealso \code{\link{TRES-deprecated}}
#' @keywords internal
NULL

#' @rdname TRES-deprecated
#' @section \code{FG_TPR}:
#' For \code{FG_TPR}, use \code{\link{TPR.fit}} with \code{method = "FG"}.
#'
#' @export
#' @importFrom stats cov

# This function gives FG estimation of tensor predictor regression
FG_TPR <- function(Xn, Yn, Gamma_init){
  .Deprecated("TPR.fit", package = "TRES")
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

    Gamma1[[i]] <- OptStiefelGBB(Gamma_init[[i]], opts=NULL, FGfun, Sigx[[i]], Uk)$Gamma

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
