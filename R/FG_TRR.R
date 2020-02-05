#' @title Envelope estimation of tensor response regression (TRR) with the full Grassmannian optimization
#' @description This function is used for envelope estimation of tensor response regression with the full Grassmannian (FG) optimization.
#'
#' @param Xn The predictor matrix of dimension \eqn{p \times n}.
#' @param Yn The response tensor instance or dimension \eqn{r_1\times r_2\times\cdots\times r_m \times n}, where \eqn{n} is the sample size.
#' @param Gamma_init The initial estimation of envelope subspace basis, can be derived from \code{TRR.fit}.
#'
#' @return
#' \describe{
#'  \item{Bhat}{The estimation of regression coefficient tensor.}
#'  \item{Gamma_hat}{The FG estimation of envelope subspace basis.}
#' }
#' @examples
#' rm(list=ls())
#'
#' r <- c(10, 10, 10)
#' m <- length(r)
#' u <- c(2, 2, 2)
#' p <- 5
#' n <- 100
#'
#' set.seed(1)
#' eta <- array(runif(prod(u)*p), c(u, p))
#' eta <- rTensor::as.tensor(eta)
#'
#' Gamma <- Gamma0 <- Omega <- Omega0 <- Sig <- Sigsqrtm <- NULL
#' for (i in 1:m){
#'   tmp <- matrix(runif(r[i]*u[i]), r[i], u[i])
#'   Gamma[[i]] <- qr.Q(qr(tmp))
#'   Gamma0[[i]] <- qr.Q(qr(Gamma[[i]]),complete=TRUE)[,(u[i]+1):r[i]]
#'
#'   A <- matrix(runif(u[i]^2), u[i], u[i])
#'   Omega[[i]] <- A %*% t(A)
#'   A <- matrix(runif((r[i]-u[i])^2), (r[i]-u[i]), (r[i]-u[i]))
#'   Omega0[[i]] <- A %*% t(A)
#'   Sig[[i]] <- Gamma[[i]] %*% Omega[[i]] %*% t(Gamma[[i]])+
#'     Gamma0[[i]] %*% Omega0[[i]] %*% t(Gamma0[[i]])
#'   Sig[[i]] <- 10*Sig[[i]]/norm(Sig[[i]], type="F")+0.01*diag(r[i])
#'   Sigsqrtm[[i]] <- pracma::sqrtm(Sig[[i]])$B
#' }
#' B <- rTensor::ttl(eta, Gamma, ms=1:m)
#' Xn <- matrix(rnorm(p*n), p, n)
#' Xn_inv <- MASS::ginv(Xn %*% t(Xn)) %*% Xn
#' Epsilon <- array(rnorm(prod(r)*n), c(r, n))
#' Epsilon <- rTensor::as.tensor(Epsilon)
#' Epsilon <- rTensor::ttl(Epsilon, Sigsqrtm, ms=1:m)
#' Yn <- Epsilon + rTensor::ttm(B, t(Xn), m+1)
#'
#' # use the result of 1D method as the initial value
#' res_1D = TRR.fit(Xn, Yn, u, method="1D")
#' \dontrun{
#'   res_FG = FG_TRR(Xn, Yn, Gamma_init=res_1D$Gamma_hat)
#' }

#' @name FG_TRR-deprecated
#' @usage FG_TRR(Xn, Yn, Gamma_init)
#' @seealso \code{\link{TRES-deprecated}}
#' @keywords internal
NULL

#' @rdname TRES-deprecated
#' @section \code{FG_TRR}:
#' For \code{FG_TRR}, use \code{\link{TRR.fit}} with \code{method = "FG"}.
#'
#' @export

FG_TRR <- function(Xn, Yn, Gamma_init) {
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
    YsnYsn <- ttt(Ysn, Ysn, ms=idx)@data*idxprod
    U <- YsnYsn - M
    Gamma1[[i]] <- OptStiefelGBB(Gamma_init[[i]], opts=NULL, FGfun, M, U)$Gamma
    PGamma[[i]] <- Gamma1[[i]] %*% t(Gamma1[[i]])
  }
  tp <- ttl(Yn, PGamma, ms=1:m)
  Bhat <- ttm(tp, Xn_inv, m+1)

  return(list(Bhat=Bhat, Gamma_hat=Gamma1))
}
