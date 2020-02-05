#' @export
EnvMU <- function(M, U, u) {
  dimM <- dim(M)
  dimU <- dim(U)
  p <- dimM[1]

  if (dimM[1] != dimM[2] & dimU[1] != dimU[2]) stop("M and U should be square matrices.")
  if (dimM[1] != dimU[1]) stop("M and U should have the same dimension.")
  if (qr(M)$rank < p) stop("M should be positive definite.")
  if (u > p & u < 0) stop("u should be between 0 and p.")
  if (u == p) {return(Gamma = diag(p))}else {
    W <- matrix(0, p, (u+1))
    for (k in 1:u) {
      Wk <- W[, 1:k]
      Ek <- M %*% Wk
      temp <- crossprod(Ek)
      QEK <- diag(p) - Ek %*% MASS::ginv(temp) %*% t(Ek)
      W[, (k+1)] <- Re(eigen(QEK %*% U %*% QEK)$vectors[, 1])
    }
    Gamma <- qr.Q(qr(W[, 2:(u+1)]))
    return(Gamma)
  }
}
