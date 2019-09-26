
#' @export
#' @importFrom stats runif


## data generation ##
MenvU_sim <- function(n, p, u, Omega=NULL, Omega0=NULL, Phi=NULL) {

  ##randomly generate a semi-orthogonal p-by-u basis matrix (Gamma) for the
  ##envelope and its orthogonal completion (Gamma0) of dimension p-by-(p-u)
  Gamma <- matrix(runif(p*u), p, u)

  ###make Gamma semi-orthogonal
  Gamma <- qr.Q(qr(Gamma))
  Gamma0 <- qr.Q(qr(Gamma),complete=TRUE)[,(u+1):p]

  ## randomly generated symmetric positive definite matrices, M and U, to have
  ## an exact u-dimensional envelope structure

  if (is.null(Phi)) {
    Phi <- matrix(runif(u^2), u, u)
    Phi <- Phi %*% t(Phi)
  }

  if (is.null(Omega)) {
      Omega <- matrix(runif(u^2), u, u)
      Omega <- Omega %*% t(Omega)
  }

  if (is.null(Omega0)) {
    Omega0 <- matrix(runif((p-u)^2),p-u,p-u)
    Omega0 <- Omega0 %*% t(Omega0)
  }
  M <- Gamma %*% Omega %*% t(Gamma) + Gamma0 %*% Omega0 %*% t(Gamma0)
  U <- Gamma %*% Phi %*% t(Gamma)

  X <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = M)
  Y <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = U)
  Mhat <- (t(X) %*% X)/n
  Uhat <- (t(Y) %*% Y)/n

  return(list(Mhat=Mhat, Uhat=Uhat, Gamma=Gamma))
}
