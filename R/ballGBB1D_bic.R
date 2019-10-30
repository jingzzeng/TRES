# GBB1D bic for selecting envelope dimension
#' @export
ballGBB1D_bic <- function(M, U, n, multiD=1, maxdim=10, opts=NULL) {
  p <- dim(M)[2]
  Mnew <- M;
  Unew <- U;
  G <- matrix(0, p, maxdim)
  G0 <- diag(p)
  phi <- rep(0, p)
  for (k in 1:maxdim) {
    if (k == p) break;

    gk <- ballGBB1D(Mnew, Unew, opts)
    phi[k] <- n*(log(t(gk) %*% Mnew %*% gk)+ log(t(gk) %*% solve(Mnew+Unew) %*% gk))+
      log(n)*multiD
    G[, k] <- G0 %*% gk
    G0 <- qr.Q(qr(G[, 1:k]), complete=TRUE)[, (k+1):p]
    Mnew <- t(G0) %*% M %*% G0
    Unew <- t(G0) %*% U %*% G0
  }

  bic_val <- rep(0, maxdim)
  for (k in 1:maxdim) {
    bic_val[k] <- sum(phi[1:k])
  }
  return(list(u=which.min(bic_val), bicval=bic_val))
}
