#' @export
#' @import rTensor
#' @import MASS
#' @importFrom stats pt
Tenv_Pval <- function(Xn, Yn, B_est){
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

  Yn <- as.tensor(tmp2)
  Xn_inv <- ginv(Xn %*% t(Xn)) %*% Xn
  Btil <- ttm(Yn, Xn_inv, m+1)
  En <- Yn - ttm(Btil, t(Xn), m+1)

  res <- kroncov(En)
  lambda <- res$lambda
  Sig <- res$S

  Sb <- ginv(cov(t(Xn)))
  if(min(dim(Sb))>1){Sb <- matrix(diag(Sb))}

  for(i in 1:m) {
    Sb <- kronecker(Sb, diag(Sig[[i]]))
  }

  Sb <- array(sqrt(lambda*Sb), c(r, p))
  Ttil <- (sqrt(n)*Btil@data)/Sb
  That <- (sqrt(n)*B_est@data)/Sb
  Ptil <- 2*(1-pt(abs(Ttil), n-1))
  Phat <- 2*(1-pt(abs(That), n-1))

  Sb <- as.tensor(Sb)
  Ptil <- as.tensor(Ptil)
  Phat <- as.tensor(Phat)

  return(list(P_OLS=Ptil, P_val=Phat, se=Sb))
}
