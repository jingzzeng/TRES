#' @export
#' @import rTensor
#' @importFrom stats pt
Tenv_Pval <- function(x, y, Bhat){
  ss <- dim(y)
  len <- length(ss)
  n <- ss[len]
  r <- ss[1:(len-1)]
  m <- length(r)
  prodr <- prod(r)
  p <- dim(x)[1]
  mux <- as.matrix(apply(x, 1, mean))
  x <- x -mux[, rep(1, n)]
  muy <- apply(y@data, c(1:m), mean)

  tmp1 <- lapply(1:n, function(x) muy)
  tmp1 <- array(unlist(tmp1), c(r, n))
  tmp2 <- y@data-tmp1

  y <- as.tensor(tmp2)
  x_inv <- chol2inv(chol(tcrossprod(x))) %*% x
  Btil <- ttm(y, x_inv, m+1)
  En <- y - ttm(Btil, t(x), m+1)

  res <- kroncov(En)
  lambda <- res$lambda
  Sig <- res$S

  Sb <- chol2inv(chol(cov(t(x))))
  if(min(dim(Sb))>1){Sb <- matrix(diag(Sb))}

  for(i in 1:m) {
    Sb <- kronecker(Sb, diag(Sig[[i]]))
  }

  Sb <- array(sqrt(lambda*Sb), c(r, p))
  Ttil <- (sqrt(n)*Btil@data)/Sb
  That <- (sqrt(n)*Bhat@data)/Sb

  Ptil <- 2*(1-pt(abs(Ttil), n-1))
  Phat <- 2*(1-pt(abs(That), n-1))

  Sb <- as.tensor(Sb)
  Ptil <- as.tensor(Ptil)
  Phat <- as.tensor(Phat)

  return(list(p_ols=Ptil, p_val=Phat, se=Sb))
}
