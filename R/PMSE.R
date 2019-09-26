#' @export
PMSE <- function(Xn, Yn, Bhat) {
  ss <- dim(Xn)
  len <- length(ss)
  n <- ss[len]
  p <- ss[1:(len-1)]
  r <- dim(Yn)[1]
  m <- length(p)
  tp1 <- matrix(Bhat, c(prod(p), r))
  tp2 <- matrix(Xn, c(prod(p), n))
  Yhat <- t(tp1) %*% tp2
  Epsilon <- Yn - Yhat
  mse <- sum(diag(Epsilon %*% t(Epsilon)))/n
  return(list(mse=mse, Yhat=Yhat))
}