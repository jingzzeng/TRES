
PMSE <- function(Xn, Yn, Bhat) {
  if(inherits(Yn, "Tensor")){
    Yn <- Yn@data
  }else if(is.matrix(Yn)){
    # ss <- dim(Yn)
  }else if(is.vector(Yn)){
    Yn <- t(as.matrix(Yn))
  }else{
    stop("Yn must be tensor, matrix or vector.")
  }
  if(inherits(Xn, "Tensor")){
    Xn <- Xn@data
  }else if(is.matrix(Xn)){
    Xn <- as.tensor(Xn)
  }else if(is.vector(Xn)){
    Xn <- as.tensor(t(as.matrix(Xn)))
  }else{
    stop("Xn must be tensor, matrix or vector.")
  }
  if(inherits(Bhat, "Tensor")){Bhat <- Bhat@data}

  ss <- dim(Xn)
  len <- length(ss)
  n <- ss[len]
  if(n != dim(Yn)[2]){stop("Unmatched dimension.")}
  p <- ss[1:(len-1)]
  r <- dim(Yn)[1]
  m <- length(p)
  tp1 <- matrix(Bhat, c(prod(p), r))
  tp2 <- matrix(Xn, c(prod(p), n))
  Yhat <- t(tp1) %*% tp2
  Epsilon <- Yn - Yhat
  mse <- sum(Epsilon^2)/n
  mse <- sum(diag(Epsilon %*% t(Epsilon)))/n
  return(list(mse=mse, Yhat=Yhat))
}
