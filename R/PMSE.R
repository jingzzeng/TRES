#' @export

PMSE <- function(Xn, Yn, Bhat){
  if(!is.matrix(Yn)){
    if(is.vector(Yn)){
      Yn <- t(as.matrix(Yn))
    }
    else stop("Yn should be vector or matrix.")
  }
  if(!inherits(Xn, "array") && !inherits(Xn, "matrix")){
    if(inherits(Xn, "Tensor")){
      Xn <- Xn@data
    }
    else stop("Xn should be matrix, array or Tensor.")
  }
  if(!inherits(Bhat, "array") && !inherits(Bhat, "matrix")){
    if(inherits(Bhat, "Tensor")){
      Bhat <- Bhat@data
    }else if(is.vector(Bhat)){
      Bhat <- t(as.matrix(Bhat))
    }
    else stop("Bhat should be vector, matrix, array or Tensor.")
  }
  ss <- dim(Xn)
  len <- length(ss)
  n <- ss[len]
  p <- ss[1:(len-1)]
  r <- dim(Yn)[1]
  m <- length(p)
  if(n!=dim(Yn)[2]){stop("Unmatched dimensions between Xn and Yn.")}
  if(dim(Bhat)[length(dim(Bhat))] != r){stop("Unmatched dimensions between Bhat and Yn.")}
  if(any(dim(Bhat)[1:(length(dim(Bhat))-1)] != p)){stop("Unmatched dimensions between Bhat and Xn.")}
  tp1 <- matrix(Bhat, c(prod(p), r))
  tp2 <- matrix(Xn, c(prod(p), n))
  Yhat <- t(tp1) %*% tp2
  Epsilon <- Yn - Yhat
  mse <- sum(Epsilon^2)/n
  return(list(mse=mse, Yhat=Yhat))
}
