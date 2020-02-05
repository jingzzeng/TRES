#' @export

PMSE <- function(x, y, Bhat){
  if(!is.matrix(y)){
    if(is.vector(y)){
      y <- t(as.matrix(y))
    }
    else stop("y should be vector or matrix.")
  }
  if(!inherits(x, "array") && !inherits(x, "matrix")){
    if(inherits(x, "Tensor")){
      x <- x@data
    }
    else stop("x should be matrix, array or Tensor.")
  }
  if(!inherits(Bhat, "array") && !inherits(Bhat, "matrix")){
    if(inherits(Bhat, "Tensor")){
      Bhat <- Bhat@data
    }else if(is.vector(Bhat)){
      Bhat <- t(as.matrix(Bhat))
    }
    else stop("Bhat should be vector, matrix, array or Tensor.")
  }
  ss <- dim(x)
  len <- length(ss)
  n <- ss[len]
  p <- ss[1:(len-1)]
  r <- dim(y)[1]
  m <- length(p)
  if(n!=dim(y)[2]){stop("Unmatched dimensions between x and y.")}
  if(dim(Bhat)[length(dim(Bhat))] != r){stop("Unmatched dimensions between Bhat and y.")}
  if(any(dim(Bhat)[1:(length(dim(Bhat))-1)] != p)){stop("Unmatched dimensions between Bhat and x.")}
  tp1 <- matrix(Bhat, c(prod(p), r))
  tp2 <- matrix(x, c(prod(p), n))
  pred <- crossprod(tp1, tp2)
  Epsilon <- y - pred
  mse <- sum(Epsilon^2)/n
  return(list(mse=mse, pred=pred))
}
