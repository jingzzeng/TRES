fitted.Tenv <- function(object, ...){
  cl <- object$call
  X <- object$X
  Bhat <- coefficients(object)
  if(cl[1] == "TRR()"){
    m <- Bhat@num_modes
    fitted <- ttm(Bhat, t(X), m)
  }else if(cl[1] == "TPR()"){
    ss <- dim(X)
    len <- length(ss)
    n <- ss[len]
    p <- ss[1:(len-1)]
    m <- length(p)
    tp1 <- matrix(Bhat@data, nrow = c(prod(p)))
    tp2 <- matrix(X@data, c(prod(p), n))
    fitted <- t(tp1) %*% tp2
  }
  fitted
}
