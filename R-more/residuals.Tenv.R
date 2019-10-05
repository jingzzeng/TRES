residuals.Tenv <- function(object, ...){
  f <- fitted(object)
  Y <- object$Y
  r <- Y-f
  r
}
