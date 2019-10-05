summary.Tenv <- function(object, ...){
  if (!inherits(object, "Tenv"))
    warning("calling summary.Tenv(<fake-Tenv-object>) ...")
  n <- dim(object$X)[length(dim(object$X))]
  if(object$call[1] == "TPR()"){
    object$mse <- sum(residuals(object)^2)/n
  }else if(object$call[1] == "TRR()"){
    object$mse <- sum(residuals(object)@data^2)/n
    object$p_val <- Tenv_Pval(object$Y, object$X, B_est = coefficients(object))$P_val
  }
  object$n <- n
  object$xdim <- dim(object$X)
  object$ydim <- dim(object$Y)
  class(object) <- "summary.Tenv"
  object
}
