#' @export

summary.Tenv <- function(object, ...){
  if (!inherits(object, "Tenv"))
    warning("calling summary.Tenv(<fake-Tenv-object>) ...")
  n <- dim(object$Xn)[length(dim(object$Xn))]
  if(object$call[1] == "TPR()"){
    object$mse <- sum(residuals(object)^2)/n
  }else if(object$call[1] == "TRR()"){
    object$mse <- sum(residuals(object)@data^2)/n
    tmp <- Tenv_Pval(object$Xn, object$Yn, B_est = coef(object))
    object$p_val <- tmp$P_val
    object$se <- tmp$se
  }
  object$n <- n
  object$xdim <- dim(object$Xn)
  object$ydim <- dim(object$Yn)
  class(object) <- "summary.Tenv"
  object
}
