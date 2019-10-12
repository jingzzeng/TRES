#' Calculate standard error for Tenv object.
#'
#' Returns the standard error for Tenv object returned from \code{\link{TRR}}. There is no \code{vcov} method for object returned
#' from \code{\link{TPR}}.
#'
#' @param object an object of class "Tenv", as from \code{\link{TPR}} or \code{\link{TRR}}.
#' @param ... arguments to be passed to or from other methods.
#' @return The standard error tensor is returned.
#' @seealso \code{\link{summary.Tenv}}
#' @examples
#' data("bat")
#' Xn <- bat$Xn
#' Yn <- bat$Yn
#' fit <- TRR(Xn, Yn, method="standard")
#' vcov(fit)
#' @export

vcov.Tenv <- function(object, ...){
  cl <- object$call
  if(cl[1] == "TPR()"){stop("there is no vcov() method for obejct from TPR().")}
  if(cl[1] == "TRR()"){
     tmp <- summary.Tenv(object)
     cat("\nThe standard error tensor:\n\n")
     tmp$se
  }
}
