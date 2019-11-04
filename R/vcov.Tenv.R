#' Calculate standard error for Tenv object.
#'
#' Returns the standard error for Tenv object returned from \code{\link{TRR.fit}}. There is no \code{vcov} method for object returned from \code{\link{TPR.fit}}.
#'
#' @param object an object of class "Tenv", as from \code{\link{TPR.fit}} or \code{\link{TRR.fit}}.
#' @param ... arguments to be passed to or from other methods.
#' @return The standard error tensor is returned.
#' @seealso \code{\link{summary.Tenv}}
#' @note Since there is no variance-covariance matrix for tensor coefficient estimator, generic function \code{vcov} for \code{Tenv} object returns the standard error elementwise which is different from what the function \code{lm()} returns.
#' @examples
#' data("bat")
#' Xn <- bat$Xn
#' Yn <- bat$Yn
#' fit <- TRR.fit(Xn, Yn, method="standard")
#' vcov(fit)
#' @export

vcov.Tenv <- function(object, ...){
  cl <- object$call
  if(cl[1] == "TPR.fit()"){stop("there is no vcov() method for obejct from TPR.fit().")}
  if(cl[1] == "TRR.fit()"){
     tmp <- summary.Tenv(object)
     cat("\nThe standard error tensor:\n\n")
     tmp$se
  }
}
