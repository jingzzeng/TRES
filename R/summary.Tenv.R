#' Summarize method for Tenv object.
#' Summary method for object returned from \code{\link{TRR}} and \code{\link{TPR}} functions.
#'
#' The p-value and standard error for coefficients are only calculated for the object returned from \code{\link{TRR}}.
#'
#' print.summary.Tenv print.summary.lm gives a more readable format of sample size, dimensions of datasets, mse,
#' and additionally gives \code{p-val} and \code{se} if \code{object} is returned from \code{\link{TRR}}.
#'
#' @param object An object of class "Tenv", as from \code{\link{TPR}} or \code{\link{TRR}}.
#' @param ... Arguments to be passed to or from other methods.
#' @param x An object of class "summary.Tenv", usually, a result of a call to summary.Tenv.
#' @name summary.Tenv
#' @return \code{object} With additional components
#' \describe{
#'  \item{n}{Sample size}
#'  \item{xdim}{Dimensions of predictor}
#'  \item{ydim}{Dimensions of response}
#'  \item{mse}{Mean squared error based on \code{residuals(object)}}
#'  \item{p_val}{Only for object returned from \code{\link{TRR}}, p-value for coefficients}
#'  \item{se}{Only for object returned from \code{\link{TRR}}, standard error for coefficients}
#' }
#' @examples
#' data("bat")
#' Xn <- bat$Xn
#' Yn <- bat$Yn
#' fit <- TRR(Xn, Yn, method="standard")
#' ##print summary
#' summary(fit)
#' @seealso \code{\link{Tenv_Pval}}

#' @rdname summary.Tenv
#' @export
#' @importFrom stats coef residuals

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

#' @rdname summary.Tenv
#' @method print summary.Tenv
#' @importFrom stats coef
print.summary.Tenv <- function(x, ...){
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n", sep = "")
  cat("\nDimensions:\n", "X:", dim(x$Xn), "\n", "Y:", dim(x$Yn), "\n\n")
  cat("Sample size:", x$n, "\n\n")
  cat("Mean squared error:", x$mse, "\n\n")
  cat("Coefficients:\n")
  print(coef(x))
  cat("\n")
  if(!is.null(x$p_val)){
    cat("p-value:\n")
    print(x$p_val)
    cat("\n")
    cat("standard error:\n")
    print(x$se)
  }
  cat("\n")
  invisible(x)
}
