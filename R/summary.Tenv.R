#' Summarize method for Tenv object.
#'
#' Summary method for object returned from \code{TRR.fit} and \code{TPR.fit} functions.
#'
#' Extract \code{call}, \code{coefficients}, \code{residuals}, \code{Gamma} from \code{object}.
#'
#' The mean squared error \code{mse} is defined as \eqn{1/n\sum_{i=1}^n\|\mathbf{Y}_i-\hat{\mathbf{Y}}_i\|_F^2}, where \eqn{\hat{\mathbf{Y}}_i} is the prediction and \eqn{\|\cdot\|_F} is the Frobenius norm of tensor.
#'
#' For the object returned from \code{TRR.fit}, return the \eqn{p}-value and the standard error of estimated coefficient. However, since \eqn{p}-value and standard error depend on \eqn{\widehat{\mathrm{cov}}^{-1}\{\mathrm{vec}(\mathbf{X})\}} which is unavailable for the ultra-high dimensional \eqn{\mathrm{vec}(\mathbf{X})} in tensor predictor regression (TPR), the two statistics are not provided for the object returned from \code{TPR.fit}.
#'
#' print.summary.Tenv gives a more readable format of the statistics contained in \code{summary.Tenv}. If \code{object} is from \code{\link{TRR.fit}}, then \code{p-val} and \code{se} are also returned.
#'
#' @param object An object of class \code{"Tenv"}, as from \code{\link{TPR.fit}} or \code{\link{TRR.fit}}.
#' @param ... Arguments to be passed to or from other methods.
#' @param x An object of class "summary.Tenv", usually, a result of a call to summary.Tenv.
#' @name summary.Tenv
#' @return Return \code{object} with additional components
#'  \item{call}{The method call}
#'  \item{n}{Sample size}
#'  \item{xdim}{Dimensions of predictor}
#'  \item{ydim}{Dimensions of response}
#'  \item{coefficients}{The tensor coefficients estimated from \code{TPR.fit} or \code{TRR.fit}}
#'  \item{residuals}{The residuals, which equals to the response minus the fitted values}
#'  \item{Gamma}{A list of envelope subspace basis}
#'  \item{mse}{Mean squared error. The mean squared Frobenius norm of the difference between each response \eqn{\mathbf{Y}_i} and fitted value \eqn{\hat{\mathbf{Y}}_i},
#'    \deqn{1/n\sum_{i=1}^n\|\mathbf{Y}_i-\hat{\mathbf{Y}}_i\|_F^2}
#'  }
#'  \item{p_val}{Only for object returned from \code{TRR.fit}, p-value for coefficients}
#'  \item{se}{Only for object returned from \code{TRR.fit}, standard error for coefficients}
#'
#' @examples
#' data("bat")
#' x <- bat$x
#' y <- bat$y
#' fit <- TRR.fit(x, y, method="standard")
#' ##print summary
#' summary(fit)
#' @seealso \code{\link{Tenv_Pval}} is used to calculate the \eqn{p}-value and standard error.
#'
#' \code{\link{PMSE}} is used to calculate mean squared error for any provided datasets and coefficient.
#'
#' Fitting functions \code{\link{TRR.fit}}, \code{\link{TPR.fit}}.

#' @rdname summary.Tenv
#' @export
#' @importFrom stats coef residuals

summary.Tenv <- function(object, ...){
  ans <- object[c("call", "coefficients", "residuals", "Gamma")]
  if (!inherits(object, "Tenv"))
    warning("calling summary.Tenv(<fake-Tenv-object>) ...")
  n <- dim(object$x)[length(dim(object$x))]
  if(object$call[1] == "TPR.fit()"){
    ans$mse <- sum(residuals(object)^2)/n
  }else if(object$call[1] == "TRR.fit()"){
    ans$mse <- sum(residuals(object)@data^2)/n
    tmp <- Tenv_Pval(object$x, object$y, Bhat = coef(object))
    ans$p_val <- tmp$p_val
    ans$se <- tmp$se
  }
  ans$n <- n
  ans$xdim <- dim(object$x)
  ans$ydim <- dim(object$y)
  class(ans) <- "summary.Tenv"
  ans
}

#' @rdname summary.Tenv
#' @method print summary.Tenv
#' @export
#' @importFrom stats coef
print.summary.Tenv <- function(x, ...){
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n", sep = "")
  cat("\nDimensions:\n", "x:", x$xdim, "\n", "y:", x$ydim, "\n\n")
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
