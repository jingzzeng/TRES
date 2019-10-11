#' @method print summary.Tenv
#' @importFrom stats coef

print.summary.Tenv <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
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
