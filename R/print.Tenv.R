#' @method print Tenv
#' @importFrom stats coef

print.Tenv <- function(x, ...){
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Coefficients:\n")
  cf <- coef(x)
  if (!is.null(cf)) {
    print(cf)
  }
  else cat("No coefficients\n")
  cat("\n")
  invisible(x)
}