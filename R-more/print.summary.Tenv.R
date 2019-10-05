print.summary.Tenv <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n", sep = "")
  cat("\nDimensions:\n", "X:", dim(x$X), "\n", "Y:", dim(x$Y), "\n\n")
  cat("Sample size:", x$n, "\n\n")
  cat("Mean squared error:", x$mse, "\n\n")
  cat("Coefficients:\n")
  print.default(format(head(coefficients(x)@data), digits = digits), print.gap = 2L,
                quote = FALSE)
  cat("\n")
  if(!is.null(x$p_val)){
    cat("p-value:\n")
    print.default(format(head(x$p_val), digits = digits), print.gap = 2L,
                  quote = FALSE)
  }
}
