print.Tenv <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cf <- coefficients(x)
  if (length(cf)) {
    cat("Numeric Tensor of", object@num_modes, "Modes\n", sep = " ")
    cat("Modes: ", object@modes, "\n", sep = " ")
    cat("Coefficients:\n")
    print.default(format(cf@data, digits = digits), print.gap = 2L,
                  quote = FALSE)
  }
  else cat("No coefficients\n")
  cat("\n")
  invisible(x)
}
