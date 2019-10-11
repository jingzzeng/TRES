#' Plot coefficients and p-value for Tenv object.
#'
#' plot method for class "Tenv"
#'
#' \code{coef(x)} must be a mode-2 tensor or a matrix. For the object return from \code{\link{TPR}},
#' only the coefficients plot is displayed. But for the object return
#' from \code{\link{TRR}}, both the coefficients plot and p-value plot is displayed.
#'
#' @param x Tenv object, typically result of \code{\link{TPR}} or \code{\link{TRR}}.
#' @param thrd significant level of p-value. Default is 0.05.
#' @param ask logical; if TRUE, the user is asked before the p-value plot. (only used for the object from \code{\link{TRR}}).
#' @param main title to each plot.
#' @param ... other parameters to be passed through to plotting functions.
#' @examples
#'  data("bat")
#'  Xn <- bat$Xn
#'  Yn <- bat$Yn
#'  fit <- TRR(Xn, Yn, method="standard")
#'  plot(fit)

#' @export
#' @importFrom stats coef
#' @importFrom grDevices devAskNewPage grey
#' @importFrom graphics image
plot.Tenv <- function(x, thrd = 0.05, ask = TRUE, main=c(paste0("Coefficient plot ", "(", x$method, ")"),
                                                         paste0("P value plot ", "(", x$method, ")")), ...){
  x <- summary(x)
  cf <- coef(x)
  data <- drop(cf@data)
  if(length(dim(data)) != 2){stop("Can only draw a 2-D plot. Check if coefficients is a 2-D matrix.")}
  devAskNewPage(FALSE)
  image(-data, axes = TRUE, col = grey(seq(0, 1, length = 256)), main=main[1], ...)
  if(!is.null(x$p_val)){
    p_val <- drop(x$p_val@data)
    if(length(dim(p_val)) != 2){stop("Can only draw a 2-D plot. Check if p-value is a 2-D matrix.")}
    devAskNewPage(ask)
    image((p_val > thrd), axes = TRUE, col = grey(seq(0, 1, length = 256)), main = main[2], ...)
  }
}