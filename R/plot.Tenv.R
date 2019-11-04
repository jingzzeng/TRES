#' Plot coefficients and p-value for Tenv object.
#'
#' Plot method for object returned from \code{\link{TRR.fit}} and \code{\link{TPR.fit}} functions.
#'
#' \code{coef(x)} must be a two-way tensor or a matrix. For the object return from \code{\link{TPR.fit}}, only the coefficients plot is displayed. But for the object return from \code{\link{TRR.fit}}, both the coefficients plot and p-value plot is displayed.
#'
#' @param x An object of class "Tenv", as from \code{\link{TPR.fit}} or \code{\link{TRR.fit}}.
#' @param thrd Significant level of p-value. Default is 0.05.
#' @param main Title to each plot.
#' @param ... Other parameters to be passed through to plotting functions.
#' @examples
#'  data("bat")
#'  Xn <- bat$Xn
#'  Yn <- bat$Yn
#'  fit <- TRR.fit(Xn, Yn, method="standard")
#'  plot(fit)

#' @export
#' @importFrom stats coef
#' @importFrom grDevices grey
#' @importFrom graphics image
plot.Tenv <- function(x, thrd = 0.05, main=c(paste0("Coefficient plot ", "(", x$method, ")"),
                                                         paste0("P value plot ", "(", x$method, ")")), ...){
  x <- summary(x)
  cf <- coef(x)
  data <- drop(cf@data)
  if(length(dim(data)) != 2){stop("Can only draw a 2-D plot. Check if coefficients is a 2-D matrix.")}
  image(-data, axes = TRUE, col = grey(seq(0, 1, length = 256)), main=main[1], ...)
  if(!is.null(x$p_val)){
    p_val <- drop(x$p_val@data)
    if(length(dim(p_val)) != 2){stop("Can only draw a 2-D plot. Check if p-value is a 2-D matrix.")}
    image((p_val > thrd), axes = TRUE, col = grey(seq(0, 1, length = 256)), main = main[2], ...)
  }
}
