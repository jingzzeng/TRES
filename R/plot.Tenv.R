#' Plot coefficients and p-value for Tenv object.
#'
#' Plot method for object returned from \code{TRR.fit} and \code{TPR.fit} functions.
#'
#' \code{coef(x)} must be a two-way tensor or a matrix. For the object return from \code{TPR.fit}, only the coefficients plot is displayed. For the object return from \code{TRR.fit}, both the coefficients plot and p-value plot are displayed.
#'
#'
#' @param x An object of class \code{"Tenv"}, as from \code{TPR.fit} or \code{TRR.fit}.
#' @param level Significant level of p-value. Default is 0.05.
#' @param main Title to coefficient plot.
#' @param main_p Title to \eqn{p}-value plot.
#' @param xticks Set tick labels of the x-axis
#' @param yticks Set tick labels of the y-axis
#' @param ... Other parameters to be passed through to plotting functions.
#'
#' @seealso \code{\link{TRR.fit}, \link{TPR.fit}}
#'
#' @examples
#'  data("bat")
#'  x <- bat$x
#'  y <- bat$y
#'  fit <- TRR.fit(x, y, method="standard")
#'  plot(fit)
#'
#'  ## Set xticks and yticks
#'  plot(fit, xticks = seq(0, 10, length.out=5), yticks = seq(0, 10, length.out=5))
#'
#'  ## Change the significant level to 0.1
#'  plot(fit, level = 0.1)

#' @export
#' @importFrom stats coef
#' @importFrom grDevices grey
#' @importFrom graphics image axis
plot.Tenv <- function(x, level = 0.05, main=paste0("Coefficient plot ", "(", x$method, ")"), main_p = paste0("P value plot ", "(", x$method, ")"), xticks, yticks, ...){
  cf <- coef(x)
  data <- drop(cf@data)
  if(length(dim(data)) != 2){stop("Can only draw a 2-D plot. Check if coefficients is a 2-D matrix.")}
  data <- t(apply(data, 2, rev))
  if(missing(xticks)){xticks <- seq(0,dim(data)[1], length.out = 5)}
  if(missing(yticks)){yticks <- seq(0,dim(data)[1], length.out = 5)}
  image(-data, axes = FALSE, col = grey(seq(0, 1, length = 256)), main=main, ...)
  axis(1, at = seq(0,1, along.with = xticks), labels = xticks, tick = TRUE, ...)
  axis(2, at = seq(0,1, along.with = yticks), labels = yticks, tick = TRUE, ...)
  axis(3, labels = FALSE, lwd.ticks = 0)
  axis(4, labels = FALSE, lwd.ticks = 0)
  if(x$call[1] == "TRR.fit()"){
    tmp <- Tenv_Pval(x$x, x$y, Bhat = coef(x))
    p_val <- drop(tmp$p_val@data)
    p_val <- t(apply(p_val, 2, rev))
    if(length(dim(p_val)) != 2){stop("Can only draw a 2-D plot. Check if p-value is a 2-D matrix.")}
    image((p_val > level), axes = FALSE, col = grey(seq(0, 1, length = 256)), main = main_p, ...)
    axis(1, at = seq(0,1, along.with = xticks), labels = xticks, tick = TRUE, ...)
    axis(2, at = seq(0,1, along.with = yticks), labels = yticks, tick = TRUE, ...)
    axis(3, labels = FALSE, lwd.ticks = 0)
    axis(4, labels = FALSE, lwd.ticks = 0)
  }
}
