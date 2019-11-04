#' Tensor Regression with Envelope Structure and Three Generic Envelope Estimation Approaches
#'
#' Provides three estimators for tensor response regression (TRR) and tensor predictor regression (TPR) models with tensor envelope structure. The three types of estimation approaches are generic and can be applied to any envelope estimation problems. The full Grassmannian (FG) optimization is often associated with likelihood-based estimation but requires heavy computation and good initialization; the one-directional optimization approaches (1D and ECD algorithms) are faster, stable and does not require carefully chosen initial values; the SIMPLS-type is motivated by the partial least squares regression and is computationally the least expensive.
#'
#' @author Wenjing Wang, Jing Zeng and Xin Zhang
#'
#' @references Cook RD, Zhang X (2016). “Algorithms for Envelope Estimation.” Journal of Computational and Graphical Statistics, 25(1), 284–300. doi:10.1080/10618600.2015.1029577.
#'
#'   Li L, Zhang X (2017). “Parsimonious Tensor Response Regression.” Journal of the American Statistical Association, 112(519), 1131–1146.
#'
#' Zhang X, Li L (2017). “Tensor Envelope Partial Least Squares Regression.” Technometrics, 59(4), 426–436.
#'
#'
#' Cook RD, Zhang X (2018). “Fast Envelope Algorithms.” Statistica Sinica, 28(3), 1179–1197.
#'
#' @docType package
#' @name TRES-package
#'
#' @examples
#' library(TRES)
#' ## Load data "bat"
#' data("bat")
#' Xn <- bat$Xn
#' Yn <- bat$Yn
#' fit <- TRR.fit(Xn, Yn, method="standard")
#'
#' ## Print cofficient
#' coef(fit)
#'
#' ## Print the summary
#' summary(fit)
#'
#' ## Make the prediction on the original dataset
#' predict(fit, Xn)
#'
#' ## Draw the plot of two-way coefficient tensor (or matrix)
#' plot(fit)
#'
#' @seealso Useful links:
#' \itemize{
#'  \item \url{https://github.com/jerryfsu3333/TRES}
#'  \item Report bugs at \url{https://github.com/jerryfsu3333/TRES/issues}
#' }
NULL
