#' Square simulated data
#'
#' Simulated data used in tensor predictor regression (TPR). The pattern of coefficient is a square.
#'
#' @docType data
#'
#' @usage data("square")
#'
#' @format A list consisting of three components:
#' \describe{
#'  \item{x}{A \eqn{32 \times 32 \times 200} tensor}
#'  \item{y}{A \eqn{1 \times 200} matrix}
#'  \item{coefficients}{A \eqn{32\times 32 \times 1} tensor with the square pattern}
#'  \item{Gamma}{Two envelope basis: \eqn{32 \times 2} matrices}
#' }
#'
#' @keywords datasets
#' @examples
#' data("square")
#' ## Coefficients
#' coeff <- square$coefficients
#' image(-coeff@data[, , 1], axes=TRUE, col = grey(seq(0, 1, length = 256)))
#' title('Coefficient matrix')
"square"
