#' Bat simulated data
#'
#' Simulated data used in Tensor Response Regression. The pattern of coefficient is a bat.
#'
#' @docType data
#'
#' @usage data(bat)
#'
#' @format A list consisting of three components:
#' \describe{
#'  \item{Xn}{A \eqn{1 \times 20} matrix}
#'  \item{Yn}{A \eqn{64\times 64\times 20} tensor}
#'  \item{coeffiicients}{A \eqn{64\times 64 \times 1} tensor with the bat pattern}
#' }
#'
#' @keywords datasets
#' @examples
#' data("bat")
#' ## Coefficients
#' coeff <- bat$coefficients
#' image(-coeff@data[, , 1], axes=TRUE, col = grey(seq(0, 1, length = 256)))
#' title('Coefficient matrix')
"bat"

