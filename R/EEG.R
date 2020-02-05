#' Electroencephalography (EEG) dataset for alcoholism study.
#'
#' EEG images data of subjects in alcoholic and control groups.
#'
#' @docType data
#'
#' @usage data("EEG")
#'
#' @format A list consisting of three components:
#' \describe{
#'  \item{x}{A \eqn{1 \times 84} matrix, denoting the subject in alcoholic and control groups as 1 and 0 respectively.}
#'  \item{y}{A \eqn{64 \times 64 \times 84} tensor, consisting 84 \emph{channels} by \emph{time} EEG images.}
#' }
#'
#' @references URL: \url{https://archive.ics.uci.edu/ml/datasets/EEG+Database}.
#'
#' Li L, Zhang X (2017). “Parsimonious Tensor Response Regression.” Journal of the American Statistical Association, 112(519), 1131–1146.
#'
#' @keywords datasets
#' @examples
#' data("EEG")
#' x <- EEG$x; y <- EEG$y
#' ## Estimate the envelope dimension, which turns out to be c(1,1).
#' # u <- TensEnv_dim(x, y)
#' u <- c(1,1)
#'
#' ## Fit the dataset with TRR.fit
#' fit_1D <- TRR.fit(x, y, u, method = "1D")
#'
#' ## The coefficient plot and p-value plot
#' plot(fit_1D, xlab = "Time", ylab = "Channels", yticks = seq(64,0, length.out=5))
#'
#' ## Uncomment to display the plots from different methods.
#' # fit_ols <- TRR.fit(x, y, method = "standard")
#' # fit_pls <- TRR.fit(x, y, u, method = "PLS")
#' # plot(fit_ols, xlab = "Time", ylab = "Channels", yticks = seq(64,0, length.out=5))
#' # plot(fit_pls, xlab = "Time", ylab = "Channels", yticks = seq(64,0, length.out=5))
#'
"EEG"
