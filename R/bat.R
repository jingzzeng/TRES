#' Bat simulated data
#'
#' Simulated data under tensor response regression (TRR) model. Each response observation is a two-dimensional image, and each binary predictor observation takes values 0 and 1, representing two groups.
#'
#' @docType data
#'
#' @usage data("bat")
#'
#' @format A list consisting of two components:
#' \describe{
#'  \item{x}{A \eqn{1 \times 20} matrix, each entry takes values 0 and 1, representing two groups.}
#'  \item{y}{A \eqn{64\times 64\times 20} tensor, each matrix \code{y@data[,,i]} represents an image.}
#' }
#'
#' @details The dataset is generated from the tensor response regression (TRR) model:
#' \deqn{\mathbf{Y}_i = \mathbf{B}X_i + \boldsymbol{\epsilon}_i, \quad i = 1,\ldots, n,}
#' where \eqn{n=20} and the regression coefficient \eqn{\mathbf{B} \in R^{64\times 64}} is a given image with rank 14, representing the mean difference of the response \eqn{\mathbf{Y}} between two groups. To make the model conform to the envelope structure, we construct the envelope basis \eqn{\boldsymbol{\Gamma}_k} and the covariance matrices \eqn{\boldsymbol{\Sigma}_k, k=1,2}, of error term as following. With the singular value decomposition of \eqn{\mathbf{B}}, namely \eqn{\mathbf{B} = \boldsymbol{\Gamma}_1 \boldsymbol{\Lambda} \boldsymbol{\Gamma}_2^{\top}}, we choose the envelope basis as \eqn{\boldsymbol{\Gamma}_k \in R^{64\times 14}, k=1,2}. Then the envelope dimensions are \eqn{u_1 =  u_2 = 14}. We generate another two matrices \eqn{\boldsymbol{\Omega}_k \in R^{14\times 14} = \mathbf{A}_k \mathbf{A}_k^{\top}} and  \eqn{\boldsymbol{\Omega}_{0k} \in R^{50\times 50} =  \mathbf{A}_{0k}\mathbf{A}_{0k}^{\top}}, where \eqn{\mathbf{A}_k \in R^{14\times 14}} and \eqn{\mathbf{A}_{0k} \in R^{50\times 50}} are randomly generated from Uniform(0,1) elementwise. Then we set the covariance matrices \eqn{\boldsymbol{\Sigma}_k = \boldsymbol{\Gamma}_k\boldsymbol{\Omega}_k\boldsymbol{\Gamma}_k^{\top} + \boldsymbol{\Gamma}_{0k}\boldsymbol{\Omega}_{0k}\boldsymbol{\Gamma}_{0k}^{\top}}, followed by normalization with their Frobenius norms. We set the first 10 predictors \eqn{X_i, i=1,\ldots, 10,} as 1 and the rest as 0. The error term is then generated from two-way tensor (matrix) normal distribution \eqn{TN(\mathbf{0}; \boldsymbol{\Sigma}_1, \boldsymbol{\Sigma}_2)}.
#'
#' @keywords datasets
#' @examples
#' data("bat")
#' ## Coefficients
#' coeff <- bat$coefficients
#' image(-coeff@data[, , 1], axes=TRUE, col = grey(seq(0, 1, length = 256)))
#' title('Coefficient matrix')
"bat"

