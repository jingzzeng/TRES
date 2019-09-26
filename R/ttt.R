#' @import rTensor
#' @export


ttt <- function(X, Y, dims) {
   s1 <- dim(X)
   s2 <- dim(Y)
   idx_1 <- which(!(1:length(s1) %in% dims))
   idx_2 <- which(!(1:length(s2) %in% dims))
   mat_1 <- rTensor::unfold(X, row_idx = idx_1, col_idx = dims)@data
   mat_2 <- rTensor::unfold(Y, row_idx = idx_2, col_idx = dims)@data
   mat <- mat_1 %*% t(mat_2)
   mat <- rTensor::as.tensor(mat)
   return(mat)
}
