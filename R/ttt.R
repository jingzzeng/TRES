#' @import rTensor
#' @export


ttt <- function(X, Y, ms) {
   s1 <- dim(X)
   s2 <- dim(Y)
   idx_1 <- which(!(1:length(s1) %in% ms))
   idx_2 <- which(!(1:length(s2) %in% ms))
   mat_1 <- rTensor::unfold(X, row_idx = idx_1, col_idx = ms)@data
   mat_2 <- rTensor::unfold(Y, row_idx = idx_2, col_idx = ms)@data
   mat <- tcrossprod(mat_1, mat_2)
   mat <- rTensor::as.tensor(mat)
   return(mat)
}
