##################################################
#         1D optimization function               #
##################################################
#' @export

fun1D <- function(W, M, U) {
  f <- log(t(W) %*% M %*% W) + log(t(W) %*% chol2inv(chol(M+U)) %*% W)
  df <- 2*(M %*% W/(as.numeric(t(W) %*% M %*% W))+
             solve(M+U) %*% W/(as.numeric(t(W) %*% chol2inv(chol(M+U)) %*% W)))

  return(list(F = f, G = df))
}
