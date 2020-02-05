#' @export
FGfun <- function (W, M, U)
{
  Fw <- log(det(t(W) %*% M %*% W)) + log(det(t(W) %*% chol2inv(chol(M+U)) %*% W))
  a <- M %*% W %*% chol2inv(chol(t(W) %*% M %*% W))
  b <- chol2inv(chol(M+U)) %*% W %*% chol2inv(chol((t(W) %*% solve(M+U) %*% W)))
  df <- 2*(a + b)
  return(list(F = Fw, G = df))
}
