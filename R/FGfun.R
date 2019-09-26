#' @export
FGfun <- function (W, M, U)
{
  Fw <- log(det(t(W) %*% M %*% W)) + log(det(t(W) %*% solve(M+U) %*% W))
  a <- M %*% W %*% solve(t(W) %*% M %*% W)
  b <- solve(M+U) %*% W %*% solve((t(W) %*% solve(M+U) %*% W))
  df <- 2*(a + b)
  return(list(F = Fw, G = df))
}