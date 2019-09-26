##################################################
#    get initial value for 1D algorithm          #
##################################################
#' @export
get_ini1D <- function(M, U) {
  
  p <- dim(U)[2]
  v1 <- eigen(M)$vectors
  v2 <- eigen(M+U)$vectors
  v <- cbind(v1, v2)
  W0 <- v[, 1]
  Fw0 <- fun1D(W0, M, U)$F
  for (i in 2:(2*p)) {
    W <- v[, i]
    Fw <- fun1D(W, M, U)$F
    if (Fw < Fw0) {
      W0 <- W
      Fw0 <- Fw
    }
  }
  return(W0)
}
