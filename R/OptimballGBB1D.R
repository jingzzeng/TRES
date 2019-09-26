
##################################################
#         1D optimization solve for gamma        #
##################################################
ballGBB1D <- function(M, U, opts=NULL) {
  W0 <- get_ini1D(M, U)
  if (is.null(opts$xtol))
    opts$xtol = 1e-8 else if (opts$xtol < 0 || opts$xtol > 1)
    opts$xtol = 1e-8


  if (is.null(opts$gtol))
    opts$gtol = 1e-8 else if (opts$gtol < 0 || opts$gtol > 1)
    opts$gtol = 1e-8

  if (is.null(opts$ftol))
    opts$ftol = 1e-12 else if (opts$ftol < 0 || opts$ftol > 1)
    opts$ftol = 1e-12

  if (is.null(opts$mxitr))
    opts$mxitr = 800

  X <- OptManiMulitBallGBB(W0, opts, fun1D, M, U)$X
  return(X)
}


##################################################
#  1D optimization solve for envelope basis      #
##################################################
#' @export
OptimballGBB1D <- function(M, U, u, opts=NULL) {

  # estimating M-envelope contains span(U)
  # where M>0 and is symmetric
  # dimension of the envelope is d
  # based on inv(M+U) and M

  if(dim(U)[1]!=dim(U)[2]){
    {U = U %*% t(U)}
  }

  if (is.null(opts$xtol))
    opts$xtol = 1e-8 else if (opts$xtol < 0 || opts$xtol > 1)
    opts$xtol = 1e-8


  if (is.null(opts$gtol))
    opts$gtol = 1e-8 else if (opts$gtol < 0 || opts$gtol > 1)
    opts$gtol = 1e-8

  if (is.null(opts$ftol))
    opts$ftol = 1e-12 else if (opts$ftol < 0 || opts$ftol > 1)
    opts$ftol = 1e-12

  if (is.null(opts$mxitr))
    opts$mxitr = 500

  p <-  dim(U)[2]

  if(u < p){
    Mnew <- M
    Unew <- U
    G <- matrix(0, p, u)
    G0 <- diag(p)
    for(k in 1:u){
      gk <- ballGBB1D(Mnew, Unew, opts)
      G[, k] <- G0 %*% gk
      G0 <- qr.Q(qr(G[, 1:k]),complete=T)[,(k+1):p]
      Mnew <- t(G0) %*% M %*% G0
      Unew <- t(G0) %*% U %*% G0
    }
    Ghat <- G
  } else { Ghat <- diag(p) }

  return(Ghat)
}
