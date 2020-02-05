##################################################
# Manifold   function 1D                         #
##################################################
#' @import ManifoldOptim
#' @import methods
mod <- Module("ManifoldOptim_module", PACKAGE = "ManifoldOptim")
mani.params <- get.manifold.params(IsCheckParams = TRUE)

fun1D_mfd <- function(M, U) {
  n <- dim(M)[2]
  mw <- function(w) { matrix(w, n, 1) }
  f <- function(w) { W <- mw(w); log(t(W) %*% M %*% W) + log(t(W) %*% chol2inv(chol(M+U)) %*% W)  }
  g <- function(w) { W <- mw(w); 2*(M %*% W/(as.numeric(t(W) %*% M %*% W))+
                                      chol2inv(chol(M+U)) %*% W/(as.numeric(t(W) %*% chol2inv(chol(M+U)) %*% W))) }

  prob <- new(mod$RProblem, f, g)
  mani.defn <- ManifoldOptim::get.stiefel.defn(n, 1)

  return(list(prob=prob, mani.defn=mani.defn))
}

##################################################
# Manifold1D optimization solve for gamma        #
##################################################

first1D <- function(M, U, opts=NULL) {

  if(is.null(opts$maxiter))
    opts$maxiter=500 else if (opts$maxiter < 0 || opts$maxiter > 2^20)
    opts$maxiter=500

  if(is.null(opts$tol))
    opts$tol=1e-08 else if (opts$tol < 0 || opts$tol > 1)
    opts$tol=1e-08

  if(is.null(opts$method))
    opts$method="RBFGS"

  if(is.null(opts$check))
    opts$check= FALSE

  res  <- fun1D_mfd(M, U)
  prob <- res$prob
  mani.defn <- res$mani.defn
  mani.params <- get.manifold.params(IsCheckParams = opts$check)
  solver.params <- get.solver.params(Max_Iteration = opts$maxiter, Tolerance=opts$tol, IsCheckParams = opts$check)

  W0 <- get_ini1D(M, U)
  gamma <- ManifoldOptim::manifold.optim(prob, mani.defn, method = opts$method,
                          mani.params = mani.params, solver.params = solver.params, x0 = W0)
  n <- dim(M)[2]
  return(matrix(gamma$xopt,n,1))
}


##################################################
#  1D optimization solve for envelope basis      #
##################################################
#' @export
manifold1D <- function(M, U, u, opts=NULL){
  # estimating M-envelope contains span(U)
  # where M>0 and is symmetric
  # dimension of the envelope is d
  # based on inv(M+U) and M
  if(is.null(opts$maxiter))
    opts$maxiter=500 else if (opts$maxiter < 0 || opts$maxiter > 2^20)
      opts$maxiter=500

  if(is.null(opts$tol))
      opts$tol=1e-08 else if (opts$tol < 0 || opts$tol > 1)
        opts$tol=1e-08

  if(is.null(opts$method))
        opts$method="RBFGS"

  if(is.null(opts$check))
        opts$check= FALSE

  if(dim(U)[1]!=dim(U)[2]) {
    U = tcrossprod(U)
  }

  p = dim(U)[2]
  if(u < p) {
    Mnew <- M
    Unew <- U
    G <- matrix(0, p, u)
    G0 <- diag(p)
    for(k in 1:u) {
      gk <- first1D(Mnew, Unew, opts)
      G[, k] <- G0 %*% gk
      G0 <- qr.Q(qr(G[, 1:k]),complete=TRUE)[,(k+1):p]
      Mnew <- t(G0) %*% M %*% G0
      Unew <- t(G0) %*% U %*% G0
    }
    Gamma <- G
  } else { Gamma <- diag(p) }
  return(Gamma)
}

