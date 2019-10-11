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
  f <- function(w) { W <- mw(w); log(t(W) %*% M %*% W) + log(t(W) %*% solve(M+U) %*% W)  }
  g <- function(w) { W <- mw(w); 2*(M %*% W/(as.numeric(t(W) %*% M %*% W))+
                                      solve(M+U) %*% W/(as.numeric(t(W) %*% solve(M+U) %*% W))) }

  prob <- new(mod$RProblem, f, g)
  mani.defn <- ManifoldOptim::get.stiefel.defn(n, 1)

  return(list(prob=prob, mani.defn=mani.defn))
}

##################################################
# Manifold1D optimization solve for gamma        #
##################################################

first1D <- function(M, U, params=NULL) {

  if(is.null(params$max_iter))
    params$max_iter=500 else if (params$max_iter < 0 || params$max_iter > 2^20)
    params$max_iter=500

  if(is.null(params$tol))
    params$tol=1e-08 else if (params$tol < 0 || params$tol > 1)
    params$tol=1e-08

  if(is.null(params$method))
    params$method="RBFGS"

  if(is.null(params$check))
    params$check= FALSE

  res  <- fun1D_mfd(M, U)
  prob <- res$prob
  mani.defn <- res$mani.defn
  mani.params <- get.manifold.params(IsCheckParams = params$check)
  solver.params <- get.solver.params(Max_Iteration = params$max_iter, Tolerance=params$tol, IsCheckParams = params$check)

  W0 <- get_ini1D(M, U)
  gamma <- ManifoldOptim::manifold.optim(prob, mani.defn, method = params$method,
                          mani.params = mani.params, solver.params = solver.params, x0 = W0)
  n <- dim(M)[2]
  return(matrix(gamma$xopt,n,1))
}


##################################################
#  1D optimization solve for envelope basis      #
##################################################
#' @export
manifold1D <- function(M, U, u, params=NULL){
  # estimating M-envelope contains span(U)
  # where M>0 and is symmetric
  # dimension of the envelope is d
  # based on inv(M+U) and M
  if(is.null(params$max_iter))
    params$max_iter=500 else if (params$max_iter < 0 || params$max_iter > 2^20)
      params$max_iter=500

  if(is.null(params$tol))
      params$tol=1e-08 else if (params$tol < 0 || params$tol > 1)
        params$tol=1e-08

  if(is.null(params$method))
        params$method="RBFGS"

  if(is.null(params$check))
        params$check= FALSE

  if(dim(U)[1]!=dim(U)[2]) {
    {U = U %*% t(U)}
  }

  p = dim(U)[2]
  if(u < p) {
    Mnew <- M
    Unew <- U
    G <- matrix(0, p, u)
    G0 <- diag(p)
    for(k in 1:u) {
      gk <- first1D(Mnew, Unew, params)
      G[, k] <- G0 %*% gk
      G0 <- qr.Q(qr(G[, 1:k]),complete=TRUE)[,(k+1):p]
      Mnew <- t(G0) %*% M %*% G0
      Unew <- t(G0) %*% U %*% G0
    }
    Ghat <- G
  } else { Ghat <- diag(p) }
  return(Ghat)
}

