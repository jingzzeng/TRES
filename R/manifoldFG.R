##################################################
# Manifold   function full                       #
##################################################
#' @import ManifoldOptim
#' @import methods
mod <- Module("ManifoldOptim_module", PACKAGE = "ManifoldOptim")
mani.params <- get.manifold.params(IsCheckParams = TRUE)

FGfun_mfd <- function(M, U, u) {
  n <- dim(M)[2]
  mw <- function(w) { matrix(w, n, u) }
  f <- function(w) { W <- mw(w); log(det(t(W) %*% M %*% W)) + log(det(t(W) %*% solve(M+U) %*% W))  }
  g <- function(w) { W <- mw(w); 2*(M %*% W %*% solve(t(W) %*% M %*% W)+
                                      solve(M + U) %*% W %*% solve((t(W) %*% solve(M + U) %*% W))) }

  prob <- new(mod$RProblem, f, g)
  mani.defn <- ManifoldOptim::get.stiefel.defn(n, u)

  return(list(prob=prob, mani.defn=mani.defn))
}


##################################################
#  FG optimization solve for envelope basis      #
##################################################
#' @export


manifoldFG <- function(M, U, u, G_ini, params=NULL) {

  # estimating M-envelope contains span(U)
  # where M>0 and is symmetric
  # dimension of the envelope is d
  # based on inv(M+U) and M

   if(is.null(params$max_iter))
    params$max_iter=500 else if (params$max_iter < 0 || params$max_iter > 2^20)
    params$max_iter=500

  if(is.null(params$tol))
    params$tol=1e-08 else if (params$tol < 0 || params$tol >+ 1)
    params$tol=1e-08

  if(is.null(params$method))
    params$method="RBFGS"

  if(is.null(params$check))
    params$check= FALSE

  if(dim(U)[1] != dim(U)[2]) { U = U %*% t(U)}
  n <- dim(M)[2]
  mani.params <- get.manifold.params(IsCheckParams = params$check)
  solver.params <- ManifoldOptim::get.solver.params(Max_Iteration = params$max_iter, Tolerance=params$tol, IsCheckParams = params$check)
  if (u < n) {
    res <- FGfun_mfd(M, U, u)
    prob <- res$prob
    mani.defn <- res$mani.defn

    gamma <- ManifoldOptim::manifold.optim(prob, mani.defn, method = params$method, mani.params = mani.params,
                            solver.params = solver.params, x0 = G_ini)
    G_hat <- matrix(gamma$xopt,n,u)
  }
  else { G_hat <- diag(n)}
  return(G_hat)
}
