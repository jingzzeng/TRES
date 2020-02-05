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
  f <- function(w) { W <- mw(w); log(det(t(W) %*% M %*% W)) + log(det(t(W) %*% chol2inv(chol(M+U)) %*% W))  }
  g <- function(w) { W <- mw(w); 2*(M %*% W %*% chol2inv(chol(t(W) %*% M %*% W))+ chol2inv(chol(M+U)) %*% W %*% chol2inv(chol((t(W) %*% chol2inv(chol(M + U)) %*% W))) ) }

  prob <- new(mod$RProblem, f, g)
  mani.defn <- ManifoldOptim::get.stiefel.defn(n, u)

  return(list(prob=prob, mani.defn=mani.defn))
}


##################################################
#  FG optimization solve for envelope basis      #
##################################################
#' @export


manifoldFG <- function(M, U, u, Gamma_init, opts=NULL) {

  # estimating M-envelope contains span(U)
  # where M>0 and is symmetric
  # dimension of the envelope is d
  # based on inv(M+U) and M

   if(is.null(opts$maxiter))
    opts$maxiter=500 else if (opts$maxiter < 0 || opts$maxiter > 2^20)
    opts$maxiter=500

  if(is.null(opts$tol))
    opts$tol=1e-08 else if (opts$tol < 0 || opts$tol >+ 1)
    opts$tol=1e-08

  if(is.null(opts$method))
    opts$method="RBFGS"

  if(is.null(opts$check))
    opts$check= FALSE

  if(dim(U)[1] != dim(U)[2]) { U = tcrossprod(U)}
  n <- dim(M)[2]
  mani.params <- get.manifold.params(IsCheckParams = opts$check)
  solver.params <- ManifoldOptim::get.solver.params(Max_Iteration = opts$maxiter, Tolerance=opts$tol, IsCheckParams = opts$check)
  if (u < n) {
    res <- FGfun_mfd(M, U, u)
    prob <- res$prob
    mani.defn <- res$mani.defn

    gamma <- ManifoldOptim::manifold.optim(prob, mani.defn, method = opts$method, mani.params = mani.params,
                            solver.params = solver.params, x0 = Gamma_init)
    Gamma <- matrix(gamma$xopt,n,u)
  }
  else { Gamma <- diag(n)}
  return(Gamma)
}
