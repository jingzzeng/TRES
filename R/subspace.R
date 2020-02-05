#' The distance between two subspaces.
#'
#' This function calculates the distance between the two subspaces with equal dimensions \eqn{\mathrm{span}(\mathbf{A})} and \eqn{\mathrm{span}(\mathbf{B})}, where \eqn{\mathbf{A}\in R^{p\times u}} and \eqn{\mathbf{B} \in R^{p\times u}} are the basis matrices of two subspaces. The distance is defined as
#' \deqn{\|\mathbf{P}_{\mathbf{A}} - \mathbf{P}_{\mathbf{B}}\|_F/(2d),}
#' where \eqn{\mathbf{P(\cdot)}} is the projection matrix onto the given subspace with the standard inner product, and \eqn{d} is the common dimension.
#'
#' @param  A A \eqn{p}-by-\eqn{u} matrix.
#' @param  B A \eqn{p}-by-\eqn{u} matrix.
#'
#' @return Returns a distance metric that is between 0 and 1
#'
#' @export

subspace <- function(A,B){
  Pa <- qr.Q(qr(A))
  Pa <- tcrossprod(Pa)
  Pb <- qr.Q(qr(B))
  Pb <- tcrossprod(Pb)
  u <- dim(A)[2]
  return(norm(Pa-Pb, type="F")/sqrt(2*u))
}
