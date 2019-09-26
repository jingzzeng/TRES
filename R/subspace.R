#' @export
subspace <- function(A,B){
  # calculate the distance between two subspaces span(A) and span(B)
  # A and B are both p-by-u matrix
  # returns a distance metric that is between 0 and 1
  Pa <- qr.Q(qr(A))
  Pa <- Pa %*% t(Pa)
  Pb <- qr.Q(qr(B))
  Pb <- Pb %*% t(Pb)
  u <- dim(A)[2]
  return(norm(Pa-Pb, type="F")/sqrt(2*u))
}
