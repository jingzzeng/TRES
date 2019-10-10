#' @export

vcov.Tenv <- function(object, ...){
  cl <- object$call
  if(cl[1] == "TPR()"){stop("there is no vcov() method for obejct from TPR().")}
  if(cl[1] == "TRR()"){
     tmp <- summary.Tenv(object)
     cat("\nThe standard error tensor:\n\n")
     tmp$se
  }
}
