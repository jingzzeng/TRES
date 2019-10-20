#' Predict method for TRR and TPR.
#'
#' Predict response for object returned from \code{\link{TRR}} and \code{\link{TPR}} functions.
#'
#' @param object An object of class "Tenv", as from \code{\link{TPR}} or \code{\link{TRR}}.
#' @param newdata The data to be used for prediction. It can be vector, matrix or tensor for \code{\link{TRR}} fits, and can be matrix
#' or tensor for \code{\link{TPR}} fits.
#' @param ... Arguments passed to or from other methods.
#' @return
#' \describe{
#'  \item{pred}{Predicted response.}
#' }
#' @examples
#' data("bat")
#' Xn <- bat$Xn
#' Yn <- bat$Yn
#' fit <- TRR(Xn, Yn, method="standard")
#' predict(fit, Xn)

#' @export
#' @importFrom stats coef

predict.Tenv <- function(object, newdata, ...){
  if(missing(newdata)){
    stop("Missing newdata.")
  }
  Bhat <- coef(object)
  if(object$call[1] == "TRR()"){
    if(!is.matrix(newdata)){
      if(is.vector(newdata)){
        newdata <- t(as.matrix(newdata))
      }
      if(inherits(newdata, "Tensor")){
        newdata <- newdata@data
      }
      else stop("newdata should be vector, matrix or tensor.")
    }
    if(is.vector(newdata)){newdata <- t(as.matrix(newdata))}
    m <- Bhat@num_modes
    pred <- rTensor::ttm(Bhat, t(newdata), m)
  }else if(object$call[1] == "TPR()"){
    if(!is.matrix(newdata)){
      if(inherits(newdata, "Tensor")){
        newdata <- newdata@data
      }
      else stop("newx should be Tensor or matrix.")
    }
    ss <- dim(newdata)
    len <- length(ss)
    n <- ss[len]
    p <- ss[1:(len-1)]
    tp1 <- matrix(Bhat, nrow = c(prod(p)))
    tp2 <- matrix(newdata, c(prod(p), n))
    pred <- t(tp1) %*% tp2
  }
  pred
}
