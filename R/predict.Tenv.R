#' Predict method for Tenv object.
#'
#' Predict response for object returned from \code{\link{TRR.fit}} and \code{\link{TPR.fit}} functions.
#'
#' @param object An object of class "Tenv", as from \code{\link{TPR.fit}} or \code{\link{TRR.fit}}.
#' @param newdata The data to be used for prediction. It can be vector, matrix or tensor for the fit returned from\code{\link{TRR.fit}}, and can be matrix or tensor for the fit returned from \code{\link{TPR.fit}}.
#' @param ... Arguments passed to or from other methods.
#' @return
#' Return the predicted response.
#'
#' @examples
#' data("bat")
#' Xn <- bat$Xn
#' Yn <- bat$Yn
#' fit <- TRR.fit(Xn, Yn, method="standard")
#' predict(fit, Xn)

#' @export
#' @importFrom stats coef

predict.Tenv <- function(object, newdata, ...){
  if(missing(newdata)){
    stop("Missing newdata.")
  }
  Bhat <- coef(object)
  if(object$call[1] == "TRR.fit()"){
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
  }else if(object$call[1] == "TPR.fit()"){
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
