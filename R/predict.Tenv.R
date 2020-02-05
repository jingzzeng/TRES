#' Predict method for Tenv object.
#'
#' Predict response for object returned from \code{\link{TRR.fit}} and \code{\link{TPR.fit}} functions.
#'
#' @param object An object of class \code{"Tenv"}, as from \code{\link{TPR.fit}} or \code{\link{TRR.fit}}.
#' @param newdata The data to be used for prediction. It can be vector, matrix or tensor for the fit returned from\code{\link{TRR.fit}}, and can be matrix or tensor for the fit returned from \code{\link{TPR.fit}}.
#' @param ... Arguments passed to or from other methods.
#' @return
#' Return the predicted response.
#' @note If \code{newdata} is missing, the fitted response from \code{object} is returned.
#'
#' @examples
#' data("bat")
#' x <- bat$x
#' y <- bat$y
#' fit <- TRR.fit(x, y, method="standard")
#' predict(fit, x)

#' @export
#' @importFrom stats coef

predict.Tenv <- function(object, newdata, ...){
  if(missing(newdata)){
    warning("newdata is missing, the fitted response from the model is returned.")
    return(object$fitted.values)
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
    pred <- crossprod(tp1, tp2)
  }
  pred
}
