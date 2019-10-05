predict.Tenv <- function(object, newdata, se.fit = FALSE,
                         interval = c("none", "confidence", "prediction"), level = 0.95, na.action = na.pass, ...){
  # type = c("response", "terms"), terms = NULL
  if(missing(newdata)){
    stop("Missing newdata.")
  }
  Bhat <- object$Bhat
  if(is.vector(newdata)){newdata <- t(as.matrix(newdata))}
  if(Bhat@modes[Bhat@num_modes] != dim(newdata)[1]){stop("The dimensions don't match.")}
  m <- Bhat@num_modes
  pred <- ttm(Bhat, t(newdata), m)
  pred
}
