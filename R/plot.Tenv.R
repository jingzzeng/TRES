#' @export

plot.Tenv <- function(object, thrd = 0.05, ask = TRUE, main=paste0("Coefficient plot ", "(", object$method, ")"),
                      main2 = paste0("P value plot ", "(", object$method, ")"), ...){
  object <- summary(object)
  cf <- coef(object)
  data <- drop(cf@data)
  if(length(dim(data)) != 2){stop("Can only draw a 2-D plot. Check if drop(coef(object)@data) is a 2-D matrix.")}
  devAskNewPage(FALSE)
  image(-data, axes = TRUE, col = grey(seq(0, 1, length = 256)), main=main, ...)
  if(!is.null(object$p_val)){
    p_val <- drop(object$p_val@data)
    if(length(dim(p_val)) != 2){stop("Can only draw a 2-D plot. Check if drop(coef(object)@data) is a 2-D matrix.")}
    devAskNewPage(ask)
    image((p_val > thrd), axes = TRUE, col = grey(seq(0, 1, length = 256)), main = main2, ...)
  }
}
