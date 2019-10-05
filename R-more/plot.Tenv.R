plot.Tenv <- function(object, thrd = 0.05, main1=NULL, main2=NULL, ask = TRUE, ...){
  object <- summary(object)
  cf <- coef(object)
  data <- drop(cf@data)
  if(length(dim(data)) != 2){stop("Can only draw a 2-D plot. Check if drop(coef(object)@data) is a 2-D matrix.")}
  if(missing(main1)){main1 <- "Coefficient plot"}
  devAskNewPage(FALSE)
  image(-data, axes = TRUE, col = grey(seq(0, 1, length = 256)), main=main1, ...)
  if(!is.null(object$p_val)){
    p_val <- drop(object$p_val)
    if(length(dim(p_val)) != 2){stop("Can only draw a 2-D plot. Check if drop(coef(object)@data) is a 2-D matrix.")}
    if(missing(main2)){main2 <- "P_val plot"}
    devAskNewPage(TRUE)
    image((p_val > thrd), axes = TRUE, col = grey(seq(0, 1, length = 256)), main = main2, ...)
  }
}
