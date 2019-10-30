#' @export
#' @import rTensor
kroncov <- function(Tn) {
  ss <- dim(Tn);
  n <- ss[length(ss)]
  r <- ss[1:(length(ss)-1)]
  m <- length(r)
  prodr <- prod(r)


  ## initialization ##
  lambda <- 1
  S <- NULL
  for (i in 1:m) {
    S[[i]] <- diag(r[i])
  }
  Sinvhalf <- S

  ## iteration ##
  if (m > 1) {
    for (isim in 1:5) {
       for (i in 1:m) {
         Si0 <- S[[i]]
         idx <- c(1:(m+1))[-i]
         len <- length(idx)
         Tsn <- rTensor::ttl(Tn, Sinvhalf[c(idx[1:(len-1)])], ms=idx[1:(len-1)])
         idxprod <- (r[i]/n)/prodr
         TsnTsn <- ttt(Tsn, Tsn, ms = idx)@data*idxprod
         S[[i]] <- TsnTsn/norm(TsnTsn, type = "F")
         Sinvhalf[[i]] <- sqrtm(S[[i]])$Binv
       }

       Tsn <- ttl(Tn, Sinvhalf, 1:m)

       lambda <- sum((Tsn@data)^2)/prod(c(r, n))
    }
  }else {
       lambda <- 1
       S[[m]] <- ttt(Tn, Tn, ms = 2)@data*(1/n)
  }
  return(list(lambda=lambda, S=S))
}
