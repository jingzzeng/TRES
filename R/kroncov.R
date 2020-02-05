#' The covariance estimation of a tensor random variable
#'
#' This function estimates the covariance of a tensor random variable. We assume the covariance of the tensor random variable has a separable Kronecker covariance structure, i.e. \eqn{\boldsymbol{\Sigma}=\boldsymbol{\Sigma}_{m}\otimes\cdots\otimes\boldsymbol{\Sigma}_{1}}. This algorithm is described in Manceur, A. M., & Dutilleul, P. (2013).
#'
#' The individual component covariance matrices \eqn{\boldsymbol{\Sigma}_i, i=1,\ldots, m} are not identifiable. So each matrix is normalized, and an overall normalizing constant \eqn{\lambda} is extract, which is defined as
#' \deqn{\boldsymbol{\Sigma}=\lambda \boldsymbol{\Sigma}_m \otimes \cdots \otimes \boldsymbol{\Sigma}_1.}
#'
#' @param Tn A \eqn{p_1\times\cdots p_m\times n} matrix, array or tensor, where \eqn{n} is the sample size.
#' @return
#' \item{lambda}{The normalizing constant.}
#' \item{S}{A matrix list, consisting of each normalized covariance matrix \eqn{\boldsymbol{\Sigma}_1,\ldots,\boldsymbol{\Sigma}_m}.}
#'
#' @references Manceur, A. M., & Dutilleul, P. (2013). Maximum likelihood estimation for the tensor normal distribution: Algorithm, minimum sample size, and empirical bias and dispersion. Journal of Computational and Applied Mathematics, 239, 37-49.
#'
#' @export
#' @import rTensor
#' @importFrom pracma sqrtm
kroncov <- function(Tn){
  ss <- dim(Tn)
  if(is.null(ss) || length(ss) <= 1) stop("The dimension of Tn should be larger than one.")
  if(!inherits(Tn, "Tensor")){
    if(is.matrix(Tn) || inherits(Tn, "array")){
      Tn <- as.tensor(Tn)
    }
    else stop("y should be matrix, array or Tensor.")
  }
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

  ## centering ##
  mu <- apply(Tn@data, seq_along(r), mean)
  tmp <- lapply(1:n, function(x) mu)
  tmp <- array(unlist(tmp), c(r, n))
  Tn <- as.tensor(Tn@data-tmp)
  ## iteration ##
  if (m > 1) {
    ## loop ##
    for (isim in 1:5) {
       for (i in 1:m) {
         Si0 <- S[[i]]
         idx <- c(1:(m+1))[-i]
         len <- length(idx)
         Tsn <- ttl(Tn, Sinvhalf[c(idx[1:(len-1)])], ms=idx[1:(len-1)])
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
