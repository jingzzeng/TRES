#' Tensor predictor regression
#'
#' This function is used for estimation of tensor predictor regression. The available method including standard OLS type estimation, PLS type of estimation as well as envelope estimation with FG, 1D and ECD approaches.
#'
#' Please refer to \strong{Details} of \code{\link{TPR_sim}} about the description of the tensor predictor regression model.
#'
#' @aliases TPR
#' @usage TPR.fit(x, y, u, method=c('standard', 'FG', '1D', 'ECD', 'PLS'), Gamma_init)
#'
#' @param x The predictor tensor instance of dimension \eqn{p_1\times p_2\times\cdots\times p_m \times n}, where \eqn{n} is the sample size. Array with the same dimensions and matrix with dimension \eqn{p\times n} are acceptable. If \code{y} is missing, \code{x} should be a list or an environment consisting of predictor and response datasets.
#' @param y The response matrix of dimension \eqn{r \times n}, where \eqn{n} is the sample size. Vector of length \eqn{n} is acceptable.
#' @param u The dimension of envelope subspace. \eqn{u=(u_1,\cdots, u_m)}. Used for methods "FG", "1D", "ECD" and "PLS". User can use \code{\link{TensPLS_cv2d3d}} to select dimension.
#' @param method The method used for estimation of tensor response regression. There are four possible choices.
#' \itemize{
#'   \item{\code{"standard"}}: The standard OLS type estimation.
#'   \item{\code{"FG"}}: Envelope estimation with full Grassmannian (FG) algorithm.
#'   \item{\code{"1D"}}: Envelope estimation with one dimensional optimization approaches by 1D algorithm.
#'   \item{\code{"ECD"}}: Envelope estimation with one dimensional optimization approaches by ECD algorithm.
#'   \item{\code{"PLS"}}: The SIMPLS-type estimation without manifold optimization.
#' }
#' @param Gamma_init A list specifying the initial envelope subspace basis for "FG" method. If missing, use the estimation from "1D" algorithm.
#'
#' @return
#'   \item{x}{The original predictor dataset.}
#'   \item{y}{The original response dataset.}
#'   \item{call}{The method call.}
#'   \item{method}{The method used.}
#'   \item{coefficients}{The estimation of regression coefficient tensor.}
#'   \item{Gamma}{The estimation of envelope subspace basis.}
#'   \item{Sigma}{A lists of estimated covariance matrices at each mode for the tensor predictors.}
#'   \item{fitted.values}{The fitted response matrix.}
#'   \item{residuals}{The residuals matrix.}
#'
#' @examples
#'
#' rm(list = ls())
#' # The dimension of predictor
#' p <- c(10, 10, 10)
#' # The envelope dimensions u.
#' u <- c(1, 1, 1)
#' # The dimension of response
#' r <- 5
#' # The sample size
#' n <- 200
#'
#' # Simulate the data with TPR_sim.
#' dat <- TPR_sim(p = p, r = r, u = u, n = n)
#' x <- dat$x
#' y <- dat$y
#' B <- dat$coefficients
#'
#' fit_std <- TPR.fit(x, y, method="standard")
#' fit_FG <- TPR.fit(x, y, u, method="FG")
#' fit_pls <- TPR.fit(x, y, u, method="PLS")
#'
#' rTensor::fnorm(B-stats::coef(fit_std))
#' rTensor::fnorm(B-stats::coef(fit_FG))
#' rTensor::fnorm(B-stats::coef(fit_pls))
#'
#' ## ----------- Pass a list or an environment to x also works ------------- ##
#' # Pass a list to x
#' l <- dat[c("x", "y")]
#' fit_std_l <- TPR.fit(l, method="standard")
#'
#' # Pass an environment to x
#' e <- new.env()
#' e$x <- dat$x
#' e$y <- dat$y
#' fit_std_e <- TPR.fit(e, method="standard")
#'
#' ## ----------- Use dataset "square" included in the package ------------- ##
#' # Note: it is time-consuming
#' \dontrun{
#'   data("square")
#'   x <- square$x
#'   y <- square$y
#'   fit_std <- TPR.fit(x, y, method="standard")
#' }
#'
#' @seealso \code{\link{summary.Tenv}} for summaries, calculating mean squared error from the prediction.
#'
#' \code{\link{plot.Tenv}}(via \code{graphics::image}) for drawing the two-dimensional coefficient plot.
#'
#' \code{\link{predict.Tenv}} for prediction.
#'
#' The generic functions \code{\link{coef}, \link{residuals}, \link{fitted}}.
#'
#' \code{\link{TensPLS_cv2d3d}} for selecting the dimension of envelope by cross-validation.
#'
#' \code{\link{TPR_sim}} for generating the simulated data used in tensor prediction regression.
#'
#' The simulated data \code{\link{square}} used in tensor predictor regression.
#'
#' @references Zhang, X., Li, L. (2017). Tensor Envelope Partial Least-Squares Regression. Technometrics, 59(4), 426-436.
#'
#' @export
#' @importFrom pracma sqrtm kron

# This function gives all the estimation of tensor predictor regression
# The tensor predictor should be 2-dimensional or 3-dimensional

TPR.fit <- function(x, y, u, method=c('standard', 'FG', '1D', 'ECD', 'PLS'), Gamma_init){
  cl <- match.call()
  method <- match.arg(method)
  if(missing(y)){
    tmp <- x
    if(is.list(tmp)){
      if(!is.null(names(tmp))){
        x <- tmp$x
        y <- tmp$y
      }
      else{
        if(length(x) < 2) stop("x or y is missing.")
        x <- tmp[[1]]
        y <- tmp[[2]]
      }
    }
    else if(is.environment(x)){
      x <- tmp$x
      y <- tmp$y
    }
    else{
      stop("y is null, x should be a list or an environment.")
    }
    if(is.null(x) || is.null(y)) stop("x or y is missing. Check names(x).")
  }
  if(!is.matrix(y)){
    if(is.vector(y)){
      y <- t(as.matrix(y))
    }
    else stop("y should be vector or matrix.")
  }
  if(!inherits(x, "Tensor")){
    if(is.matrix(x) || inherits(x, "array")){
      x <- as.tensor(x)
    }
    else stop("x should be matrix, array or Tensor.")
  }
  x_old <- x
  y_old <- y
  ss <- dim(x)
  len <- length(ss)
  n <- ss[len]
  if(n != dim(y)[2]){stop("Unmatched dimension.")}
  p <- ss[1:(len-1)]
  m <- length(p)
  r <- dim(y)[1]
  ##center the data
  muy <- as.matrix(apply(y, 1, mean))
  y <- y - muy[, rep(1, n)]
  mux <- apply(x@data, c(1:m), mean)
  ttmp <- lapply(1:n, function(x) mux)
  ttmp <- array(unlist(ttmp), c(p, n))
  ttmp2 <- x@data - ttmp

  x <- rTensor::as.tensor(ttmp2)
  vecx <- matrix(x@data, prod(p), n)
  res <- kroncov(x)
  lambda <- res$lambda
  Sigx <- res$S
  Sigx[[1]] <- lambda*Sigx[[1]]


  if(method == "standard") {
    if(length(dim(x))==4){
      tmp6 <- kron(chol2inv(chol(Sigx[[3]])), chol2inv(chol(Sigx[[2]])))
      Btil <- kron(tmp6, chol2inv(chol(Sigx[[1]]))) %*% tcrossprod(vecx, y)/n
    }else if(length(dim(x))==3) {
      tmp6 <- kron(chol2inv(chol(Sigx[[2]])), chol2inv(chol(Sigx[[1]])))
      Btil <- tmp6 %*% tcrossprod(vecx, y)/n
    }else if(length(dim(x))==2) {
      Btil <- chol2inv(chol(Sigx[[1]])) %*% tcrossprod(vecx, y)/n
    }
    Btil <- array(Btil, c(p, r))
    Bhat <- Btil
    Gamma1 <- NULL
  }else{
    if(missing(u)){stop("A user-defined u is required.")}
    if(method == "PLS") {
      res_PLS <- TensPLS_fit(x, y, Sigx, u)
      Gamma1 <- res_PLS$Gamma; PGamma <- res_PLS$PGamma

      if(length(dim(x))==4) {
        tmp7 <- kron(PGamma[[2]], PGamma[[1]])
        Bhat_pls <- kron(PGamma[[3]], tmp7) %*% tcrossprod(vecx, y)/n
      }else if(length(dim(x))==3) {
        tmp7 <- kron(PGamma[[2]], PGamma[[1]])
        Bhat_pls <- tmp7 %*% tcrossprod(vecx, y)/n
      }else if(length(dim(x))==2) {
        Bhat_pls <- PGamma[[1]] %*% tcrossprod(vecx, y)/n
      }
      Bhat <- array(Bhat_pls, c(p, r))
    }

    if(method == "1D") {
      Sinvhalf <- NULL
      for (i in 1:m) {
        Sinvhalf[[i]] <- sqrtm(Sigx[[i]])$Binv
      }
      SigY <- (n-1)*cov(t(y))/n
      Sinvhalf[[m+1]] <- sqrtm(SigY)$Binv

      C <- ttm(x, y, m+1)/n
      Gamma1 <- PGamma <- NULL
      for (i in 1:m) {
        idx <- c(1:(m+1))[-i]
        Ck <- ttl(C, Sinvhalf[idx], ms = idx)
        U <- unfold(Ck, row_idx = i, col_idx = idx)@data
        Uk <- tcrossprod(U)
        Gamma1[[i]] <- OptimballGBB1D(Sigx[[i]], Uk, u[i])
        tmp8 <- t(Gamma1[[i]]) %*% Sigx[[i]] %*% Gamma1[[i]]
        PGamma[[i]] <- Gamma1[[i]] %*% chol2inv(chol(tmp8)) %*% t(Gamma1[[i]]) %*% Sigx[[i]]
      }

      if(length(dim(x))==4) {
        tmp9 <- kron(PGamma[[2]], PGamma[[1]])
        Bhat_env <- kron(PGamma[[3]], tmp9) %*% tcrossprod(vecx, y)/n
      }else if(length(dim(x))==3) {
        tmp9 <- kron(PGamma[[2]], PGamma[[1]])
        Bhat_env <- tmp9 %*% tcrossprod(vecx, y)/n
      }else if(length(dim(x))==2) {
        Bhat_env <- PGamma[[1]] %*% tcrossprod(vecx, y)/n
      }
      Bhat <- array(Bhat_env, c(p, r))
    }

    if(method == "ECD") {
      Sinvhalf <- NULL
      for (i in 1:m) {
        Sinvhalf[[i]] <- sqrtm(Sigx[[i]])$Binv
      }
      SigY <- (n-1)*cov(t(y))/n
      Sinvhalf[[m+1]] <- sqrtm(SigY)$Binv

      C <- ttm(x, y, m+1)/n
      Gamma1 <- PGamma <- NULL
      for (i in 1:m) {
        idx <- c(1:(m+1))[-i]
        Ck <- ttl(C, Sinvhalf[idx], ms = idx)
        U <- unfold(Ck, row_idx = i, col_idx = idx)@data
        Uk <- tcrossprod(U)
        Gamma1[[i]] <- ECD(Sigx[[i]], Uk, u[i])
        tmp8 <- t(Gamma1[[i]]) %*% Sigx[[i]] %*% Gamma1[[i]]
        PGamma[[i]] <- Gamma1[[i]] %*% chol2inv(chol(tmp8)) %*% t(Gamma1[[i]]) %*% Sigx[[i]]
      }

      if(length(dim(x))==4) {
        tmp9 <- kron(PGamma[[2]], PGamma[[1]])
        Bhat_env <- kron(PGamma[[3]], tmp9) %*% tcrossprod(vecx, y)/n
      }else if(length(dim(x))==3) {
        tmp9 <- kron(PGamma[[2]], PGamma[[1]])
        Bhat_env <- tmp9 %*% tcrossprod(vecx, y)/n
      }else if(length(dim(x))==2) {
        Bhat_env <- PGamma[[1]] %*% tcrossprod(vecx, y)/n
      }
      Bhat <- array(Bhat_env, c(p, r))
    }

    if(method=='FG'){
      Sinvhalf <- NULL
      for (i in 1:m) {
        Sinvhalf[[i]] <- sqrtm(Sigx[[i]])$Binv
      }
      SigY <- (n-1)*cov(t(y))/n
      Sinvhalf[[m+1]] <- sqrtm(SigY)$Binv
      C <- ttm(x, y, m+1)/n
      Gamma1 <- PGamma <- NULL
      for (i in 1:m) {
        idx <- c(1:(m+1))[-i]
        Ck <- ttl(C, Sinvhalf[idx], ms = idx)
        U <- unfold(Ck, row_idx = i, col_idx = idx)@data
        idxprod <- (p[i]/r)/prod(p)
        Uk <- idxprod * tcrossprod(U)
        if(missing(Gamma_init)){
          init <-  OptimballGBB1D(Sigx[[i]], Uk, u[i], opts=NULL)
        }else{
          init <- Gamma_init[[i]]
        }
        Gamma1[[i]] <- OptStiefelGBB(init, opts=NULL, FGfun, Sigx[[i]], Uk)$Gamma
        tmp8 <- t(Gamma1[[i]]) %*% Sigx[[i]] %*% Gamma1[[i]]
        PGamma[[i]] <- Gamma1[[i]] %*% chol2inv(chol(tmp8)) %*% t(Gamma1[[i]]) %*% Sigx[[i]]
      }

      if(length(dim(x))==4) {
        tmp9 <- kron(PGamma[[2]], PGamma[[1]])
        Bhat_env <- kron(PGamma[[3]], tmp9) %*% tcrossprod(vecx, y)/n
      }else if(length(dim(x))==3) {
        tmp9 <- kron(PGamma[[2]], PGamma[[1]])
        Bhat_env <- tmp9 %*% tcrossprod(vecx, y)/n
      }else if(length(dim(x))==2) {
        Bhat_env <- PGamma[[1]] %*% tcrossprod(vecx, y)/n
      }
      Bhat <- array(Bhat_env, c(p, r))
    }
  }

  Bhat <- as.tensor(Bhat)
  tp1 <- matrix(Bhat@data, nrow = c(prod(p)))
  tp2 <- matrix(x_old@data, c(prod(p), n))
  fitted.values <- crossprod(tp1, tp2)
  residuals <- y_old - fitted.values
  output <- list(x = x_old, y = y_old, call = cl, method = method, coefficients=Bhat, Gamma=Gamma1, Sigma=Sigx, fitted.values = fitted.values, residuals=residuals)
  class(output) <- "Tenv"
  output
}
