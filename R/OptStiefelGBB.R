#' @export
OptStiefelGBB <- function (X, opts=NULL, fun, ...)
{
  #
  # curvilinear search algorithm for optimization on Stiefel manifold
  #
  #   min F(X), S.t., X'*X = I_k, where X \in R^{n,k}
  #
  #   H = [G, X]*[X -G]'
  #   U = 0.5*tau*[G, X];    V = [X -G]
  #   X(tau) = X - 2*U * inv( I + V'*U ) * V'*X
  #
  #   -------------------------------------
  #   U = -[G,X];  V = [X -G];  VU = V'*U;
  #   X(tau) = X - tau*U * inv( I + 0.5*tau*VU ) * V'*X
  #
  #
  # Input:
  #           X --- n by k matrix such that X'*X = I
  #         fun --- objective function and its gradient:
  #                 fun(X,  data1, data2)
  #                 data1, data2 are addtional data, and can be more
  #               Calling syntax:
  #                 OptStiefelGBB(X0, fun, opts, data1, data2);
  #
  #        opts --- option structure with fields:
  #                 record = 0, no print out
  #                 mxitr       max number of iterations
  #                 xtol        stop control for ||X_k - X_{k-1}||
  #                 gtol        stop control for the projected gradient
  #                 ftol        stop control for |F_k - F_{k-1}|/(1+|F_{k-1}|)
  #                             usually, max{xtol, gtol} > ftol
  #
  # Output:
  #           X --- solution
  #         Out --- output information
  #
  # -------------------------------------
  # For example, consider the eigenvalue problem F(X) = -0.5*Tr(X'*A*X);
  #
  # function demo
  #
  # fun <- function(X,  A) {
  #   G = -(A %*% X)
  #   F = 0.5*sum(G*X)
  #   return(list(F = F, G = G))
  # }
  #
  #
  # n = 1000; k = 6;
  # A = matrix(rnorm(n^2), n, n); A = t(A) %*% A
  #
  # opts=c()
  # opts$record = 0;
  # opts$mxitr  = 1000;
  # opts$xtol = 1e-5;
  # opts$gtol = 1e-5;
  # opts$ftol = 1e-8;
  #
  # X0 = matrix(rnorm(n*k), n, k);    X0 = qr.Q(qr(X0));
  #
  # eva <- OptStiefelGBB(X0, opts, fun, A)
  # X <- eva$X
  # out <- eva$out
  # out$fval = -2*out$fval; # convert the function value to the sum of eigenvalues
  # sprintf('\nOptM: obj: %7.6f, itr: %f, nfe: %f, norm(XT*X-I): %3.2f \n',
  #            out$fval, out$itr, out$nfe, norm(t(X) %*% X - diag(k), type = "F") )
  #
  # -------------------------------------
  #
  # Reference:
  #  Z. Wen and W. Yin
  #  A feasible method for optimization with orthogonality constraints
  #
  # Author: Zaiwen Wen, Wotao Yin
  #   Version 1.0 .... 2010/10


  ##Size information
  X <- as.matrix(X)
  if (length(X) == 0)
    print("input X is an empty matrix") else {
    n <- dim(X)[1]; k <- dim(X)[2] }

  if (is.null(opts$xtol))
    opts$xtol = 1e-8 else if (opts$xtol < 0 || opts$xtol > 1)
    opts$xtol = 1e-8


  if (is.null(opts$gtol))
     opts$gtol = 1e-8 else if (opts$gtol < 0 || opts$gtol > 1)
     opts$gtol = 1e-8

  if (is.null(opts$ftol))
    opts$ftol = 1e-12 else if (opts$ftol < 0 || opts$ftol > 1)
    opts$ftol = 1e-12

  # parameters for control the linear approximation in line search
  if (is.null(opts$rho))
     opts$rho = 1e-04 else if (opts$rho < 0 || opts$rho > 1)
     opts$rho = 1e-04

  # factor for decreasing the step size in the backtracking line search
  if (is.null(opts$eta))
     opts$eta = 0.2 else if (opts$eta < 0 || opts$eta > 1)
     opts$eta = 0.1


  # parameters for updating C by HongChao, Zhang
  if (is.null(opts$gamma))
     opts$gamma = 0.85 else if (opts$gamma < 0 || opts$gamma > 1)
     opts$gamma = 0.85


  if (is.null(opts$tau))
      opts$tau = 1e-03  else if (opts$tau < 0 || opts$tau > 1)
      opts$tau = 1e-03


  #parameters for the  nonmontone line search by Raydan
  if (is.null(opts$STPEPS))
     opts$STPEPS = 1e-10


  if (is.null(opts$projG))
     opts$projG = 1  else if (opts$projG != 1 && opts$projG != 2)
     opts$projG = 1



  if (is.null(opts$iscomplex))
     opts$iscomplex = 0 else if (opts$iscomplex !=0 && opts$iscomplex != 1)
     opts$iscomplex = 0

  if (is.null(opts$mxitr))
     opts$mxitr = 500 else if (opts$mxitr < 0 || opts$mxitr > 2^20)
     opts$mxitr = 500

  if (is.null(opts$nt))
     opts$nt = 5 else if (opts$nt < 0 || opts$nt > 100)
     opts$nt = 5

  if (is.null(opts$record))
     opts$record = 0

  # copy parameters
  xtol <- opts$xtol
  gtol <- opts$gtol
  ftol <- opts$ftol
  rho  <- opts$rho
  STPEPS <- opts$STPEPS
  eta <- opts$eta
  gamma <- opts$gamma
  iscomplex <- opts$iscomplex
  record <- opts$record

  nt <- opts$nt
  crit <- matrix(1, opts$mxitr, 3)

  invH <- TRUE; if (k < n/2) { invH <- FALSE; eye2k <- diag(2*k) }

  ## Initial function value and gradient
  # prepare for iterations
  args = list(X, ...)
  eva <- do.call(fun, args)
  F <- eva$F; G <- eva$G
  out <- c()
  out$nfe <- 1
  #GX = t(G) %*% X
  GX <- crossprod(G, X)

  if (invH == TRUE) {
    #GXT <- G %*% t(X)
    GXT <- tcrossprod(G, X)
    H <- 0.5*(GXT - t(GXT))
    RX <- H %*% X
  } else {
    if (opts$projG == 1) {
      U <- cbind(G, X)
      V <- cbind(X, -G)
      #VU <- t(V) %*% U
      VU <- crossprod(V, U)
    }else if (opts$projG == 2){
      GB <- G - 0.5*X %*% (t(X) %*% G)
      U <- cbind(GB, X)
      V <- cbind(X, -GB)
      #VU <- t(V) %*% U
      VU <- crossprod(V, U)
    }
      # VX <- t(V) %*% X
      VX <- crossprod(V, X)
  }
  dtX <- G - X %*% GX
  nrmG <- norm(dtX, type = "F")

  Q <- 1; Cval <- F; tau <- opts$tau

  ## print iteration header if debug == 1
  if (opts$record == 1){
    cat(paste('------ Gradient Method with Line search ----- ',"\n"),
     sprintf("%4s %8s %8s %10s %10s", 'Iter', 'tau', 'F(X)', 'nrmG', 'XDiff')
    )
  }

  ## main iteration
  for (itr in 1:opts$mxitr) {
    XP <- X; FP <- F; GP <- G; dtXP <- dtX

    #scale step size
    nls <- 1; deriv <- rho*nrmG^2
    while (1 > 0) {
      if (invH == TRUE) {
        X <- solve(diag(n) + tau*H, XP - tau*RX) } else {
          aa <- solve(eye2k + (0.5*tau)*VU, VX)
          X <- XP - U %*% (tau*aa)
        }
      if (iscomplex == 0 && is.complex(X) == 1 )
        cat("error: X is complex")

      args <- list(X, ...)
      eva <- do.call(fun, args)
      F <- eva$F; G <- eva$G
      out$nfe <- out$nfe + 1

      if (F <= Cval - tau*deriv || nls >= 5)
        break;

      tau <- eta*tau; nls <- nls + 1
    }

    #GX <- t(G) %*% X
    GX <- crossprod(G, X)
    if (invH == TRUE) {
      #GXT <- G %*% t(X)
      GXT <- tcrossprod(G, X)
      H <- 0.5*(GXT - t(GXT))
      RX <- H %*% X
    } else {
      if (opts$projG == 1) {
        U = cbind(G, X); V = cbind(X, -G); #VU = t(V) %*% U;
        VU = crossprod(V, U)
      } else if (opts$projG == 2) {
        GB = G - 0.5 * X %*% (t(X) %*% G)
        U = cbind(GB, X); V = cbind(X, -GB); #VU = t(V) %*% U
        VU = crossprod(V, U)
      }
      # VX = t(V) %*% X
      VX = crossprod(V, X)
    }
    dtX = G - X %*% GX; nrmG = norm(dtX, type = "F")
    S = X - XP; XDiff = norm(S, type = "F")/sqrt(n)
    tau = opts$tau; FDiff = abs(FP - F)/(abs(FP) + 1)

    if (iscomplex == 1) {
      Y = dtX - dtXP; SY = abs(sum(Conj(S)*Y))
      if (itr %% 2 == 0) {
        tau = sum(Conj(S)*S)/SY } else {
          tau = SY/sum(Conj(Y)*Y)
        }
    } else {
      Y = dtX - dtXP; SY = abs(sum(S*Y))
      if (itr %% 2 == 0) {
        tau = sum(S*S)/SY } else {
          tau = SY/sum(Y*Y)
        }
    }

    if (is.na(tau))
      tau = opts$tau

    tau <- max(min(tau, 1e+20), 1e-20)

    if (record >= 1)
      sprintf("%4s  %3.2s  %4.3s  %3.2s  %3.2s  %3.2s  %2s",
              'iter', 'tau', 'F', 'nrmG', 'XDiff','FDiff', 'nls')

    crit[itr,] <- cbind(nrmG, XDiff, FDiff)
    r <- length((itr - min(nt, itr)+1):itr)
    tmp <- matrix(crit[(itr - min(nt, itr)+1):itr,], r, 3)
    mcrit <- colMeans(tmp)

    if ((XDiff < xtol && FDiff < ftol) || nrmG < gtol
        || all(mcrit[2:3] < 10 * c(xtol, ftol))) {
      if (itr <= 2) {
        ftol <- 0.1*ftol;
        xtol <- 0.1*xtol;
        gtol <- 0.1*gtol;
      } else {
        out$msg = "converge"
        break
      }
    }

    Qp <- Q; Q <- gamma*Qp + 1; Cval <- (gamma*Qp*Cval + F)/Q
  }

  if (itr >= opts$mxitr)
    out$msg = "exceed max iteration"

  #out$feasi <- norm((t(X) %*% X - diag(k)), type = "F")
  out$feasi <- norm((crossprod(X)- diag(k)), type = "F")
  if (out$feasi > 1e-13) {
    X <- qr.Q(qr(X))
    args <- list(X, ...)
    eva <- do.call(fun, args)
    F <- eva$F; G <- eva$G
    out$nfe <- out$nfe + 1
    #out$feasi <- norm((t(X) %*% X - diag(k)), type = "F")
    out$feasi <- norm((crossprod(X)- diag(k)), type = "F")
  }

  out$nrmG <- nrmG
  out$fval <- F
  out$itr <- itr

  return(list(X = X, out = out))
}
