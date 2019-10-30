#' @export
OptManiMulitBallGBB <- function(X, opts=NULL, fun, ...) {
  # Line search algorithm for optimization on manifold:
  #
  #    min f(X), s.t., ||X_i||_2 = 1, where X \in R^{n,p}
  #        g(X) = grad f(X)
  #    X = [X_1, X_2, ..., X_p]
  #
  #
  #  each column of X lies on a unit sphere
  #  Input:
  #          X --- ||X_i||_2 = 1, each column of X lies on a unit sphere
  #      fun --- objective function and its gradient:
  #      [F, G] = fun(X,  data1, data2)
  #      F, G are the objective function value and gradient, repectively
  #      data1, data2 are addtional data, and can be more
  #      Calling syntax:
  #      [X, out]= OptManiMulitBallGBB(X0, @fun, opts, data1, data2);
  #
  #      opts --- option structure with fields:
  #             record = 0, no print out
  #             maxiter       max number of iterations
  #             xtol        stop control for ||X_k - X_{k-1}||
  #             gtol        stop control for the projected gradient
  #             ftol        stop control for |F_k - F_{k-1}|/(1+|F_{k-1}|)
  #                         usually, max{xtol, gtol} > ftol
  #
  #        Output:
  #             x --- solution
  #             g --- gradient of x
  #             Out --- output information
  #
  #
  #          For example, consider the maxcut SDP:
  #          X is n by n matrix
  #          max Tr(C*X), s.t., X_ii = 1, X psd
  #
  #          low rank model is:
  #           X = V'*V, V = [V_1, ..., V_n], V is a p by n matrix
  #           max Tr(C*V'*V), s.t., ||V_i|| = 1,
  #
  #           function [f, g] = maxcut_quad(V, C)
  #           g = 2*(V*C);
  #           f = sum(dot(g,V))/2;
  #          end
  #
  #         [x, g, out]= OptManiMulitBallGBB(x0, @maxcut_quad, opts, C);
  #
  #
  #
  #         Reference:
  #             Z. Wen and W. Yin
  #         A feasible method for optimization with orthogonality constraints
  #

  ##Size information
  X <- as.matrix(X)
  if (length(X) == 0)
    print("input X is an empty") else {
      n <- dim(X)[1]; k <- dim(X)[2]}

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


  # parameters for the  nonmontone line search by Raydan
  if (is.null(opts$m))
    opts$m = 10

  if (is.null(opts$STPEPS))
    opts$STPEPS = 1e-10

  if (is.null(opts$maxiter))
    opts$maxiter = 500 else if (opts$maxiter < 0 || opts$maxiter > 2^20)
    opts$maxiter = 500

  if (is.null(opts$nt))
    opts$nt = 5  else if (opts$nt < 0 || opts$nt > 100)
    opts$nt = 5

  if (is.null(opts$record))
    opts$record = 0

  # copy parameters
    xtol <- opts$xtol
    gtol <- opts$gtol
    ftol <- opts$ftol
    rho  <- opts$rho
    m <- opts$m
    STPEPS <- opts$STPEPS
    eta <- opts$eta
    gamma <- opts$gamma
    record <- opts$record

    nt <- opts$nt

    crit <- matrix(1, opts$maxiter, 3)

    # normalize x so that ||x||_2 = 1
    nrmX <- colSums(X*X)
    nrmX <- matrix(nrmX, 1, k)
    if (norm((nrmX-1), type = "F") > 1e-8) {
      X <- sweep(X, 2, sqrt(nrmX),"/")
     }
    args = list(X, ...)
    eva <- do.call(fun, args)
    f <- eva$F; g <- as.matrix(eva$G)
    out <- c()
    out$nfe <- 1


    Xtg <- colSums(X*g)
    Xtg <- matrix(Xtg, 1, k)
    gg <- colSums(g*g)
    gg <- matrix(gg, 1, k)
    XX <- colSums(X*X)
    XX <- matrix(XX, 1, k)
    XXgg <- XX*gg
    temp <- sweep(X, 2, Xtg, "*")
    dtX <- matrix(temp, n, k) - g
    nrmG <- norm(dtX, type = "F")

    Q <- 1; Cval <- f; tau <- opts$tau

    ## print iteration header if debug == 1
    if (record >= 1) {
        cat(paste('------ Gradient Method with Line search ----- ',"\n"),
            sprintf("%4s %8s %8s %10s %10s", 'Iter', 'tau', 'F(X)', 'nrmG', 'XDiff'))

        }
    if (record == 10) out$fvec = f

    ##main iteration
    for (itr in 1 : opts$maxiter) {
        Xp <- X; fp <- f; gp <- g; dtXP <- dtX

        nls <- 1; deriv = rho*nrmG^2

        while (1 > 0) {
           ## calculate g, f
           tau2 <- tau/2
           beta <- (1 + (tau2^2)*(-(Xtg^2) + XXgg))
           a1 <- ((1 + tau2*Xtg)^2 - (tau2^2)*XXgg)/beta
           a2 <- -tau*XX/beta
           X <- sweep(Xp, 2, a1, "*") + sweep(gp, 2, a2, "*")


           args = list(X, ...)
           eva <- do.call(fun, args)
           f <- eva$F; g <- as.matrix(eva$G)
           out$nfe <- out$nfe + 1

           if (f <= Cval - tau*deriv || nls >= 5)
                  break
               tau <- eta * tau
               nls <- nls + 1
        }

        if (record == 10)
          out$fvec <- rbind(out$fvec,f)

        #Xtg <- sapply(1:k, function(d) sum(X[,d]*g[,d]))
          Xtg <- colSums(X*g)
          Xtg <- matrix(Xtg, 1, k)
          gg <- colSums(g*g)
          gg <- matrix(gg, 1, k)
          XX <- colSums(X*X)
          XX <- matrix(XX, 1, k)
          XXgg <- XX*gg
          temp <- sweep(X, 2, Xtg, "*")
          dtX <- matrix(temp, n, k) - g

          nrmG <- norm(dtX, type = "F")
          s <- X - Xp
          XDiff <- norm(s, type = "F")/sqrt(n)
          FDiff <- abs(fp - f)/(abs(fp) + 1)

          if (record >= 1)
              cat(paste('------ Gradient Method with Line search ----- ',"\n"),
                            sprintf("%4s %8s %8s %10s %10s %10s", 'Iter', 'tau', 'F(X)', 'nrmG', 'XDiff','nls'))


          crit[itr,] <- cbind(nrmG, XDiff, FDiff)
          r <- length((itr - min(nt, itr)+1):itr)
          temp1 <- matrix(crit[(itr - min(nt, itr)+1):itr,], r, 3)
          mcrit <- colMeans(temp1)

          if ((XDiff < xtol && FDiff < ftol) || nrmG < gtol
                || all(mcrit[2:3] < 10 * c(xtol, ftol))) {
                   out$msg = "converge"
                   break
              }

           y <- dtX - dtXP
           sy <- sum(s*y)
           tau <- opts$tau
           sy <- abs(sy)

           if (sy > 0) {
              if (itr %% 2 == 0)
                  tau <- sum(s*s)/sy else { tau <- sy/sum(y*y)}

                  tau <- max(min(tau, 1e+20), 1e-20)
            }
            Qp <- Q; Q <- gamma*Qp + 1
            Cval <- (gamma*Qp*Cval + f)/Q
    }

            if (itr >= opts$maxiter)
            out$msg <- "exceed max iteration"


            Xn <- colSums(X*X)
            Xn <- matrix(Xn, 1, k)
            out$feasi <- norm((Xn - 1),type = "2")

            if (out$feasi > 1e-14) {
                nrmX <- colSums(X*X)
                X <- sweep(X, 2, sqrt(nrmX),"/")
                args = list(X, ...)
                eva <- do.call(fun, args)
                f <- eva$F; g <- as.matrix(eva$G)
                out$nfe <- out$nfe + 1
                nrmX.n <- colSums(X*X)
                out$feasi <- norm((nrmX.n - 1), type = "2" )
            }

            out$nrmG <- nrmG
            out$fval <- f
            out$itr <- itr

            return(list(X = X, g = g, out = out))
}
