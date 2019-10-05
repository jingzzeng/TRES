# Section 4.4
# Dimension selection with different sample size
library(TRES)
set.seed(1)

p <- 50
u <- 5

# Generate M and U
tmp <- matrix(runif(p*u), p, u)
Gamma <- qr.Q(qr(tmp))
Gamma0 <- qr.Q(qr(Gamma),complete=T)[,(u+1):p]

A <- matrix(runif(u^2), u, u)
Omega <- A%*%t(A)

A <- matrix(runif((p-u)^2), p-u, p-u)
Omega0 <- A%*%t(A)

A <- matrix(runif(u^2), u, u)
Phi <- A%*%t(A)

# Model 1
U <- Gamma%*%Phi%*%t(Gamma)
M <- Gamma%*%Omega%*%t(Gamma) + Gamma0%*%Omega0%*%t(Gamma0)
M <- M + 0.00001*diag(1,p,p)

n0 <- c(50,70,100, 200, 400, 800)
u_est3 <- NULL
for (n in n0){
  # Generate Wishart sample Mhat and Uhat
  X <- mvrnorm(n, rep(0,p), M)
  Mhat <- (t(X)%*%X)/n
  X <- mvrnorm(n, rep(0,p), U)
  Uhat <- (t(X)%*%X)/n
  output <- ballGBB1D_bic(Mhat, Uhat, n, bic_max = p/2)
  u_est3 <- c(u_est3, output$u)
}

u_est3