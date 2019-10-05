# Section 4.3
# Compare the performance of the FG algorithm with initial value as 1D estimator and 
# randomly generated matrix
library(TRES)

set.seed(1)

p <- 20
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


# Model M1
U <- Gamma%*%Phi%*%t(Gamma)
M <- Gamma%*%Omega%*%t(Gamma) + Gamma0%*%Omega0%*%t(Gamma0)
M <- M + 0.00001*diag(1,p,p)

# The 1D estimator
G3 <- manifold1D(M, U, u)
G4 <- OptimballGBB1D(M, U, u)

G5 <- manifoldFG(M, U, u, G3)
d1 <- subspace(G5, Gamma)
G6 <- OptStiefelGBB(G4, opts=NULL, FGfun, M, U)$X
d2 <- subspace(G6, Gamma)

# The randomly generated matrix
A_tmp <- matrix(runif(p*u), p, u)

G7 <- manifoldFG(M, U, u, A_tmp)
d3 <- subspace(G7, Gamma)
G8 <- OptStiefelGBB(A_tmp, opts=NULL, FGfun, M, U)$X
d4 <- subspace(G8, Gamma)

output3 <- c(d1,d2,d3,d4)
output3