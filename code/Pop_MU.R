# Section 4.3, Table 3.
# Compare the execution time and estimation error for each functions
library(TRES)

set.seed(1)
times <- 50
exe_time <- matrix(0,times,6, dimnames = list(NULL,  c('PLS', 'ECD', '1D_Mani', '1D_Feasi', 'FG_Mani', 'FG_Feasi')))
dist <- matrix(0,times,6, dimnames = list(NULL,  c('PLS', 'ECD', '1D_Mani', '1D_Feasi', 'FG_Mani', 'FG_Feasi')))

for(i in 1:times){
  cat('Time', i, '\n')
  p <- 20
  u <- 5
  
  # Construct U and M
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

  # # Model M2
  # U <- Gamma%*%Phi%*%t(Gamma)
  # M <- Gamma%*%t(Gamma) + 0.01*Gamma0%*%t(Gamma0)

  # # Model M3
  # U <- Gamma%*%Phi%*%t(Gamma)
  # M <- 0.01*Gamma%*%t(Gamma) + Gamma0%*%t(Gamma0)

  start_time <- Sys.time()
  Ghat_pls <- EnvMU(M, U, u)
  end_time <- Sys.time()
  exe_time[i, 1] <- difftime(end_time, start_time, units = 'secs')
  dist[i,1] <- subspace(Ghat_pls, Gamma)

  start_time <- Sys.time()
  Ghat_ecd <- ECD(M, U, u)
  end_time <- Sys.time()
  exe_time[i, 2] <- difftime(end_time, start_time, units = 'secs')
  dist[i,2] <- subspace(Ghat_ecd, Gamma)

  start_time <- Sys.time()
  Ghat_mani1D <- manifold1D(M, U, u)
  end_time <- Sys.time()
  exe_time[i, 3] <- difftime(end_time, start_time, units = 'secs')
  dist[i,3] <- subspace(Ghat_mani1D, Gamma)

  start_time <- Sys.time()
  Ghat_feasi1D <- OptimballGBB1D(M, U, u)
  end_time <- Sys.time()
  exe_time[i, 4] <- difftime(end_time, start_time, units = 'secs')
  dist[i,4] <- subspace(Ghat_feasi1D, Gamma)

  start_time <- Sys.time()
  Ghat_maniFG <- manifoldFG(M, U, u, Ghat_mani1D)
  end_time <- Sys.time()
  exe_time[i, 5] <- difftime(end_time, start_time, units = 'secs')
  dist[i,5] <- subspace(Ghat_maniFG, Gamma)

  start_time <- Sys.time()
  Ghat_feasiFG <- OptStiefelGBB(Ghat_feasi1D, opts=NULL, FGfun, M, U)$X
  end_time <- Sys.time()
  exe_time[i, 6] <- difftime(end_time, start_time, units = 'secs')
  dist[i,6] <- subspace(Ghat_feasiFG, Gamma)

}

# Average execution time and standard error for each method
mean_time <- apply(exe_time, 2, mean)
sd_time <- apply(exe_time, 2, sd)/sqrt(times)

# Average subspace distance and standard error for each method
mean_dist <- apply(dist, 2, mean)
sd_dist <- apply(dist, 2, sd)/sqrt(times)

mean_time
max(sd_time)
max(mean_dist)
max(sd_dist)