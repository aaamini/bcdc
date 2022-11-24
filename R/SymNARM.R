multrnd_histc <- function(n, p) {
  require(pracma)
  
  p <- as.vector(p)
  
  edges <- c(0, cumsum(p))
  r <- pracma::histc(runif(n) * edges[length(edges)], edges)$cnt
  
  return(r[1:(length(r)-1)])
}

crt <- function(nobs, theta) {
  if (nobs==0) return(0)
  # Draws from the Chinese Restaurant Table (CRT) distribution 
  # N ~ sum_{t=1}^nobs Ber(theta / (theta + t-1))
  t = 1:nobs
  p = theta / (theta + t - 1 + 1e-5)
  
  sum(rbinom(nobs, 1, p))
}

crt_vec <- function(n_vec, theta_vec){
  if (length(n_vec) != length(theta_vec)) {
    stop("nvec and theta_vec should have the same length")
  }
  mapply(crt, n_vec, theta_vec)
  
  # sapply(seq_along(n_vec), function(i) crt(n_vec[i], theta_vec[i]))
  
  # idx = n_vec > 0
  # out = n_vec
  # out[idx] = mapply(crt, n_vec[idx], theta_vec[idx])
}

SymNARM <- function(B, X, K = 100, idx_train, idx_test
                    , Burnin = 1500, Collections = 1500, realmin = 1e-30) {
  # "Leveraging Node Attributes for Incomplete Relational Data"
  # https://github.com/ethanhezhao/NARM/blob/master/SymNARM/SymNARM.m
  
  # B           : upper triangular of adjacency matrix
  # X           : Covariates
  # K           : number of atoms
  # Burnin      : length of burn in
  # Collections : number of iterations after burn
  
  require(extraDistr)
  
  iterMax <- Burnin + Collections
  N <- nrow(B)
  
  BTrain_Mask <- matrix(0, N, N)
  BTrain_Mask[idx_train] <- 1
  BTrain_Mask <- BTrain_Mask + t(BTrain_Mask)
  
  BTrain <- B
  BTrain[idx_test] <- 0
  
  find <- which(BTrain > 0, arr.ind = TRUE)
  ii <- find[, 1]
  jj <- find[, 2]
  idx <- which(BTrain > 0)
  
  ProbSamples <- matrix(0, N, N)
  count <- 0
  EPS <- 0
  
  links <- which(B[idx_test] > 0)
    
  Epsilon <- 1
  beta1 <- 1
  beta2 <- 1
  c0 <- 1
  e_0 <- 1
  f_0 <- 1
  gamma0 <- 1
  
  Phi <- matrix(rgamma(N*K, 1, 1), N, K)
  r_k <- matrix(1/K, K, 1)
  Lambda_KK <- tcrossprod(r_k)
  diag(Lambda_KK) <- Epsilon * r_k
  c_i <- matrix(1, N, 1)
  
  X <- cbind(X, 1)
  
  L <- ncol(X)
  
  h <- matrix(1, L, K)
  
  g <- exp(X %*% log(h))
  
  mu_0 <- 1
  
  nu_0 <- 1 / mu_0
  
  active_nodes <- vector("list", length = L)
  
  for (l in 1:L) {
    active_nodes[[l]] <- which(X[, l] == 1, arr.ind = TRUE)
  }
  
  for (iter in 1:iterMax) {
    
    # EPM: Draw latent count for each edge
    Rate <- rowSums(Phi[ii, ] %*% Lambda_KK * Phi[jj, ])
    M <- extraDistr::rtpois(length(Rate), Rate)
    
    m_i_k_dot_dot <- matrix(0, K, N)
    m_dot_k_k_dot <- matrix(0, K, K)
    for (ij in 1:length(idx)) {
      pmf <- t(Phi[ii[ij], , drop = FALSE]) %*% Phi[jj[ij], ] * Lambda_KK
      mij_kk <- matrix(multrnd_histc(M[ij], pmf), K, K)
      m_i_k_dot_dot[, ii[ij]] <- m_i_k_dot_dot[, ii[ij]] + colSums(mij_kk)
      m_i_k_dot_dot[, jj[ij]] <- m_i_k_dot_dot[, jj[ij]] + t(rowSums(mij_kk))
      m_dot_k_k_dot <- m_dot_k_k_dot + mij_kk + t(mij_kk)
    } # end ij loop
    
    diag(m_dot_k_k_dot) <- diag(m_dot_k_k_dot) / 2
    
    Phi_times_Lambda_KK <- Phi %*% Lambda_KK
    
    # Node attributes
    p <- BTrain_Mask %*% Phi_times_Lambda_KK
    
    # Sample tables by CRP
    t <- matrix(0, N, K)
    t[t(m_i_k_dot_dot) > 0] <- 1
    for (i in 1:N) {
      for (k in 1:K) {
        if(m_i_k_dot_dot[k, i] < 2) next
        for (j in 1:(m_i_k_dot_dot[k, i] - 1)) {
          t[i, k] <- t[i, k] + runif(1) < g[i, k] / (g[i, k] + j)
        }
      }
    }
    
    # According to this: https://docs.tibco.com/pub/enterprise-runtime-for-R/5.0.0/doc/html/Language_Reference/base/sweep.html
    # sweep creates an array with the same dimensions as x by copying the data in STATS 
    # across the dimensions in x *not* specified by MARGIN.
    # log_p <- log(sweep(sweep(p, 2, c_i, "+"), 2, c_i, "/"))
    log_p <- log(sweep(sweep(p, 1, c_i, "+"), 1, c_i, "/"))
    new_h <- matrix(rgamma(L*K, scale = 1, shape = mu_0 + crossprod(X, t)), L, K)
    for (l in 1:L) {
      an_l <- active_nodes[[l]]
      new_h_l <- new_h[l, ] / (colSums(g[an_l, , drop = FALSE] * log_p[an_l, , drop = FALSE]) + nu_0 * h[l, ])
      
      g[an_l, ] <- g[an_l, ] * new_h_l
      
      h[l, ] <- new_h_l * h[l, ]
    }
    
    # EPM: Sample phi_ik
    # MATLAB: Y = randg(A) returns a matrix of random values chosen from gamma distributions with unit scale. Y is the same size as A, and randg generates each element of Y using a shape parameter equal to the corresponding element of A
    # MATLAB: Phi_temp = randg(g + m_i_k_dot_dot')
    Phi_temp <- matrix(rgamma(N*K, shape = g + t(m_i_k_dot_dot)), N, K)
    
    for (i in sample(N)) {
      Phi[i, ] <- Phi_temp[i, ] / (c_i[i] + BTrain_Mask[i, ] %*% Phi_times_Lambda_KK)
      Phi_times_Lambda_KK[i, ] <- Phi[i, ] %*% Lambda_KK
    }
    
   # EPM: Sample c_i
   # MATLAB: gamrnd(a,b) in MATLAB has shape parameter a and scale parameter b
   # MATLAB: c_i = gamrnd(1e-0 + sum(g,2), 1./(1e-0 +  sum(Phi,2)));  
   # c_i <- rgamma(1 + rowSums(g, 2), 1 / (1 +  rowSums(Phi)))
   c_i <- matrix(rgamma(N, shape = 1 + rowSums(g), rate = 1 +  rowSums(Phi)), ncol=1)
   
   
   Phi_KK <- crossprod(Phi, BTrain_Mask) %*% Phi
   
   diag(Phi_KK) <- diag(Phi_KK) / 2
   
   triu1dex <- triu(matrix(1, K, K), 1)
   diagdex <- Diagonal(K) # sparse(1:K,1:K,true);
   
   # EPM: Sample r_k
   L_KK = matrix(0, K, K)
   temp_p_tilde_k= matrix(0, K, 1)
   p_kk_prime_one_minus = matrix(0, K, K)
   for (k in sample(K)) {
    R_KK = t(r_k)
    R_KK[k] = Epsilon
    beta3 = beta2*matrix(1,1,K)
    beta3[k] = beta1
    p_kk_prime_one_minus[k,] = beta3 / (beta3 + Phi_KK[k,])
        
    # L_KK[k,] = CRT_sum_mex_matrix(sparse(m_dot_k_k_dot(k,:)),r_k(k)*R_KK);
    L_KK[k,] = crt_vec(m_dot_k_k_dot[k,], r_k[k]*R_KK)
    temp_p_tilde_k[k] = -sum(R_KK*log(max(p_kk_prime_one_minus[k,], realmin)))
    r_k[k] = rgamma(1, gamma0/K+sum(L_KK[k,])) / (c0+temp_p_tilde_k[k])
   }
  
   #  EPM: Sample gamma0 with independence chain M-H
   ell_tilde = crt_vec(as.numeric(rowSums(L_KK)), gamma0/K*rep(1,K))  # This is just ell_tilde_k ~ CRT(sum_{k_2} l_{k k_2}, gamma_0/K)
     
  } # end iter loop
  L_KK # for test only, we get the zero matrix
}