get_upper_tri = function(M, diag = F) {
  M[which(upper.tri(M, diag = diag))]
}

# log multivariate Beta function
lmbeta = function(a) {
  sum(lgamma(a)) - lgamma(sum(a))
}

get_sbm_bic = function(A, z, type = "exact", eps = 1e-20) {
  nn = table(z)
  K = length(nn)
  n = length(z)
  out = comp_blk_sums_and_sizes(A, z-1, K)
  M = out$lambda
  N = out$NN
  
  if (type == "exact") {
    # cat("exact calc \n")
    temp = lbeta(M + 1, N - M + 1)
    term1 = sum(get_upper_tri(temp, diag = T))
    term2 = lmbeta(nn+1)
    bic = -2*(term1 + term2)
    
  } else { # approximate BIC calculation
    # cat("approx. calc\n")
    etah = M / (N + eps)
    pih = nn / n
    
    temp = M * log(etah + eps) + (N-M) * log(1-etah + eps) 
    term1 = sum(get_upper_tri(temp, diag = T))
    term2 = sum(nn * log(pih + eps)) 
    cK = (K-1) + K*(K+1)/2 
    bic = -2*(term1 + term2) + c(K)*log(n*(n+1)/2)  
  }
  
  bic
}

get_bcdc_marg_loglike = function(A, z, eps = 1e-20) {
  nn = table(z)
  K = max(z) # length(nn)
  n = length(z)
  out = comp_blk_sums_and_sizes(A, z-1, K)
  M = out$lambda
  N = out$NN
  
  temp = lbeta(M + 1, N - M + 1)
  term1 = sum(get_upper_tri(temp, diag = T))
  # term2 = lmbeta(nn+1)
  term1 
}

