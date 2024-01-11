get_sbm_waic <- function(A, z) {
  nn <- table(z)
  K <- length(nn)
  n <- length(z)
  
  out <- comp_blk_sums_and_sizes(A, z-1, K)
  M <- out$lambda
  N <- out$NN
  
  waic <- matrix(0, n, n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      
    if (A[i, j] == 1) {
      E_post <- A[i,j] * (log(M[z[i], z[j]] + 1) - log(N[z[i], z[j]] + 2))
      var_post <- trigamma(M[z[i], z[j]] + 1) - trigamma(N[z[i], z[j]] + 2)
    } else {
      E_post <- (1-A[i,j]) * (log(N[z[i], z[j]] - M[z[i], z[j]] + 1) - log(N[z[i], z[j]] + 2))
      var_post <- trigamma(N[z[i], z[j]] - M[z[i], z[j]] + 1) - trigamma(N[z[i], z[j]] + 2)
    }
      
      waic[i,j] <- E_post + var_post
    }
  }
  
  -sum(waic)
}