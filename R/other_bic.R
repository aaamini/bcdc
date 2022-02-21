count_edges <- function(A, xi, L) {
  # counts edges and non-edges between clusters
  # adopted from: https://github.com/danieledurante/TESTsbm
  #
  # A     : adjacency matrix A | z = k, xi ~ SBM(eta_k, xi)
  # xi    : cluster assignments
  # L     : number of clusters
  
  require(Matrix)
  
  n <- nrow(A)
  xi_mat <- Matrix(0, n, L, sparse = TRUE)
  
  for (i in seq_len(n)){
    xi_mat[i, xi[i]] <- 1
  }
  
  tmp   <- A %*% xi_mat
  edges <- crossprod(xi_mat, tmp) - diag(0.5*colSums(tmp*xi_mat), ncol(xi_mat))
  
  A_c <- 1 - A
  diag(A_c) <- 0
  tmp   <- A_c %*% xi_mat
  non_edges <- crossprod(xi_mat, tmp) - diag(0.5*colSums(tmp*xi_mat), ncol(xi_mat))
  
  return(m = list(edges = as.matrix(edges), non_edges = as.matrix(non_edges)))
}

sbm_bic <- function(A, Z) {
  n <- nrow(A)
  K <- length(unique(Z))
  
  m <- count_edges(A, Z, K)
  edges <- m$edges
  non_edges <- m$non_edges
  
  log.lik <- sum(lbeta(1 + gdata::lowerTriangle(edges, diag = TRUE)
                       , 1 + gdata::lowerTriangle(non_edges, diag = TRUE)))
  
  (K+1)*K/2 * log(n*(n-1)/2) - 2*log.lik
  
}