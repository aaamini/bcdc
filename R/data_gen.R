scale_sim_data <- function(n, p, r, K) {
  # n       : number of nodes
  # p       : prob of edge within cluster
  # r       : p*r = prob of edge between clusters
  # K       : number of communities
  
  require(igraph)
  require(abind)
  
  P <- matrix(p*r, K, K)
  diag(P) <- p
  
  # n_i <- rep(n/K, K)
  z_true = as.integer( gl(K, ceiling(n/K))[1:n] )
  A <- get.adjacency(sample_sbm(n, P, table(z_true)))
  
  Xc = matrix(rnorm(n, 2*(z_true %% 2) -1), ncol=1)
  Xd = ifelse(z_true <= K/2 , 0, 1)
  # feature <- cbind(c(rep(0:1, c(n/2, n/2)))
  #                  , rnorm(n, rep(2*(1:K)%%2-1, each = n/K)))
  # 
  return(list(A = A, Xc = Xc, Xd = Xd, z_true = z_true))
}