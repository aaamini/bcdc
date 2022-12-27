cat1_sim_data <- function(n, p, r) {
  # n       : number of nodes
  # p       : prob of edge within cluster
  # r       : p*r = prob of edge between clusters
  
  require(igraph)
  require(abind)
  
  # sample SBM ----
  n1 <- n2 <- floor(n*1/3)
  n3 <- n - n1 - n2
  
  z_true <- rep(1:K, c(n1, n2, n3))
  
  P <- matrix(p*r, 3, 3)
  diag(P) <- p
  A <- get.adjacency(sample_sbm(n, P, c(n1, n2, n3)))
  
  # sample features ----
  feature <- cbind(c(rep(0:2, c(n1, n2, n3)))      # signal
                   , sample(0:2, n, replace = TRUE)) # noise
  
  return(list(A = A, Xd = feature, z_true = z_true))
}

cat2_sim_data <- function(n, p, r) {
  # n       : number of nodes
  # p       : prob of edge within cluster
  # r       : p*r = prob of edge between clusters
  
  require(igraph)
  require(abind)
  
  require(igraph)
  require(MCMCpack)
  require(abind)
  
  # sample SBM ----
  n1 <- floor(n*1/3)
  n2 <- n - n1
  
  z_true <- rep(1:K, c(n1, n2))
  
  P <- matrix(p*r, 2, 2)
  diag(P) <- p
  A <- get.adjacency(sample_sbm(n, P, c(n1, n2)))
  
  # sample features ----
  dir_p1 <- MCMCpack::rdirichlet(1, rep(1, 4))
  dir_p2 <- MCMCpack::rdirichlet(1, rep(1, 4))
  feature <- cbind(c(sample(0:3, n1, prob = dir_p1, replace = TRUE)
                     , sample(0:3, n2, prob = dir_p2, replace = TRUE))
                   , sample(0:3, n, prob = rep(1/4, 4), replace = TRUE))
  
  return(list(A = A, Xd = feature, z_true = z_true))
}

scale_sim_data <- function(n, p, r, K) {
  # n       : number of nodes
  # p       : prob of edge within cluster
  # r       : p*r = prob of edge between clusters
  # K       : number of communities
  
  require(igraph)
  require(abind)
  
  P <- matrix(p*r, K, K)
  diag(P) <- p
  
  z_true = as.integer( gl(K, ceiling(n/K))[1:n] )
  A <- get.adjacency(sample_sbm(n, P, table(z_true)))
  
  Xc = matrix(rnorm(n, 2*(z_true %% 2) -1), ncol=1)
  Xd = matrix(ifelse(z_true <= K/2 , 0, 1), ncol=1)
  
  return(list(A = A, Xc = Xc, Xd = Xd, z_true = z_true))
}

cont_sim_data <- function(n, p, r, mu) {
  # n       : number of nodes
  # p       : prob of edge within cluster
  # r       : p*r = prob of edge between clusters
  # mu      : continuous feature signal
  
  z_true <- sample(1:2, n, replace = T, prob = c(1, 2))
  
  eta <- matrix(c(p, r*p, r*p, p), nrow = 2)
  
  A <- nett::fast_sbm(z_true, eta)
  
  Mu <- rbind(c(mu, 0), c(-mu,0))
  X <- do.call(cbind, lapply(1:ncol(Mu), function(j) {
    rnorm(n, mean = Mu[z_true,j])
  }))
  
  return(list(A = A, Xc = X, z_true = z_true))
}

sparse_sim_data <- function(n) {
  # n       : number of nodes

  require(igraph)
  require(abind)
  
  require(igraph)
  require(abind)
  
  # sample SBM ----
  n1 <- floor(n*3/12)
  n2 <- floor(n*4/12)
  n3 <- n - n1 - n2
  
  z_true <- rep(1:3, c(n1, n2, n3))
  
  P <- cbind(c(1.6, 1.2, 0.16)
             , c(1.2, 1.6, 0.02)
             , c(0.16, 0.02, 1.2)) * 0.01
  A <- get.adjacency(sample_sbm(n, P, c(n1, n2, n3)))
  
  # sample features ----
  feature <- cbind(c(rnorm(n1, 0, 1), rnorm(n2, -1, 1), rnorm(n3, 1, 1))
                   , c(rnorm(n1, 2, 1), rnorm(n2, -0.8, 1), rnorm(n3, -0.8, 1))
                   , matrix(rnorm(n*98, 0, 1), n, 98)) # noise
  
  return(list(A = A, Xc = feature, z_true = z_true))
}

homophily_sim_data <- function(n, p, r, beta) {
  # n       : number of nodes
  # p       : prob of edge within cluster
  # r       : p*r = prob of edge between clusters
  # beta.   : homophily parameter
  
  n1 <- n2 <- floor(n*1/3)
  n3 <- n - n1 - n2
  
  community <- rep(0:2, c(n1, n2, n3))
  feature <- sample(0:1, n, replace = TRUE)
  
  A <- matrix(0, n, n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      p_ij <- ifelse(community[i] == community[j], p, p*r) +
        beta * (feature[i] == feature[j])
      
      A[i, j] <- sample(c(0, 1), 1, prob = c(1-p_ij, p_ij))
    }
  }
  
  A[lower.tri(A)] <- t(A)[lower.tri(A)]
  
  num <- data.frame(community, feature) %>%
    mutate(cluster = group_indices(., community, feature)) %>%
    select(cluster)
  
  return(list(A = as(A, "dgCMatrix"), Xd = cbind(feature), z_true = num$cluster))
  
}