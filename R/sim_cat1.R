# Libraries ----
library(pbmcapply)
library(NMI)
library(tidyverse)

# Functions ----
Rcpp::sourceCpp("./src/SBM.cpp", verbose = T)
Rcpp::sourceCpp("./src/CovarSBM.cpp", verbose = T)
source("./R/inference.R")
source("./R/CASC/cascTuningMethods.R")

simdata <- function(n, p, r) {
  # n       : number of nodes
  # p       : prob of edge within cluster
  # r       : p*r = prob of edge between clusters
  
  require(igraph)
  require(abind)
  
  # sample SBM ----
  n1 <- n2 <- floor(n*1/3)
  n3 <- n - n1 - n2
  
  P <- matrix(p*r, 3, 3)
  diag(P) <- p
  A <- get.adjacency(sample_sbm(n, P, c(n1, n2, n3)))
  
  # sample features ----
  feature <- cbind(c(rep(1:3, c(n1, n2, n3)))      # signal
                   , sample(3, n, replace = TRUE)) # noise
  
  return(list(A = A, Cov = feature, num = c(n1, n2, n3)))
}

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

bcdc <- function(A, Cov
                 , alpha = 10, beta = 1
                 , omega = 1
                 , burn, iter) {
  # A       : adjacency matrix
  # Cov     : covariates
  # alpha   : clusters ~ CRP(alpa)
  # beta    : pi ~ Beta(beta, beta)
  # omega
  # burn    : length of burn in
  # iter    : number of iterations
  
  require(abind)
  require(MCMCpack)
  
  fea <- dim(Cov)[2]
  
  # initialize ----
  
  # initialize clusters ~ CRP(alpha)
  cluster <- rep(NA, n)
  cluster[1] <- 1
  cluster_counts <- L <- 1
  for(i in 2:n) {
    p <- c(cluster_counts, alpha)
    cluster_star <- sample(seq_len(L+1), 1, prob = p)
    cluster[i] <- cluster_star
    if (cluster_star == (L+1)) {
      cluster_counts <- c(cluster_counts, 1)
      L <- L+1
    } else {
      cluster_counts[cluster_star] <- cluster_counts[cluster_star] + 1
    }
  }
  
  # initialize eta
  n_cat1 <- length(unique(Cov[, 1]))
  cd1 <- t(sapply(seq_len(L), function(i) {
    prop.table(tabulate(Cov[which(cluster == i), 1, drop = FALSE], n_cat1) + 1)
  }))
  
  n_cat2 <- length(unique(Cov[, 2]))
  cd2 <- t(sapply(seq_len(L), function(i) {
    prop.table(tabulate(Cov[which(cluster == i), 2, drop = FALSE], n_cat2) + 1)
  }))
  
  # initialize pi ~ Beta(beta, beta)
  pi <- matrix(nrow = L, ncol = L)
  pi[lower.tri(pi, diag = TRUE)] <- rbeta(L*(L+1)/2, beta, beta)
  pi[upper.tri(pi)] <- t(pi)[upper.tri(pi)]
  
  # MCMC
  sims <- matrix(NA, iter, n)
  
  for(k in seq_len(iter)) {
    
    # update clusters ----
    for(i in seq_len(n)) {
      
      pre <- cluster[i]
      cluster_counts[pre] <- cluster_counts[pre] - 1
      
      if (cluster_counts[pre] == 0) { # was only node in cluster
        removed_cluster <- TRUE
        
        L <- L-1
        cd1 <- cd1[-pre, , drop = FALSE]
        cd2 <- cd2[-pre, , drop = FALSE]
        pi <- pi[-pre, -pre, drop = FALSE]
        
        cluster <- ifelse(cluster > pre, cluster - 1, cluster)
        cluster_counts <- cluster_counts[-pre]
        
      } else {
        removed_cluster <- FALSE
      }
      
      # probability node is is assigned old cluster
      dn <- sapply(seq_len(L), function(j) log(cd1[j, Cov[i, 1]]) + log(cd2[j, Cov[i, 2]]))
      log_p <- sapply(seq_len(L), function(j) {
        A[i, -i] %*% log(pi[j, cluster[-i]]) +
          (1 - A[i, -i]) %*% log((1 - pi[j, cluster[-i]]))
      })
      p <- log(cluster_counts) + omega * dn + log_p
      
      # probability node is is assigned new cluster
      cd1_star <- MCMCpack::rdirichlet(1, rep(1, n_cat1))
      cd2_star <- MCMCpack::rdirichlet(1, rep(1, n_cat2))
      pi_star <- matrix(nrow = L+1, ncol = L+1)
      pi_star[seq_len(L), seq_len(L)] = pi
      pi_star[L+1, ] <- rbeta(L+1, beta, beta)
      pi_star[, L+1] <- pi_star[L+1, ]
      
      dn <- log(cd1_star[Cov[i, 1]]) + log(cd2_star[Cov[i, 2]])
      log_p <- A[i, -i] %*% log(pi_star[L+1, cluster[-i]]) +
        (1 - A[i, -i]) %*% log((1 - pi_star[L+1, cluster[-i]]))
      p <- c(p, log(alpha) + omega*dn + log_p)
      
      # sample cluster
      new <- sample(seq_len(L+1), 1, prob = exp(p - max(p)))
      cluster[i] <- new
      
      # if new cluster, update L, eta, and pi
      if (new == L+1) {
        L <- L+1
        cd1 <- rbind(cd1, cd1_star)
        cd2 <- rbind(cd2, cd2_star)
        pi <- pi_star
        
        cluster_counts <- c(cluster_counts, 1)
        
      } else {
        cluster_counts[new] <- cluster_counts[new] + 1
      }
      
    }
    
    # update eta ----
    cd1 <- t(sapply(seq_len(L), function(i) {
      prop.table(tabulate(Cov[which(cluster == i), 1, drop = FALSE], n_cat1) + 1)
    }))
    
    cd2 <- t(sapply(seq_len(L), function(i) {
      prop.table(tabulate(Cov[which(cluster == i), 2, drop = FALSE], n_cat2) + 1)
    }))
    
    # update pi ----
    if (L == 1) {
      A_n <- sum(A) / 2
      B_n <- n*(n - 1)/2 - A_n
      pi <- matrix(rbeta(1, beta + A_n, beta + B_n))
      
    } else {
      pi <- matrix(nrow = L, ncol = L)
      
      AB <- sapply(seq_len(L), function(l) {
        A_n <- sum(A[cluster == l, cluster == l]) / 2
        B_n <- (sum(cluster == l)^2 - sum(cluster == l)) / 2 - A_n
        return(c(A_n, B_n))
      })
      diag(pi) <- rbeta(L, beta + AB[1, ], beta + AB[2, ])
      
      if (L == 2) {
        A_n <- sum(A[cluster == 1, cluster == 2])
        B_n <- sum(cluster == 1) * sum(cluster == 2) - A_n
        pi[1, 2] <- rbeta(1, beta + A_n, beta + B_n)
        pi[2, 1] <- pi[1, 2]
        
      } else {
        AB <- do.call(cbind, sapply(seq_len(L-1), function(i) {
          sapply((i+1):L, function(j) {
            A_n <- sum(A[cluster == i, cluster == j])
            B_n <- sum(cluster == i) * sum(cluster == j) - A_n
            return(c(A_n, B_n))
          })
        }))
        
        pi[lower.tri(pi)] <- rbeta(L * ((L+1)/2 - 1)
                                   , beta + AB[1, ]
                                   , beta + AB[2, ])
        pi[upper.tri(pi)] <- t(pi)[upper.tri(pi)]  
      }
    } 
    
    # save clusters ----
    sims[k, ] <- cluster
  }
  
  mode_result <- apply(sims[(burn+1):iter, ], 2, getmode)
  
  return(mode_result)
  
}

# Simulation ----
K <- 3
n <- 150
p <- 0.1
r <- seq(0.1, 0.8, 0.05)
n_iter <- 1500
n_reps <- n_cores <- detectCores()

grid <- expand.grid(n, p, r)
names(grid) <- c("n", "p", "r")

results <- as.data.frame(matrix(nrow = nrow(grid) * n_reps, ncol = 8))
names(results) <- c("n", "p", "r"
                    , "bcdc", "casc", "kmeans", "SC", "bsbm")

for (tt in seq_len(nrow(grid))) {
  
  print(tt)
  
  n <- grid[tt, "n"]
  p <- grid[tt, "p"]
  r <- grid[tt, "r"]
  
  res <- simplify2array(
    pbmclapply(seq_len(n_reps), mc.cores = n_cores, function(i)
    {
      set.seed(i)
      
      sim <- simdata(n, p, r)
      A <- sim$A
      Cov <- sim$Cov
      simnum <- sim$num
      
      Z_true <- rep(1:K, simnum)
      
      # BCDC
      Z_bcdc <- bcdc(as.matrix(A), Cov, burn = .1*n_iter, iter = n_iter)
      
      # CASC
      Z_casc <- kmeans(getCascAutoSvd(A, Cov, 2, enhancedTuning = TRUE)$singVec
                       , centers = K, nstart = 20)$cluster
      
      # k-means
      Z_kmeans <- kmeans(Cov, centers = K, nstart = 20)$cluster
      
      # SC
      Z_SC <- nett::spec_clust(A, K)
      
      # BSBM
      bsbm <- new(SBM, A, K, alpha = 1, beta = 1)
      Z_bsbm <- get_map_labels(bsbm$run_gibbs(n_iter))$z_map
      
      return(c(nmi_wrapper(Z_true, Z_bcdc)
               , nmi_wrapper(Z_true, Z_casc)
               , nmi_wrapper(Z_true, Z_kmeans)
               , nmi_wrapper(Z_true, Z_SC)
               , nmi_wrapper(Z_true, Z_bsbm)))
      
    })
  )
  
  ind <- seq_len(n_reps) + n_reps * (tt - 1)
  
  results[ind, "n"] <- n
  results[ind, "p"] <- p
  results[ind, "r"] <- r
  results[ind, "bcdc"] <- unlist(res[1, ])
  results[ind, "casc"] <- unlist(res[2, ])
  results[ind, "kmeans"] <- unlist(res[3, ])
  results[ind, "SC"] <- unlist(res[4, ])
  results[ind, "bsbm"] <- unlist(res[5, ])
  
}

# Visualize ----
res <- tibble(results) %>% 
  pivot_longer(-c(n,p,r), names_to = "method", values_to = "nmi") %>% 
  mutate(method = factor(method)) %>% 
  mutate(r = signif(r, 1))

levels(res$method) <- list(BCDC = "bcdc", BSBM = "bsbm" , CASC = "casc", SC = "SC",  `k-means` = "kmeans")

res %>%
  ggplot(aes(x = factor(r), y = nmi, fill = method)) +
  geom_boxplot() +
  ylab("NMI") + xlab("r") +
  theme(legend.position = "none") +
  theme_minimal(base_size = 15)
