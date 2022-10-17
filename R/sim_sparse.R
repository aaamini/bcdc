# Libraries ----
library(pbmcapply)
library(NMI)
library(tidyverse)
# devtools::install_github("aaamini/bcdc/bcdc_package")
library(bcdc)

# Functions ----
source("./R/inference.R")
source("./R/CASC/cascTuningMethods.R")

simdata <- function(n, p, r, beta) {
  # n       : number of nodes
  # p       : prob of edge within cluster
  # r       : p*r = prob of edge between clusters
  # beta    : homophily parameter
  
  require(igraph)
  require(abind)
  
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
  
  return(list(A = as(A, "dgCMatrix"), Cov = cbind(feature), num = num$cluster))
}

# Simulation ----
K <- 3
n <- 150
p <- 0.3
r <- .35
beta <- seq(-.1, .1, length.out = 15)
n_iter <- 1500
n_reps <- n_cores <- detectCores()

grid <- expand.grid(n, p, r, beta)
names(grid) <- c("n", "p", "r", "beta")

results <- as.data.frame(matrix(nrow = nrow(grid) * n_reps, ncol = 8))
names(results) <- c("n", "p", "r", "beta"
                    , "bcdc", "casc", "SC", "bsbm")

for (tt in seq_len(nrow(grid))) {
  
  print(tt)
  
  n <- grid[tt, "n"]
  p <- grid[tt, "p"]
  r <- grid[tt, "r"]
  beta <- grid[tt, "beta"]
  
  res <- simplify2array(
    pbmclapply(seq_len(n_reps), mc.cores = n_cores, function(i)
    {
      set.seed(i)
      
      sim <- simdata(n, p, r, beta)
      A <- sim$A
      Cov <- sim$Cov
      Z_true <- sim$num
      
      # BCDC
      bcdc_model <- new(CovarSBM, A, alpha = 1, beta = 1, dp_concent = 10)
      bcdc_model$set_discrete_features(Cov)
      Z_bcdc <- get_map_labels(bcdc_model$run_gibbs(n_iter))$z_map
      
      # CASC
      Z_casc <- kmeans(getCascAutoSvd(A, Cov, 1, enhancedTuning = TRUE)$singVec
                       , centers = 2*K, nstart = 20)$cluster
      
      # SC
      Z_SC <- nett::spec_clust(A, 2*K)
      
      # BSBM
      bsbm <- new(SBM, A, 2*K, alpha = 1, beta = 1)
      Z_bsbm <- get_map_labels(bsbm$run_gibbs(n_iter))$z_map
      
      return(c(nmi_wrapper(Z_true, Z_bcdc)
               , nmi_wrapper(Z_true, Z_casc)
               , nmi_wrapper(Z_true, Z_SC)
               , nmi_wrapper(Z_true, Z_bsbm)))
      
    })
  )
  
  ind <- seq_len(n_reps) + n_reps * (tt - 1)
  
  results[ind, "n"] <- n
  results[ind, "p"] <- p
  results[ind, "r"] <- r
  results[ind, "beta"] <- beta
  results[ind, "bcdc"] <- unlist(res[1, ])
  results[ind, "casc"] <- unlist(res[2, ])
  results[ind, "SC"] <- unlist(res[3, ])
  results[ind, "bsbm"] <- unlist(res[4, ])
  
}

# Visualize ----
res <- tibble(results) %>% 
  pivot_longer(-c(n,p,r,beta), names_to = "method", values_to = "nmi") %>% 
  mutate(method = factor(method)) %>% 
  mutate(r = signif(r, 1))

levels(res$method) <- list(BCDC = "bcdc", BSBM = "bsbm" , CASC = "casc", SC = "SC",  `k-means` = "kmeans")

res %>%
  ggplot(aes(x = factor(round(beta, 2)), y = nmi, fill = method)) +
  geom_boxplot() +
  ylab("NMI") + xlab("beta") +
  theme(legend.position = "none") +
  theme_minimal(base_size = 15)
