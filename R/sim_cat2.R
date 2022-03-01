# Libraries ----
library(pbmcapply)
library(NMI)
library(tidyverse)

# Functions ----
Rcpp::sourceCpp("./src/SBM.cpp", verbose = T)
Rcpp::sourceCpp("./src/CovarSBM.cpp", verbose = T)
source("./R/inference.R")
source("./R/CASC/cascTuningMethods.R")

simdata <- function(n) {
  # n       : number of nodes
  
  require(igraph)
  require(abind)
  
  # sample SBM ----
  n1 <- floor(n*3/12)
  n2 <- floor(n*4/12)
  n3 <- n - n1 - n2
  
  P <- cbind(c(1.6, 1.2, 0.16)
             , c(1.2, 1.6, 0.02)
             , c(0.16, 0.02, 1.2)) * 0.01
  A <- get.adjacency(sample_sbm(n, P, c(n1, n2, n3)))
  
  # sample features ----
  feature <- cbind(c(rnorm(n1, 0, 1), rnorm(n2, -1, 1), rnorm(n3, 1, 1))
                   , c(rnorm(n1, 2, 1), rnorm(n2, -0.8, 1), rnorm(n3, -0.8, 1))
                   , matrix(rnorm(n*98, 0, 1), n, 98)) # noise
  
  return(list(A = A, Cov = feature, num = c(n1, n2, n3)))
}

# Simulation ----
K <- 3
n <- 800
n_iter <- 1000
n_reps <- n_cores <- detectCores()

results <- as.data.frame(matrix(nrow = n_reps, ncol = 5))
names(results) <- c("bcdc", "casc", "kmeans", "SC", "bsbm")

res <- simplify2array(
  pbmclapply(seq_len(n_reps), mc.cores = n_cores, function(i)
  {
    set.seed(i)
    
    sim <- simdata(n)
    A <- sim$A
    Cov <- sim$Cov
    simnum <- sim$num
    
    Z_true <- rep(1:K, simnum)
    
    # BCDC
    bcdc <- new(CovarSBM, A, alpha = 1, beta = 1, dp_concent = 10)
    bcdc$set_continuous_features(Cov)
    Z_bcdc <- get_map_labels(bcdc$run_gibbs(n_iter))$z_map
    
    # CASC
    Z_casc <- kmeans(getCascAutoSvd(A, Cov, K, enhancedTuning = TRUE)$singVec
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

results[, "bcdc"] <- unlist(res[1, ])
results[, "casc"] <- unlist(res[2, ])
results[, "kmeans"] <- unlist(res[3, ])
results[, "SC"] <- unlist(res[4, ])
results[, "bsbm"] <- unlist(res[5, ])

# Visualize ----
res <- tibble(results) %>% 
  pivot_longer(everything(), names_to = "method", values_to = "nmi") %>% 
  mutate(method = factor(method))

levels(res$method) <- list(BCDC = "bcdc", BSBM = "bsbm" , CASC = "casc", SC = "SC",  `k-means` = "kmeans")

res %>% 
  ggplot(aes(x = method, y = nmi, fill = method)) +
  geom_boxplot() +
  ylab("NMI") + theme(legend.position = "none") + xlab("") +
  theme_minimal(base_size = 15)
