# Libraries ----
library(pbmcapply)
library(NMI)
library(tidyverse)
# devtools::install_github("aaamini/bcdc/bcdc_package")
library(bcdc)

# Functions ----
source("./R/inference.R")
source("./R/CASC/cascTuningMethods.R")

simdata <- function(n, p, r, K) {
  # n       : number of nodes
  # p       : prob of edge within cluster
  # r       : p*r = prob of edge between clusters
  # K       : number of communities
  
  require(igraph)
  require(abind)
  
  P <- matrix(p*r, K, K)
  diag(P) <- p
  
  n_i <- rep(n/K, K)
  A <- get.adjacency(sample_sbm(n, P, n_i))
  
  feature <- cbind(c(rep(0:1, c(n/2, n/2)))
                   , rnorm(n, rep(2*(1:K)%%2-1, each = n/K)))
  
  return(list(A = A, Cov = feature, num = n_i))
}

# Simulation ----
n <- seq(300, 1000, 100)
p <- 0.3
r <- .35
n_iter <- 1500
n_reps <- n_cores <- detectCores()

grid <- expand.grid(n, p, r)
names(grid) <- c("n", "p", "r")
grid$K <- grid$n / 50

results <- as.data.frame(matrix(nrow = nrow(grid) * n_reps, ncol = 14))
res_names <- c("bcdc", "casc", "kmeans", "SC", "bsbm")
names(results) <- c("n", "p", "r", "K"
                    , res_names, paste0("time_", res_names))

for (tt in seq_len(nrow(grid))) {
  
  print(tt)
  
  n <- grid[tt, "n"]
  p <- grid[tt, "p"]
  r <- grid[tt, "r"]
  K <- grid[tt, "K"]
  
  res <- simplify2array(
    pbmclapply(seq_len(n_reps), mc.cores = n_cores, function(i)
    {
      set.seed(i)
      
      sim <- simdata(n, p, r, K)
      A <- sim$A
      Cov <- sim$Cov
      simnum <- sim$num
      
      Z_true <- rep(1:K, simnum)
      
      # BCDC
      start.time <- Sys.time()
      bcdc_model <- new(CovarSBM, A, alpha = 10, beta = 1, dp_concent = 10)
      bcdc_model$set_discrete_features(Cov[, 1, drop = FALSE])
      bcdc_model$set_continuous_features(Cov[, 2, drop = FALSE])
      Z_bcdc <- get_map_labels(bcdc_model$run_gibbs(n_iter))$z_map
      time_bcdc <- Sys.time() - start.time
      
      # CASC
      start.time <- Sys.time()
      Z_casc <- kmeans(getCascAutoSvd(A, Cov, 1, enhancedTuning = TRUE)$singVec
                       , centers = K, nstart = 20)$cluster
      time_casc <- Sys.time() - start.time
      
      # k-means
      start.time <- Sys.time()
      Z_kmeans <- kmeans(Cov, centers = K, nstart = 20)$cluster
      time_kmeans <- Sys.time() - start.time
      
      # SC
      start.time <- Sys.time()
      Z_SC <- nett::spec_clust(A, K)
      time_SC <- Sys.time() - start.time
      
      # BSBM
      start.time <- Sys.time()
      bsbm <- new(SBM, A, K, alpha = 1, beta = 1)
      Z_bsbm <- get_map_labels(bsbm$run_gibbs(n_iter))$z_map
      time_bsbm <- Sys.time() - start.time
      
      return(c(nmi_wrapper(Z_true, Z_bcdc)
               , nmi_wrapper(Z_true, Z_casc)
               , nmi_wrapper(Z_true, Z_kmeans)
               , nmi_wrapper(Z_true, Z_SC)
               , nmi_wrapper(Z_true, Z_bsbm)
               , time_bcdc
               , time_casc
               , time_kmeans
               , time_SC
               , time_bsbm))
      
    })
  )
  
  ind <- seq_len(n_reps) + n_reps * (tt - 1)
  
  results[ind, "n"] <- n
  results[ind, "p"] <- p
  results[ind, "r"] <- r
  results[ind, "K"] <- K
  results[ind, "bcdc"] <- unlist(res[1, ])
  results[ind, "casc"] <- unlist(res[2, ])
  results[ind, "kmeans"] <- unlist(res[3, ])
  results[ind, "SC"] <- unlist(res[4, ])
  results[ind, "bsbm"] <- unlist(res[5, ])
  results[ind, "time_bcdc"] <- unlist(res[6, ])
  results[ind, "time_casc"] <- unlist(res[7, ])
  results[ind, "time_kmeans"] <- unlist(res[8, ])
  results[ind, "time_SC"] <- unlist(res[9, ])
  results[ind, "time_bsbm"] <- unlist(res[10, ])
  
}

# Visualize ----
res <- tibble(results[, 1:9]) %>% 
  pivot_longer(-c(n,p,r,K), names_to = "method", values_to = "nmi") %>% 
  mutate(method = factor(method)) %>% 
  mutate(r = signif(r, 1))

levels(res$method) <- list(BCDC = "bcdc", BSBM = "bsbm" , CASC = "casc", SC = "SC",  `k-means` = "kmeans")

res %>%
  ggplot(aes(x = factor(n), y = nmi, fill = method)) +
  geom_boxplot() +
  ylab("NMI") + xlab("n") +
  guides(fill = "none") +
  theme_minimal(base_size = 15)

res_time <- tibble(results[, c(1:4, 10:14)]) %>% 
  pivot_longer(-c(n,p,r,K), names_to = "method", values_to = "time") %>% 
  mutate(method = factor(method)) %>% 
  mutate(r = signif(r, 1))

levels(res_time$method) <- list(BCDC = "time_bcdc", BSBM = "time_bsbm" , CASC = "time_casc", SC = "time_SC",  `k-means` = "time_kmeans")

res_time %>% 
  ggplot(aes(x = factor(n), y = time, fill = method)) +
  geom_boxplot() +
  ylab("Seconds") + xlab("n") +
  guides(fill = "none") +
  theme_minimal(base_size = 15)
