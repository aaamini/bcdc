# Libraries ----
library(pbmcapply)
library(NMI)
library(tidyverse)
# devtools::install_github("aaamini/bcdc/bcdc_package")
library(bcdc)

# Functions ----
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

results <- as.data.frame(matrix(nrow = n_reps, ncol = 10))
res_names <- c("bcdc", "casc", "kmeans", "SC", "bsbm")
names(results) <- c(res_names, paste0("time_", res_names))

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
    start.time <- Sys.time()
    bcdc <- new(CovarSBM, A, alpha = 1, beta = 1, dp_concent = 10)
    bcdc$set_continuous_features(Cov)
    Z_bcdc <- get_map_labels(bcdc$run_gibbs(n_iter))$z_map
    time_bcdc <- Sys.time() - start.time
    
    # CASC
    start.time <- Sys.time()
    Z_casc <- kmeans(getCascAutoSvd(A, Cov, K, enhancedTuning = TRUE)$singVec
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

results[, "bcdc"] <- unlist(res[1, ])
results[, "casc"] <- unlist(res[2, ])
results[, "kmeans"] <- unlist(res[3, ])
results[, "SC"] <- unlist(res[4, ])
results[, "bsbm"] <- unlist(res[5, ])
results[, "time_bcdc"] <- unlist(res[6, ])
results[, "time_casc"] <- unlist(res[7, ])
results[, "time_kmeans"] <- unlist(res[8, ])
results[, "time_SC"] <- unlist(res[9, ])
results[, "time_bsbm"] <- unlist(res[10, ])

# Visualize ----
res <- tibble(results[, 1:5]) %>% 
  pivot_longer(everything(), names_to = "method", values_to = "nmi") %>% 
  mutate(method = factor(method))

levels(res$method) <- list(BCDC = "bcdc", BSBM = "bsbm" , CASC = "casc", SC = "SC",  `k-means` = "kmeans")

res %>% 
  ggplot(aes(x = method, y = nmi, fill = method)) +
  geom_boxplot() +
  ylab("NMI") + theme(legend.position = "none") + xlab("") +
  guides(fill = "none") +
  theme(text = element_text(size = 20))

res_time <- tibble(results[, 6:10]) %>% 
  pivot_longer(everything(), names_to = "method", values_to = "time") %>% 
  mutate(method = factor(method))

levels(res_time$method) <- list(BCDC = "time_bcdc", BSBM = "time_bsbm" , CASC = "time_casc", SC = "time_SC",  `k-means` = "time_kmeans")

res_time %>% 
  ggplot(aes(x = method, y = time, fill = method)) +
  geom_boxplot() +
  ylab("Seconds") + theme(legend.position = "none") + xlab("") +
  guides(fill = "none") +
  theme(text = element_text(size = 20))
