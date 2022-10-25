# Libraries ----
library(pbmcapply)
library(NMI)
library(tidyverse)
# devtools::install_github("aaamini/bcdc/bcdc_package")
library(bcdc)

# Functions ----
source("./R/inference.R")
source("./R/CASC/cascTuningMethods.R")

simdata <- function(n, p, r) {
  # n       : number of nodes
  # p       : prob of edge within cluster
  # r       : p*r = prob of edge between clusters
  
  require(igraph)
  require(MCMCpack)
  require(abind)
  
  # sample SBM ----
  n1 <- floor(n*1/3)
  n2 <- n - n1
  
  P <- matrix(p*r, 2, 2)
  diag(P) <- p
  A <- get.adjacency(sample_sbm(n, P, c(n1, n2)))
  
  # sample features ----
  dir_p1 <- MCMCpack::rdirichlet(1, rep(1, 4))
  dir_p2 <- MCMCpack::rdirichlet(1, rep(1, 4))
  feature <- cbind(c(sample(0:3, n1, prob = dir_p1, replace = TRUE)
                     , sample(0:3, n2, prob = dir_p2, replace = TRUE))
                   , sample(0:3, n, prob = rep(1/4, 4), replace = TRUE))
  
  return(list(A = A, Cov = feature, num = c(n1, n2)))
}

# Simulation ----
K <- 2
n <- 150
p <- 0.1
r <- seq(0.1, 0.8, 0.1)
n_iter <- 1500
n_reps <- n_cores <- detectCores()

grid <- expand.grid(n, p, r)
names(grid) <- c("n", "p", "r")

results <- as.data.frame(matrix(nrow = nrow(grid) * n_reps, ncol = 7))
names(results) <- c("n", "p", "r"
                    , "bcdc", "casc", "SC", "bsbm")

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
      bcdc_model <- new(CovarSBM, A, alpha = 1, beta = 1, dp_concent = 10)
      bcdc_model$set_discrete_features(Cov)
      Z_bcdc <- get_map_labels(bcdc_model$run_gibbs(n_iter))$z_map
      
      # CASC
      Z_casc <- kmeans(getCascAutoSvd(A, Cov, K, enhancedTuning = TRUE)$singVec
                       , centers = K, nstart = 20)$cluster
      
      # SC
      Z_SC <- nett::spec_clust(A, K)
      
      # BSBM
      bsbm <- new(SBM, A, K, alpha = 1, beta = 1)
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
  results[ind, "bcdc"] <- unlist(res[1, ])
  results[ind, "casc"] <- unlist(res[2, ])
  results[ind, "SC"] <- unlist(res[3, ])
  results[ind, "bsbm"] <- unlist(res[4, ])
  
}

# Visualize ----
res <- tibble(results) %>% 
  pivot_longer(-c(n,p,r), names_to = "method", values_to = "nmi") %>% 
  mutate(method = factor(method)) %>% 
  mutate(r = signif(r, 1))

levels(res$method) <- list(BCDC = "bcdc", BSBM = "bsbm" , CASC = "casc", SC = "SC")

res %>% 
  ggplot(aes(x = factor(r), y = nmi, fill = method)) +
  geom_boxplot() +
  ylab("NMI") + xlab("r") +
  guides(fill = "none") +
  theme_minimal(base_size = 15)
