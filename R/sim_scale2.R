# Libraries ----
library(pbmcapply)
library(NMI)
library(tidyverse)
library(tictoc)
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



methods <- list()

methods[["BSBM"]] <- function(A, X, K) {
  model <- new(SBM, A, K, alpha = 1, beta = 1)
  get_map_labels(model$run_gibbs(n_iter))$z_map
}

methods[["k-means"]] <- function(A, X, K) {
  kmeans(X, centers = K, nstart = 20)$cluster
}

methods[["SC"]] <- function(A, X, K) {
  nett::spec_clust(A, K)
}

methods[["BCDC"]] <- function(A, X, K) {
  bcdc_model <- new(CovarSBM, A, alpha = 10, beta = 1, dp_concent = 10)
  bcdc_model$set_discrete_features(X[, 1, drop = FALSE])
  bcdc_model$set_continuous_features(X[, 2, drop = FALSE])
  get_map_labels(bcdc_model$run_gibbs(n_iter))$z_map
}

# methods[["CASC"]] <- function(A, X, K) {
#   kmeans(getCascAutoSvd(A, scale(X), K, enhancedTuning = TRUE)$singVec
#          , centers = K, nstart = 20)$cluster
# }

mtd_names <- names(methods)


runs <- expand.grid(n = n, p = p, r = r, rep = seq_len(n_reps))
runs$K <- runs$n / 50

res <- do.call(rbind, pbmclapply(seq_len(nrow(runs)), function(ri) {r
  set.seed(ri)
  
  rep <- runs[ri,"rep"]
  n <- runs[ri, "n"]
  p <- runs[ri, "p"]
  r <- runs[ri, "r"]
  K <- runs[ri, "K"]
  
  sim <- simdata(n, p, r, K)
  A <- sim$A
  X <- sim$Cov
  simnum <- sim$num
  z_true <- rep(1:K, simnum)

  do.call(rbind, lapply(seq_along(methods), function(j) { 
    start_time = Sys.time()
    zh <- as.vector(methods[[j]](A, X, K))
    data.frame(
      time = as.numeric(Sys.time() - start_time)
      , rep = rep
      , n = n
      , p = p
      , r = r
      , K = K
      , nmi = nmi_wrapper(z_true, zh)
      , method = mtd_names[j])
  }))
}, mc.cores = n_cores))

res <- res %>%
  mutate(method = factor(method
                         , levels = c("BCDC", "BSBM", "CASC", "SC", "k-means")))


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

res %>% 
  ggplot(aes(x = factor(n), y = time, fill = method)) +
  geom_boxplot() +
  ylab("Seconds") + xlab("n") +
  guides(fill = "none") +
  theme_minimal(base_size = 15)
