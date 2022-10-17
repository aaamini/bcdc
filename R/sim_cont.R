# Libraries ----
library(pbmcapply)
library(tictoc)
library(dplyr)
library(ggplot2)
# devtools::install_github("aaamini/bcdc/bcdc_package")
library(bcdc)

# Functions ----
source("./R/inference.R")
source("./R/CASC/cascTuningMethods.R")

# Simulation ----
n <- 150
K <- 2
alpha <- beta <- 1
dp_concent <- 10

niter <- 1000
nreps <- n_cores <- detectCores()

r <- 0.3
p <- 0.1
eta <- matrix(c(p, r*p, r*p, p), nrow = 2)

methods <- list()

methods[["BSBM"]] <- function(A, X) {
  model <- new(SBM, A, K, alpha, beta)
  get_map_labels(model$run_gibbs(niter))$z_map
}

methods[["k-means"]] <- function(A, X) {
  kmeans(X, K, nstart = 20)$cluster
}

methods[["SC"]] <- function(A, X) {
  nett::spec_clust(A, K)
}

methods[["BCDC"]] <- function(A, X) {
  bcdc <- new(CovarSBM, A, alpha, beta, dp_concent)
  bcdc$set_continuous_features(X)
  get_map_labels(bcdc$run_gibbs(niter))$z_map
}

methods[["CASC"]] <- function(A, X) {
  kmeans(getCascAutoSvd(A, X, K, enhancedTuning = TRUE)$singVec
         , centers = K, nstart = 20)$cluster
}

mtd_names <- names(methods)

runs <- expand.grid(mu = seq(0, 2, by = 0.25), rep = seq_len(nreps))

res <- do.call(rbind, pbmclapply(seq_len(nrow(runs)), function(ri) {r
  set.seed(ri)
  
  rep <- runs[ri,"rep"]
  mu <- runs[ri,"mu"]
  
  z_tru <- sample(1:K, n, replace = T, prob = c(1, 2))
  A <- nett::fast_sbm(z_tru, eta)
  
  Mu <- rbind(c(mu, 0), c(-mu,0))
  X <- do.call(cbind, lapply(1:ncol(Mu), function(j) {
    rnorm(n, mean = Mu[z_tru,j])
  }))
  
  do.call(rbind, lapply(seq_along(methods), function(j) { 
    tic()
    zh <- as.vector(methods[[j]](A, X)) + 1
    dt <- toc(quiet = T)
    data.frame(
      dt = as.numeric(dt$toc - dt$tic)
      , rep = rep
      , mu = mu
      , nmi = nmi_wrapper(z_tru, zh)
      , method = mtd_names[j])
  }))
}, mc.cores = n_cores))

res <- res %>%
  mutate(method = factor(method
                         , levels = c("BCDC", "BSBM", "CASC", "SC", "k-means")))

res %>%
  ggplot(aes(x = as.factor(mu), y = nmi, fill = method)) + 
  geom_boxplot() +
  theme(legend.background = element_blank()
        , legend.title = element_blank()
        , legend.position = c(0.2, 0.4)) +
  guides(colour = guide_legend(keywidth = 2, keyheight = .75)) +
  ylab("NMI") + xlab("mu") +
  labs(title = sprintf("r = %2.1f", r)) +
  theme_minimal(base_size = 15)
