# Libraries ----
library(pbmcapply)
library(NMI)
library(tidyverse)
# devtools::install_github("aaamini/bcdc/bcdc_package")
library(bcdc)


# Functions ----
source("./R/inference.R")
source("./R/CASC/cascTuningMethods.R")

source("./R/data_gen.R")
source("./R/setup_methods.R")

# Simulation ----
n <- seq(300, 1000, 200)
p <- 0.3
r <- .35
n_iter <- 1500
n_reps <- n_cores <- 3 # detectCores()


runs <- expand.grid(n = n, p = p, r = r, rep = seq_len(n_reps))
runs$K <- runs$n / 50

res <- do.call(rbind, mclapply(seq_len(nrow(runs)), function(ri) {r
  set.seed(ri)
  cat('.')
  
  rep <- runs[ri,"rep"]
  n <- runs[ri, "n"]
  p <- runs[ri, "p"]
  r <- runs[ri, "r"]
  K <- runs[ri, "K"]
  
  sim <- scale_sim_data(n, p, r, K)
  A <- sim$A
  Xc = sim$Xc
  Xd = sim$Xd
  z_true <- sim$z_true

  do.call(rbind, lapply(seq_along(methods), function(j) { 
    start_time = Sys.time()
    zh <- as.vector(methods[[j]](A, Xc, Xd, K))
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
