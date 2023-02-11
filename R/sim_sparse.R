# Libraries ----
library(parallel)
library(doParallel)
library(NMI)
library(tidyverse)
# devtools::install_github("aaamini/bcdc/bcdc_package")
library(bcdc)
library(CASCORE)

# Functions ----
source("./R/inference.R")
source("./R/CASC/cascTuningMethods.R")

source("./R/data_gen.R")

scaled = FALSE
source("./R/setup_methods.R")

# Simulation ----
K <- 2
n_iter <- 1000
n_reps <- 100
n_cores <- 45 # detectCores()

runs <- expand.grid(n = 800, rep = seq_len(n_reps))

progress_file = "progress.log"
system(sprintf("> %s", progress_file))

res <- do.call(rbind, mclapply(seq_len(nrow(runs)), function(ri) {
  set.seed(ri)
  
  rep <- runs[ri,"rep"]
  n <- runs[ri, "n"]
  
  sim <- sparse_sim_data(n)
  A <- sim$A
  Xc <- sim$Xc
  Xd = NULL
  z_true <- sim$z_true
  
  out = do.call(rbind, lapply(seq_along(methods), function(j) {
    start_time = Sys.time()
    zh <- as.vector(methods[[j]](A, Xc, Xd, K))
    end_time = as.numeric(Sys.time() - start_time)
    cat(sprintf("\n%10s (n = %4d) done in %3.2g (s)", mtd_names[j], n, end_time))
    data.frame(
      time = end_time
      , rep = rep
      , n = n
      , nmi = nmi_wrapper(z_true, zh)
      , method = mtd_names[j])
  }))
  system(sprintf('echo -n "." >> %s', progress_file))
  
  out 
}, mc.cores = n_cores))

res <- res %>%
  mutate(method = factor(method
                         , levels = c("BCDC", "BSBM", "CASC", "CASCORE", "SC", "k-means")))

save(res, file = sprintf("sparse_results_nreps_%d.RData", n_reps))

# Visualize ----
n_reps = length(unique(res$rep)) # redefine n_reps based on "res"
res %>%
  ggplot(aes(x = method, y = nmi, fill = method)) +
  geom_boxplot() +
  ylab("NMI") + xlab("") +
  guides(fill = "none") +
  scale_y_continuous(breaks = round(seq(0.3, 0.9, by = 0.1),1)) +
  theme_minimal(base_size = 15)

ggsave(sprintf("sim_sparse_nmi_nreps_%d.pdf", n_reps), width = 6, height = 6)

res %>%
  ggplot(aes(x = method, y = time, fill = method)) +
  geom_boxplot() +
  ylab("Seconds") + xlab("") +
  guides(fill = "none") +
  theme_minimal(base_size = 15)

ggsave(sprintf("sim_sparse_time_nreps_%d.pdf", n_reps), width = 6, height = 6)

res %>% 
  group_by(method) %>% 
  summarise(mean_nmi = mean(nmi))
