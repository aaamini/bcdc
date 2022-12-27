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
source("./R/setup_methods.R")

methods[["k-means"]] <- NULL
mtd_names <- names(methods)

# Simulation ----
K <- 3
n <- 150
p <- 0.3
r <- .7
n_iter <- 1500
n_reps <- n_cores <- detectCores()

runs <- expand.grid(beta = seq(-.2, .2, by = .05), rep = seq_len(n_reps))

progress_file = "progress.log"
system(sprintf("> %s", progress_file))

res <- do.call(rbind, mclapply(seq_len(nrow(runs)), function(ri) {
  set.seed(ri)
  
  rep <- runs[ri,"rep"]
  beta <- runs[ri, "beta"]
  
  sim <- homophily_sim_data(n, p, r, beta)
  A <- sim$A
  Xc <- NULL
  Xd = sim$Xd
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
      , p = p
      , r = r
      , beta = beta
      , nmi = nmi_wrapper(z_true, zh)
      , method = mtd_names[j])
  }))
  system(sprintf('echo -n "." >> %s', progress_file))
  
  out 
}, mc.cores = n_cores))

res <- res %>%
  mutate(method = factor(method
                         , levels = c("BCDC", "BSBM", "CASC", "CASCORE", "SC")))

save(res, file = "homophily_results.RData")

# Visualize ----
res %>%
  ggplot(aes(x = factor(round(beta, 2)), y = nmi, fill = method, color = method)) +
  geom_boxplot() +
  ylab("NMI") + xlab(expression(beta)) +
  guides(fill = "none") +
  theme_minimal(base_size = 15)
