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

# Simulation ----
n <- seq(300, 1000, 100)
p <- 0.3
r <- .35
n_iter <- 1500
n_reps <- 500
n_cores <- detectCores()

runs <- expand.grid(n = n, p = p, r = r, rep = seq_len(n_reps))
runs$K <- runs$n / 50

progress_file = "progress.log"
system(sprintf("> %s", progress_file))

res <- do.call(rbind, mclapply(seq_len(nrow(runs)), function(ri) {
  set.seed(ri)
  
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
      , K = K
      , nmi = nmi_wrapper(z_true, zh)
      , method = mtd_names[j])
  }))
  system(sprintf('echo -n "." >> %s', progress_file))
  
  out 
}, mc.cores = n_cores))

res <- res %>%
  mutate(method = factor(method
                         , levels = c("BCDC", "BSBM", "CASC", "CASCORE", "SC", "k-means")))
save(res, file = "scale_results.RData")

# Visualize ----
mean_res =  res %>% 
  group_by(method, n) %>% 
  summarise(mean_nmi = mean(nmi), lower=quantile(nmi,.25), upper=quantile(nmi,.75), .groups="drop")

mean_res %>% 
  ggplot(aes(x = n, y = mean_nmi, color = method)) +
  geom_line(size = 1.2) +
  theme_minimal() +
  ggplot2::theme(
    legend.background = ggplot2::element_blank(),
    legend.title = ggplot2::element_blank(),
    legend.position = c(0.8, 0.825),
    # legend.text = ggplot2::element_text(size=18),
  ) +
  ggplot2::guides(colour = ggplot2::guide_legend(keywidth = 2, keyheight = .75)) +
  geom_ribbon(aes(ymin = lower, ymax=upper, fill= method), alpha= 0.1, linetype = "blank") +
  ylab("NMI") + xlab("r")

ggsave("sim_scale.pdf", width = 6, height = 6)

# Tried trimmed mean to avoid outliers, but there are too many outlies for CASC, etc.
mean_res =  res %>% 
  group_by(method, n) %>% 
  summarise(mean_time = mean(time, trim = 0.0), lower=quantile(time,.25), upper=quantile(time,.75), .groups="drop")

# "sqrt" scale for the y-axis to improve visualization
mean_res %>% 
  ggplot(aes(x = n, y = mean_time, color = method)) +
  geom_line(size = 1.2) +
  theme_minimal() +
  theme(legend.position="none") +
  ggplot2::guides(colour = ggplot2::guide_legend(keywidth = 2, keyheight = .75)) +
  geom_ribbon(aes(ymin = lower, ymax=upper, fill= method), alpha= 0.1, linetype = "blank") +
  ylab("Seconds") + xlab("n") +
  scale_y_continuous(trans="sqrt")

ggsave("sim_scale_time_sqrt_nreps_500.pdf", width = 6, height = 6)
