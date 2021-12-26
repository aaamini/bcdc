load("BCDC/sim_cont.RData")

library(tidyverse)
res = tibble(results) %>% 
  pivot_longer(-c(n,p,r,mu), names_to = "method", values_to = "nmi") %>% 
  mutate(method = factor(method)) %>% 
  mutate(r = signif(r, 1))

# Rename and reorder levels
levels(res$method) = list(BCDC = "bcdc", BSBM = "bsbm" , CASC = "casc", SC = "SC",  `k-means` = "kmeans")
res$method

# rvec = unique(res$r)[-1]
rvec = c(0.3, 0.4, 0.5, 0.7)
for (curr_r in rvec) {
  plt = res %>% 
    filter(p == 0.1 & r == curr_r) %>%
    mutate(r = sprintf("r = %2.1f", r)) %>% 
    # filter(p == 0.1) %>%
    # pivot_longer(-c(n,p,r,mu), names_to = "method", values_to = "nmi") %>% 
    ggplot(aes(x = factor(mu), y = nmi, fill = method)) +
    geom_boxplot() +
    theme_minimal(base_size = 15) +
    ylab("NMI") + xlab(expression(mu)) +
    facet_wrap(~r)
  
  # if (curr_r == rvec[1]) {
  if (curr_r == .5) {
    plt = plt  + theme(
      legend.background = ggplot2::element_blank(),
      legend.title = ggplot2::element_blank(),
      # legend.position = c(0.9, 0.3),
      legend.position = c(0.2, 0.8),
      # legend.text = ggplot2::element_text(size=18),
    )
  } else {
    plt = plt + theme(legend.position = "none")
  }
  plt
  #labs(title = sprintf("r = %2.1f",  curr_r))
  ggsave(sprintf("sim_cont_%1d.pdf", round(10*curr_r)), width = 5, height = 5)
}



### Sparse sim
load("BCDC/sim_sparse.RData")
res = tibble(results) %>% 
  pivot_longer(everything(), names_to = "method", values_to = "nmi") %>% 
  mutate(method = factor(method))

# Rename and reorder levels
levels(res$method) = list(BCDC = "bcdc", BSBM = "bsbm" , CASC = "casc", SC = "SC",  `k-means` = "kmeans")
res$method

plt = res %>% 
  ggplot(aes(x = method, y = nmi, fill = method)) +
  geom_boxplot() +
  theme_minimal(base_size = 15) +
  ylab("NMI") + theme(legend.position = "none") + xlab("")
plt
ggsave("sim_sparse.pdf", width = 5, height = 5)
       