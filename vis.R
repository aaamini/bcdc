library(tidyverse)

fig_dir = "Images"

if (!exists("res")) {
  # stop("Need a `res` variable loaded in the workspace containing the results.")
  load("Data/cat1_results.RData")
}
if (!dir.exists(fig_dir)) {
  dir.create(fig_dir)
}

base_name = file.path(fig_dir, gsub("\\.", "p", sprintf("cat1_%d_%2.1f", res$n[1], res$p[1])))
fig1_name = paste0(base_name, ".pdf")
fig2_name = paste0(base_name, "_mean.pdf")

res %>% 
  # group_by(method, mu) %>% summarise(nmi = mean(nmi)) %>% 
  ggplot(aes(x = as.factor(r), y = nmi, fill = method)) + 
  geom_boxplot() +
  theme_minimal() +
  ggplot2::theme(
    legend.background = ggplot2::element_blank(),
    legend.title = ggplot2::element_blank(),
    legend.position = c(0.9, 0.8),
    # legend.text = ggplot2::element_text(size=18),
  ) + 
  ggplot2::guides(colour = ggplot2::guide_legend(keywidth = 2, keyheight = .75)) +
  ylab("NMI") + xlab("r") 
  # labs(title = sprintf("n = %d,  p = %2.1f,  r = %2.1f", n, p, r))
ggsave(fig1_name, width = 6, height = 6)


res %>% 
  group_by(method, r) %>% summarise(nmi = mean(nmi)) %>% 
  ggplot(aes(x = r, y = nmi, color = method)) + 
  geom_line(size = 1.2) +
  theme_minimal() +
  ggplot2::theme(
    legend.background = ggplot2::element_blank(),
    legend.title = ggplot2::element_blank(),
    legend.position = c(0.9, 0.8),
    # legend.text = ggplot2::element_text(size=18),
  ) + 
  ggplot2::guides(colour = ggplot2::guide_legend(keywidth = 2, keyheight = .75)) +
  ylab("NMI") + xlab("r") 
# labs(title = sprintf("n = %d,  p = %2.1f,  r = %2.1f", n, p, r))
ggsave(fig2_name, width = 6, height = 6)

res %>% 
  group_by(method) %>% 
  summarise(dt = mean(time)) %>% 
  knitr::kable(digits = 3)

