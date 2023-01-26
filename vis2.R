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
fig1_name = paste0(base_name, "_boxplot_plus_mean.pdf")
fig2_name = paste0(base_name, "_bands_plus_mean.pdf")

mean_res =  res %>% 
  group_by(method, r) %>% 
  summarise(mean_nmi = mean(nmi), lower=quantile(nmi,.25), upper=quantile(nmi,.75), .groups="drop")

res %>% 
  # group_by(method, mu) %>% summarise(nmi = mean(nmi)) %>% 
  ggplot(aes(x = r, y = nmi, group = interaction(r, method), fill = method)) + 
  geom_boxplot(alpha=1) +
  theme_minimal() +
  ggplot2::theme(
    legend.background = ggplot2::element_blank(),
    legend.title = ggplot2::element_blank(),
    legend.position = c(0.9, 0.8),
    # legend.text = ggplot2::element_text(size=18),
  ) + 
  ggplot2::guides(colour = ggplot2::guide_legend(keywidth = 2, keyheight = .75)) +
  ylab("NMI") + xlab("r") +
  geom_path(aes(x = r, y = mean_nmi, group = method, color = method), 
            data=mean_res,
            alpha=1,
            size = 1.3)
   # labs(title = sprintf("n = %d,  p = %2.1f,  r = %2.1f", n, p, r))
ggsave(fig1_name, width = 6, height = 6)


mean_res %>% 
  ggplot(aes(x = r, y = mean_nmi, color = method)) +
  geom_line(size = 1.2) +
  theme_minimal() +
  ggplot2::theme(
    legend.background = ggplot2::element_blank(),
    legend.title = ggplot2::element_blank(),
    legend.position = c(0.9, 0.8),
    # legend.text = ggplot2::element_text(size=18),
  ) +
  ggplot2::guides(colour = ggplot2::guide_legend(keywidth = 2, keyheight = .75)) +
  geom_ribbon(aes(ymin = lower, ymax=upper, fill= method), alpha= 0.1, linetype = "blank") +
  ylab("NMI") + xlab("r")
ggsave(fig2_name, width = 6, height = 6)

res %>% 
  group_by(method) %>% 
  summarise(dt = mean(time)) %>% 
  knitr::kable(digits = 3)

