library(igraph)
library(nett)
library(tidyverse)

# df = read.delim("data/temp/511145.protein.links.v11.5.txt", sep=" ")
df = read.delim("data/pi_nets//9606.protein.links.v11.5.txt.gz", sep=" ")
# df = readr::read_table("data/pi_nets//9606.protein.links.v11.5.txt.gz")
df = data.table::fread("data/pi_nets//9606.protein.links.v11.5.txt.gz", sep = " ")
hist(df$combined_score)
df = tibble(df)

g = graph_from_data_frame(df %>% filter(combined_score > 900), directed = F)
g2 = induced_subgraph(g, coreness(g) > 3)
g2 = igraph::simplify(g2)
g2
plot_net(g2)

A = as_adj(g2)

library(parallel)
nreps = 2
niter = 200
zhs = do.call(cbind, mclapply(1:nreps, function(rep) {
  model = new(CovarSBM, A, 1, 1, 1)
  # model$set_node_features(U) 
  tic()
  zout = model$run_gibbs(niter) # 
  toc()
  get_map_labels(zout)$z_map
}))

zhs
  
# model = new(CovarSBM, A, 1, 1, 1)
# # model$set_node_features(U) 
# tic()
# zout = model$run_gibbs(niter) # 
# toc()
# # nett::compute_mutual_info(z, model$z+1)
# zh = get_map_labels(zout)$z_map
# # nett::compute_mutual_info(z, zh)
# table(model$z)
# table(zh)

# zh1 = zh
# zh2 = zh
# compute_mutual_info(zh, zh2)

V(g2)$comm = zhs[,1]
out = nett::plot_net(g2, extract_lcc = 1, community = V(g2)$comm, vertex_alpha = .3, niter = 1500)


data.frame(itr=1:niter, nmi=sapply(1:niter, function(itr) nett::compute_mutual_info(zout[,itr], zout[,itr+1]))) %>% 
  ggplot(aes(itr, nmi)) + 
  geom_line() +
  xlab("Iteration") + ylab("NMI")

# cls = data.table::fread("data/pi_nets//9606.clusters.proteins.v11.5.txt.gz", sep = "\t")
# cls = tibble(cls) %>% 
#   janitor::clean_names() %>% 
#   mutate_all(factor)
# cls$cluster_id
# unique(cls$cluster_id)
# 
# ctree = data.table::fread("data/pi_nets//9606.clusters.tree.v11.5.txt.gz", sep = "\t")
# ctree = ctree[,3:2]
# tree = graph_from_data_frame(ctree)
# plot_net(tree, coord = layout_as_tree(tree),  edge.arrow.size = .2)
