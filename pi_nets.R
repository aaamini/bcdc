library(igraph)
library(nett)
library(tidyverse)

# # df = read.delim("data/temp/511145.protein.links.v11.5.txt", sep=" ")
# df = read.delim("data/pi_nets//9606.protein.links.v11.5.txt.gz", sep=" ")
# # df = readr::read_table("data/pi_nets//9606.protein.links.v11.5.txt.gz")
df = data.table::fread("data/pi_nets//9606.protein.links.v11.5.txt.gz", sep = " ")
hist(df$combined_score)
df = tibble(df)

# g = graph_from_data_frame(df %>% filter(combined_score > 900), directed = F)
g = graph_from_data_frame(df %>% filter(combined_score > 975), directed = F)
g2 = induced_subgraph(g, coreness(g) > 3)
g2 = extract_largest_cc(g2)
g2 = igraph::simplify(g2)
g2
plot_net(g2, niter = 2000)

# hist(degree(g2))
plot_deg_dist(g2, logx=F)
A = as_adj(g2)

library(parallel)
library(foreach)
library(doParallel)
nreps = 32
niter = 1000
n_cores = 32
cl <- makeCluster(n_cores)
registerDoParallel(cl)

zout_list = foreach(ri = 1:nreps, .combine = 'c') %dopar% {
# zhs = do.call(cbind, mclapply(1:nreps, function(rep) {
  Rcpp::sourceCpp("src/CovarSBM.cpp")
  model = new(CovarSBM, A, 1, 1, 1)
  # model$set_node_features(U) 
  # tic()
  zout = model$run_gibbs(niter) # 
  # toc()
  # get_map_labels(zout)$z_map
  list(zout)
}
stopCluster(cl) 
# }, mc.cores = 32))



sort_labels = function(z) {
  z_new = z
  sorted_comms = sort(table(z), dec = T)
  old_labels = as.integer(names(sorted_comms))
  old_labels
  for (i in seq_along(old_labels)) {
    z_new[z == old_labels[i]] = i
  }
  z_new
}

zhs = sapply(zout_list, function(zout) get_map_labels(zout)$z_map)
zhs2 = apply(zhs, 2, sort_labels)

apply(zhs2, 2, function(zh) table(zh))

nmi_mat = apply(zhs, 2, function(y) apply(zhs, 2, function(z) nett::compute_mutual_info(z, y)))
nmi_seq = nmi_mat[upper.tri(nmi_mat)]
png("nmi_hist.png")
hist(nmi_seq)
dev.off()


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
# table(zhs[,3])

# mar = par()$mar
for (rep in 1:nreps) {
  png(sprintf("figs/pi_net_%d.png",rep), width = 1000, height = 1000)
  par(mar = rep(0,4))
  V(g2)$comm = zhs2[, rep]
  if (rep == 1) {
    out = nett::plot_net(g2, extract_lcc = 1, community = V(g2)$comm, vertex_alpha = .3, niter = 2000)  
    coord = out$coord
  } else {
    out = nett::plot_net(g2, extract_lcc = 1, community = V(g2)$comm, vertex_alpha = .3, niter = 2000, coord = coord)  
  }
  dev.off()
}

get_seq_nmi = function(zout) {
  sapply(1:(ncol(zout)-1), function(itr) nett::compute_mutual_info(zout[,itr], zout[,itr+1]))
}

cl <- makeForkCluster(n_cores)
registerDoParallel(cl)
seq_nmi = foreach(ri = 1:nreps, .combine = 'rbind') %dopar% {
 data.frame(itr=1:niter, nmi = get_seq_nmi(zout_list[[ri]]), rep = ri)
}
stopCluster(cl) 

#zout = zout_list[[1]]
#data.frame(itr=1:niter, nmi=get_seq_nmi(zout)) %>% 
seq_nmi %>% 
 #  filter(rep == 3) %>% 
  ggplot(aes(itr, nmi)) + 
  geom_line(aes(group = rep), alpha = .2) + 
  xlab("Iteration") + ylab("NMI")
ggsave("figs/pi_nets_seq_nmi.png", width = 6, height = 5)

zh_sc = sort_labels(spec_clust(A, 8))
table(zh_sc)
png("figs/pi_net_SC.png", width = 1000, height = 1000)
V(g2)$comm = zh_sc
out = nett::plot_net(g2, extract_lcc = 1, community = V(g2)$comm, vertex_alpha = .3, niter = 2000, coord = coord)  
dev.off()

# tstat = snac_resample(A, nrep = 10, ncores = 3)
# plot_smooth_profile(tstat)

(like = apply(zhs2, 2, function(z) nett::eval_dcsbm_like(A, z)))
bic = apply(zhs2, 2, function(z) nett::eval_dcsbm_bic(A, z, K = length(unique(z)), poi = T))
(max_like_idx = which.max(like))
which.max(bic)

nett::eval_dcsbm_like(A, zh_sc)
like[max_like_idx]
