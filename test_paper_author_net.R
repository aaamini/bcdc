library(igraph)
library(tictoc)
library(ggplot2)
Rcpp::sourceCpp("src/CovarSBM.cpp", verbose = T)
source("R/inference.R")

net_name = "cs_author_paper"
load(sprintf("data/%s.RData", net_name))
g = graph_from_incidence_matrix(A)

out = bipartite_projection(g)
papers_g =out$proj2
A = as_adj(papers_g) 
# g2 = nett::extract_largest_cc(papers_g)

V(papers_g)$comm = z 
g2 = induced_subgraph(papers_g, coreness(papers_g) > 5)
out = nett::plot_net(g2, extract_lcc = 1, community = V(g2)$comm, vertex_alpha = .3, niter = 1500)
coord = out$coord

# Spectral clustering
zh = nett::spec_clust(A, 5)
nett::compute_mutual_info(zh, z)
V(papers_g)$comm = zh
g2 = induced_subgraph(papers_g, coreness(papers_g) > 5)
out = nett::plot_net(g2, extract_lcc = 1, community = V(g2)$comm, vertex_alpha = .3, niter = 1500, coord = coord)

# nett::plot_deg_dist(papers_g)
# summary(degree(papers_g))

U = RSpectra::svds(X, 20)$u

niter = 1000
model = new(CovarSBM, A, 1, 1, .5)
model$set_node_features(U) 
tic()
zout = model$run_gibbs(niter) # 
toc()
nett::compute_mutual_info(z, model$z+1)
zh = get_map_labels(zout)$z_map
nett::compute_mutual_info(z, zh)
table(model$z)
table(zh)

V(papers_g)$comm = zh
g2 = induced_subgraph(papers_g, coreness(papers_g) > 5)
out = nett::plot_net(g2, extract_lcc = 1, community = V(g2)$comm, vertex_alpha = .3, niter = 1500, coord = coord)

df = data.frame(itr=1:niter, nmi=sapply(1:niter, function(itr) nett::compute_mutual_info(z, zout[,itr])))
df %>% ggplot(aes(itr, nmi)) + 
  geom_line() +
  xlab("Iteration") + ylab("NMI")


