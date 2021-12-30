library(igraph)
library(nett)
library(tidyverse)
library(Matrix)
library(tictoc)
source("R/inference.R")

Rcpp::sourceCpp("src/CovarSBM.cpp", verbose = T)

# df = read_delim("data/food_web/full_data.csv", delim = ";")
df = data.table::fread("data/food_web/full_data.csv", sep = ";")
df = tibble(df)


# df$BodyWeight
# df$Mobility %>% as.integer()
df = df %>% 
  mutate(across(c(BodyWeight, TrophicLevel), as.double)) %>% 
  mutate(across(Mobility, as.integer)) %>% 
  filter(Network == "weddell")

table(df$FeedingType)
df %>% filter(FeedingType != "NULL")

df = df %>% 
  mutate(across(c(Group, Kingdom, Phylum, Class, Order, Family, Genus, FeedingMode, FeedingType, MetabolicCategory), 
                  as.factor)) 

df2 = data.table::fread("data/food_web/weddell_spnames.csv", sep = ",") %>% as_tibble()
# cbind(df2[,1], colnames(df2[,-1]))ll
A = Matrix(df2 %>% column_to_rownames(var = "V1") %>% as.matrix())


g = graph_from_adjacency_matrix(A, "undirected") #  %>% extract_largest_cc()
g2 = g

# g = graph_from_adjacency_matrix(A) #  %>% extract_largest_cc()
# A = as_adj(g)
# M = as(A %*% t(A), "dgTMatrix")
# idx = M@x >= 5
# M@x = rep(1, sum(idx))# M@x[idx]
# M@i = M@i[idx]
# M@j = M@j[idx]
# g2 = graph_from_adjacency_matrix(M)

plot_net(g2)
plot_deg_dist(g2, logx = F)
# sdegree(g2) %>% hist(50)
A = as_adj(g2)
diag(A) = 0 # important!

# reorder df according to the order of vertices of the graph
idx = sapply(seq_along(V(g2)), function(i) match(V(g2)[i]$name, df$Species))
df = df[idx,]
# cbind(df$Species, V(g2)$name)

levels(df$FeedingType)
#df = df %>% mutate(FeedingType = fct_collapse(FeedingType, carnivorous = c("carnivorous", "carnivorous/necrovorous", "detrivorous")))
#levels(df$FeedingType)

z_tru = as.integer(df$FeedingType)
n = length(z_tru)
K_tru = length(unique(z_tru))


out = nett::plot_net(g2, extract_lcc = 1, community = z_tru, vertex_alpha = .5, niter = 2000, asp=1)
coord = out$coord
TeachingDemos::zoomplot(xlim = c(-.75,.75), ylim = c(-.75,.75))
dev.copy(png, 'tru.png', width = 1000, height=1000)
dev.off()

# nco = norm_coords(coord)
# min_co = apply(nco, 2, min)
# max_co = apply(nco, 2, max)
# axis(1)
# axis(2)

# Apply SC
zh_sc = sort_labels(spec_clust(A, K_tru))
table(zh_sc)
compute_mutual_info(z_tru, zh_sc)
NMI::NMI(cbind(1:n, z_tru), cbind(1:n, zh_sc))

# bench::mark(v1 = compute_mutual_info(z_tru, zh_sc), v2 = NMI::NMI(cbind(1:n, z_tru), cbind(1:n, zh_sc)), check = F)

#png("figs/pi_net_SC.png", width = 1000, height = 1000)
out = nett::plot_net(g2, extract_lcc = 1, community = zh_sc, vertex_alpha = .5, niter = 2000, asp=1, coord = coord)
TeachingDemos::zoomplot(xlim = c(-.75,.75), ylim = c(-.75,.75))  # coordinates are normalized to [-1,1]^2
dev.copy(png, 'sc.png', width = 1000, height=1000)
dev.off()

# X = scale(as.matrix(log(df$BodyWeight), ncol=1))
X = as.matrix(log(df$BodyWeight), ncol=1)

niter = 500
model = new(CovarSBM, A, 1, 1, 1)
model$set_continuous_features(X)
tic()
zout = model$run_gibbs(niter) # 
toc()
zh_bcdc = sort_labels(get_map_labels(zout)$z_map)
table(zh_bcdc)
nett::compute_mutual_info(z_tru, zh_bcdc)

out = nett::plot_net(g2,extract_lcc = 1, community = zh_bcdc, vertex_alpha = .5, niter = 2000, asp=1, coord=coord)
TeachingDemos::zoomplot(xlim = c(-.75,.75), ylim = c(-.75,.75))
dev.copy(png, 'bcdc.png', width = 1000, height=1000)
dev.off()


seq_nmi = data.frame(iter = 1:niter, nmi = get_seq_nmi(zout))
seq_nmi %>% 
  #  filter(rep == 3) %>% 
  ggplot(aes(iter, nmi)) + 
  geom_line() + 
  xlab("Iteration") + ylab("NMI")

