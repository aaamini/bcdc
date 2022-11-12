# Libraries ----
library(bcdc)
library(ggplot2)
library(igraph)
library(dplyr)
library(R.matlab)
library(pROC)
library(RSpectra)
source("R/inference.R")

train_ratio = 0.25
out = readMat(sprintf("Data/facebook-ego-%2.2f.mat", train_ratio))
# out = readMat("../Data/nips12.mat")
# out = readMat("../Data/lazega-cowork.mat")
A = out$A
M = out$mask # training mask
Xd = out$features
test_idx <- which(triu(!M,1))
image(Xd)

cat(sprintf("Traning density = %2.2f\n", sum(M) / prod(dim(M))))

g = graph_from_adjacency_matrix(A, "undirected")
nett::plot_net(g)

image(A * M)
image(M)
image(A)


out = RSpectra::svds(Xd, 20)
plot(out$d)
Xc = out$u[,10, drop=F]


set.seed(123)
mod = new(CovarPSBM, A, M, alpha = 1, beta = 10, dp_concent = 10)
mod$s2 = .1
# mod = new(CovarSBM, A, alpha = 1, beta = 1, dp_concent = 10)
# mod$set_discrete_features(Xd)
mod$set_continuous_features(Xc)
n_iter = 500
zmat = mod$run_gibbs(n_iter)
plot(get_seq_nmi(zmat))

zh = get_map_labels(zmat)$z_map

(class_freqs = table(zh))
classes = as.integer(names(class_freqs))
# eta = mod$eta[classes, classes]
# nett::plot_net(g, community = zh)


Zh = nett::label_vec2mat(zh, nrow(mod$eta))
EA = Zh %*% mod$eta %*% t(Zh)


roc1 = roc(A[test_idx], EA[test_idx], quiet = T)


zh = sort_labels(zh)
eta = nett::estim_dcsbm(A, sort_labels(zh))$B
Zh = nett::label_vec2mat(zh)
EA = Zh %*% eta %*% t(Zh)
roc2 = roc(A[test_idx], EA[test_idx], quiet = T)

cat(sprintf("auc1 = %2.2f, auc2 = %2.2f", roc1$auc, roc2$auc))
# plot(roc)

PRROC::pr.curve(EA[test_idx], weights.class0 = A[test_idx])

# This is from SymNARM from MATLAB code
out = readMat(sprintf("Data/facebook-ego-%2.2f-result.mat", train_ratio))
roc3 = roc(A[test_idx], out$ProbAve[test_idx], quiet = T)
roc3$auc
PRROC::pr.curve( out$ProbAve[test_idx], weights.class0 = A[test_idx])
# #yardstick::pr_auc_vec(factor(A[test_idx]), out$ProbAve[test_idx])


# mod = new(BasicPSBM, A, M, 3)
# mod$z
# mod$set_z_to_random_labels()
# j = 20
# mod$col_adj_compress(j)
# mod$col_mask_compress(j)
# z = mod$z+1
# aggregate(mod$mask[j+1,], list(z), sum)
#
# mod$update_blk_sums_and_sizes()
# temp = nett::compute_block_sums(mod$A, z)
# diag(temp) = diag(temp)/2
# temp
# mod$M
#
# temp = nett::compute_block_sums(mod$mask, z)
# diag(temp) = diag(temp)/2
# temp
# mod$N
