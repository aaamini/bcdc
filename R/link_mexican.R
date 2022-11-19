# Libraries ----
library(bcdc)
library(ggplot2)
library(igraph)
library(dplyr)
library(R.matlab)
library(pROC)
library(RSpectra)
source("R/inference.R")
library(Matrix)
library(network)

set.seed(123)

weddell <- read.csv("./Data/weddell_net.csv", row.names = 1)
full <- read.csv("./Data/weddell_cov.csv", sep = ";")

ind <- match(rownames(weddell), full$Species)
ind.null <- which(full[ind, "FeedingType"] == "NULL")
ind <- ind[-ind.null]

G <- as.matrix(weddell[-ind.null, -ind.null])
A <- as(ifelse(tcrossprod(G) >= 5, 1, 0), "dgCMatrix")
diag(A) <- 0
n <- nrow(A)

Cov_unscaled <- log(as.numeric(full[ind, "BodyWeight"])) # for plot
Cov <- scale(Cov_unscaled, center = TRUE, scale = TRUE)

train_ratio = 0.1
M <- matrix(0, n, n)
no_pairs <- choose(n, 2)
no_observe <- as.integer(train_ratio * no_pairs)
no_latent <- no_pairs - no_observe
M[upper.tri(M)] <- sample(c(rep(1, no_observe)
                            , rep(0, no_latent)))
M[lower.tri(M)] <- t(M)[lower.tri(M)]
M <- Matrix(M, sparse = TRUE)

test_idx <- which(triu(!M,1))

# w/ covariates
mod = new(CovarPSBM, A, M, alpha = 1, beta = 1, dp_concent = 1)
mod$set_continuous_features(Cov)
n_iter = 200
zmat = mod$run_gibbs(n_iter)
plot(get_seq_nmi(zmat))

zh = get_map_labels(zmat)$z_map

(class_freqs = table(zh))
classes = as.integer(names(class_freqs))

Zh = nett::label_vec2mat(zh, nrow(mod$eta))
EA = Zh %*% mod$eta %*% t(Zh)

roc(A[test_idx], EA[test_idx], quiet = T)

# w/out covariates
mod = new(CovarPSBM, A, M, alpha = 1, beta = 1, dp_concent = 1)
zmat = mod$run_gibbs(n_iter)
plot(get_seq_nmi(zmat))
zh = get_map_labels(zmat)$z_map
(class_freqs = table(zh))
classes = as.integer(names(class_freqs))
Zh = nett::label_vec2mat(zh, nrow(mod$eta))
EA = Zh %*% mod$eta %*% t(Zh)
roc(A[test_idx], EA[test_idx], quiet = T)