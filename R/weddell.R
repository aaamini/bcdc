# Libraries ----
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(igraph)
library(NMI)
library(scales)
library(dplyr)
# devtools::install_github("aaamini/bcdc/bcdc_package")
library(bcdc)
library(CASCORE)

# Functions ----
source("./R/inference.R")
source("./R/CASC/cascTuningMethods.R")
source("R/bic.R")

# Data ----
weddell <- read.csv("./Data/weddell_net.csv", row.names = 1)
full <- read.csv("./Data/weddell_cov.csv", sep = ";")

ind <- match(rownames(weddell), full$Species)
ind.null <- which(full[ind, "FeedingType"] == "NULL")
ind <- ind[-ind.null]

# network
G <- as.matrix(weddell[-ind.null, -ind.null])
A <- as(ifelse(tcrossprod(G) >= 5, 1, 0), "dgCMatrix")
diag(A) <- 0
n <- nrow(A)

# covariates
Cov_unscaled <- log(as.numeric(full[ind, "BodyWeight"])) # for plot
Cov <- scale(Cov_unscaled, center = TRUE, scale = TRUE)
fea <- dim(Cov)[2]

# clusters
Z_true <- full[ind, "FeedingType"]
Z_true <- factor(Z_true
                 , levels = c("primaryproducer"
                              , "herbivorous/detrivorous"
                              , "detrivorous", "carnivorous", "carnivorous/necrovorous"
                              , "omnivorous")
                 , labels = c("Producer"
                              , "Herbivore"
                              , "Carnivore", "Carnivore", "Carnivore"
                              , "Omnivore"))
K <- length(unique(Z_true))

# Plot covariates ----
df <- data.frame(y = Cov_unscaled
                 , z = Z_true)
df <- df[order(df$z), ]
df$x <- seq_len(n)

ggplot(df) +
  geom_point(aes(x = x, y = y, color = z, shape = z)
             , size = 2) +
  labs(x = "", y = "Log of body weight", color = "", shape = "") +
  theme(axis.title.x = element_blank()
        , axis.text.x = element_blank()
        , axis.ticks.x = element_blank()) +
  theme_minimal(base_size = 15)

# BCDC ----
set.seed(575)
n_iter <- 200

bcdc <- new(CovarSBM, A, alpha = 1, beta = 1, dp_concent = 20)
bcdc$set_continuous_features(Cov)
Z_bcdc <- get_map_labels(bcdc$run_gibbs(n_iter))$z_map

# CASC ----
Z_casc <- kmeans(getCascAutoSvd(A, scale(Cov)
                                , K, enhancedTuning = TRUE)$singVec
                 , centers = K, nstart = 20)$cluster

# CASCORE ----
Z_cascore <- CASCORE(as.matrix(A), Cov, K)

# k-means ----
Z_kmeans <- kmeans(Cov, centers = K, nstart = 20)$cluster

# SC ----
Z_SC <- nett::spec_clust(A, K)

# BSBM ----
bsbm <- new(SBM, A, K, alpha = 1, beta = 1)
Z_bsbm <- get_map_labels(bsbm$run_gibbs(n_iter))$z_map

# Results ----

# NMI
nmi_wrapper(Z_true, Z_bcdc)
nmi_wrapper(Z_true, Z_casc)
nmi_wrapper(Z_true, Z_cascore)
nmi_wrapper(Z_true, Z_kmeans)
nmi_wrapper(Z_true, Z_SC)
nmi_wrapper(Z_true, Z_bsbm)

# BIC
get_sbm_bic(A, as.integer(Z_true))
get_sbm_bic(A, sort_labels(Z_bcdc))
get_sbm_bic(A, Z_cascore)
get_sbm_bic(A, Z_casc)
get_sbm_bic(A, Z_kmeans)
get_sbm_bic(A, Z_SC)
get_sbm_bic(A, Z_bsbm)

# Visualize network ----
Z_sort <- c()
for (type in levels(Z_true)) {
  ind <- which(Z_true == type)
  Z_sort <- c(Z_sort
              , ind[hclust(dist(A[ind, ind], method = "binary"))$order])
}

pheatmap(A[Z_sort, Z_sort]
         , color = colorRampPalette(brewer.pal(9,"Greys")[c(1,8)])(30)
         , cluster_cols = F, cluster_rows = F
         , gaps_row = cumsum(tabulate(Z_true))
         , gaps_col = cumsum(tabulate(Z_true))
         , annotation_names_row = F
         , show_rownames = F, show_colnames = F
         , legend = F)

Z_bcdc_sorted <- sort_labels(Z_bcdc)
idx <- order(Z_bcdc_sorted)

pheatmap(A[idx, idx]
         , color = colorRampPalette(brewer.pal(9,"Greys")[c(1,8)])(30)
         , cluster_cols = F, cluster_rows = F
         , gaps_row = cumsum(tabulate(Z_bcdc_sorted))
         , gaps_col = cumsum(tabulate(Z_bcdc_sorted))
         , annotation_names_row = F
         , show_rownames = F, show_colnames = F
         , annotation_legend = F, legend = F)