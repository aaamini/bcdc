# Libraries ----
library(network)
library(igraph)
library(NMI)
library(ggplot2)
library(scales)
# devtools::install_github("aaamini/bcdc/bcdc_package")
library(bcdc)
library(CASCORE)

# Functions ----
source("./R/inference.R")
source("./R/CASC/cascTuningMethods.R")
source("R/bic.R")
source("R/waic.R")

mytriangle <- function(coords, v = NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size <- 1/175 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  
  symbols(x = coords[, 1], y = coords[, 2], bg = vertex.color
          , stars = cbind(vertex.size, vertex.size, vertex.size)
          , add = TRUE, inches = FALSE)
}
add_shape("triangle", clip = shapes("circle")$clip, plot = mytriangle)

# Data ----
mexican <- read.paj("./Data/mexican_power.paj")

# network
A <- as(as.matrix.network(mexican$networks$mexican_power.net)
        , "dgCMatrix")
n <- nrow(A)

# covariates
Cov <- as.matrix(mexican$partitions$mexican_year.clu)
Cov <- scale(Cov, center = TRUE, scale = TRUE)
fea <- dim(Cov)[2]

# clusters
Z_true <- mexican$partitions$mexican_military.clu
K <- length(unique(Z_true))

# Visualization ----
df <- data.frame(y = mexican$partitions$mexican_year.clu
                 , z = Z_true)
df <- df[order(df$y), ]
df$x <- seq_len(n)
df$z <- factor(df$z, levels = 1:2, labels = c("Military", "Civilian"))

ggplot(df) +
  geom_point(aes(x = x, y = y, color = z, shape = z)
             , size = 2) +
  labs(x = "", y = "Years since 1990", color = "", shape = "") +
  theme(axis.title.x = element_blank()
        , axis.text.x = element_blank()
        , axis.ticks.x = element_blank()) +
  scale_shape_manual(values = ifelse(df$z == "Military", "circle", "square")) +
  theme_minimal(base_size = 15)

# BCDC ----
set.seed(575)
n_iter <- 1000

bcdc <- new(CovarSBM, A, alpha = 10, beta = 10, dp_concent = 1)
bcdc$set_continuous_features(Cov)
zout <- get_map_labels(bcdc$run_gibbs(n_iter))
Z_bcdc <- zout$z_map

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
z2out <- get_map_labels(bsbm$run_gibbs(n_iter))
Z_bsbm <- z2out$z_map

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

# WAIC
get_sbm_waic(A, as.integer(Z_true))
mean(sapply(1:501, function(i) get_sbm_waic(A, sort_labels(zout$z_mat[, i]))))
get_sbm_waic(A, Z_cascore)
get_sbm_waic(A, Z_casc)
get_sbm_waic(A, Z_kmeans)
get_sbm_waic(A, Z_SC)
mean(sapply(1:501, function(i) get_sbm_waic(A, sort_labels(z2out$z_mat[, i]))))

G <- graph_from_adjacency_matrix(A, mode = "undirected")

par(mfrow = c(1, 2))
layout <- layout_with_kk(G)
plot(G
     , edge.width = 2
     , vertex.size = 20
     , vertex.color = ifelse(Z_true == 1, hue_pal()(2)[1], hue_pal()(2)[2])
     , vertex.shape = ifelse(Z_true == 1, "circle", "square")
     , vertex.label = 1:n
     , vertex.font.color = "black"
     , layout = layout, asp = 0)

plot(G
     , edge.width = 2
     , vertex.size = 20
     , vertex.color = ifelse(Z_bcdc == 1, hue_pal()(2)[1], ifelse(Z_bcdc == 2, hue_pal()(2)[2], hue_pal()(4)[4]))
     , vertex.shape = ifelse(Z_bcdc == 1, "circle", ifelse(Z_bcdc == 2, "square", "triangle"))
     , vertex.label = 1:n
     , vertex.font.color = "black"
     , layout = layout, asp = 0)