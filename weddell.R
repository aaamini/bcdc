# Libraries ----
library(ggplot2)
library(NMI)
library(igraph)
library(scales)

# Functions ----
Rcpp::sourceCpp("src/SBM.cpp", verbose = T)
Rcpp::sourceCpp("src/CovarSBM.cpp", verbose = T)
source("R/inference.R")
# source("R/cascTuningMethods.R")
# Rcpp::sourceCpp("./JCDC/JCDC.cpp", verbose = T)

# Data ----
weddell <- read.csv("data/food_web/weddell_spnames.csv", row.names = 1)
full <- read.csv("data/food_web/full_data.csv", sep = ";")

ind <- match(rownames(weddell), full$Species)
ind.null <- which(full[ind, "FeedingType"] == "NULL")
ind <- ind[-ind.null]

# network
G <- as.matrix(weddell[-ind.null, -ind.null])
A <- as(ifelse(tcrossprod(G) >= 5, 1, 0), "dgCMatrix")
diag(A) <- 0
n <- nrow(A)
Matrix::image(A)

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

sort_labels(as.integer(Z_true))
table(Z_true)

# Visualization ----
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

ggsave(filename = "./weddell.pdf", width = 8, height = 5)

# Discrete covariates
da = tibble(full[ind, ])
# glimpse(da)
db = da %>% mutate(across(c(Group, Kingdom, Phylum, Class, Order, Family, Genus, FeedingMode, MetabolicCategory, TrophicLevel, Mobility, Environment), ~ as.integer(factor(.x))))
# glimpse(db)
Xd = db %>% select(Mobility, FeedingMode) %>% as.matrix() - 1
apply(Xd, 2, function(x) length(unique(x))) # number of categories for each discrete cov.

# Try uncommenting to see the Producer vs. Non-producer effect
# levels(Z_true) = c("p","np","np","np")

# BCDC ----
library(tictoc)
library(tibble)
library(tidyverse)
set.seed(575)
n_iter <- 200

bcdc <- new(CovarSBM, A, alpha = 1, beta = 1, dp_concent = 20)

bcdc$set_discrete_features(Xd)
bcdc$set_continuous_features(Cov)

zout = bcdc$run_gibbs(n_iter)
Z_bcdc <- get_map_labels(zout)$z_map

# bench::mark(v1 = bcdc$run_gibbs(n_iter))

nmi_wrapper(Z_true, Z_bcdc)
NMI(cbind(seq_len(n), Z_true), cbind(seq_len(n), Z_bcdc))
data.frame(itr=1:n_iter, nmi = get_seq_nmi(zout)) %>% 
  ggplot(aes(itr, nmi)) + 
  geom_line() + xlab("Iteration") + ylab("Seqential NMI")

Z_bcdc_sorted = sort_labels(Z_bcdc)
table(Z_bcdc_sorted)
idx = order(Z_bcdc_sorted)
Matrix::image(A[idx,idx])
Matrix::image(A)



# CASC ----
Z_casc <- kmeans(getCascSvd(graphMat = A
                            , covariates = Cov
                            , hTuningParam = mean(colSums(A))
                            , nBlocks = K)$singVec
                 , centers = K, nstart = 20)$cluster



# k-means ----
Z_kmeans <- kmeans(Cov, centers = K, nstart = 20)$cluster

# SC ----
Z_SC <- nett::spec_clust(A, K)

# BSBM ----
bsbm <- new(SBM, A, K, alpha = 1, beta = 1)
Z_bsbm <- get_map_labels(bsbm$run_gibbs(n_iter))$z_map

# JCDC ----
RUNJCDC <- function(A, phi, K) {
  
  n <- nrow(A)
  
  D.inv <- diag(1 / (sqrt(apply(A, 1, sum)) + 1e-7))
  Laplacian <- D.inv %*% A %*% D.inv
  L.svd <- svd(Laplacian)
  U.K <- L.svd$u[, 1:K]
  spec.cluster <- kmeans(U.K, K, nstart = 20)$cluster
  G.fit <- array(0, c(n, K))
  for(k in 1:K) G.fit[spec.cluster == k, k] <- 1
  
  if (length(dim(phi)) == 2) {
    p <- 1
  } else {
    p <- dim(phi)[3]    
  }
  
  W_max <- 5
  result <- JCDC(A, phi, p, G.fit, 1, K, 20, 30, 20, W_max, 2, 0, 1)
  
  return(as.vector(result$G.fit %*% seq_len(K)))
  
}

Z_jcdc <- RUNJCDC(as.matrix(A)
                  , as.matrix(dist(Cov, method = "manhattan"))
                  , K)

# Results ----
nmi_wrapper(Z_true, Z_bcdc)
nmi_wrapper(Z_true, Z_casc)
nmi_wrapper(Z_true, Z_kmeans)
nmi_wrapper(Z_true, Z_SC)
nmi_wrapper(Z_true, Z_bsbm)
nmi_wrapper(Z_true, Z_jcdc)

# NMI(cbind(seq_len(n), Z_true), cbind(seq_len(n), Z_bcdc))
# NMI(cbind(seq_len(n), Z_true), cbind(seq_len(n), Z_casc))
# NMI(cbind(seq_len(n), Z_true), cbind(seq_len(n), Z_kmeans))
# NMI(cbind(seq_len(n), Z_true), cbind(seq_len(n), Z_SC))
# NMI(cbind(seq_len(n), Z_true), cbind(seq_len(n), Z_bsbm))
# NMI(cbind(seq_len(n), Z_true), cbind(seq_len(n), Z_jcdc))

# Stability ----
library(pbmcapply)
n_reps <- 10
n_cores <- parallel::detectCores()

res <- simplify2array(
  pbmclapply(seq_len(n_reps), mc.cores = n_cores, function(i)
  {
    set.seed(i)
    bcdc <- new(CovarSBM, A, alpha = 1, beta = 1, dp_concent = 20)
    # bcdc$set_node_features(Cov)
    
    bcdc$set_discrete_features(Xd)
    bcdc$set_continuous_features(Cov)

    return(get_map_labels(bcdc$run_gibbs(n_iter))$z_map)
    
    # bsbm <- new(SBM, A, K, alpha = 1, beta = 1)
    # return(get_map_labels(bsbm$run_gibbs(n_iter))$z_map)
    
  }))

pw_NMI <- unlist(sapply(seq_len(n_reps - 1), function (i) {
  sapply((i+1):n_reps, function (j) {
    NMI(cbind(1:n, res[, i]), cbind(1:n, res[, j]))
  })
}))

mean(pw_NMI)
