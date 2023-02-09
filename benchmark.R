# Libraries ----
library(parallel)
library(doParallel)
library(NMI)
library(tidyverse)
# devtools::install_github("aaamini/bcdc/bcdc_package")
library(bcdc)
library(CASCORE)

# Functions ----
source("./R/inference.R")
source("./R/CASC/cascTuningMethods.R")

source("./R/data_gen.R")

# source("./R/setup_methods.R")

methods <- list()

methods[["BSBM"]] <- function(A, Xc, Xd, K) {
  model <- new(SBM, A, K, alpha = 1, beta = 1)
  get_map_labels(model$run_gibbs(n_iter))$z_map
}

methods[["k-means"]] <- function(A, Xc, Xd, K) {
  kmeans(scale(cbind(Xc, Xd)), centers = K, nstart = 20)$cluster # added scale
}

methods[["SC"]] <- function(A, Xc, Xd, K) {
  nett::spec_clust(A, K)
}

methods[["BCDC"]] <- function(A, Xc, Xd, K) {
  bcdc_model <- new(CovarSBM, A, alpha = 1, beta = 1, dp_concent = 10)
  if (!is.null(Xd)) bcdc_model$set_discrete_features(Xd)
  if (!is.null(Xc)) bcdc_model$set_continuous_features(Xc)
  
  get_map_labels(bcdc_model$run_gibbs(n_iter))$z_map
}


methods[["k-means (unscaled)"]] <- function(A, Xc, Xd, K) {
  kmeans(cbind(Xc, Xd), centers = K, nstart = 20)$cluster
}

# methods[["CASC (unscaled)"]] <- function(A, Xc, Xd, K) {
#   kmeans(getCascAutoSvd(A, cbind(Xc, Xd), K, enhancedTuning = TRUE)$singVec
#          , centers = K, nstart = 20)$cluster
# }

# methods[["CASC"]] <- function(A, Xc, Xd, K) {
#   kmeans(getCascAutoSvd(A, scale(cbind(Xc, Xd)), K, enhancedTuning = TRUE)$singVec
#          , centers = K, nstart = 20)$cluster
# }
# 
# methods[["CASCORE"]] <- function(A, Xc, Xd, K) {
#   CASCORE(as.matrix(A), cbind(Xc, Xd+1), K)
# }


mtd_names <- names(methods)

n <- 800
K <- 2
n_iter <- 1000

set.seed(123)
sim <- sparse_sim_data(n)
A <- sim$A
Xc <- sim$Xc
Xd = NULL
z_true <- sim$z_true

res = do.call(rbind, lapply(seq_along(methods), function(j) {
  print(mtd_names[j])
  start_time = Sys.time()
  zh <- as.vector(methods[[j]](A, Xc, Xd, K))
  end_time = as.numeric(Sys.time() - start_time)
  data.frame(method = mtd_names[j], 
             nmi_wrapper(z_true, zh),
             time = end_time)
}))

print( knitr::kable(res, digits = 4, format = "pipe") )
