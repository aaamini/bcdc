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
  bcdc_model <- new(CovarSBM, A, alpha = 10, beta = 1, dp_concent = 10)
  bcdc_model$set_discrete_features(Xd)
  bcdc_model$set_continuous_features(Xc)
  get_map_labels(bcdc_model$run_gibbs(n_iter))$z_map
}

methods[["CASC"]] <- function(A, Xc, Xd, K) {
  kmeans(getCascAutoSvd(A, scale(cbind(Xc, Xd)), K, enhancedTuning = TRUE)$singVec
         , centers = K, nstart = 20)$cluster
}

mtd_names <- names(methods)