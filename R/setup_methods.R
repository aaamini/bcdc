if (!exists("scaled")) {
  scaled = FALSE
}
if (!exists("scx")) {
  scx = FALSE
}

methods <- list()

methods[["BSBM"]] <- function(A, Xc, Xd, K) {
  model <- new(SBM, A, K, alpha = 1, beta = 1)
  get_map_labels(model$run_gibbs(n_iter))$z_map
}

if (scaled) {
  methods[["k-means"]] <- function(A, Xc, Xd, K) {
    kmeans(scale(cbind(Xc, Xd)), centers = K, nstart = 20)$cluster # added scale
  }
} else {
  methods[["k-means"]] <- function(A, Xc, Xd, K) {
    kmeans(cbind(Xc, Xd), centers = K, nstart = 20)$cluster
  }
}

if (scx) {
  methods[["SC-X"]] <- function(A, Xc, Xd, K) {
    X = cbind(Xc, Xd)
    # Perform truncated SVD on X
    svd_res = svd(X, nu = K)
    # Use the first K singular vectors as the new X
    # X = svd_res$u[, 1:K]
    X = svd_res$u
    
    kmeans(X, centers = K, nstart = 20)$cluster
  }
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

if (scaled) {
  methods[["CASC"]] <- function(A, Xc, Xd, K) {
    kmeans(getCascAutoSvd(A, scale(cbind(Xc, Xd)), K, enhancedTuning = TRUE)$singVec
           , centers = K, nstart = 20)$cluster
  }
} else {
  methods[["CASC"]] <- function(A, Xc, Xd, K) {
    kmeans(getCascAutoSvd(A, cbind(Xc, Xd), K, enhancedTuning = TRUE)$singVec
           , centers = K, nstart = 20)$cluster
  }
}
  
  
methods[["CASCORE"]] <- function(A, Xc, Xd, K) {
  CASCORE(as.matrix(A), cbind(Xc, Xd+1), K)
}

mtd_names <- names(methods)