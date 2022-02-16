###--------------------------------------------------------------------
# functions to run various spectral clustering methods
###--------------------------------------------------------------------


# ---------------------------------------------------------------------
# load any necessary libraries and scripts
# ---------------------------------------------------------------------
if (!require(Matrix)) {
  install.packages('Matrix', dependencies = T)
  require(Matrix)
}

if (!require(irlba)) {
  install.packages('irlba', dependencies = T)
  require(irlba)
}

source("./R/CASC/irlbaMod.R")

# ---------------------------------------------------------------------
# MAIN FUNCTIONS
# ---------------------------------------------------------------------

# ---------------------------------------------------------------------
# returns cluster memberships for CCA based clustering
# ---------------------------------------------------------------------
getCcaClusters = function(adjacencyMat, covariates, nBlocks,
                          method = "regLaplacian") {
  
  randStarts = 10 #number of random starts for kmeans
  
  gilSingVec = getGilSvd(getGraphMatrix(adjacencyMat, method), covariates,
                         nBlocks)$singVec
  
  return( kmeans(gilSingVec, nBlocks, nstart = randStarts)$cluster )
}

# ---------------------------------------------------------------------
# returns cluster memberships for CASC based clustering
# ---------------------------------------------------------------------
getCascClusters = function(adjacencyMat, covariates, hTuningParam,
                           nBlocks, method = "regLaplacian") {
  
  randStarts = 10 #number of random starts for kmeans
  
  cascSingVec = getCascSvd(getGraphMatrix(adjacencyMat, method), covariates,
                           hTuningParam, nBlocks)$singVec
  
  return( kmeans(cascSingVec, nBlocks, nstart = randStarts)$cluster )
  
}

# ---------------------------------------------------------------------
# returns cluster memberships for CASC based clustering takes graphMat
# ---------------------------------------------------------------------
getCascResults = function(graphMat, covariates, hTuningParam,
                          nBlocks) {
  
  randStarts = 10 #number of random starts for kmeans
  
  cascSvd = getCascSvd(graphMat, covariates, hTuningParam, nBlocks)
  
  ortho = getOrtho(graphMat, covariates, cascSvd$singVec, cascSvd$singVal,
                   hTuningParam, nBlocks)
  
  kmeansResults = kmeans(cascSvd$singVec, nBlocks, nstart = randStarts)
  
  return( list(cluster = kmeansResults$cluster,
               wcss = kmeansResults$tot.withinss,
               singGap = cascSvd$singVal[nBlocks] -
                 cascSvd$singVal[nBlocks + 1],
               orthoL = ortho$orthoL,
               orthoX = ortho$orthoX,
               singVecK = cascSvd$singVec[, nBlocks],
               singVecKPlus = cascSvd$singVecKPlus) )    
}

# ---------------------------------------------------------------------
# returns cluster memberships for SC based graph clustering
# ---------------------------------------------------------------------
getGraphScClusters = function(adjacencyMat, nBlocks,
                              method = "regLaplacian") {
  
  randStarts = 10 #number of random starts for kmeans
  
  scSingVec = getGraphScSvd(getGraphMatrix(adjacencyMat, method),
                            nBlocks)$singVec
  
  return( kmeans(scSingVec, nBlocks, nstart = randStarts)$cluster )
}

# ---------------------------------------------------------------------
# returns cluster memberships for SC based covariate clustering
# ---------------------------------------------------------------------
getCovScClusters = function(covariates, nBlocks) {
  
  randStarts = 10 #number of random starts for kmeans
  
  scSingVec = getCovScSvd(covariates, nBlocks)$singVec
  
  return( kmeans(scSingVec, nBlocks, nstart = randStarts)$cluster )
}


# ---------------------------------------------------------------------
# HELPER FUNCTIONS
# ---------------------------------------------------------------------

# ---------------------------------------------------------------------
# returns left singular vectors and values for GIL based clustering
# ---------------------------------------------------------------------
getGilSvd = function(graphMat, covariates, nBlocks) {
  
  #insure irlba internal representation is large enough
  if(nBlocks > 10) {
    internalDim = 2 * nBlocks
  }
  else {
    internalDim = 20
  }
  
  #approximate the generalized inverse of L using top eigenvectors
  svdL = irlba(graphMat, nu = 2*nBlocks, m_b = internalDim)
  graphMatP = svdL$u %*% Diagonal(1/svdL$d) %*% svdL$u^T
  
  #define a custom matrix vector multiply function
  matrixMulti = function(aList, aVector, transposeBool) {
    return( as.vector(aList$graphMat %*%
                        (aList$covariates %*% (aList$covariates^T
                                               %*% aVector))) )
  } 
  
  singDecomp = irlbaMod(list(graphMat = graphMatP, covariates = covariates,
                             hTuningParam = hTuningParam), nu = nBlocks + 1, nv = 0,
                        m_b = internalDim, matmul = matrixMulti) 
  
  return( list(singVec = singDecomp$u[, 1:nBlocks],
               singVal = singDecomp$d) ) 
}

# ---------------------------------------------------------------------
# returns left singular vectors and values for CCA based clustering
# ---------------------------------------------------------------------
getCcaSvd = function(graphMat, covariates, nBlocks) {
  
  singDecomp = svd(graphMat %*% covariates, nu = nBlocks, nv = 0)
  
  return( list(singVec = singDecomp$u[, 1:nBlocks],
               singVal = singDecomp$d) ) 
}

# ---------------------------------------------------------------------
# returns left singular vectors and values for CASC based clustering
# ---------------------------------------------------------------------
getCascSvd = function(graphMat, covariates, hTuningParam, nBlocks) {
  
  #insure irlba internal representation is large enough
  internalDim = max(2*nBlocks, 20)
  
  #define a custom matrix vector multiply function
  matrixMulti = function(aList, aVector, transposeBool) {
    return( as.vector(aList$graphMat %*% aVector +
                        aList$hTuningParam * aList$covariates %*%
                        crossprod(aList$covariates, aVector) ))
  } 
  
  singDecomp = irlbaMod(list(graphMat = graphMat, covariates = covariates,
                             hTuningParam = hTuningParam), nu = nBlocks + 1, nv = 0,
                        m_b = internalDim, matmul = matrixMulti)
  
  return( list(singVec = singDecomp$u[, 1:nBlocks],
               singVal = singDecomp$d,
               singVecKPlus = singDecomp$u[, nBlocks+1]) ) 
}

# ---------------------------------------------------------------------
# returns left singular vectors and values for graph SC based clustering
# ---------------------------------------------------------------------
getGraphScSvd = function(graphMat, nBlocks) {
  
  #insure irlba internal representation is large enough
  if(nBlocks > 10) {
    internalDim = 2 * nBlocks
  }
  else {
    internalDim = 20
  }
  
  singDecomp = irlba(graphMat, nu = nBlocks, nv = 0, m_b = internalDim)
  
  return( list(singVec = singDecomp$u[, 1:nBlocks],
               singVal = singDecomp$d) ) 
}

# ---------------------------------------------------------------------
# returns left singular vectors and values for covariate SC based clustering
# ---------------------------------------------------------------------
getCovScSvd = function(covMat, nBlocks) {
  
  singDecomp = svd(covMat, nu = nBlocks, nv = 0)
  
  return( list(singVec = singDecomp$u[, 1:nBlocks],
               singVal = singDecomp$d) ) 
}

# ---------------------------------------------------------------------
# returns the proportion of the eigenvalues due to X in the top eigenspace
# ---------------------------------------------------------------------
getOrtho <- function(graphMat, covariates, cascSvdSingVec, cascSvdSingVal,
                     h, nBlocks) {
  orthoL <- as.numeric((t(cascSvdSingVec[, nBlocks])%*%graphMat%*%
                          cascSvdSingVec[, nBlocks])/cascSvdSingVal[nBlocks])
  orthoX <- as.numeric(h*(t(cascSvdSingVec[, nBlocks])%*%covariates%*%
                            t(covariates)%*%cascSvdSingVec[, nBlocks])/
                         cascSvdSingVal[nBlocks])
  return( list(orthoL = orthoL/(orthoL + orthoX),
               orthoX = orthoX/(orthoL + orthoX)) )
}

# ---------------------------------------------------------------------
# returns the graph matrix corresponding to the given method
# ---------------------------------------------------------------------
getGraphMatrix = function(adjacencyMat, method) {
  
  if(method == "regLaplacian") {
    rSums = rowSums(adjacencyMat)
    tau = mean(rSums)
    normMat = Diagonal(length(rSums), 1/sqrt(rSums + tau))
    return(forceSymmetric(normMat %*% adjacencyMat %*% normMat))
  }
  else if(method == "laplacian") {
    rSums = rowSums(adjacencyMat)
    normMat = Diagonal(length(rSums), 1/sqrt(rSums))
    return(forceSymmetric(normMat %*% adjacencyMat %*% normMat))
  }
  else if(method == "adjacency"){
    return(adjacencyMat)
  }
  else {
    stop("Method given not valid.")
  }
  
  return(-1)
}