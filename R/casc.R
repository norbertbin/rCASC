# casc
#' Covariate Assisted Spectral Clustering
#' 
#' @param adjMat An adjacency matrix
#' @param covMat A covariate matrix
#' @param K The number of clusters
#' @param method The form of the adjacency matrix to be used.
#' @param rowNorm True if row normalization should be
#' done before running kmeans.
#' @param nIter Number of iterations to find  the optimal tuning
#' parameter.
#'
#' @export
#' @return A list with node cluster assignments, the
#' the value of the tuning parameter used, the within
#' cluster sum of squares.
#'
#' @keywords spectral clustering
casc <- function(adjMat, covMat, K, method = "regLaplacian",
                 rowNorm = TRUE, nIter = 10) {
    
    adjMat <- getGraphMatrix(adjMat, method)
    covMat <- scale(covMat)
    
}

# ---------------------------------------------------------------------
# Helper methods
# ---------------------------------------------------------------------

# ---------------------------------------------------------------------
# returns the graph matrix corresponding to the given method
# ---------------------------------------------------------------------
getGraphMatrix = function(adjacencyMat, method) {

    if(method == "regLaplacian") {
        rSums = rowSums(adjacencyMat)
        tau = mean(rSums)
        normMat = Diagonal(length(rSums), 1/sqrt(rSums + tau))
        return(normMat %*% adjacencyMat %*% normMat)
    }
    else if(method == "laplacian") {
        rSums = rowSums(adjacencyMat)
        normMat = Diagonal(length(rSums), 1/sqrt(rSums))
        return(normMat %*% adjacencyMat %*% normMat)
    }
    else if(method == "adjacency"){
        return(adjacencyMat)
    }
    else {
        stop(paste("Error: method =", method, "Not valid"))
    }
}

# ---------------------------------------------------------------------
# gets a good range for the tuning parameter in CASC
# ---------------------------------------------------------------------
getTuningRange = function(graphMatrix, covariates, nBlocks) {
    
    #ensure irlba internal representation is large enough
    if(nBlocks > 10) {
        internalDim = 2 * nBlocks
    }
    else {
        internalDim = 20
    }
    
    singValGraph = irlba(graphMatrix, nu = nBlocks + 1, nv = 0, m_b =
        internalDim)$d
    singValCov = svd(covariates, nu = nBlocks)$d

    hmax = singValGraph[1]/singValCov[nBlocks]^2

    hmin = (singValGraph[nBlocks] - singValGraph[nBlocks + 1])/
        singValCov[1]^2


    return( list( hmax = hmax, hmin = hmin ) )
}
