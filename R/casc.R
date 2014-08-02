# casc
#' Covariate Assisted Spectral Clustering
#' 
#' @param adjMat An adjacency matrix
#' @param covMat A covariate matrix
#' @param nBlocks The number of clusters
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
casc <- function(adjMat, covMat, nBlocks, nIter = 10,
                 method = "regLaplacian", rowNorm = TRUE) {
    
    adjMat <- getGraphMatrix(adjMat, method)
    covMat <- scale(covMat)

    return( getCascClusters(adjMat, covMat, nBlocks, nIter, rowNorm) )    
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
# returns CASC optimal h tuning parameter SVD
# ---------------------------------------------------------------------
getCascClusters = function(graphMat, covariates, nBlocks,
    nIter, rowNorm) {

    rangehTuning = getTuningRange(graphMat, covariates, nBlocks)

    hTuningSeq = seq(rangehTuning$hmin, rangehTuning$hmax,
        length.out = nIter)

    clusterMat = matrix(0, ncol=dim(graphMat)[1], nrow = nIter)
    wcssVec = vector(length = nIter)
    gapVec = vector(length = nIter)
    
    for(i in 1:nIter) {
        cascResults = getCascResults(graphMat, covariates, hTuningSeq[i],
            nBlocks, rowNorm)
        clusterMat[i, ] = cascResults$cluster
        wcssVec[i] = cascResults$wcss
        gapVec[i] = cascResults$eGap
    }

    # restrict possible values to those past any phase transition
    if(min(gapVec) < .9*min(gapVec[1], gapVec[nPoints]) &
       min(gapVec) < .02) {
        starth = match(min(gapVec), gapVec) + 1
        warning("A potential phase transition found. Restricting range
                 of tuning parameters.")
    }
    else {
        starth = 1
    }
        
    minWcssI = match(min(wcssVec[starth:nPoints]),
        wcssVec[starth:nPoints]) + starth - 1

    return(list(clusters = clusterMat[minWcssI, ],
                hOpt = hTuningSeq[minWcssI],
                hSeq = hTuningSeq,
                wcss = wcssVec,
                eGap = gapVec))
}

# ---------------------------------------------------------------------
# returns cluster memberships for CASC based clustering takes graphMat
# ---------------------------------------------------------------------
getCascResults = function(graphMat, covariates, hTuningParam,
    nBlocks, rowNorm) {

    randStarts = 10 #number of random starts for kmeans
    
    cascSvd = getCascSvd(graphMat, covariates, hTuningParam, nBlocks)

    if(rowNorm == TRUE) {
        cascSingVec = cascSvd$singVec/sqrt(rowSums(cascSvd$singVec^2))
    } else {
        cascSingVec = cascSvd$singVec
    }    
    
    kmeansResults = kmeans(cascSingVec, nBlocks, nstart = randStarts)
    
    return( list(cluster = kmeansResults$cluster,
                 wcss = kmeansResults$tot.withinss,
                 eGap = cascSvd$singVal[nBlocks] -
                 cascSvd$singVal[nBlocks + 1]) )
    
}

# ---------------------------------------------------------------------
# returns left singular vectors and values for CASC based clustering
# ---------------------------------------------------------------------
getCascSvd = function(graphMat, covariates, hTuningParam, nBlocks) {

    #ensure irlba internal representation is large enough
    internalDim = max(2*nBlocks, 20)

    #define a custom matrix vector multiply function
    matrixMulti = function(aList, aVector, transposeBool) {
        return( as.vector(aList$graphMat %*% aVector +
                          aList$hTuningParam * aList$covariates %*%
                          crossprod(aList$covariates, aVector) ))
    } 

    singDecomp = irlbaMod(list(graphMat = graphMat,
        covariates = covariates,
        hTuningParam = hTuningParam), nu = nBlocks + 1, nv = 0,
        m_b = internalDim, matmul = matrixMulti)

    return( list(singVec = singDecomp$u[, 1:nBlocks],
                 singVal = singDecomp$d) ) 
}


# ---------------------------------------------------------------------
# gets a good range for the tuning parameter in CASC
# ---------------------------------------------------------------------
getTuningRange = function(graphMatrix, covariates, nBlocks) {
    
    #ensure irlba internal representation is large enough
    internalDim = max(2*nBlocks, 20)
    
    singValGraph = irlba(graphMatrix, nu = nBlocks + 1, nv = 0, m_b =
        internalDim)$d
    singValCov = svd(covariates, nu = nBlocks)$d

    hmax = singValGraph[1]/singValCov[nBlocks]^2

    hmin = (singValGraph[nBlocks] - singValGraph[nBlocks + 1])/
        singValCov[1]^2


    return( list( hmax = hmax, hmin = hmin ) )
}
