# casc
#' Covariate Assisted Spectral Clustering
#' 
#' @param adjMat An adjacency matrix
#' @param covMat A covariate matrix
#' @param nBlocks The number of clusters
#' @param nIter Number of iterations to find  the optimal tuning
#' parameter.
#' @param method The form of the adjacency matrix to be used.
#' @param rowNorm True if row normalization should be
#' done before running kmeans.
#' @param enhancedTuning If true, then the enhanced tuning procedure is used.
#' @param center A boolean indicating if the covariate matrix columns
#' should be centered.
#' 
#'
#' @export
#' @return A list with node cluster assignments, the
#' the value of the tuning parameter used, the within
#' cluster sum of squares, and the eigengap.
#'
#' @keywords spectral clustering
casc <- function(adjMat, covMat, nBlocks, nIter = 30,
                 method = "regLaplacian", rowNorm = F,
                 enhancedTuning = F, center = F) {

    # Matrix has Namespace problems when using dsCMatrix
    adjMat = as(adjMat, "dgCMatrix")
    
    adjMat <- getGraphMatrix(adjMat, method)
    covMat <- scale(covMat, center = center,
                    scale = sqrt(Matrix::colSums(covMat^2)))

    return( getCascClusters(adjMat, covMat, nBlocks, nIter,
                            rowNorm, enhancedTuning) )    
}

# ---------------------------------------------------------------------
# Helper methods
# ---------------------------------------------------------------------

# ---------------------------------------------------------------------
# returns the graph matrix corresponding to the given method
# ---------------------------------------------------------------------
getGraphMatrix = function(adjacencyMat, method) {

    if(method == "regLaplacian") {
        rSums = Matrix::rowSums(adjacencyMat)
        tau = mean(rSums)
        normMat = Diagonal(length(rSums), 1/sqrt(rSums + tau))
        return(normMat %*% adjacencyMat %*% normMat)
    }
    else if(method == "laplacian") {
        rSums = Matrix::rowSums(adjacencyMat)
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
    nPoints, rowNorm, enhancedTuning) {

    # value for detecting a transition
    epsilon = .05
    
    rangehTuning = getTuningRange(graphMat, covariates, nBlocks)

    hTuningSeq = seq(rangehTuning$hmin, rangehTuning$hmax,
        length.out = nPoints)
    wcssVec = vector(length = nPoints)
    gapVec = vector(length = nPoints)
    orthoX = vector(length = nPoints)
    orthoL = vector(length = nPoints)
    
    for(i in 1:nPoints) {
        cascResults = getCascResults(graphMat, covariates, hTuningSeq[i],
            nBlocks, rowNorm)
        orthoX[i] = cascResults$orthoX
        orthoL[i] = cascResults$orthoL
        wcssVec[i] = cascResults$wcss
        gapVec[i] = cascResults$singGap
    }

# get transition points of static eigenvectors
    subspaces = getSubspaces(orthoX, orthoL, nPoints, epsilon)
    nSubspaces = length(subspaces$subintervalStart)    

    if((enhancedTuning == T) & (nSubspaces > 1)) {

        subMinIndex = vector(length = nSubspaces)
        subMaxIndex = vector(length = nSubspaces)
        for(i in 1:nSubspaces) {
             subMinIndex[i] = which.min(wcssVec[
                            subspaces$subintervalStart[i]:
                                subspaces$subintervalEnd[i]]) +
                                    subspaces$subintervalStart[i] - 1
            subMaxIndex[i] = which.max(wcssVec[
                           subspaces$subintervalStart[i]:
                               subspaces$subintervalEnd[i]]) +
                                   subspaces$subintervalStart[i] - 1
        }

        # keep only those intervals that are not dominated in terms of wcss
         includeVec = (rowSums(outer(wcssVec[subMinIndex], wcssVec[subMaxIndex],
                       function(x, y) {x > y})) == 0)
        
        minCountSubspaces = ((1:nSubspaces)[includeVec == 1])[
                             which.min(subspaces$orthoCounts[includeVec == 1])]

        # min WCSS on most overlapping set of subspaces
        startIndex = subspaces$subintervalStart[minCountSubspaces]
        endIndex = subspaces$subintervalEnd[minCountSubspaces]
        minInterval = unlist(apply(cbind(startIndex, endIndex), 1, function(x)
            {x[1]:x[2]}))
        minWcssSubindex = which.min(wcssVec[minInterval])
        hOpt = (hTuningSeq[minInterval])[minWcssSubindex]
    } else {
        hOpt = hTuningSeq[which.min(wcssVec)]
    }
    
    hOptResults = getCascResults(graphMat, covariates, hOpt, nBlocks, rowNorm)
    
    return( list(cluster = hOptResults$cluster,
                 h = hOpt,
                 wcss = hOptResults$wcss,
                 eigenGap = hOptResults$eigenGap) )
}

# ---------------------------------------------------------------------
# returns cluster memberships for CASC based clustering takes graphMat
# ---------------------------------------------------------------------
getCascResults = function(graphMat, covariates, hTuningParam,
    nBlocks, rowNorm) {

    randStarts = 10 #number of random starts for kmeans
    
    cascSvd = getCascSvd(graphMat, covariates, hTuningParam, nBlocks)

    ortho = getOrtho(graphMat, covariates, cascSvd$singVec, cascSvd$singVal,
        hTuningParam, nBlocks)

    if(rowNorm == T) {
        cascSvd$singVec = cascSvd$singVec/sqrt(colSums(cascSvd$singVec^2))
    }
    
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

    nCov = dim(covariates)[2]
    
    singValGraph = irlba(graphMatrix, nu = nBlocks + 1, nv = 0, m_b =
        internalDim)$d
    singValCov = svd(covariates, nu = min(nBlocks, nCov))$d

    hmax = singValGraph[1]/singValCov[min(nBlocks, nCov)]^2

    hmin = (singValGraph[nBlocks] - singValGraph[nBlocks + 1])/singValCov[1]^2

    return( list( hmax = hmax, hmin = hmin ) )
}

# ---------------------------------------------------------------------
# Finds leading subspace discontinuities.
# Returns the start and end of a continuous interval and
# the number of orthogonal components in the leading subspace
# on the interval.
# ---------------------------------------------------------------------
getSubspaces = function(orthoX, orthoL, nPoints, epsilon) {

    indicatorOut = vector(length = nPoints)
    indicatorIn = vector(length = nPoints)
    
    for(i in 1:(nPoints - 1)) {
        if((orthoX[i] < epsilon) & (orthoX[i+1] > epsilon)) {
            indicatorOut[i+1] = 1
        }
        else if((orthoL[i+1] < epsilon) & (orthoL[i] > epsilon)) {
            indicatorIn[i+1] = 1
        }
    }

    orthoCounts = cumsum(indicatorIn) - cumsum(indicatorOut) +
        max(cumsum(indicatorOut))
    subintervalStart = unique(c(which(indicatorIn == 1),
        which(indicatorOut == 1)))
    subintervalEnd = sort(c(subintervalStart-1, nPoints))
    subintervalStart = sort(c(1, subintervalStart))
    orthoCounts = orthoCounts[subintervalStart]

    return( list(orthoCounts = orthoCounts,
                 subintervalStart = subintervalStart,
                 subintervalEnd = subintervalEnd) )
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
