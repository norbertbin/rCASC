rCASC
=====

A simple R package for doing covariate assisted spectral clustering using a modified version of the **irlba** package for partial SVD computation. 

Installation
===

The **rCASC** package can be installed in R directly from GitHub by using **devtools**.

```r
library(devtools)
install_github("norbertbin/rCASC")
```

Basic Usage
===
The required input for the *casc* function includes an adjacency matrix, *adjMat*, a node covariate matrix, *covMat*, and the number of blocks to be recovered, *nBlocks*. For more details see the documentation. 

```r
casc(adjMat, covMat, nBlocks)
```