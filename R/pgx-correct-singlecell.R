# This file is part of the Omics Playground project.
# Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
# Single-cell batch integration through Harmony, fastMNN, and BBKNN.
# Keep graph and embedding workflows outside the canonical bulk correction API.

#' @export
runHarmony <- function(X, batch) {
  library(Seurat)
  library(harmony)
  nx <- ncol(X)
  cn <- colnames(X)
  X1 <- X
  if (ncol(X) < 10) {
    X1 <- cbind(X, X, X)
    batch <- rep(batch, 3)
  }
  colnames(X1) <- paste0("col", 1:ncol(X1), "_", colnames(X1))
  M <- data.frame(batch = batch)
  rownames(M) <- colnames(X1)
  if (is.null(rownames(X1))) rownames(X1) <- paste0("row", 1:nrow(X1))
  obj <- CreateSeuratObject(2**X1, meta.data = M)
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj, verbose = FALSE)
  sX <- t(scale(t(X1)))
  hvg <- VariableFeatures(obj)
  sX <- sX[hvg, ]
  obj[["RNA"]]@layers[["scale.data"]] <- sX
  npcs <- min(30L, ncol(X1) / 2)
  npcs
  nn <- npcs
  obj <- RunPCA(obj, npcs = npcs, verbose = FALSE)
  obj <- RunUMAP(obj, reduction = "pca", dims = 1:npcs, n.neighbors = nn)
  pos0 <- obj@reductions[["umap"]]@cell.embeddings
  ## DimPlot(obj, reduction="umap", group.by='batch',pt.size=5)
  hm.obj <- harmony::RunHarmony(obj, "batch", verbose = TRUE, sigma = 0.1)
  hm.obj <- RunUMAP(hm.obj, reduction = "harmony", dims = 1:npcs, n.neighbors = nn)
  ## DimPlot(hm.obj, reduction="umap", group.by="batch", pt.size=4)
  pos1 <- hm.obj@reductions[["umap"]]@cell.embeddings
  ## reconstruct corrected matrix
  u2 <- Loadings(hm.obj, "harmony")
  u2[is.nan(u2) | is.na(u2)] <- 0 ## sometimes...
  u2 <- t(t(u2) / sqrt(1e-4 + colSums(u2**2)))
  v2 <- Embeddings(hm.obj, "harmony")
  X2 <- as.matrix(u2 %*% t(v2))

  X2 <- X2[, 1:nx]
  pos0 <- pos0[1:nx, ]
  pos1 <- pos1[1:nx, ]
  colnames(X2) <- rownames(pos0) <- rownames(pos0) <- cn
  list(corrected = X2, umap = pos1, umap.orig = pos0)
}

#' @export
MNNcorrect <- function(X, batch, controls = NULL) {
  if (is.null(controls)) {
    hvg <- head(rownames(X)[order(-apply(X, 1, sd, na.rm = TRUE))], 1000)
    res <- batchelor::fastMNN(X, batch = batch, subset.row = hvg, correct.all = TRUE)
  } else {
    xx <- tapply(1:ncol(X), batch, function(i) X[, i, drop = FALSE])
    rr <- tapply(1:ncol(X), batch, function(i) which(controls[i]))
    res <- batchelor::fastMNN(xx, restrict = rr)
  }
  cX <- as.matrix(SummarizedExperiment::assay(res))
  cX
}

##' @name bbknn State Matrix
##' @rdname bbknn
##'
##' @title Run bbknn clustering algorithm
##'
##' @description Implements the bbknn clustering algorithm in R using
##'   reticulate to run the Python version. Requires the python
##'   "bbknn" and "igraph" modules to be installed. Returns a vector
##'   of partition indices. From: https://rdrr.io/github/TomKellyGenetics/bbknn/src/R/bbknn.R
##'
##' @param data_matrix A matrix (genes x samples or cells) for expression data
##' @param batch An integer vector of batches to correct for (converts factors or numeric vectors)
##' @param pca whether to compute pca (defaults to TRUE) or apply correction to the raw matrix (FALSE)
##' @param compute_pca whether to compute PCA in Python (defaults to TRUE, requires scanpy library) or with R functions (FALSE)
##' @param nPcs number of principal components to compute (defaults to 50 if more than 50 genes)
##'
##' @return returns a list with the following components
##'   \item{corrected matrix}{matrix of data corrected by the BBKNN
##'   (batch based K nearest neighbours)}\item{pca}{principal
##'   components(matrix with row for every sample and column for each
##'   component)}\item{tsne}{t-distributed stochastic neighbour
##'   embedding (matrix with row for every sample)}\item{umap}{uniform
##'   manifold approximation and projection (matrix with row for every
##'   sample)}
##'
##' @keywords graph network igraph mvtnorm simulation
##' @import reticulate
##' @importFrom stats prcomp
##' @export
bbknn <- function(data_matrix, batch, pca = TRUE, compute_pca = "python", nPcs = NULL) {
  ## reticulate::py_install("anndata")
  ## reticulate::py_install("bbknn")

  # import python modules with reticulate
  if (!is.matrix(data_matrix)) {
    warning("matrix expected for data_matrix")
    data_matrix <- as.matrix(data_matrix)
  }
  if (is.null(nPcs)) {
    nPcs <- min(50, nrow(data_matrix), ncol(data_matrix))
  }
  if (nPcs > nrow(data_matrix)) {
    warning("number of genes less than nPcs")
    print(paste("using", nrow(data_matrix), "components"))
    ## nPcs <- nrow(data_matrix)
    nPcs <- min(nPcs, dim(data_matrix) - 1)
  }
  # reticulate::use_python("/usr/local/bin/python3")
  ##  reticulate::py_install("anndata")
  ##  reticulate::py_install("bbknn")
  ##  reticulate::py_install("scanpy")

  anndata <- reticulate::import("anndata", convert = FALSE)
  bbknn <- reticulate::import("bbknn", convert = FALSE)
  ##  sc <- reticulate::import("scanpy.api",convert=FALSE)
  sc <- reticulate::import("scanpy", convert = FALSE)

  # set up annotation data for batches
  if (is.character(batch)) batch <- as.factor(batch)
  if (is.factor(batch)) batch <- as.numeric(batch)
  if (is.numeric(batch)) batch <- as.integer(batch)

  # perform PCA
  if (pca) {
    # sc$tl$pca(adata)
    if (compute_pca == "python") {
      # use PCA computed in Python
      pca <- sc$pp$pca(t(data_matrix))
    } else if (compute_pca != "python") {
      # use PCA computed in R
      print("test")
      pca <- reticulate::r_to_py(t(prcomp(data_matrix)$x[1:nPcs, ]))
    }
    adata <- anndata$AnnData(X = pca, obs = batch)
    sc$tl$pca(adata)
    adata$obsm$X_pca <- pca
  } else {
    # use full matrix
    adata <- anndata$AnnData(X = t(data_matrix), obs = batch)
    sc$tl$pca(adata)
  }
  # perform BBKNN to derive corrected components
  bbknn$bbknn(adata, batch_key = 0)
  corrected_matrix <- t(reticulate::py_to_r(adata$X))
  sc$tl$pca(adata)
  pca <- reticulate::py_to_r(adata$obsm["X_pca"])
  sc$tl$tsne(adata)
  tsne <- reticulate::py_to_r(adata$obsm["X_tsne"])
  sc$tl$umap(adata)
  umap <- reticulate::py_to_r(adata$obsm["X_umap"])
  output <- list(corrected_matrix = corrected_matrix, pca = pca, tsne = tsne, umap = umap)
  return(output)
}
