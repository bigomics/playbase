##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##
##


#' @title Cluster genes
#'
#' @param pgx A PGX object containing gene expression data
#' @param method Character vector of clustering methods to apply (choose between "pca", "tsne", and "umap")
#' @param dims Dimensions for dimensionality reduction (by default 2,and 3)
#' @param reduce.pca Number of PCs to reduce to for PCA
#' @param perplexity Perplexity parameter for tSNE
#' @param level Data level to cluster ("gene", "exon", "junction")
#' @param rank.tf Rank transform expression first? Logical
#' @param center.rows Center expression rows? Logical
#' @param scale.rows Scale expression rows? Logical
#' @param X Expression matrix to use instead of pgx$X
#' @param umap.pkg UMAP package to use ("uwot" or "umap")
#'
#' @return The updated PGX object with gene clusters added
#'
#' @description
#' Clusters the genes in a PGX object using dimensionality reduction methods like PCA, tSNE and UMAP.
#'
#' @details
#' This function takes a PGX object and performs dimensionality reduction on the gene expression data,
#' followed by clustering in the reduced space. Methods like PCA, tSNE and UMAP are used to reduce the dimensions.
#'
#' The gene expression matrix is extracted from the PGX object, optionally transformed, then embedded.
#' The resulting coordinates are clustered using k-means to assign genes to clusters.
#'
#' The gene clusters are added to the PGX object for downstream use.
#'
#' @export
pgx.clusterGenes <- function(pgx, methods = c("pca", "tsne", "umap"), dims = c(2, 3),
                             reduce.pca = 50, perplexity = 30, level = "gene",
                             rank.tf = FALSE, center.rows = TRUE, scale.rows = FALSE,
                             find.clusters = FALSE, X = NULL, umap.pkg = "uwot") {
  if (!is.null(X)) {
    message("using provided X matrix...")
  } else if (!is.null(pgx$X) && level == "gene") {
    message("using expression gene X matrix...")
    X <- pgx$X
  } else if (!is.null(pgx$gsetX) && level == "geneset") {
    message("using expression geneset X matrix...")
    X <- pgx$gsetX
    X <- X[complete.cases(X), , drop = FALSE]
    if (nrow(X) == 0) {
      message("WARNING:: pgx$gsetX has 0 complete cases. Returning pgx.")
      return(pgx)
    }
  } else {
    message("WARNING:: could not find matrix X")
    return(pgx)
  }

  if (center.rows) {
    X <- X - rowMeans(X, na.rm = TRUE)
  }
  if (scale.rows) {
    X <- X / (1e-6 + matrixStats::rowSds(X, na.rm = TRUE))
  }
  if (rank.tf) {
    X <- scale(apply(X, 2, rank))
  }

  ## Do dimensionality reduction
  message("[pgx.clusterGenes] computing dimensionality reductions...")
  clust <- pgx.clusterBigMatrix(
    t(X),
    methods = methods,
    dims = dims,
    perplexity = perplexity,
    center.features = FALSE,
    scale.features = FALSE,
    reduce.sd = -1,
    reduce.pca = reduce.pca,
    find.clusters = FALSE,
    umap.pkg = umap.pkg
  )

  ## Find clusters
  clust$membership <- NULL
  clust.index <- NULL
  if(find.clusters) {
    message("[pgx.clusterGenes] finding clusters using Louvain...")
    posx <- do.call(cbind, clust)
    posx <- apply(posx,2,scale)
    rownames(posx) <- rownames(X)
    dim(posx)
    clust.index <- pgx.findLouvainClusters(
      posx, graph.method="knn", knn=200, level=1,
      prefix="M", small.zero=0.01
    )
    names(clust.index) <- rownames(posx)    
  }
  
  ## remove empty space in tSNE/UMAP
  ii <- grep("tsne|umap", names(clust))
  if (length(ii) > 0) {
    clust[ii] <- lapply(clust[ii], pos.compact)
  }

  ## put in slot 'gene cluster'
  if (level == "gene") {
    pgx$cluster.genes <- NULL
    pgx$cluster.genes$pos <- clust
    pgx$cluster.genes$index <- clust.index
  }
  if (level == "geneset") {
    pgx$cluster.gsets <- NULL
    pgx$cluster.gsets$pos <- clust
    pgx$cluster.gsets$index <- clust.index
  }

  message("[pgx.clusterGenes] done!")
  pgx
}


#' Cluster samples in a PGX object
#'
#' @param pgx A PGX object
#' @param method Character vector of clustering methods to apply (choose between "pca", "tsne", and "umap")
#' @param dims Dimensions for dimensionality reduction (by default 2,and 3)
#' @param reduce.sd Perform dimensionality reduction if original dims > reduce.sd
#' @param reduce.pca For PCA, take first reduce.pca components
#' @param perplexity Perplexity parameter for tSNE and UMAP
#' @param center.rows Center rows before clustering?
#' @param scale.rows Scale rows to unit variance before clustering?
#' @param X Use this matrix instead of pgx$X
#' @param umap.pkg UMAP package to use (umap or uwot)
#' @param replace.orig Replace original pgx$X with reduced data?
#'
#' @return The updated PGX object with cluster assignments in pgx$cluster$sample
#'
#' @description
#' Performs dimensionality reduction and clustering on the samples in a PGX object.
#'
#' @details
#' This function takes a PGX object and performs PCA, tSNE, and/or UMAP on the
#' samples. It clusters the samples in the reduced space and stores the cluster
#' assignments in the pgx$cluster$sample slot.
#'
#' Dimensionality reduction is performed on pgx$X or a provided matrix X. If the
#' original number of dimensions is large, PCA is used to reduce to reduce.pca dimensions
#' first. tSNE and UMAP are then applied to further reduce to 2 or 3 dimensions.
#'
#' Multiple parameter options are available to control the clustering, including
#' choice of methods, perplexity, centering/scaling samples first, etc.
#'
#' By default the original pgx$X matrix is replaced with the reduced data. Set
#' replace.orig=FALSE to retain the original.
#'
#' @export
pgx.clusterSamples <- function(pgx, methods = c("pca", "tsne", "umap"),
                               dims = c(2, 3),
                               reduce.sd = 1000, reduce.pca = 50, perplexity = 30,
                               center.rows = TRUE, scale.rows = FALSE,
                               X = NULL, umap.pkg = "uwot", replace.orig = TRUE) {
  if (!is.null(X)) {
    message("using provided X matrix...")
  } else if (!is.null(pgx$impX)) {
    message("using pgx$impX matrix...")
    X <- pgx$impX
  } else if (!is.null(pgx$X)) {
    message("using pgx$X matrix...")
    X <- pgx$X
  } else {
    message("using logCPM(pgx$counts)...")
    X <- logCPM(pgx$counts, total = NULL)
  }

  if (any(is.na(X))) {
    ## NEED RETHINK: We should use impX here if available. Some
    ## datasets have missing values on all rows!!!
    ## X <- X[complete.cases(X), , drop = FALSE]
    X <- svdImpute2(X) ## IK
  }

  clust.pos <- pgx.clusterBigMatrix(
    X,
    methods = methods,
    dims = dims,
    perplexity = perplexity,
    center.features = center.rows,
    scale.features = scale.rows,
    reduce.sd = reduce.sd,
    reduce.pca = reduce.pca,
    find.clusters = FALSE,
    umap.pkg = umap.pkg
  )

  clust.index <- clust.pos$membership
  clust.pos$membership <- NULL

  if (replace.orig) {
    message("[pgx.clusterSamples2] update tsne2d/tsne3d and 'cluster' pheno...")
    pgx$samples$cluster <- clust.index
    pgx$tsne2d <- clust.pos[["tsne2d"]]
    pgx$tsne3d <- clust.pos[["tsne3d"]]
  } else {
    message("[pgx.clusterSamples2] skipping tsne2d/tsne3d update...")
  }

  pgx$cluster <- NULL
  pgx$cluster$pos <- clust.pos
  pgx$cluster$index <- clust.index

  message("[pgx.clusterSamples2] done!")
  pgx
}

#' @export
pgx.clusterSamples2 <- function(...) pgx.clusterSamples(...)


#' Find optimal number of sample clusters
#'
#' @title Find optimal sample clusters
#'
#' @param X Numeric matrix of expression data
#' @param method Clustering methods to evaluate (kmeans, hclust, louvain, meta)
#' @param top.sd Number of top variable features to use
#' @param npca Number of PCA components for dimensionality reduction
#'
#' @return Optimal number of clusters
#'
#' @description Finds the optimal number of clusters in the samples using multiple methods.
#'
#' @details This function takes a gene expression matrix and evaluates multiple clustering algorithms
#' (kmeans, hierarchical, Louvain) across a range of cluster numbers to determine the optimal number
#' of clusters.
#'
#' The data is first reduced to the top variable features and PCA dimensions. Clustering is then
#' performed for a range of k clusters. The optimal k is determined by the minimum average silhouette
#' width across methods.
#' @export
pgx.FindClusters <- function(X, method = c("kmeans", "hclust", "louvain", "meta"),
                             top.sd = 1000, npca = 50, scale=TRUE) {
  message("[FindClusters] called...")

  km.sizes <- c(2, 3, 4, 5, 7, 10, 15, 20, 25, 50, 100)
  km.sizes <- km.sizes[km.sizes < ncol(X)]
  km.sizes
  if (length(method) == 1 && method[1] == "meta") {
    method <- c("kmeans", "hclust", "louvain", "meta")
  }

  ## reduce dimensions
  ## X <- Matrix::head(X[order(apply(X, 1, stats::sd, na.rm = TRUE)), ], top.sd)
  X <- Matrix::head(X[order(matrixStats::rowSds(X, na.rm = TRUE)), ], top.sd)
  if(scale) {
    X <- t(scale(t(X))) ## scale features??
  }

  if (nrow(X) > npca) {
    npca <- min(npca, dim(X) - 1)
    suppressMessages(suppressWarnings(
      out <- irlba::irlba(X, nv = npca)
    ))
    X <- t(out$v)
  }

  index <- list()

  ## perform K-means
  if ("kmeans" %in% method) {
    message("perform K-means...")
    km <- lapply(km.sizes, function(k) stats::kmeans(t(X), k, iter.max = 10))
    km.idx <- do.call(cbind, lapply(km, function(r) r$cluster))
    colnames(km.idx) <- paste0("kmeans.", km.sizes)
    index[["kmeans"]] <- km.idx
  }

  ## perform hclust (on positions)
  if ("hclust" %in% method) {
    message("perform hclust...")
    hc <- fastcluster::hclust(stats::dist(t(X)), method = "ward.D")
    hc.idx <- lapply(km.sizes, function(k) stats::cutree(hc, k))
    hc.idx <- do.call(cbind, hc.idx)
    colnames(hc.idx) <- paste0("hclust.", km.sizes)
    index[["hclust"]] <- hc.idx
  }

  ## perform Louvain clustering
  if ("louvain" %in% method) {
    message("perform Louvain clustering...")
    gr <- scran::buildSNNGraph(X)
    gr.idx <- hclustGraph(gr, k = 4) ## iterative cluster until level3
    rownames(gr.idx) <- colnames(X)
    nc <- apply(gr.idx, 2, function(x) length(unique(x)))
    colnames(gr.idx) <- paste0("louvain.", nc)
    index[["louvain"]] <- gr.idx
  }

  ## find meta-index
  if ("meta" %in% method && length(index) > 1) {
    message("perform meta clustering...")
    K <- do.call(cbind, index)
    k.rows <- split(K, row(K))
    d1 <- outer(k.rows, k.rows, Vectorize(function(x, y) sum(x != y)))
    rownames(d1) <- colnames(d1) <- rownames(K)
    hc <- fastcluster::hclust(stats::as.dist(d1))
    meta.idx <- do.call(cbind, lapply(km.sizes, function(k) stats::cutree(hc, k)))
    colnames(meta.idx) <- paste0("meta.", km.sizes)
    rownames(meta.idx) <- rownames(K)
    index[["meta"]] <- meta.idx
  }

  ## sort cluster index from big to small clusters
  relevelBig2small <- function(idx) as.integer(factor(idx, levels = names(sort(-table(idx)))))
  for (i in 1:length(index)) {
    index[[i]] <- apply(index[[i]], 2, relevelBig2small)
  }

  return(index)
}


#' Cluster large matrices
#'
#' @title Cluster large matrices
#'
#' @param X Numeric matrix of data to cluster
#' @param methods Character vector of methods to use (pca, tsne, umap)
#' @param dims Output dimensions for dimensionality reduction
#' @param perplexity Perplexity parameter for tSNE and UMAP
#' @param reduce.sd Perform PCA if original dimensions > reduce.sd
#' @param reduce.pca Number of PCA components for initial reduction
#' @param center.features Center the rows before clustering?
#' @param scale.features Scale the rows before clustering?
#' @param find.clusters Detect optimal number of clusters?
#' @param umap.pkg UMAP package to use (umap or uwot)
#'
#' @return List of cluster assignments
#'
#' @description
#' Performs dimensionality reduction and clustering on large data matrices
#'
#' @details
#' This function takes a large data matrix and reduces the dimensions using
#' PCA and tSNE/UMAP before clustering the rows. It is designed to handle
#' matrices that are too large to embed directly with tSNE or UMAP.
#'
#' PCA is first used to reduce the dimensions if the original size is over
#' reduce.sd rows. tSNE or UMAP is then applied to further reduce the
#' dimensions for clustering. k-means clustering is performed on the
#' reduced data.
#'
#' Multiple parameters allow customization of the methods, perplexity,
#' dimensions, etc. By default it tries to detect the optimal number of
#' clusters using the elbow method.
#'
#' @export
pgx.clusterMatrix <- function(X,
                              ## methods = c("pca", "tsne", "umap", "pacmap"),
                              methods = c("pca", "tsne", "umap"),
                              dims = c(2, 3),
                              perplexity = 30, reduce.sd = 1000, reduce.pca = 50,
                              center.features = TRUE, scale.features = FALSE,
                              find.clusters = FALSE, umap.pkg = "uwot") {
  methods <- intersect(methods, c("pca", "tsne", "umap", "pacmap"))
  if (length(methods) == 0) methods <- "pca"

  ## Reduce dimensions by SD
  dimx <- dim(X) ## original dimensions
  namesx <- colnames(X)
  if (reduce.sd > 0 && nrow(X) > reduce.sd) {
    sdx <- matrixStats::rowSds(X, na.rm = TRUE)
    is.constant <- all(abs(sdx - mean(sdx, na.rm = TRUE)) < 1e-8, na.rm = TRUE)
    if (is.constant) {
      message("WARNING:: SD is constant. Skipping SD reduction...\n")
    } else {
      message("Reducing to ", reduce.sd, " max SD features...\n")
      X <- X[utils::head(order(-sdx), reduce.sd), ]
    }
  }

  ## scale and augment if few samples
  if (center.features) {
    X <- X - rowMeans(X, na.rm = TRUE) ## do??
  }
  if (scale.features) {
    sdx <- matrixStats::rowSds(X, na.rm = TRUE)
    X <- X / sdx
  }

  ## impute on row median
  ## if (any(is.na(X))) {
  ##    X <- imputeMedian(X)
  ## }

  if (ncol(X) <= 6) X <- cbind(X, X, X, X, X, X)
  if (nrow(X) <= 3) X <- rbind(X, X, X, X)

  ## add small variation...
  X <- X + 1e-3 * matrix(stats::rnorm(length(X)), nrow(X), ncol(X))

  ## Further pre-reduce dimensions using SVD
  res.svd <- NULL
  if (reduce.pca > 0) {
    reduce.pca <- max(3, min(c(reduce.pca, dim(X) - 1)))
    message("Reducing to ", reduce.pca, " PCA dimensions...")
    cnx <- colnames(X)
    suppressMessages(suppressWarnings(
      res.svd <- irlba::irlba(X, nv = reduce.pca)
    ))
    X <- t(res.svd$v) * res.svd$d ## really weight with D??
    colnames(X) <- cnx
  }

  all.pos <- list()

  if ("pca" %in% methods && 2 %in% dims) {
    message("calculating PCA 2D/3D...")
    if (is.null(res.svd)) {
      suppressMessages(suppressWarnings(
        res.svd <- irlba::irlba(X, nv = 3)
      ))
    }
    pos <- res.svd$v[, 1:2]
    rownames(pos) <- colnames(X)
    pos <- pos[1:dimx[2], ] ## if augmented
    colnames(pos) <- paste0("PC-", c("x", "y"))
    all.pos[["pca2d"]] <- pos
  }

  if ("pca" %in% methods && 3 %in% dims) {
    if (is.null(res.svd)) {
      suppressMessages(suppressWarnings(
        res.svd <- irlba::irlba(X, nv = 3)
      ))
    }
    pos <- res.svd$v[, 1:3]
    rownames(pos) <- colnames(X)
    pos <- pos[1:dimx[2], ] ## if augmented
    colnames(pos) <- paste0("PC-", c("x", "y", "z"))
    all.pos[["pca3d"]] <- pos
  }

  if ("tsne" %in% methods && 2 %in% dims) {
    message("calculating t-SNE 2D...")
    perplexity <- pmax(min(ncol(X) / 4, perplexity), 2)
    perplexity
    res1 <- Rtsne::Rtsne(t(X),
      dims = 2,
      is_distance = FALSE,
      check_duplicates = FALSE,
      perplexity = perplexity,
      num_threads = 1 ## NOTE: multi-threading may have MEM problems...
    )
    pos <- res1$Y
    rownames(pos) <- colnames(X)
    pos <- pos[1:dimx[2], ] ## if augmented
    colnames(pos) <- paste0("tSNE-", c("x", "y"))
    all.pos[["tsne2d"]] <- pos
    remove(res1)
  }

  if ("tsne" %in% methods && 3 %in% dims) {
    message("calculating t-SNE 3D...")

    perplexity <- pmax(min(dimx[2] / 4, perplexity), 2)
    perplexity
    pos <- Rtsne::Rtsne(t(X),
      dims = 3,
      is_distance = FALSE,
      check_duplicates = FALSE,
      perplexity = perplexity,
      num_threads = 1 ## NOTE: multi-threading may have MEM problems...
    )$Y
    rownames(pos) <- colnames(X)
    pos <- pos[1:dimx[2], ] ## if augmented
    colnames(pos) <- paste0("tSNE-", c("x", "y", "z"))
    all.pos[["tsne3d"]] <- pos
  }

  if ("umap" %in% methods && 2 %in% dims) {
    message("calculating UMAP 2D...")
    if (umap.pkg == "uwot") {
      nb <- ceiling(pmax(min(dimx[2] / 4, perplexity), 2))
      pos <- uwot::tumap(t(X),
        n_components = 2,
        n_neighbors = nb,
        local_connectivity = ceiling(nb / 15)
      )
    } else {
      custom.config <- umap.defaults
      custom.config$n_components <- 2
      custom.config$n_neighbors <- pmax(min(dimx[2] / 4, perplexity), 2)
      pos <- umap::umap(t(X), custom.config)$layout
    }
    rownames(pos) <- colnames(X)
    pos <- pos[1:dimx[2], ] ## if augmented
    colnames(pos) <- paste0("UMAP-", c("x", "y"))
    all.pos[["umap2d"]] <- pos
  }

  if ("umap" %in% methods && 3 %in% dims) {
    message("calculating UMAP 3D...")
    if (umap.pkg == "uwot") {
      nb <- ceiling(pmax(min(dimx[2] / 4, perplexity), 2))
      pos <- uwot::tumap(t(X),
        n_components = 3,
        n_neighbors = nb,
        local_connectivity = ceiling(nb / 15)
      )
    } else {
      custom.config <- umap.defaults
      custom.config$n_components <- 3
      custom.config$n_neighbors <- pmax(min(dimx[2] / 4, perplexity), 2)
      pos <- umap::umap(t(X), custom.config)$layout
    }
    rownames(pos) <- colnames(X)
    pos <- pos[1:dimx[2], ] ## if augmented
    colnames(pos) <- paste0("UMAP-", c("x", "y", "z"))
    all.pos[["umap3d"]] <- pos
  }

  ## NOTE: this pacmap package seem to require GB memory for importing
  ## the python package and calling. Need more testing!!!
  if (FALSE && "pacmap" %in% methods) {
    has.pacmap <- reticulate::py_module_available("pacmap")
    if (has.pacmap && 2 %in% dims) {
      message("calculating PACMAP 2D...")
      ## reticulate::py_install("pacmap")
      pacmap <- reticulate::import("pacmap")
      reducer <- pacmap$PaCMAP(n_components = 2L)
      pos <- reducer$fit_transform(t(X))
      rownames(pos) <- colnames(X)
      pos <- pos[1:dimx[2], ] ## if augmented
      colnames(pos) <- paste0("PACMAP-", c("x", "y"))
      all.pos[["pacmap2d"]] <- pos
    }
    if (has.pacmap && 3 %in% dims) {
      message("calculating PACMAP 3D...")
      pacmap <- reticulate::import("pacmap")
      reducer <- pacmap$PaCMAP(n_components = 3L)
      pos <- reducer$fit_transform(t(X))
      rownames(pos) <- colnames(X)
      pos <- pos[1:dimx[2], ] ## if augmented
      colnames(pos) <- paste0("PACMAP-", c("x", "y", "z"))
      all.pos[["pacmap3d"]] <- pos
    }
  }

  all.pos$membership <- NULL
  if (find.clusters) {
    message("*** DEPRECATED *** please call seperately")
    message("calculating Louvain memberships (from reduced X)...")
    idx <- pgx.findLouvainClusters(t(X), level = 1, prefix = "C", small.zero = 0.01)
    all.pos$membership <- idx[1:dimx[2]]
  }

  gc()
  return(all.pos)
}

#' @export
pgx.clusterBigMatrix <- function(...) pgx.clusterMatrix(...)


#' Find Louvain clusters
#'
#' @title Find Louvain clusters
#'
#' @description Find clusters in a dataset using the Louvain algorithm.
#'
#' @param X A data matrix, data frame, or object that can be coerced to a matrix.
#' @param graph.method Method for building graph from data. Either "dist" or "snn".
#' @param level Number of levels to build hierarchy.
#' @param prefix Prefix for cluster labels.
#' @param gamma Exponent for distance matrix.
#' @param small.zero Proportion of smallest clusters to reassign to 0.
#'
#' @details Builds a graph from the input data and runs the Louvain algorithm to find clusters.
#' The number of levels in the hierarchy can be specified. Small clusters can optionally be
#' reassigned to cluster 0.
#'
#' @return Vector of cluster assignments
#'
#' @export
pgx.findLouvainClusters <- function(X, graph.method = "dist", level = 1, prefix = "C",
                                    knn=100, gamma = 1, small.zero = 0.01) {
  ## find clusters from t-SNE positions
  idx <- NULL
  message("Finding clusters using Louvain...\n")
  knn <- min(round(nrow(X)/5),knn)

  if (graph.method == "dist") {
    dist <- stats::as.dist(stats::dist(scale(X)))
    adjmatrix <- 1.0 / dist**gamma  ## 'power' like WGCNA
    gr <- igraph::graph_from_adjacency_matrix(
      as.matrix(adjmatrix),
      diag = FALSE,
      mode = "undirected"
    )
  } else if (graph.method == "snn") {
    suppressMessages(suppressWarnings(gr <- scran::buildSNNGraph(t(X), d = knn)))
  } else if (graph.method == "knn") {
    suppressMessages(suppressWarnings(gr <- graph_from_knn(X, k = knn)))
  } else {
    stop("FATAL: unknown graph method ", graph.method)
  }

  ## should we iteratively cluster (louvain)???
  hc <- hclustGraph(gr, k = level)
  idx <- hc[, min(level, ncol(hc))]

  ## set small cluster to 'zero'
  if (!is.null(idx) && small.zero > 0) {
    ## ------------ zap small clusters to "0"
    sort(table(idx))
    min.size <- pmax(3, small.zero * length(idx))
    small.clusters <- names(which(table(idx) < min.size))
    idx[which(idx %in% small.clusters)] <- "0"
  }

  ## rename levels with largest cluster first
  idx <- factor(idx, levels = names(sort(-table(idx))))
  levels(idx) <- paste0(prefix, 1:length(levels(idx)))
  idx <- as.character(idx)
  message("Found ", length(unique(idx)), " clusters...")
  return(idx)
}

#' @export
pacmap <- function(X, n_components = 2L, ...) {
  has.pacmap <- reticulate::py_module_available("pacmap")
  if (!has.pacmap) {
    stop("pacmap python module not installed")
  }
  ## reticulate::py_install("pacmap")
  pacmap <- reticulate::import("pacmap")
  reducer <- pacmap$PaCMAP(n_components = as.integer(n_components))
  pos <- reducer$fit_transform(X)
  rownames(pos) <- rownames(X)
  colnames(pos) <- paste0("pacmap_", 1:ncol(pos))
  pos
}

#' Create KNN igraph from position matrix.
#'
graph_from_knn <- function(pos, k = 10) {
  ## return: edge weight are distances
  if (is.null(rownames(pos))) stop("pos must have rownames")
#  if (ncol(pos) > 3 || NCOL(pos) == 1) {
#    stop("positions must be 2 or 3 columns\n")
#  }
  ## use fast KNN package
  res <- FNN::get.knn(pos, k = k)
  idx <- res$nn.index
  xval <- res$nn.dist
  xval <- as.vector(unlist(xval))
  sp.idx <- do.call(rbind, lapply(1:nrow(idx), function(i) cbind(i, idx[i, ])))
  sp <- Matrix::sparseMatrix(i = sp.idx[, 1], j = sp.idx[, 2], x = xval, dims = c(nrow(pos), nrow(pos)))
  sp <- (sp + Matrix::t(sp)) / 2
  rownames(sp) <- colnames(sp) <- rownames(pos)
  g <- igraph::graph_from_adjacency_matrix(sp, mode = "undirected", diag = FALSE, weighted = TRUE)
#  g$layout <- pos
  return(g)
}
