##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


#' @title Convert Seurat to PGX
#'
#' @param obj Seurat object to convert
#' @param do.cluster Logical indicating whether to cluster samples. Default is FALSE.
#'
#' @return PGX object
#'
#' @description Converts a Seurat single-cell RNA-seq object into a PGX object
#'
#' @details This function takes a Seurat object containing single-cell RNA-seq data and converts it into a PGX object.
#' The count matrix, normalized expression matrix, and sample metadata are extracted from the Seurat object.
#' Gene annotations are added using the gene symbols.
#'
#' If do.cluster=TRUE, dimensionality reduction and clustering of samples is performed.
#' Any existing tsne/umap embeddings and cluster assignments are copied over from the Seurat object.
#'
#' @export
seurat2pgx <- function(obj, species, do.cluster = FALSE) {
  ## Convert a Seurat object to a minimal PGX object.
  message("[createPGX.10X] creating PGX object...")
  pgx <- list()
  pgx$name <- "SeuratProject"
  pgx$description <- "Seurat object converted using seurat2pgx"
  pgx$date <- Sys.Date()
  pgx$datatype <- "scRNA-seq"
  pgx$counts <- obj[["RNA"]]@counts
  pgx$X <- obj[["RNA"]]@data
  pgx$samples <- obj@meta.data

  mart <- biomaRt::useMart(biomart = "ensembl", dataset = species)
  pgx$genes <- ngs.getGeneAnnotation(rownames(pgx$counts), mart = mart, organism = pgx$organism)

  if (do.cluster) {
    message("[seurat2pgx] clustering samples")
    pgx <- pgx.clusterSamples2(
      pgx,
      dims = c(2, 3), methods = c("pca", "tsne", "umap")
    )
    names(pgx$cluster$pos)
  }

  ## copy clustering from Seurat
  if ("tsne" %in% names(obj@reductions)) {
    pos <- obj@reductions[["tsne"]]@cell.embeddings[, 1:2]
    pgx$cluster$pos[["tsne2d"]] <- pos
    pgx$tsne2d <- pos
    pgx$tsne3d <- cbind(pos, 0)
  }
  if ("umap" %in% names(obj@reductions)) {
    pgx$cluster$pos[["umap2d"]] <- obj@reductions[["umap"]]@cell.embeddings[, 1:2]
  }
  pgx$samples$cluster <- obj@meta.data[, "seurat_clusters"]
  names(pgx)
  return(pgx)
}


#' @title Pool single cell counts into pseudo-bulk groups
#'
#' @description This function pools single cell RNA-seq count matrices into
#' pseudo-bulk groups. Cells are assigned to groups based on metadata or
#' clustering. Counts within each group are aggregated using summation or other statistics.
#'
#' @param counts Single cell count matrix with cells as columns.
#' @param ncells Number of cells to sample per group. If NULL, takes all cells in each group.
#' @param groups Group assignments for each cell. If NULL, a single group is used.
#' @param stats Aggregation method for pooling counts within each group (sum, mean, median, etc)
#' @param clust.method Clustering method used to assign groups if groups=NULL (umap, pca, etc)
#' @param prior Per-group prior used for pooling. Default is 1 (no prior).
#' @param X Optional log-expression matrix for clustering if groups=NULL.
#' @param meta Optional cell metadata for clustering if groups=NULL.
#' @param verbose Print progress messages?
#'
#' @details This function takes a single cell count matrix and pools the counts into pseudo-bulk groups.
#' If groups are provided, cells are assigned to these groups. Otherwise clustering is performed using umap/pca
#' on the log-expression matrix X to infer groups. Within each group, cell counts are aggregated into a pseudo-bulk profile
#' using summation or other statistics. A per-group prior can be used to normalize pooling.
#' The output is a pooled count matrix with groups as columns.
#'
#' @return Matrix of pooled counts for each group
#'
#' @export
pgx.poolCells <- function(counts, ncells, groups = NULL, stats = "sum",
                          clust.method = "umap", prior = 1, X = NULL,
                          meta = NULL, verbose = TRUE) {
  if (is.null(groups)) {
    groups <- rep("grp1", ncol(counts))
  }
  groups <- as.integer(factor(as.character(groups)))

  pool.counts <- c()
  cluster.id <- c()

  if (is.null(X)) {
    ## compute log-expression matrix from counts
    X <- counts
    X <- 1e6 * t(t(X) / (1e-8 + Matrix::colSums(X))) ## CPM
    nz <- Matrix::which(X != 0, arr.ind = TRUE)
    X[nz] <- log2(1 + X[nz])
  }

  clusterX <- function(X1, k, method, topsd = 1000, nv = 50) {
    if (ncol(X1) <= k) {
      cluster <- paste0("c", 1:ncol(X1))
      names(cluster) <- colnames(X1)
      return(cluster)
    }
    avgx <- Matrix::rowMeans(X1)

    sdx <- (Matrix::rowMeans(X1**2) - avgx**2)**0.5
    wt <- sdx * avgx
    X1 <- X1[utils::head(order(-wt), topsd), , drop = FALSE]

    X1 <- X1 - Matrix::rowMeans(X1) ## center features
    nv <- min(nv, ncol(X1) - 1)
    sv <- irlba::irlba(X1, nv = nv)
    V <- sv$v %*% diag(sv$d**0.5)
    V <- t(scale(t(V))) ## rowscale
    cluster <- NULL
    if (method == "tsne") {
      px <- min(nrow(V) / 4, 30)
      pos <- Rtsne::Rtsne(V, perplexity = px)$Y
      cluster <- stats::kmeans(pos, k, iter.max = 100)$cluster
    } else if (method == "umap") {
      pos <- uwot::umap(V) ## umap is tighter and faster than t-SNE
      cluster <- stats::kmeans(pos, k, iter.max = 100)$cluster
    } else if (method == "kmeans") {
      cluster <- stats::kmeans(V, k, iter.max = 100)$cluster
    } else if (method == "hclust") {
      cluster <- stats::cutree(fastcluster::hclust(stats::dist(V)), k)
    } else {
      stop("ERROR:: unknown clustering method")
    }
    ## cluster <- stats::kmeans(po s,k)$cluster
    cluster <- paste0("c", cluster)
    names(cluster) <- colnames(X1)
    cluster
  }

  ngroup <- length(unique(groups))
  if (verbose) {
    message("[pgx.poolCells] ", ngroup, " groups defined")
    message("[pgx.poolCells] clustering method = ", clust.method)
  }

  cluster.id <- rep(NA, ncol(counts))
  names(cluster.id) <- colnames(counts)
  all.groups <- sort(unique(groups))
  g <- groups[1]
  for (g in all.groups) {
    if (verbose) message("[pgx.poolCells] clustering group", g, " ...")
    ## do quick clustering
    sel <- which(groups == g)
    X1 <- X[, sel, drop = FALSE]
    k <- ceiling(ncells / ngroup) ## equalize between groups
    k
    ngroup
    cluster <- NULL
    if (ncol(X1) <= k) {
      cluster <- paste0("g", g, "-c", 1:ncol(X1))
    } else {
      cluster <- clusterX(X1, k = k, method = clust.method)
      if (!is.null(groups)) cluster <- paste0("g", g, "-", cluster)
    }
    names(cluster) <- colnames(X1)
    cluster.id[colnames(X1)] <- cluster
  }

  if (verbose) message("[pgx.poolCells] pooling cells...")
  cluster.id <- cluster.id[colnames(counts)]
  if (stats == "mean") {
    pool.counts <- tapply(1:ncol(counts), cluster.id, function(ii) {
      Matrix::rowMeans(counts[, ii, drop = FALSE])
    })
  }
  if (stats == "sum") {
    pool.counts <- tapply(1:ncol(counts), cluster.id, function(ii) {
      Matrix::rowSums(counts[, ii, drop = FALSE])
    })
  }
  pool.counts <- do.call(cbind, pool.counts)

  new.meta <- NULL
  if (!is.null(meta)) {
    max.class <- function(x) names(which.max(table(x)))
    new.meta <- apply(meta, 2, function(x) {
      tapply(x, cluster.id, max.class)
    })
    new.meta <- data.frame(new.meta, check.names = FALSE)
  }

  res <- list()
  res$cluster.id <- cluster.id
  res$counts <- pool.counts
  res$meta <- new.meta
  return(res)
}


#' @title Integrate single-cell data across batches
#'
#' @param X Numeric matrix of expression values, cells as columns
#' @param batch Factor specifying batch for each cell
#' @param method Methods to use for batch integration. Options are "ComBat", "limma", "CCA", "MNN", "Harmony", "liger"
#'
#' @return List containing batch-integrated expression matrices by method
#'
#' @description Integrate single-cell RNA-seq data from multiple batches using various batch correction methods.
#'
#' @details This function takes a single-cell expression matrix \code{X} and a batch vector \code{batch} as input.
#' It applies different batch correction methods to integrate the data across batches.
#'
#' The following methods can be selected via the \code{method} parameter:
#' \itemize{
#' \item ComBat: Apply ComBat batch correction from sva package
#' \item limma: Apply removeBatchEffect from limma package
#' \item CCA: Apply canonical correlation analysis using Seurat package
#' \item MNN: Apply mutual nearest neighbors correction
#' \item Harmony: Apply Harmony integration using Harmony package
#' \item liger: Apply integration via liger package
#' }
#'
#' The batch-integrated expression matrices are returned as a list by method name.
#'
#' @export
pgx.scBatchIntegrate <- function(X, batch,
                                 method = c("ComBat", "limma", "CCA", "MNN", "Harmony", "liger")) {
  res <- list()

  if ("ComBat" %in% method) {
    message("[pgx.scBatchIntegrate] single-cell batch correction using ComBat...")

    ## ComBat correction
    res[["ComBat"]] <- sva::ComBat(X, batch = batch, par.prior = TRUE)
  }
  if ("limma" %in% method) {
    message("[pgx.scBatchIntegrate] single-cell batch correction using LIMMA...")

    ## LIMMA correction
    res[["limma"]] <- limma::removeBatchEffect(X, batch = batch)
  }
  if ("CCA" %in% method) {
    message("[pgx.scBatchIntegrate] single-cell batch correction using CCA (Seurat)...")
    ## Seurat CCA correction
    counts <- pmax(2**X - 1, 0)
    try(res[["CCA"]] <- pgx.SeuratBatchIntegrate(counts, batch = batch))
  }
  if ("MNN" %in% method) {
    message("[pgx.scBatchIntegrate] single-cell batch correction using MNN...")
    ## MNN correction
    try(mnn <- batchelor::mnnCorrect(X, batch = batch, cos.norm.in = TRUE, cos.norm.out = FALSE))
    res[["MNN"]] <- MultiAssayExperiment::assays(mnn)[["corrected"]]
  }
  if ("Harmony" %in% method) {
    message("[pgx.scBatchIntegrate] single-cell batch correction using Harmony...")

    ## Harmony corrections
    nv <- min(floor(dim(X) * 0.8), 30)
    out <- irlba::irlba(X, nu = nv, nv = nv)
    V <- t(out$v)
    meta_data <- data.frame(batch = batch)
    try(hm <- harmony::HarmonyMatrix(
      V, meta_data, "batch",
      do_pca = FALSE, npcs = nv,
      return_object = TRUE
    ))
    ## corrected PCA embeddings
    hX <- (out$u %*% diag(out$d) %*% hm$Z_corr)
    dimnames(hX) <- dimnames(X)
    res[["Harmony"]] <- hX
  }
  if ("liger" %in% method) {
    message("[pgx.scBatchIntegrate] single-cell batch correction using LIGER...")

    xlist <- tapply(1:ncol(X), batch, function(i) pmax(2**X[, i] - 1, 0))
    liger <- rliger::createLiger(xlist, take.gene.union = TRUE)
    liger <- rliger::normalize(liger)


    liger@var.genes <- Matrix::head(rownames(X)[order(-apply(X, 1, stats::sd))], 100)
    liger <- rliger::scaleNotCenter(liger)
    vg <- liger@var.genes
    xdim <- sapply(xlist, ncol)
    k <- 15
    k <- round(min(30, length(vg) / 3, stats::median(xdim / 2)))
    ## OFTEN GIVES ERROR!!!!!
    liger <- try(rliger::optimizeALS(liger, k = k))


    if (inherits(liger, "try-error")) {

    } else {
      liger <- rliger::quantile_norm(liger)
      cX <- t(liger@H.norm %*% liger@W)
      cat("[pgx.scBatchIntegrate] WARNING:: LIGER returns smaller matrix")
      res[["liger"]] <- cX
    }
  }
  if (length(res) == 1) {
    res <- res[[1]]
  }
  return(res)
}


#' @title Integrate single-cell data with Seurat
#'
#' @param counts Single-cell count matrix
#' @param batch Batch vector assigning batches to cells
#' @param qc.filter Logical indicating whether to filter cells by QC metrics. Default is FALSE.
#' @param nanchors Number of anchor points to use for integration. Default is -1 to auto-determine.
#' @param sct Logical indicating whether to use SCTransform normalization. Default is FALSE.
#'
#' @return Seurat object containing integrated data
#'
#' @description Integrate single-cell RNA-seq data from multiple batches using canonical correlation analysis via the Seurat package.
#'
#' @details This function takes a single-cell count matrix \code{counts} and a \code{batch} vector as input.
#' It sets up a Seurat object for each batch and integrates them using FindIntegrationAnchors/IntegrateData functions.
#'
#' Cells can be filtered by QC metrics like mitochondrial content if \code{qc.filter=TRUE}.
#' The \code{nanchors} parameter controls the number of anchor points used for integration.
#' Normalization and scaling can be done using SCTransform if \code{sct=TRUE}.
#'
#' The integrated Seurat object is returned containing the corrected expression matrix.
#'
#' @export
pgx.SeuratBatchIntegrate <- function(counts, batch, qc.filter = FALSE,
                                     nanchors = -1, sct = FALSE) {
  ## From Seurat vignette: Integration/batch correction using
  ## CCA. Note there is no QC filtering for samples on ribo/mito
  ## content. You need to do that before.

  nbatch <- length(unique(batch))
  message("[pgx.SeuratBatchIntegrate] Processing ", nbatch, " batches...")
  obj.list <- list()
  i <- 1
  b <- batch[1]
  batches <- unique(batch)
  for (i in 1:length(batches)) {
    sel <- which(batch == batches[i])

    counts1 <- counts[, sel]
    obj <- Seurat::CreateSeuratObject(counts1)
    if (sct) {
      obj <- Seurat::SCTransform(obj, vars.to.regress = NULL, verbose = FALSE)
    } else {
      obj <- Seurat::NormalizeData(obj,
        normalization.method = "LogNormalize",
        scale.factor = 10000, verbose = FALSE
      )
      obj <- Seurat::FindVariableFeatures(obj, selection.method = "vst", verbose = FALSE)
    }
    obj$batch <- b
    obj.list[[i]] <- obj
  }

  anchor.features <- NULL
  if (sct) {
    ## See: https://satijalab.org/seurat/v3.0/integration.html
    nfeatures <- nrows(counts)
    if (nanchors > 0) nfeatures <- nanchors
    anchor.features <- Seurat::SelectIntegrationFeatures(
      object.list = obj.list, nfeatures = nanchors
    )
    obj.list <- Seurat::PrepSCTIntegration(
      object.list = obj.list,
      anchor.features = anchor.features,
      verbose = FALSE
    )
  } else {
    anchor.features <- nrow(counts) ## really?
    if (nanchors > 0) anchor.features <- nanchors ## really?
  }

  message("[pgx.SeuratBatchIntegrate] Finding anchors...")
  options(future.globals.maxSize = 8 * 1024^3) ## set to 8GB

  NUM.CC <- max(min(20, min(table(batch)) - 1), 1)
  bdims <- sapply(obj.list, ncol)
  kmax <- max(min(bdims) - 1, 1)
  normalization.method <- ifelse(sct, "SCT", "LogNormalize")
  message("[pgx.SeuratBatchIntegrate] NUM.CC = ", NUM.CC)
  message("[pgx.SeuratBatchIntegrate] normalization.method = ", normalization.method)

  anchors <- Seurat::FindIntegrationAnchors(
    obj.list,
    dims = 1:NUM.CC,
    k.filter = min(200, kmax),
    k.anchor = min(5, kmax),
    k.score = min(30, kmax),
    anchor.features = anchor.features,
    normalization.method = normalization.method,
    verbose = FALSE
  )

  message("[pgx.SeuratBatchIntegrate] Integrating data...")
  len.anchors <- length(anchors@anchor.features)
  message("[pgx.SeuratBatchIntegrate] number of anchors = ", len.anchors)

  integrated <- Seurat::IntegrateData(anchorset = anchors)
  integrated <- Seurat::IntegrateData(
    anchorset = anchors,
    k.weight = min(100, kmax), ##  troublesome...
    dims = 1:NUM.CC,
    normalization.method = normalization.method,
    verbose = FALSE
  )

  key <- ifelse(sct, "SCT", "integrated")
  mat.integrated <- as.matrix(integrated[[key]]@data)
  mat.integrated <- exp(mat.integrated) ## natural log!!
  mat.integrated <- mat.integrated[, colnames(counts)]

  ## set previously zero counts to zero again
  zc <- Matrix::which(counts[rownames(mat.integrated), ] == 0, arr.ind = TRUE)
  mat.integrated[zc] <- 0

  ## set missing to zero...

  mat.integrated[is.na(mat.integrated)] <- 0

  if (!nrow(mat.integrated) == nrow(counts)) {
    cat("WARNING:: number of rows of integrated matrix has changed!")
  }

  return(mat.integrated)
}

#' @export
pgx.read_singlecell_counts <- function(filename) {
  counts <- NULL
  if (grepl("[.]csv$", filename)) {
    counts <- as.matrix(data.table::fread(filename, header = TRUE), row.names = 1)
  }
  if (grepl("[.]mtx$", filename)) {
    dir <- dirname(filename)
    barcode.file <- file.path(dir, "barcodes.tsv")
    genes.file <- file.path(dir, "genes.tsv")
    if (!file.exists(filename)) stop("could not find counts matrix: ", filename)
    if (!file.exists(barcode.file)) stop("could not find barcode file: ", barcode.file)
    if (!file.exists(genes.file)) stop("could not find genes file: ", genes.file)
    counts <- Matrix::readMM(filename)
    bc <- read.csv(barcode.file, header = FALSE, sep = "\t")
    gn <- read.csv(genes.file, header = FALSE, sep = "\t")
    rownames(counts) <- gc[, 2] ## gene names?
    colnames(counts) <- bc[, 1]
  }
  if (grepl("[.]h5$", filename)) {
    counts <- Seurat::Read10X_h5(filename, use.names = TRUE, unique.features = TRUE)
  }
  dim(counts)
  counts
}


#' @export
pgx.supercell <- function(counts, meta, group = NULL, gamma = 20) {
  ## require(SuperCell)
  
  X <- log2(1 + edgeR::cpm(counts) / 100)
  if (is.null(group) && "group" %in% colnames(meta)) {
    cat("using group detected in meta\n")
    group <- meta[, "group"]
  }

  SC <- SuperCell::SCimplify(X,
    gamma = gamma,
    n.var.genes = 1000,
    ## cell.annotation = group,
    cell.split.condition = group
  )

  dsel <- which(sapply(meta, class) %in% c("factor", "character", "logical"))
  group.argmax <- function(x) tapply(x, SC$membership, function(x) names(which.max(table(x))))
  dmeta <- apply(meta[, dsel, drop = FALSE], 2, function(x) as.character(group.argmax(x)))
  rownames(dmeta) <- sort(unique(SC$membership))
  
  csel <- which(sapply(meta, class) %in% c("numeric", "integer"))
  group.mean <- function(x) tapply(x, SC$membership, function(x) mean(x, na.rm = TRUE))
  cmeta <- apply(meta[, csel, drop = FALSE], 2, function(x) group.mean(x))

  sc.meta <- data.frame(dmeta)
  if (length(csel) > 0) sc.meta <- cbind(sc.meta, cmeta)
  ii <- setdiff(match(colnames(meta), colnames(sc.meta)), NA)
  sc.meta <- sc.meta[, ii]

  ## Compute metacall expression as sum of counts
  sc.counts <- SuperCell::supercell_GE(counts, mode = "sum", groups = SC$membership)
  sc.membership <- paste0("mc", SC$membership)
  colnames(sc.counts) <- paste0("mc", 1:ncol(sc.counts))
  rownames(sc.meta) <- colnames(sc.counts)

  list(counts = sc.counts, meta = sc.meta, membership = sc.membership)
}

#' @export
pgx.supercellX <- function(X, meta, group = NULL, gamma = 20) {
  ## require(SuperCell)
  
  if (is.null(group) && "group" %in% colnames(meta)) {
    cat("using group detected in meta\n")
    group <- meta[, "group"]
  }

  SC <- SuperCell::SCimplify(X,
    gamma = gamma,
    n.var.genes = 1000,
    ## cell.annotation = group,
    cell.split.condition = group
  )

  dsel <- which(sapply(meta, class) %in% c("factor", "character", "logical"))
  idx <- SC$membership
  group.argmax <- function(x) tapply(x, idx, function(x) names(which.max(table(x))))
  dmeta <- apply(meta[, dsel, drop = FALSE], 2, function(x) as.character(group.argmax(x)))
  rownames(dmeta) <- sort(unique(idx))

  csel <- which(sapply(meta, class) %in% c("numeric", "integer"))
  group.mean <- function(x) tapply(x, idx, function(x) mean(x, na.rm = TRUE))
  cmeta <- apply(meta[, csel, drop = FALSE], 2, function(x) group.mean(x))

  sc.meta <- data.frame(dmeta)
  if (length(csel) > 0) sc.meta <- cbind(sc.meta, cmeta)
  ii <- setdiff(match(colnames(meta), colnames(sc.meta)), NA)
  sc.meta <- sc.meta[, ii]

  ## Compute metacall expression as sum of counts
  sc.X <- SuperCell::supercell_GE(X, mode = "average", groups = idx)
  colnames(sc.X) <- paste0("mc", 1:ncol(sc.X))
  rownames(sc.meta) <- colnames(sc.X)

  list(X = sc.X, meta = sc.meta, membership = sc.membership)
}

pgx.supercell_BIG <- function(counts, meta, group = NULL, gamma = 20, batch.size = 1e5) {
  ## require(SuperCell)
  if (is.null(group) && "group" %in% colnames(meta)) {
    cat("using group detected in meta\n")
    group <- meta[, "group"]
  }

  nbatch <- ceiling(ncol(counts) / batch.size)
  idx <- head(as.vector(sapply(1:nbatch, rep, batch.size)), ncol(counts))
  ## idx <- sample(idx)
  table(idx)

  i <- 1
  sc.counts <- c()
  sc.membership <- c()
  gamma

  i <- 1
  for (i in sort(unique(idx))) {
    cat("processing index", i, "\n")
    ii <- which(idx == i)

    ## https://stackoverflow.com/questions/39284774/column-rescaling-for-a-very-large-sparse-matrix-in-r
    ## X <- log2(1 + edgeR::cpm(X) / 100)
    X1 <- as(counts[, ii], "sparseMatrix")
    X1@x <- X1@x / rep.int(colSums(X1), diff(X1@p)) * 1e4
    X1@x <- log2(1 + X1@x)

    grp <- NULL
    if (!is.null(group)) grp <- group[ii]

    SC <- SuperCell::SCimplify(X1,
      gamma = gamma, n.var.genes = 1000,
      ## cell.annotation = Y$Mouse[ii]
      cell.split.condition = grp
    )

    ## Compute metacall expression as sum of counts
    n0 <- length(unique(unlist(sc.membership)))
    sc1 <- SuperCell::supercell_GE(counts[, ii], mode = "sum", groups = SC$membership)
    colnames(sc1) <- paste0("mc", n0 + 1:ncol(sc1))

    sc.counts <- cbind(sc.counts, sc1)
    sc.membership <- c(sc.membership, paste0("mc", n0 + SC$membership))
    ## sc.membership[[i]] <- paste0("mc.",n0 + SC$membership)
  }

  dim(sc.counts)
  table(sc.membership)

  meta <- meta[, colMeans(is.na(meta)) < 1, drop = FALSE]

  dsel <- which(sapply(meta, class) %in% c("factor", "character", "logical"))
  group.argmax <- function(x) tapply(x, sc.membership, function(x) names(which.max(table(x))))
  dmeta <- apply(meta[, dsel, drop = FALSE], 2, function(x) group.argmax(x))

  csel <- which(sapply(meta, class) %in% c("numeric", "integer"))
  group.mean <- function(x) tapply(x, sc.membership, function(x) mean(x, na.rm = TRUE))
  cmeta <- apply(meta[, csel, drop = FALSE], 2, function(x) group.mean(x))

  sc.meta <- data.frame(dmeta)
  if (length(csel) > 0) sc.meta <- cbind(sc.meta, cmeta)
  ii <- setdiff(match(colnames(meta), colnames(sc.meta)), NA)
  sc.meta <- sc.meta[, ii]
  mc.name <- paste0("mc", 1:ncol(sc.counts))
  colnames(sc.counts) <- rownames(sc.meta) <- mc.name

  list(counts = sc.counts, meta = sc.meta, membership = sc.membership)
}

pgx.createSeuratObject <- function(counts, samples, sct = TRUE) {
  obj <- Seurat::CreateSeuratObject(counts = counts, meta.data = samples)
  obj <- Seurat::NormalizeData(obj)
  obj <- Seurat::PercentageFeatureSet(obj, pattern = "^MT-|^Mt-", col.name = "percent.mt")
  obj <- Seurat::PercentageFeatureSet(obj, pattern = "^RP[LS]|^Rp[ls]", col.name = "percent.ribo")
  if (sct) {
    sct <- Seurat::SCTransform(obj, method = "glmGamPoi", verbose = FALSE)
    sct <- Seurat::RunPCA(sct, verbose = FALSE)
    sct <- Seurat::FindNeighbors(sct, dims = 1:30, verbose = FALSE)
    sct <- Seurat::FindClusters(sct, verbose = FALSE, resolution = 0.8) ## smaller resolution
    sct <- Seurat::RunUMAP(sct, dims = 1:30, verbose = FALSE)
    obj <- sct
  }
  obj
}

pgx.runAzimuth <- function(counts, reference = "pbmcref") {
  obj <- pgx.createSeuratObject(counts, samples = NULL, sct = FALSE)
  obj1 <- Azimuth::RunAzimuth(obj, reference = reference, verbose = FALSE)
  k1 <- !(colnames(obj1@meta.data) %in% colnames(obj@meta.data))
  k2 <- !grepl("score$|refAssay$", colnames(obj1@meta.data))
  obj1@meta.data[, (k1 & k2)]
}




## =====================================================================================
## =========================== END OF FILE =============================================
## =====================================================================================
