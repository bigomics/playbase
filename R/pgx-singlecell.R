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
seurat2pgx <- function(obj, do.cluster = FALSE, organism) {
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

  probes <- rownames(pgx$counts)
  pgx$genes <- ngs.getGeneAnnotation(rownames(pgx$counts), organism = organism)

  if (do.cluster) {
    message("[seurat2pgx] clustering samples")
    pgx <- pgx.clusterSamples(
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
      Matrix::rowSums(counts[, ii, drop = FALSE], na.rm = TRUE)
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
    res[["ComBat"]] <- sva::ComBat(X, batch = batch, par.prior = TRUE)
  }
  if ("limma" %in% method) {
    message("[pgx.scBatchIntegrate] single-cell batch correction using LIMMA...")
    res[["limma"]] <- limma::removeBatchEffect(X, batch = batch)
  }
  if ("CCA" %in% method) {
    message("[pgx.scBatchIntegrate] single-cell batch correction using CCA (Seurat)...")
    counts <- pmax(2**X - 1, 0)
    CCA <- try(pgx.SeuratBatchIntegrate(counts, batch = batch))    
    if (!inherits(CCA, "try-error")) {
      res[["CCA"]] <- CCA
    }
  }
  if ("MNN" %in% method) {
    message("[pgx.scBatchIntegrate] single-cell batch correction using MNN...")
    ##mnn <- try(batchelor::mnnCorrect(X, batch = batch, cos.norm.in = TRUE, cos.norm.out = FALSE))
    mnn <- try(batchelor::mnnCorrect(X, batch = batch))
    if (!inherits(mnn, "try-error")) {
      res[["MNN"]] <- MultiAssayExperiment::assays(mnn)[["corrected"]]
    }
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
    if (!inherits(liger, "try-error")) {
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

#' @title SuperCell down sampling. Uniform down samplsing using gamma = 20.
#' 
#' @export
pgx.supercell <- function(counts, meta, group = NULL, gamma = 20) {
  ## require(SuperCell)

  ## supercell uses log2 matrix
  ##X <- log2(1 + edgeR::cpm(counts) / 100)
  X <- logCPM(counts, total=1e4)
  if (!is.null(group) && "group" %in% colnames(meta)) {
    cat("using group detected in meta\n")
    group <- meta[, "group"]
  }
  
  SC <- SuperCell::SCimplify(X,
    gamma = gamma,
    n.var.genes = 1000,
    ## cell.annotation = group,
    cell.split.condition = group
  )

  meta <- as.data.frame(meta)
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
  sc.meta <- sc.meta[, ii, drop=FALSE]
  
  ## Compute metacall expression as sum of counts
  sc.counts <- SuperCell::supercell_GE(counts, mode = "sum", groups = SC$membership)
  sc.membership <- paste0("mc", SC$membership)
  colnames(sc.counts) <- paste0("mc", 1:ncol(sc.counts))
  rownames(sc.meta) <- colnames(sc.counts)
  
  list(counts = sc.counts, meta = sc.meta, membership = sc.membership)
}


#' @title SuperCell downsampling by group targeting equal sizes per group.
#' 
#' @export
pgx.supercell2 <- function(counts, meta, group, target_n = 20) {

  ## metacell equal across celltype and phenotype
  g <- group[1]
  sc.membership <- rep(NA,length(group))
  sc.counts <- c()
  sc.meta <- c()
  for(g in unique(group)) {
    sel <- which(group == g)
    k <- max(1, length(sel) / target_n)
    sc <- pgx.supercell( counts[,sel], samples[sel,], gamma=k)
    names(sc)
    colnames(sc$counts) <- paste0(g,".",colnames(sc$counts))
    rownames(sc$meta) <- colnames(sc$counts)
    sc$meta$group <- g
    sc$membership <- paste0(g,".",sc$membership)
    sc.membership[sel] <- sc$membership
    sc.counts <- cbind(sc.counts, sc$counts)
    sc.meta <- rbind(sc.meta, sc$meta)
  }
  dim(sc.counts)

  list(counts = sc.counts, meta = sc.meta, membership = sc.membership)  
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

pgx.sc_anchors <- function( counts, sc.counts, sc.membership,
                           sc.group = colnames(sc.counts) ) {
  
  ## determine closest sample to metacell reference (aka anchor)
  sc.ref <- rep(NA,ncol(sc.counts))
  X1 <- logCPM(counts,1e4)
  X2 <- as.matrix(logCPM(sc.counts,1e4))
  i=1
  for(i in 1:ncol(sc.counts)) {
    mc <- sc.group[i]
    sel <- which(sc.membership == mc)
    x1 <- as.matrix(X1[,sel,drop=FALSE])
    colnames(x1) <- colnames(X1)[sel]
    rmax <- which.max(cor(x1, X2[,i])[,1])
    sc.ref[i] <- colnames(x1)[rmax]
  }
  sc.ref
}

#' @export
pgx.justSeuratObject <- function(counts, samples) {
  options(Seurat.object.assay.calcn = TRUE)
  getOption("Seurat.object.assay.calcn")
  Seurat::CreateSeuratObject(counts = counts, meta.data = samples)
}

#' @export
pgx.createSeuratObject <- function(counts,
                                   samples,
                                   batch,
                                   cellcyclescores = TRUE,
                                   filter = TRUE,
                                   preprocess = TRUE,
                                   method = "Harmony") {

  options(Seurat.object.assay.calcn = TRUE)
  getOption("Seurat.object.assay.calcn")
        
  if (is.null(samples)) {
    samples <- data.frame(row.names = colnames(counts))
  }

  ## samples$batch <- batch
  rownames(samples) <- colnames(counts)
  samples <- as.data.frame(samples)
  samples$nCount_RNA <- Matrix::colSums(counts, na.rm = TRUE)
  samples$nFeature_RNA <- Matrix::colSums(counts>0, na.rm = TRUE)  
  obj <- Seurat::CreateSeuratObject(counts = counts, meta.data = samples)
  obj <- Seurat::PercentageFeatureSet(obj, pattern = "^MT-|^Mt-", col.name = "percent.mt")
  obj <- Seurat::PercentageFeatureSet(obj, pattern = "^RP[LS]|^Rp[ls]", col.name = "percent.ribo")
  obj <- Seurat::PercentageFeatureSet(obj, pattern = "^HB|^Hb", col.name = "percent.hb")

  if (cellcyclescores) {
    ## CellCycleScoring needs normalized data
    counts <- obj@assays$RNA$counts
    totcounts <- Matrix::colSums(counts, na.rm = TRUE)
    cpm <- sweep(counts, 2, totcounts, FUN = "/") * 1e4
    obj@assays$RNA$data <- log1p(cpm)    
    g2m.ff <- Seurat::cc.genes$g2m.genes
    s.ff <- Seurat::cc.genes$s.genes
    obj <- Seurat::CellCycleScoring(obj, g2m.features = g2m.ff, s.features = s.ff)
  }
  
  if(filter) {
    ## Filter on number of counts/features and mitochondrial gene content.
    message("[pgx.createIntegratedSeuratObject] filtering cells")    
    qN <- quantile( obj$nCount_RNA, probs = c(0.01,0.99))
    qF <- quantile( obj$nFeature_RNA, probs = c(0.01,0.99))
    obj <- subset(obj, subset =
                         nCount_RNA > qN[1] & nCount_RNA < qN[2] &                       
                         nFeature_RNA > qF[1] & nFeature_RNA < qF[2] &
                         percent.mt < 5)
    dim(obj)
  }

  if(preprocess) {
    message("[pgx.createIntegratedSeuratObject] preprocessing Seurat object")   
    if(!is.null(batch)) {
      obj <- seurat.integrate(obj, 'batch', sct = TRUE, method = "Harmony") 
    } else {
      obj <- seurat.preprocess(obj, sct = TRUE)
    }
  }

  obj  

}

#' @export
seurat.preprocess <- function(obj,
                              sct = FALSE,
                              tsne = TRUE,
                              umap = TRUE) {

  options(future.globals.maxSize= 4*1024^4)
  
  if(sct) {
    message("[seurat.preprocess] normalization method = SCT")
    obj <- Seurat::SCTransform(obj, method = "glmGamPoi", verbose = FALSE)
  } else {
    # pre-process dataset (without integration)
    ## obj <- Seurat::NormalizeData(obj) ## causes "integer overflow issue"
    counts <- obj@assays$RNA$counts
    totcounts <- Matrix::colSums(counts, na.rm = TRUE)
    cpm <- sweep(counts, 2, totcounts, FUN = "/") * 1e4
    obj@assays$RNA$data <- log1p(cpm)
    obj <- Seurat::FindVariableFeatures(obj)
    obj <- Seurat::ScaleData(obj)
    message("[seurat.preprocess] norm, findVarFeatures and scaling completed.")
  }
  
  message("[seurat.preprocess] running PCA")
  npcs <- min(30, ncol(obj)/2)
  obj <- Seurat::RunPCA(obj, npcs = npcs, verbose = FALSE)
  
  # clustering using integrated data or original pca
  dr <- "pca"
  nn <- min(30L, ncol(obj)/5)
  message("[seurat.preprocess] running FindNeighbors & Clusters")
  obj <- Seurat::FindNeighbors(obj, dims=1:npcs, reduction = dr, verbose = FALSE)
  obj <- Seurat::FindClusters(obj, resolution = 1, verbose = FALSE) 

  if (tsne) {
    message("[seurat.preprocess] running tSNE")
    obj <- Seurat::RunTSNE(obj, dims=1:npcs, reduction = dr, verbose = FALSE)
  }

  if (umap) {
    message("[seurat.preprocess] running UMAP")
    obj <- Seurat::RunUMAP(obj, dims = 1:npcs, n.neighbors = nn,
      reduction = dr, verbose = FALSE)
  }
  obj
}

#' @export
seurat.integrate <- function(obj, batch, sct = TRUE, method = "Harmony") {

  obj[["RNA"]] <- split( obj[["RNA"]], f = obj@meta.data[,batch] )

  if(sct) {
    message("[seurat.integrate] normalization method = SCT")
    obj <- Seurat::SCTransform(obj, method = "glmGamPoi", verbose = FALSE)
  } else {
    # pre-process dataset (without integration)
    obj <- Seurat::NormalizeData(obj)
    obj <- Seurat::FindVariableFeatures(obj)
    obj <- Seurat::ScaleData(obj)
  }
  
  message("[seurat.integrate] running PCA")  
  obj <- Seurat::RunPCA(obj, npcs = 30, verbose = FALSE)

  sel.method <- Seurat::HarmonyIntegration
  sel.method <- switch(
    method,
    CCA = Seurat::CCAIntegration,
    RPCA = Seurat::RPCAIntegration,
    Harmony = Seurat::HarmonyIntegration,
    JointPCA = Seurat::JointPCAIntegration        
  )

  message("[seurat.integrate] integration method = ",method)
  obj <- Seurat::IntegrateLayers(
    object = obj,
    method = sel.method,
    normalization.method = ifelse(sct, "SCT", "LogNormalize"),
    new.reduction = "integrated.dr",
    verbose = TRUE
  )

  # re-join layers after integration
  obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
  
  # clustering using integrated data or original pca
  dr <- "integrated.dr"
  obj <- Seurat::FindNeighbors(obj, dims = 1:30, reduction = dr, verbose = FALSE)
  obj <- Seurat::FindClusters(obj, resolution = 1, verbose = FALSE, ) 
  obj <- Seurat::RunUMAP(obj, dims=1:30, reduction = dr, verbose = FALSE)
  ## obj <- Seurat::RunTSNE(obj, dims=1:30, reduction = dr, verbose = FALSE)
  obj
}


#' @export
pgx.runAzimuth <- function(counts, k.weight = NULL, reference = NULL) {
  require(Seurat)
  ## options(future.globals.maxSize= 4*1024^3)  ## needed  
  options(future.globals.maxSize= 4*1024^4) 
  obj <- pgx.justSeuratObject(counts, samples = NULL)
  if (is.null(k.weight)) {
    k.weight <- 20
    ## k.weight <- round(min(50, ncol(obj)/10))
  }
  dbg("[pgx.runAzimuth] k.weight = ", k.weight)
  if (is.null(reference)) {
    reference <- "pbmcref"
  }
  obj1 <- try(Azimuth::RunAzimuth(
    obj,
    reference = reference,
    k.weight = k.weight, 
    verbose = FALSE
  ))
  if(!"try-error" %in% class(obj1)) {
    k1 <- !(colnames(obj1@meta.data) %in% colnames(obj@meta.data))
    k2 <- !grepl("score$|refAssay$", colnames(obj1@meta.data))
    meta1 <- obj1@meta.data[, (k1 & k2)]
    return(meta1)
  } else {
    dbg("[pgx.runAzimuth] Azimuth failed: ref.atlas might be incorrect.")
    return(NULL)
  }
}


#' @title Random downsampling by group targeting equal sizes per group.
#' 
#' @export
pgx.downsample <- function(counts, meta, group, target_n = 20) {
  ## metacell equal across celltype and phenotype
  g <- group[1]
  sc.membership <- rep(NA,length(group))
  sc.counts <- c()
  sc.meta <- c()
  sel <- tapply( 1:ncol(counts), group, function(ii)
    head( sample(ii), target_n) )
  sel <- unlist(sel)
  sc.counts  <- counts[,sel]
  sc.meta    <- meta[sel,]
  sc.meta$group <- group[sel]
  list(counts = sc.counts, meta = sc.meta, membership=NULL)  
}

seurat.downsample <- function(obj, target_n = 1000, target_g = 20, group = NULL) {
  if( ncol(obj) <= target_n ) return(obj)
  if(is.null(group)) {
    sub <- obj[ , sample(ncol(obj), target_n)]
    return(sub)
  }
  obj@meta.data$downsample.group <- group
  sel <- tapply(1:ncol(obj), group, function(ii) head(sample(ii), target_g))
  sub <- obj[ , unlist(sel)]
  return(sub)
}


#' Transfer labels from reference to query (Seurat objects)
#' 
seurat.transferLabels <- function(ref, query, labels, sct=TRUE) {

  
  ## find anchors
  npcs <- ncol(ref@reductions[[1]])
  nn <- min(30L, ncol(ref)/5)
  anchors <- Seurat::FindTransferAnchors(
    reference = ref, query = query, dims = 1:npcs, k.score = nn,
    normalization.method = ifelse(sct, "SCT", "LogNormalize"),
    reference.reduction = "pca")

  ## redudant??? MapQuery does the same??
  if(length(labels)==1 && labels %in% colnames(ref@meta.dat)) {
    refdata <- ref@meta.data[,labels]
  } else {
    refdata <- labels      
  }
  predictions <- Seurat::TransferData(anchorset = anchors, refdata = refdata, dims = 1:npcs)
  query <- AddMetaData(query, metadata = predictions[,"predicted.id",drop=FALSE])
  #table(query$predicted.id)
  
  ## UMAP projection  
  ref <- Seurat::RunUMAP(ref, dims = 1:npcs, n.neighbors=nn,
                         reduction = "pca", return.model = TRUE)
  n.anchors <- nrow(anchors@anchors)
  n.anchors
  query <- Seurat::MapQuery(
    anchorset = anchors, reference = ref, query = query,
    refdata = list( label = labels ),
    reference.reduction = "pca", reduction.model = "umap",
    transferdata.args = list( k.weight = min(n.anchors/2,50) )
  )
  
  do.plot = FALSE
  if(do.plot) {
    library(ggplot2)
    ref@meta.data$label <- labels
    p1 <- DimPlot(ref, reduction = "umap", group.by = "label", label = TRUE,
                  label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
    p2 <- DimPlot(query, reduction = "ref.umap", group.by = "predicted.label",
                  label = TRUE, label.size = 3, pt.size=0.3, repel = TRUE) + NoLegend() +
      ggtitle("Query transferred labels")
    p3 <- DimPlot(query, reduction = "umap", group.by = c("predicted.label","stim"),
                  label = TRUE, label.size = 3, pt.size=0.3, repel = TRUE) + NoLegend() +
      ggtitle("Query transferred labels")
    (p1 + p2) / p3 
  }

  ## yes. this is redundant but which one is better?
  predicted <- query@meta.data[,c("predicted.id","predicted.label")] 

  list(ref=ref, query=query, predicted=predicted)
}

#' Transfer labels from reference to query (matrices)  
#'
transferLabels <- function(ref.mat, query.mat, labels) {

  ref.meta <- data.frame(label=labels)
  rownames(ref.meta) <- colnames(ref.mat)
  ref   <-  pgx.justSeuratObject(ref.mat, samples=ref.meta)
  query <-  pgx.justSeuratObject(query.mat, samples=NULL)
  sct=FALSE
  ref <- seurat.preprocess(ref, sct=sct)
  query <- seurat.preprocess(query, sct=sct)  
  tf <- seurat.transferLabels(ref, query, ref$label, sct=sct)
  names(tf)
  
  do.plot=FALSE
  if(do.plot) {
    names(tf$ref@reductions)
    names(tf$query@reductions)  
    ref.pos <- tf$ref@reductions[['umap']]@cell.embeddings
    query.pos <- tf$query@reductions[['ref.umap']]@cell.embeddings
    query.pos2 <- tf$query@reductions[['umap']]@cell.embeddings  

    par(mfrow=c(2,2))
    plot(ref.pos, col = factor(tf$ref$label), main="REFERENCE")
    plot(query.pos, col = factor(tf$query$predicted.label), main="TRANSFERRED (ref.umap)")
    plot(query.pos2, col = factor(tf$query$predicted.label), main="TRANSFERRED (umap)")
  }
   
  list(
    ref = tf$ref,
    query = tf$query,
    predicted = tf$predicted
  )
}


pgx.createSingleCellPGX_DEPRECATED <- function(counts,
                                               samples,
                                               pheno,
                                               batch, 
                                               azimuth.reference = "pbmcref") {

  ## if 'celltype' is not in samples then we do Azimuth
  if("celltype" %in% colnames(samples)) {
    message("[pgx.createSingleCellPGX] using 'celltype' column from sample info")
  } else {
    message("[pgx.createSingleCellPGX] running Azimuth for celltype ...")
    azm <- pgx.runAzimuth(counts, reference = azimuth.reference)
    colnames(azm)
    azm <- azm[,grep("predicted",colnames(azm))]
    ntype <- apply(azm, 2, function(a) length(unique(a)))
    ntype
    ## select smallest level, or highest with at most 10 celltypes
    sel <- ifelse(min(ntype) > 10, which.min(ntype), tail(which(ntype <= 10),1))
    sel
    samples$celltype <- azm[,sel]
    table(samples$celltype)
  }

  dim(counts)
  sc.membership <- NULL
  if(ncol(counts) > 2000) {
    group <- paste0(samples[,"celltype"],":",samples[,pheno])
    if(!is.null(batch)) {
      group <- paste0(group,":",samples[,batch])
    }
    table(group)
    q10 <- quantile(table(group), probs=0.25)
    nb <- round( ncol(counts) / 2000 )
    nb <- ceiling(round( q10 / 20 ))
    nb
    message("[pgx.createSingleCellPGX] running SuperCell. nb = ", nb)    
    sc <- pgx.supercell(counts, samples, group = group, gamma = nb)
    message("[pgx.createSingleCellPGX] SuperCell: ", ncol(counts)," -> ",ncol(sc$counts))    
    counts <- sc$counts
    samples <- sc$meta
    sc.membership <- sc$membership
    dim(counts)
    remove(sc)
  }

  ## Create full Seurat object. Optionally integrate by batch.
  table(samples$celltype)
  batch.vec <- NULL
  message("[pgx.createSingleCellPGX] Creating Seurat object ...")
  if(!is.null(batch)) {
    message("[pgx.createSingleCellPGX] Integrating by batch = ", batch)     
    batch.vec <- samples[,batch]
  }
  obj <- pgx.createSeuratObject(counts, samples, batch = batch.vec,
                                filter=TRUE, method="Harmony") 

  message("[pgx.createSingleCellPGX] Addding Seurat clustering ...")
  r <- "pca"
  names(obj@reductions)
  if(!is.null(batch)) r <- "integrated.dr"  
  obj <- Seurat::RunTSNE(obj, dims=1:30, reduction = r, verbose = FALSE)
  obj <- Seurat::RunTSNE(obj, dim.embed = 3L, dims=1:30, reduction = r,
                         reduction.name ="tsne.3d", reduction.key ="tsne3d_",
                         verbose = FALSE)
  obj <- Seurat::RunUMAP(obj, n.components = 3L, dims=1:30, reduction = r,
                         reduction.name ="umap.3d", reduction.key = "umap3d_", 
                         verbose = FALSE)
  names(obj@reductions)

  ## create balanced down-sampled object. We target about n=20 cells
  ## per statistical condition, per celltype.
  message("[pgx.createSingleCellPGX] Down-sampling Seurat object ...")  
  meta <- obj@meta.data
  group <- paste0(meta$celltype, ":", meta[,pheno])
  table(group)
  length(table(group))
  sub <- seurat.downsample(obj, target_g = 20, group = group) 
  dim(sub)
  table(sub$downsample.group)
  dim(sub)
    
  do.plot = FALSE
  if(do.plot) {

    Seurat::DimPlot(obj, group.by = c("celltype",pheno))
    Seurat::DimPlot(sub, group.by = c("celltype",pheno))  
    names(sub@reductions)

    ## 2D plot
    dim(pos.full)
    pos.full <- obj@reductions[['tsne']]@cell.embeddings
    pos.sub  <- sub@reductions[['tsne']]@cell.embeddings      
    plot(pos.full, pch=20, cex=0.5, col='grey90')
    cc <- factor(sub$celltype)
    points(pos.sub, pch=20, cex=0.8, col=cc)  
    
    ## 3D plot
    pos3d <- obj@reductions[['tsne.3d']]@cell.embeddings
    dim(pos3d)
    pos3d <- uscale(pos3d)  
    colnames(pos3d) <- c("x","y","z")
    pos3d <- data.frame(pos3d)
    pos3d$celltype <- obj@meta.data$celltype
    pos3d$pheno <- obj@meta.data[,pheno]
    fig <- plotly::plot_ly(pos3d, x = ~x, y = ~y, z = ~z, color = ~celltype)
    fig %>% plotly::add_markers( marker = list(size=3) )
  }

  ## results for pgxCreate
  message("[pgx.createSingleCellPGX] Creating PGX object ...")
  counts = sub[['RNA']]$counts
  samples = sub@meta.data
  df <- samples[,c(pheno,"celltype")]
  contrasts <- pgx.makeAutoContrastsStratified(
    df, strata.var = "celltype", mingrp = 3, max.level = 99,
    ref = NULL, slen = 20, fix.degenerate = FALSE, skip.hidden = TRUE) 
  colnames(contrasts) <- sub(".*@","",colnames(contrasts))
  colnames(contrasts) <- gsub("[ ]","_",colnames(contrasts))

  ## single-cell specific normalization (10k)
  X <- logCPM( counts, total = 1e4, prior=1 )
  
  pgx <- pgx.createPGX(
    counts,
    samples,
    contrasts,
    organism = "Human",
    custom.geneset = NULL,
    annot_table = NULL,
    max.genesets = 5000,
    name = "Data set",
    datatype = "scRNA-seq",
    probe_type = NULL,
    creator = "unknown",
    description = "No description provided.",
    X = X,
    impX = NULL,
    norm_method = "CPM",
    is.logx = FALSE,
    batch.correct = FALSE,
    auto.scale = TRUE,
    filter.genes = TRUE,
    prune.samples = FALSE,
    only.known = TRUE,
    only.hugo = TRUE,
    convert.hugo = TRUE,
    only.proteincoding = TRUE,
    remove.xxl = TRUE,
    remove.outliers = TRUE) 

  dim(counts)
  dim(pgx$X)

  ## We take the clusterings from Seurat because these are
  ## 'integrated' (batch corrected).
  cluster = list(
    pca2d = sub@reductions[['pca']]@cell.embeddings[,1:2],
    pca3d = sub@reductions[['pca']]@cell.embeddings[,1:3],
    tsne2d = sub@reductions[['tsne']]@cell.embeddings,
    tsne3d = sub@reductions[['tsne.3d']]@cell.embeddings,
    umap2d = sub@reductions[['umap']]@cell.embeddings,
    umap3d = sub@reductions[['umap.3d']]@cell.embeddings            
  )

  cluster.full = list(
    pca2d = obj@reductions[['pca']]@cell.embeddings[,1:2],
    pca3d = obj@reductions[['pca']]@cell.embeddings[,1:3],
    tsne2d = obj@reductions[['tsne']]@cell.embeddings,
    tsne3d = obj@reductions[['tsne.3d']]@cell.embeddings,
    umap2d = obj@reductions[['umap']]@cell.embeddings,
    umap3d = obj@reductions[['umap.3d']]@cell.embeddings            
  )

  pgx$cluster$pos <- cluster
  pgx$cluster$pos.full <- cluster.full    ## 'full' set
  dim(pgx$cluster$pos[[1]])
  dim(pgx$cluster$pos.full[[1]])
  
  message("[pgx.createSingleCellPGX] done!")  
  return(pgx)
}


#' @export
pgx.createSingleCellPGX <- function(counts,
                                    samples,
                                    contrasts,
                                    organism,
                                    scrnaseq_pheno,
                                    batch = NULL, 
                                    azimuth_ref,
                                    sc_pheno) {

  message("[pgx.createSingleCellPGX]==========================================")
  message("[pgx.createSingleCellPGX]======= pgx.createSingleCellPGX ==========")
  message("[pgx.createSingleCellPGX]==========================================")

  if (!is.null(counts)) {
    message("[createSingleCellPGX] dim.counts: ", dim(counts)[1], ",", dim(counts)[2])
    message("[createSingleCellPGX] class.counts: ", class(counts))
  } else {
    stop("[createPGX] FATAL: counts must be provided")
  }

  if (is.null(samples)) {
    stop("[createSingleCellPGX] FATAL: counts must be provided")
  }
  
  if (!all(colnames(counts) == rownames(samples))) {
    stop("colnames of counts and rownames of samples do not match\n")
  }

  if (is.null(organism)) {
    stop("[createSingleCellPGX] FATAL: organism must be provided")
  }

  if (is.null(azimuth_ref)) {
    azimuth_ref <- "pbmcref"
    message("[createSingleCellPGX] Ref. atlas for Azimuth not specified. Default to pbmcref.")
  } else {
    message("[createSingleCellPGX] Azimuth ref. atlas: ", azimuth_ref)
  }

  if (is.null(sc_pheno)) {
    stop("[createSingleCellPGX] FATAL: phenotype of interest must be provided")
  } else {
    pheno <- sc_pheno
    if (!pheno %in% colnames(samples)) {
      stop("[createSingleCellPGX] FATAL: phenotype of interest not found in sample file.")
    }
  }
  
  cl <- c("celltype", "cell_type", "cell.type", "CELLTYPE", "CELL_TYPE", "CellType")
  cl <- unique(c(cl, toupper(cl), tolower(cl)))
  kk <- intersect(cl, colnames(samples))
  if (length(kk)>0) {
    message("[pgx.createSingleCellPGX] using 'celltype' column from sample info")
    colnames(samples)[match(kk, colnames(samples))] <- "celltype"
  } else {
    message("[pgx.createSingleCellPGX] Inferring cell types with Azimuth.")
    message("[pgx.createSingleCellPGX] Using ", azimuth_ref, " as reference atlas.")
    azm <- pgx.runAzimuth(counts, reference = azimuth_ref)
    azm <- azm[, grep("predicted",colnames(azm))]
    ntype <- apply(azm, 2, function(a) length(unique(a)))
    sel <- ifelse(min(ntype) > 10, which.min(ntype), tail(which(ntype <= 10),1))
    samples$celltype <- azm[, sel]
  }

  sc.membership <- NULL
  if(ncol(counts) > 2000000000) {
    group <- paste0(samples[, "celltype"], ":", samples[, pheno])
    if(!is.null(batch)) {
      group <- paste0(group, ":", samples[, batch])
    }
    q10 <- quantile(table(group), probs=0.25)
    nb <- round( ncol(counts) / 2000 )
    ## nb <- ceiling(round( q10 / 20 ))
    message("[pgx.createSingleCellPGX]=======================================")
    message("[pgx.createSingleCellPGX] running SuperCell. nb = ", nb)    
    sc <- pgx.supercell(counts, samples, group = group, gamma = nb)
    message("[pgx.createSingleCellPGX] SuperCell done: ", ncol(counts), " -> ", ncol(sc$counts))
    message("[pgx.createSingleCellPGX]=======================================")
    counts <- sc$counts
    samples <- sc$meta
    sc.membership <- sc$membership
    remove(sc)
  }

  ## Create full Seurat object. Optionally integrate by batch.
  batch.vec <- NULL
  message("[pgx.createSingleCellPGX] Creating Seurat object ...")
  if (!is.null(batch)) {
    message("[pgx.createSingleCellPGX] Integrating by batch = ", batch)     
    batch.vec <- samples[, batch, drop = FALSE]
  }
  obj <- pgx.createSeuratObject(
    counts,
    samples,
    batch = batch.vec,
    filter = TRUE,
    method = "Harmony"
  ) 

  message("[pgx.createSingleCellPGX] Perform Seurat PCA, t-SNE, UMAP.")
  r <- "pca"
  if (!is.null(batch)) { r <- "integrated.dr" }
  obj <- Seurat::RunTSNE(obj, dims=1:30, reduction = r, verbose = FALSE)
  obj <- Seurat::RunTSNE(obj, dim.embed = 3L, dims=1:30,
    reduction = r, reduction.name ="tsne.3d",
    reduction.key ="tsne3d_", verbose = FALSE)
  obj <- Seurat::RunUMAP(obj, n.components = 3L, dims=1:30,
    reduction = r, reduction.name ="umap.3d",
    reduction.key = "umap3d_",  verbose = FALSE)
  ##  names(obj@reductions)

  message("[pgx.createSingleCellPGX] dim(seurat counts): ", paste0(dim(obj@assays$RNA$counts), collapse=","))
  message("[pgx.createSingleCellPGX] dim(seurat X): ", paste0(dim(obj@assays$RNA$data), collapse=","))
  message("[pgx.createSingleCellPGX] dim(meta.data): ", paste0(dim(obj@meta.data), collapse=","))
  
  ## create balanced down-sampled object. We target about n=20 cells
  ## per statistical condition, per celltype.
  downsample = FALSE
  if (downsample) {
    message("[pgx.createSingleCellPGX] Down-sampling Seurat object ...")  
    group <- paste0(obj@meta.data$celltype, ":", obj@meta.data[, pheno])
    sub <- seurat.downsample(obj, target_g = 20, group = group)
  } else {
    sub <- obj
  }
  
  ## results for pgxCreate
  message("[pgx.createSingleCellPGX] Creating PGX object ...")
  counts = sub[['RNA']]$counts
  samples = sub@meta.data
  df <- samples[,c(pheno,"celltype")]
  contrasts <- pgx.makeAutoContrastsStratified(
    df, strata.var = "celltype", mingrp = 3, max.level = 99,
    ref = NULL, slen = 20, fix.degenerate = FALSE, skip.hidden = TRUE) 
  colnames(contrasts) <- sub(".*@","",colnames(contrasts))
  colnames(contrasts) <- gsub("[ ]","_",colnames(contrasts))

  ## single-cell specific normalization (10k)
  X <- playbase::logCPM(counts, total = 1e4, prior = 1)

  message("[pgx.createSingleCellPGX] dim(counts): ", paste0(dim(counts), collapse=","))
  message("[pgx.createSingleCellPGX] dim(X): ", paste0(dim(X), collapse=","))
  message("[pgx.createSingleCellPGX] dim(samples): ", paste0(dim(samples), collapse=","))

  pgx <- pgx.createPGX(
    counts = counts,
    samples = samples,
    contrasts = contrasts,
    organism = organism,
    custom.geneset = NULL,
    annot_table = NULL,
    max.genesets = 5000,
    name = "Data set",
    datatype = "scRNAseq", ## hack from scRNA-seq.
    azimuth_ref = azimuth_ref, ## NEW AZ
    probe_type = NULL,
    creator = "unknown",
    description = paste0(azimuth_ref, "_scRNAseq_dataset"), ## "No description provided.",
    X = X,
    impX = NULL,
    norm_method = "CPM",
    is.logx = FALSE,
    batch.correct = FALSE,
    auto.scale = TRUE,
    filter.genes = TRUE,
    prune.samples = FALSE,
    only.known = TRUE,
    only.hugo = TRUE,
    convert.hugo = TRUE,
    only.proteincoding = TRUE,
    remove.xxl = TRUE,
    remove.outliers = TRUE
  ) 

  message("[pgx.createSingleCellPGX] dim(pgx$counts): ", paste0(dim(pgx$counts), collapse=","))
  message("[pgx.createSingleCellPGX] dim(pgx$X): ", paste0(dim(pgx$X), collapse=","))
  message("[pgx.createSingleCellPGX] dim(pgx$samples): ", paste0(dim(pgx$samples), collapse=","))

  ## We take the clusterings from Seurat because these are
  ## 'integrated' (batch corrected).
  cluster = list(
    pca2d = sub@reductions[['pca']]@cell.embeddings[,1:2],
    pca3d = sub@reductions[['pca']]@cell.embeddings[,1:3],
    tsne2d = sub@reductions[['tsne']]@cell.embeddings,
    tsne3d = sub@reductions[['tsne.3d']]@cell.embeddings,
    umap2d = sub@reductions[['umap']]@cell.embeddings,
    umap3d = sub@reductions[['umap.3d']]@cell.embeddings            
  )

  cluster.full = list(
    pca2d = obj@reductions[['pca']]@cell.embeddings[,1:2],
    pca3d = obj@reductions[['pca']]@cell.embeddings[,1:3],
    tsne2d = obj@reductions[['tsne']]@cell.embeddings,
    tsne3d = obj@reductions[['tsne.3d']]@cell.embeddings,
    umap2d = obj@reductions[['umap']]@cell.embeddings,
    umap3d = obj@reductions[['umap.3d']]@cell.embeddings            
  )

  pgx$cluster$pos <- cluster
  pgx$cluster$pos.full <- cluster.full
  dim(pgx$cluster$pos[[1]])
  dim(pgx$cluster$pos.full[[1]])

  message("\n\n")
  message("[pgx.createSingleCellPGX]==========================================")
  message("[pgx.createSingleCellPGX]===== pgx.createSingleCellPGX: DONE! =====")
  message("[pgx.createSingleCellPGX]==========================================")
  message("\n\n")
  return(pgx)

}



## =====================================================================================
## =========================== END OF FILE =============================================
## =====================================================================================
