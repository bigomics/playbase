## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.

#' @title Convert Seurat to PGX
#' @param obj Seurat object to convert
#' @param do.cluster Logical indicating whether to cluster samples. Default is FALSE.
#' @return PGX object
#' @description Converts a Seurat single-cell RNA-seq object into a PGX object
#' @details This function takes a Seurat object containing single-cell RNA-seq data and converts it into a PGX object.
#' The count matrix, normalized expression matrix, and sample metadata are extracted from the Seurat object.
#' Gene annotations are added using the gene symbols.
#' If do.cluster=TRUE, dimensionality reduction and clustering of samples is performed.
#' Any existing tsne/umap embeddings and cluster assignments are copied over from the Seurat object.
#' @export
seurat2pgx <- function(obj, do.cluster = FALSE, organism) {

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

  return(pgx)

}

#' @title Integrate single-cell data across batches
#' @param X Numeric matrix of expression values, cells as columns
#' @param batch Factor specifying batch for each cell
#' @param method Methods to use for batch integration. Options are "ComBat", "limma", "CCA", "MNN", "Harmony", "liger"
#' @return List containing batch-integrated expression matrices by method
#' @description Integrate single-cell RNA-seq data from multiple batches using various batch correction methods.
#' @details This function takes a single-cell expression matrix \code{X} and a batch vector \code{batch} as input.
#' It applies different batch correction methods to integrate the data across batches.
#' The following methods can be selected via the \code{method} parameter:
#' \itemize{
#' \item ComBat: Apply ComBat batch correction from sva package
#' \item limma: Apply removeBatchEffect from limma package
#' \item CCA: Apply canonical correlation analysis using Seurat package
#' \item MNN: Apply mutual nearest neighbors correction
#' \item Harmony: Apply Harmony integration using Harmony package
#' \item liger: Apply integration via liger package
#' }
#' The batch-integrated expression matrices are returned as a list by method name.
#' @export
pgx.scBatchIntegrate <- function(X,
                                 batch,
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
    ## mnn <- try(batchelor::mnnCorrect(X, batch = batch, cos.norm.in = TRUE, cos.norm.out = FALSE))
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
    ## liger@var.genes <- Matrix::head(rownames(X)[order(-apply(X, 1, stats::sd))], 100)
    liger@var.genes <- Matrix::head(rownames(X)[order(-matrixStats::rowSds(X, na.rm = TRUE))], 100)
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
#' @param counts Single-cell count matrix
#' @param batch Batch vector assigning batches to cells
#' @param qc.filter Logical indicating whether to filter cells by QC metrics. Default is FALSE.
#' @param nanchors Number of anchor points to use for integration. Default is -1 to auto-determine.
#' @param sct Logical indicating whether to use SCTransform normalization. Default is FALSE.
#' @return Seurat object containing integrated data
#' @description Integrate single-cell RNA-seq data from multiple batches using canonical correlation analysis via the Seurat package.
#' @details This function takes a single-cell count matrix \code{counts} and a \code{batch} vector as input.
#' It sets up a Seurat object for each batch and integrates them using FindIntegrationAnchors/IntegrateData functions.
#' Cells can be filtered by QC metrics like mitochondrial content if \code{qc.filter=TRUE}.
#' The \code{nanchors} parameter controls the number of anchor points used for integration.
#' Normalization and scaling can be done using SCTransform if \code{sct=TRUE}.
#' The integrated Seurat object is returned containing the corrected expression matrix.
#' @export
pgx.SeuratBatchIntegrate <- function(counts,
                                     batch,
                                     qc.filter = FALSE,
                                     nanchors = -1,
                                     sct = FALSE) {

  ## From Seurat vignette: Integration/batch correction using
  ## CCA. Note there is no QC filtering for samples on ribo/mito
  ## content. You need to do that before.

  nbatch <- length(unique(batch))
  message("[pgx.SeuratBatchIntegrate] Processing ", nbatch, " batches...")
  obj.list <- list()
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
    anchor.features <- nrow(counts)
    if (nanchors > 0) anchor.features <- nanchors
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

  counts

}

#' @title SuperCell down sampling. Uniform down samplsing using gamma = 20.
#' @export
pgx.supercell <- function(counts,
                          meta,
                          group = NULL,
                          gamma = 20,
                          nvargenes = 1000,
                          log.transform = TRUE) {

  if (log.transform) { ## supercell uses log2 matrix
    X <- logCPM(counts, total = 1e4)
  } else {
    X <- counts
  }

  if (is.null(group) && "group" %in% colnames(meta)) {
    message("using group column detected in meta\n")
    group <- meta[, "group"]
  }

  if (!is.null(group) && any(group %in% colnames(meta))) {
    group <- intersect(group, colnames(meta))
    message("using groups: ", paste(group, collapse = "."))
    group <- meta[, group]
    if (NCOL(group) > 1) group <- apply(group, 1, paste, collapse = ".")
  }

  SC <- SuperCell::SCimplify(X,
    gamma = gamma,
    n.var.genes = nvargenes,
    cell.split.condition = group
  )
  message("[pgx.supercell] SuperCell::SCimplify completed")

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
  sc.meta <- sc.meta[, ii, drop = FALSE]

  ## Compute metacall expression as sum of counts
  counts <- as.matrix(counts)
  if (log.transform) {
    sc.counts <- SuperCell::supercell_GE(
      counts,
      mode = "sum",
      groups = SC$membership
    )
  } else {
    sc.counts <- SuperCell::supercell_GE(
      counts,
      mode = "average",
      groups = SC$membership
    )
  }

  message("[pgx.supercell] SuperCell::supercell_GE completed")
  sc.membership <- paste0("mc", SC$membership)
  colnames(sc.counts) <- paste0("mc", 1:ncol(sc.counts))
  rownames(sc.meta) <- colnames(sc.counts)

  list(counts = sc.counts, meta = sc.meta, membership = sc.membership)

}


#' @title SuperCell downsampling by group targeting equal sizes per group.
#' @export
pgx.supercell2 <- function(counts, meta, group, target_n = 20) {
  ## metacell equal across celltype and phenotype
  g <- group[1]
  sc.membership <- rep(NA, length(group))
  sc.counts <- c()
  sc.meta <- c()
  for (g in unique(group)) {
    sel <- which(group == g)
    k <- max(1, length(sel) / target_n)
    sc <- pgx.supercell(counts[, sel], samples[sel, ], gamma = k)
    colnames(sc$counts) <- paste0(g, ".", colnames(sc$counts))
    rownames(sc$meta) <- colnames(sc$counts)
    sc$meta$group <- g
    sc$membership <- paste0(g, ".", sc$membership)
    sc.membership[sel] <- sc$membership
    sc.counts <- cbind(sc.counts, sc$counts)
    sc.meta <- rbind(sc.meta, sc$meta)
  }

  list(counts = sc.counts, meta = sc.meta, membership = sc.membership)

}

## pgx.sc_anchors <- function(counts,
##                            sc.counts, sc.membership,
##                            sc.group = colnames(sc.counts)) {

##   ## determine closest sample to metacell reference (aka anchor)
##   sc.ref <- rep(NA, ncol(sc.counts))
##   X1 <- logCPM(counts, 1e4)
##   X2 <- as.matrix(logCPM(sc.counts, 1e4))
##   for (i in 1:ncol(sc.counts)) {
##     mc <- sc.group[i]
##     sel <- which(sc.membership == mc)
##     x1 <- as.matrix(X1[, sel, drop = FALSE])
##     colnames(x1) <- colnames(X1)[sel]
##     rmax <- which.max(cor(x1, X2[, i])[, 1])
##     sc.ref[i] <- colnames(x1)[rmax]
##   }

##   return(sc.ref)

## }

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
                                   sc_compute_settings = list(),
                                   preprocess = TRUE,
                                   method = "Harmony") {

  message("[pgx.createSingleCellPGX] Creating Seurat object ...")

  options(Seurat.object.assay.calcn = TRUE)
  getOption("Seurat.object.assay.calcn")

  if (is.null(samples)) {
    samples <- data.frame(row.names = colnames(counts))
  }

  samples <- as.data.frame(samples)
  rownames(samples) <- colnames(counts)
  if (!is.null(batch)) samples$batch <- batch
  samples$nCount_RNA <- Matrix::colSums(counts, na.rm = TRUE)
  samples$nFeature_RNA <- Matrix::colSums(counts > 0, na.rm = TRUE)
  obj <- Seurat::CreateSeuratObject(counts = counts, meta.data = samples)
  obj <- Seurat::PercentageFeatureSet(obj, pattern = "^MT-|^Mt-", col.name = "percent.mt")
  obj <- Seurat::PercentageFeatureSet(obj, pattern = "^RP[LS]|^Rp[ls]", col.name = "percent.ribo")
  obj <- Seurat::PercentageFeatureSet(obj, pattern = "^HB|^Hb", col.name = "percent.hb")
  obj@meta.data$percent.mt <- obj@meta.data$percent.mt + 1e-5
  obj@meta.data$percent.ribo <- obj@meta.data$percent.ribo + 1e-5
  obj@meta.data$percent.hb <- obj@meta.data$percent.hb + 1e-5

  if (cellcyclescores) { ## needs normalized data
    counts <- obj@assays$RNA$counts
    totcounts <- Matrix::colSums(counts, na.rm = TRUE)
    cpm <- sweep(counts, 2, totcounts, FUN = "/") * 1e4
    obj@assays$RNA$data <- log1p(cpm)
    g2m.ff <- Seurat::cc.genes$g2m.genes
    s.ff <- Seurat::cc.genes$s.genes
    obj <- Seurat::CellCycleScoring(obj, s.features = s.ff, g2m.features = g2m.ff)
  }

  ncells0 <- ncol(obj)

  if (filter) {
    message("[pgx.createSeuratObject] Filtering cells")
    flts <- unlist(lapply(sc_compute_settings, function(x) isTRUE(x)))

    if (any(flts)) {
      flts <- names(flts)[which(flts)]

      nfeat_thr <- c(0, 1e200)
      if ("nfeature_threshold" %in% flts) {
        nfeat_thr <- sc_compute_settings[["nfeature_threshold"]]
      }

      mt_thr <- 100
      if ("mt_threshold" %in% flts) {
        mt_thr <- sc_compute_settings[["mt_threshold"]]
      }

      hb_thr <- 100
      if ("hb_threshold" %in% flts) {
        hb_thr <- sc_compute_settings[["hb_threshold"]]
      }

      obj <- subset(obj,
        subset =
          nFeature_RNA > nfeat_thr[1] &
            nFeature_RNA < nfeat_thr[2] &
            percent.mt < mt_thr &
            percent.hb < hb_thr
      )
      ncells1 <- ncol(obj)
      message("[pgx.createSeuratObject] Filtering cells: ", ncells0, " --> ", ncells1)
    } else {
      message("[pgx.createSeuratObject]: No user-selected filters.")
      ncells1 <- ncol(obj)
      message("[pgx.createSeuratObject] N. of cells: ", ncells0, " --> ", ncells1)
    }
  }

  if (preprocess) {
    message("[pgx.createIntegratedSeuratObject] Preprocessing Seurat object")
    has.batch <- "batch" %in% tolower(colnames(obj@meta.data))
    if (!is.null(batch)) {
      obj <- seurat.integrate(obj, "batch", sct = TRUE, method = "Harmony")
    } else {
      obj <- seurat.preprocess(obj, sc_compute_settings, sct = TRUE)
    }
  }

  return(obj)

}

#' @export
seurat.preprocess <- function(obj,
                              sc_compute_settings = list(),
                              sct = FALSE,
                              tsne = TRUE,
                              umap = TRUE) {

  options(future.globals.maxSize = 4 * 1024^4)

  vars.to.regress <- NULL

  if (length(sc_compute_settings)) {
    if (sc_compute_settings[["regress_mt"]]) {
      vars.to.regress <- c(vars.to.regress, "percent.mt")
    }
    if (sc_compute_settings[["regress_hb"]]) {
      vars.to.regress <- c(vars.to.regress, "percent.hb")
    }
    if (sc_compute_settings[["regress_ribo"]]) {
      vars.to.regress <- c(vars.to.regress, "percent.ribo")
    }
    if (sc_compute_settings[["regress_ccs"]]) {
      vars.to.regress <- c(vars.to.regress, "S.Score", "G2M.Score")
    }
  }

  if (sct) {
    message("[seurat.preprocess] Performing Seurat SCT normalization")
    obj <- Seurat::SCTransform(obj,
      method = "glmGamPoi",
      vars.to.regress = vars.to.regress, verbose = FALSE
    )
  } else {
    message("[seurat.preprocess] Performing Seurat standard Normalization, findVarFeatures, scaling")
    ## obj <- Seurat::NormalizeData(obj) causes "integer overflow" issue
    counts <- obj@assays$RNA$counts
    totcounts <- Matrix::colSums(counts, na.rm = TRUE)
    cpm <- sweep(counts, 2, totcounts, FUN = "/") * 1e4
    obj@assays$RNA$data <- log1p(cpm)
    obj <- Seurat::FindVariableFeatures(obj)
    if (!is.null(vars.to.regress)) {
      message("[seurat.preprocess] Regressing out ", paste0(vars.to.regress, collapse = "; "))
    }
    obj <- Seurat::ScaleData(obj, vars.to.regress = vars.to.regress)
  }

  message("[seurat.preprocess] Performing Seurat PCA")
  npcs <- min(30, ncol(obj) / 2)
  obj <- Seurat::RunPCA(obj, npcs = npcs, verbose = FALSE)

  # clustering using integrated data or original pca
  dr <- "pca"
  nn <- min(30L, ncol(obj) / 5)
  message("[seurat.preprocess] Performing Seurat FindNeighbors & Clusters")
  obj <- Seurat::FindNeighbors(obj, dims = 1:npcs, reduction = dr, verbose = FALSE)
  obj <- Seurat::FindClusters(obj, resolution = 1, verbose = FALSE)

  if (tsne) {
    message("[seurat.preprocess] Performing Seurat 2D tSNE")
    obj <- Seurat::RunTSNE(obj, dims = 1:npcs, reduction = dr, verbose = FALSE)
  }

  if (umap) {
    message("[seurat.preprocess] Performing Seurat 2D UMAP")
    obj <- Seurat::RunUMAP(obj,
      dims = 1:npcs,
      n.neighbors = nn, reduction = dr, verbose = FALSE
    )
  }

  return(obj)

}

#' @export
seurat.integrate <- function(obj,
                             batch,
                             sct = TRUE,
                             method = "Harmony") {

  obj[["RNA"]] <- split(obj[["RNA"]], f = obj@meta.data[, batch])

  if (sct) {
    message("[seurat.integrate] normalization method = SCT")
    obj <- Seurat::SCTransform(obj, method = "glmGamPoi", verbose = FALSE)
  } else {
    obj <- Seurat::NormalizeData(obj)
    obj <- Seurat::FindVariableFeatures(obj)
    obj <- Seurat::ScaleData(obj)
  }

  message("[seurat.integrate] Performing Seurat PCA")
  obj <- Seurat::RunPCA(obj, npcs = 30, verbose = FALSE)

  sel.method <- Seurat::HarmonyIntegration
  sel.method <- switch(method,
    CCA = Seurat::CCAIntegration,
    RPCA = Seurat::RPCAIntegration,
    Harmony = Seurat::HarmonyIntegration,
    JointPCA = Seurat::JointPCAIntegration
  )

  message("[seurat.integrate] integration method = ", method)
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
  obj <- Seurat::RunUMAP(obj, dims = 1:30, reduction = dr, verbose = FALSE)

  return(obj)

}


#' @export
pgx.runAzimuth <- function(counts, k.weight = NULL, reference = NULL) {

  options(future.globals.maxSize = 4 * 1024^4)

  obj <- pgx.justSeuratObject(counts, samples = NULL)

  if (is.null(k.weight)) k.weight <- 20

  message("[pgx.runAzimuth] k.weight = ", k.weight)

  if (is.null(reference)) reference <- "pbmcref"

  obj1 <- try(Azimuth::RunAzimuth(
    obj,
    reference = reference,
    k.weight = k.weight,
    verbose = FALSE
  ))

  if (!"try-error" %in% class(obj1)) {
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
#' @export
pgx.downsample <- function(counts, meta, group, target_n = 20) {
  ## metacell equal across celltype and phenotype
  g <- group[1]
  sc.membership <- rep(NA, length(group))
  sc.counts <- c()
  sc.meta <- c()
  sel <- tapply(1:ncol(counts), group, function(ii) {
    head(sample(ii), target_n)
  })
  sel <- unlist(sel)
  sc.counts <- counts[, sel]
  sc.meta <- meta[sel, ]
  sc.meta$group <- group[sel]
  list(counts = sc.counts, meta = sc.meta, membership = NULL)
}

seurat.downsample <- function(obj, target_n = 1000, target_g = 20, group = NULL) {
  if (ncol(obj) <= target_n) {
    return(obj)
  }
  if (is.null(group)) {
    sub <- obj[, sample(ncol(obj), target_n)]
    return(sub)
  }
  obj@meta.data$downsample.group <- group
  sel <- tapply(1:ncol(obj), group, function(ii) head(sample(ii), target_g))
  sub <- obj[, unlist(sel)]
  return(sub)
}


#' Transfer labels from reference to query (Seurat objects)
seurat.transferLabels <- function(ref, query, labels, sct = TRUE) {
  ## find anchors
  npcs <- ncol(ref@reductions[[1]])
  nn <- min(30L, ncol(ref) / 5)
  anchors <- Seurat::FindTransferAnchors(
    reference = ref, query = query, dims = 1:npcs, k.score = nn,
    normalization.method = ifelse(sct, "SCT", "LogNormalize"),
    reference.reduction = "pca"
  )

  ## redudant??? MapQuery does the same??
  if (length(labels) == 1 && labels %in% colnames(ref@meta.dat)) {
    refdata <- ref@meta.data[, labels]
  } else {
    refdata <- labels
  }
  predictions <- Seurat::TransferData(anchorset = anchors, refdata = refdata, dims = 1:npcs)
  query <- AddMetaData(query, metadata = predictions[, "predicted.id", drop = FALSE])
  # table(query$predicted.id)

  ## UMAP projection
  ref <- Seurat::RunUMAP(ref,
    dims = 1:npcs, n.neighbors = nn,
    reduction = "pca", return.model = TRUE
  )
  n.anchors <- nrow(anchors@anchors)
  n.anchors
  query <- Seurat::MapQuery(
    anchorset = anchors, reference = ref, query = query,
    refdata = list(label = labels),
    reference.reduction = "pca", reduction.model = "umap",
    transferdata.args = list(k.weight = min(n.anchors / 2, 50))
  )

  do.plot <- FALSE
  if (do.plot) {
    library(ggplot2)
    ref@meta.data$label <- labels
    p1 <- DimPlot(ref,
      reduction = "umap", group.by = "label", label = TRUE,
      label.size = 3, repel = TRUE
    ) + NoLegend() + ggtitle("Reference annotations")
    p2 <- DimPlot(query,
      reduction = "ref.umap", group.by = "predicted.label",
      label = TRUE, label.size = 3, pt.size = 0.3, repel = TRUE
    ) + NoLegend() +
      ggtitle("Query transferred labels")
    p3 <- DimPlot(query,
      reduction = "umap", group.by = c("predicted.label", "stim"),
      label = TRUE, label.size = 3, pt.size = 0.3, repel = TRUE
    ) + NoLegend() +
      ggtitle("Query transferred labels")
    (p1 + p2) / p3
  }

  ## yes. this is redundant but which one is better?
  predicted <- query@meta.data[, c("predicted.id", "predicted.label")]

  list(ref = ref, query = query, predicted = predicted)
}

#' Transfer labels from reference to query (matrices)
transferLabels <- function(ref.mat, query.mat, labels) {
  ref.meta <- data.frame(label = labels)
  rownames(ref.meta) <- colnames(ref.mat)
  ref <- pgx.justSeuratObject(ref.mat, samples = ref.meta)
  query <- pgx.justSeuratObject(query.mat, samples = NULL)
  sct <- FALSE
  ref <- seurat.preprocess(ref, sct = sct)
  query <- seurat.preprocess(query, sct = sct)
  tf <- seurat.transferLabels(ref, query, ref$label, sct = sct)
  do.plot <- FALSE
  if (do.plot) {
    ref.pos <- tf$ref@reductions[["umap"]]@cell.embeddings
    query.pos <- tf$query@reductions[["ref.umap"]]@cell.embeddings
    query.pos2 <- tf$query@reductions[["umap"]]@cell.embeddings
    par(mfrow = c(2, 2))
    plot(ref.pos, col = factor(tf$ref$label), main = "REFERENCE")
    plot(query.pos, col = factor(tf$query$predicted.label), main = "TRANSFERRED (ref.umap)")
    plot(query.pos2, col = factor(tf$query$predicted.label), main = "TRANSFERRED (umap)")
  }
  list(ref = tf$ref, query = tf$query, predicted = tf$predicted)
}


#' @export
pgx.createSingleCellPGX <- function(counts,
                                    samples,
                                    contrasts,
                                    organism,
                                    azimuth_ref,
                                    batch = NULL,
                                    sc_compute_settings = list()) {
  message("[pgx.createSingleCellPGX]==========================================")
  message("[pgx.createSingleCellPGX]======= pgx.createSingleCellPGX ==========")
  message("[pgx.createSingleCellPGX]==========================================")

  if (!is.null(counts)) {
    message("[createSingleCellPGX] dim.counts: ", nrow(counts), ",", ncol(counts))
  } else {
    stop("[createPGX] FATAL: counts must be provided")
  }

  if (is.null(samples)) {
    stop("[createSingleCellPGX] FATAL: samples must be provided")
  }

  if (!all(colnames(counts) == rownames(samples))) {
    stop("colnames of counts and rownames of samples do not match\n")
  }

  if (is.null(organism)) {
    stop("[createSingleCellPGX] FATAL: organism must be provided")
  }

  if (is.null(azimuth_ref)) {
    azimuth_ref <- "pbmcref"
    message("[createSingleCellPGX] Azimuth ref. atlas not specified. Defaulting to pbmcref.")
  } else {
    message("[createSingleCellPGX] Azimuth ref. atlas: ", azimuth_ref)
  }

  sc_params <- names(sc_compute_settings)
  if (!is.null(sc_compute_settings) && length(sc_compute_settings) > 0) {
    info("[createSingleCellPGX] scRNAseq parameters: ", sc_params)
  } else {
    sc_compute_settings <- list()
  }

  message("[pgx.createSingleCellPGX] Azimuth reference atlas: ", azimuth_ref)
  library(Seurat) # DO NOT REMOVE THIS LINE (Azimuth fails if not sourcing Seurat)
  azm <- pgx.runAzimuth(counts, reference = azimuth_ref)
  azm <- azm[, grep("predicted", colnames(azm))]
  ntype <- apply(azm, 2, function(a) length(unique(a)))
  sel <- ifelse(min(ntype) > 10, which.min(ntype), tail(which(ntype <= 10), 1))
  samples <- as.data.frame(samples)
  samples$celltype <- azm[, sel]

  ## because the samples get downsamples, also the sample and
  ## contrasts matrices get downsampled (by majority label).
  samplesx <- cbind(samples, contrasts)

  sc.membership <- NULL
  do.supercells <- sc_compute_settings[["compute_supercells"]]
  if (is.null(do.supercells)) {
    do.supercells <- FALSE
  }

  if (do.supercells || ncol(counts) > 10000) {
    if (do.supercells) {
      message("[pgx.createSingleCellPGX] User choice: performing SuperCell")
    }
    if (ncol(counts) > 10000) {
      message("[pgx.createSingleCellPGX] >10K cells: performing SuperCell")
    }

    ct <- samplesx[, "celltype"]
    group <- paste0(ct, ":", apply(contrasts, 1, paste, collapse = "_"))
    ##  group <- paste0(samples[, "celltype"], ":", samples[, pheno])
    if (!is.null(batch)) group <- paste0(group, ":", samples[, batch])

    q10 <- quantile(table(group), probs = 0.25)
    if (ncol(counts) > 2000) {
      nb <- round(ncol(counts) / 2000)
    } else {
      d <- round(ncol(counts) / 8, 1)
      nb <- round(ncol(counts) / d) ## temporary...
    }
    ## nb <- ceiling(round( q10 / 20 ))
    message("[pgx.createSingleCellPGX] running SuperCell. nb = ", nb)
    sc <- pgx.supercell(counts, samplesx, group = group, gamma = nb)
    message("[pgx.createSingleCellPGX] SuperCell done: ", ncol(counts), " -> ", ncol(sc$counts))
    counts <- sc$counts
    samplesx <- sc$meta
    sc.membership <- sc$membership
    remove(sc)
  }

  ## Create full Seurat object. Optionally integrate by batch.
  batch.vec <- NULL
  if (!is.null(batch)) {
    message("[pgx.createSingleCellPGX] Integrating by batch = ", batch)
    batch.vec <- as.vector(unlist(samplesx[, batch]))
  }

  message("[pgx.createSingleCellPGX] Creating Seurat object ...")
  obj <- pgx.createSeuratObject(
    counts = counts,
    samples = samplesx,
    batch = batch.vec,
    cellcyclescores = TRUE,
    filter = TRUE,
    sc_compute_settings = sc_compute_settings,
    preprocess = TRUE,
    method = "Harmony"
  )

  r <- "pca"
  if (!is.null(batch)) r <- "integrated.dr"

  message("[pgx.createSingleCellPGX] Performing Seurat 2D and 3D t-SNE")
  obj <- Seurat::RunTSNE(obj, dims = 1:30, reduction = r, verbose = FALSE)
  obj <- Seurat::RunTSNE(obj,
    dim.embed = 3L, dims = 1:30,
    reduction = r, reduction.name = "tsne.3d",
    reduction.key = "tsne3d_", verbose = FALSE
  )

  message("[pgx.createSingleCellPGX] Performing Seurat 3D UMAP")
  obj <- Seurat::RunUMAP(obj,
    n.components = 3L, dims = 1:30,
    reduction = r, reduction.name = "umap.3d",
    reduction.key = "umap3d_", verbose = FALSE
  )

  ## Create balanced down-sampled object.
  ## Target about n=20 cells per statistical condition, per celltype.
  downsample <- FALSE
  if (!downsample) sub <- obj

  ## results for pgx.createPGX
  message("[pgx.createSingleCellPGX] Creating PGX object ...")
  counts2 <- sub[["RNA"]]$counts
  kk <- setdiff(colnames(sub@meta.data), colnames(contrasts))
  samples2 <- sub@meta.data[, kk]

  ## stratify contrast matrix by celltype
  contrasts2 <- sub@meta.data[, colnames(contrasts), drop = FALSE]
  contrasts2 <- stratifyContrasts(contrasts2, samples2$celltype)
  colnames(contrasts2) <- sub(".*@", "", colnames(contrasts2))
  colnames(contrasts2) <- gsub("[ ]", "_", colnames(contrasts2))

  ## single-cell specific normalization (10k)
  X <- logCPM(counts2, total = 1e4, prior = 1)

  pgx <- pgx.createPGX(
    counts = counts2,
    samples = samples2,
    contrasts = contrasts2,
    organism = organism,
    custom.geneset = NULL,
    annot_table = NULL,
    max.genesets = 5000,
    name = "Data set",
    datatype = "scRNAseq", ## hack from scRNA-seq.
    azimuth_ref = azimuth_ref,
    probe_type = NULL,
    creator = "unknown",
    description = "No description provided.",
    X = X,
    norm_method = "CPM",
    is.logx = FALSE,
    # batch.correct = FALSE,
    batch.correct.method = "no_batch_correct", ## new
    batch.pars = "<autodetect>", ## new
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

  ## We take the clusterings from Seurat because
  ## these are 'integrated' (batch corrected).
  cluster <- list(
    pca2d = sub@reductions[["pca"]]@cell.embeddings[, 1:2],
    pca3d = sub@reductions[["pca"]]@cell.embeddings[, 1:3],
    tsne2d = sub@reductions[["tsne"]]@cell.embeddings,
    tsne3d = sub@reductions[["tsne.3d"]]@cell.embeddings,
    umap2d = sub@reductions[["umap"]]@cell.embeddings,
    umap3d = sub@reductions[["umap.3d"]]@cell.embeddings
  )

  if (!downsample) {
    cluster.full <- cluster
  } else {
    cluster.full <- list(
      pca2d = obj@reductions[["pca"]]@cell.embeddings[, 1:2],
      pca3d = obj@reductions[["pca"]]@cell.embeddings[, 1:3],
      tsne2d = obj@reductions[["tsne"]]@cell.embeddings,
      tsne3d = obj@reductions[["tsne.3d"]]@cell.embeddings,
      umap2d = obj@reductions[["umap"]]@cell.embeddings,
      umap3d = obj@reductions[["umap.3d"]]@cell.embeddings
    )
  }

  pgx$cluster$pos <- cluster
  pgx$cluster$pos.full <- cluster.full

  message("\n\n")
  message("[pgx.createSingleCellPGX]============================================")
  message("[pgx.createSingleCellPGX]======= pgx.createSingleCellPGX: DONE! =====")
  message("[pgx.createSingleCellPGX]============================================")
  message("\n\n")

  return(pgx)

}
## =====================================================================================
## =========================== END OF FILE =============================================
## =====================================================================================
