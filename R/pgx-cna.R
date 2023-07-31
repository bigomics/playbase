##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


#' Infer copy number variations from single cell RNA-seq data
#'
#' @param ngs An NGS object containing single cell RNA-seq data
#' @param refgroup Name of the reference sample group for normalization
#' @param progress Show progress bar?
#'
#' @return Updated NGS object with inferred CNV segments
#'
#' @description Infers copy number variations from single cell RNA-seq data using
#' InferCNV.
#'
#' @details This function takes an NGS object containing single cell RNA-seq
#' data and runs the InferCNV algorithm to infer copy number variations.
#'
#' It extracts the expression matrix and gene annotations from the NGS object.
#' The expression of each gene is normalized to the reference group to estimate
#' relative copy number levels.
#'
#' InferCNV is then applied to segment the genome and identify focal/broad CNV regions.
#' The output CNV segments are added to the NGS object for downstream analysis and plotting.
#'
#' @export
pgx.inferCNV <- function(ngs, refgroup = NULL, progress = NULL) {
  ## InferCNV: Inferring copy number alterations from tumor single
  ## cell RNA-Seq data
  ##
  ## https://github.com/broadinstitute/inferCNV/wiki
  ##

  symbol <- as.vector(as.list(org.Hs.eg.db::org.Hs.egSYMBOL))
  chrloc <- as.list(org.Hs.eg.db::org.Hs.egCHRLOC)
  chr <- as.vector(sapply(chrloc, function(x) names(x)[1]))
  pos <- abs(as.integer(as.vector(sapply(chrloc, function(x) x[1]))))
  chr[sapply(chr, is.null)] <- NA
  chr <- as.character(unlist(chr))
  chr <- chr[match(ngs$genes$gene_name, symbol)]
  pos <- pos[match(ngs$genes$gene_name, symbol)]
  genes <- data.frame(
    chr = paste0("chr", chr),
    start = pos - 1000,
    stop = pos - 1000
  ) ## fake start/stop
  rownames(genes) <- ngs$genes$gene_name

  ## filter known genes
  jj <- which(genes$chr %in% paste0("chr", c(1:22, "X", "Y")) &
    !is.na(genes$start) & !is.na(genes$stop))
  genes <- genes[jj, ]

  ## prepare data objects
  gg <- intersect(rownames(genes), rownames(ngs$counts))
  data <- ngs$counts[gg, ]
  genes <- genes[gg, ]
  annots <- ngs$samples[, "group", drop = FALSE]

  if (FALSE && is.null(refgroup)) {
    ## if no reference group is given, we create a reference by
    ## random sampling of genes.
    ref <- t(apply(data, 1, function(x) sample(x, 50, replace = TRUE)))
    colnames(ref) <- paste0("random.", 1:ncol(ref))
    data <- cbind(data, ref)
    annots <- matrix(c(annots[, 1], rep("random", ncol(ref))), ncol = 1)
    rownames(annots) <- colnames(data)
    colnames(annots) <- "group"
    refgroup <- c("random")
  }

  ## take out tiny groups
  selgrp <- names(which(table(annots[, 1]) >= 2))
  kk <- which(annots[, 1] %in% selgrp)
  data <- data[, kk]
  annots <- annots[colnames(data), , drop = FALSE]

  ## From inferCNV vignette
  infercnv_obj <- infercnv::CreateInfercnvObject(
    raw_counts_matrix = data,
    gene_order_file = genes,
    annotations_file = annots,
    ref_group_names = refgroup
  )

  out_dir <- "/tmp/Rtmpn8rPtL/file19b68b27f09/"
  out_dir <- tempfile()
  cat("DBG pgx.inferCNV:: setting out_dir=", out_dir, "\n")

  infercnv_obj <- infercnv::run(infercnv_obj,
    cutoff = 1,
    out_dir = out_dir,
    cluster_by_groups = TRUE,
    num_threads = 4,
    no_plot = FALSE
  )


  img.file <- paste0(out_dir, "/infercnv.png")
  suppressWarnings(cnv <- data.table::fread(file.path(out_dir, "expr.infercnv.dat"), check.names = FALSE))
  symbol <- cnv[[1]]
  cnv <- as.data.frame(cnv, check.names = FALSE)[2:ncol(cnv)]
  cnv <- as.matrix(cnv)
  rownames(cnv) <- symbol

  genes <- genes[rownames(cnv), ]
  pos <- (genes$start + genes$stop) / 2
  ichr <- as.integer(sub("X", 23, sub("Y", 24, sub("chr", "", genes$chr))))
  jj <- order(ichr, pos)
  pos <- pos[jj]
  chr <- as.character(genes$chr)[jj]
  logcnv <- log2(cnv[jj, ] / mean(cnv, na.rm = TRUE)) ## logarithmic



  img <- png::readPNG(img.file)

  res <- list(cna = logcnv, chr = chr, pos = pos, png = img)

  ## clean up folder??
  unlink(out_dir, recursive = TRUE)

  return(res)
}


#' Estimate copy number from gene expression
#'
#' @param ngs An NGS object containing gene expression data
#' @param nsmooth Smoothing window size for copy number estimation. Default 40.
#'
#' @return A list with estimated copy number values, chromosome, position,
#' and ideogram image for each gene.
#'
#' @description Estimates copy number variation from gene expression data by smoothing
#' relative expression values within genomic windows.
#'
#' @details This function takes an NGS object containing normalized gene expression data.
#' It calculates the relative expression level of each gene compared to the mean expression.
#' These relative levels are then smoothed within sliding genomic windows of size \code{nsmooth}
#' genes to estimate regional copy number variation. The smoothed values are returned along
#' with the chromosome, genomic position, and an ideogram image for data visualization.
#'
#' @export
pgx.CNAfromExpression <- function(ngs, nsmooth = 40) {
  ## This estimates CNV by local smoothing of relative expression
  ## values.
  symbol <- as.vector(as.list(org.Hs.eg.db::org.Hs.egSYMBOL))
  chrloc <- as.list(org.Hs.eg.db::org.Hs.egCHRLOC)
  chr <- as.vector(sapply(chrloc, function(x) names(x)[1]))
  pos <- abs(as.integer(as.vector(sapply(chrloc, function(x) x[1]))))
  chr[sapply(chr, is.null)] <- NA
  chr <- as.character(unlist(chr))
  chr <- chr[match(ngs$genes$gene_name, symbol)]
  pos <- pos[match(ngs$genes$gene_name, symbol)]
  genes <- data.frame(chr = chr, pos = pos)
  rownames(genes) <- ngs$genes$gene_name

  sel <- which(!is.na(genes$chr) & !is.na(genes$pos))
  genes <- genes[sel, ]

  if (!is.null(ngs$counts)) {
    cna <- log2(100 + edgeR::cpm(ngs$counts)) ## moderated log2
  } else {
    cna <- ngs$X
  }
  gg <- intersect(rownames(genes), rownames(cna))
  cna <- cna[gg, ]
  genes <- genes[gg, ]

  ## ---------------------------------------------------------------------
  ## order genes and matrix according genomic position
  ## ---------------------------------------------------------------------
  jj <- which(genes$chr %in% c(1:22, "X", "Y"))
  genes <- genes[jj, ]
  genes$chr <- factor(genes$chr, levels = c(1:22, "X", "Y"))
  jj <- order(genes$chr, genes$pos)
  genes <- genes[jj, ]
  cna <- cna[rownames(genes), ]
  cna0 <- cna

  ## ---------------------------------------------------------------------
  ## apply 'crude' moving average filter (THIS SHOULD BE IMPROVED!)
  ## ---------------------------------------------------------------------
  mavg <- function(x, n = nsmooth) {
    stats::filter(x, rep(1 / n, n), sides = 2, circular = TRUE)
  }
  cna <- t(scale(t(cna), center = FALSE)) ## z-score
  cna <- apply(cna, 2, mavg)
  cna <- cna - apply(cna, 1, median, na.rm = TRUE)
  rownames(cna) <- rownames(cna0)

  res <- list(cna = cna, chr = genes$chr, pos = genes$pos)
  return(res)
}


#' Plot copy number alteration heatmap
#'
#' @title Plot CNA heatmap
#'
#' @param ngs An NGS object containing copy number data
#' @param res CNA segmentation results from pgx.segmentCN()
#' @param annot Data frame with sample annotations
#' @param pca.filter Filter samples by PCA clustering (-1 to disable)
#' @param lwd Line width for chromosome lines
#' @param downsample Downsample CNA matrix for plotting (integer factor)
#' @param order.by Order samples by "clust" or "annot"
#' @param clip Clip copy number ratio values (0 to disable)
#' @param lab.cex Label cex size
#'
#' @return None. Plot is produced as a side-effect.
#'
#' @description
#' Generates a heatmap to visualize copy number alterations across samples.
#'
#' @details
#' This function takes CNA segmentation results from pgx.segmentCN() and
#' plots a heatmap of the copy number ratios. Samples can be ordered by
#' annotations or clustering. The full CNA matrix can be downsampled for
#' easier visualization. Chromosome lines and sample annotations are added.
#'
#' @export
pgx.plotCNAHeatmap <- function(ngs, res, annot = NA, pca.filter = -1, lwd = 1,
                               downsample = 10,
                               order.by = "clust", clip = 0, lab.cex = 0.6) {
  cna <- res$cna
  chr <- res$chr
  chr <- as.character(chr)
  pos <- res$pos

  ## ---------------------------------------------------------------------
  ## Downsample if needed
  ## ---------------------------------------------------------------------
  if (downsample > 1) {
    ## Downsample
    cat("downsampling CNA matrix...\n")
    n <- downsample
    jj <- as.vector(sapply(1:((nrow(cna) + n) / n), rep, n))[1:nrow(cna)]
    cna <- apply(cna, 2, function(x) tapply(x, jj, mean))
    gg <- tapply(rownames(res$cna), jj, paste, collapse = ",")
    rownames(cna) <- gg
    j1 <- which(!duplicated(jj))
    chr <- chr[j1]
    pos <- tapply(pos, jj, mean)
  }

  ## ---------------------------------------------------------------------
  ## take out small groups/chromsomes
  ## ---------------------------------------------------------------------
  ii <- which(chr %in% names(which(table(chr) > 3)))
  cna <- cna[ii, ]
  pos <- pos[ii]
  chr <- chr[ii]

  ## ensure order on chrpos
  ichr <- as.integer(sub("X", 23, sub("Y", 24, sub("chr", "", chr))))
  jj <- order(ichr, pos)
  cna <- cna[jj, ]
  chr <- chr[jj]
  pos <- pos[jj]

  cna <- cna - rowMeans(cna, na.rm = TRUE)
  cna <- cna / max(abs(cna), na.rm = TRUE)
  cna <- tanh(1.3 * cna)
  cna <- t(t(cna) - apply(cna, 2, median))

  if (pca.filter > 0) {
    k <- 20
    k <- pca.filter
    k <- ceiling(min(0.33 * ncol(cna), k))
    sv <- irlba::irlba(cna, nv = k)
    cna2 <- sv$u[, 1:k] %*% diag(sv$d[1:k]) %*% t(sv$v[, 1:k])
    colnames(cna2) <- colnames(cna)
    rownames(cna2) <- rownames(cna)
    cna <- cna2
  }

  hc <- NULL
  sv1 <- NULL
  if (order.by == "pc1") {
    ## by default order on SV1
    sv1 <- irlba::irlba(cna, nv = 1)$v[, 1]
    jj <- order(sv1)
    cna <- cna[, jj]
    sv1 <- sv1[jj]
  } else {
    ## order by hierarchical clustering
    jj <- Matrix::head(order(-apply(cna, 1, sd)), 1000)
    hc <- fastcluster::hclust(dist(t(cna[jj, ])), method = "ward.D2")
    cna <- cna[, hc$order]
  }

  ## create annation matrix
  ann.mat <- NULL
  if (!is.null(annot)) {
    if (is.na(annot)) {
      k <- c(grep("cell.type|tissue|cluster|group",
        colnames(ngs$samples),
        ignore.case = TRUE
      ), 1)[1]
    } else {
      k <- match(annot, colnames(ngs$samples))
    }
    k
    y <- as.character(ngs$samples[colnames(cna), k])
    ny <- length(setdiff(unique(y), NA))
    if (ny >= 2) {
      y[is.na(y)] <- "_"
      ann.mat <- model.matrix(~ 0 + y)
      colnames(ann.mat) <- sub("^y", "", colnames(ann.mat))
      rownames(ann.mat) <- colnames(cna)
    }
  }

  BLUERED2 <- colorRampPalette(c("blue3", "white", "red3"))

  ## ---------- do plotting ------------

  par(mgp = c(0.8, 0.4, 0))
  wa <- 0.1
  if (!is.null(ann.mat)) wa <- 0.05 + 0.016 * ncol(ann.mat)
  plotly::layout(matrix(1:3, 1, 3), widths = c(0.2, 0.7, wa))

  if (!is.null(hc)) {
    par(mar = c(8, 2, 12, 0))
    plot(as.dendrogram(hc),
      horiz = TRUE, leaflab = "none",
      yaxs = "i", xaxt = "n", yaxt = "n"
    )
  } else if (!is.null(sv1)) {
    par(mar = c(8, 3, 12, 0.3))
    barplot(sv1,
      horiz = TRUE, border = NA, col = "grey50", width = 0.1,
      space = 0, yaxs = "i", xaxt = "n"
    )
    mtext("PC1", side = 2, cex = 0.8)
  } else {
    frame()
  }

  ## main heatmap
  par(mar = c(8, 0.2, 12, 0))
  cna0 <- cna
  cna0 <- tanh(3 * cna0)
  cna0[which(abs(cna0) < clip)] <- NA
  Matrix::image(1:nrow(cna), 1:ncol(cna), cna0[, ],
    col = BLUERED2(16),
    ylab = "samples", xlab = "DNA copy number  (log2R)",
    yaxt = "n", yaxs = "i", xaxt = "n", xaxs = "i",
    zlim = c(-1, 1) * 1.0
  )

  ichr <- as.integer(sub("X", 23, sub("Y", 24, sub("chr", "", chr))))
  chrbrk <- which(diff(ichr) != 0)
  chrmid <- c(0, chrbrk) + diff(c(0, chrbrk, nrow(cna))) / 2
  abline(v = chrbrk, col = "grey50", lty = 1, lwd = lwd)
  chrlen <- length(unique(chr))
  j0 <- seq(1, chrlen, 2)
  j1 <- seq(2, chrlen, 2)
  mtext(unique(chr)[j0], side = 3, at = chrmid[j0], cex = lab.cex, line = 0.25)
  mtext(unique(chr)[j1], side = 3, at = chrmid[j1], cex = lab.cex, line = 0.9)

  if (!is.null(ann.mat)) {
    par(mar = c(8, 0.5, 12, 2))
    Matrix::image(1:ncol(ann.mat), 1:nrow(ann.mat), t(ann.mat),
      col = rev(grey.colors(2)), xlab = "", ylab = "",
      yaxt = "n", yaxs = "i", xaxt = "n", xaxs = "i"
    )
    mtext(colnames(ann.mat),
      side = 3, at = 1:ncol(ann.mat),
      las = 3, cex = lab.cex, line = 0.25
    )
  } else {
    frame()
  }

  ## done plotting
}
