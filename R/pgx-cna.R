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
