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
pgx.inferCNV <- function(pgx, refgroup = NULL, progress = NULL) {
  ## InferCNV: Inferring copy number alterations from tumor single
  ## cell RNA-Seq data
  ##
  ## https://github.com/broadinstitute/inferCNV/wiki
  ##

  # Prepare input data
  annots <- pgx$samples[, "group", drop = FALSE]

  #Get gene symbol and chromosome location info
  genes_info <- data.table::data.table(pgx$genes[rownames(data), , drop = FALSE])

  # Get approx start and stop
  genes_info[, c("chr", "start", "stop") :=
               .(paste0("chr", chr),
                 pos,
                 pos  + tx_len)]

  ## filter known genes
  jj <- which(genes_info[["chr"]] %in% paste0("chr", c(1:22, "X", "Y")) &
                !is.na(genes_info[["start"]]) &
                !is.na(genes_info[["stop"]]))
  genes_info <- genes_info[jj, ]

  ## Reshape objects for infercnv
  gg <- intersect(genes_info[[1]], rownames(ngs$counts))
  data <- data[gg, ]
  genes_info <- genes_info[, .SD, .SDcols = c("chr", "start", "stop")]
  genes_info <- as.data.frame(genes_info)
  rownames(genes_info) <- rownames(data)

  ## take out tiny groups
  selgrp <- names(which(table(annots[, 1]) >= 2))
  kk <- which(annots[, 1] %in% selgrp)
  data <- data[, kk]
  annots <- annots[colnames(data), , drop = FALSE]

  ## Run inferCNV
  infercnv_obj <- infercnv::CreateInfercnvObject(
    raw_counts_matrix = data,
    gene_order_file = genes_info,
    annotations_file = annots,
    ref_group_names = refgroup
  )

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

  genes_info <- genes_info[rownames(cnv), ]
  pos <- (genes_info$start + genes_info$stop) / 2
  ichr <- as.integer(sub("X", 23, sub("Y", 24, sub("chr", "", genes_info$chr))))
  jj <- order(ichr, pos)
  pos <- pos[jj]
  chr <- as.character(genes_info$chr)[jj]
  logcnv <- log2(cnv[jj, ] / mean(cnv, na.rm = TRUE)) ## logarithmic



  img <- png::readPNG(img.file)
  res <- list(cna = logcnv, chr = chr, pos = pos, png = img)

  ## clean up folder??
  unlink(out_dir, recursive = TRUE)

  return(res)
}
