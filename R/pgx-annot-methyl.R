##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

## ================================================================================
## =========================== METHYLATION FUNCTIONS ==============================
## ================================================================================

#' Handle entire annotation for 450K+850K methylomics.
#' @param organism Strictly human (for now).
#' @param probes 450K or 850K array methylation probes. 850K is the EPIC array.
#' @return Ann. dframe: feature;symbol;gene_name;gene_title;chr;source;position;uniprot;
#' @export
annotate_methylomics <- function(organism = "Human", probes = probes, meth_type = "450K array") {

  msg <- function(...) message("[playbase::annotate_methylomics] ", ...)

  if (!organism %in% c("Human", "Homo sapiens")) {
    msg("Error: annotation of methylomics probes limited to only Human.")
    return(NULL)
  }

  msg("Annotating methylomics data...")

  if (!meth_type %in% c("450K array", "EPIC array")) {
    msg("meth_type not recognized. Must be '450K array' or 'EPIC array'")
    msg("Defaulting to '450K array'")
  }
  pkg <- "IlluminaHumanMethylation450kanno.ilmn12.hg19"
  if (meth_type == "EPIC array") pkg <- "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"
  require(pkg, character.only = TRUE)

  annot <- as.data.frame(minfi::getAnnotation(get(pkg)))
  kk <- c("chr", "pos", "strand", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "Relation_to_Island")
  annot <- annot[probes, intersect(kk, colnames(annot)), drop = FALSE]    
  colnames(annot)[colnames(annot) == "UCSC_RefGene_Name"] <- "symbol"
  colnames(annot)[colnames(annot) == "UCSC_RefGene_Group"] <- "genomic_location"    

  kk <- c("symbol", "genomic_location")
  for(i in 1:ncol(annot[, kk, drop = FALSE])) {
    ff <- lapply(annot[, kk[i]], function(x) strsplit(x, ";")[[1]])
    annot[, kk[i]] <- unlist(lapply(ff, function(x) paste0(unique(x),collapse=";")))
  }
  
  ff <- ifelse(is.na(annot$symbol) | annot$symbol == "", "probe", annot$symbol)
  ff <- playbase::make_unique(ff)
  genes <- getGeneAnnotation("Human", ff)
  genes <- genes[match(ff, genes$feature), ]
  rownames(genes) <- genes$feature <- rownames(annot)
  genes$symbol <- genes$gene_name <- annot$symbol
  kk <- setdiff(colnames(annot), colnames(genes))
  if (length(kk)) genes <- cbind(genes, annot[, kk, drop = FALSE])
  rm(annot, ff)

  msg("Annotation completed\n")
  return(genes)
  
}

## ADD meth_type to mergeCpG!!

#' @title mergeCpG
#' @description Collapse CpG into genes. Use average beta value.
#' @param X Matrix of methyation values. Can be M or Beta. Beta values will be used.
#' @param genes Annotation matrix: pgx$genes. if NULL, automatically retrieved.
#' @param collapse.by collapse.by variable. Can be later flexibly user-defined.
#' @return List of 2 elements: [[1]] collapsed CpG matrix; [[2]] updated annotation dataframe.
#' @export
mergeCpG <- function(data, genes = NULL, collapse.by = "gene") {

  msg <- function(...) message("[playbase::mergeCpG] ", ...)
  msg("Methylomics: collapsing CpG into genes.")

  if (!collapse.by %in% c("gene")) {
    msg("'collapse.by' must be one of gene,... Returning input matrix.\n")
    return(X)
  }

  X <- mToBeta(data)

  if (is.null(genes)) {
    genes <- annotate_methylomics(probes = rownames(X), meth_type = "450K array")
  }
  
  jj <- which(!as.character(genes$symbol) %in% c("", "NA", NA))
  genes <- genes[jj, , drop = FALSE]
  if ("genomic_location" %in% colnames(genes)) {
    hh <- grep("Body|3'UTR", genes$genomic_location, ignore.case = TRUE)
    if (length(hh) > 0) genes <- genes[-hh, , drop = FALSE]
  }

  kk <- intersect(rownames(X), rownames(genes))
  if (length(kk) == 0) {
    msg("No shared features between X and annotation. Returning input matrix.\n")
    return(X)
  }
  X <- X[kk, , drop = FALSE]
  genes <- genes[kk, , drop = FALSE]

  ff <- unique(as.character(genes$symbol))

  ## Collapse CpGs
  LL=list()
  for(i in 1:length(ff)) {
    jj <- which(genes$symbol == ff[i])
    LL[[ff[i]]] <- t(as.matrix(colMeans(X[jj, , drop = FALSE], na.rm = TRUE)))
  }
  X <- do.call(rbind, LL)
  rownames(X) <- ff

  ## Collapse annotation
  LL=list()
  for(i in 1:length(ff)) {
    jj <- which(genes$symbol == ff[i])
    genes1 <- genes[jj[1], , drop = FALSE]
    rownames(genes1) <- genes1$feature <- ff[i]
    genes1$cpg_probe <- paste0(unique(genes$feature[jj]), collapse=";")
    uniprots <- unlist(lapply(genes$uniprot[jj], function(x) strsplit(x, split=";")[[1]]))
    uniprots <- uniprots[!is.na(uniprots)]
    genes1$uniprot <- paste0(unique(uniprots), collapse=";")
    kk <- c("pos", "strand", "genomic_location", "Relation_to_Island")
    kk <- intersect(kk, colnames(genes))
    for(k in 1:length(kk)) {
      genes1[, kk[k]] <- paste0(unique(genes[jj, kk[k]]), collapse=";")
    }
    LL[[ff[i]]] <- genes1
  }
  genes <- do.call(rbind, LL)
  rm(LL)

  ## Safety alignment
  kk <- intersect(rownames(X), rownames(genes))
  X <- X[kk, , drop = FALSE]
  genes <- genes[kk, , drop = FALSE]
  
  msg("Mapping completed. Final matrix:", nrow(X), " regions.\n")
  return(list(data = X, genes = genes))

}

#' @title infer_sex_methyl
#' @description Infer biological sex from array methylomics data
#' @param X Matrix of methyation values. Can be M or Beta. Beta values will be used.
#' @param genes Annotation matrix: pgx$genes. if NULL, automatically retrieved.
#' @param meth_type Array type. '450K array' or 'EPIC array'
#' @return list of 3 elements: named vector of sex, chrX and chrY median beta values.
#' @export
infer_sex_methyl <- function(data, genes = NULL, meth_type = "450K array") {

  msg <- function(...) message("[playbase::infer_sex_methyl] ", ...)
  msg("Methylomics: infer biological sex using sex-linked CpG methylation profiles.")

  X <- mToBeta(data)

  if (is.null(genes)) {
    c1 <- is.null(meth_type)
    c2 <- !meth_type %in% c("450K array", "EPIC array")
    if (c1 | c2) meth_type = "450K array"
    pkg <- "IlluminaHumanMethylation450kanno.ilmn12.hg19"
    if (meth_type == "EPIC array") pkg <- "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"
    require(pkg, character.only = TRUE)
    genes <- as.data.frame(minfi::getAnnotation(get(pkg)))
  }
  genes <- genes[grep("X|chrX|Y|chrY", rownames(genes)), , drop = FALSE]
  kk <- intersect(rownames(X), rownames(genes))
  if (length(kk) == 0) {
    msg("No X- or Y-linked CpG probes detected in input matrix. Exiting")
    return(NULL)
  }
  X <- X[kk, , drop = FALSE]
  genes <- genes[kk, , drop = FALSE]

  chrX <- rownames(genes)[grep("X|chrX", rownames(genes))]
  chrY <- rownames(genes)[grep("Y|chrY", rownames(genes))]

  x_med <- y_med <- pred_sex <- NULL
  if (length(chrX) > 0) x_med <- colMedians(X[chrX, , drop = FALSE ], na.rm = TRUE)
  if (length(chrY) > 0) y_med <- colMedians(X[chrY, , drop = FALSE ], na.rm = TRUE)

  if (!is.null(y_med)) {
    pred_sex <- ifelse(y_med > 0.1, "M", "F")
  } else if (!is.null(x_med)) {
    pred_sex <- ifelse(x_med > 0.4, "M", "F")
  }
  names(pred_sex) <- colnames(X)
  
  return(list(pred_sex = pred_sex, x_med = x_med, y_med = y_med))

}

