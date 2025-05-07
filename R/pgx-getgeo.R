##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


## -------------------------------------------------------------------------------------
## Query GEO
## -------------------------------------------------------------------------------------

#' @title pgx.getGEOseries
#' @description Download and process GEO dataset
#'
#' @param id GEO series ID to download
#' @param archs.h5 Path to ARCHS4 HDF5 file containing GEO data
#' @param convert.hugo Logical, convert symbols to HUGO if TRUE
#'
#' @return List containing processed counts, samples, and genes
#'
#' @details This function downloads and processes a GEO dataset specified by the ID.
#' It first checks if the data is in the ARCHS4 file. If not, it retrieves it from GEO using geoquery.
#' The data matrices are subset to intersecting samples, converted to HUGO symbols if specified,
#' and duplicate symbols are summed.
#'
#' The result is a list containing the expression matrix, sample metadata, and gene symbols.
#'
#' @export
pgx.getGEOseries <- function(id, archs.h5 = "human_matrix.h5", convert.hugo = TRUE) {
  ## Highly automagic download of GEO datasets from different
  ## sources with automatic probe/gene conversion and creating
  ## autocontrasts. The GEO series is first searched in a locally
  ## stored ARCHS4 H5 file, then if it is available at the recount
  ## database, if not it is retrieved from GEO using geoquery.

  is.valid.id <- is.GEO.id.valid(id) 
  if (!is.valid.id) stop("[pgx.getGEOseries] FATAL: ID is invalid. Exiting.")
  id <- as.character(id)

  ## get data/pheno matrices
  geo <- pgx.getGEOcounts(id, archs.h5 = archs.h5)
  counts <- geo$expr

  ## get sample info
  meta <- pgx.getGeoMetadata(id)

  ## conform matrices
  samples <- intersect(rownames(meta), colnames(counts))
  meta <- meta[samples, , drop = FALSE]
  counts <- counts[, samples, drop = FALSE]

  ## convert to latest official HUGO???
  if (convert.hugo) {
    symbol <- alias2hugo(rownames(counts)) ## auto-detect mouse/human
    rownames(counts) <- symbol
  }

  ## sum up values duplicated symbols
  ndup <- sum(duplicated(rownames(counts)))
  if (ndup > 0) {
    counts1 <- tapply(1:nrow(counts), symbol, function(i) Matrix::colSums(counts[i, , drop = FALSE]))
    counts <- do.call(rbind, counts1)
    remove(counts1)
  }

  ## get annotation
  genes <- ngs.getGeneAnnotation(rownames(counts))

  ## get categorical phenotypes
  meta1 <- apply(meta, 2, trimsame)
  rownames(meta1) <- rownames(meta)
  sampleinfo <- pgx.discretizePhenotypeMatrix(meta1, min.ncat = 2,
    max.ncat = 20, remove.dup = TRUE)
  sampleinfo <- data.frame(sampleinfo, stringsAsFactors = FALSE, check.names = FALSE)

  ## automagically create contrast matrix
  contrasts <- NULL
  if (NCOL(sampleinfo) > 0) {
    mingrp <- 3
    slen <- 15
    ref <- NA
    ct <- pgx.makeAutoContrasts(sampleinfo, mingrp = 3, slen = 20, ref = NA)
    if (is.null(ct)) {
      ct <- pgx.makeAutoContrasts(sampleinfo, mingrp = 2, slen = 20, ref = NA)
    }
    if (!is.null(ct$exp.matrix)) {
      contrasts <- ct$exp.matrix
    } else {
      contrasts <- ct$design %*% ct$contr.matrix
    }
  }

  info <- pgx.getGeoExperimentInfo(id)

  out <- list(
    counts = counts,
    genes = genes,
    samples = sampleinfo,
    contrasts = contrasts,
    meta = meta,
    info = info,
    source = geo$source
  )

  return(out)

}


#' @describeIn pgx.getGEOcounts Download count data from GEO. First check
#' if the GEO ID is in archs5, then in recount. If not, try to get from GEO. 
#' @return List of counts matrix and source.
#' @param id GEO accession ID.
#' @param archs.h5 Path to archs.h5 dataset.
#' @export
pgx.getGEOcounts <- function(id, archs.h5) {

  is.valid.id <- is.GEO.id.valid(id) 
  if (!is.valid.id) stop("[pgx.getGEOcounts] FATAL: ID is invalid. Exiting.")
  id <- as.character(id)

  expr=NULL; src=""
  
  if (!is.null(archs.h5) && is.null(expr)) {
    message("[pgx.getGEOcounts]: pgx.getGEOcounts.archs4....")
    expr <- pgx.getGEOcounts.archs4(id, archs.h5)
    if (!is.null(expr)) src <- "ARCHS4"
  }

  if (is.null(expr)) {
    message("[pgx.getGEOcounts]: pgx.getGEOcounts.recount....")
    expr <- pgx.getGEOcounts.recount(id)
    if (!is.null(expr)) src <- "recount"
  }

  if (is.null(expr)) {
    message("[pgx.getGEOcounts]: pgx.getGEOcounts.GEOquery....")
    expr <- pgx.getGEOcounts.GEOquery(id)
    if (!is.null(expr)) src <- "GEO"
  }

  if (is.null(expr)) {
    cat("WARNING:: Could not get GEO expression. please download manually.\n")
    return(NULL)
  }
  
  list(expr = expr, source = src)

}

#' @describeIn pgx.getGEOmetadata Download and extract the metadata from a GEO ID.
#' It attemtps without GSEMatrix first, and then with GSEMatrix.
#' @export
pgx.getGEOmetadata <- function(id) {

  is.valid.id <- is.GEO.id.valid(id) 
  if (!is.valid.id) stop("[pgx.getGEOmetadata] FATAL: ID is invalid. Exiting.")
  id <- as.character(id)

  meta <- pgx.getGEOmetadata.fromGSM(id) ## no GSEMatrix 
  if (is.null(meta)) meta <- pgx.getGEOmetadata.fromEset(id) ## with GSEMatrix
  
  ## Sometimes the phenotype is coded in the title string
  #if ("title" %in% colnames(meta) && NCOL(meta) == 0) {
  #  px <- title2pheno(meta$title, split = NULL, trim = TRUE, summarize = TRUE)
  #  if (!is.null(px) && NCOL(px) > 0 && is.null(meta))
  #    pheno <- px
  #  if (!is.null(px) && NCOL(px) > 0 && !is.null(pheno))
  #    pheno <- cbind(pheno, px)
  #}

  return(meta)

}


## -------------------------------------------------------------------------------------
## Query GEO expression
## -------------------------------------------------------------------------------------

#' @describeIn pgx.getGEOcounts.archs4 Downloads and extracts gene expression count
#' data from a GEO series stored in an HDF5 file (if available). It searches the
#' HDF5 file metadata to find samples matching the input GEO series ID, and returns
#' the count matrix for those samples.
#' @export
pgx.getGEOcounts.archs4 <- function(id, h5.file) {

  is.valid.id <- is.GEO.id.valid(id) 
  if (!is.valid.id) stop("[pgx.getGEOcounts.archs4] FATAL: ID is invalid. Exiting.\n")
  id <- as.character(id)

  if (is.null(h5.file) || h5.file == "")
    stop("[pgx.getGEOcounts.archs4] FATAL: invalid path to h5.file ID. Exiting.\n")

  sample.series <- rhdf5::h5read(h5.file, "meta/Sample_series_id")
  sample.series <- strsplit(as.character(sample.series), split = "Xx-xX")
  idx <- which(sapply(sample.series, function(s) id %in% s))
  if (!id %in% sample.series) {
    message("[pgx.getGEOcounts.archs4] WARNING: series ", id, " not in ARCHS4. Exiting.\n")
    return(NULL)
  }
  message("[pgx.getGEOcounts.archs4] Series ", id, " found in ARCHS4.")

  ## get matrix
  counts <- rhdf5::h5read(h5.file, "data/expression", index = list(NULL, idx))
  sample.acc <- rhdf5::h5read(h5.file, "meta/Sample_geo_accession")
  gene_name <- rhdf5::h5read(h5.file, "meta/genes")
  colnames(counts) <- sample.acc[idx]
  rownames(counts) <- gene_name

  ## ensure counts
  qx <- as.numeric(stats::quantile(counts, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = T))
  is.count <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2) ## from GEO2R script
  if (!is.count) counts <- 2 ** counts
  
  ## rm missing genes and sum linear intensities
  jj <- which(!is.na(gene_name) & gene_name != "")
  counts <- counts[jj, ]
  gene_name <- gene_name[jj]
  counts <- tapply(1:nrow(counts), gene_name, function(ii) {
    Matrix::colSums(counts[ii, , drop = FALSE], na.rm = TRUE)
  })
  counts <- do.call(rbind, counts)
  message("[pgx.getGEOcounts.archs4] Success!")

  return(counts)

}


#' @describeIn pgx.getGEOcounts.recount Downloads and processes gene-level count data
#' for a GEO series from the recount database. It takes a GEO ID, searches recount,
#' downloads the RangedSummarizedExperiment object, and returns the count matrix.
#' Vignette recount-quickstart.html
#' @export
pgx.getGEOcounts.recount <- function(id) {

  is.valid.id <- is.GEO.id.valid(id) 
  if (!is.valid.id) stop("[pgx.getGEOcounts.recount] FATAL: ID is invalid. Exiting.")
  id <- as.character(id)

  project_info <- recount::abstract_search(id)
  pid <- project_info$project
  if (length(pid) == 0) {
    message("[pgx.getGEOcounts.recount] WARNING: series ", id, " not in recount. Exiting.\n")
    return(NULL)
  }
  message("[pgx.getGEOcounts.recount] Series ", id, " found in recount database.")
  
  ## Download the gene-level RangedSummarizedExperiment data
  outdir <- file.path(tempdir(), pid)
  cc <- try(recount::download_study(pid, outdir = outdir), silent = TRUE)
  if (inherits(cc, "try-error")) {
    message("[pgx.getGEOcounts.recount] Error: could not retrieve ", id, ". Exiting.\n")
    return(NULL)
  }

  ## Load the data
  load(file.path(outdir, "rse_gene.Rdata"))

  ## Scale counts by taking into account the total coverage per sample
  rse <- recount::scale_counts(rse_gene)
  counts <- MultiAssayExperiment::assay(rse)

  ## ensure counts
  qx <- as.numeric(stats::quantile(counts, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = T))
  is.count <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2) ## from GEO2R script
  if (!is.count) counts <- 2 ** counts
  message("[pgx.getGEOcounts.recount] Success!")
  
  return(counts)

}


#' @describeIn  pgx.getGEOcounts.GEOquery retrieves expression count data for a
#' GEO accession ID using the GEOquery package. It downloads the series matrix data,
#' platform metadata, and probe annotations from GEO into R objects.
#' @export
pgx.getGEOcounts.GEOquery <- function(id) {
  
  is.valid.id <- is.GEO.id.valid(id) 
  if (!is.valid.id) stop("[pgx.getGEOcounts.GEOquery] FATAL: ID is invalid. Exiting.")
  id <- as.character(id)

  gse <- try(GEOquery::getGEO(id, GSEMatrix = TRUE, getGPL = TRUE), silent = TRUE)
  if (inherits(gse, "try-error")) {
    message("[pgx.getGEOcounts.GEOquery] Error: getGEO Efailed to retrieve ", id, "\n")
    return(NULL)
  }

  require(Biobase)
  has.expr <- sapply(gse, function(x) nrow(Biobase::exprs(x)) > 0)
  
  if (any(has.expr)) {
    gse <- gse[which(has.expr)]    
  } else {
    message("[pgx.getGEOcounts.GEOquery] WARNING: no data found in ", id, " from GEO.\n")
    supp_file <- sapply(gse, function(g) g@experimentData@other$supplementary_file)
    if (class(supp_file) == "character")
      message("Supplementary file available: ", paste(supp_file, collapse = " "), "\n")
    return(NULL)
  }

  ## select preferred platform is multiple exists
  k=1; expr.list=list()
  for (k in 1:length(gse)) {

    eset <- gse[[k]]
    ex <- exprs(eset)
    if (ncol(ex) <= 3) {
      message("[pgx.getGEOcounts.GEOquery] WARNING: ", id, " contains <= 3 samples. Skipping.\n")
      next()
    }

    ## perform linear transformation (unlog) if required
    qx <- as.numeric(stats::quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = T))
    is.count <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0) ||
      (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2) ## from GEO2R script
    if (!is.count) ex <- 2 ** ex

    ## featuredata
    has.fdata <- !is.null(Biobase::fData(eset)) && NCOL(Biobase::fData(eset)) > 0
    if (has.fdata) {
      fdata <- Biobase::fData(eset)
    } else {
      gpl.annot <- GEOquery::getGEO(eset@annotation)
      fdata <- GEOquery::Table(gpl.annot)
    }
    if ("ID" %in% colnames(fdata)) {
      cm <- intersect(rownames(ex), as.character(fdata$ID))
      jj <- match(cm, as.character(fdata$ID))
      fdata <- fdata[jj, , drop = FALSE]
      rownames(fdata) <- cm
      fdata <- fdata[, colnames(fdata) != "ID"]
      ex <- ex[cm, , drop = FALSE]
    } else {
      cm <- intersect(rownames(ex), rownames(fdata))
      fdata <- fdata[cm, , drop = FALSE]
      ex <- ex[cm, , drop = FALSE]
    }
    
    ## get symbol from featuredata; clean and sum linear intensities
    fsymbol <- pgx.getSymbolFromFeatureData(fdata)    
    jj <- which(!is.na(fsymbol) & fsymbol != "")
    ex <- ex[jj, ]
    fsymbol <- fsymbol[jj]
    fsymbol <- gsub(" /// ", ";", fsymbol)
    ex2 <- tapply(1:nrow(ex), fsymbol, function(ii) {
      Matrix::colSums(ex[ii, , drop = FALSE], na.rm = TRUE) ## not log!!
    })
    ex2 <- do.call(rbind, ex2)
    expr.list[[names(gse)[k]]] <- ex2

  }

  if (length(expr.list) == 0) return(NULL)
  
  if (length(expr.list) > 1) {
    ## merge/join all expressions
    probes <- sort(unique(unlist(lapply(expr.list, rownames))))
    samples <- sort(unique(unlist(lapply(expr.list, colnames))))
    expr.list2 <- lapply(expr.list, function(x) {
      x[match(probes, rownames(x)), match(samples, colnames(x))]
    })
    expr.list2 <- lapply(expr.list2, function(x) { x[is.na(x)] <- 0; x })
    expr <- Reduce("+", expr.list2)
    colnames(expr) <- samples
    rownames(expr) <- probes
  } else {
    expr <- expr.list[[1]]
  }

  return(expr) ## linear intensities

}


## -------------------------------------------------------------------------------------
## Query GEO metadata
## -------------------------------------------------------------------------------------

#' @describeIn pgx.getGEOexperimentInfo Retrieves GEO accession ID info using GEOquery.
#' @param id GEO accession ID
#' @export
pgx.getGEOexperimentInfo <- function(id) {

  is.valid.id <- is.GEO.id.valid(id) 
  if (!is.valid.id) stop("[pgx.getGEOexperimentInfo] FATAL: ID is invalid. Exiting.")
  id <- as.character(id)

  suppressMessages(gse <- try(GEOquery::getGEO(id, GSEMatrix = FALSE, getGPL = FALSE), silent = TRUE))
  if (inherits(gse, "try-error")) {
    message("[pgx.getGEOexperimentInfo] Error: GEOquery::getGEO failed to get", id, ".\n")
    return(NULL)
  }

  return(gse@header)

}


#' @describeIn pgx.getGeoMetadata.fromGSM Retrieves metadata for a GEO accession ID.
#' It attemps without GSE Series Matrix files from GEO.
#' @param id GEO accession ID.
#' @export
pgx.getGEOmetadata.fromGSM <- function(id) {

  is.valid.id <- is.GEO.id.valid(id)
  if (!is.valid.id) stop("[pgx.getGEOmetadata.fromGSM] FATAL: ID is invalid. Exiting.")
  id <- as.character(id)

  message("[pgx.getGEOmetadata.fromGSM] Attempt to download metadata without GSEMatrix...")
  suppressMessages(
    gse <- try(GEOquery::getGEO(id, GSEMatrix = FALSE, getGPL = FALSE), silent = TRUE)
  )
  if (inherits(gse, "try-error")) {
    message("[pgx.getGEOmetadata.fromGSM] Error: getGEO failed to retrieve metadata for ", id, "\n")
    return(NULL)
  }

  if (length(gse@gsms) == 0) {
    message("[pgx.getGEOmetadata.fromGSM] WARNING: no GSM information in object. Exiting. \n")
    return(NULL)
  }

  ## get metadata
  #summary <- gse@header$summary
  gsm.title <- sapply(gse@gsms, function(g) g@header$title)
  gsm.source <- sapply(gse@gsms, function(g) g@header$source_name_ch1)
  gsm.gpl <- sapply(gse@gsms, function(g) g@header$platform_id)
  gsm.samples <- gse@header$sample_id
  meta <- data.frame(
    GPL = gsm.gpl,
    GSM = gsm.samples,
    title = gsm.title,
    source = gsm.source,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  ch1_info <- lapply(gse@gsms, function(g) g@header$characteristics_ch1)
  cm <- intersect(names(ch1_info), gsm.samples)
  if (!is.null(ch1_info) && length(cm)) {
    gsm.samples <- gsm.samples[match(cm, gsm.samples)]
    ch1_info <- ch1_info[match(cm, names(ch1_info))]
    meta <- meta[match(cm, meta$GSM), ]
    ch1_info <- lapply(ch1_info, function(x) sub("^Clinical info: ", "", x))
    ch1_vars <- unique(unlist(lapply(ch1_info, function(x) trimws(sub("[:=].*", "", x)))))
    ch1_info <- lapply(ch1_info, function(x) {
      xvar <- trimws(sub("[:=].*", "", x))
      x <- trimws(sub(".*[:=] ", "", x))
      names(x)=xvar; x=x[match(ch1_vars, names(x))]
      return(x)
    })
    ch1_info <- do.call(rbind, ch1_info)
    colnames(ch1_info) <- ch1_vars
    meta <- data.frame(cbind(meta, ch1_info), stringsAsFactors = FALSE, check.names = FALSE)
    colnames(meta) <- gsub("[ ]", "_", colnames(meta))

    message("[pgx.getGEOmetadata.fromGSM] Success!") 
    return(meta)
  } else {
    message("[pgx.getGEOmetadata.fromGSM] WARNING: no shared samples between GSM & ch1_info. Exiting.\n")
    return(NULL)
  }

  #is.underscored <- length(gsm.title) && all(grepl("_", gsm.title))
  #title_info <- NULL
  ## NEED RETHINK!!!!!!!!!!!!!!!!!!
  #if (FALSE && is.underscored) {
  #  title2 <- trimws(gsm.title)
  #  title_info <- eset.parsePhenoFromTitle(title2, split = "_")
  #}
  #if (!is.null(title_info)) sample_info <- cbind(sample_info, title_info)

}

#' @describeIn pgx.getGEOmetadata.fromEset Retrieves sample metadata from Eset.
#' It downloads the GSE Series Matrix files from GEO and extract sample metadata.
#' @param id GEO accession ID 
#' @export
pgx.getGEOmetadata.fromEset <- function(id) {

  is.valid.id <- is.GEO.id.valid(id)
  if (!is.valid.id) stop("[pgx.getGEOmetadata.fromEset] FATAL: ID is invalid. Exiting.")
  id <- as.character(id)

  message("[pgx.getGEOmetadata.fromEset] Attempt to download metadata with GSEMatrix...")
  suppressMessages(
    gse <- try(GEOquery::getGEO(id, GSEMatrix = TRUE, getGPL = FALSE), silent = TRUE)
  )
  if (inherits(gse, "try-error")) {
    message("[pgx.getGEOmetadata.fromEset] Error: getGEO failed to retrieve metadata for ", id, "\n")
    return(NULL)
  }

  nsamples <- sapply(gse, function(s) nrow(Biobase::pData(Biobase::phenoData(s))))
  gse <- gse[which(nsamples >= 3)]
  meta.list <- lapply(gse, function(x) pgx.getGEOmetadata.fromEset.helper(x))
  meta <- do.call(rbind, meta.list)
  meta <- data.frame(meta, stringsAsFactors = FALSE, check.names = FALSE)
  colnames(meta) <- gsub("[ ]", "_", colnames(meta))
  message("[pgx.getGEOmetadata.fromEset] Success!")
  rm(meta.list)

  return(meta)

}


#' @describeIn pgx.getGEOmetadata.fromEset.helper Extracts phenotype data from a Eset object.
#' @param eset Eset object
#' @export
pgx.getGEOmetadata.fromEset.helper <- function(eset) {

  if (!class(eset) %in% "ExpressionSet") {
    message("[pgx.getGEOmetadata.fromEset.helper] Error: eset must be of class 'ExpressionSet'")
    return(NULL)
  }

  ## get metadata
  meta0 <- Biobase::pData(Biobase::phenoData(eset))
  gsm.gpl <- as.character(meta0$platform_id)
  gsm.samples <- as.character(meta0$geo_accession)
  gsm.title <- as.character(meta0$title)
  gsm.source <- as.character(meta0$source_name_ch1)
  meta <- data.frame(
    GPL = gsm.gpl,
    GSM = gsm.samples,
    title = gsm.title,
    source = gsm.source,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  ch1_info <- NULL
  ch1_sel <- grepl("characteristics_ch1", colnames(meta0))
  if (any(ch1_sel)) {
    ch1_info <- meta0[, ch1_sel, drop = FALSE]
    colnames(ch1_info) <- paste0("characteristics_", 1:ncol(ch1_info))
    ch1_info <- apply(ch1_info, 2, function(x) sub("^Clinical info: ", "", x))
    ch_vars <- apply(ch1_info, 2, function(x) return(unique(sub("[:=].*", "", x))))
    colnames(ch1_info) <- unname(ch_vars)
    ch1_info <- apply(ch1_info, 2, function(x) trimws(sub(".*[:=] ", "", x)))
    meta <- cbind(meta, ch1_info)
    meta <- data.frame(meta, stringsAsFactors = FALSE, check.names = FALSE)    
  }
  rm(meta0)

  ## Base sample_info from characteristics (ch1) column
  #ch1_info <- eset.getCH1(eset)
  # We can get extra information from title
  #is.underscored <- length(gsm.title) && all(grepl("_", gsm.title))
  #title_info <- NULL
  #if (FALSE && is.underscored) {
  #  title2 <- trimws(gsm.title)
  #  title_info <- eset.parsePhenoFromTitle(title2, split = "_")
  #}
  ## All sample_info: from characterisctis_ch1 and title
  #if (!is.null(ch1_info)) sample_info <- cbind(sample_info, ch1_info)
  #if (!is.null(title_info)) sample_info <- cbind(sample_info, title_info)

  return(meta)

}


## -------------------------------------------------------------------------------------
## HELPER functions
## -------------------------------------------------------------------------------------

#' @param id GEO accession ID
#' @return Boolean
#' @details Checks whether GEO accession ID is alphanumeric as per convention. 
#' @export
is.GEO.id.valid <- function(id) {
  if (is.null(id) || id == "" || !grepl("[A-Za-z]",id) || !grepl("[0-9]",id)) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}


#' Extract phenotype data field from ExpressionSet
#' @param eset ExpressionSet object
#' @param field Character string specifying phenotype data field name
#' @return Vector of values for the specified field
#' @details Extracts a specific phenotype data field from an ExpressionSet object.
#' The \code{field} parameter specifies the phenotype data column name to extract.
#' The phenotype data pData slot is extracted from the ExpressionSet; the specified
#' field is returned as a vector.
#' @export
eset.getPhenoData <- function(eset, field) { pData(phenoData(eset))[, field] }


#' Get pheno from title [??]
#' @param title character vector with the titles
#' @param split delimiter character to split on `c(",", ";", "\\|", "_", " ")`.
#' @param trim trim leading and trailing whitespace from each term
#' @param summarize summarize the terms by counting the number of occurrences of each term
#' @export
title2pheno <- function(title, split = NULL, trim = TRUE, summarize = TRUE) {
  if (is.null(split)) {
    split.chars <- c(",", ";", "\\|", "_", " ")
    i=1; ss=c()
    for (i in 1:length(title)) {
      ns <- sapply(split.chars, function(s) sum(gregexpr(s, title[i])[[1]] > 0))
      split0 <- names(ns)[which.max(ns)]
      ns1 <- ns[setdiff(names(ns), " ")]
      if (split0 == " " && any(ns1 > 0)) {
        split0 <- names(ns1)[which.max(ns1)]
      }
      ss[i] <- split0
    }
    split <- names(which.max(table(ss)))
  }

  ## Check if all titles have equal splitted parts (e.g. nicely formatted)
  nsplit <- sapply(title, function(tt) sum(gregexpr(split, tt)[[1]] > 0))
  nsplit.equal <- all(nsplit == nsplit[1])
  nsplit.equal
  if (!nsplit.equal) {
    cat("splitted title terms not equal lengths\n")
    return(NULL)
  }

  ## split
  tt <- as.character(sapply(as.character(title), function(s) trimws(s)))
  ff <- sapply(as.character(tt), strsplit, split = split)

  ## cleanup
  ff <- lapply(ff, trimws) ## trim whitespace
  ff <- lapply(ff, function(s) gsub("hours$|hour$|hrs$|hr$", "h", s)) ## hours
  ff <- lapply(ff, function(s) gsub("[ ][ ]*", " ", s)) ## double space

  ## make dataframe
  F1 <- do.call(rbind, ff)
  F1[is.na(F1)] <- NA
  rownames(F1) <- NULL

  ## Guess names
  getmax.term <- function(s) {
    tt <- table(unlist(strsplit(s, split = "[ _]")))
    vip.tt <- grep("hour|repl|hr|time|treat|infec|pati|sampl", names(tt))
    if (length(vip.tt)) tt[vip.tt] <- 1.1 * tt[vip.tt] ## boost known keywords
    names(which.max(tt))
  }
  maxterm <- apply(F1, 2, function(s) getmax.term(s))
  maxterm <- paste0("_", maxterm)
  colnames(F1) <- maxterm

  ## trims same words/characters on both ends
  if (trim) F1 <- apply(F1, 2, trimsame, summarize = summarize)
  F1
}


#' @describeIn eset.getPhenoData Phenotype data from the title of an
#' ExpressionSet object by splitting the title on a specified delimiter and guessing column names.
#' @export
eset.parsePhenoFromTitle <- function(title, split = NULL) {

  if (!all(grepl(split, title))) return(NULL)
  
  tt <- as.character(sapply(as.character(title), function(s) trimws(s)))
  tt <- sapply(tt, function(s) gsub("[ ]*hours|[ ]*hour|[ ]*hrs|[ ]*hr", "h", s)) ## hours
  tt <- sapply(tt, function(s) gsub("([0-9]*)[ _]h", "\\1h", s))
  tt <- trimsame(tt, split = split, ends = TRUE)
  tt <- gsub(paste0(split, split, split), split, tt)
  tt <- gsub(paste0(split, split), split, tt)
  ff <- sapply(as.character(tt), strsplit, split = split)
  nf <- max(sapply(ff, length))
  ff <- lapply(ff, function(x) Matrix::head(c(x, rep(NA, nf)), nf))

  ## cleanup
  ff <- lapply(ff, trimws) ## trim whitespace
  ff <- lapply(ff, function(s) gsub("hours$|hour$|hrs$|hr$", "h", s)) ## hours
  ff <- lapply(ff, function(s) gsub("[ ][ ]*", " ", s)) ## double space

  F1 <- do.call(rbind, ff)
  F1[is.na(F1)] <- NA
  lapply(apply(F1, 2, table), sort, decreasing = TRUE)

  AA <- setdiff(unique(GENETIC_CODE), "*")

  i=1; G=list()
  for (i in 1:(ncol(F1) - 1)) {
    k <- min(ncol(F1), (i + 1))
    a2 <- factor(as.vector(F1[, i:k]))
    aa.dict <- levels(a2)
    names(aa.dict) <- AA[1:length(levels(a2))]
    levels(a2) <- AA[1:length(levels(a2))]
    F2 <- matrix(a2, nrow(F1))
    F2[is.na(F2)] <- "-"

    ff <- apply(F2, 1, paste, collapse = "")
    names(ff) <- paste0("tt", 1:nrow(F2))
    aln <- as.character(msa::msa(ff, type = "protein"))
    aln <- aln[names(ff)]
    F.aln <- do.call(rbind, sapply(aln, strsplit, split = ""))
    F.aln2 <- apply(F.aln, 2, function(x) aa.dict[x])

    if (ncol(F.aln2) > ncol(F2) && i < (ncol(F1) - 1)) {
      G[[i]] <- F.aln2[, 1:(ncol(F.aln) - 1)]
    } else if (i == (ncol(F1) - 1)) {
      G[[i]] <- F.aln2
    } else {
      G[[i]] <- F.aln2[, 1]
    }
  }
  G <- do.call(cbind, G)
  G <- G[, colMeans(is.na(G)) < 1, drop = FALSE]
  rownames(G) <- NULL
  colnames(G) <- paste0("V", 1:ncol(G))

  return(G)

}


#' @title Extract gene symbols from GEO feature data
#' @param fdata The featureData table from a GEOquery GEO dataset object.
#' @description Extracts official gene symbols from feature datatable of a GEO dataset downloaded with GEOquery.
#' It first looks for a column containing gene symbols by matching against the org.Hs.egSYMBOL database.
#' If no direct symbol column is found, it looks for an ENTREZ identifier and maps to symbols using org.Hs.egSYMBOL.
#' Then it looks at REFSEQ identifiers. Then it looks at ENSEMBLE IDs. If no approach works, it returns NULL.
#' @return A character vector of gene symbols, or NULL if symbols could not be extracted.
#' @export
pgx.getSymbolFromFeatureData <- function(fdata) {

  symbol <- NULL
  
  ## SYMBOL column
  SYMBOL <- as.character(unlist(as.list(org.Hs.eg.db::org.Hs.egSYMBOL)))
  symbol.col <- grep("symbol|gene|hugo", colnames(fdata), ignore.case = TRUE)
  if (any(symbol.col)) {
    ok.symbol <- apply(
      fdata[, symbol.col, drop = FALSE], 2,
      function(g) mean(toupper(g[!is.na(g)]) %in% SYMBOL)
    )
    if (any(ok.symbol > 0.5)) {
      k <- which.max(ok.symbol)
      symbol <- fdata[, symbol.col[k]]
      message("[pgx.getSymbolFromFeatureData] SYMBOL column found. Returning gene symbols...")
      return(symbol)
    }
  }

  ## ENTREZ column
  ENTREZ <- biomaRt::keys(org.Hs.eg.db::org.Hs.egSYMBOL)
  entrez.col <- grep("entrez", colnames(fdata), ignore.case = TRUE)
  if (any(entrez.col)) {
    entrez.match <- apply(
      fdata[, entrez.col, drop = FALSE], 2,
      function(g) mean(g[!is.na(g)] %in% ENTREZ)
    )
    entrez.ok <- length(entrez.col) && entrez.match > 0.5
    if (entrez.ok) {
      k <- entrez.col[which.max(entrez.match)]
      probes <- as.character(fdata[, k])
      symbol <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, probes, "SYMBOL", "ENTREZID")
      message("[pgx.getSymbolFromFeatureData] ENTREZ column found. Returning gene symbols...")
      return(symbol)
    }
  }

  ## REFSEQ column
  REFSEQ <- unlist(as.list(org.Hs.eg.db::org.Hs.egREFSEQ))
  refseq.col <- grep("refseq", colnames(fdata), ignore.case = TRUE)
  if (any(refseq.col)) {
    refseq.match <- apply(
      fdata[, refseq.col, drop = FALSE], 2,
      function(g) mean(sub("[.].*", "", g[!is.na(g)]) %in% REFSEQ)
    )
    refseq.ok <- length(refseq.col) && refseq.match > 0.5
    if (refseq.ok) {
      k <- refseq.col[which.max(refseq.match)]
      probes <- sub("[.].*", "", as.character(fdata[, k]))
      symbol <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, probes, "SYMBOL", "REFSEQ")
      message("[pgx.getSymbolFromFeatureData] REFSEQ column found. Returning gene symbols...")
      return(symbol)
    }
  }

  ## EnsembleID column
  gene.column <- grep("gene|mrna|transcript", colnames(fdata), ignore.case = TRUE)
  has.ens <- apply(fdata[, gene.column, drop = FALSE], 2, function(s) mean(grepl("ENS", s)))
  if (any(has.ens > 0.3)) {
    ens.col <- ifelse(max(has.ens) > 0, names(which.max(has.ens)), NA)
    ens.ann <- lapply(fdata[, ens.col], function(a) trimws(strsplit(a, split = "//|///")[[1]]))
    ens.probes <- sapply(ens.ann, function(s) Matrix::head(grep("^ENS", s, value = TRUE), 1))
    ens.probes[sapply(ens.probes, length) == 0] <- NA
    ens.probes <- unlist(ens.probes)
    symbol <- probe2symbol(ens.probes)
    message("[pgx.getSymbolFromFeatureData] ENSEMBLE ID column found. Returning gene symbols...")
    return(symbol)
  }

  message("WARNING:: could not parse symbol information from featureData!")
  return(NULL)

}
