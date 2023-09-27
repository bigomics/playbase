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
  ##
  ## Highly automagic download of GEO datasets from different
  ## sources with automatic probe/gene conversion and creating
  ## autocontrasts. The GEO series is first searched in a locally
  ## stored ARCHS4 H5 file, then if it is available at the recount
  ## database, if not it is retrieved from GEO using geoquery.
  ##
  ## id:      GEO id
  ## return:  object with counts, samples, genes.
  ##


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
  sampleinfo <- pgx.discretizePhenotypeMatrix(
    meta1,
    min.ncat = 2, max.ncat = 20, remove.dup = TRUE
  )
  sampleinfo <- data.frame(sampleinfo, stringsAsFactors = FALSE, check.names = FALSE)

  ## automagically create contrast matrix
  contrasts <- NULL
  if (NCOL(sampleinfo) > 0) {
    mingrp <- 3
    slen <- 15
    ref <- NA
    ct <- pgx.makeAutoContrasts(sampleinfo, mingrp = 3, slen = 20, ref = NA)
    is.null(ct)
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


#' @describeIn pgx.getGEOseries Download and process count data from GEO into a PGX object.
#' It takes a GEO accession as input, downloads the count matrix and sample metadata, processes
#' it into a PGX object with gene annotation, sample info, and contrasts.
#' @export
pgx.getGEOcounts <- function(id, archs.h5) {
  expr <- NULL
  src <- ""

  if (!is.null(archs.h5) && is.null(expr)) {
    expr <- pgx.getGEOcounts.archs4(id, archs.h5)
    if (!is.null(expr)) src <- "ARCHS4"
  }

  if (is.null(expr)) {
    expr <- pgx.getGEOcounts.recount(id)
    if (!is.null(expr)) src <- "recount"
  }

  if (is.null(expr)) {
    ## Try with GEOquery
    expr <- pgx.getGEOcounts.GEOquery(id)
    if (!is.null(expr)) src <- "GEO"
  }

  is.null(expr)
  if (is.null(expr)) {
    cat("WARNING:: could not get GEO expression. please download manually.\n")
    return(NULL)
  }
  list(expr = expr, source = src)
}


#' @describeIn pgx.getGEOseries Download and extract the metadata from a GEO series.
#' @export
pgx.getGeoMetadata <- function(id) {
  ##
  ## load series and platform data from GEO
  ##

  id
  ## First try without downloading the GSEMatrix
  pheno <- NULL
  pheno <- pgx.getGeoMetadata.fromGSM(id)
  is.null(pheno)
  if (is.null(pheno)) {
    ## try from Eset
    pheno <- pgx.getGeoMetadata.fromEset(id)
  }

  ## Sometimes the phenotype is coded in the title string
  has.title <- "title" %in% colnames(pheno)
  if (has.title && NCOL(pheno) == 0) {
    px <- title2pheno(pheno$title, split = NULL, trim = TRUE, summarize = TRUE)
    if (!is.null(px) && NCOL(px) > 0 && is.null(pheno)) {
      pheno <- px
    }
    if (!is.null(px) && NCOL(px) > 0 && !is.null(pheno)) {
      pheno <- cbind(pheno, px)
    }
  }

  colnames(pheno) <- gsub("[ ]", "_", colnames(pheno)) ## no spaces???

  return(pheno)
}


## -------------------------------------------------------------------------------------
## Query GEO expression
## -------------------------------------------------------------------------------------


#' @describeIn pgx.getGEOseries Downloads and extracts gene expression count
#' data from a GEO series stored in an HDF5 file. It searches the HDF5 file
#' metadata to find samples matching the input GEO series ID, and returns the
#' count matrix for those samples.
#' @export
pgx.getGEOcounts.archs4 <- function(id, h5.file) {
  rhdf5::h5ls(h5.file)
  sample.series <- rhdf5::h5read(h5.file, "meta/Sample_series_id")
  sample.series <- strsplit(as.character(sample.series), split = "Xx-xX")

  gse.series <- sort(unique(unlist(sample.series)))
  id %in% gse.series

  idx <- which(sapply(sample.series, function(s) id %in% s))

  if (length(idx) == 0) {
    cat("WARNING: series", id, "not in ARCHS4 matrix file\n")
    return(NULL)
  }
  X <- rhdf5::h5read(h5.file, "data/expression", index = list(NULL, idx))

  sample.acc <- rhdf5::h5read(h5.file, "meta/Sample_geo_accession")
  gene_name <- rhdf5::h5read(h5.file, "meta/genes")
  colnames(X) <- sample.acc[idx]
  rownames(X) <- gene_name

  ## collapse by symbol
  jj <- !is.na(gene_name) & gene_name != ""
  X <- X[jj, ]
  gene_name <- gene_name[jj]
  ## sum intensities (linear)
  X2 <- tapply(1:nrow(X), gene_name, function(ii) {
    Matrix::colSums(X[ii, , drop = FALSE], na.rm = TRUE) ## not log!!
  })
  X2 <- do.call(rbind, X2)

  return(X2)
}


#' @describeIn pgx.getGEOseries Downloads and processes gene-level count data for a GEO
#' series from the recount database. It takes a GEO accession ID, searches
#' recount, downloads the RangedSummarizedExperiment object, and returns the count matrix.
#' @export
pgx.getGEOcounts.recount <- function(id) {
  ## Vignette recount-quickstart.html
  ## Load library


  ## Find a project of interest
  project_info <- recount::abstract_search(id)
  project_info$project

  if (length(project_info$project) == 0) {
    cat("could not find", id, "in recount database\n")
    return(NULL)
  }

  ## Download the gene-level RangedSummarizedExperiment data
  outdir <- file.path(tempdir(), project_info$project)
  recount::download_study(project_info$project, outdir = outdir)

  ## Load the data
  load(file.path(outdir, "rse_gene.Rdata"))

  ## Scale counts by taking into account the total coverage per sample
  rse <- recount::scale_counts(rse_gene)

  counts <- MultiAssayExperiment::assay(rse)

  return(counts)
}


#' @describeIn  pgx.getGEOseries Retrieves gene expression count data for a
#' specified GEO accession ID using the GEOquery package. It downloads the
#' series matrix data, platform metadata, and probe annotations from GEO into
#' R objects.
#' @export
pgx.getGEOcounts.GEOquery <- function(id) {
  ## Retrieve expression matrix, phenotype and probe annotation
  ## matrices for a certain GEO id.
  ##


  ## load series and platform data from GEO
  id

  gse <- try(GEOquery::getGEO(id, GSEMatrix = TRUE, getGPL = TRUE))

  class(gse)
  if (inherits(gse, "try-error")) {
    cat("ERROR: GEOquery::getGEO() error\n")
    return(NULL)
  }

  # The f exprs comes from unknown source, but it is used
  has.expr <- sapply(gse, function(x) nrow(exprs(x)) > 0)
  has.expr

  if (!any(has.expr)) {
    cat("WARNING: dataset has no included expression data\n")
    supp_file <- sapply(gse, function(g) g@experimentData@other$supplementary_file)
    supp_file
    if (any(nchar(supp_file)) > 5) {
      cat("Supplementary file available: ", paste(supp_file, collapse = " "), "\n")
    }
    return(NULL)
  }

  ## select which has expression
  gse <- gse[which(has.expr)]

  ## select preferred platform is multiple exists
  expr.list <- list()
  k <- 1
  for (k in 1:length(gse)) {
    ## get expression
    eset <- gse[[k]]
    ex <- exprs(eset)

    if (ncol(ex) <= 3) {
      ## too small dataset
      next()
    }

    ## perform linear transformation (unlog) if required
    qx <- as.numeric(stats::quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = T))
    qx
    is.count <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0) ||
      (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2) ## from GEO2R script
    is.count
    if (!is.count) {
      ex <- 2**ex
    }

    ## featuredata
    has.fdata <- !is.null(fData(eset)) && NCOL(fData(eset)) > 0
    has.fdata
    if (has.fdata) {
      fdata <- fData(eset)
    } else {
      eset@annotation
      gpl.annot <- GEOquery::getGEO(eset@annotation)
      fdata <- GEOquery::Table(gpl.annot)
      fdata <- fdata[match(rownames(ex), fdata$ID), ]
    }

    ## get symbol from featuredata
    fsymbol <- pgx.getSymbolFromFeatureData(fdata)

    ## collapse by symbol
    jj <- which(!is.na(fsymbol) & fsymbol != "")
    ex <- ex[jj, ]
    fsymbol <- fsymbol[jj]
    ## sum intensities (linear)
    ex2 <- tapply(1:nrow(ex), fsymbol, function(ii) {
      Matrix::colSums(ex[ii, , drop = FALSE], na.rm = TRUE) ## not log!!
    })
    ex2 <- do.call(rbind, ex2)
    expr.list[[names(gse)[k]]] <- ex2
  }

  if (length(expr.list) == 0) {
    return(NULL)
  }

  if (length(expr.list) > 1) {
    ## merge/join all expressions
    probes <- sort(unique(unlist(lapply(expr.list, rownames))))
    samples <- sort(unique(unlist(lapply(expr.list, colnames))))
    expr.list2 <- lapply(expr.list, function(x) {
      x[match(probes, rownames(x)), match(samples, colnames(x))]
    })
    expr.list2 <- lapply(expr.list2, function(x) {
      x[is.na(x)] <- 0
      x
    })
    expr <- Reduce("+", expr.list2)
    colnames(expr) <- samples
    rownames(expr) <- probes
  } else {
    expr <- expr.list[[1]]
  }

  return(expr) ## return always linear intensities
}


## -------------------------------------------------------------------------------------
## Query GEO metadata
## -------------------------------------------------------------------------------------


#' @describeIn pgx.getGEOseries Retrieves basic experiment metadata for a GEO accession ID using GEOquery.
#' @export
pgx.getGeoExperimentInfo <- function(id) {
  suppressMessages(gse <- try(GEOquery::getGEO(id, GSEMatrix = FALSE, getGPL = FALSE)))
  info <- gse@header

  return(info)
}


#' @describeIn pgx.getGEOseries retrieves sample metadata from a GEO dataset by first downloading the
#' GSEMatrix and then extracting the phenotype data.
#' @export
pgx.getGeoMetadata.fromEset <- function(id) {
  ## If not succesful, try with downloading the GSEMatrix
  suppressMessages(gse <- try(GEOquery::getGEO(id, GSEMatrix = TRUE, getGPL = FALSE)))

  if (inherits(gse, "try-error")) {
    res <- list(error = "ERROR: pgx.getGeoMetadata.fromEset() error")
    return(res)
  }

  nsamples <- sapply(gse, function(s) nrow(pData(phenoData(s))))

  gse <- gse[nsamples >= 3]
  eset <- gse[[1]]
  pheno.list <- lapply(gse, pgx.getGeoMetadata.fromEset1)
  pheno.list <- lapply(pheno.list, function(m) {
    rn <- rownames(m)
    m1 <- as.matrix(apply(m, 2, as.character))
    rownames(m1) <- rn
    m1
  })

  vars <- sort(unique(unlist(lapply(pheno.list, colnames))))
  pheno.list <- lapply(pheno.list, function(m) {
    m1 <- m[, match(vars, colnames(m))]
    colnames(m1) <- vars
    m1
  })

  gpl <- sub("_series.*", "", sub(".*-GPL", "GPL", names(gse)))
  pheno.list <- lapply(1:length(pheno.list), function(i) {
    pheno.list[[i]] <- cbind(GPL = gpl[i], pheno.list[[i]])
  })

  pheno <- do.call(rbind, pheno.list)
  pheno <- data.frame(pheno, stringsAsFactors = FALSE, check.names = FALSE)

  return(pheno)
}


#' @describeIn pgx.getGEOseries helper function that extracts the phenotype data
#' from a GEOquery GSEMatrix object.
#' @param eset Eset object
#' @export
pgx.getGeoMetadata.fromEset1 <- function(eset) {
  ##
  ## load series and platform data from GEO
  ##

  ## Get summary
  summary <- experimentData(eset)@abstract

  ## pdata object
  pdata <- pData(phenoData(eset))

  gsm.title <- as.character(pdata$title)
  gsm.source <- as.character(pdata$source_name_ch1)
  gsm.samples <- as.character(pdata$geo_accession)

  ## Base sample_info from characteristics (ch1) column
  ch1_info <- eset.getCH1(eset)

  ## We can get extra information from title
  is.underscored <- length(gsm.title) && all(grepl("_", gsm.title))

  title_info <- NULL
  if (FALSE && is.underscored) {
    title2 <- trimws(gsm.title)
    title_info <- eset.parsePhenoFromTitle(title2, split = "_")
  }

  ## All sample_info: from characterisctis_ch1 and title
  sample_info <- data.frame(
    GSM = gsm.samples, title = gsm.title,
    source = gsm.source,
    stringsAsFactors = FALSE, check.names = FALSE
  )
  if (!is.null(ch1_info)) sample_info <- cbind(sample_info, ch1_info)
  if (!is.null(title_info)) sample_info <- cbind(sample_info, title_info)

  sample_info <- data.frame(sample_info, stringsAsFactors = FALSE, check.names = FALSE)

  sample_info
}


#' @describeIn pgx.getGEOseries retrieves sample metadata for a GEO sample accession
#' ID by downloading the full GEO series and extracting the metadata for that sample.
#' @export
pgx.getGeoMetadata.fromGSM <- function(id) {
  ##
  ## load series and platform data from GEO
  ##
  id
  suppressMessages(gse <- try(GEOquery::getGEO(id, GSEMatrix = FALSE, getGPL = FALSE)))

  if (inherits(gse, "try-error")) {
    res <- list(error = "ERROR: GEOquery::getGEO() error")
    return(res)
  }


  if (length(gse@gsms) == 0) {
    cat("WARNING:: no GSM information in object\n")
    return(NULL)
  }

  ## Summary and sample names
  summary <- gse@header$summary
  gsm.samples <- gse@header$sample_id

  ## Get sample_info from characteristics (ch1) column
  ch1_info <- lapply(gse@gsms, function(g) g@header$characteristics_ch1)

  is.null(ch1_info)
  if (!is.null(ch1_info)) {
    ch1_info <- lapply(ch1_info, function(x) sub("^Clinical info: ", "", x))
    ch1_vars <- unique(unlist(lapply(ch1_info, function(x) trimws(sub("[:=].*", "", x)))))
    ch1_info <- lapply(ch1_info, function(x) {
      xvar <- trimws(sub("[:=].*", "", x))
      x <- trimws(sub(".*[:=] ", "", x))
      names(x) <- xvar
      x <- x[match(ch1_vars, names(x))]
      x
    })
    ch1_info <- do.call(rbind, ch1_info)
    colnames(ch1_info) <- ch1_vars
  }

  ## We can get more information from title??
  gsm.title <- sapply(gse@gsms, function(g) g@header$title)
  gsm.source <- sapply(gse@gsms, function(g) g@header$source_name_ch1)
  gsm.gpl <- sapply(gse@gsms, function(g) g@header$platform_id)
  is.underscored <- length(gsm.title) && all(grepl("_", gsm.title))
  title_info <- NULL
  ## NEED RETHINK!!!!!!!!!!!!!!!!!!
  if (FALSE && is.underscored) {
    title2 <- trimws(gsm.title)
    title_info <- eset.parsePhenoFromTitle(title2, split = "_")
  }

  ## All sample_info: from characterisctis_ch1 and title
  sample_info <- data.frame(
    GPL = gsm.gpl, GSM = gsm.samples, title = gsm.title,
    source = gsm.source, stringsAsFactors = FALSE
  )
  if (!is.null(ch1_info)) sample_info <- cbind(sample_info, ch1_info)
  if (!is.null(title_info)) sample_info <- cbind(sample_info, title_info)

  sample_info <- data.frame(sample_info, stringsAsFactors = FALSE, check.names = FALSE)

  return(sample_info)
}




## -------------------------------------------------------------------------------------
## HELPER functions
## -------------------------------------------------------------------------------------


#' Extract phenotype data field from ExpressionSet
#'
#' @param eset ExpressionSet object
#' @param field Character string specifying phenotype data field name
#'
#' @return Vector of values for the specified field
#'
#' @details This function extracts a specific phenotype data field from an ExpressionSet object.
#' The \code{field} parameter specifies the phenotype data column name to extract.
#'
#' The phenotype data pData slot is extracted from the ExpressionSet, and the specified
#' field is returned as a vector.
#'
#' @export
eset.getPhenoData <- function(eset, field) {
  pData(phenoData(eset))[, field]
}




#' @describeIn eset.getPhenoData Extracts sample phenotype terms from GEO dataset titles by splitting on a delimiter.
#' @param title character vector with the titles
#' @param split delimiter character to split on `c(",", ";", "\\|", "_", " ")`.
#' @param trim trim leading and trailing whitespace from each term
#' @param summarize summarize the terms by counting the number of occurrences of each term
#' @export
title2pheno <- function(title, split = NULL, trim = TRUE, summarize = TRUE) {
  ## determine the split character
  if (is.null(split)) {
    split.chars <- c(",", ";", "\\|", "_", " ")
    ss <- c()
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
  split

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
  F1
  F1[is.na(F1)] <- NA
  rownames(F1) <- NULL

  ## Guess names
  getmax.term <- function(s) {
    tt <- table(unlist(strsplit(s, split = "[ _]")))
    vip.tt <- grep("hour|repl|hr|time|treat|infec|pati|sampl", names(tt))
    vip.tt
    if (length(vip.tt)) tt[vip.tt] <- 1.1 * tt[vip.tt] ## boost known keywords
    names(which.max(tt))
  }
  maxterm <- apply(F1, 2, function(s) getmax.term(s))
  maxterm <- paste0("_", maxterm)
  maxterm
  colnames(F1) <- maxterm

  ## trims same words/characters on both ends
  if (trim) F1 <- apply(F1, 2, trimsame, summarize = summarize)
  F1
}


#' @describeIn eset.getPhenoData Phenotype data from the title of an
#' ExpressionSet object by splitting the title on a specified delimiter and guessing column names.
#' @export
eset.parsePhenoFromTitle <- function(title, split = NULL) {
  if (!all(grepl(split, title))) {
    return(NULL)
  }

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
  F1
  F1[is.na(F1)] <- NA
  lapply(apply(F1, 2, table), sort, decreasing = TRUE)

  AA <- setdiff(unique(GENETIC_CODE), "*")

  G <- list()
  i <- 1
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
