##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

#' Create a pgx object
#'
#' This function creates a pgx object from files, which is the core object in the
#' OmicsPlayground. It then runs the specified differential expression methods.
#'
#' @param counts.file Path to counts data file. Rows are genes, columns are samples.
#' @param samples.file Path to samples data file. Rows are samples, columns are sample info.
#' @param contrasts.file (optional) Path to contrasts file. Rows and columns define contrasts.
#' @param gxmethods a string with the gene-level methods to use. The default value is \code{"trend.limma,edger.qlf,deseq2.wald"}
#' @param gsetmethods a string with the gene-set methods to use. The default value is \code{"fisher,gsva,fgsea"}
#' @param extra a string with the extra modules to use. The default value is \code{"meta.go,deconv,infer,drugs,wordcloud"}
#'
#' @return list. represents a pgx object. It contains the data and analysis results.
#' @examples
#' \dontrun{
#'
#' library(playbase)
#' counts <- system.file("extdata", "counts.csv", package = "playbase")
#' contrasts <- system.file("extdata", "contrasts.csv", package = "playbase")
#' samples <- system.file("extdata", "samples.csv", package = "playbase")
#'
#' mypgx <- pgx.createFromFiles(counts, samples, contrasts)
#' }
#' @export
pgx.createFromFiles <- function(counts.file,
                                samples.file,
                                contrasts.file = NULL,
                                gxmethods = "trend.limma,edger.qlf,deseq2.wald",
                                gsetmethods = "fisher,gsva,fgsea",
                                extra = "meta.go,deconv,infer,drugs,wordcloud",
                                pgx.dir = "./data",
                                libx.dir = "./libx") {
  ## read counts table (allow dup rownames)
  counts <- read.as_matrix(counts.file)

  ## compile sample table
  samples <- read.as_matrix(samples.file)
  samples <- data.frame(samples, check.names = FALSE)

  ## parse requested phenotypes
  if (!is.null(contrasts.file) && file.exists(contrasts.file)) {
    message("reading contrasts file ", contrasts.file)
    contrasts <- read.as_matrix(contrasts.file)
  } else {
    ## take first (not-dotted) column in samples as phenotype vector
    group.col <- head(grep("group|condition", colnames(samples), ignore.case = TRUE), 1)
    if (length(group.col) == 0) {
      group.col <- head(grep("^.*", colnames(samples), invert = TRUE), 1)
    }
    if (length(group.col) == 0) {
      group.col <- colnames(samples)[1]
    }
    Y <- samples[, group.col, drop = FALSE]
    ## automatically guess contrasts
    contr <- pgx.makeAutoContrasts(Y, mingrp = 3, slen = 20, ref = NA)
    contrasts <- contrastAsLabels(contr$exp.matrix)
  }

  ## other params
  gx.methods <- strsplit(gxmethods, split = ",")[[1]]
  gset.methods <- strsplit(gsetmethods, split = ",")[[1]]
  extra.methods <- strsplit(extra, split = ",")[[1]]

  ## create initial PGX object
  pgx <- pgx.createPGX(
    counts,
    samples = samples,
    contrasts = contrasts,
    X = NULL,
    is.logx = NULL,
    dotimeseries = FALSE,
    batch.correct.method = "no_batch_correct",
    batch.pars = "<autodetect>",
    covariates = NULL,
    auto.scale = TRUE,
    filter.genes = TRUE,
    prune.samples = FALSE,
    only.known = TRUE,
    average.duplicated = FALSE,
    only.hugo = TRUE,
    convert.hugo = TRUE,
    only.proteincoding = TRUE,
    max.genesets = 10000
  )

  ## start computing PGX object
  pgx <- pgx.computePGX(
    pgx,
    max.genes = 40000,
    gx.methods = gx.methods,
    gset.methods = gset.methods,
    extra.methods = extra.methods,
    cluster.contrasts = FALSE,
    do.clustergenes = TRUE,
    do.clustergenesets = TRUE,
    do.cluster = TRUE,
    use.design = FALSE,
    prune.samples = TRUE,
    pgx.dir = pgx.dir,
    libx.dir = libx.dir,
    progress = NULL
  )

  ## save
  pgx
}

#' Create a PGX object
#' This function creates a pgx object, which is the core object in the
#' OmicsPlayground.
#' @param counts Matrix of count data with genes as rows and samples as columns.
#' @param samples Data frame containing sample information.
#' @param organism Default "Human", it indicates the species used
#' for the gene annotation table and the probe to symbol conversion.
#' @param contrasts Data frame defining sample contrasts.
#' @param X (Optional) Matrix of normalized expression data. If NULL, will be calculated from counts.
#' @param is.logx Logical indicating if count matrix is already log-transformed. If NULL, guessed automatically.
#' @param dotimeseries Logical indicating if timeseries analysis has been activated by the user at upload
#' @param batch.correct.method BC method. Default is "no_batch_correct" (meaning no batch correction).
#' @param batch.pars BC variable. Default "autodetect" as per QC/BC tab in upload.
#' @param covariates variables to regress out. Valid only for linear model-based tests.
#' @param auto.scale Logical indicating whether to automatically scale/center genes. Default is TRUE.
#' @param filter.genes Logical indicating whether to filter lowly expressed genes. Default is TRUE.
#' @param prune.samples Logical indicating whether to remove samples without contrasts. Default is FALSE.
#' @param only.known Logical indicating whether to keep only known genes. Default is TRUE.
#' @param average.duplicated Logical whether average duplicated features (if any). Default FALSE (thus keep all by making unique).
#' @param only.hugo Logical indicating whether to convert symbols to HUGO names. Default is TRUE.
#' @param convert.hugo Logical indicating whether to convert symbols to HUGO names. Default is TRUE.
#' @param only.proteincoding Logical indicating whether to keep only protein-coding genes. Default is TRUE.
#' @param custom.geneset Custom gene sets to test, as a named list with gmt and info elements.
#' @param max.genesets Maximum number of gene sets to test. Default is 5000.
#'
#' @details
#' pgx.createPGX creates a pgx object with the following slots:
#'
#' - `name`: Name of the dataset
#' - `organism`: Organism for the dataset
#' - `version`: Dataset version
#' - `date`: Date the dataset was created
#' - `creator`: Creator of the dataset
#' - `datatype`: Type of data (e.g. RNA-seq, microarray)
#' - `description`: Description of the dataset
#' - `samples`: Sample metadata
#' - `counts`: Raw count matrix
#' - `contrasts`: Contrast matrix
#' - `X`: Normalized expression matrix
#' - `total_counts`: Total counts per sample
#' - `counts_multiplier`: Counts multiplier for each sample
#' - `genes`: Gene annotation data.frame (initially NULL)
#' - `all_genes`: Full list of genes
#' - `probe_type`: Probe type according to biomaRt classification(e.g. ensemble_id)
#' - `GMT`: Gene set matrix
#' @import data.table
#' @return List. PGX object containing input data and parameters.
#'
#' @export
pgx.createPGX <- function(counts,
                          samples,
                          contrasts,
                          organism,
                          custom.geneset = NULL,
                          annot_table = NULL,
                          max.genesets = 5000,
                          name = "Data set",
                          datatype = "RNA-seq",
                          azimuth_ref = "pbmcref",
                          probe_type = NULL,
                          creator = "unknown",
                          description = "No description provided.",
                          X = NULL,
                          norm_method = "CPM",
                          is.logx = NULL,
                          dotimeseries = FALSE,
                          batch.correct.method = "no_batch_correct",
                          batch.pars = "<autodetect>",
                          covariates = NULL, ## new
                          auto.scale = TRUE,
                          filter.genes = TRUE,
                          exclude.genes = NULL,
                          prune.samples = FALSE,
                          only.known = TRUE,
                          average.duplicated = FALSE,
                          only.hugo = TRUE, ## DEPRECATED
                          convert.hugo = FALSE,
                          only.proteincoding = TRUE,
                          remove.xxl = TRUE, ## DEPRECATED
                          remove.outliers = TRUE, ## DEPRECATED
                          add.gmt = TRUE,
                          settings = list(),
                          sc_compute_settings = list()) {

  message("[pgx.createPGX]===========================================")
  message("[pgx.createPGX]=========== pgx.createPGX =================")
  message("[pgx.createPGX]===========================================")
  message("\n")
  message("[pgx.createPGX] datatype = ", datatype, "\n")

  if (is.null(counts)) stop("[pgx.createPGX] FATAL: counts must be provided")
  if (is.null(samples)) stop("[pgx.createPGX] FATAL: samples must be provided")
  if (is.null(organism)) stop("[pgx.createPGX] FATAL: organism must be provided")

  message("[pgx.createPGX] dim.counts: ", dim(counts)[1], " x ", dim(counts)[2])
  message("[pgx.createPGX] class.counts: ", class(counts))
  message("[pgx.createPGX] counts has ", sum(is.na(counts)), " missing values")
  
  ndup <- sum(duplicated(rownames(counts)))
  if (ndup > 0) {
    if (average.duplicated) {
      message("[pgx.createPGX] ", ndup, " duplicated feature(s) detected. Averaging....")
      counts <- playbase::counts.mergeDuplicateFeatures(counts, is.counts = TRUE)
      if (!is.null(X)) X <- playbase::counts.mergeDuplicateFeatures(X, is.counts = FALSE)
    } else {
      message("[pgx.createPGX] ", ndup, " duplicated feature(s) detected. Making unique to keep all...")
      rownames(counts) <- playbase::make_unique(rownames(counts))
      if (!is.null(X)) rownames(X) <- playbase::make_unique(rownames(X))
      if (!is.null(annot_table)) rownames(annot_table) <- rownames(counts)
    }
  }

  if (datatype == "scRNA-seq") {
    pgx <- pgx.createSingleCellPGX(
      counts = counts,
      samples = samples,
      contrasts = contrasts,
      organism = organism,
      azimuth_ref = azimuth_ref,
      batch = NULL,
      sc_compute_settings = sc_compute_settings
    )
    return(pgx)
  }

  if (is.null(X)) {
    min.nz <- min(counts[counts > 0], na.rm = TRUE)
    prior <- ifelse(grepl("CPM|TMM|TPM", norm_method), 1, min.nz)
    message("[pgx.createPGX] creating X as log2(counts+p) with p = ", prior)
    X <- log2(counts + prior)
  }

  if (!is.null(annot_table)) {
    message("[pgx.createPGX] dim(annot_table) = ", nrow(annot_table), " x ", ncol(annot_table))
    ndiff <- sum(rownames(annot_table) != rownames(counts))
    message("[pgx.createPGX] WARNING: annot_table has ", ndiff, " different rownames as counts")
    ndups <- sum(duplicated(rownames(annot_table)))
    message("[pgx.createPGX] annot_table has ", ndups, " duplicated rows")
    if (nrow(annot_table) != nrow(counts)) {
      message("[pgx.createPGX] WARNING: annot_table has different nrows. forcing dimensions.")
      ii <- match(rownames(counts), rownames(annot_table))
      annot_table <- annot_table[ii, ]
      rownames(annot_table) <- rownames(counts)
    }
  }

  if (sum(is.na(X)) > 0) {
    message("[pgx.createPGX] X has ", sum(is.na(X)), " missing values")
  }

  if (!is.null(X) && !all(dim(counts) == dim(X))) {
    stop("[pgx.createPGX] dimension of counts and X do not match\n")
  }

  if (!all(rownames(counts) == rownames(X))) {
    stop("rownames of counts and X do not match\n")
  }

  if (datatype == "multi-omics") {
    has.colons <- mean(grepl("[:]", rownames(counts)), na.rm = TRUE) > 0.9
    if (!has.colons) stop("[pgx.createPGX] FATAL: features must have multi-omics prefix\n")
  }

  ## -------------------------------------------------------------------
  ## clean up input files
  ## -------------------------------------------------------------------
  samples <- as.data.frame(samples, drop = FALSE)
  counts <- as.matrix(counts)
  X <- as.matrix(X)
  if (is.null(contrasts)) contrasts <- samples[, 0]

  ## convert old-style contrast matrix to sample-wise labeled contrasts
  contrasts <- contrasts.convertToLabelMatrix(contrasts, samples)
  contrasts <- fixContrastMatrix(contrasts)

  ## ---------------------------------------------------------------------
  ## Time series conducted if user checked the box during upload
  ## ---------------------------------------------------------------------
  if (dotimeseries) contrasts <- contrasts.addTimeInteraction(contrasts, samples)
  
  ## -------------------------------------------------------------------
  ## Auto-scaling (scale down huge values, often in proteomics)
  ## -------------------------------------------------------------------
  # res <- counts.autoScaling(counts)
  # counts <- res$counts
  # counts_multiplier <- res$counts_multiplier
  counts_multiplier <- Inf
  # remove(res)

  ## -------------------------------------------------------------------
  ## conform all matrices
  ## -------------------------------------------------------------------
  message("[createPGX] conforming matrices...")

  ## prune unused samples
  contrasts[contrasts %in% c("", " ", "NA")] <- NA
  used.samples <- names(which(rowSums(!is.na(contrasts)) > 0))
  if (prune.samples && length(used.samples) < ncol(counts)) {
    counts <- counts[, used.samples, drop = FALSE]
    samples <- samples[used.samples, , drop = FALSE]
    contrasts <- contrasts[used.samples, , drop = FALSE] ## sample-based!!!
  }

  ## align samples
  kk <- intersect(colnames(counts), rownames(samples))
  kk <- intersect(kk, colnames(X))
  counts <- counts[, kk, drop = FALSE]
  X <- X[, kk, drop = FALSE]
  samples <- samples[kk, , drop = FALSE]
  samples <- utils::type.convert(samples, as.is = TRUE) ## automatic type conversion
  if (all(kk %in% rownames(contrasts))) {
    contrasts <- contrasts[kk, , drop = FALSE]
  }

  ## sanity checks
  if (ncol(X) == 0) {
    info("[createPGX] FATAL. ncol(X) == 0")
    return(NULL)
  }

  ## -------------------------------------------------------------------
  ## Special case for PTM phospho-proteomics.
  ## -------------------------------------------------------------------
  is.phospho <- annotate_phospho_residue(rownames(counts), detect.only = TRUE)
  if (datatype == "proteomics" && is.phospho) {
    info("[createPGX] annotating rownames with phospho residue...")
    newnames <- annotate_phospho_residue(rownames(counts))
    newnames <- make_unique(newnames)
    rownames(counts) <- newnames
    rownames(X) <- newnames
    if (!is.null(annot_table)) {
      rownames(annot_table) <- newnames
      pos.col <- grep("site|position|phosho", colnames(annot_table), ignore.case = TRUE)
      phosphosite <- sub(".*_|[.].*", "", newnames)
      if (length(pos.col)) {
        i <- pos.col[1]
        annot_table[, i] <- phosphosite
      } else {
        annot_table$site_position <- phosphosite
      }
    }
  }

  ## -------------------------------------------------------------------
  ## create pgx object
  ## -------------------------------------------------------------------
  message("[createPGX] creating pgx object...")

  ## remove special characters from description (other columns too??)
  description <- gsub("[\"\']", " ", description) ## remove quotes (important!!)
  description <- gsub("[\n]", ". ", description) ## replace newline
  description <- trimws(gsub("[ ]+", " ", description)) ## remove ws

  ## add to setting info
  settings$filter.genes <- filter.genes
  settings$exclude.genes <- exclude.genes
  settings$only.known <- only.known
  settings$only.proteincoding <- only.proteincoding
  settings$convert.hugo <- convert.hugo
  settings$custom.geneset <- !is.null(custom.geneset)

  ## add versions info
  versions <- list()
  versions$playbase_version <- packageVersion("playbase")
  versions$playdata_version <- packageVersion("playdata")

  pgx <- list(
    name = name,
    organism = organism,
    version = packageVersion("playbase"), # useless, just keep for back compatibility
    date = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    creator = creator,
    datatype = datatype,
    description = description,
    samples = data.frame(samples, check.names = FALSE),
    counts = as.matrix(counts),
    contrasts = contrasts,
    X = X,
    norm_method = norm_method,
    total_counts = Matrix::colSums(counts, na.rm = TRUE),
    counts_multiplier = counts_multiplier,
    covariates = covariates,
    settings = settings,
    versions = versions,
    sc_compute_settings = sc_compute_settings
  )

  ## -------------------------------------------------------------------
  ## create gene annotation table
  ## -------------------------------------------------------------------
  pgx$genes <- NULL
  pgx$probe_type <- probe_type

  message("[createPGX] annotating genes")
  pgx$genes <- getProbeAnnotation(
    organism = pgx$organism,
    probes = rownames(pgx$counts),
    datatype = pgx$datatype,
    probetype = pgx$probe_type,
    annot_table = annot_table
  )

  if (is.null(pgx$genes)) stop("[pgx.createPGX] FATAL: Could not build gene annotation")

  if (!"symbol" %in% colnames(pgx$genes) && "gene_name" %in% colnames(pgx$genes)) {
    dbg("[pgx.createPGX] WARNING! no symbol column. copying deprecated gene_name column as symbol")
    pgx$genes$symbol <- pgx$genes$gene_name
  }

  if (all(is.na(pgx$genes$symbol))) {
    dbg("[pgx.createPGX] WARNING! all symbol NA. copying rownames as symbol")
    pgx$genes$symbol <- gsub(".*:|[.].*", "", rownames(pgx$genes))
  }

  ## -------------------------------------------------------------------
  ## Filter out not-expressed
  ## -------------------------------------------------------------------
  if (filter.genes) {
    nexpr <- sum(rowSums(pgx$counts, na.rm = TRUE) == 0)
    message("[pgx.createPGX] Filtering out ", nexpr, " not-expressed genes...")
    pgx <- pgx.filterZeroCounts(pgx)
    ii <- match(rownames(pgx$counts), rownames(pgx$genes))
    pgx$genes <- pgx$genes[ii, , drop = FALSE]
  }

  ## -------------------------------------------------------------------
  ## Filter genes
  ## -------------------------------------------------------------------
  do.filter <- (only.known || only.proteincoding || !is.null(exclude.genes))
  if (do.filter) {
    if (only.known) {
      message("[pgx.createPGX] Removing genes without symbol...")
      no.symbol <- (is.na(pgx$genes$symbol) | pgx$genes$symbol %in% c("", "-"))
      pgx$genes <- pgx$genes[which(!no.symbol), ]
    }

    if (only.proteincoding) {
      message("[pgx.createPGX] Removing Rik/ORF/LOC genes...")
      is.unknown <- grepl("^rik|^loc|^orf", tolower(pgx$genes$symbol))
      is.unknown <- is.unknown & !is.na(pgx$genes$symbol)
      pgx$genes <- pgx$genes[which(!is.unknown), ]
    }

    if (!is.null(exclude.genes)) {
      message("[pgx.createPGX] excluding genes: ", exclude.genes)
      exstr <- strsplit(tolower(exclude.genes), split = "[ ,]")[[1]]
      exexpr <- paste(c(paste0("^", exstr), paste0(exstr, "$")), collapse = "|")
      exgene <- grepl(exexpr, tolower(pgx$genes$symbol))
      if (sum(exgene)) {
        pgx$genes <- pgx$genes[which(!exgene), ]
      }
    }

    ## conform
    pgx$counts <- pgx$counts[rownames(pgx$genes), , drop = FALSE]
    pgx$X <- pgx$X[rownames(pgx$genes), , drop = FALSE]
  }

  ## -------------------------------------------------------------------
  ## collapse probe-IDs to gene symbol and aggregate duplicates
  ## -------------------------------------------------------------------
  ## if feature/rownames are not symbol, we paste symbol to row name.
  pp <- sub("^[a-zA-Z]+:", "", rownames(pgx$genes))
  mean_feature_is_symbol <- mean(pp == pgx$genes$symbol, na.rm = TRUE)

  ## NOTE: this was old chunk to convert rownames to HUGO gene
  ## symbol. It now serves to append symbol to rownames/feature names.
  if (convert.hugo && mean_feature_is_symbol < 0.10) {
    symbol <- pgx$genes$symbol
    symbol[is.na(symbol)] <- ""
    feature_is_symbol <- (sub("^[a-zA-Z]+:", "", rownames(pgx$genes)) == symbol)

    new.names <- combine_feature_names(pgx$genes, target = c("rownames", "_", "symbol"))
    new.names <- ifelse(feature_is_symbol, rownames(pgx$genes), new.names)
    new.names <- make_unique(new.names)

    rownames(pgx$genes) <- new.names
    pgx$genes$gene_name <- new.names ## gene_name should also be renamed??
    pgx$genes$feature <- new.names ## feature should also be renamed??
    rownames(pgx$counts) <- new.names
    rownames(pgx$X) <- new.names
  }

  ## -------------------------------------------------------------------
  ## Infer cell cycle/gender here (before any batchcorrection)
  ## -------------------------------------------------------------------
  info("[createPGX] infer cell cycle")
  pgx <- compute_cellcycle_gender(pgx, pgx$counts)

  ## -------------------------------------------------------------------
  ## Add GMT
  ## -------------------------------------------------------------------
  ## If no organism, no custom annotation table and no custom geneset,
  ## then create empty GMT
  unknown.organism <- (pgx$organism %in% c("No organism", "custom", "unkown"))
  unknown.datatype <- (pgx$datatype %in% c("custom", "unkown"))
  no3 <- unknown.organism && is.null(annot_table) && is.null(custom.geneset)
  if (no3 || unknown.datatype || !add.gmt) {
    message("[pgx.createPGX] WARNING: empty GMT matrix. No gene sets. ")
    pgx$GMT <- Matrix::Matrix(0, nrow = 0, ncol = 0, sparse = TRUE)
  } else {
    pgx <- pgx.add_GMT(pgx = pgx, custom.geneset = custom.geneset, max.genesets = max.genesets)
  }

  ## --------------------------------
  ## rm NA contrasts
  ## --------------------------------
  if (ncol(pgx$samples) > 1) {
    pgx$samples <- pgx$samples[, colMeans(is.na(pgx$samples)) < 1, drop = FALSE]
  }

  ## -------------------------------------------------------------------
  ## Batch correction if user-selected
  ## -------------------------------------------------------------------
  if (batch.correct.method != "no_batch_correct" && nrow(pgx$samples) > 2) {
    batch <- NULL
    mm <- batch.correct.method[1]
    if (length(batch.pars) == 0) batch.pars <- "<autodetect>"
    X <- pgx$X
    samples <- pgx$samples
    contrasts <- pgx$contrasts

    message("[pgx.createPGX] batch.correct.method=", batch.correct.method)
    message("[pgx.createPGX] batch.pars=", batch.pars)

    pars <- playbase::get_model_parameters(X, samples, pheno = NULL, contrasts)
    if (any(grepl("<autodetect>", batch.pars))) batch.pars <- pars$batch.pars
    if (any(grepl("<none>", batch.pars))) batch.pars <- ""
    batch.pars <- intersect(batch.pars, colnames(samples))
    if (length(batch.pars)) batch <- samples[, batch.pars, drop = FALSE]
    pheno <- pars$pheno

    message("[pgx.createPGX] Batch correction using ", mm)
    if (sum(is.na(X)) == 0) {
      xlist <- playbase::runBatchCorrectionMethods(X, batch, pheno, methods = mm, ntop = Inf)
      cX <- xlist[[mm]]
    } else {
      is.mox <- is.multiomics(rownames(X))
      if (is.mox) {
        impX <- imputeMissing.mox(X, method = "SVD2")
      } else {
        impX <- imputeMissing(X, method = "SVD2")
      }
      xlist <- playbase::runBatchCorrectionMethods(impX, batch, pheno, methods = mm, ntop = Inf)
      cX <- xlist[[mm]]
      jj <- which(is.na(X), arr.ind = TRUE)
      cX[jj] <- NA ## Batch corrected X; original NAs restored
    }

    ## Compute correctedCounts from corrected X.
    counts <- pgx$counts ## same as the one originally uploaded by user.
    jj <- which(rownames(counts) %in% rownames(cX))
    kk <- which(colnames(counts) %in% colnames(cX))
    counts <- counts[jj, kk]
    tc.counts <- colSums(counts, na.rm = TRUE)

    prior <- 0
    if (min(counts, na.rm = TRUE) == 0 || any(is.na(counts))) {
      prior <- min(counts[counts > 0], na.rm = TRUE)
    }
    if (grepl("CPM|TMM", norm_method)) prior <- 1
    rc.counts <- pmax(2**cX - prior, 0) # recomputed counts
    tc.rc.counts <- colSums(rc.counts, na.rm = TRUE)

    ## Put back to original total counts.
    corrected.counts <- t(t(rc.counts) / tc.rc.counts * tc.counts)

    ## Restore original NAs in correctedCounts.
    jj <- which(is.na(counts), arr.ind = TRUE)
    if (any(jj)) corrected.counts[jj] <- NA

    message("[pgx.createPGX] Batch correction completed\n")

    pgx$X <- cX
    pgx$counts <- corrected.counts

    rm(xlist, cX, counts, corrected.counts)
  }

  rm(counts, X, samples, contrasts)

  message("\n\n")
  message("[pgx.createPGX]======================================")
  message("[pgx.createPGX]======== pgx.createPGX: DONE! ========")
  message("[pgx.createPGX]======================================")
  message("\n\n")

  return(pgx)
}


#' @title Compute PGX
#' @description Main function to populate pgx with results. The function computes the analysis on a pgx object
#'
#' @param pgx A pgx object containing the input data
#' @param max.genes Maximum number of genes to test. Default is 19999.
#' @param gx.methods Methods for differential expression analysis at the gene level. Default is c("ttest.welch", "trend.limma", "edger.qlf").
#' @param gset.methods Methods for differential analysis at the gene set level. Default is c("fisher", "gsva", "fgsea").
#' @param do.cluster Logical indicating whether to run sample clustering. Default is TRUE.
#' @param do.clustergenesets Logical indicating whether to cluster gene sets.
#' @param do.clustergenes Logical indicating whether to cluster genes. Default is TRUE.
#' @param use.design Whether to use model design matrix for testing. Default is FALSE.
#' @param prune.samples Whether to remove samples without valid contrasts. Default is TRUE.
#' @param time Whether perform time series analysis or not. Default FALSE
#' @param extra.methods Additional analysis methods to run. Default is c("meta.go", "infer", "deconv", "drugs", "wordcloud", "wgcna")[c(1, 2)].
#' @param libx.dir Directory containing custom analysis modules.
#' @param progress A progress object for tracking status.
#'
#' @details
#' The slots created by pgx.computePGX are the following:
#'
#' - `tsne2d`: 2D tSNE coordinates matrix
#' - `tsne3d`: 3D tSNE coordinates matrix
#' - `cluster`: List containing sample clustering results
#' - `cluster.genes`: List containing gene clustering results
#' - `model.parameters`: Model parameters from normalization
#' - `timings`: Matrix of timings for computations
#' - `gx.meta`: Gene metadata data.frame
#' - `gset.meta`: Gene set metadata data.frame
#' - `gsetX`: Gene set scores matrix
#' - `cluster.gsets`: List of gene set clustering results
#' - `meta.go`: GO graph and metadata
#'
#' @return An updated pgx object containing analysis results.
#'
#' @export
pgx.computePGX <- function(pgx,
                           max.genes = 19999,
                           gx.methods = c("trend.limma", "edger.qlf", "deseq2.wald"),
                           gset.methods = c("fisher", "gsva", "fgsea"),
                           custom.geneset = list(gmt = NULL, info = NULL),
                           custom_fc = NULL,
                           do.cluster = TRUE,
                           cluster.contrasts = FALSE,
                           do.clustergenesets = TRUE,
                           do.clustergenes = TRUE,
                           use.design = FALSE,
                           prune.samples = TRUE,
                           extra.methods = c(
                             "meta.go", "infer", "deconv", "drugs",
                             "connectivity", "wordcloud", "wgcna",
                             "mofa"
                           )[c(1, 2)],
                           pgx.dir = NULL,
                           libx.dir = NULL,
                           progress = NULL,
                           user_input_dir = getwd()) {

  message("[pgx.computePGX]===========================================")
  message("[pgx.computePGX]========== pgx.computePGX =================")
  message("[pgx.computePGX]===========================================")
  message("\n")
  message("[pgx.computePGX] Starting pgx.computePGX")
  message("\n")

  if (!"contrasts" %in% names(pgx)) {
    stop("[pgx.computePGX] FATAL:: no contrasts in object")
  }

  if (!all(grepl("_vs_", colnames(pgx$contrasts)))) {
    stop("[pgx.computePGX] FATAL:: all contrast names must include _vs_")
  }

  ## -----------------------------------------------------------------------------
  ## Time series: check methods
  ## -----------------------------------------------------------------------------
  timeseries <- any(grepl("IA:*", colnames(pgx$contrasts)))
  if (timeseries) {
    ts.mm <- c("trend.limma", "deseq2.lrt", "deseq2.wald", "edger.lrt", "edger.qlf")
    cm <- intersect(gx.methods, ts.mm)
    if (length(cm) == 0) {
      message(
        "[pgx.computePGX] For time series analysis, gx.methods must be among ",
        paste0(ts.mm, collapse = "; "), " Skipping time series analysis."
      )
      hh <- grep("IA:*", colnames(pgx$contrasts))
      pgx$contrasts <- pgx$contrasts[, -hh, drop = FALSE]
    } else {
      gx.methods <- cm
    }
  }

  contr.matrix <- contrasts.convertToLabelMatrix(pgx$contrasts, pgx$samples)
  contr.matrix <- makeContrastsFromLabelMatrix(contr.matrix)
  contr.matrix <- sign(contr.matrix) ## sign is fine

  ## sanity check
  if (NCOL(contr.matrix) == 0) {
    message("[pgx.computePGX] WARNING: FATAL ERROR. zero contrasts")
    return(pgx)
  }

  ## select valid contrasts
  sel <- Matrix::colSums(contr.matrix == -1) > 0 & Matrix::colSums(contr.matrix == 1) > 0
  contr.matrix <- contr.matrix[, sel, drop = FALSE]
 
  ## -------------------------------------------------------------------
  ## Clustering
  ## -------------------------------------------------------------------
  ## Cluster by sample
  if (do.cluster || cluster.contrasts) {
    message("[pgx.computePGX] clustering samples...")
    mm <- c("pca", "tsne", "umap")
    pgx <- pgx.clusterSamples(pgx, dims = c(2, 3), perplexity = NULL, X = NULL, methods = mm)
  }

  ## Make contrasts by cluster
  if (cluster.contrasts) {
    ## NEED RETHINK: for the moment we use combination of t-SNE/UMAP
    posx <- cbind(pgx$cluster$pos[["umap2d"]], pgx$cluster$pos[["tsne2d"]])
    posx <- scale(posx)
    idx <- pgx.findLouvainClusters(posx, level = 1, prefix = "c", small.zero = 0.0)
    if (length(unique(idx)) == 1) {
      ## try again with finer settings if single cluster...
      idx <- pgx.findLouvainClusters(posx, level = 2, prefix = "c", small.zero = 0.01)
    }
    pgx$samples$.cluster <- idx ## really add??

    ## Add cluster contrasts
    message("[pgx.computePGX] adding cluster contrasts...")
    Y <- pgx$samples[, ".cluster", drop = FALSE]
    if (length(unique(Y[, 1])) < 2) {
      message("[pgx.computePGX] warning: only one cluster.")
    } else {
      ct <- makeDirectContrasts(Y, ref = "others")
      ctx <- contrastAsLabels(ct$exp.matrix)
      if (ncol(pgx$contrasts) == 0) {
        pgx$contrasts <- ctx
      } else {
        pgx$contrasts <- cbind(pgx$contrasts, ctx)
      }
    }
  }

  ## Cluster by genes
  if (do.clustergenes) {
    message("[pgx.computePGX] clustering genes...")
    mm <- "umap"
    if (pgx$datatype == "scRNAseq") mm <- c("pca", "tsne", "umap")
    pgx <- pgx.clusterGenes(pgx, methods = mm, level = "gene")
  }

  ## -----------------------------------------------------------------------------
  ## Filter genes (previously in compute_testGenesSingleOmics). NEED
  ## RETHINK?? MOVE TO PGXCREATE??
  ## -----------------------------------------------------------------------------

  ## Shrink number of genes (highest SD/var)
  if (max.genes > 0 && nrow(pgx$counts) > max.genes) {
    message("shrinking data matrices: n= ", max.genes)
    logcpm <- logCPM(pgx$counts, total = NULL)
    sdx <- matrixStats::rowSds(logcpm, na.rm = TRUE)
    jj <- Matrix::head(order(-sdx), max.genes) ## how many genes?
    jj0 <- setdiff(seq_len(nrow(pgx$counts)), jj)
    pgx$filtered[["low.variance"]] <- paste(rownames(pgx$counts)[jj0], collapse = ";")
    pgx$counts <- pgx$counts[jj, ]
  }

  gg <- intersect(rownames(pgx$counts), rownames(pgx$X))
  pgx$counts <- pgx$counts[gg, ]
  pgx$X <- pgx$X[gg, ]

  ## ======================================================================
  ## ================= Run tests ==========================================
  ## ======================================================================

  pgx$timings <- c()
  GENETEST.METHODS <- c(
    "ttest", "ttest.welch", "ttest.rank",
    "voom.limma", "trend.limma", "notrend.limma",
    "edger.qlf", "edger.lrt", "deseq2.wald", "deseq2.lrt"
  )
  GENESETTEST.METHODS <- c(
    "fisher", "gsva", "ssgsea", "spearman",
    "camera", "fry", "fgsea"
  ) ## no GSEA, too slow...

  ## ------------------ gene level tests ---------------------
  if (!is.null(progress)) progress$inc(0.1, detail = "testing genes")

  timeseries <- any(grepl("^IA:*", colnames(pgx$contrasts)))
  
  message("[pgx.computePGX] testing genes...")
  pgx <- compute_testGenes(
    pgx,
    contr.matrix,
    max.features = max.genes,
    test.methods = gx.methods,
    custom_fc = custom_fc,
    use.design = use.design,
    prune.samples = prune.samples,
    timeseries = timeseries,
    remove.outputs = TRUE
  )

  ## ------------------ gene set tests -----------------------
  if (!is.null(progress)) progress$inc(0.2, detail = "testing gene sets")

  if ((pgx$organism != "No organism" && !is.null(pgx$GMT) && nrow(pgx$GMT) > 0) ||
      (pgx$organism == "No organism" && !is.null(custom.geneset$gmt))) {
    message("[pgx.computePGX] testing genesets...")

    pgx <- compute_testGenesets(
      pgx = pgx,
      custom.geneset = custom.geneset,
      test.methods = gset.methods,
      use.replaid = TRUE
    )

    ## Cluster by genes
    if (do.clustergenesets) {
      message("[pgx.computePGX] clustering genesets...")
      pgx <- pgx.clusterGenes(pgx, methods = "umap", X = NULL, level = "geneset")
    }
  } else {
    message("[pgx.computePGX] Skipping genesets test")
  }

  ## ------------------ extra analyses ---------------------
  if (!is.null(progress)) progress$inc(0.3, detail = "extra modules")
  message("[pgx.computePGX] computing extra modules: ", paste0(extra.methods, collapse = "; "))
  pgx <- compute_extra(
    pgx,
    extra = extra.methods,
    pgx.dir = pgx.dir,
    libx.dir = libx.dir,
    user_input_dir = user_input_dir
  )

  info("[pgx.computePGX] DONE")
  return(pgx)
}


## ===================================================================
## =================== UTILITY FUNCTIONS =============================
## ===================================================================

counts.removeSampleOutliers <- function(counts) {
  ## remove samples with 1000x more or 1000x less total counts (than median)
  totcounts <- colSums(counts, na.rm = TRUE)
  mx <- median(log10(totcounts))
  ex <- (log10(totcounts) - mx)
  sel <- which(abs(ex) > 3 | totcounts < 1) ## allowed: 0.001x - 1000x
  sel
  if (length(sel)) {
    message("[createPGX] WARNING: bad samples. Removing samples: ", paste(sel, collapse = " "))
    counts <- counts[, -sel, drop = FALSE]
  }
  counts
}


counts.removeXXLvalues <- function(counts, xxl.val = NA, zsd = 10) {
  ## remove extra-large and infinite values
  ## X <- log2(1 + counts)
  X <- logCPM(counts)
  which.xxl <- which(is.xxl(X), arr.ind = TRUE)
  nxxl <- nrow(which.xxl)
  if (nxxl > 0) {
    message("[createPGX] WARNING: setting ", nxxl, " XXL values to NA")
    counts[which.xxl] <- xxl.val
  } else {
    message("[createPGX] no XXL values detected")
  }
  counts
}

counts.imputeMissing <- function(counts, method = "SVD2") {
  epsx <- min(counts[counts > 0], na.rm = TRUE)
  X <- log2(epsx + counts)
  is.mox <- is.multiomics(rownames(X))
  if (is.mox) {
    impX <- imputeMissing.mox(X, method = method)
  } else {
    impX <- imputeMissing(X, method = method)
  }
  pmax(2**impX - epsx, 0)
}

#' @export
counts.autoScaling <- function(counts) {
  message("[createPGX] scaling counts...")
  counts_multiplier <- 1

  ## If the difference in total counts is too large, we need to
  ## euqalize them because the thresholds can become strange. Here
  ## we decide if normalizing is necessary (WARNING changes total
  ## counts!!!)
  totcounts <- Matrix::colSums(counts, na.rm = TRUE)
  totratio <- log10(max(1 + totcounts, na.rm = TRUE) / min(1 + totcounts, na.rm = TRUE))
  totratio

  if (totratio > 6) {
    message("[createPGX:autoscale] WARNING: too large total counts ratio. forcing normalization.")
    meancounts <- exp(mean(log(1 + totcounts), na.rm = TRUE))
    counts <- t(t(counts) / totcounts) * meancounts
  }

  ## Check if too big (more than billion reads). This is important
  ## for some proteomics intensity signals that are in billions of
  ## units.
  mean.counts <- mean(Matrix::colSums(counts, na.rm = TRUE))
  is.toobig <- log10(mean.counts) > 9
  if (is.toobig) {
    ## scale to about 10 million reads
    message("[createPGX:autoscale] WARNING: too large total counts. Scaling down to 10e6 reads.")
    unit <- 10**(round(log10(mean.counts)) - 7)
    unit
    counts <- counts / unit
    counts_multiplier <- unit
  }
  counts_multiplier
  message("[createPGX:autoscale] count_multiplier= ", counts_multiplier)

  list(counts = counts, counts_multiplier = counts_multiplier)
}

#' @export
counts.mergeDuplicateFeatures <- function(counts, is.counts = TRUE) {
  counts <- counts[rownames(counts) != "", ]
  counts[which(is.nan(counts))] <- NA
  ndup <- sum(duplicated(rownames(counts)))
  if (ndup > 0) {
    if (!is.counts) counts <- 2**counts
    message("[mergeDuplicateFeatures] ", ndup, " duplicated rownames: averaging rows (in counts).")
    counts <- playbase::rowmean(counts, group = rownames(counts), reorder = TRUE)
    counts[which(is.nan(counts))] <- NA
    if (!is.counts) counts <- log2(counts)
  }
  counts
}

#' @export
pgx.filterZeroCounts <- function(pgx) {
  ## There is second filter in the statistics computation. This
  ## first filter is primarily to reduce the counts table.
  ## AZ: added na.rm=TRUE to avoid introducing NAs and edit to keep NAs.
  keep <- (Matrix::rowMeans(pgx$counts > 0, na.rm = TRUE) > 0) ## at least in one...

  nas <- which(is.na(keep))
  jj <- which(keep)
  if (is.null(nas)) {
    keep <- names(keep)[jj]
  } else {
    keep <- names(keep)[c(nas, jj)]
  }

  pgx$counts <- pgx$counts[keep, , drop = FALSE]
  pgx$X <- pgx$X[keep, , drop = FALSE]
  pgx$genes <- pgx$genes[keep, , drop = FALSE]

  pgx
}

#' @export
pgx.filterLowExpressed <- function(pgx, prior.cpm = 1) {
  AT.LEAST <- ceiling(pmax(2, 0.01 * ncol(pgx$counts)))
  message("filtering for low-expressed genes: > ", prior.cpm, " CPM in >= ", AT.LEAST, " samples")
  keep <- (rowSums(edgeR::cpm(pgx$counts) > prior.cpm, na.rm = TRUE) >= AT.LEAST)
  pgx$filtered <- NULL
  pgx$filtered[["low.expressed"]] <- paste(rownames(pgx$counts)[which(!keep)], collapse = ";")
  pgx$counts <- pgx$counts[keep, , drop = FALSE]
  message("filtering out ", sum(!keep), " low-expressed genes")
  message("keeping ", sum(keep), " expressed genes")
  if (!is.null(pgx$X)) {
    ## WARNING: counts and X should match dimensions.
    pgx$X <- pgx$X[which(keep), , drop = FALSE]
  }
  pgx
}


#' Internal use: append gmt (list) to a sparse gene set matrix
#' 
.append_gmt_to_matrix <- function(gmt, G, all_genes, minsize, maxsize) {

  if(is.null(all_genes)) {
    all_genes <- unique(unlist(gmt))
    if(!is.null(G)) all_genes <- unique(c(rownames(G),all_genes))
  }

  gmt <- lapply(gmt, function(s) intersect(s, all_genes))
  gmt.size <- sapply(gmt, length)
  if( sum(gmt.size >= minsize & gmt.size <= maxsize) == 0) {
    message("[.append_gmt_to_matrix] warning no valid gmt to add")
    return(G)
  }
  
  add_gmt <- createSparseGenesetMatrix(
    gmt.all = gmt,
    min.geneset.size = minsize,
    max.geneset.size = maxsize,
    min_gene_frequency = 1,
    all_genes = all_genes,
    annot = NULL,
    filter_genes = FALSE
  )

  # G and custom_gmt have to be SYMBOL alligned
  if (!is.null(add_gmt) && ncol(add_gmt) > 0) {
    # only run this code if custom_gmt has columns (genes)
    ## colnames(custom_gmt) <- probe2symbol(
    ##   colnames(custom_gmt), pgx$genes, "symbol",
    ##   fill_na = TRUE
    ## )
    G <- merge_sparse_matrix(G, Matrix::t(add_gmt) )
    remove(add_gmt)
  }
  return(G)
}

pgx.add_GMT <- function(pgx, custom.geneset = NULL, max.genesets = 20000) {

  if (!"symbol" %in% colnames(pgx$genes)) {
    message(paste(
      "[pgx.add_GMT] ERROR: could not find 'symbol' column.",
      "Is this an old gene annotation?"
    ))
    return(pgx)
  }

  ## -----------------------------------------------------------
  ## Load Geneset matrix and filter genes by gene or homologous
  ## -----------------------------------------------------------
  message("[pgx.add_GMT] Creating GMT matrix... ")

  # Load geneset matrix from playdata. add metabolomics if data.type
  # is metabolomics
  target <- c("human_ortholog", "symbol", "gene_name", "rownames")
  ortho.col <- intersect(target, colnames(pgx$genes))
  if (length(ortho.col) == 0) {
    symbol <- toupper(pgx$genes$symbol)
  } else {
    symbol <- pgx$genes[, ortho.col[1]]  ## human symbol!
  }

  ## check if we have genes/proteins
  symbol <- sub(".*:","",symbol) ## strip prefix
  sum.px <- sum(symbol %in% colnames(playdata::GSETxGENE),na.rm=TRUE) 
  has.px <- sum.px >= 10

  ## check if we have metabolites/lipids
  has.mx1 <- grepl("metabolomics|lipidomics",pgx$datatype,ignore.case=TRUE)
  has.mx2 <- pgx$datatype=="multi-omics" && any(grepl("mx|metabolomics|lipidomics",pgx$genes$data_type))
  has.mx3 <- pgx$datatype=="multi-omics" && !all(grepl("mx|metabolomics|lipidomics",pgx$genes$data_type))  

  has.mx  <- has.mx1 || has.mx2
  has.px2 <- !has.mx1 || has.mx3
  
  dbg("[pgx.add_GMT] 1: has.px = ", has.px)
  dbg("[pgx.add_GMT] 1: has.px2 = ", has.px2)  
  dbg("[pgx.add_GMT] 1: has.mx = ", has.mx)
  
  ## Note!!!: Rownames of G must be in species symbol (not anymore
  ## human ortholog).
  G <- NULL
  
  ## add metabolomic gene sets
  if (has.mx) {
    info("[pgx.add_GMT] Retrieving metabolomics genesets")
    G <- mx.create_metabolite_sets(
      annot = pgx$genes,
      gmin = 0,
      metmin = 5,
      as_matrix = TRUE
    )
  }

  ## add SYMBOL (classic) gene sets
  if (has.px) {
    info("[pgx.add_GMT] Retrieving transcriptomics/proteomics genesets")
    G1 <- Matrix::t(playdata::GSETxGENE)
    G1 <- rename_by2(G1, pgx$genes, new_id="symbol")  ## symbol!
    G <- merge_sparse_matrix(G, G1)
  }

  # create a feature list that will be used to filter and reduce dimensions of G
  full_feature_list <- c(
    pgx$genes$human_ortholog, pgx$genes$symbol,
    rownames(pgx$genes)
  )
  full_feature_list <- full_feature_list[!is.na(full_feature_list)]
  full_feature_list <- full_feature_list[full_feature_list != ""]
  full_feature_list <- unique(full_feature_list)

  if (!is.null(G)) {
    G <- G[rownames(G) %in% full_feature_list, , drop = FALSE]
    G <- G[, Matrix::colSums(G != 0) > 0, drop = FALSE]
    if (nrow(G) == 0 || ncol(G) == 0) G <- NULL
  }

  
  ## Add organism specific GO gene sets. This is species gene symbol.
  if (has.px2) {
    ## add species GO genesets from AnnotationHub
    go.genesets <- NULL
    info("[pgx.add_GMT] Adding species GO for organism", pgx$organism)
    go.genesets <- tryCatch(
      {
        getOrganismGO(pgx$organism)
      },
      error = function(e) {
        message("Error in getOrganismsGO:", e)
      }
    )

    if (!is.null(go.genesets)) {
      dbg("[pgx.add_GMT] got", length(go.genesets), "GO genesets")
      all_genes <- unique(pgx$genes$symbol)
      go_genes <- unique(unlist(go.genesets))
      go_genes2 <- paste0("SYMBOL:",unique(unlist(go.genesets)))
      if( sum(go_genes2 %in% all_genes) > sum(go_genes %in% all_genes) ) {
        go.genesets <- lapply(go.genesets, function(m) paste0("SYMBOL:",m))
      }
      G <- .append_gmt_to_matrix(go.genesets, G, all_genes, minsize=15, maxsize=400)
    } ## end-if go.genesets
  } ## end-if !metabolics

  ## Add custom gene sets if provided
  if (!is.null(custom.geneset$gmt)) {
    message("[pgx.add_GMT] Adding custom genesets...")
    ## convert gmt standard to SPARSE matrix: gset in rows, genes in columns.
    custom_gmt <- custom.geneset$gmt
    customG <- lapply(custom_gmt, function(s) probe2symbol(s, pgx$genes, "symbol"))
    all_genes <- unique(pgx$genes$symbol)
    G <- .append_gmt_to_matrix(customG, G, all_genes, minsize=3, maxsize=9999)
  }

  ## -----------------------------------------------------------
  ##  Prioritize gene sets by fast rank-correlation
  ## -----------------------------------------------------------
  ## NEED RETHINK!! IK. Probably not needed anymore with generalized
  ## features. Generally G and X are not aligned anymore.
  ## !!!!!!!!!!!!
  ## NOTE: this can be replace by PLAID??
  
  if (is.null(max.genesets)) max.genesets <- 20000
  if (max.genesets < 0) max.genesets <- 20000
  if (!is.null(G) && ncol(G) > max.genesets) {
    message("[pgx.add_GMT] Matching gene set matrix...")
    # we use SYMBOL as rownames
    gX <- pgx$X
    if (!all(rownames(gX) %in% pgx$genes$symbol)) {
      gX <- rename_by(gX, pgx$genes, "symbol", unique = TRUE)
    }

    ## if reduced samples
    ss <- rownames(pgx$model.parameters$exp.matrix)
    if (!is.null(ss)) {
      gX <- gX[, ss, drop = FALSE]
    }

    ## Align the GENESETxGENE matrix with genes in X_geneset
    gg <- rownames(gX)
    ii <- intersect(gg, rownames(G))
    G <- G[ii, , drop = FALSE]
    ## gX <- gX[ii, , drop = FALSE]
    xx <- setdiff(gg, rownames(G))
    matX <- Matrix::Matrix(0, nrow = length(xx), ncol = ncol(G), sparse = TRUE)
    rownames(matX) <- xx
    colnames(matX) <- colnames(G)
    G <- rbind(G, matX)
    G <- G[match(gg, rownames(G)), , drop = FALSE]
    rownames(G) <- rownames(gX) ## must be symbol

    ## Prioritize gene sets by fast rank-correlation
    message("[pgx.add_GMT] Reducing gene set matrix... ")
    ## Reduce gene sets by selecting top varying genesets. We use the
    ## very fast sparse rank-correlation for approximate single sample
    ## geneset activation.
    cX <- gX - rowMeans(gX, na.rm = TRUE) ## center!
    cX <- t(matrixStats::colRanks(cX))
    if (ncol(cX) <= 5000) {
      gsetX <- qlcMatrix::corSparse(G, cX)
    } else { ## split into chuncks. faster & needs less memory.
      index <- unique(c(seq(1, ncol(cX), by = round(ncol(cX) / 10, 0)), ncol(cX)))
      i <- 1
      LL.cor <- list()
      for (i in 1:(length(index) - 1)) {
        if (index[i] == 1) jj <- 1:index[i + 1]
        if (index[i] > 1) jj <- (index[i] + 1):index[i + 1]
        LL.cor[[i]] <- qlcMatrix::corSparse(G, cX[, jj, drop = FALSE])
      }
      gsetX <- do.call(cbind, LL.cor)
      rm(index, LL.cor)
    }
    grp <- pgx$model.parameters$group
    gsetX.bygroup <- NULL
    ## If groups/conditions are present we calculate the SD by group
    if (!is.null(grp)) {
      gsetX.bygroup <- tapply(1:ncol(gsetX), grp, function(i) rowMeans(gsetX[, i, drop = FALSE], na.rm = TRUE))
      gsetX.bygroup <- do.call(cbind, gsetX.bygroup)
      ## sdx <- apply(gsetX.bygroup, 1, stats::sd, na.rm = TRUE)
      sdx <- matrixStats::rowSds(gsetX.bygroup, na.rm = TRUE)
    } else {
      sdx <- matrixStats::rowSds(gsetX, na.rm = TRUE)
    }
    names(sdx) <- colnames(G)
    jj <- Matrix::head(order(-sdx), max.genesets)
    must.include <- "hallmark|kegg|^go|^celltype|^pathway|^custom|^metabo"
    jj <- unique(c(jj, grep(must.include, colnames(G), ignore.case = TRUE)))
    jj <- jj[order(colnames(G)[jj])] ## sort alphabetically
    G <- G[, jj, drop = FALSE]
    rm(gsetX.bygroup, gsetX)
  }

  ## -----------------------------------------------------------------------
  ## Clean up and return pgx object
  ## -----------------------------------------------------------------------

  # final check: drop genesets in G based on geneset size
  if (!is.null(G)) {
    gmt.size <- Matrix::colSums(G != 0)
    has.metabolites <- sum(grepl("^[0-9]+$|CHEBI|LIPID", rownames(G))) >= 10
    has.metabolites
    ## if (pgx$datatype %in% c("metabolomics","multi-omics")) {
    if (has.metabolites) {
      # metabolomics genesets are MUCH smaller than transcriptomics,
      # metabolomics have usually less features, so we need to reduce
      # the min size
      size.ok <- which(gmt.size >= 3 & gmt.size <= 400)
    } else {
      size.ok <- which(gmt.size >= 10 & gmt.size <= 400)
    }

    # add all custom genesets to size.ok
    idx_custom_gmt <- grep("CUSTOM", colnames(G))
    # make sure we dont miss CUSTOM genesets due to size.ok exclusion
    if (length(idx_custom_gmt) > 0) {
      names(idx_custom_gmt) <- colnames(G)[idx_custom_gmt]
      size.ok <- union(size.ok, idx_custom_gmt)
    }
    G <- G[, size.ok, drop = FALSE]
  }

  # add random genesets if G is too small
  if (is.null(G) || ncol(G) < 30 || nrow(G) < 3) {
    add.gmt <- NULL
    rr <- sample(3:400, 50)
    gg <- pgx$genes$symbol
    random.gmt <- lapply(rr, function(n) head(sample(gg), min(n, length(gg) / 2)))
    names(random.gmt) <- paste0("TEST:random_geneset.", 1:length(random.gmt))
    # Extreme low feature count control, avoids crash
    if (all(lapply(random.gmt, length) |> unlist() < 3)) {
      min.geneset.size <- 1
    } else {
      min.geneset.size <- 3
    }

    G <- .append_gmt_to_matrix(
      random.gmt, G,
      all_genes = unique(pgx$genes$symbol),
      minsize = min.geneset.size,
      maxsize = 400
    )

  }

  # normalize columns (required for some methods downstream)log2foldchange
  G <- normalize_cols(G)

  pgx$GMT <- G
  pgx$custom.geneset <- custom.geneset
  message(glue::glue("[pgx.add_GMT] Final GMT: {nrow(G)} x {ncol(G)}"))
  rm(G)

  gc()
  return(pgx)
}



## ----------------------------------------------------------------------
## -------------------------- end of file -------------------------------
## ----------------------------------------------------------------------
