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
pgx.createFromFiles <- function(counts.file, samples.file, contrasts.file = NULL,
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
    cat("reading contrasts file", contrasts.file, "\n")
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
    batch.correct = TRUE,
    auto.scale = TRUE,
    filter.genes = TRUE,
    prune.samples = FALSE,
    only.known = TRUE,
    only.hugo = TRUE,
    convert.hugo = TRUE,
    do.cluster = TRUE,
    cluster.contrasts = FALSE,
    do.clustergenes = TRUE,
    only.proteincoding = TRUE
  )


  ## start computing PGX object
  pgx <- pgx.computePGX(
    pgx,
    max.genes = 40000,
    max.genesets = 10000,
    gx.methods = gx.methods,
    gset.methods = gset.methods,
    extra.methods = extra.methods,
    do.cluster = TRUE,
    use.design = TRUE,
    prune.samples = FALSE,
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
#' @param batch.correct Logical indicating whether to perform batch correction. Default is TRUE.
#' @param auto.scale Logical indicating whether to automatically scale/center genes. Default is TRUE.
#' @param filter.genes Logical indicating whether to filter lowly expressed genes. Default is TRUE.
#' @param prune.samples Logical indicating whether to remove samples without contrasts. Default is FALSE.
#' @param only.known Logical indicating whether to keep only known genes. Default is TRUE.
#' @param only.hugo Logical indicating whether to convert symbols to HUGO names. Default is TRUE.
#' @param convert.hugo Logical indicating whether to convert symbols to HUGO names. Default is TRUE.
#' @param only.proteincoding Logical indicating whether to keep only protein-coding genes. Default is TRUE.
#' @param custom.geneset Custom gene sets to test, as a named list with gmt and info elements.
#' @param max.genesets Maximum number of gene sets to test. Default is 5000.
#' 
#' @import data.table
#' @return List. PGX object containing input data and parameters.
#'
#' @export
pgx.createPGX <- function(counts,
                          samples,
                          contrasts,
                          organism = "Human",
                          custom.geneset = NULL,
                          max.genesets = 5000,
                          name = "Data set",
                          datatype = "unknown",
                          creator = "unknown",
                          description = "No description provided.",
                          X = NULL,
                          is.logx = NULL,
                          batch.correct = TRUE,
                          auto.scale = TRUE,
                          filter.genes = TRUE,
                          prune.samples = FALSE,
                          only.known = TRUE,
                          only.hugo = TRUE,
                          convert.hugo = TRUE,
                          only.proteincoding = TRUE,
                          remove.xxl = TRUE,
                          remove.outliers = TRUE,
                          normalize = TRUE,
                          use_biomart = NA) {
  if (!is.null(X) && !all(dim(counts) == dim(X))) {
    stop("[createPGX] dimension of counts and X do not match\n")
  }

  ## -------------------------------------------------------------------
  ## clean up input files
  ## -------------------------------------------------------------------
  samples <- data.frame(samples, drop = FALSE)
  counts <- as.matrix(counts)
  if (is.null(contrasts)) contrasts <- samples[, 0]

  message("[createPGX] input: dim(counts) = ", paste(dim(counts), collapse = "x"))
  message("[createPGX] input: dim(samples) = ", paste(dim(samples), collapse = "x"))
  message("[createPGX] input: dim(contrasts) = ", paste(dim(contrasts), collapse = "x"))


  ## convert old-style contrast matrix to sample-wise labeled contrasts
  contrasts <- playbase::contrasts.convertToLabelMatrix(contrasts, samples)

  # prune unused samples
  contrasts[contrasts %in% c("", " ", "NA")] <- NA
  used.samples <- names(which(rowSums(!is.na(contrasts)) > 0))
  if (prune.samples && length(used.samples) < ncol(counts)) {
    counts <- counts[, used.samples, drop = FALSE]
    samples <- samples[used.samples, , drop = FALSE]
    contrasts <- contrasts[used.samples, , drop = FALSE] ## sample-based!!!
  }

  ## -------------------------------------------------------------------
  ## check counts: linear or logarithm?
  ## -------------------------------------------------------------------
  message("[createPGX] check logarithm/linear...")
  guess.log <- (min(counts, na.rm = TRUE) < 0 || max(counts, na.rm = TRUE) < 100)
  guess.log <- guess.log && (is.null(is.logx) || is.logx == TRUE)
  if (is.null(is.logx)) {
    is.logx <- guess.log
  }
  if (is.logx) {
    message("[createPGX] input assumed logarithm: undo-ing logarithm")
    counts <- pmax(2**counts - 1, 0) ## undo logarithm
  } else {
    message("[createPGX] input assumed counts (not logarithm)")
  }

  ## -------------------------------------------------------------------
  ## How to deal with missing or infinite values??
  ## -------------------------------------------------------------------

  ## remove XXL/Infinite values and set to NA
  if (remove.xxl) {
    zsd <- 10 ## default value
    if (is.numeric(remove.xxl)) zsd <- remove.xxl
    counts <- counts.removeXXLvalues(counts, xxl.val = NA, zsd = zsd)
  }

  ## impute missing values
  if (any(is.na(counts))) {
    nmissing <- sum(is.na(counts))
    message("[createPGX] WARNING: data has ", nmissing, " missing values.")
    impute.method <- "SVD2"
    message("[createPGX] Imputing missing values using ", impute.method)
    counts <- counts.imputeMissing(counts, method = impute.method)
  }

  ## -------------------------------------------------------------------
  ## Check bad samples (in total counts, after imputation)
  ## -------------------------------------------------------------------

  ## remove samples from counts matrix with extreme (1000x more or
  ## 1000x less) total counts (than median).
  if (remove.outliers) {
    message("[createPGX] removing outliers samples ")
    counts <- counts.removeOutliers(counts)
  }

  ## -------------------------------------------------------------------
  ## Auto-scaling (scale down huge values, often in proteomics)
  ## -------------------------------------------------------------------
  res <- counts.autoScaling(counts)
  counts <- res$counts
  counts_multiplier <- res$counts_multiplier
  remove(res)

  ## -------------------------------------------------------------------
  ## conform all matrices (after filtering)
  ## -------------------------------------------------------------------
  message("[createPGX] conforming matrices...")
  kk <- intersect(colnames(counts), rownames(samples))
  counts <- counts[, kk, drop = FALSE]
  samples <- samples[kk, , drop = FALSE]
  samples <- utils::type.convert(samples, as.is = TRUE) ## automatic type conversion
  if (!is.null(X)) X <- X[, kk, drop = FALSE]
  if (all(kk %in% rownames(contrasts))) {
    contrasts <- contrasts[kk, , drop = FALSE]
  }

  message("[createPGX] final: dim(counts) = ", paste(dim(counts), collapse = "x"))
  message("[createPGX] final: dim(samples) = ", paste(dim(samples), collapse = "x"))
  message("[createPGX] final: dim(contrasts) = ", paste(dim(contrasts), collapse = "x"))

  ## -------------------------------------------------------------------
  ## COMPUTE LOG NORMALIZE EXPRESSION (if not given)
  ## -------------------------------------------------------------------
  if (is.null(X)) {
    message("[createPGX] creating log-expression matrix X...")
    X <- log2(1 + counts)
  } else {
    message("[createPGX] using passed log-expression matrix X...")
  }

  if (normalize) {
    message("[createPGX] NORMALIZING log-expression matrix X...")
    X <- playbase::logCPM(pmax(2**X - 1, 0), total = 1e6, prior = 1)
    X <- limma::normalizeQuantiles(X) ## in log space
  } else {
    message("[createPGX] SKIPPING NORMALIZATION!")
  }

  ## -------------------------------------------------------------------
  ## create pgx object
  ## -------------------------------------------------------------------
  message("[createPGX] creating pgx object...")

  pgx <- list(
    name = name,
    organism = organism,
    version = packageVersion("playbase"),
    date = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    creator = creator,
    datatype = datatype,
    description = description,
    samples = data.frame(samples, check.names = FALSE),
    counts = as.matrix(counts),
    contrasts = contrasts,
    X = X,
    total_counts = Matrix::colSums(counts, na.rm = TRUE), # input normalized log-expression (can be NULL)
    counts_multiplier = counts_multiplier
  )

  ## -------------------------------------------------------------------
  ## create gene annotation table
  ## -------------------------------------------------------------------
  pgx$genes <- NULL
  if (is.na(use_biomart)) {
    use_biomart <- !(organism %in% c("Mouse", "Human", "Rat"))
  }
  if (!use_biomart && organism %in% c("Mouse", "Human", "Rat")) {
    message("[createPGX] annotating genes using R libraries")
    probe_type <- detect_probe_DEPRECATED(
      probes = rownames(pgx$counts),
      organism = organism
    )
    probe_type
    pgx$genes <- ngs.getGeneAnnotation_DEPRECATED(
      probes = rownames(pgx$counts),
      probe_type = probe_type,
      organism = organism
    )
  }

  if (!use_biomart && !(organism %in% c("Mouse", "Human", "Rat"))) {
    message("ERROR: organism '", organism, "' not supported using R libraries.")
    stop("ERROR: you must set 'use_biomart=TRUE' for organism: ", organism)
  }

  if (use_biomart) {
    message("[createPGX] annotating genes using BiomaRt")
    pgx <- playbase::pgx.gene_table(pgx, organism = organism)
  }

  if (is.null(pgx$genes)) {
    stop("[createPGX] FATAL: Could not get gene annotation!")
  } else {
    message("[createPGX] Gene annotation table created: ", 
            paste(dim(pgx$genes), collapse = "x"))
  }


  ## -------------------------------------------------------------------
  ## convert probe-IDs to gene symbol and aggregate duplicates
  ## -------------------------------------------------------------------
  
  if (convert.hugo) {
    message("[createPGX] converting probes to symbol...")
    symbol <- pgx$genes[rownames(pgx$counts), "symbol"]
    mapped_symbols <- !is.na(symbol) & symbol != ""
    probes_with_symbol <- pgx$genes[mapped_symbols, "feature"]

    ## Update counts and genes
    pgx$counts <- pgx$counts[probes_with_symbol, , drop = FALSE]
    pgx$genes <- pgx$genes[probes_with_symbol, , drop = FALSE]
    pgx$genes$gene_name <- symbol[mapped_symbols]

    # Sum columns of rows with the same gene symbol
    selected_symbols <- symbol[mapped_symbols]
    rownames(pgx$counts) <- selected_symbols
    if (sum(duplicated(selected_symbols)) > 0) {
      message("[createPGX:autoscale] duplicated rownames detected: summing up rows (counts).")
      pgx$counts <- rowsum(pgx$counts, selected_symbols)
    }
    if (!is.null(pgx$X)) {
      # For X, sum the 2^X values of rows with the same gene symbol
      # And then take log2 again.
      pgx$X <- pgx$X[probes_with_symbol, , drop = FALSE]
      rownames(pgx$X) <- selected_symbols
      pgx$X <- log2(rowsum(2**pgx$X, selected_symbols))
    }

    # Collapse feature as a comma-separated elements
    # if multiple rows match to the same gene, then collapse them

    features_collapsed_by_symbol <- aggregate(feature ~ symbol, data = pgx$genes, function(x) paste(unique(x), collapse = "; "))
    pgx$genes <- pgx$genes[!duplicated(pgx$genes$symbol), , drop = FALSE]

    # merge by symbol (we need to remove feature, as the new feature is collapsed)
    pgx$genes$feature <- NULL

    # merge features_collapsde_by_symbol with pgx$genes by the column symbol
    pgx$genes <- merge(pgx$genes, features_collapsed_by_symbol, by = "symbol")
    rownames(pgx$genes) <- pgx$genes$symbol
    pgx$counts <- pgx$counts[rownames(pgx$genes), , drop = FALSE]
  }


  ## -------------------------------------------------------------------
  ## Filter out not-expressed
  ## -------------------------------------------------------------------
  
  if (filter.genes) {
    ## There is second filter in the statistics computation. This
    ## first filter is primarily to reduce the counts table.
    message("[createPGX] filtering out not-expressed genes (zero counts)...")
    pgx <- pgx.filterZeroCounts(pgx)

    ## prefiltering for low-expressed genes (recommended for edgeR and
    ## DEseq2). Require at least in 2 or 1% of total. Specify the
    ## PRIOR CPM amount to regularize the counts and filter genes
    pgx <- pgx.filterLowExpressed(pgx, prior.cpm = 1)

    # Conform gene table
    keep <- intersect(unique(pgx$genes$gene_name), rownames(pgx$counts))
    pgx$genes <- pgx$genes[keep, , drop = FALSE]
  }


  ## -------------------------------------------------------------------
  ## Filter genes?
  ## -------------------------------------------------------------------

  do.filter <- (only.hugo | only.known | only.proteincoding)
  if (do.filter) {
    pgx$genes <- pgx$genes[!is.na(pgx$genes$symbol) | pgx$genes$symbol == "", ]
    if (only.proteincoding) {
      pgx$genes <- pgx$genes[pgx$genes$gene_biotype %in% c("protein_coding", "protein-coding"), ]
    }
    keep <- intersect(unique(pgx$genes$gene_name), rownames(pgx$counts))
    pgx$counts <- pgx$counts[keep, , drop = FALSE]
    if (!is.null(pgx$X)) {
      pgx$X <- pgx$X[keep, , drop = FALSE]
    }
  }


  ## -------------------------------------------------------------------
  ## Infer cell cycle/gender here (before any batchcorrection)
  ## -------------------------------------------------------------------
  pgx <- playbase::compute_cellcycle_gender(pgx, pgx$counts)


  ## -------------------------------------------------------------------
  ## Batch-correction (if requested. WARNING: changes counts )
  ## -------------------------------------------------------------------
  
  has.batchpar <- any(grepl("^batch|^batch2", colnames(pgx$samples), ignore.case = TRUE))
  if (batch.correct && has.batchpar) {
    b <- "batch"
    bb <- grep("^batch|^batch2", colnames(pgx$samples), ignore.case = TRUE, value = TRUE)
    for (b in bb) {
      message("[createPGX] batch correcting for parameter '", b, "'\n")
      zz <- which(pgx$counts == 0, arr.ind = TRUE)
      cX <- log2(1 + pgx$counts)
      bx <- pgx$sample[, b]
      if (length(unique(bx[!is.na(bx)])) > 1) {
        message("[createPGX] batch correcting for counts using LIMMA\n")
        cX <- limma::removeBatchEffect(cX, batch = bx) ## in log-space
        cX <- pmax(2**cX - 1, 0)
        cX[zz] <- 0
        pgx$counts <- pmax(cX, 0) ## batch corrected counts...

        if (!is.null(pgx$X)) {
          message("[createPGX] batch correcting for logX using LIMMA\n")
          pgx$X <- limma::removeBatchEffect(pgx$X, batch = bx) ## in log-space
          pgx$X[zz] <- 0
        }
      } else {
        message("createPGX] invalid batch paramater")
      }
    }
    remove(cX)
  }


  ## -------------------------------------------------------------------
  ## Add GMT
  ## -------------------------------------------------------------------

  pgx <- pgx.add_GMT(pgx, custom.geneset = custom.geneset, max.genesets = max.genesets)

  ### done
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
#' @param cluster.contrasts Logical indicating whether to cluster contrasts. Default is FALSE.
#' @param do.clustergenes Logical indicating whether to cluster genes. Default is TRUE.
#' @param use.design Whether to use model design matrix for testing. Default is TRUE.
#' @param prune.samples Whether to remove samples without valid contrasts. Default is FALSE.
#' @param extra.methods Additional analysis methods to run. Default is c("meta.go", "infer", "deconv", "drugs", "wordcloud", "wgcna")[c(1, 2)].
#' @param libx.dir Directory containing custom analysis modules.
#' @param progress A progress object for tracking status.
#'
#' @return An updated pgx object containing analysis results.
#'
#' @export
pgx.computePGX <- function(pgx,
                           max.genes = 19999,
                           max.genesets = 5000,
                           gx.methods = c("trend.limma", "edger.qlf", "deseq2.wald"),
                           gset.methods = c("fisher", "gsva", "fgsea"),
                           custom.geneset = list(gmt = NULL, info = NULL),
                           do.cluster = TRUE,
                           cluster.contrasts = TRUE,
                           do.clustergenesets = TRUE,
                           do.clustergenes = TRUE,
                           use.design = TRUE,
                           prune.samples = FALSE,
                           extra.methods = c(
                             "meta.go", "infer", "deconv", "drugs",
                             "connectivity", "wordcloud", "wgcna"
                           )[c(1, 2)],
                           pgx.dir = NULL,
                           libx.dir = NULL,
                           progress = NULL) {
  ## ======================================================================
  ## ======================================================================
  ## ======================================================================

  if (!"contrasts" %in% names(pgx)) {
    stop("[pgx.computePGX] FATAL:: no contrasts in object")
  }
  if (!all(grepl("_vs_", colnames(pgx$contrasts)))) {
    stop("[pgx.computePGX] FATAL:: all contrast names must include _vs_")
  }

  contr.matrix <- playbase::contrasts.convertToLabelMatrix(pgx$contrasts, pgx$samples)
  contr.matrix <- playbase::makeContrastsFromLabelMatrix(contr.matrix)
  contr.matrix <- sign(contr.matrix) ## sign is fine

  ## select valid contrasts
  sel <- Matrix::colSums(contr.matrix == -1) > 0 & Matrix::colSums(contr.matrix == 1) > 0
  contr.matrix <- contr.matrix[, sel, drop = FALSE]


  ## -------------------------------------------------------------------
  ## Clustering
  ## -------------------------------------------------------------------
  
  # Cluster by sample
  if (do.cluster || cluster.contrasts) {
    message("[createPGX] clustering samples...")
    pgx <- playbase::pgx.clusterSamples2(
      pgx,
      dims = c(2, 3),
      perplexity = NULL,
      methods = c("pca", "tsne", "umap")
    )

    ## NEED RETHINK: for the moment we use combination of t-SNE/UMAP
    posx <- scale(cbind(pgx$cluster$pos[["umap2d"]], pgx$cluster$pos[["tsne2d"]]))
    idx <- playbase::pgx.findLouvainClusters(posx, level = 1, prefix = "c", small.zero = 0.0)
    if (length(unique(idx)) == 1) {
      ## try again with finer settings if single cluster...
      idx <- pgx.findLouvainClusters(posx, level = 2, prefix = "c", small.zero = 0.01)
    }
    pgx$samples$cluster <- idx
  }

  # Cluster by contrasts
  if (cluster.contrasts) {
    ## Add cluster contrasts
    message("[createPGX] adding cluster contrasts...")
    Y <- pgx$samples[, "cluster", drop = FALSE]
    if (length(unique(Y[, 1])) < 2) {
      message("[createPGX] warning: only one cluster.")
    } else {
      ct <- playbase::makeDirectContrasts(Y, ref = "others")
      ctx <- playbase::contrastAsLabels(ct$exp.matrix)
      if (ncol(pgx$contrasts) == 0) {
        pgx$contrasts <- ctx
      } else {
        pgx$contrasts <- cbind(pgx$contrasts, ctx)
      }
    }
  }

  # Cluster by genes
  if (do.clustergenes) {
    message("[createPGX] clustering genes...")
    pgx <- playbase::pgx.clusterGenes(pgx, methods = "umap", dims = c(2, 3), level = "gene")
  }


  ## -----------------------------------------------------------------------------
  ## Filter genes (previously in compute_testGenesSingleOmics). NEED
  ## RETHINK?? MOVE TO PGXCREATE??
  ## -----------------------------------------------------------------------------
  
  ## Shrink number of genes (highest SD/var)
  if (max.genes > 0 && nrow(pgx$counts) > max.genes) {
    cat("shrinking data matrices: n=", max.genes, "\n")
    logcpm <- logCPM(pgx$counts, total = NULL)
    sdx <- apply(logcpm, 1, stats::sd)
    jj <- Matrix::head(order(-sdx), max.genes) ## how many genes?
    jj0 <- setdiff(seq_len(nrow(pgx$counts)), jj)
    pgx$filtered[["low.variance"]] <- paste(rownames(pgx$counts)[jj0], collapse = ";")
    pgx$counts <- pgx$counts[jj, ]
  }

  if (!is.null(pgx$X)) {
    gg <- intersect(rownames(pgx$counts), rownames(pgx$X))
    pgx$X <- pgx$X[gg, ]
  }

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
  message("[pgx.computePGX] testing genes...")

  pgx <- playbase::compute_testGenes(
    pgx, contr.matrix,
    max.features = max.genes,
    test.methods = gx.methods,
    use.design = use.design,
    prune.samples = prune.samples
  )

  ## ------------------ gene set tests -----------------------
  if (!is.null(progress)) progress$inc(0.2, detail = "testing gene sets")

  message("[pgx.computePGX] testing genesets...")
  pgx <- compute_testGenesets(
    pgx = pgx,
    custom.geneset = custom.geneset,
    test.methods = gset.methods
  )

  # Cluster by genes
  if (do.clustergenes) {
    message("[createPGX] clustering genesets...")
    pgx <- playbase::pgx.clusterGenes(pgx, methods = "umap", dims = c(2, 3), level = "geneset")
  }

  ## ------------------ extra analyses ---------------------
  if (!is.null(progress)) progress$inc(0.3, detail = "extra modules")
  message("[pgx.computePGX] computing extra modules...")
  pgx <- compute_extra(
    pgx,
    extra = extra.methods,
    pgx.dir = pgx.dir,
    libx.dir = libx.dir
  )

  message("[pgx.computePGX] done!")

  return(pgx)
}



## ===================================================================
## =================== UTILITY FUNCTIONS =============================
## ===================================================================

counts.removeOutliers <- function(counts) {
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
  X <- playbase::logCPM(counts)
  ## sdx <- apply(X, 1, function(x) mad(x[x > 0], na.rm = TRUE))
  sdx <- matrixStats::rowSds(X, na.rm = TRUE)
  sdx[is.na(sdx)] <- 0
  sdx0 <- 0.8 * sdx + 0.2 * mean(sdx, na.rm = TRUE) ## moderated SD
  mx <- rowMeans(X, na.rm = TRUE)
  z <- (X - mx) / sdx0
  ## table(abs(z)>10)
  which.xxl <- which(abs(z) > zsd, arr.ind = TRUE)
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
  X <- log2(1 + counts)
  table(is.na(X))
  impX <- playbase::imputeMissing(X, method = method)
  pmax(2**impX - 1, 0)
}

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
    meancounts <- exp(mean(log(1 + totcounts)))
    counts <- t(t(counts) / totcounts) * meancounts
  }

  ## Check if too big (more than billion reads). This is important
  ## for some proteomics intensity signals that are in billions of
  ## units.
  mean.counts <- mean(Matrix::colSums(counts, na.rm = TRUE))
  is.toobig <- log10(mean.counts) > 9
  if (is.toobig) {

    ## scale to about 10 million reads
    message("[createPGX:autoscale] WARNING: too large total counts. Scaling down to 10e6 reads.\n")
    unit <- 10**(round(log10(mean.counts)) - 7)
    unit
    counts <- counts / unit
    counts_multiplier <- unit
  }
  counts_multiplier
  message("[createPGX:autoscale] count_multiplier= ", counts_multiplier, "\n")

  list(counts = counts, counts_multiplier = counts_multiplier)
}

normalizeCounts <- function(M, method = c("TMM", "TMMwsp", "RLE", "upperquartile", "none")) {
  method <- method[1]
  dge <- edgeR::DGEList(M)
  dge <- edgeR::calcNormFactors(dge, method = method)
  logCPM <- edgeR::cpm(dge, log = TRUE)
  logCPM
}


counts.mergeDuplicateFeatures <- function(counts) {
  ## take only first gene as rowname, retain others as alias
  gene0 <- rownames(counts)
  gene1 <- sapply(gene0, function(s) strsplit(s, split = "[;,\\|]")[[1]][1])
  ndup <- sum(duplicated(gene1))
  ndup
  if (ndup > 0) {
    message("[mergeDuplicateFeatures] ", ndup, " duplicated rownames: summing rows (in counts).")
    counts <- base::rowsum(counts, gene1, na.rm = TRUE)
  }
  counts
}

pgx.filterZeroCounts <- function(pgx) {
  ## There is second filter in the statistics computation. This
  ## first filter is primarily to reduce the counts table.
  message("[createPGX] filtering out not-expressed genes (zero counts)...")
  keep <- (Matrix::rowMeans(pgx$counts > 0) > 0) ## at least in one...
  pgx$counts <- pgx$counts[keep, , drop = FALSE]
  if (!is.null(pgx$X)) {
    pgx$X <- pgx$X[keep, , drop = FALSE]
  }
  pgx
}

pgx.filterLowExpressed <- function(pgx, prior.cpm = 1) {
  AT.LEAST <- ceiling(pmax(2, 0.01 * ncol(pgx$counts)))
  cat("filtering for low-expressed genes: >", prior.cpm, "CPM in >=", AT.LEAST, "samples\n")
  keep <- (rowSums(edgeR::cpm(pgx$counts) > prior.cpm, na.rm = TRUE) >= AT.LEAST)
  pgx$filtered <- NULL
  pgx$filtered[["low.expressed"]] <- paste(rownames(pgx$counts)[which(!keep)], collapse = ";")
  pgx$counts <- pgx$counts[keep, , drop = FALSE]
  cat("filtering out", sum(!keep), "low-expressed genes\n")
  cat("keeping", sum(keep), "expressed genes\n")
  if (!is.null(pgx$X)) {
    ## WARNING: counts and X should match dimensions.
    pgx$X <- pgx$X[keep, , drop = FALSE]
  }
  pgx
}

pgx.add_GMT <- function(pgx, custom.geneset = NULL, max.genesets = 20000) {

  if (is.null(pgx$genes$human_ortholog)) {
    # this is needed in case the species is human, and we dont have the homolog column or if we have an old pgx
    # which will ensure consistency between old and new pgx
    pgx$genes$human_ortholog <- NA
  }

  # Load geneset matrix
  G <- playdata::GSETxGENE


  ## -----------------------------------------------------------
  ## Filter genes by gene or homologous, if it exists
  ## -----------------------------------------------------------

  if (pgx$organism != "Human") {
    human_genes <- ifelse(!is.na(pgx$genes$human_ortholog),
      pgx$genes$human_ortholog,
      pgx$genes$symbol
    )
  } else {
    human_genes <- pgx$genes$symbol
  }
  genes <- pgx$genes$symbol

  # Change HUMAN gene names to species symbols if NOT human and human_ortholog column is NOT all NA
  G <- G[, colnames(G) %in% human_genes]

  if (pgx$organism != "Human" && !all(is.na(pgx$genes$human_ortholog))) {
    colnames(G) <- pgx$genes$symbol[match(colnames(G), pgx$genes$human_ortholog)]
  }

  # Normalize G after removal of genes
  G <- playbase::normalize_matrix_by_row(G)


  ## -----------------------------------------------------------
  ## Filter gene sets on size
  ## -----------------------------------------------------------

  message("[pgx.add_GMT] Filtering gene sets on size...")
  gmt.size <- Matrix::rowSums(G != 0)
  size.ok <- (gmt.size >= 15 & gmt.size <= 400)

  # If dataset is too small that size.ok == 0, then select top 100
  if (sum(size.ok) == 0) {
    top_100gs <- utils::head(sort(gmt.size, decreasing = TRUE), 100)
    size.ok <- names(gmt.size) %in% names(top_100gs)
  }

  G <- G[which(size.ok), ]
  G <- Matrix::t(G)


  ## -----------------------------------------------------------
  ## Add custom gene sets if provided
  ## -----------------------------------------------------------

  if (!is.null(custom.geneset$gmt)) {
    message("[pgx.add_GMT] Adding custom genesets...")
    # convert gmt standard to SPARSE matrix
    custom_gmt <- playbase::createSparseGenesetMatrix(
      gmt.all = custom.geneset$gmt,
      min.geneset.size = 3,
      max.geneset.size = 9999,
      min_gene_frequency = 1,
      all_genes = pgx$all_genes,
      annot = pgx$genes,
      filter_genes = FALSE
    )
    custom_gmt <- custom_gmt[, colnames(custom_gmt) %in% genes, drop = FALSE]
    custom_gmt <- playbase::normalize_matrix_by_row(custom_gmt)
    G <- playbase::merge_sparse_matrix(m1 = G, m2 = Matrix::t(custom_gmt))
    remove(custom_gmt)
  }


  ## -----------------------------------------------------------
  ## create the GENESETxGENE matrix
  ## -----------------------------------------------------------
    
  message("[pgx.add_GMT] Matching gene set matrix...")

  gg <- pgx$genes$symbol
  ii <- intersect(gg, rownames(G))
  G <- G[ii, , drop = FALSE]
  xx <- setdiff(gg, rownames(G))
  matX <- Matrix::Matrix(0, nrow = length(xx), ncol = ncol(G), sparse = TRUE)
  rownames(matX) <- xx
  colnames(matX) <- colnames(G)
  G <- rbind(G, matX)
  G <- G[match(gg, rownames(G)), , drop = FALSE]


  ## -----------------------------------------------------------
  ## create the full GENE matrix (always collapsed by gene)
  ## -----------------------------------------------------------

  X <- pgx$X
  single.omics <- TRUE ## !!! for now...
  if (single.omics) {
    ## normalized matrix
    X <- X
  } else {
    data.type <- gsub("\\[|\\].*", "", rownames(pgx$counts))
    jj <- which(data.type %in% c("gx", "mrna"))
    if (length(jj) == 0) {
      stop("[pgx.add_GMT] FATAL. could not find gx/mrna values.")
    }
    X <- X[jj, ]
  }

  if (!all(rownames(X) %in% pgx$genes$symbol)) {
    X <- rename_by(X, pgx$genes, "symbol")
    X <- X[!rownames(X) == "", , drop = FALSE]
    if (any(duplicated(rownames(X)))) {
      X <- log2(rowsum(2**X, rownames(X)))
    }
  }


  ## -----------------------------------------------------------
  ## Prioritize gene sets by fast rank-correlation
  ## -----------------------------------------------------------
  
  if (is.null(max.genesets)) max.genesets <- 20000
  if (max.genesets < 0) max.genesets <- 20000

  if (max.genesets > 0) {
    message("[pgx.add_GMT] Reducing gene set matrix... ")
    ## Reduce gene sets by selecting top varying genesets. We use the
    ## very fast sparse rank-correlation for approximate single sample
    ## geneset activation.
    cX <- X - rowMeans(X, na.rm = TRUE) ## center!
    cX <- apply(cX, 2, rank)
    gsetX <- qlcMatrix::corSparse(G, cX) ## slow!
    grp <- pgx$model.parameters$group
    gsetX.bygroup <- NULL
    ## If groups/conditions are present we calculate the SD by group
    if (!is.null(grp)) {
      gsetX.bygroup <- tapply(1:ncol(gsetX), grp, function(i) rowMeans(gsetX[, i, drop = FALSE]))
      gsetX.bygroup <- do.call(cbind, gsetX.bygroup)
      sdx <- apply(gsetX.bygroup, 1, stats::sd)
    } else {
      sdx <- matrixStats::rowSds(gsetX)
    }

    names(sdx) <- colnames(G)
    jj <- Matrix::head(order(-sdx), max.genesets)
    must.include <- "hallmark|kegg|^go|^celltype|^pathway|^custom"
    jj <- unique(c(jj, grep(must.include, colnames(G), ignore.case = TRUE)))
    jj <- jj[order(colnames(G)[jj])] ## sort alphabetically
    G <- G[, jj, drop = FALSE]
  }


  ## -----------------------------------------------------------------------
  ## Clean up and return pgx object
  ## -----------------------------------------------------------------------

  pgx$GMT <- G
  message(glue::glue("[pgx.add_GMT] Final GMT: {nrow(G)}x{ncol(G)}"))
  rm(gsetX.bygroup, gsetX, G)
  gc()
  return(pgx)
}


## ----------------------------------------------------------------------
## -------------------------- end of file -------------------------------
## ----------------------------------------------------------------------
