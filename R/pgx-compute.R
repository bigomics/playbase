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
    impX = NULL,
    is.logx = NULL,
    batch.correct = TRUE,
    auto.scale = TRUE,
    filter.genes = TRUE,
    prune.samples = FALSE,
    only.known = TRUE,
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
                          organism = NULL,
                          custom.geneset = NULL,
                          annot_table = NULL,
                          max.genesets = 5000,
                          name = "Data set",
                          datatype = "unknown",
                          creator = "unknown",
                          description = "No description provided.",
                          X = NULL,
                          impX = NULL,
                          norm_method = "CPM",
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
                          remove.outliers = TRUE
                          ) {

    message("[createPGX] datatype = ", datatype)

    if(!is.null(counts)) {
        message("[createPGX] dim.counts: ", dim(counts)[1], ",", dim(counts)[2])
        message("[createPGX] class.counts: ", class(counts))
        nmissing <- sum(is.na(counts))
        message("[createPGX] counts has ", nmissing, " missing values")
    } else {
        stop("[createPGX] FATAL: counts must be provided")
    }

    if(!is.null(X)) {
        message("[createPGX] class.X: ", class(X))
        message("[createPGX] dim.X: ", dim(X)[1], ", ", dim(X)[2])
        message("[createPGX] Normalization method:" , norm_method)
        nmissing <- sum(is.na(X))
        message("[createPGX] X has ", nmissing, " missing values")
    } else {
        stop("[createPGX] FATAL: X must be provided")
    }

    if(!is.null(impX)) {
        message("[createPGX] dim.impX: impX matrix also provided.")
        message("[createPGX] dim.impX: ", dim(impX)[1], ", ", dim(impX)[2])
    }

    if (!is.null(X) && !all(dim(counts) == dim(X))) {
        stop("[createPGX] dimension of counts and X do not match\n")
    }

    if (!all(rownames(counts) == rownames(X))) {
        stop("rownames of counts and X do not match\n")
    }

    if (is.null(organism)) {
        stop("[createPGX] FATAL: organism must be provided")
    }

    ## -------------------------------------------------------------------
    ## clean up input files
    ## -------------------------------------------------------------------
    samples <- as.data.frame(samples, drop = FALSE)
    counts <- as.matrix(counts)
    if (is.null(contrasts)) contrasts <- samples[, 0]

    ## convert old-style contrast matrix to sample-wise labeled contrasts
    contrasts <- playbase::contrasts.convertToLabelMatrix(contrasts, samples)

    ## prune unused samples
    contrasts[contrasts %in% c("", " ", "NA")] <- NA
    used.samples <- names(which(rowSums(!is.na(contrasts)) > 0))
    if (prune.samples && length(used.samples) < ncol(counts)) {
        counts <- counts[, used.samples, drop = FALSE]
        samples <- samples[used.samples, , drop = FALSE]
        contrasts <- contrasts[used.samples, , drop = FALSE] ## sample-based!!!
    }
    
    ## -------------------------------------------------------------------
    ## Auto-scaling (scale down huge values, often in proteomics)
    ## -------------------------------------------------------------------
    res <- counts.autoScaling(counts)
    counts <- res$counts
    counts_multiplier <- res$counts_multiplier
    remove(res)

    ## -------------------------------------------------------------------
    ## conform all matrices
    ## -------------------------------------------------------------------
    message("[createPGX] conforming matrices...")
    kk <- intersect(colnames(counts), rownames(samples))
    kk <- intersect(kk, colnames(X))
    if(!is.null(impX)) {
        kk <- intersect(kk, colnames(impX))
    }
    counts <- counts[, kk, drop = FALSE]
    X <- X[, kk, drop = FALSE]
    samples <- samples[kk, , drop = FALSE]
    if(!is.null(impX)) {
        impX <- impX[, kk, drop = FALSE]
    }
    samples <- utils::type.convert(samples, as.is = TRUE) ## automatic type conversion
    if (all(kk %in% rownames(contrasts))) {
        contrasts <- contrasts[kk, , drop = FALSE]
    }

    ## -------------------------------------------------------------------
    ## create pgx object
    ## -------------------------------------------------------------------
    message("[createPGX] creating pgx object...")

    ## remove special characters from description (other columns too??)
    description <- gsub("[\"\']", " ", description) ## remove quotes (important!!)
    description <- gsub("[\n]", ". ", description) ## replace newline
    description <- trimws(gsub("[ ]+", " ", description)) ## remove ws

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
        impX = impX,
        norm_method = norm_method,
        total_counts = Matrix::colSums(counts, na.rm = TRUE),
        counts_multiplier = counts_multiplier
    )

    ## -------------------------------------------------------------------
    ## create gene annotation table
    ## -------------------------------------------------------------------
    pgx$genes <- NULL

    message("[createPGX] annotating genes")
    pgx <- pgx.addGeneAnnotation(pgx, organism = organism, annot_table = annot_table)
    if (is.null(pgx$genes)) {
        stop("[createPGX] FATAL: Could not build gene annotation")
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
        ## AZ (16.6.24): it crashes in presence of NAs
        ## NEED RETHINK - add datatype as argument.
        if(datatype != "proteomics") {
            message("[createPGX] filtering out lowly expressed genes (zero counts)...")
            pgx <- pgx.filterLowExpressed(pgx, prior.cpm = 1)
        }
        ## Conform gene table
        ii <- match(rownames(pgx$counts), rownames(pgx$genes))
        pgx$genes <- pgx$genes[ii, , drop = FALSE]
    }

    ## -------------------------------------------------------------------
    ## Filter genes?
    ## -------------------------------------------------------------------
    do.filter <- (only.known | only.proteincoding)
    if (do.filter) {
        if (only.known) {
            has.symbol <- (!is.na(pgx$genes$symbol) & pgx$genes$symbol != "")
            pgx$genes <- pgx$genes[which(has.symbol), ]
        }

        is.proteincoding <- grepl("protein.coding", pgx$genes$gene_biotype)
        if (only.proteincoding && any(is.proteincoding)) {
            pgx$genes <- pgx$genes[which(is.proteincoding), , drop = FALSE]
        }

        keep <- match(rownames(pgx$genes), rownames(pgx$counts))
        pgx$counts <- pgx$counts[keep, , drop = FALSE]
        keep <- match(rownames(pgx$genes), rownames(pgx$X))
        pgx$X <- pgx$X[keep, , drop = FALSE]
        if (!is.null(pgx$impX)) {
            keep <- match(rownames(pgx$genes), rownames(pgx$impX))
            pgx$impX <- pgx$impX[keep, , drop = FALSE]
        }
    }

    ## -------------------------------------------------------------------
    ## collapse probe-IDs to gene symbol and aggregate duplicates
    ## -------------------------------------------------------------------
    if (convert.hugo) {
        ## Average duplicated rows if any
        group <- rownames(pgx$counts)
        pgx$counts <- playbase::rowmean(pgx$counts, group = group, reorder = TRUE)
        pgx$counts <- pgx$counts[rownames(pgx$counts) != "", , drop = FALSE]
        prior <- ifelse(norm_method == "CPM", 1, 1e-3)
        pgx$X <- playbase::rowmean(2 ** pgx$X - prior, group = group, reorder = TRUE)
        pgx$X <- log2(pgx$X + prior)
        pgx$X <- pgx$X[rownames(pgx$X) != "", , drop = FALSE]
        if (!is.null(pgx$impX)) {
            pgx$impX <- playbase::rowmean(2 ** pgx$impX - prior, group = group, reorder = TRUE)
            pgx$impX <- log2(pgx$impX + prior)
            pgx$impX <- pgx$impX[rownames(pgx$impX) != "", , drop = FALSE]
        }
        
        ## Collapse features as a comma-separated elements
        message("[createPGX] collapsing probes by SYMBOL")
        symbol <- pgx$genes[rownames(pgx$counts), "symbol"]
        symbol[is.na(symbol)] <- "" ## avoids warning

        agg_features <- aggregate(feature ~ symbol, data = pgx$genes,
                                  function(x) paste(unique(x), collapse = "; "))

        ## Merge by symbol, replace features by collapsed features
        ## pgx$genes <- pgx$genes[match(rownames(pgx$counts), symbol), ]
        ## rownames(pgx$genes) <- rownames(pgx$counts)
        jj <- match(rownames(pgx$counts), rownames(pgx$genes))
        pgx$genes <- pgx$genes[jj, ]
        jj <- match(pgx$genes$symbol, agg_features$symbol)
        pgx$genes$feature <- agg_features[jj, "feature"]

        ## Rename gene_name with new rownames
        pgx$genes$gene_name <- rownames(pgx$counts)
    }

    ## -------------------------------------------------------------------
    ## Infer cell cycle/gender here (before any batchcorrection)
    ## -------------------------------------------------------------------
    ## AZ: for proteomics, this needs gene names.
    pgx <- playbase::compute_cellcycle_gender(pgx, pgx$counts)

    ## -------------------------------------------------------------------
    ## Add GMT
    ## -------------------------------------------------------------------
    ## If no organism, no custom annotation table and no custom geneset, then create empty GMT
    if (pgx$organism == "No organism" && is.null(annot_table) && is.null(custom.geneset)) {
        pgx$GMT <- Matrix::Matrix(0, nrow = 0, ncol = 0, sparse = TRUE)
    } else {
        pgx <- pgx.add_GMT(pgx = pgx, custom.geneset = custom.geneset, max.genesets = max.genesets)
    }

    message("[createPGX] done! PGX being returned")
    message("[createPGX] dim.pgx$counts: ", dim(pgx$counts)[1], ",", dim(pgx$counts)[2])
    message("[createPGX] dim.pgx$X: ", dim(pgx$X)[1], ",", dim(pgx$X)[2])
    message("[createPGX] pgx$counts has ", sum(is.na(pgx$counts)), " missing values")
    message("[createPGX] pgx$X has ", sum(is.na(pgx$X)), " missing values")

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
#' @param use.design Whether to use model design matrix for testing. Default is TRUE.
#' @param prune.samples Whether to remove samples without valid contrasts. Default is FALSE.
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
                           progress = NULL,
                           user_input_dir = getwd()) {

    message("[pgx.computePGX] Starting pgx.computePGX")
    
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

    ## Cluster by sample
    if (do.cluster || cluster.contrasts) {
        message("[pgx.computePGX] clustering samples...")
        mm <- c("pca", "tsne", "umap")
        pgx <- playbase::pgx.clusterSamples2(pgx, dims = c(2, 3), perplexity = NULL, X = NULL, methods = mm)

        ## NEED RETHINK: for the moment we use combination of t-SNE/UMAP
        posx <- scale(cbind(pgx$cluster$pos[["umap2d"]], pgx$cluster$pos[["tsne2d"]]))
        idx <- playbase::pgx.findLouvainClusters(posx, level = 1, prefix = "c", small.zero = 0.0)
        if (length(unique(idx)) == 1) {
            ## try again with finer settings if single cluster...
            idx <- pgx.findLouvainClusters(posx, level = 2, prefix = "c", small.zero = 0.01)
        }
        pgx$samples$cluster <- idx
    }

    ## Cluster by contrasts
    if (cluster.contrasts) {
        ## Add cluster contrasts
        message("[pgx.computePGX] adding cluster contrasts...")
        Y <- pgx$samples[, "cluster", drop = FALSE]
        if (length(unique(Y[, 1])) < 2) {
            message("[pgx.computePGX] warning: only one cluster.")
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

    ## Cluster by genes
    if (do.clustergenes) {
        message("[pgx.computePGX] clustering genes...")
        pgx <- playbase::pgx.clusterGenes(pgx, methods = "umap", dims = c(2, 3), X = pgx$impX, level = "gene")
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
        ## sdx <- apply(logcpm, 1, stats::sd, na.rm = TRUE)
        jj <- Matrix::head(order(-sdx), max.genes) ## how many genes?
        jj0 <- setdiff(seq_len(nrow(pgx$counts)), jj)
        pgx$filtered[["low.variance"]] <- paste(rownames(pgx$counts)[jj0], collapse = ";")
        pgx$counts <- pgx$counts[jj, ]
    }

    gg <- intersect(rownames(pgx$counts), rownames(pgx$X))
    pgx$counts <- pgx$counts[gg, ]
    pgx$X <- pgx$X[gg, ]
    if( !is.null(pgx$impX)) {
        pgx$impX <- pgx$impX[gg, ]
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
    pgx <- playbase::compute_testGenes(pgx, contr.matrix,
                                       max.features = max.genes,
                                       test.methods = gx.methods,
                                       use.design = use.design,
                                       prune.samples = prune.samples)


    ## ------------------ gene set tests -----------------------
    if (!is.null(progress)) progress$inc(0.2, detail = "testing gene sets")

    if (pgx$organism != "No organism" && !is.null(pgx$GMT) && nrow(pgx$GMT) > 0) {
        message("[pgx.computePGX] testing genesets...")

        pgx <- compute_testGenesets(pgx = pgx,
                                    custom.geneset = custom.geneset,
                                    test.methods = gset.methods)

        ## Cluster by genes
        if (do.clustergenesets) {
            message("[pgx.computePGX] clustering genesets...")
            pgx <- playbase::pgx.clusterGenes(pgx, methods = "umap",
                                              dims = c(2, 3), X = NULL, level = "geneset")
        }
    } else {
        message("[pgx.computePGX] Skipping genesets test")
    }

    ## ------------------ extra analyses ---------------------
    if (!is.null(progress)) progress$inc(0.3, detail = "extra modules")
    message("[pgx.computePGX] computing extra modules...")

    pgx <- playbase::compute_extra(pgx, extra = extra.methods,
                                   pgx.dir = pgx.dir,
                                   libx.dir = libx.dir,
                                   user_input_dir = user_input_dir)

    message("[pgx.computePGX] done!")

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
  X <- playbase::logCPM(counts)
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
    genes <- sapply(rownames(counts), function(s) strsplit(s, split = "[;,\\|]")[[1]][1])
    rownames(counts) <- unname(genes)
    ndup <- sum(duplicated(rownames(counts)))
    if (ndup > 0) {
        if (!is.counts) counts <- 2 ** counts
        message("[mergeDuplicateFeatures] ", ndup, " duplicated rownames: averaging rows (in counts).")
        counts <- playbase::rowmean(counts, group=rownames(counts), reorder=TRUE)
        counts[which(is.nan(counts))] <- NA
        if (!is.counts) counts <- log2(counts)
    }
    counts
}

counts.mergeDuplicateFeatures.OLD <- function(counts, is.counts = TRUE, keep.NA = FALSE) { ## Needs review.
    counts <- counts[rownames(counts) != "", ]
    counts[which(is.nan(counts))] <- NA
    gene0 <- rownames(counts)
    gene1 <- sapply(gene0, function(s) strsplit(s, split = "[;,\\|]")[[1]][1])
    ndup <- sum(duplicated(gene1))
    if (ndup > 0) {
        if (!is.counts) counts <- 2 ** counts - 1
        if (!keep.NA) {
            message("[mergeDuplicateFeatures] ", ndup, " duplicated rownames: summing rows (in counts).")
            counts <- base::rowsum(counts, gene1, na.rm = TRUE)
        } else {
            message("[mergeDuplicateFeatures] ", ndup, " duplicated rownames: summing rows (in counts) and keeping NAs.")
            dups <- unique(gene0[duplicated(gene0)])
            i=1
            for(i in 1:length(dups)) {
                jj <- which(rownames(counts) == dups[i])
                mv <- apply(counts[jj,], 2, function(x) sum(is.na(x)))
                ss <- names(mv[mv==length(jj)])
                new_row <- colMeans(counts[jj, ], na.rm = TRUE) ## was colSums()
                if(length(ss)>0) new_row[ss] <- NA
                counts <- rbind(new_row, counts[-jj, ])
                rownames(counts)[1] <- dups[i]
            }
            rownames(counts) <- unname(sapply(rownames(counts), function(s) strsplit(s, split = "[;,\\|]")[[1]][1]))
        }
        if (!is.counts) counts <- log2(counts + 1)
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
    if(is.null(nas)) {
        keep <- names(keep)[jj]
    } else {
        keep <- names(keep)[c(nas,jj)]
    }
    
    pgx$counts <- pgx$counts[keep, , drop = FALSE]
    pgx$X <- pgx$X[keep, , drop = FALSE]
    if (!is.null(pgx$impX)) {
        pgx$impX <- pgx$impX[keep, , drop = FALSE]
    }
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

#' Get GO gene sets for organism directly from
#' AnnotationHub/OrgDB. Restrict to genes as background.
#'
#' export
getOrganismGO <- function(organism, genes, ah = NULL) {
  if (tolower(organism) == "human") organism <- "Homo sapiens"
  if (tolower(organism) == "mouse") organism <- "Mus musculus"
  if (tolower(organism) == "rat") organism <- "Rattus norvegicus"

  ## Load the annotation resource.
  if (is.null(ah)) ah <- AnnotationHub::AnnotationHub()
  cat("querying AnnotationHub for", organism, "\n")

  ahDb <- AnnotationHub::query(ah, pattern = c(organism, "OrgDb"))

  ## select on exact organism name
  ahDb <- ahDb[which(tolower(ahDb$species) == tolower(organism))]
  if (length(ahDb) == 0) {
    return(list())
  }

  ## select latest/last
  ahDb$species
  k <- length(ahDb)
  cat("selecting database for", ahDb$species[k], "\n")
  orgdb <- ahDb[[k]] ## last one, newest version

  go.gmt <- list()
  AnnotationDbi::keytypes(orgdb)
  if (!"GOALL" %in% AnnotationDbi::keytypes(orgdb)) {
    cat("WARNING:: missing GO annotation in database!\n")
  } else {
    ## create GO annotets
    cat("\nCreating GO annotation from AnnotationHub...\n")
    ont_classes <- c("BP", "CC", "MF")
    k <- "BP"
    for (k in ont_classes) {
      go_id <- AnnotationDbi::mapIds(orgdb,
        keys = k, keytype = "ONTOLOGY",
        column = "GO", multiVals = "list"
      )[[1]]
      go_id <- unique(go_id)
      sets <- AnnotationDbi::mapIds(orgdb,
        keys = go_id, keytype = "GOALL",
        column = "SYMBOL", multiVals = "list"
      )

      # intersect with genes should be done later, as we need to count the gset size before
      # sets <- lapply(sets, function(s) intersect(s, genes))

      ## get GO title
      go <- sapply(GO.db::GOTERM[names(sets)], Term)
      new_names <- paste0("GO_", k, ":", go, " (", sub("GO:", "GO_", names(sets)), ")")
      names(sets) <- new_names

      ## add to list

      go.gmt <- c(go.gmt, sets)
    }
  }
  go.gmt
}


pgx.add_GMT <- function(pgx, custom.geneset = NULL, max.genesets = 20000) {
  if (!"symbol" %in% colnames(pgx$genes)) {
    message("[pgx.add_GMT] ERROR: could not find 'symbol' column. Is this an old gene annotation?")
    return(pgx)
  }

  is.human <- (tolower(pgx$organism) %in% c("human", "homo sapiens"))
  is.mouse <- (tolower(pgx$organism) %in% c("mouse", "mus musculus"))
  is.rat <- (tolower(pgx$organism) %in% c("rat", "rattus norvegicus"))

  ## -----------------------------------------------------------
  ## Load Geneset matrix and filter genes by gene or homologous
  ## -----------------------------------------------------------
  message("[pgx.add_GMT] Creating GMT matrix... ")
  # Load geneset matrix
  G <- Matrix::t(playdata::GSETxGENE)
  if (!is.human && !is.null(pgx$genes$human_ortholog)) {
    human_genes <- ifelse(
      !is.na(pgx$genes$human_ortholog),
      pgx$genes$human_ortholog,
      pgx$genes$symbol
    )
  } else {
    human_genes <- pgx$genes$symbol
  }
  G <- G[rownames(G) %in% human_genes, , drop = FALSE]
  dim(G)

  if (nrow(G) == 0) {
    message("[pgx.add_GMT] WARNING : no overlapping genes. no GMT added.")
    return(pgx)
  }

  # Change HUMAN gene names to species symbols if NOT human and human_ortholog column is NOT NULL
  if (!is.human && !is.null(pgx$genes$human_ortholog)) {
    rownames(G) <- pgx$genes$symbol[match(rownames(G), pgx$genes$human_ortholog)]
  }

  ## -----------------------------------------------------------
  ## Filter gene sets on size
  ## -----------------------------------------------------------

  message("[pgx.add_GMT] Filtering gene sets on size...")
  gmt.size <- Matrix::colSums(G != 0)
  size.ok <- which(gmt.size >= 15 & gmt.size <= 400)
  if (length(size.ok) < 100) {
    ## take at least 100 gene sets, adding random selected
    jj <- head(unique(c(size.ok, sample(1:ncol(G)))), 100)
    G <- G[, jj, drop = FALSE]
  } else {
    G <- G[, size.ok, drop = FALSE]
  }

  ## -----------------------------------------------------------
  ## Add custom gene sets if provided
  ## -----------------------------------------------------------

  ## add species GO genesets from AnnotationHub
  dbg("[pgx.add_GMT] Adding species GO for organism", pgx$organism)
  go.genesets <- NULL

  go.genesets <- tryCatch(
    {
      getOrganismGO(pgx$organism, rownames(G), ah = NULL)
    },
    error = function(e) {
      message("Error in getOrganismsGO:", e)
    }
  )

  if (!is.null(go.genesets)) {
    dbg("[pgx.add_GMT] got", length(go.genesets), "genesets")
    go.size <- sapply(go.genesets, length)
    size.ok <- which(go.size >= 15 & go.size <= 400)
    go.genesets <- go.genesets[size.ok]

    # add CU
    custom.geneset$gmt <- c(custom.geneset$gmt, go.genesets)

    # get the length of go.genesets and add to gmt info
    go.size <- sapply(go.genesets, length)
    custom.geneset$info$GSET_SIZE <- c(custom.geneset$info$GSET_SIZE, go.size)
  }
  if (!is.null(custom.geneset$gmt)) {
    message("[pgx.add_GMT] Adding custom genesets...")
    ## convert gmt standard to SPARSE matrix: gset in rows, genes in
    ## columns.

    custom_gmt <- playbase::createSparseGenesetMatrix(
      gmt.all = custom.geneset$gmt,
      min.geneset.size = 3,
      max.geneset.size = 9999,
      min_gene_frequency = 1,
      all_genes = pgx$all_genes,
      annot = pgx$genes,
      filter_genes = FALSE
    )

    custom_gmt <- custom_gmt[, colnames(custom_gmt) %in% pgx$genes$symbol, drop = FALSE]
    ## merge_sparse_matrix removes duplicated genesets
    G <- playbase::merge_sparse_matrix(
      m1 = G,
      m2 = Matrix::t(custom_gmt)
    )
    G <- G[rownames(G) %in% pgx$genes$symbol, , drop = FALSE]
    remove(custom_gmt)
  }

  ## -----------------------------------------------------------
  ## create the full GENE matrix (always collapsed by gene)
  ## -----------------------------------------------------------

  X <- pgx$X
  if (!all(rownames(X) %in% pgx$genes$symbol)) {
    X <- rename_by(X, pgx$genes, "symbol") ## pgx-functions.R
    X <- X[!rownames(X) == "", , drop = FALSE]
    if (any(duplicated(rownames(X)))) {
      X <- log2(rowsum(2**X, rownames(X)))
    }
  }

  ## if reduced samples
  ss <- rownames(pgx$model.parameters$exp.matrix)
  if (!is.null(ss)) {
    X <- X[, ss, drop = FALSE]
  }

  ## -----------------------------------------------------------
  ## Align the GENESETxGENE matrix with genes in X
  ## -----------------------------------------------------------
  message("[pgx.add_GMT] Matching gene set matrix...")
  gg <- rownames(X)
  ii <- intersect(gg, rownames(G))
  G <- G[ii, , drop = FALSE]
  xx <- setdiff(gg, rownames(G))
  matX <- Matrix::Matrix(0, nrow = length(xx), ncol = ncol(G), sparse = TRUE)
  rownames(matX) <- xx
  colnames(matX) <- colnames(G)
  G <- rbind(G, matX)
  G <- G[match(gg, rownames(G)), , drop = FALSE]
  rownames(G) <- rownames(X) ## original name (e.g. mouse)

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
      gsetX.bygroup <- tapply(1:ncol(gsetX), grp, function(i) rowMeans(gsetX[, i, drop = FALSE], na.rm = TRUE))
      gsetX.bygroup <- do.call(cbind, gsetX.bygroup)
      sdx <- apply(gsetX.bygroup, 1, stats::sd, na.rm = TRUE)
    } else {
      sdx <- matrixStats::rowSds(gsetX, na.rm = TRUE)
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

  G <- playbase::normalize_cols(G)

  pgx$GMT <- G
  pgx$custom.geneset <- custom.geneset
  message(glue::glue("[pgx.add_GMT] Final GMT: {nrow(G)}x{ncol(G)}"))
  rm(gsetX.bygroup, gsetX, G)
  gc()
  return(pgx)
}


## ----------------------------------------------------------------------
## -------------------------- end of file -------------------------------
## ----------------------------------------------------------------------
