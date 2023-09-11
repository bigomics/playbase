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
                                extra = "meta.go,deconv,infer,drugs,wordcloud") {
  ## compile sample table
  samples <- data.table::fread(samples.file, header = TRUE)
  samples <- data.frame(samples, check.names = FALSE, row.names = 1)

  ## read counts table (allow dup rownames)
  counts <- data.table::fread(counts.file)
  counts.rownames <- counts[[1]]
  counts <- as.matrix(counts[, -1])
  rownames(counts) <- counts.rownames

  ## undo logarithm if necessary
  if (max(counts) < 100) {
    cat("assuming counts were log2 values. undoing logarithm...")
    counts <- 2**counts
  }

  ## match sample table and counts
  kk <- sort(intersect(colnames(counts), rownames(samples)))
  counts <- counts[, kk, drop = FALSE]
  samples <- samples[kk, ]

  ## parse requested phenotypes
  if (!is.null(contrasts.file) && file.exists(contrasts.file)) {
    cat("reading contrasts file", contrasts.file, "\n")
    contrasts <- data.table::fread(contrasts.file, header = TRUE)
    contrasts <- data.frame(contrasts, check.names = FALSE, row.names = 1)
  } else {
    ## take first (not-dotted) column as phenotype vector
    pheno <- utils::head(grep("^[.]", colnames(samples), value = TRUE, invert = TRUE), 1)
    pheno <- intersect(pheno, colnames(samples))
    Y <- samples[, pheno, drop = FALSE]
    ## automatically guess contrasts
    ac <- pgx.makeAutoContrasts(Y, mingrp = 3, slen = 20, ref = NA)
    contrasts <- contrastAsLabels(ac$exp.matrix)
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
#' @param do.cluster Logical indicating whether to run sample clustering. Default is TRUE.
#' @param cluster.contrasts Logical indicating whether to cluster contrasts. Default is FALSE.
#' @param do.clustergenes Logical indicating whether to cluster genes. Default is TRUE.
#' @param only.proteincoding Logical indicating whether to keep only protein-coding genes. Default is TRUE.
#'
#' @return List. PGX object containing input data and parameters.
#'
#' @export
pgx.createPGX <- function(counts, samples, contrasts, X = NULL, ## genes,
                          is.logx = NULL, batch.correct = TRUE,
                          auto.scale = TRUE, filter.genes = TRUE, prune.samples = FALSE,
                          only.known = TRUE, only.hugo = TRUE, convert.hugo = TRUE,
                          do.cluster = TRUE, cluster.contrasts = FALSE, do.clustergenes = TRUE,
                          only.proteincoding = TRUE) {
  if (!is.null(X) && !all(dim(counts) == dim(X))) {
    stop("dimension of counts and X do not match\n")
  }

  ## -------------------------------------------------------------------
  ## clean up input files
  ## -------------------------------------------------------------------
  samples <- data.frame(samples)
  counts <- as.matrix(counts)
  if (is.null(contrasts)) contrasts <- samples[, 0]

  message("[createPGX] input: dim(counts) = ", paste(dim(counts), collapse = "x"))
  message("[createPGX] input: dim(samples) = ", paste(dim(samples), collapse = "x"))
  message("[createPGX] input: dim(contrasts) = ", paste(dim(contrasts), collapse = "x"))

  ## contrast matrix
  is.numbered <- all(unique(as.vector(contrasts)) %in% c(-1, 0, 1))
  is.numbered <- all(sapply(utils::type.convert(data.frame(contrasts), as.is = TRUE), class) %in% c("numeric", "integer"))
  ct.type <- c("labeled (new style)", "numbered (old style)")[1 + 1 * is.numbered]
  is.numbered
  if (is.numbered && ncol(contrasts) > 0) {
    contrasts <- contrastAsLabels(contrasts)
  }

  ## convert group-wise contrast to sample-wise
  grp.idx <- grep("group|condition", tolower(colnames(samples)))[1]
  if (any(!is.na(grp.idx))) {
    # only run the code below if we identify at least one group

    is.group.contrast <- all(rownames(contrasts) %in% samples[, grp.idx])
    is.group.contrast
    if (is.group.contrast && nrow(contrasts) < nrow(samples)) {
      ## group
      grp <- as.character(samples[, grp.idx])
      contrasts.new <- contrasts[grp, , drop = FALSE]
      rownames(contrasts.new) <- rownames(samples)
      contrasts <- contrasts.new
    }
  }

  # prune unused samples
  contrasts[contrasts == ""] <- NA
  used.samples <- names(which(rowSums(!is.na(contrasts)) > 0))
  if (prune.samples && length(used.samples) < ncol(counts)) {
    counts <- counts[, used.samples, drop = FALSE]
    samples <- samples[used.samples, , drop = FALSE]
    contrasts <- contrasts[used.samples, , drop = FALSE] ## sample-based!!!
  }

  ## -------------------------------------------------------------------
  ## conform
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
  ## check counts
  ## -------------------------------------------------------------------
  message("[createPGX] check logarithm/linear...")
  guess.log <- (min(counts, na.rm = TRUE) < 0 || max(counts, na.rm = TRUE) < 100)
  guess.log <- guess.log && is.null(X) && (is.null(is.logx) || is.logx == TRUE)
  guess.log
  if (is.null(is.logx)) {
    is.logx <- guess.log
  }
  is.logx
  if (is.logx) {
    cat("[createPGX] input assumed log-expression (logarithm)\n")
    cat("[createPGX] ...undo-ing logarithm\n")
    counts <- pmax(2**counts - 1, 0) ## undo logarithm
  } else {
    cat("[createPGX] input assumed counts (not logarithm)\n")
  }

  ## -------------------------------------------------------------------
  ## How to deal with missing values??
  ## -------------------------------------------------------------------
  if (any(is.na(counts)) || any(is.infinite(counts))) {
    message("[createPGX] setting missing values to zero")
    counts[is.na(counts) | is.infinite(counts)] <- 0
  }

  ## -------------------------------------------------------------------
  ## global scaling (no need for CPM yet)
  ## -------------------------------------------------------------------
  message("[createPGX] scaling counts...")
  counts_multiplier <- 1
  totcounts <- Matrix::colSums(counts, na.rm = TRUE)
  totcounts
  if (auto.scale) {
    ## If the difference in total counts is too large, we need to
    ## euqalize them because the thresholds can become
    ## strange. Here we decide if normalizing is necessary (WARNING
    ## changes total counts!!!)
    totratio <- log10(max(totcounts, na.rm = TRUE) / min(totcounts, na.rm = TRUE))
    totratio
    if (totratio > 6) {
      cat("[createPGX:autoscale] WARNING: too large differences in total counts. forcing normalization.")
      meancounts <- exp(mean(log(totcounts)))
      meancounts
      counts <- t(t(counts) / totcounts) * meancounts
    }

    ## check if too big (more than billion reads)
    mean.counts <- mean(Matrix::colSums(counts, na.rm = TRUE))
    mean.counts
    is.toobig <- log10(mean.counts) > 9
    is.toobig
    if (is.toobig) {
      ## scale to about 10 million reads
      #
      cat("[createPGX:autoscale] WARNING: too large total counts. Scaling down to 10e6 reads.\n")
      unit <- 10**(round(log10(mean.counts)) - 7)
      unit
      counts <- counts / unit
      counts_multiplier <- unit
    }
    counts_multiplier
    cat("[createPGX:autoscale] count_multiplier= ", counts_multiplier, "\n")
  }

  ## -------------------------------------------------------------------
  ## convert probe-IDs to gene symbol (do not translate yet to HUGO)
  ## -------------------------------------------------------------------
  message("[createPGX] converting probes to symbol...")
  symbol <- probe2symbol(rownames(counts), type = NULL) ## auto-convert function
  if (mean(rownames(counts) == symbol, na.rm = TRUE) < 0.5) { ## why??
    jj <- which(!is.na(symbol))
    counts <- as.matrix(counts[jj, ])
    rownames(counts) <- symbol[jj]
    if (!is.null(X)) {
      rownames(X) <- rownames(counts)
    }
  }

  ## -------------------------------------------------------------------
  ## create pgx object
  ## -------------------------------------------------------------------
  message("[createPGX] creating pgx object...")

  pgx <- list() ## empty object
  pgx$name <- "data set"
  this.date <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  pgx$date <- this.date
  pgx$datatype <- "unknown"
  pgx$description <- "data set"

  pgx$samples <- data.frame(samples, check.names = FALSE)
  pgx$counts <- as.matrix(counts)
  pgx$contrasts <- contrasts
  pgx$X <- X ## input normalized log-expression (can be NULL)

  pgx$total_counts <- totcounts
  pgx$counts_multiplier <- counts_multiplier

  ## -------------------------------------------------------------------
  ## collapse multiple row for genes by summing up counts
  ## -------------------------------------------------------------------
  ## take only first gene as rowname, retain others as alias
  gene0 <- rownames(pgx$counts)
  gene1 <- gene0
  gene1 <- sapply(gene0, function(s) strsplit(s, split = "[;,\\|]")[[1]][1])

  if (convert.hugo) {
    message("[createPGX] converting to HUGO symbols...")
    gene1 <- alias2hugo(gene1) ## convert to latest HUGO
  } else {
    message("[createPGX] skip conversion to HUGO symbols")
  }
  ndup <- sum(duplicated(gene1))
  ndup
  if (ndup > 0) {
    message("[createPGX:autoscale] duplicated rownames detected: summing up rows (counts).")
    x1 <- tapply(1:nrow(pgx$counts), gene1, function(i) {
      Matrix::colSums(pgx$counts[i, , drop = FALSE])
    })
    if (ncol(pgx$counts) == 1) {
      x1 <- matrix(x1, ncol = 1, dimnames = list(names(x1), colnames(pgx$counts)[1]))
    } else {
      x1 <- do.call(rbind, x1)
    }
    pgx$counts <- x1
    remove(x1)
  }
  if (ndup > 0 && !is.null(pgx$X)) {
    x1 <- tapply(1:nrow(pgx$X), gene1, function(i) {
      log2(Matrix::colSums(2**pgx$X[i, , drop = FALSE]))
    })
    x1 <- do.call(rbind, x1)
    pgx$X <- x1
    remove(x1)
  }

  ## -------------------------------------------------------------------
  ## create gene annotation if not given (no HUGO conversion)
  ## -------------------------------------------------------------------
  message("[createPGX] annotating genes...")
  pgx$genes <- ngs.getGeneAnnotation(genes = rownames(pgx$counts))
  rownames(pgx$genes) <- rownames(pgx$counts)
  pgx$genes[is.na(pgx$genes)] <- ""

  ## -------------------------------------------------------------------
  ## Filter out not-expressed
  ## -------------------------------------------------------------------
  if (filter.genes) {
    ## There is second filter in the statistics computation. This
    ## first filter is primarily to reduce the counts table.
    message("[createPGX] filtering out not-expressed genes...")
    keep <- (Matrix::rowMeans(pgx$counts > 0) > 0) ## at least in one...
    pgx$counts <- pgx$counts[keep, , drop = FALSE]
    pgx$genes <- pgx$genes[keep, , drop = FALSE]
    if (!is.null(pgx$X)) {
      pgx$X <- pgx$X[keep, , drop = FALSE]
    }
  }

  ## -------------------------------------------------------------------
  ## Filter genes?
  ## -------------------------------------------------------------------
  organism <- pgx.getOrganism(pgx$counts)
  message("[createPGX] detected organism: ", organism, "")

  do.filter <- (only.hugo | only.known | only.proteincoding)
  if (do.filter && organism %in% c("mouse", "rat")) {
    message("[createPGX] filtering genes...")
    SYMBOL <- unlist(as.list(org.Mm.eg.db::org.Mm.egSYMBOL))
    is.hugo <- is.known <- is.protcoding <- TRUE
    if (only.hugo) is.hugo <- (pgx$genes$gene_name %in% SYMBOL)
    #
    if (only.known) {
      is.known <- !grepl("Rik|^Orf|^Loc", pgx$genes$gene_name) ## ???
    }
    if (only.proteincoding) {
      is.protcoding <- pgx$genes$gene_biotype %in% c("protein_coding")
    }
    keep <- (is.known & is.hugo & is.protcoding)
    pgx$counts <- pgx$counts[keep, , drop = FALSE]
    pgx$genes <- pgx$genes[keep, , drop = FALSE]
    if (!is.null(pgx$X)) pgx$X <- pgx$X[keep, , drop = FALSE]
  }
  if (do.filter && organism == "human") {
    message("[createPGX] filtering genes...")
    SYMBOL <- unlist(as.list(org.Hs.eg.db::org.Hs.egSYMBOL))
    is.hugo <- is.protcoding <- is.known <- TRUE
    if (only.hugo) is.hugo <- (pgx$genes$gene_name %in% SYMBOL)
    if (only.known) {
      is.known <- !grepl("^ORF|^LOC", pgx$genes$gene_name) ## ???
    }
    if (only.proteincoding) {
      is.protcoding <- pgx$genes$gene_biotype %in% c("protein_coding")
    }
    keep <- (is.known & is.hugo & is.protcoding)
    pgx$counts <- pgx$counts[keep, , drop = FALSE]
    pgx$genes <- pgx$genes[keep, , drop = FALSE]
    if (!is.null(pgx$X)) pgx$X <- pgx$X[keep, , drop = FALSE]
  }

  ## -------------------------------------------------------------------
  ## Check bad samples...
  ## -------------------------------------------------------------------
  min.counts <- 1e-3 * mean(colSums(pgx$counts, na.rm = TRUE))
  sel <- which(colSums(pgx$counts, na.rm = TRUE) < pmax(min.counts, 1))
  if (length(sel)) {
    message("[createPGX] *WARNING* bad samples. Removing samples: ", paste(sel, collapse = " "))
    pgx$counts <- pgx$counts[, -sel, drop = FALSE]
    pgx$samples <- pgx$samples[-sel, , drop = FALSE]
    pgx$contrasts <- pgx$contrasts[-sel, , drop = FALSE]
  }

  ## -------------------------------------------------------------------
  ## Infer cell cycle/gender here (before any batchcorrection)
  ## -------------------------------------------------------------------
  pgx <- compute_cellcycle_gender(pgx)

  ## -------------------------------------------------------------------
  ## Batch-correction (if requested. WARNING: changes counts )
  ## -------------------------------------------------------------------
  batch.par <- c("batch", "batch2")
  has.batchpar <- any(batch.par %in% colnames(pgx$samples))
  if (batch.correct && has.batchpar) {
    b <- "batch"
    bb <- intersect(colnames(pgx$samples), batch.par)
    for (b in bb) {
      message("[createPGX] batch correcting for parameter '", b, "'\n")
      batch <- pgx$samples$batch
      zz <- which(pgx$counts == 0, arr.ind = TRUE)
      cX <- log2(1 + pgx$counts)
      bx <- pgx$sample[, b]

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
    }
    remove(cX)
  }

  ## -------------------------------------------------------------------
  ## Pre-calculate t-SNE for and get clusters early so we can use it
  ## for doing differential analysis.
  ## -------------------------------------------------------------------
  if (do.cluster || cluster.contrasts) {
    message("[createPGX] clustering samples...")
    pgx <- pgx.clusterSamples2(
      pgx,
      dims = c(2, 3),
      perplexity = NULL,
      methods = c("pca", "tsne", "umap")
    )

    ## NEED RETHINK: for the moment we use combination of t-SNE/UMAP
    posx <- scale(cbind(pgx$cluster$pos[["umap2d"]], pgx$cluster$pos[["tsne2d"]]))
    idx <- pgx.findLouvainClusters(posx, level = 1, prefix = "c", small.zero = 0.0)
    if (length(unique(idx)) == 1) {
      ## try again with finer settings if single cluster...
      idx <- pgx.findLouvainClusters(posx, level = 2, prefix = "c", small.zero = 0.01)
    }
    pgx$samples$cluster <- idx
  }

  if (cluster.contrasts) {
    ## Add cluster contrasts
    message("[createPGX] adding cluster contrasts...")
    Y <- pgx$samples[, "cluster", drop = FALSE]
    if (length(unique(Y[, 1])) < 2) {
      message("[createPGX] warning: only one cluster.")
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

  ## -------------------------------------------------------------------
  ## Add normalized log-expression
  ## -------------------------------------------------------------------
  if (is.null(pgx$X)) {
    message("[createPGX] calculating log-expression matrix X...")
    pgx$X <- logCPM(pgx$counts, total = 1e6, prior = 1)
  } else {
    message("[createPGX] using passed log-expression X...")
  }

  if (!all(dim(pgx$X) == dim(pgx$counts))) {
    stop("[createPGX] dimensions of X and counts do not match\n")
  }

  if (do.clustergenes) {
    message("[createPGX] clustering genes...")
    pgx <- pgx.clusterGenes(pgx, methods = "umap", dims = c(2, 3), level = "gene")
  }

  return(pgx)
}


#' @title Compute PGX
#' @description Main function to populate pgx with results. The function computes the analysis on a pgx object
#'
#' @param pgx A pgx object containing the input data
#' @param max.genes Maximum number of genes to test. Default is 19999.
#' @param max.genesets Maximum number of gene sets to test. Default is 5000.
#' @param gx.methods Methods for differential expression analysis at the gene level. Default is c("ttest.welch", "trend.limma", "edger.qlf").
#' @param gset.methods Methods for differential analysis at the gene set level. Default is c("fisher", "gsva", "fgsea").
#' @param custom.geneset Custom gene sets to test, as a named list with gmt and info elements.
#' @param do.cluster Whether to run sample clustering. Default is TRUE.
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
                           gx.methods = c("ttest.welch", "trend.limma", "edger.qlf"),
                           gset.methods = c("fisher", "gsva", "fgsea"),
                           custom.geneset = c(gmt = NULL, info = NULL),
                           do.cluster = TRUE,
                           use.design = TRUE,
                           prune.samples = FALSE,
                           extra.methods = c(
                             "meta.go", "infer", "deconv", "drugs",
                             "connectivity", "wordcloud", "wgcna"
                           )[c(1, 2)],
                           libx.dir = NULL,
                           progress = NULL) {
  ## ======================================================================
  ## ======================================================================
  ## ======================================================================

  if (!"contrasts" %in% names(pgx)) {
    stop("[pgx.computePGX] FATAL:: no contrasts in object")
  }

  ## make proper contrast matrix
  contr.matrix <- pgx$contrasts
  contr.values <- unique(as.vector(contr.matrix))
  is.numcontrast <- all(contr.values %in% c(NA, -1, 0, 1))
  is.numcontrast <- is.numcontrast && (-1 %in% contr.values) && (1 %in% contr.values)
  if (!is.numcontrast) {
    contr.matrix <- makeContrastsFromLabelMatrix(contr.matrix)
    contr.matrix <- sign(contr.matrix) ## sign is fine
  }

  ## select valid contrasts
  sel <- Matrix::colSums(contr.matrix == -1) > 0 & Matrix::colSums(contr.matrix == 1) > 0
  contr.matrix <- contr.matrix[, sel, drop = FALSE]

  ## ======================================================================
  ## ======================================================================
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

  pgx <- compute_testGenes(
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
    pgx,
    custom.geneset = custom.geneset,
    max.features = max.genesets,
    test.methods = gset.methods
  )

  if (do.cluster) {
    message("[pgx.computePGX] clustering genes...")
    pgx <- pgx.clusterGenes(pgx, methods = "umap", dims = c(2, 3), level = "geneset") ## gsetX not ready!!
  }


  ## ------------------ extra analyses ---------------------
  if (!is.null(progress)) progress$inc(0.3, detail = "extra modules")
  message("[pgx.computePGX] computing extra modules...")
  pgx <- compute_extra(pgx, extra = extra.methods, libx.dir = libx.dir)

  message("[pgx.computePGX] done!")

  return(pgx)
}
