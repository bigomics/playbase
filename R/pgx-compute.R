##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

#' Create a pgx object
#'
#' This function creates a pgx object from files, which is the core object in the
#' OmicsPlayground.
#'
#' @param counts.file a playbase::COUNTS file
#' @param samples.file a playbase::SAMPLES file
#' @param contrasts.file a playbase::CONTRASTS file
#'
#' @return list. represents a pgx object
#' @export
#' @examples
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
    pheno <- head(grep("^[.]", colnames(samples), value = TRUE, invert = TRUE), 1)
    pheno <- intersect(pheno, colnames(samples))
    Y <- samples[, pheno, drop = FALSE]
    ## automatically guess contrasts
    ac <- pgx.makeAutoContrasts(Y, mingrp = 3, slen = 20, ref = NA)
    contrasts <- contrastAsLabels(ac$exp.matrix)
  }

  head(contrasts)

  ## reduce sample table? Only phenotypes in contrasts
  vv <- unique(sub("[:].*", "", colnames(contrasts)))
  vv
  samples <- data.frame(samples[, vv, drop = FALSE])

  ## other params
  gx.methods <- strsplit(gxmethods, split = ",")[[1]]
  gset.methods <- strsplit(gsetmethods, split = ",")[[1]]
  extra.methods <- strsplit(extra, split = ",")[[1]]

  gx.methods
  gset.methods
  extra.methods

  ## create initial PGX object
  pgx <- pgx.createPGX(
    counts,
    samples = samples,
    contrasts = contrasts,
    X = NULL, ## genes,
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
    max.genesets = 20000,
    gx.methods = gx.methods,
    gset.methods = gset.methods,
    extra.methods = extra.methods,
    do.cluster = TRUE,
    use.design = TRUE,
    prune.samples = FALSE,
    progress = NULL
  )

  ## save
  names(pgx)
  pgx
}


#' Create a pgx object
#'
#' This function creates a pgx object, which is the core object in the
#' OmicsPlayground.
#'
#' @param counts.file a playbase::COUNTS file
#' @param samples.file a playbase::SAMPLES file
#' @param contrasts.file a playbase::CONTRASTS file
#' @param X dataframe. value
#' @param is.logx boolean. value
#' @param batch.correct boolean. value
#' @param auto.scale boolean. value
#' @param filter.genes boolean. value
#' @param prune.samples boolean. value
#' @param only.known boolean. value
#' @param only.hugo boolean. value
#' @param convert.hugo boolean. value
#' @param do.cluster boolean. value
#' @param cluster.contrasts boolean. value
#' @param do.clustergenes boolean. value
#' @param only.proteincoding boolean. value
#'
#' @return list. represents a pgx object
#' @export
#'
#' @examples
#' # first step is to create pgx
#' pgx <- playbase::pgx.createPGX(
#'   counts = playbase::COUNTS,
#'   samples = playbase::SAMPLES,
#'   contrasts = playbase::CONTRASTS
#' )
#'
#' # once pgx is created, we can compute the modules
#' pgx <- playbase::pgx.computePGX(
#'   pgx = pgx
#' )
#'
#' # if you want a more minimal (and quick) example for testing, use the settings below
#'
#' pgx <- playbase::pgx.createPGX(
#'   counts = playbase::COUNTS,
#'   samples = playbase::SAMPLES,
#'   contrasts = playbase::CONTRASTS[1]
#' )
#'
#' pgx <- playbase::pgx.computePGX(
#'   pgx = pgx,
#'   max.genes = 10000,
#'   max.genesets = 1000,
#'   gx.methods = c("ttest.welch"),
#'   gset.methods = c("fisher")
#' )
pgx.createPGX <- function(counts, samples, contrasts, X = NULL, ## genes,
                          is.logx = NULL, batch.correct = TRUE,
                          auto.scale = TRUE, filter.genes = TRUE, prune.samples = FALSE,
                          only.known = TRUE, only.hugo = TRUE, convert.hugo = TRUE,
                          do.cluster = TRUE, cluster.contrasts = FALSE, do.clustergenes = TRUE,
                          only.proteincoding = TRUE) {
  if (0 && !"group" %in% colnames(samples)) {
    stop("samples information must have 'group' column\n")
    return(NULL)
  }

  if (!is.null(X) && !all(dim(counts) == dim(X))) {
    stop("dimension of counts and X do not match\n")
  }

  ## -------------------------------------------------------------------
  ## clean up input files
  ## -------------------------------------------------------------------
  samples <- data.frame(samples)
  counts <- as.matrix(counts)
  if (is.null(contrasts)) contrasts <- samples[, 0]


  ## contrast matrix
  colnames(contrasts)
  is.numbered <- all(unique(as.vector(contrasts)) %in% c(-1, 0, 1))
  is.numbered <- all(sapply(type.convert(data.frame(contrasts), as.is = TRUE), class) %in% c("numeric", "integer"))
  ct.type <- c("labeled (new style)", "numbered (old style)")[1 + 1 * is.numbered]
  is.numbered
  if (is.numbered) {
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

  ## sanity check...
  if (!all(rownames(contrasts) == rownames(samples)) &&
    !all(rownames(contrasts) == colnames(counts))) {
    stop("[createPGX] FATAL :: matrices do not match")
  }



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
  samples <- type.convert(samples, as.is = TRUE) ## automatic type conversion
  if (!is.null(X)) X <- X[, kk, drop = FALSE]
  if (all(kk %in% rownames(contrasts))) {
    contrasts <- contrasts[kk, , drop = FALSE]
  }

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

      cat("[createPGX:autoscale] WARNING: too large total counts. Scaling down to 10e6 reads.\n")
      unit <- 10**(round(log10(mean.counts)) - 7)
      unit
      counts <- counts / unit
      counts_multiplier <- unit
    }
    counts_multiplier
    cat("[createPGX:autoscale] count_multiplier= ", counts_multiplier, "\n")
  }

  if (0 && auto.scale) {
    ## auto-scale down billions of counts like sometimes for proteomics
    q10 <- quantile(counts[counts > 0.25], probs = 0.10)
    q10
    if (q10 > 100) {
      counts <- counts / q10
      counts_multiplier <- q10
    }
    cat("[createPGX:autoscale] count_multiplier= ", counts_multiplier, "\n")
  }

  ## -------------------------------------------------------------------
  ## convert probe-IDs to gene symbol (do not translate yet to HUGO)
  ## -------------------------------------------------------------------
  message("[createPGX] converting probes to symbol...")
  symbol <- playbase::probe2symbol(rownames(counts), type = NULL) ## auto-convert function
  if (mean(rownames(counts) == symbol, na.rm = TRUE) < 0.5) { ## why??
    jj <- which(!is.na(symbol))
    counts <- as.matrix(counts[jj, ])
    rownames(counts) <- symbol[jj]
    if (!is.null(X)) {
      rownames(X) <- rownames(counts)
    }
  }
  dim(counts)

  ## -------------------------------------------------------------------
  ## create ngs object
  ## -------------------------------------------------------------------
  message("[createPGX] creating pgx object...")


  ngs <- list() ## empty object
  ngs$name <- "data set"
  this.date <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

  ngs$date <- this.date
  ngs$datatype <- "unknown"
  ngs$description <- "data set"

  ngs$samples <- data.frame(samples, check.names = FALSE)
  ngs$counts <- as.matrix(counts)
  ngs$contrasts <- contrasts
  ngs$X <- X ## input normalized log-expression (can be NULL)

  ngs$total_counts <- totcounts
  ngs$counts_multiplier <- counts_multiplier

  ## -------------------------------------------------------------------
  ## collapse multiple row for genes by summing up counts
  ## -------------------------------------------------------------------
  ## take only first gene as rowname, retain others as alias
  gene0 <- rownames(ngs$counts)
  gene1 <- gene0
  gene1 <- sapply(gene0, function(s) strsplit(s, split = "[;,\\|]")[[1]][1])

  if (convert.hugo) {
    message("[createPGX] converting to HUGO symbols...")
    gene1 <- playbase::alias2hugo(gene1) ## convert to latest HUGO
  } else {
    message("[createPGX] skip conversion to HUGO symbols")
  }
  ndup <- sum(duplicated(gene1))
  ndup
  if (ndup > 0) {
    message("[createPGX:autoscale] duplicated rownames detected: summing up rows (counts).")
    x1 <- tapply(1:nrow(ngs$counts), gene1, function(i) {
      Matrix::colSums(ngs$counts[i, , drop = FALSE])
    })
    x1 <- do.call(rbind, x1)
    ngs$counts <- x1
    remove(x1)
  }
  if (ndup > 0 && !is.null(ngs$X)) {
    x1 <- tapply(1:nrow(ngs$X), gene1, function(i) {
      log2(Matrix::colSums(2**ngs$X[i, , drop = FALSE]))
    })
    x1 <- do.call(rbind, x1)
    ngs$X <- x1
    remove(x1)
  }

  ## -------------------------------------------------------------------
  ## create gene annotation if not given (no HUGO conversion)
  ## -------------------------------------------------------------------
  message("[createPGX] annotating genes...")
  ngs$genes <- playbase::ngs.getGeneAnnotation(genes = rownames(ngs$counts))
  rownames(ngs$genes) <- rownames(ngs$counts)
  ngs$genes[is.na(ngs$genes)] <- ""

  ## -------------------------------------------------------------------
  ## Filter out not-expressed
  ## -------------------------------------------------------------------
  if (filter.genes) {
    ## There is second filter in the statistics computation. This
    ## first filter is primarily to reduce the counts table.
    message("[createPGX] filtering out not-expressed genes...")
    keep <- (Matrix::rowMeans(ngs$counts > 0) > 0) ## at least in one...
    ngs$counts <- ngs$counts[keep, ]
    ngs$genes <- ngs$genes[keep, , drop = FALSE]
    if (!is.null(ngs$X)) {
      ngs$X <- ngs$X[keep, ]
    }
  }

  ## -------------------------------------------------------------------
  ## Filter genes?
  ## -------------------------------------------------------------------
  cap.fraction <- mean(grepl("^[A-Z][a-z]+", rownames(ngs$counts)), na.rm = TRUE)
  message("[createPGX: filter genes] rownames.ngs.counts = ", head(rownames(ngs$counts)))
  message("[createPGX: filter genes] cap.frac = ", cap.fraction)

  is.mouse <- (cap.fraction > 0.9)
  org <- ifelse(is.mouse, "mouse", "human")
  message("[createPGX] detected organism: ", org, "")

  do.filter <- (only.hugo | only.known | only.proteincoding)
  if (do.filter && org == "mouse") {
    message("[createPGX] filtering genes...")
    SYMBOL <- unlist(as.list(org.Mm.eg.db::org.Mm.egSYMBOL))
    is.hugo <- is.known <- is.protcoding <- TRUE
    if (only.hugo) is.hugo <- (ngs$genes$gene_name %in% SYMBOL)

    if (only.known) {
      is.known <- !grepl("Rik|^Orf|^Loc", ngs$genes$gene_name) ## ???
    }
    if (only.proteincoding) {
      is.protcoding <- ngs$genes$gene_biotype %in% c("protein_coding")
    }
    keep <- (is.known & is.hugo & is.protcoding)
    table(keep)
    ngs$counts <- ngs$counts[keep, ]
    ngs$genes <- ngs$genes[keep, ]
    if (!is.null(ngs$X)) ngs$X <- ngs$X[keep, ]
  }
  if (do.filter && org == "human") {
    message("[createPGX] filtering genes...")
    SYMBOL <- unlist(as.list(org.Hs.eg.db::org.Hs.egSYMBOL))
    is.hugo <- is.protcoding <- is.known <- TRUE
    if (only.hugo) is.hugo <- (ngs$genes$gene_name %in% SYMBOL)
    if (only.known) {
      is.known <- !grepl("^ORF|^LOC", ngs$genes$gene_name) ## ???
    }
    if (only.proteincoding) {
      is.protcoding <- ngs$genes$gene_biotype %in% c("protein_coding")
    }
    keep <- (is.known & is.hugo & is.protcoding)
    table(keep)
    ngs$counts <- ngs$counts[keep, ]
    ngs$genes <- ngs$genes[keep, ]
    if (!is.null(ngs$X)) ngs$X <- ngs$X[keep, ]
  }

  ## -------------------------------------------------------------------
  ## Infer cell cycle/gender here (before any batchcorrection)
  ## -------------------------------------------------------------------
  ngs <- playbase::compute_cellcycle_gender(ngs)
  Matrix::head(ngs$samples)

  ## -------------------------------------------------------------------
  ## Batch-correction (if requested. WARNING: changes counts )
  ## -------------------------------------------------------------------
  batch.par <- c("batch", "batch2")
  has.batchpar <- any(batch.par %in% colnames(ngs$samples))
  if (batch.correct && has.batchpar) {
    b <- "batch"
    bb <- intersect(colnames(ngs$samples), batch.par)
    bb
    for (b in bb) {
      message("[createPGX] batch correcting for parameter '", b, "'\n")
      batch <- ngs$samples$batch
      zz <- which(ngs$counts == 0, arr.ind = TRUE)
      cX <- log2(1 + ngs$counts)
      bx <- ngs$sample[, b]

      message("[createPGX] batch correcting for counts using LIMMA\n")
      cX <- limma::removeBatchEffect(cX, batch = bx) ## in log-space
      cX <- pmax(2**cX - 1, 0)
      cX[zz] <- 0
      ngs$counts <- pmax(cX, 0) ## batch corrected counts...

      if (!is.null(ngs$X)) {
        message("[createPGX] batch correcting for logX using LIMMA\n")
        ngs$X <- limma::removeBatchEffect(ngs$X, batch = bx) ## in log-space
        ngs$X[zz] <- 0
      }
    }
    remove(cX)
  }


  ## -------------------------------------------------------------------
  ## Pre-calculate t-SNE for and get clusters early so we can use it
  ## for doing differential analysis.
  ## -------------------------------------------------------------------
  if (do.cluster) {
    message("[createPGX] clustering samples...")


    ngs <- playbase::pgx.clusterSamples2(
      ngs,
      dims = c(2, 3),
      perplexity = NULL,
      methods = c("pca", "tsne", "umap")
    )

    ## NEED RETHINK: for the moment we use combination of t-SNE/UMAP

    posx <- scale(cbind(ngs$cluster$pos[["umap2d"]], ngs$cluster$pos[["tsne2d"]]))

    idx <- playbase::pgx.findLouvainClusters(posx, level = 1, prefix = "c", small.zero = 0.0)
    table(idx)
    if (length(unique(idx)) == 1) {
      ## try again with finer settings if single cluster...
      idx <- playbase::pgx.findLouvainClusters(posx, level = 2, prefix = "c", small.zero = 0.01)
    }
    ngs$samples$cluster <- idx
    Matrix::head(ngs$samples)
    table(ngs$samples$cluster)
  }

  if (cluster.contrasts) {
    ## Add cluster contrasts
    message("[createPGX] adding cluster contrasts...")
    Y <- ngs$samples[, "cluster", drop = FALSE]
    if (length(unique(Y)) < 2) {
      message("[createPGX] warning: only one cluster.")
    } else {
      ct <- playbase::makeDirectContrasts(Y, ref = "others")
      ctx <- playbase::contrastAsLabels(ct$exp.matrix)
      if (ncol(ngs$contrasts) == 0) {
        ngs$contrasts <- ctx
      } else {
        ngs$contrasts <- cbind(ngs$contrasts, ctx)
      }
    }
  }

  ## -------------------------------------------------------------------
  ## Add normalized log-expression
  ## -------------------------------------------------------------------
  if (is.null(ngs$X)) {
    message("[createPGX] calculating log-expression matrix X...")

    ngs$X <- playbase::logCPM(ngs$counts, total = 1e6, prior = 1)

    dim(ngs$X)
  } else {
    message("[createPGX] using passed log-expression X...")
  }

  if (!all(dim(ngs$X) == dim(ngs$counts))) {
    stop("[createPGX] dimensions of X and counts do not match\n")
  }

  if (do.clustergenes) {
    message("[createPGX] clustering genes...")
    ngs <- playbase::pgx.clusterGenes(ngs, methods = "umap", dims = c(2, 3), level = "gene")
  }

  return(ngs)
}


#' Main function to populate pgx with results
#'
#' @param pgx value
#' @param max.genes value
#' @param max.genesets value
#' @param gx.methods value
#' @param gset.methods value
#' @param do.cluster value
#' @param use.design value
#' @param prune.samples value
#' @param extra.methods value
#' @param progress value
#'
#' @return
#' @export
#'
#' @examples
#' # first step is to create pgx
#' pgx <- playbase::pgx.createPGX(
#'   counts = playbase::COUNTS,
#'   samples = playbase::SAMPLES,
#'   contrasts = playbase::CONTRASTS
#' )
#'
#' # once pgx is created, we can compute the modules
#' pgx <- playbase::pgx.computePGX(
#'   pgx = pgx
#' )
#'
#' # if you want a more minimal (and quick) example for testing, use the settings below
#'
#' pgx <- playbase::pgx.createPGX(
#'   counts = playbase::COUNTS,
#'   samples = playbase::SAMPLES,
#'   contrasts = playbase::CONTRASTS[1]
#' )
#'
#' pgx <- playbase::pgx.computePGX(
#'   pgx = pgx,
#'   max.genes = 10000,
#'   max.genesets = 1000,
#'   gx.methods = c("ttest.welch"),
#'   gset.methods = c("fisher")
#' )
pgx.computePGX <- function(pgx,
                           max.genes = 19999,
                           max.genesets = 5000,
                           gx.methods = c("ttest.welch", "trend.limma", "edger.qlf"),
                           gset.methods = c("fisher", "gsva", "fgsea"),
                           custom.geneset = c(gmt = NULL, info = NULL),
                           do.cluster = TRUE,
                           use.design = TRUE,
                           prune.samples = FALSE,
                           extra.methods = c("meta.go", "infer", "deconv", "drugs", "wordcloud", "wgcna")[c(1, 2)],
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
    contr.matrix <- playbase::makeContrastsFromLabelMatrix(contr.matrix)
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

  pgx <- playbase::compute_testGenes(
    pgx, contr.matrix,
    max.features = max.genes,
    test.methods = gx.methods,
    use.design = use.design,
    prune.samples = prune.samples
  )
  Matrix::head(pgx$gx.meta$meta[[1]])

  ## ------------------ gene set tests -----------------------
  if (!is.null(progress)) progress$inc(0.2, detail = "testing gene sets")

  message("[pgx.computePGX] testing genesets...")
  pgx <- compute_testGenesets(
    pgx,
    custom.geneset = custom.geneset,
    max.features = max.genesets,
    test.methods = gset.methods
  )
  Matrix::head(pgx$gset.meta$meta[[1]])


  if (do.cluster) {
    message("[pgx.computePGX] clustering genes...")
    pgx <- pgx.clusterGenes(pgx, methods = "umap", dims = c(2, 3), level = "geneset") ## gsetX not ready!!
  }


  ## ------------------ extra analyses ---------------------
  if (!is.null(progress)) progress$inc(0.3, detail = "extra modules")
  message("[pgx.computePGX] computing extra modules...")

  pgx <- compute_extra(pgx, extra = extra.methods, libx.dir = libx.dir)

  message("[pgx.computePGX] done!")
  pgx$timings
  return(pgx)
}
