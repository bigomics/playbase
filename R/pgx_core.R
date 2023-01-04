#' Create a pgx object
#'
#' This function creates a pgx object (the core object
#' used in playbase) from an input of counts, samples, and
#' contrasts.
#'
#' The pgx object is a list containing the following elements:
#'  - name
#'  - date
#'  - datatype
#'  - description
#'  - samples
#'  - counts
#'  - contrasts
#'  - X
#'  - total_counts
#'  - counts_multiplier
#'  - genes
#'  - and others that are dynamic / undocumented:
#'    - tsne2d
#'    - tsne3d
#'    - cluster
#'    - cluster.genes
#'
#' This function performs the following actions in order on the data:
#'  - process inputs (counts/samples/contrasts)
#'  -
#'
#' Notes:
#'  - consider making "pgx_options" function/object to contain
#'    all of the possible options / parameters that can be passed
#'    to the "create_pgx" and "compute_pgx" functions
#'
#' @param counts dataframe. The counts of the data
#' @param samples dataframe. The samples of the data
#' @param contrasts dataframe. The contrasts of the data
#'
#' @return a list (pgx object)
#' @export
#'
#' @examples
#' counts <- read_counts(example_file("counts.csv"))
#' samples <- read_samples(example_file("samples.csv"))
#' contrasts <- read_contrasts(example_file("contrasts.csv"))
#' ## takes too long
#' # my_pgx <- create_pgx(counts, samples, contrasts)
create_pgx <- function(counts,
                       samples,
                       contrasts,
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

  ## contrast matrix
  is.numbered <- all(sapply(type.convert(data.frame(contrasts)), class) %in% c("numeric", "integer"))
  ct.type <- c("labeled (new style)", "numbered (old style)")[1 + 1 * is.numbered]
  if (is.numbered) {
    contrasts <- contrast_as_labels(contrasts)
  }

  ## convert group-wise contrast to sample-wise
  is.group1 <- FALSE
  is.group2 <- FALSE
  if ("group" %in% colnames(samples)) is.group1 <- all(rownames(contrasts) %in% samples$group)
  if ("condition" %in% colnames(samples)) is.group2 <- all(rownames(contrasts) %in% samples$condition)
  is.group.contrast <- (is.group1 || is.group2)
  if (is.group.contrast && nrow(contrasts) < nrow(samples)) {
    ## group
    if (is.group1) grp <- as.character(samples$group)
    if (is.group2) grp <- as.character(samples$condition)
    contrasts.new <- contrasts[grp, , drop = FALSE]
    rownames(contrasts.new) <- rownames(samples)
    contrasts <- contrasts.new
  }

  ## sanity check...
  if (!all(rownames(contrasts) == rownames(samples)) &&
    !all(rownames(contrasts) == colnames(counts))) {
    stop("matrices do not match")
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
  kk <- intersect(colnames(counts), rownames(samples))
  counts <- counts[, kk, drop = FALSE]
  samples <- samples[kk, , drop = FALSE]
  samples <- type.convert(samples) ## automatic type conversion
  if (!is.null(X)) X <- X[, kk, drop = FALSE]
  if (all(kk %in% rownames(contrasts))) {
    contrasts <- contrasts[kk, , drop = FALSE]
  }

  ## -------------------------------------------------------------------
  ## check counts
  ## -------------------------------------------------------------------
  guess.log <- (min(counts, na.rm = TRUE) < 0 || max(counts, na.rm = TRUE) < 100)
  guess.log <- guess.log && is.null(X) && (is.null(is.logx) || is.logx == TRUE)
  if (is.null(is.logx)) is.logx <- guess.log
  if (is.logx) {
    counts <- pmax(2**counts - 1, 0) ## undo logarithm
  }

  ## -------------------------------------------------------------------
  ## How to deal with missing values??
  ## -------------------------------------------------------------------
  if (any(is.na(counts)) || any(is.infinite(counts))) {
    counts[is.na(counts) | is.infinite(counts)] <- 0
  }

  ## -------------------------------------------------------------------
  ## global scaling (no need for CPM yet)
  ## -------------------------------------------------------------------
  counts_multiplier <- 1
  totcounts <- Matrix::colSums(counts, na.rm = TRUE)
  if (auto.scale) {
    ## If the difference in total counts is too large, we need to
    ## euqalize them because the thresholds can become
    ## strange. Here we decide if normalizing is necessary (WARNING
    ## changes total counts!!!)
    totratio <- log10(max(totcounts, na.rm = TRUE) / min(totcounts, na.rm = TRUE))
    if (totratio > 6) {
      meancounts <- exp(mean(log(totcounts)))
      counts <- t(t(counts) / totcounts) * meancounts
    }

    ## check if too big (more than billion reads)
    mean.counts <- mean(Matrix::colSums(counts, na.rm = TRUE))
    is.toobig <- log10(mean.counts) > 9
    if (is.toobig) {
      ## scale to about 10 million reads
      unit <- 10**(round(log10(mean.counts)) - 7)
      counts <- counts / unit
      counts_multiplier <- unit
    }
  }

  ## -------------------------------------------------------------------
  ## convert probe-IDs to gene symbol (do not translate yet to HUGO)
  ## -------------------------------------------------------------------
  symbol <- probe_to_symbol(rownames(counts), type = NULL) ## auto-convert function
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
    gene1 <- alias_to_hugo(gene1) ## convert to latest HUGO
  }
  ndup <- sum(duplicated(gene1))

  if (ndup > 0) {
    x1 <- tapply(1:nrow(pgx$counts), gene1, function(i) {
      Matrix::colSums(pgx$counts[i, , drop = FALSE])
    })
    x1 <- do.call(rbind, x1)
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
  pgx$genes <- get_gene_annotation(genes = rownames(pgx$counts))
  rownames(pgx$genes) <- rownames(pgx$counts)
  pgx$genes[is.na(pgx$genes)] <- ""

  ## -------------------------------------------------------------------
  ## Filter out not-expressed
  ## -------------------------------------------------------------------
  if (filter.genes) {
    ## There is second filter in the statistics computation. This
    ## first filter is primarily to reduce the counts table.
    keep <- (Matrix::rowMeans(pgx$counts > 0) > 0) ## at least in one...
    pgx$counts <- pgx$counts[keep, ]
    pgx$genes <- pgx$genes[keep, , drop = FALSE]
    if (!is.null(pgx$X)) {
      pgx$X <- pgx$X[keep, ]
    }
  }

  ## -------------------------------------------------------------------
  ## Filter genes?
  ## -------------------------------------------------------------------
  is.mouse <- (mean(grepl("[a-z]", rownames(pgx$counts))) > 0.9)
  org <- ifelse(is.mouse, "mouse", "human")
  do.filter <- (only.hugo | only.known | only.proteincoding)
  if (do.filter && org == "mouse") {
    SYMBOL <- unlist(as.list(org.Mm.eg.db::org.Mm.egSYMBOL))
    is.hugo <- is.known <- is.protcoding <- TRUE
    if (only.hugo) is.hugo <- (pgx$genes$gene_name %in% SYMBOL)
    if (only.known) {
      is.known <- !grepl("Rik|^Orf|^Loc", pgx$genes$gene_name) ## ???
    }
    if (only.proteincoding) {
      is.protcoding <- pgx$genes$gene_biotype %in% c("protein_coding")
    }
    keep <- (is.known & is.hugo & is.protcoding)
    pgx$counts <- pgx$counts[keep, ]
    pgx$genes <- pgx$genes[keep, ]
    if (!is.null(pgx$X)) pgx$X <- pgx$X[keep, ]
  }
  if (do.filter && org == "human") {
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
    pgx$counts <- pgx$counts[keep, ]
    pgx$genes <- pgx$genes[keep, ]
    if (!is.null(pgx$X)) pgx$X <- pgx$X[keep, ]
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
      batch <- pgx$samples$batch
      zz <- which(pgx$counts == 0, arr.ind = TRUE)
      cX <- log2(1 + pgx$counts)
      bx <- pgx$sample[, b]

      cX <- limma::removeBatchEffect(cX, batch = bx) ## in log-space
      cX <- pmax(2**cX - 1, 0)
      cX[zz] <- 0
      pgx$counts <- cX ## batch corrected counts...

      if (!is.null(pgx$X)) {
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
  if (do.cluster) {
    pgx <- cluster_samples2(
      pgx,
      dims = c(2, 3),
      perplexity = NULL,
      methods = c("pca", "tsne", "umap")
    )

    ## NEED RETHINK: for the moment we use combination of t-SNE/UMAP
    posx <- scale(cbind(pgx$cluster$pos[["umap2d"]], pgx$cluster$pos[["tsne2d"]]))
    idx <- find_louvain_clusters(posx, level = 1, prefix = "c", small.zero = 0.0)
    if (length(unique(idx)) == 1) {
      ## try again with finer settings if single cluster...
      idx <- find_louvain_clusters(posx, level = 2, prefix = "c", small.zero = 0.01)
    }
    pgx$samples$cluster <- idx
  }

  if (cluster.contrasts) {
    ## Add cluster contrasts
    Y <- pgx$samples[, "cluster", drop = FALSE]
    if (length(unique(Y)) > 1) {
      ct <- make_direct_contrasts(Y, ref = "others")
      ctx <- contrast_as_labels(ct$exp.matrix)
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
    pgx$X <- log_cpm(pgx$counts, total = 1e6, prior = 1)
  }

  if (!all(dim(pgx$X) == dim(pgx$counts))) {
    stop("dimensions of X and counts do not match\n")
  }

  if (do.clustergenes) {
    pgx <- cluster_genes(pgx, methods = "umap", dims = c(2, 3), level = "gene")
  }

  return(pgx)
}



#' Compute the pgx object
#'
#' This function works on a pgx object that has already been created from
#' counts/samples/contrasts data using the `pgx_create()` function.
#'
#' @param ngs datatype
#' @param max.genes datatype
#' @param max.genesets datatype
#' @param gx.methods datatype
#' @param gset.methods datatype
#' @param do.cluster datatype
#' @param use.design datatype
#' @param prune.samples datatype
#' @param extra.methods datatype
#' @param lib.dir datatype
#' @param progress datatype
#'
#' @return a list (representing a pgx object)
#' @export
#'
#' @examples
#' x <- 1
compute_pgx <- function(ngs,
                        max.genes = 19999,
                        max.genesets = 9999,
                        gx.methods = c("ttest.welch", "trend.limma", "edger.qlf"),
                        gset.methods = c("fisher", "gsva", "fgsea"),
                        do.cluster = TRUE,
                        use.design = TRUE,
                        prune.samples = FALSE,
                        extra.methods = c(
                          "meta.go", "deconv", "infer",
                          "drugs", "wordcloud"
                        ),
                        lib.dir = "../lib",
                        progress = NULL) {
  ## ======================================================================
  ## ======================================================================
  ## ======================================================================

  if (!"contrasts" %in% names(ngs)) {
    stop("No contrasts found in pgx object")
  }

  ## make proper contrast matrix
  contr.matrix <- ngs$contrasts
  contr.values <- unique(as.vector(contr.matrix))
  is.numcontrast <- all(contr.values %in% c(NA, -1, 0, 1))
  is.numcontrast <- is.numcontrast && (-1 %in% contr.values) && (1 %in% contr.values)
  if (!is.numcontrast) {
    contr.matrix <- make_contrasts_from_label_matrix(contr.matrix)
    contr.matrix <- sign(contr.matrix) ## sign is fine
  }

  ## select valid contrasts
  sel <- Matrix::colSums(contr.matrix == -1) > 0 & Matrix::colSums(contr.matrix == 1) > 0
  contr.matrix <- contr.matrix[, sel, drop = FALSE]

  ## ======================================================================
  ngs$timings <- c()
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
  ngs <- test_genes(
    ngs, contr.matrix,
    max.features = max.genes,
    test.methods = gx.methods,
    use.design = use.design,
    prune.samples = prune.samples
  )

  ## ------------------ gene set tests -----------------------
  ngs <- test_genesets(
    ngs,
    max.features = max.genesets,
    test.methods = gset.methods,
    lib.dir = lib.dir
  )

  if (do.cluster) {
    ## gsetX not ready!!
    ngs <- cluster_genes(
      ngs,
      methods = "umap",
      dims = c(2, 3),
      level = "geneset"
    )
  }

  ## ------------------ extra analyses ---------------------
  ngs <- compute_extra(ngs, extra = extra.methods, lib.dir = lib.dir)

  return(ngs)
}


#' Alter the pgx options list
#'
#' This function creates the options list which will be used (eventually) to
#' control how the pgx object is created and computed.
#'
#' @param is.logx boolean.
#' @param batch.correct boolean.
#' @param auto.scale boolean.
#' @param filter.genes boolean.
#' @param prune.samples boolean.
#' @param only.known boolean.
#' @param only.hugo boolean.
#' @param convert.hugo boolean.
#' @param do.cluster boolean.
#' @param cluster.contrasts boolean.
#' @param do.clustergenes boolean.
#' @param only.proteincoding boolean.
#'
#' @return list of options
#' @export
#'
#' @examples
#' opt <- pgx_options()
pgx_options <- function(is.logx = NULL,
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
                        only.proteincoding = TRUE) {
  opts <- list(
    'is.logx' = is.logx,
    'batch.correct' = batch.correct,
    'auto.scale' = auto.scale,
    'filter.genes' = filter.genes,
    'prune.samples' = prune.samples,
    'only.known' = only.known,
    'only.hugo' = only.hugo,
    'convert.hugo' = convert.hugo,
    'do.cluster' = do.cluster,
    'cluster.contrasts' = cluster.contrasts,
    'do.clustergenes' = do.clustergenes,
    'only.proteincoding' = only.proteincoding
  )
  return(opts)
}
