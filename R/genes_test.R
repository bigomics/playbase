
# compute.testGenes
test_genes <- function(pgx, contr.matrix, max.features = 1000,
                       test.methods = c("trend.limma", "deseq2.wald", "edger.qlf"),
                       use.design = TRUE, prune.samples = FALSE,
                       remove.outputs = TRUE) {
  single.omics <- mean(grepl("\\[", rownames(pgx$counts))) < 0.1
  data.types <- unique(gsub("\\[|\\].*", "", rownames(pgx$counts)))
  if (single.omics || length(data.types) == 1) {
    ## single-omics, no missing values
    pgx <- test_genes_singleomics(
      pgx = pgx,
      contr.matrix = contr.matrix,
      max.features = max.features,
      test.methods = test.methods,
      use.design = use.design,
      prune.samples = prune.samples,
      remove.outputs = remove.outputs
    )
  } else {
    ## multi-omics, missing values allowed
    pgx <- test_genes_multiomics(
      pgx = pgx, ## type is inferred
      contr.matrix = contr.matrix,
      max.features = max.features,
      test.methods = test.methods,
      use.design = use.design,
      remove.outputs = remove.outputs
    )
  }
  return(pgx)
}

# compute.testGenesSingleOmics
test_genes_singleomics <- function(pgx, contr.matrix, max.features = 1000,
                                   filter.low = TRUE, remove.outputs = TRUE,
                                   use.design = TRUE, prune.samples = FALSE,
                                   test.methods = c("trend.limma", "deseq2.wald", "edger.qlf")) {
  contr.matrix0 <- contr.matrix ## SAVE

  ## -----------------------------------------------------------------------------
  ## Check parameters, decide group level
  ## -----------------------------------------------------------------------------
  if (!("counts" %in% names(pgx))) {
    stop("[test_genes_singleomics] FATAL: cannot find counts in pgx object")
  }
  if (!("X" %in% names(pgx))) {
    stop("[test_genes_singleomics] FATAL: cannot find normalized expression X in pgx object")
  }

  is.expmatrix <- all(rownames(contr.matrix) %in% rownames(pgx$samples))
  if (!is.expmatrix) {
    stop("[test_genes_singleomics] FATAL: contrast must be sample-wise")
  }

  stat.group <- NULL
  if (use.design) {
    stat.group <- get_conditions(contr.matrix, nmax = 0) ## !!!
    names(stat.group) <- rownames(contr.matrix)
    nlev <- length(unique(stat.group))
    if (nlev >= nrow(contr.matrix)) {
      use.design <- FALSE
    }
  }

  if (use.design) {
    ## stat.group = rownames(pgx$samples)
    ## convert sample-wise contrasts to group-wise contrasts
    stat0 <- sort(unique(stat.group))
    contr.matrix <- contr.matrix[match(stat0, stat.group), , drop = FALSE]
    rownames(contr.matrix) <- stat0
  } else if (!use.design) {
    stat.group <- rownames(contr.matrix)
    names(stat.group) <- rownames(contr.matrix)
  }

  ## take out any empty comparisons
  sel <- which(Matrix::colSums(contr.matrix > 0) & Matrix::colSums(contr.matrix < 0))
  contr.matrix <- contr.matrix[, sel, drop = FALSE]
  contr.matrix[is.na(contr.matrix)] <- 0

  ## -----------------------------------------------------------------------------
  ## normalize contrast matrix to zero mean and signed sums to one
  ## -----------------------------------------------------------------------------
  ## normalize?? why??
  for (i in 1:ncol(contr.matrix)) {
    m <- contr.matrix[, i]
    m[is.na(m)] <- 0
    contr.matrix[, i] <- 1 * (m > 0) / sum(m > 0) - 1 * (m < 0) / sum(m < 0)
  }

  ## -----------------------------------------------------------------------------
  ## create design matrix from defined contrasts (group or clusters)
  ## -----------------------------------------------------------------------------

  no.design <- all(stat.group %in% rownames(pgx$samples)) ## sample-wise design
  design <- NULL

  if (no.design || !use.design) {
    ## SAMPLE-WISE DESIGN
    design <- NULL
    exp.matrix <- contr.matrix
  } else {
    ## GROUP DESIGN
    notk <- which(!stat.group %in% rownames(contr.matrix))
    if (length(notk)) {
      stat.group[notk] <- "_"
    }
    design <- model.matrix(~ 0 + stat.group) ## clean design no batch effects...
    colnames(design) <- sub("^stat.group", "", colnames(design))
    if (is.null(names(stat.group))) {
      stop("[test_genes_singleomics] FATAL:: stat.group must have names")
    }
    rownames(design) <- names(stat.group)

    ## make sure matrix align and compute experiment matrix
    design <- design[, match(rownames(contr.matrix), colnames(design)), drop = FALSE]
    colnames(design) <- rownames(contr.matrix)
    exp.matrix <- (design %*% contr.matrix)

    ## check contrasts for sample sizes (at least 2 in each group) and
    ## remove otherwise
    keep <- rep(TRUE, ncol(contr.matrix))
    keep <- (Matrix::colSums(exp.matrix > 0) >= 1 & Matrix::colSums(exp.matrix < 0) >= 1)
    contr.matrix <- contr.matrix[, keep, drop = FALSE]
    exp.matrix <- exp.matrix[, keep, drop = FALSE]
  }

  model.parameters <- list(
    design = design,
    contr.matrix = contr.matrix,
    exp.matrix = exp.matrix,
    group = stat.group
  )
  pgx$model.parameters <- model.parameters

  ## -----------------------------------------------------------------------------
  ## Filter genes
  ## -----------------------------------------------------------------------------
  if (is.null(names(stat.group))) {
    stop("[compute.testGenesSingleOmics] FATAL2:: stat.group must have names")
  }

  ## notice original counts will not be affected
  ss <- names(stat.group)
  gg <- rownames(pgx$counts)
  if (!is.null(pgx$X)) gg <- intersect(gg, rownames(pgx$X))
  counts <- pgx$counts[gg, ss, drop = FALSE]
  genes <- pgx$genes[gg, ]
  samples <- pgx$samples[ss, ]

  ## Rescale if too low. Often EdgeR/DeSeq can give errors of total counts
  ## are too low. Happens often with single-cell (10x?). We rescale
  ## to a minimum of 1 million counts (CPM)

  mean.counts <- mean(Matrix::colSums(counts, na.rm = TRUE))
  if (mean.counts < 1e6) {
    warning("low total counts = ", mean.counts, "\n")
    counts <- counts * 1e6 / mean.counts
  }
  mean(Matrix::colSums(counts, na.rm = TRUE))

  ## prefiltering for low-expressed genes (recommended for edgeR and
  ## DEseq2). Require at least in 2 or 1% of total. Specify the
  ## PRIOR CPM amount to regularize the counts and filter genes
  PRIOR.CPM <- 1
  if (filter.low) {
    PRIOR.CPM <- 1
    AT.LEAST <- ceiling(pmax(2, 0.01 * ncol(counts)))
    keep <- (rowSums(edgeR::cpm(counts) > PRIOR.CPM, na.rm = TRUE) >= AT.LEAST)
    pgx$filtered <- NULL
    pgx$filtered[["low.expressed"]] <- paste(rownames(counts)[which(!keep)], collapse = ";")
    counts <- counts[which(keep), , drop = FALSE]
    genes <- genes[which(keep), , drop = FALSE]
  }

  ## -----------------------------------------------------------------------------
  ## Shrink number of genes before testing (highest SD/var)
  ## -----------------------------------------------------------------------------
  if (is.null(max.features)) max.features <- -1
  if (max.features > 0 && nrow(counts) > max.features) {
    logcpm <- log_cpm(counts, total = NULL)
    sdx <- apply(logcpm, 1, sd)
    jj <- Matrix::head(order(-sdx), max.features) ## how many genes?
    jj0 <- setdiff(1:nrow(counts), jj)
    pgx$filtered[["low.variance"]] <- paste(rownames(counts)[jj0], collapse = ";")
    counts <- counts[jj, ]
    genes <- genes[jj, ]
  }
  genes <- genes[, c("gene_name", "gene_title")]

  ## -----------------------------------------------------------------------------
  ## Do the fitting
  ## -----------------------------------------------------------------------------
  methods <- test.methods

  ## Run all test methods
  ##
  X <- pgx$X[rownames(counts), colnames(counts)]

  gx.meta <- fit_contrasts_with_all_methods(
    counts = counts,
    X = X, ## type = type,
    samples = samples,
    genes = NULL, ## genes=genes,
    methods = methods,
    design = design,
    contr.matrix = contr.matrix,
    prune.samples = prune.samples,
    prior.cpm = PRIOR.CPM, ## prior count regularization
    ## quantile.normalize = TRUE,  ## only for logCPM???
    remove.batch = FALSE, ## we do explicit batch correction instead
    conform.output = TRUE,
    do.filter = FALSE,
    correct.AveExpr = TRUE,
    custom = NULL, custom.name = NULL
  )

  ## --------------------------------------------------------------------------------
  ## set default matrices
  ## --------------------------------------------------------------------------------

  rownames(gx.meta$timings) <- paste0("[test.genes]", rownames(gx.meta$timings))
  pgx$timings <- rbind(pgx$timings, gx.meta$timings)
  gx.meta$timings <- NULL
  gx.meta$X <- NULL
  pgx$model.parameters <- model.parameters
  pgx$gx.meta <- gx.meta
  pgx$X <- X ## replace with filtered

  ## remove large outputs... (uncomment if needed!!!)
  if (remove.outputs) {
    pgx$gx.meta$outputs <- NULL
  }

  return(pgx)
}

# compute.testGenesMultiOmics
test_genes_multiomics <- function(pgx, contr.matrix, max.features = 1000,
                                  test.methods = c("trend.limma", "deseq2.wald", "edger.qlf"),
                                  use.design = TRUE, prune.samples = FALSE,
                                  remove.outputs = TRUE) {
  pgx$gx.meta <- NULL
  pgx$model.parameters <- NULL
  pgx$gx.meta$meta <- vector("list", ncol(contr.matrix))
  pgx$X <- c()
  pgx$timings <- c()
  for (j in 1:4) {
    nk <- ncol(contr.matrix)
    pgx$gx.meta$sig.counts[[j]] <- vector("list", nk)
  }

  data.type <- gsub("\\[|\\].*", "", rownames(pgx$counts))
  data.types <- unique(data.type)
  dt <- data.types[1]
  for (dt in data.types) {
    ## get data block
    pgx1 <- pgx
    jj <- which(data.type == dt)
    pgx1$counts <- pgx1$counts[jj, ]
    pgx1$genes <- pgx1$genes[jj, ]

    ## determine if datatype are counts or not
    type <- "not.counts"
    if (min(pgx1$counts, na.rm = TRUE) >= 0 &&
      max(pgx1$counts, na.rm = TRUE) >= 50) {
      type <- "counts"
    }

    ## do test
    pgx1 <- test_genes_singleomics(
      pgx = pgx1, type = type,
      contr.matrix = contr.matrix,
      max.features = max.features,
      test.methods = test.methods
    )

    ## copy results
    pgx$model.parameters <- pgx1$model.parameters
    for (k in 1:ncol(contr.matrix)) {
      pgx$gx.meta$meta[[k]] <- rbind(
        pgx$gx.meta$meta[[k]],
        pgx1$gx.meta$meta[[k]]
      )
    }
    names(pgx$gx.meta$meta) <- names(pgx1$gx.meta$meta)
    for (j in 1:4) {
      nk <- ncol(contr.matrix)
      for (k in 1:nk) {
        cnt1 <- pgx1$gx.meta$sig.counts[[j]][[k]]
        cnt0 <- pgx$gx.meta$sig.counts[[j]][[k]]
        rownames(cnt1) <- paste0("[", dt, "]", rownames(cnt1))
        pgx$gx.meta$sig.counts[[j]][[k]] <- rbind(cnt0, cnt1)
      }
      names(pgx$gx.meta$sig.counts[[j]]) <- names(pgx1$gx.meta$sig.counts[[j]])
    }
    names(pgx$gx.meta$sig.counts) <- names(pgx1$gx.meta$sig.counts)
    pgx$timings <- rbind(pgx$timings, pgx1$timings)
    pgx$X <- rbind(pgx$X, pgx1$X)
  }

  gg <- rownames(pgx$counts)
  pgx$X <- pgx$X[match(gg, rownames(pgx$X)), ]
  pgx$model.parameters <- pgx1$model.parameters
  return(pgx)
}
