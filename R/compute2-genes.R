##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

#' Compute Test Genes
#'
#' This function computes gene tests based on the input data and contrast matrix.
#' It performs different test methods depending on whether the data is single-omics or multi-omics.
#' For single-omics data or when there is only one data type, it uses the \code{compute_testGenesSingleOmics} function.
#' For multi-omics data with the possibility of missing values, it uses the \code{compute_testGenesMultiOmics} function.
#'
#' @param pgx An object containing the input data for analysis.
#' @param contr.matrix The contrast matrix for the gene tests.
#' @param max.features The maximum number of features to consider in the gene tests.
#' @param test.methods The test methods to use for gene testing.
#' @param use.design A logical value indicating whether to use the design matrix in the analysis.
#' @param prune.samples A logical value indicating whether to prune samples with missing data.
#' @param remove.outputs A logical value indicating whether to remove intermediate outputs.
#' @return An updated object with gene test results.
#' @export
compute_testGenes <- function(pgx, contr.matrix, max.features = 1000,
                              test.methods = c("trend.limma", "deseq2.wald", "edger.qlf"),
                              use.design = TRUE, prune.samples = FALSE,
                              remove.outputs = TRUE) {
  single.omics <- mean(grepl("\\[", rownames(pgx$counts))) < 0.1
  single.omics
  data.types <- unique(gsub("\\[|\\].*", "", rownames(pgx$counts)))
  ## data.types
  if (single.omics || length(data.types) == 1) {
    ## single-omics, no missing values
    cat(">>> computing gene tests for SINGLE-OMICS\n")
    pgx <- compute_testGenesSingleOmics(
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
    cat(">>> computing gene tests for MULTI-OMICS\n")
    pgx <- compute_testGenesMultiOmics(
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

#' Compute Test Genes for Single-Omics Data
#'
#' This function computes gene tests for single-omics data based on the input data and contrast matrix.
#' It performs various steps such as normalization, creating a design matrix, filtering genes, and fitting the data using different test methods.
#'
#' @param pgx An object containing the input data for analysis.
#' @param contr.matrix The contrast matrix for the gene tests.
#' @param max.features The maximum number of features to consider in the gene tests.
#' @param filter.low A logical value indicating whether to filter low-expressed genes.
#' @param remove.outputs A logical value indicating whether to remove intermediate outputs.
#' @param use.design A logical value indicating whether to use the design matrix in the analysis.
#' @param prune.samples A logical value indicating whether to prune samples with missing data.
#' @param test.methods The test methods to use for gene testing.
#' @return An updated object with gene test results for single-omics data.
#' @export
compute_testGenesSingleOmics <- function(pgx, contr.matrix, max.features = 1000,
                                         filter.low = TRUE, remove.outputs = TRUE,
                                         use.design = TRUE, prune.samples = FALSE,
                                         test.methods = c("trend.limma", "deseq2.wald", "edger.qlf")) {
  ## -----------------------------------------------------------------------------
  ## Check parameters, decide group level
  ## -----------------------------------------------------------------------------
  if (!("counts" %in% names(pgx))) {
    stop("[compute_testGenesSingleOmics] FATAL: cannot find counts in pgx object")
  }
  if (!("X" %in% names(pgx))) {
    stop("[compute_testGenesSingleOmics] FATAL: cannot find normalized expression X in pgx object")
  }

  is.expmatrix <- all(rownames(contr.matrix) %in% rownames(pgx$samples))
  is.expmatrix
  if (!is.expmatrix) {
    stop("[compute_testGenesSingleOmics] FATAL: contrast must be sample-wise")
  }

  stat.group <- NULL
  if (use.design) {
    message("[compute_testGenesSingleOmics] detecting stat groups...")
    stat.group <- pgx.getConditions(contr.matrix, nmax = 0) ## !!!
    names(stat.group) <- rownames(contr.matrix)
    nlev <- length(unique(stat.group))
    nlev
    if (nlev >= nrow(contr.matrix)) {
      message("[compute_testGenesSingleOmics] cannot use groups, switching to no design")
      use.design <- FALSE
    }
  }

  if (use.design) {
    message("[compute_testGenesSingleOmics] contrasts on groups (use design)")
    ## convert sample-wise contrasts to group-wise contrasts
    message("replacing contrast matrix...")
    stat0 <- sort(unique(stat.group))
    contr.matrix <- contr.matrix[match(stat0, stat.group), , drop = FALSE]
    rownames(contr.matrix) <- stat0
  } else if (!use.design) {
    message("[compute_testGenesSingleOmics] contrasts on samples (no design)")
    stat.group <- rownames(contr.matrix)
    names(stat.group) <- rownames(contr.matrix)
  }

  message("[compute_testGenesSingleOmics] pruning unused contrasts")
  ## take out any empty comparisons
  sel <- which(Matrix::colSums(contr.matrix > 0) & Matrix::colSums(contr.matrix < 0))
  contr.matrix <- contr.matrix[, sel, drop = FALSE]
  contr.matrix[is.na(contr.matrix)] <- 0

  ## -----------------------------------------------------------------------------
  ## normalize contrast matrix to zero mean and signed sums to one
  ## -----------------------------------------------------------------------------
  ## normalize?? why??
  message("[compute_testGenesSingleOmics] normalizing contrasts")
  for (i in seq_len(ncol(contr.matrix))) {
    m <- contr.matrix[, i]
    m[is.na(m)] <- 0
    contr.matrix[, i] <- 1 * (m > 0) / sum(m > 0) - 1 * (m < 0) / sum(m < 0)
  }

  ## -----------------------------------------------------------------------------
  ## create design matrix from defined contrasts (group or clusters)
  ## -----------------------------------------------------------------------------

  no.design <- all(stat.group %in% rownames(pgx$samples)) ## sample-wise design
  design <- NULL
  no.design

  if (no.design || !use.design) {
    message("[compute_testGenesSingleOmics] 6 : no design matrix ")
    ## SAMPLE-WISE DESIGN
    design <- NULL
    exp.matrix <- contr.matrix
  } else {
    message("[compute_testGenesSingleOmics] 6 : creating model design matrix ")

    ## GROUP DESIGN
    notk <- which(!stat.group %in% rownames(contr.matrix))
    if (length(notk)) {
      stat.group[notk] <- "_"
    }
    design <- model.matrix(~ 0 + stat.group) ## clean design no batch effects...
    colnames(design) <- sub("^stat.group", "", colnames(design))
    if (is.null(names(stat.group))) {
      stop("[compute_testGenesSingleOmics] FATAL:: stat.group must have names")
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
    stop("[compute_testGenesSingleOmics] FATAL2:: stat.group must have names")
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
  mean.counts
  if (mean.counts < 1e6) {
    cat("[compute_testGenesSingleOmics] WARNING:: low total counts = ", mean.counts, "\n")
    cat("[compute_testGenesSingleOmics] applying global mean scaling to 1e6...\n")
    counts <- counts * 1e6 / mean.counts
  }

  ## prefiltering for low-expressed genes (recommended for edgeR and
  ## DEseq2). Require at least in 2 or 1% of total. Specify the
  ## PRIOR CPM amount to regularize the counts and filter genes
  PRIOR.CPM <- 1
  if (filter.low) {
    PRIOR.CPM <- 0.25
    PRIOR.CPM <- 1
    PRIOR.CPM
    AT.LEAST <- ceiling(pmax(2, 0.01 * ncol(counts)))
    cat("filtering for low-expressed genes: >", PRIOR.CPM, "CPM in >=", AT.LEAST, "samples\n")
    keep <- (rowSums(edgeR::cpm(counts) > PRIOR.CPM, na.rm = TRUE) >= AT.LEAST)
    pgx$filtered <- NULL
    pgx$filtered[["low.expressed"]] <-
      paste(rownames(counts)[which(!keep)], collapse = ";")
    counts <- counts[which(keep), , drop = FALSE]
    genes <- genes[which(keep), , drop = FALSE]
    cat("filtering out", sum(!keep), "low-expressed genes\n")
    cat("keeping", sum(keep), "expressed genes\n")
  }

  ## -----------------------------------------------------------------------------
  ## Shrink number of genes before testing (highest SD/var)
  ## -----------------------------------------------------------------------------
  if (is.null(max.features)) max.features <- -1
  if (max.features > 0 && nrow(counts) > max.features) {
    cat("shrinking data matrices: n=", max.features, "\n")
    logcpm <- logCPM(counts, total = NULL)
    sdx <- apply(logcpm, 1, sd)
    jj <- Matrix::head(order(-sdx), max.features) ## how many genes?
    jj0 <- setdiff(seq_len(nrow(counts)), jj)
    pgx$filtered[["low.variance"]] <- paste(rownames(counts)[jj0], collapse = ";")
    counts <- counts[jj, ]
    genes <- genes[jj, ]
  }
  genes <- genes[, c("gene_name", "gene_title")]

  ## -----------------------------------------------------------------------------
  ## Do the fitting
  ## -----------------------------------------------------------------------------
  methods <- test.methods
  methods
  cat(">>> Testing differential expressed genes (DEG) with methods:", methods, "\n")

  ## Run all test methods
  ##
  X <- pgx$X[rownames(counts), colnames(counts)]

  message("[compute_testGenesSingleOmics] 12 : start fitting... ")

  gx.meta <- ngs.fitContrastsWithAllMethods(
    counts = counts,
    X = X,
    samples = samples,
    genes = NULL,
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

  message("[compute_testGenesSingleOmics] 13 : fitting done!")

  names(gx.meta)
  names(gx.meta$outputs)
  print(gx.meta$timings)

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

  ## remove large outputs.
  if (remove.outputs) {
    pgx$gx.meta$outputs <- NULL
  }

  message("[compute_testGenesSingleOmics] done!")

  return(pgx)
}


#' Compute Multi-Omics Test Genes
#'
#' Computes test genes for multi-omics data.
#'
#' @param pgx A data object representing multi-omics data.
#' @param contr.matrix A contrast matrix specifying the experimental design.
#' @param max.features Maximum number of features to consider.
#' @param test.methods Methods to use for testing.
#' @param use.design Logical indicating whether to use the experimental design.
#' @param prune.samples Logical indicating whether to prune samples.
#' @param remove.outputs Logical indicating whether to remove outputs.
#'
#' @return The updated \code{pgx} object with computed test genes.
#'
#' @export
compute_testGenesMultiOmics <- function(pgx, contr.matrix, max.features = 1000,
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
  data.types
  dt <- "cn"
  dt <- "gx"
  dt <- data.types[1]
  dt
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
    dt
    type

    ## do test
    pgx1 <- compute_testGenesSingleOmics(
      pgx = pgx1,
      contr.matrix = contr.matrix,
      max.features = max.features,
      test.methods = test.methods
    )

    ## copy results
    pgx$model.parameters <- pgx1$model.parameters
    names(pgx1$gx.meta)
    for (k in seq_len(ncol(contr.matrix))) {
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
