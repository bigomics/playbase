##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

#' Compute Test Genes
#'
#' This function computes gene tests based on the input data and contrast matrix.
#' It performs different test methods depending on whether the data is single-omics or multi-omics.
#'
#' @param pgx An object containing the input data for analysis.
#' @param contr.matrix The contrast matrix for the gene tests.
#' @param max.features The maximum number of features to consider in the gene tests.
#' @param test.methods The test methods to use for gene testing.
#' @param prune.samples A logical value indicating whether to prune samples with missing data.
#' @param remove.outputs A logical value indicating whether to remove intermediate outputs.
#' @return An updated object with gene test results.
#' @export
compute_testGenes <- function(pgx,
                              contr.matrix,
                              max.features = 1000,
                              test.methods = c("trend.limma", "deseq2.wald", "edger.qlf"),
                              custom_fc = NULL,
                              prune.samples = TRUE,
                              remove.outputs = TRUE,
                              timeseries = FALSE) {

  message("[compute_testGenes] detecting stat groups...")

  ## -----------------------------------------------------------------------------
  ## Check parameters, decide group level
  ## -----------------------------------------------------------------------------
  if (!("counts" %in% names(pgx))) {
    stop("[compute_testGenes] FATAL: cannot find counts in pgx object")
  }

  if (!("X" %in% names(pgx))) {
    stop("[compute_testGenes] FATAL: cannot find normalized expression X in pgx object")
  }

  is.expmatrix <- all(rownames(contr.matrix) %in% rownames(pgx$samples))
  if (!is.expmatrix) {
    stop("[compute_testGenes] FATAL: contrast must be sample-wise")
  }

  if (NCOL(contr.matrix) == 0) {
    stop("[compute_testGenes] FATAL: zero contrasts")
  }

  stat.group <- NULL
  message("[compute_testGenes] contrasts on samples (no design)")
  stat.group <- rownames(contr.matrix)
  names(stat.group) <- rownames(contr.matrix)

  message("[compute_testGenes] pruning unused contrasts")
  ## take out any empty comparisons
  sel <- which(Matrix::colSums(contr.matrix > 0) & Matrix::colSums(contr.matrix < 0))
  contr.matrix <- contr.matrix[, sel, drop = FALSE]
  contr.matrix[is.na(contr.matrix)] <- 0

  ## -----------------------------------------------------------------------------
  ## normalize contrast matrix to zero mean and signed sums to one
  ## -----------------------------------------------------------------------------
  ## normalize?? why??
  message("[compute_testGenes] normalizing contrasts")
  for (i in seq_len(ncol(contr.matrix))) {
    m <- contr.matrix[, i]
    m[is.na(m)] <- 0
    contr.matrix[, i] <- 1 * (m > 0) / sum(m > 0) - 1 * (m < 0) / sum(m < 0)
  }

  ## -----------------------------------------------------------------------------
  ## create design matrix from defined contrasts (group or clusters)
  ## -----------------------------------------------------------------------------
  message("[compute_testGenes] no design matrix ")
  ## SAMPLE-WISE DESIGN
  design <- NULL
  exp.matrix <- contr.matrix

  model.parameters <- list(
    design = design,
    contr.matrix = contr.matrix,
    exp.matrix = exp.matrix,
    group = stat.group
  )
  pgx$model.parameters <- model.parameters

  if (is.null(names(stat.group))) {
    stop("[compute_testGenes] FATAL2:: stat.group must have names")
  }

  ## -----------------------------------------------------------------------------
  ## Conform data matrices
  ## -----------------------------------------------------------------------------
  ## notice original counts will not be affected
  ss <- names(stat.group)
  gg <- intersect(rownames(pgx$X), rownames(pgx$counts))
  counts <- pgx$counts[gg, ss, drop = FALSE]
  samples <- pgx$samples[ss, ]
  X <- pgx$X[gg, ss, drop = FALSE]

  ## -----------------------------------------------------------------------------
  ## Do the fitting
  ## -----------------------------------------------------------------------------
  methods <- test.methods
  message("Testing differential expression methods: ", paste(methods, collapse = ", "))
  PRIOR.CPM <- 1
  
  ## Run all test methods
  message("[compute_testGenes] start fitting... ")
  gx.meta <- ngs.fitContrastsWithAllMethods(
    counts = counts,
    X = X,
    samples = samples,
    genes = NULL,
    methods = methods,
    covariates = pgx$covariates,
    contr.matrix = contr.matrix,
    prune.samples = prune.samples,
    prior.cpm = PRIOR.CPM, ## prior count regularization
    remove.batch = FALSE, ## we do explicit batch correction instead
    conform.output = TRUE,
    do.filter = FALSE,
    correct.AveExpr = TRUE,
    custom = custom_fc,
    custom.name = NULL,
    timeseries = timeseries
  )

  message("[compute_testGenes]: fitting completed!")

  
  ## --------------------------------------------------------------------------------
  ## set default matrices
  ## --------------------------------------------------------------------------------
  rownames(gx.meta$timings) <- paste0("[test.genes]", rownames(gx.meta$timings))
  pgx$timings <- rbind(pgx$timings, gx.meta$timings)
  gx.meta$timings <- NULL
  gx.meta$X <- NULL
  pgx$model.parameters <- model.parameters
  pgx$X <- X ## adopt
  pgx$gx.meta <- gx.meta

  ## remove large outputs.
  if (remove.outputs) pgx$gx.meta$outputs <- NULL
  
  message("[compute_testGenes] done!")

  return(pgx)

}
