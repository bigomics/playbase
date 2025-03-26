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
#' @param use.design A logical value indicating whether to use the design matrix in the analysis.
#' @param prune.samples A logical value indicating whether to prune samples with missing data.
#' @param remove.outputs A logical value indicating whether to remove intermediate outputs. 
#' @return An updated object with gene test results.
#' @export
compute_testGenes <- function(pgx,
                              contr.matrix,
                              max.features = 1000,
                              test.methods = c("trend.limma", "deseq2.wald", "edger.qlf"),
                              custom_fc = NULL,
                              use.design = FALSE,
                              prune.samples = TRUE,
                              remove.outputs = TRUE,
                              timeseries = FALSE
                              ) {

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

  ## sanity check
  if(NCOL(contr.matrix)==0) {
    stop("[compute_testGenes] FATAL: zero contrasts")
  }

  stat.group <- NULL
  if (use.design) {
    message("[compute_testGenes] detecting stat groups...")
    stat.group <- pgx.getConditions(contr.matrix, nmax = 0) ## !!!
    names(stat.group) <- rownames(contr.matrix)
    nlev <- length(unique(stat.group))
    if (nlev >= nrow(contr.matrix)) {
      message("[compute_testGenes] cannot use groups, switching to no design")
      use.design <- FALSE
    }
  }

  if (use.design) {
    message("[compute_testGenes] contrasts on groups (use design)")
    ## convert sample-wise contrasts to group-wise contrasts
    stat0 <- sort(unique(stat.group))
    contr.matrix <- contr.matrix[match(stat0, stat.group), , drop = FALSE]
    rownames(contr.matrix) <- stat0
  } else {
    message("[compute_testGenes] contrasts on samples (no design)")
    stat.group <- rownames(contr.matrix)
    names(stat.group) <- rownames(contr.matrix)
  }
  
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
  no.design <- all(stat.group %in% rownames(pgx$samples)) ## sample-wise design
  design <- NULL

  if (no.design || !use.design) {
    message("[compute_testGenes] 6 : no design matrix ")
    ## SAMPLE-WISE DESIGN
    design <- NULL
    exp.matrix <- contr.matrix
  } else {
    message("[compute_testGenes] 6 : creating model design matrix ")

    ## GROUP DESIGN
    notk <- which(!stat.group %in% rownames(contr.matrix))
    if (length(notk)) stat.group[notk] <- "_"
    design <- stats::model.matrix(~ 0 + stat.group) ## clean design no batch effects...
    colnames(design) <- sub("^stat.group", "", colnames(design))
    if (is.null(names(stat.group))) {
      stop("[compute_testGenes] FATAL:: stat.group must have names")
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
  message("Testing differential expresssion methods: ", paste(methods, collapse = ", "))
  PRIOR.CPM <- 1

  ## Run all test methods
  message("[compute_testGenes] 12 : start fitting... ")
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
  #remove.outputs = FALSE
  if (remove.outputs) pgx$gx.meta$outputs <- NULL

  message("[compute_testGenes] done!")

  return(pgx)
}


