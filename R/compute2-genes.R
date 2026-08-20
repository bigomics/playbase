##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

#' Compute Test Genes
#' This function computes gene tests based on the input data and contrast matrix.
#' It performs different test methods depending on whether the data is single-omics or multi-omics.
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

  ## Check parameters, decide group level
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

  ## Normalize contrast matrix to zero mean and signed sums to one
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
    group = stat.group,
    ## imputation method recorded upstream (createPGX), surfaced for the
    ## deterministic AI-report methods block; NULL when no imputation ran.
    impute_method = pgx$impute_method %||% "none (no imputation applied)"
  )
  pgx$model.parameters <- model.parameters

  if (is.null(names(stat.group))) {
    stop("[compute_testGenes] FATAL2:: stat.group must have names")
  }

  ## Conform data matrices
  ## notice original counts will not be affected
  ss <- names(stat.group)
  gg <- intersect(rownames(pgx$X), rownames(pgx$counts))
  counts <- pgx$counts[gg, ss, drop = FALSE]
  samples <- pgx$samples[ss, ]
  X <- pgx$X[gg, ss, drop = FALSE]

  methods <- test.methods
  message("Testing differential expression methods: ", paste(methods, collapse = ", "))
  PRIOR.CPM <- 1

  ## Methylomics: fill gx.meta rather than fit it.
  ##
  ## The slot is required - pgx.checkObject() lists it, pgxinfo.updateDatasetFolder()
  ## skips any pgx failing that check, and opg_server.R reads gx.meta$meta[[1]]$fc
  ## unguarded - but nothing renders it. Methylomics opens the standalone Methylome
  ## app, which refits limma itself from pgx$X with the user's own covariates,
  ## masking and SVA. Fitting here computes a model no screen displays, and the
  ## expensive part is not limma: it is betaToM() over the whole probe matrix,
  ## which on an 850K EPIC array is two more copies of the data.
  ##
  ## The fill is neutral, not random: zero effect, q = 1. This dataset therefore
  ## contributes no fold changes to the cross-dataset FC index, which is honest -
  ## no differential test was run - where invented values would not be.
  if (!is.null(pgx$datatype) && pgx$datatype == "methylomics") {
    message("[compute_testGenes] methylomics: filling gx.meta, not fitting")
    ctd <- colnames(contr.matrix)
    n <- nrow(X)
    one_col <- function(v) {
      m <- matrix(v, nrow = n, ncol = 1,
                  dimnames = list(rownames(X), "not.fitted"))
      I(m)
    }
    stub <- lapply(ctd, function(k) {
      data.frame(
        meta.fx = rep(0, n), meta.p = rep(1, n), meta.q = rep(1, n),
        avg.0 = rep(0, n), avg.1 = rep(0, n),
        fc = one_col(0), p = one_col(1), q = one_col(1),
        row.names = rownames(X), check.names = FALSE
      )
    })
    names(stub) <- ctd
    pgx$model.parameters <- model.parameters
    ## X stays on the beta scale here. compute_testGenes normally hands back
    ## M-values and pgx.computePGX converts them; that conversion now checks the
    ## range, so leaving beta is safe.
    pgx$X <- X
    pgx$gx.meta <- list(
      meta = stub,
      meta.covs = NULL,
      sig.counts = NULL
    )
    message("[compute_testGenes] done (methylomics, not fitted)")
    return(pgx)
  }

  if (!is.null(pgx$datatype) & pgx$datatype == "methylomics") {

    require_epigenetics()

    if ("Differentially methylated regions" %in% pgx$dma) {

      message("[playbase::compute_testGenes] Methylomics: DMRs...")

      vv <- range(counts, na.rm = TRUE)
      is.beta <- (vv[1] >= 0 & vv[2] <= 1) ## original counts
      MG <- playbase.epigenetics::mergeCpG(data = counts, genes = pgx$genes)
      counts <- playbase.epigenetics::betaToM(MG$data)
      if (is.beta) pgx$counts=MG$data else pgx$counts=counts ## restore as original (beta or m)
      rm(MG); gc()

      MG <- playbase.epigenetics::mergeCpG(data = X, genes = pgx$genes)
      X <- playbase.epigenetics::betaToM(MG$data)
      pgx$genes <- MG$genes
      rm(MG); gc()
      
    } else {
      counts <- X <- playbase.epigenetics::betaToM(counts)
    }
  }

  ## Run methods
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
    prior.cpm = PRIOR.CPM,
    remove.batch = FALSE,
    conform.output = TRUE,
    do.filter = FALSE,
    correct.AveExpr = TRUE,
    custom = custom_fc,
    custom.name = NULL,
    timeseries = timeseries
  )

  message("[compute_testGenes]: fitting completed!")
  
  ## Set default matrices
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
