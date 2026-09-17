## A pgx object in which counts retain the upload while X follows outlier removal.
##
## The removal is the REAL one -- playbase.preprocess drops two samples of
## playbase::COUNTS at threshold 2 -- so the split is produced rather than
## asserted into existence. Built once and reused; nothing mutates it.
outlier_pgx <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) {
      return(cached)
    }
    samples <- playbase::SAMPLES
    pp <- suppressMessages(playbase.preprocess::pgx.preprocess(
      counts = as.matrix(playbase::COUNTS),
      samples = samples,
      contrasts = playbase::CONTRASTS,
      options = list(
        norm_method = "CPM",
        remove_outliers = TRUE,
        outlier_threshold = 2,
        impute = FALSE
      )
    ))

    ## the design, cut to X's samples where pgx.computePGX() cuts it
    labels <- suppressMessages(
      playbase::contrasts.convertToLabelMatrix(playbase::CONTRASTS, samples)
    )
    contr <- sign(playbase::makeContrastsFromLabelMatrix(labels))
    contr <- contr[intersect(rownames(contr), colnames(pp$X)), , drop = FALSE]

    pgx <- list(
      counts = pp$counts,
      X = pp$X,
      samples = samples,
      contrasts = labels,
      organism = "Human",
      datatype = "RNA-seq",
      model.parameters = list(exp.matrix = contr, contr.matrix = contr)
    )
    ## gene sets live on X's samples, as compute_testGenesets() leaves them
    set.seed(1)
    pgx$gsetX <- matrix(
      stats::rnorm(2 * ncol(pp$X)),
      nrow = 2,
      dimnames = list(c("GS1", "GS2"), colnames(pp$X))
    )
    cached <<- suppressMessages(
      playbase::pgx.clusterSamples(pgx, methods = c("pca", "tsne"), dims = 2)
    )
    cached
  }
})

## The same split on a multi-omics object carrying MOFA factorizations, which is
## what `pgx.compute_importance(multiomics = 1)` needs to reach its own branch.
##
## Two real omics blocks are cut out of playbase::COUNTS and prefixed the way
## mofa.split_data() reads them. The removal is again the REAL one, here at the
## SHIPPED DEFAULT outlier_threshold = 3, which drops one sample. The
## factorizations are real too -- `pca` and `nmf2` are the two cheapest kernels,
## and two are needed because mofa.compute_meta_importance() drops to a vector
## on a single one. Built once and reused; nothing mutates it.
outlier_mofa_pgx <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) {
      return(cached)
    }
    counts <- as.matrix(playbase::COUNTS)
    samples <- playbase::SAMPLES
    half <- floor(nrow(counts) / 2)
    gx <- counts[1:half, , drop = FALSE]
    px <- counts[(half + 1):nrow(counts), , drop = FALSE]
    rownames(gx) <- paste0("gx:", rownames(gx))
    rownames(px) <- paste0("px:", rownames(px))

    pp <- suppressMessages(playbase.preprocess::pgx.preprocess(
      counts = rbind(gx, px),
      samples = samples,
      contrasts = playbase::CONTRASTS,
      options = list(
        norm_method = "CPM",
        remove_outliers = TRUE,
        outlier_threshold = 3,
        impute = FALSE
      )
    ))

    labels <- suppressMessages(
      playbase::contrasts.convertToLabelMatrix(playbase::CONTRASTS, samples)
    )
    contr <- sign(playbase::makeContrastsFromLabelMatrix(labels))
    contr <- contr[intersect(rownames(contr), colnames(pp$X)), , drop = FALSE]

    pgx <- list(
      counts = pp$counts,
      X = pp$X,
      samples = samples,
      contrasts = labels,
      organism = "Human",
      datatype = "multi-omics",
      model.parameters = list(exp.matrix = contr, contr.matrix = contr)
    )
    ## as pgx.computePGX() leaves them on a multi-omics object
    pgx$mofa <- list(
      factorizations = suppressMessages(suppressWarnings(
        mofa.compute_factorizations(
          ## internal
          playbase::mofa.split_data(pgx$X),
          samples[colnames(pgx$X), , drop = FALSE],
          labels[colnames(pgx$X), , drop = FALSE],
          numfactors = 4,
          kernels = c("pca", "nmf2")
        )
      ))
    )
    cached <<- pgx
    cached
  }
})
