## A pgx object shaped the way pgx.createPGX() + pgx.computePGX() leave one when
## outlier removal drops a sample: `counts`, `samples` and `contrasts` span the
## whole upload, `X` and the design span whatever the removal left (D-24).
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
    pp <- suppressMessages(playbase::pgx.preprocess(
      counts = as.matrix(playbase::COUNTS),
      samples = samples,
      contrasts = playbase::CONTRASTS,
      options = list(
        datatype = "RNA-seq", norm_method = "CPM",
        remove_outliers = TRUE, outlier_threshold = 2, impute = FALSE
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
      nrow = 2, dimnames = list(c("GS1", "GS2"), colnames(pp$X))
    )
    cached <<- suppressMessages(
      playbase::pgx.clusterSamples(pgx, methods = c("pca", "tsne"), dims = 2)
    )
    cached
  }
})
