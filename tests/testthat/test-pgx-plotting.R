## Regression cover for the plotting functions against the counts/X sample split.
##
## `pgx.createPGX()` stopped forcing the sample sets of `counts` and `X` equal
## (D-24): `counts`, `samples` and `contrasts` span the whole upload, `X` and the
## design span whatever outlier removal left. Both functions below indexed an
## X-shaped quantity by `rownames(pgx$samples)` and went out of bounds.
##
## The fixture runs the REAL removal -- playbase.preprocess drops two samples of
## playbase::COUNTS at threshold 2 -- so the split is produced, not asserted into
## existence, and every test asserts the drop happened before it tests anything.
## The fixture itself is helper-outlier-pgx.R's `outlier_pgx()`.

## the plotly data behind a scatter plot, one row per point, named by sample
scatter_points <- function(...) {
  suppressMessages(playbase::pgx.scatterPlot(..., plotlib = "plotly"))$x$visdat[[2]]()
}

## the plotly traces behind an expression plot
expression_traces <- function(...) {
  p <- suppressMessages(playbase::pgx.plotExpression(..., plotlib = "plotly"))
  plotly::plotly_build(p)$x$data
}

test_that("the fixture really drops a sample", {
  pgx <- outlier_pgx()
  expect_length(setdiff(rownames(pgx$samples), colnames(pgx$X)), 2)
  expect_equal(ncol(pgx$counts), nrow(pgx$samples))
  expect_lt(ncol(pgx$X), nrow(pgx$samples))
  expect_equal(nrow(pgx$model.parameters$exp.matrix), ncol(pgx$X))
})

test_that("pgx.scatterPlot draws X's samples when an outlier was dropped", {
  pgx <- outlier_pgx()
  dropped <- setdiff(rownames(pgx$samples), colnames(pgx$X))
  expect_gt(length(dropped), 0)
  gene <- rownames(pgx$X)[1]
  comp <- colnames(pgx$model.parameters$exp.matrix)[1]

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(suppressMessages(playbase::pgx.scatterPlot(pgx, pheno = "group")))
  expect_no_error(suppressMessages(playbase::pgx.scatterPlot(pgx, gene = gene)))
  expect_no_error(suppressMessages(playbase::pgx.scatterPlot(pgx, contrast = comp)))
  expect_no_error(suppressMessages(playbase::pgx.scatterPlot(pgx, geneset = "GS1")))

  pts <- scatter_points(pgx, pheno = "group")
  expect_setequal(rownames(pts), colnames(pgx$X))
  expect_false(any(dropped %in% rownames(pts)))
  expect_equal(as.character(pts$value), as.character(pgx$samples[rownames(pts), "group"]))

  pts <- scatter_points(pgx, gene = gene)
  expect_equal(as.numeric(pts$value), as.numeric(pgx$X[gene, rownames(pts)]),
    tolerance = 1e-3 ## the hover value is rounded for display
  )
})

test_that("pgx.scatterPlot rejects positions that are not sample-shaped", {
  pgx <- outlier_pgx()
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  gene.pos <- matrix(0, nrow = 4, ncol = 2, dimnames = list(rownames(pgx$X)[1:4], NULL))
  expect_error(
    suppressMessages(playbase::pgx.scatterPlot(pgx, pos = gene.pos, pheno = "group")),
    "dimension mismatch"
  )
  expect_error(
    suppressMessages(playbase::pgx.scatterPlot(pgx, pos = unname(pgx$tsne2d), pheno = "group")),
    "dimension mismatch"
  )
})

test_that("pgx.plotExpression draws X's samples when an outlier was dropped", {
  pgx <- outlier_pgx()
  dropped <- setdiff(rownames(pgx$samples), colnames(pgx$X))
  expect_gt(length(dropped), 0)
  gene <- rownames(pgx$X)[1]
  comp <- colnames(pgx$model.parameters$exp.matrix)[1]

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(suppressMessages(playbase::pgx.plotExpression(pgx, gene, comp)))
  expect_no_error(suppressMessages(playbase::pgx.plotExpression(pgx, gene, comp, grouped = TRUE)))
  expect_no_error(suppressMessages(playbase::pgx.plotExpression(pgx, gene, comp, showothers = FALSE)))
  expect_no_error(suppressMessages(
    playbase::pgx.plotExpression(pgx, "GS1", comp, level = "geneset")
  ))

  bars <- expression_traces(pgx, probe = gene, comp = comp)[[1]]
  drawn <- as.character(bars$x)
  expect_setequal(drawn, colnames(pgx$X))
  expect_false(any(dropped %in% drawn))
  expect_equal(as.numeric(bars$y), as.numeric(pgx$X[gene, drawn]))

  bars <- expression_traces(pgx, probe = "GS1", comp = comp, level = "geneset")[[1]]
  expect_equal(as.numeric(bars$y), as.numeric(pgx$gsetX["GS1", as.character(bars$x)]))
})

test_that("pgx.plotExpression reads group names off the contrast, not a recycled index", {
  ## pgx$contrasts spans the whole upload and the index into it is built on X's
  ## samples, so leaving it unsubset recycles silently instead of erroring.
  pgx <- outlier_pgx()
  expect_gt(nrow(pgx$contrasts), nrow(pgx$model.parameters$exp.matrix))
  gene <- rownames(pgx$X)[1]
  comp <- colnames(pgx$model.parameters$exp.matrix)[1]
  ct <- pgx$model.parameters$exp.matrix[, comp]
  labels <- playbase::contrastAsLabels(pgx$contrasts)[, comp]

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  bars <- expression_traces(pgx, probe = gene, comp = comp, grouped = TRUE)[[1]]
  groups <- as.character(bars$x)

  ## "other" is the bucket for samples this contrast does not compare
  expect_setequal(setdiff(groups, "other"), c(
    unique(labels[names(which(ct < 0))]),
    unique(labels[names(which(ct > 0))])
  ))
  ## and each bar is the mean over the samples that contrast side actually names
  members <- list(
    names(which(ct < 0)), names(which(ct > 0)), names(which(ct == 0))
  )
  names(members) <- c(
    unique(labels[names(which(ct < 0))]),
    unique(labels[names(which(ct > 0))]),
    "other"
  )
  expect_equal(
    as.numeric(bars$y),
    vapply(members[groups], function(s) mean(pgx$X[gene, s]), numeric(1)),
    ignore_attr = TRUE
  )
})
