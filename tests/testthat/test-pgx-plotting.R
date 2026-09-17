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

## Review UNCONFIRMED cx5, site 1. `pgx.contrastScatter()` read the sides of a
## contrast as POSITIONS in the design's rows and then used them as POSITIONS in
## pgx$X's columns. The two axes agreed only because pgx.computePGX() built the
## design in the sample table's order and the sample table happened to be in
## X's order -- a coincidence nothing documented or tested. The design is now
## built in X's order, and the reader names its samples instead of counting.

test_that("pgx.contrastScatter names the contrast's samples instead of counting", {
  pgx <- outlier_pgx()
  comp <- colnames(pgx$model.parameters$exp.matrix)[1]
  ct <- pgx$model.parameters$exp.matrix[, comp]
  grp0 <- names(which(ct < 0))
  grp1 <- names(which(ct > 0))
  expect_gt(length(grp0), 0)
  expect_gt(length(grp1), 0)

  ## meta stats are all pgx.contrastScatter needs beyond the design
  gg <- rownames(pgx$X)
  pgx$gx.meta <- list(meta = stats::setNames(list(data.frame(
    meta.fx = seq_along(gg) / length(gg),
    meta.q = rep(0.01, length(gg)),
    row.names = gg
  )), comp))

  expected <- cbind(
    x0 = rowMeans(pgx$X[, grp0, drop = FALSE], na.rm = TRUE),
    x1 = rowMeans(pgx$X[, grp1, drop = FALSE], na.rm = TRUE)
  )

  ## the order the object happens to arrive in
  xy <- suppressMessages(playbase::pgx.contrastScatter(pgx, comp, data = TRUE))
  expect_equal(as.numeric(xy), as.numeric(expected))

  ## and a design whose rows are X's columns REVERSED, which is what the old
  ## positional read could not survive: same samples, different order, and the
  ## drawn values must not move
  flip <- pgx
  rev.order <- rev(rownames(pgx$model.parameters$exp.matrix))
  flip$model.parameters$exp.matrix <-
    pgx$model.parameters$exp.matrix[rev.order, , drop = FALSE]
  xy.flip <- suppressMessages(playbase::pgx.contrastScatter(flip, comp, data = TRUE))
  expect_equal(as.numeric(xy.flip), as.numeric(expected))
})

## Review W10. The scatterPlot guard was weakened from "pos has one row per
## sample" to "pos names samples". This measures what the weakening costs: a
## `pos` that names one sample of many now draws, and what it draws is that
## sample's own value -- the guard's job was never to stop under-drawing, and
## every read of `pos` below it is by name.

test_that("pgx.scatterPlot draws a short pos correctly rather than refusing it", {
  pgx <- outlier_pgx()
  gene <- rownames(pgx$X)[1]
  one <- pgx$tsne2d[1, , drop = FALSE]
  expect_lt(nrow(one), nrow(pgx$samples))

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  pts <- scatter_points(pgx, pos = one, gene = gene)
  expect_setequal(rownames(pts), rownames(one))
  expect_equal(as.numeric(pts$value), as.numeric(pgx$X[gene, rownames(pts)]),
    tolerance = 1e-3 ## the hover value is rounded for display
  )

  ## a wider one is still drawn sample-for-sample, not shifted
  three <- pgx$tsne2d[2:4, , drop = FALSE]
  pts <- scatter_points(pgx, pos = three, gene = gene)
  expect_setequal(rownames(pts), rownames(three))
  expect_equal(as.numeric(pts$value), as.numeric(pgx$X[gene, rownames(pts)]),
    tolerance = 1e-3
  )
})
