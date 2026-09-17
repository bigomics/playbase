#' tiledb.prepareData() hands tiledb.writeData() two matrices that it walks with
#' ONE index `j`, so they have to be aligned on both axes. Rows were aligned by
#' gene id already; samples were not, and under D-24 outlier removal leaves `X`
#' with fewer samples than `counts`. Unaligned, every sample after the first
#' dropped one is written with ANOTHER sample's z-scores, and the last ones run
#' off the end of the matrix.
test_that("tiledb.prepareData aligns z-scores to counts on both axes", {
  counts <- matrix(1:24 * 10,
    nrow = 4, dimnames = list(paste0("G", 1:4), paste0("S", 1:6))
  )
  ## X is what preprocessing leaves: G4 removed, S3 and S6 removed as outliers.
  X <- log2(counts[1:3, c("S1", "S2", "S4", "S5")] + 1)
  genes <- data.frame(
    human_ortholog = paste0("H", 1:4), symbol = paste0("H", 1:4),
    row.names = rownames(counts)
  )

  d <- playbase:::tiledb.prepareData(list(counts = counts, X = X, genes = genes))

  ## the two matrices writeData() indexes together are one shape
  expect_identical(dim(d$zscores), dim(d$counts))
  expect_identical(length(d$genes), nrow(d$counts))
  expect_identical(length(d$samples), ncol(d$counts))
  expect_identical(d$samples, colnames(counts))

  ## the samples X does not have carry no z-score rather than someone else's
  expect_true(all(is.na(d$zscores[, c(3, 6)])))
  expect_false(anyNA(d$zscores[1:3, c(1, 2, 4, 5)]))

  ## and the ones it does have carry their own
  z <- (X - rowMeans(X)) / apply(X, 1, stats::sd)
  expect_equal(unname(d$zscores[1:3, c(1, 2, 4, 5)]), unname(z))

  ## the feature X does not have is dropped with its counts row, as before
  expect_identical(d$genes, paste0("H", 1:4))
  expect_true(all(is.na(d$zscores[4, ])))
})
