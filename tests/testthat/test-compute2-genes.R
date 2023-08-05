#' Test for compute_testGenes
pgx_data <- playbase:::get_mini_example_data()

pgx <- playbase::pgx.createPGX(
  samples = pgx_data$samples,
  counts = pgx_data$counts,
  contrasts = pgx_data$contrast
)
contr.matrix <- pgx$contrasts
contr.matrix <- makeContrastsFromLabelMatrix(contr.matrix)
contr.matrix <- sign(contr.matrix) ## sign is fine

sel <- Matrix::colSums(contr.matrix == -1) > 0 & Matrix::colSums(contr.matrix == 1) > 0
contr.matrix <- contr.matrix[, sel, drop = FALSE]


test_that("compute_testGenes returns correct number of genes", {
  test_genes <- compute_testGenes(
    pgx, contr.matrix,
    test.methods = c("ttest.welch", "trend.limma", "edger.qlf")
  )

  # Check that all the test have ran
  expect_equal(length(test_genes$gx.meta$meta), 2)
  expect_equal(length(test_genes$gx.meta$sig.counts), 4)
  expect_equal(dim(test_genes$gx.meta$meta$act12h_vs_notact)[2], 8)
  expect_equal(dim(test_genes$gx.meta$meta$act24h_vs_notact)[2], 8)

  # Check static value of in tables
  sum_vals <- sapply(test_genes$gx.meta$meta, sum)
  expected_vals <- c(1929.020, 2137.772)
  names(expected_vals) <- names(test_genes$gx.meta$meta)
  expect_equal(sum_vals, expected_vals, tolerance = 0.1)
})

#' Test for compute_testGenesSingleOmics
test_that("compute_testGenesSingleOmics runs without errors", {
  result <- compute_testGenesSingleOmics(pgx, contr.matrix)
  expected_slots <- c("meta", "sig.counts")

  expect_equal(names(result$gx.meta), expected_slots)
  # Check that all the test have ran
  expect_equal(length(result$gx.meta$meta), 2)
  expect_equal(length(result$gx.meta$sig.counts), 4)
  expect_equal(dim(result$gx.meta$meta$act12h_vs_notact)[2], 8)
  expect_equal(dim(result$gx.meta$meta$act24h_vs_notact)[2], 8)

  # Check static value of in tables
  sum_vals <- sapply(result$gx.meta$meta, sum)
  expected_vals <- c(1929.020, 2137.772)
  names(expected_vals) <- names(result$gx.meta$meta)
  expect_equal(sum_vals, expected_vals, tolerance = 0.1)
})
