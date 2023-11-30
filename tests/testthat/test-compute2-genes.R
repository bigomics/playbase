#' Test for compute_testGenes
test_that("compute_testGenes returns correct number of genes", {
  
  # Call mock data 
  pgx <- playbase::PGX_CREATE
  ## make proper contrast matrix
  contr.matrix <- pgx$contrasts
  contr.values <- unique(as.vector(contr.matrix))
  is.numcontrast <- all(contr.values %in% c(NA, -1, 0, 1))
  is.numcontrast <- is.numcontrast && (-1 %in% contr.values) && (1 %in% contr.values)
  if (!is.numcontrast) {
    contr.matrix <- playbase::makeContrastsFromLabelMatrix(contr.matrix)
    contr.matrix <- sign(contr.matrix) ## sign is fine
  }

  ## select valid contrasts
  sel <- Matrix::colSums(contr.matrix == -1) > 0 & Matrix::colSums(contr.matrix == 1) > 0
  contr.matrix <- contr.matrix[, sel, drop = FALSE]

  suppressWarnings(
  test_genes <- playbase::compute_testGenes(
    pgx, contr.matrix,
    test.methods = c("ttest.welch", "trend.limma", "edger.qlf")
  ))

  # Check that all the test have ran
  expect_equal(length(test_genes$gx.meta$meta), 2)
  expect_equal(length(test_genes$gx.meta$sig.counts), 4)
  expect_equal(dim(test_genes$gx.meta$meta$act12h_vs_notact)[2], 8)
  expect_equal(dim(test_genes$gx.meta$meta$act24h_vs_notact)[2], 8)

  # Check static value of in tables
  sum_vals <- sapply(test_genes$gx.meta$meta, sum)
  expected_vals <- c(8333, 8523)
  names(expected_vals) <- names(test_genes$gx.meta$meta)
  expect_equal(sum_vals, expected_vals, tolerance = 0.1)
})

#' Test for compute_testGenesSingleOmics
test_that("compute_testGenesSingleOmics runs without errors", {
  # Call mock data 
  pgx <- playbase::PGX_CREATE
  ## make proper contrast matrix
  contr.matrix <- pgx$contrasts
  contr.values <- unique(as.vector(contr.matrix))
  is.numcontrast <- all(contr.values %in% c(NA, -1, 0, 1))
  is.numcontrast <- is.numcontrast && (-1 %in% contr.values) && (1 %in% contr.values)
  if (!is.numcontrast) {
    contr.matrix <- playbase::makeContrastsFromLabelMatrix(contr.matrix)
    contr.matrix <- sign(contr.matrix) ## sign is fine
  }

  ## select valid contrasts
  sel <- Matrix::colSums(contr.matrix == -1) > 0 & Matrix::colSums(contr.matrix == 1) > 0
  contr.matrix <- contr.matrix[, sel, drop = FALSE]

  suppressWarnings(
    result <- playbase::compute_testGenesSingleOmics(pgx, contr.matrix)
  )
  expected_slots <- c("meta", "sig.counts")

  expect_equal(names(result$gx.meta), expected_slots)
  # Check that all the test have ran
  expect_equal(length(result$gx.meta$meta), 2)
  expect_equal(length(result$gx.meta$sig.counts), 4)
  expect_equal(dim(result$gx.meta$meta$act12h_vs_notact)[2], 8)
  expect_equal(dim(result$gx.meta$meta$act24h_vs_notact)[2], 8)

  # Check static value of in tables
  sum_vals <- sapply(result$gx.meta$meta, sum)
  expected_vals <- c(8244, 7820)
  names(expected_vals) <- names(result$gx.meta$meta)
  expect_equal(sum_vals, expected_vals, tolerance = 0.1)
})
