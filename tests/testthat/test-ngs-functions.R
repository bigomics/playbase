#' Test for ngs.getGeneAnnotation
test_that("ngs.getGeneAnnotation returns annotation for genes", {
  d <- get_mini_example_data()
  genes <- sample(rownames(d$counts), 10)

  result <- playbase::ngs.getGeneAnnotation(genes)

  # Check class
  expect_s3_class(result, "data.frame")

  # Check number of rows
  expect_equal(nrow(result), length(genes))

  # Check presence of expected columns
  expected_cols <- c("gene_name", "gene_title", "gene_biotype", "chr", "pos", "tx_len", "map")
  expect_true(all(expected_cols %in% colnames(result)))

  # Check gene names match
  expect_equal(result$gene_name, genes)
})
