#' Test for pgx.getFamilies
test_that("pgx.getFamilies returns expected families", {
  # Create mock data
  pgx <- list(
    samples = playbase::SAMPLES,
    counts = as.matrix(playbase::COUNTS),
    genes = playbase::GENES_TABLE,
    organism = "Human"
  )
  pgx$counts <- pgx$counts[!duplicated(rownames(pgx$counts)), , drop = FALSE]
  pgx$genes <- pgx$genes[!pgx$genes$symbol == "", , drop = FALSE]
  pgx$counts <- pgx$counts[rownames(pgx$counts) %in% pgx$genes$symbol, , drop = FALSE]

  # Test with default parameters
  families <- playbase::pgx.getFamilies(pgx)

  # Expect 65 families returned
  expect_equal(length(families), 65)

  # Expect <all> family returned
  expect_equal(families[1], "<all>")
  expect_equal(families[65], "FAMILY:Zinc fingers C2H2-type (HGNC)")
})
