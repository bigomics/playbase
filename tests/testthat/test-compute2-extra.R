#' Test for compute_cellcycle_gender
test_that("compute_cellcycle_gender adds cell cycle and gender data", {
  
  # Create mock data
  pgx <- list(
    samples = playbase::SAMPLES,
    counts = as.matrix(playbase::COUNTS),
    genes = playbase::GENES_TABLE,
    organism = "Human"
  )
  pgx$counts <- pgx$counts[!duplicated(rownames(pgx$counts)), , drop = FALSE]
  pgx$genes <- pgx$genes[!pgx$genes$symbol == "", ,drop = FALSE]
  pgx$counts <- pgx$counts[rownames(pgx$counts) %in% pgx$genes$symbol, , drop = FALSE]

  result <- playbase::compute_cellcycle_gender(pgx, pgx$counts)

  # Check cell cycle and gender added to sample df
  expect_true(".cell_cycle" %in% colnames(result$samples))
  expect_true(".gender" %in% colnames(result$samples))
})

test_that("compute_cellcycle_gender adds cell cycle and gender data", {
  
  # Create mock data
  pgx <- list(
    samples = playbase::SAMPLES,
    counts = as.matrix(playbase::COUNTS),
    genes = playbase::GENES_TABLE,
    organism = "Human"
  )
  pgx$counts <- pgx$counts[!duplicated(rownames(pgx$counts)), , drop = FALSE]
  pgx$genes <- pgx$genes[!pgx$genes$symbol == "", ,drop = FALSE]
  pgx$counts <- pgx$counts[rownames(pgx$counts) %in% pgx$genes$symbol, , drop = FALSE]
  
  result <- playbase::compute_cellcycle_gender(pgx)
  expect_cell_cycle_stages <- c("G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", 
  "G1", "S", "S", "S", "S", "S", "S", "S")
  # Check cell cycle and gender added to sample df
  expect_equal(result$samples$.cell_cycle, expect_cell_cycle_stages)
})
