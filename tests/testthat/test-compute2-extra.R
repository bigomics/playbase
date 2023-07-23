
#' Test for compute_cellcycle_gender
test_that("compute_cellcycle_gender adds cell cycle and gender data", {

  now = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  ngs <- list(
    name = "data set",
    this.date = now,
    date = now,
    datatype = "unknown",
    description = "data set",
    samples = data.frame(playbase::SAMPLES, check.names = FALSE),
    counts = as.matrix(playbase::COUNTS),
    contrasts = playbase::CONTRASTS,
    X = NULL,
    total_counts = Matrix::colSums(d$counts, na.rm = TRUE),
    counts_multiplier = 1
  )
  ngs$genes <- playbase::ngs.getGeneAnnotation(genes = rownames(ngs$counts))
  ngs <- playbase::compute_cellcycle_gender(ngs)
  
  result <- playbase::compute_cellcycle_gender(ngs)

  # Check cell cycle and gender added to sample df
  expect_true(".cell_cycle" %in% colnames(result$samples))
  expect_true(".gender" %in% colnames(result$samples))
  
})
