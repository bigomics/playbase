
#' Test for pgx.clusterSamples2
test_that("pgx.clusterSamples2 performs dimensionality reduction", {

  d <- get_mini_example_data()

  now = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  ngs <- list(
      name = "data set",
      this.date = now,
      date = now,
      datatype = "unknown",
      description = "data set",
      samples = data.frame(d$samples, check.names = FALSE),
      counts = as.matrix(d$counts),
      contrasts = d$contrasts,
      X = NULL,
      total_counts = Matrix::colSums(d$counts, na.rm = TRUE),
      counts_multiplier = 1
    )
  
  clustering_tests <- c("pca", "tsne", "umap")
  ngs <- playbase::pgx.clusterSamples2(
        ngs,
        dims = c(2),
        perplexity = NULL,
        methods = clustering_tests
      )

  # Check clustering performed
  expect_true("cluster" %in% names(ngs))

  # Expect number of test performed
  expect_equal(length(ngs$cluster$pos), length(clustering_tests))
  
  # Expect dimensions
  expect_equal(ncol(ngs$cluster$pos$pca2d), 2)
  expect_equal(ncol(ngs$cluster$pos$tsne2d), 2)
  expect_equal(ncol(ngs$cluster$pos$umap2d), 2)
  
})

#' Test for pgx.findLouvainClusters
test_that("pgx.findLouvainClusters runs correctly", {

  d <- get_mini_example_data()

  now = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  ngs <- list(
      name = "data set",
      this.date = now,
      date = now,
      datatype = "unknown",
      description = "data set",
      samples = data.frame(d$samples, check.names = FALSE),
      counts = as.matrix(d$counts),
      contrasts = d$contrasts,
      X = NULL,
      total_counts = Matrix::colSums(d$counts, na.rm = TRUE),
      counts_multiplier = 1
    )
  
  clustering_tests <- c("pca", "tsne", "umap")
  ngs <- playbase::pgx.clusterSamples2(
        ngs,
        dims = c(2),
        perplexity = NULL,
        methods = clustering_tests
      )
  posx <- scale(cbind(ngs$cluster$pos[["umap2d"]], ngs$cluster$pos[["tsne2d"]]))
  idx <- playbase::pgx.findLouvainClusters(posx, level = 1, prefix = "c", small.zero = 0.0)

  # Expect number of test performed
  expect_equal(length(ngs$cluster$pos), length(clustering_tests))
  
  # Expect dimensions
  expect_equal(length(idx), nrow(posx))
  expect_type(idx, "character")
  
})
