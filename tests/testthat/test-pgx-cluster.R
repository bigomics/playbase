# Generate data for all functions
d <- playbase::get_mini_example_data()

now <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
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

#' Test for pgx.clusterSamples2
#' Run clusters
clustering_tests <- c("pca", "tsne", "umap")
ngs <- playbase::pgx.clusterSamples(
  ngs,
  dims = c(2),
  perplexity = NULL,
  methods = clustering_tests
)

test_that("pgx.clusterSamples2 performs dimensionality reduction", {
  # Check clustering performed
  expect_true("cluster" %in% names(ngs))

  # Expect number of test performed
  expect_equal(length(ngs$cluster$pos), length(clustering_tests) + 1)

  # Expect dimensions
  expected_clusts <- paste0(clustering_tests, 2, "d")
  expect_equal(ncol(ngs$cluster$pos[[expected_clusts[1]]]), 2)
  expect_equal(ncol(ngs$cluster$pos[[expected_clusts[2]]]), 2)
  expect_equal(ncol(ngs$cluster$pos[[expected_clusts[3]]]), 2)
})

#' Test for pgx.findLouvainClusters
#' Run Louvain clusters
posx <- scale(cbind(ngs$cluster$pos[["umap2d"]], ngs$cluster$pos[["tsne2d"]]))
idx <- playbase::pgx.findLouvainClusters(posx, level = 1, prefix = "c", small.zero = 0.0)

test_that("pgx.findLouvainClusters runs correctly", {
  # Expect number of test performed
  expect_equal(length(ngs$cluster$pos), length(clustering_tests) + 1)

  # Expect dimensions
  expect_equal(length(idx), nrow(posx))
  expect_type(idx, "character")
})

#' Test for pgx.clusterGenes
# Run logCPM normalisation and clustering for genes
ngs$X <- playbase::logCPM(ngs$counts, total = 1e6, prior = 1)
# Run pgx.clusterGenes for only one clustering method (umap)
ngs <- playbase::pgx.clusterGenes(ngs, methods = "umap", dims = c(3), level = "gene")

test_that("pgx.findLouvainClusters runs correctly", {
  # Expect number of test performed
  expect_true("cluster.genes" %in% names(ngs))
  expect_equal(length(ngs$cluster.genes), 1)

  # Expect dimensions
  expected_gene_clusts <- paste0("umap", 2, "d")
  expect_equal(ncol(ngs$cluster.genes$pos[[expected_gene_clusts[1]]]), 2)
})
