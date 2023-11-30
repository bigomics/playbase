# Generate sample data
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

clustering_tests <- c("pca", "tsne", "umap")
ngs <- playbase::pgx.clusterSamples2(
  ngs,
  dims = c(2),
  perplexity = NULL,
  methods = clustering_tests
)
posx <- scale(cbind(ngs$cluster$pos[["umap2d"]], ngs$cluster$pos[["tsne2d"]]))
idx <- playbase::pgx.findLouvainClusters(posx, level = 1, prefix = "c", small.zero = 0.0)



#' Test for makeDirectContrasts
ngs$samples$cluster <- idx
Y <- ngs$samples[, "cluster", drop = FALSE]
ct <- playbase::makeDirectContrasts(Y, ref = "others")
test_that("makeDirectContrasts runs correctly", {
  # Check makeDirectContrasts performed
  expect_equal(names(ct), c("contr.matrix", "group", "exp.matrix"))
  expect_true(all(ct$contr.matrix %in% c(1, -1)))
  expect_true(all(ct$exp.matrix %in% c(1, -1)))
})

#' Test for contrastAsLabels
# Run contrastAsLabels on ct
ctx <- playbase::contrastAsLabels(ct$exp.matrix)

test_that("contrastAsLabels runs on ct", {
  # Check makeDirectContrasts performed
  expected <- rep("character", ncol(ctx))
  names(expected) <- colnames(ctx)
  expect_equal(apply(ctx, 2, class), expected)
  # validate first column for ctx
  string <- colnames(ctx)[1]
  string <- gsub("^[^:]*:", "", string)

  # Generate the expected labels
  char_labels <- strsplit(string, "_vs_")[[1]]
  expected <- ifelse(ct$exp.matrix[, 1] == 1, char_labels[1], char_labels[2])
  names(expected) <- NULL
  expect_equal(ctx[, 1], expected)
})

# Run contrastAsLabels on contrast
test_that("contrastAsLabels runs on input contrast", {
  contrast <- playbase::CONTRASTS

  # Run function
  result <- playbase::contrastAsLabels(contrast)

  # Check class
  expect_s3_class(result, "data.frame")

  # Check dimensions
  expect_equal(dim(result), dim(contrast))

  # Check preserve rownames
  expect_true(all(rownames(contrast) == rownames(result)))
  expect_true(all(colnames(contrast) == colnames(result)))

  # Check 0 become NAs
  expect_equal(sum(contrast == 0), sum(is.na(result)))
})
