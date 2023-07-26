#' Test for gx.markermap
#'
#'

#' Test for gx.PCAheatmap
#'
#'

#' Test for gx.PCAcomponents
#'
#'

#' Test for gx.imagemap
#'
#'


#' Test for gx.splitmap
test_that("gx.heatmap runs correctly with defaults", {
  set.seed(123)
  gx <- matrix(rnorm(250)^2, nrow=25, ncol=10) 
  rownames(gx) <- sample(LETTERS, 25)
  colnames(gx) <- sample(letters, 10)
  # Create a grouping variable

  playbase::gx.splitmap(
    gx, 
    scale = "row",
    cluster_rows = FALSE,
    show_key = FALSE,
    cluster_columns = FALSE
  )

  # Test that an object has been created
  p <- dev.capture()
  expect_true(is.null(p))

})


#' Test for gx.heatmap
test_that("gx.heatmap runs correctly with defaults", {
  set.seed(123)
  gx <- matrix(rnorm(250)^2, nrow=25, ncol=10) 
  rownames(gx) <- sample(LETTERS, 25)
  colnames(gx) <- sample(letters, 10)
  # Create a grouping variable

  clust <- playbase::gx.heatmap(gx)

  # Check the clusters as outputs
  expect_type(clust, "list")

})

#' Test for clustermap
#'
#'


#' Test for frozenmap
test_that("frozenmap returns matrix with selected dim", {
  # Generate test data
  mat <- matrix(rnorm(20), 5, 4)
  rownames(mat) <- replicate(5, paste(sample(LETTERS, 4, replace = TRUE), collapse = ""))
  mat

  # Run function
  z <- frozenmap(mat, m = 2, n = 2)

  # Test class
  expect_equal(class(z), c("matrix", "array"))

  # Test dim
  expect_equal(dim(z), c(2, 2))
})

#' Test for multi.dist

test_that("multi.dist runs correctly", {
  data <- c(
    1, 2, 3, 4, 5,
    5, 4, 3, 2, 1,
    5, 5, 5, 5, 5,
    2, 2, 2, 2, 2
  )
  mat <- matrix(data, 5, 4)
  rownames(mat) <- replicate(5, paste(sample(LETTERS, 4, replace = TRUE), collapse = ""))
  D <- multi.dist(mat)
  D2 <- multi.dist(mat, p = 20)
  # Expected result
  expected <- c(0.3061224, 0.6326531, 0.8775510, 1.0000000, 0.3877551, 0.5918367, 0.8775510, 0.3877551, 0.6326531, 0.3061224)

  # Test object type
  expect_s3_class(D, "dist")

  # Test that rownames pased as labels
  expect_equal(attr(D, "Labels"), rownames(mat))

  # Test that the distance matrix is correct
  expect_equal(as.vector(D), expected, tolerance = 1e-3)
})
