#' Test for imputeMedian
test_that("imputeMedian replaces NA with median", {
  # Test matrix input
  X <- matrix(c(1, NA, 3, NA, 2, NA), ncol = 2)
  expected <- matrix(c(1, 2, 3, 1, 2, 3), ncol = 2)

  imputed <- playbase::imputeMedian(X)

  # Compare result to expected
  expect_equal(imputed, expected)
  expect_true(sum(is.na(imputed)) == 0)
})

#' Test for averageByGroup
mat <- matrix(rep(1:10, 10), ncol = 10, byrow = TRUE)
groups <- gl(2, 5)

test_that("averageByGroup returns a numeric matrix", {
  means <- playbase::averageByGroup(mat, groups)

  expect_type(means, "double")
  expect_equal(dim(means), c(10, 2))
})

test_that("averageByGroup calculates correct group means", {
  means <- playbase::averageByGroup(mat, groups)
  expected_means <- matrix(c(
    rep(3, 10),
    rep(8, 10)
  ), ncol = 2)
  colnames(expected_means) <- c(1, 2)
  expect_equal(means, expected_means)
})

#' Test for gmean
test_that("gmean works as expected", {
  # Test for running as expected
  expect_equal(playbase::gmean(c(1, 2, 3)), 1.8171205928321396)

  # Test for handling 0
  expect_equal(playbase::gmean(c(0, 1, 2)), 0)

  # Test for handling NA
  expect_equal( round(playbase::gmean(c(NA, 1, 2)),4), 1.4142)
})


#' Test for mat2hugo
#'
#'

#' Test for gx.hist
#'
#'

#' Test for val2col
test_that("val2col returns colour vector of lenght n", {
  # Generate Data
  set.seed(10)
  z <- rnorm(100)
  cols <- playbase::val2col(z, zlim = c(-3, 3))

  expect_type(cols, "character")
  expect_equal(length(cols), length(z))
})
