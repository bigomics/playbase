## Covers detectOutlierSamples() method selection and the isoforest path,
## which pgx.preprocess() enables via opt$outlier_methods.

mk_X <- function(n) {
  set.seed(1)
  X <- matrix(rnorm(500 * n), 500, n,
    dimnames = list(paste0("g", 1:500), paste0("s", 1:n)))
  X[, n] <- X[, n] + 4 ## last sample is the outlier
  X
}

test_that("detectOutlierSamples returns one Z column per requested method", {
  X <- mk_X(8)
  res <- playbase::detectOutlierSamples(X, plot = FALSE)
  expect_equal(colnames(res$Z), c("z.correlation", "z.distance", "z.features"))
  expect_equal(unname(which.max(res$z.outlier)), 8L)

  ## a single method must stay a matrix, else rowMeans() fails
  res1 <- playbase::detectOutlierSamples(X, methods = "z.distance", plot = FALSE)
  expect_true(is.matrix(res1$Z))
  expect_equal(colnames(res1$Z), "z.distance")
  expect_length(res1$z.outlier, ncol(X))
})

test_that("detectOutlierSamples rejects an unknown method", {
  ## used to return NaN z-scores silently, reading as "no outliers"
  expect_error(
    playbase::detectOutlierSamples(mk_X(8), methods = "z.correlaton", plot = FALSE),
    "should be one of"
  )
})

test_that("isoforest is off by default and scores the outlier when asked", {
  skip_if_not_installed("isotree")
  X <- mk_X(8)
  expect_false("z.isoforest" %in% colnames(playbase::detectOutlierSamples(X, plot = FALSE)$Z))

  set.seed(42)
  res <- playbase::detectOutlierSamples(X, methods = "z.isoforest", plot = FALSE)
  expect_equal(colnames(res$Z), "z.isoforest")
  expect_equal(unname(which.max(res$z.outlier)), 8L)

  ## used to abort in irlba(nv = 3) with fewer than 4 samples
  small <- playbase::detectOutlierSamples(mk_X(3), methods = "z.isoforest", plot = FALSE)
  expect_length(small$z.outlier, 3L)
})

test_that("pgx.preprocess plumbs outlier_methods to detectOutlierSamples", {
  X <- mk_X(8)
  counts <- 2^X
  samples <- data.frame(group = rep(c("a", "b"), each = 4), row.names = colnames(X))
  run <- function(...) playbase::pgx.preprocess(counts, samples, contrasts = NULL,
    options = list(remove_outliers = TRUE, outlier_threshold = 3, ...))

  ## the option must reach detectOutlierSamples()' validation, not be ignored
  expect_error(run(outlier_methods = "nonsense"), "should be one of")
  expect_no_error(run(outlier_methods = "z.distance"))
})
