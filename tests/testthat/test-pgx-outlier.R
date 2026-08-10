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

test_that("isoforest is off by default and runs when asked", {
  skip_if_not_installed("isotree")
  X <- mk_X(8)
  expect_false("z.isoforest" %in% colnames(playbase::detectOutlierSamples(X, plot = FALSE)$Z))

  ## NB: structure only. The isoforest's score scale and sign are under review
  ## (it is two-sided and sd-normalized, unlike the mad-based methods), so this
  ## deliberately does not assert which sample scores highest.
  res <- playbase::detectOutlierSamples(X, methods = "z.isoforest", plot = FALSE)
  expect_equal(colnames(res$Z), "z.isoforest")
  expect_true(all(is.finite(res$z.outlier)))

  ## used to abort in irlba(nv = 3) with fewer than 4 samples
  small <- playbase::detectOutlierSamples(mk_X(3), methods = "z.isoforest", plot = FALSE)
  expect_length(small$z.outlier, 3L)
})

test_that("the mad-based methods score the injected outlier well above threshold", {
  res <- playbase::detectOutlierSamples(mk_X(8), plot = FALSE)
  expect_gt(res$z.outlier[["s8"]], 3)
  expect_true(all(res$z.outlier[paste0("s", 1:7)] < 3))
})

test_that("tied samples give finite scores instead of NaN", {
  ## mad() is 0 when over half the values tie, which used to yield 0/0 for
  ## every method and crash pgx.preprocess with "missing value where
  ## TRUE/FALSE needed". Same sample uploaded 3x, twice.
  set.seed(1)
  X <- matrix(rnorm(500 * 6), 500, 6,
    dimnames = list(paste0("g", 1:500), paste0("s", 1:6)))
  X[, 2] <- X[, 1]; X[, 3] <- X[, 1]
  X[, 5] <- X[, 4]; X[, 6] <- X[, 4]

  res <- playbase::detectOutlierSamples(X, plot = FALSE)
  expect_true(all(is.finite(res$Z)))
  expect_true(all(is.finite(res$z.outlier)))

  samples <- data.frame(group = rep(c("a", "b"), each = 3), row.names = colnames(X))
  out <- playbase::pgx.preprocess(2^X, samples, contrasts = NULL,
    options = list(remove_outliers = TRUE, outlier_threshold = 3))
  expect_equal(ncol(out$X), 6L) ## identical samples, none is an outlier
})

test_that("plotting works on degenerate input and restores par", {
  pdf(NULL)
  on.exit({ dev.off(); unlink("Rplots.pdf") }, add = TRUE)
  before <- graphics::par("mfrow")

  ## non-finite Z used to abort with "need finite 'ylim' values"
  res <- playbase::detectOutlierSamples(mk_X(3), plot = TRUE)
  expect_length(res$z.outlier, 3L)
  expect_equal(graphics::par("mfrow"), before)

  ## all four methods, 6 panels
  skip_if_not_installed("isotree")
  expect_length(playbase::detectOutlierSamples(mk_X(8), methods = NULL, plot = TRUE)$z.outlier, 8L)
  expect_equal(graphics::par("mfrow"), before)
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
