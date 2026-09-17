## Playbase integration with the direct preprocessing leaf.
##
## These tests cover metadata reduction and the explicit final alignment.
## Numerical family behavior is tested in playbase.preprocess itself.

make_preprocess_fixture <- function() {
  counts <- matrix(
    c(
      10,
      11,
      30,
      31,
      0,
      0,
      8,
      9,
      5,
      6,
      7,
      8
    ),
    nrow = 3L,
    byrow = TRUE,
    dimnames = list(c("a", "b", "c"), paste0("s", 1:4))
  )
  samples <- data.frame(
    group = c("control", "control", "case", "case"),
    batch = c("one", "two", "one", "two"),
    row.names = colnames(counts)
  )
  contrasts <- matrix(
    samples$group,
    ncol = 1L,
    dimnames = list(colnames(counts), "case_vs_control")
  )
  list(counts = counts, samples = samples, contrasts = contrasts)
}

test_that("leaf reduces the contrast design without changing source counts", {
  fixture <- make_preprocess_fixture()
  result <- playbase.preprocess::pgx.preprocess(
    counts = fixture$counts,
    samples = fixture$samples,
    contrasts = fixture$contrasts,
    options = list(
      normalize = FALSE,
      filter_missing = FALSE,
      remove_outliers = FALSE
    )
  )

  expect_identical(result$counts, fixture$counts)
  expect_named(
    result,
    c("counts", "X", "annot", "prior", "space", "alignment", "options")
  )
  expect_identical(result$alignment$cols, seq_len(ncol(fixture$counts)))
})

test_that("leaf derives plain target and batch metadata", {
  fixture <- make_preprocess_fixture()
  options <- list(
    normalize = FALSE,
    batch.correct.method = "limma"
  )
  result <- playbase.preprocess::pgx.preprocess(
    counts = fixture$counts,
    samples = fixture$samples,
    contrasts = fixture$contrasts,
    options = options
  )

  expect_identical(
    unname(result$options$target),
    unname(apply(fixture$contrasts, 1L, paste, collapse = "_"))
  )
  expect_identical(
    unname(result$options$batch$batch),
    fixture$samples$batch
  )
})

test_that("unknown legacy options fail at the matrix boundary", {
  fixture <- make_preprocess_fixture()
  expect_error(
    playbase.preprocess::pgx.preprocess(
      counts = fixture$counts,
      samples = fixture$samples,
      contrasts = fixture$contrasts,
      options = list(datatype = "RNA-seq")
    ),
    "unknown option"
  )
})

test_that("private count-scale policy supplies explicit leaf inputs", {
  fixture <- make_preprocess_fixture()
  result <- playbase.preprocess::pgx.preprocess(
    counts = fixture$counts,
    samples = fixture$samples,
    contrasts = fixture$contrasts,
    options = list(normalize = FALSE)
  )
  pgx <- list(
    counts = result$counts,
    X = result$X,
    settings = list(preprocess = .pgx_preprocess_metadata(result))
  )
  count_scale <- .pgx_count_scale_matrix(pgx)

  expect_identical(dimnames(count_scale), dimnames(result$X))
  expect_equal(colSums(count_scale), colSums(fixture$counts))
})

test_that("configurable CPM defines empty dense and sparse libraries", {
  counts <- matrix(
    c(0, 0, 1, 2),
    nrow = 2L,
    dimnames = list(c("g1", "g2"), c("empty", "observed"))
  )
  dense <- .pgx_log_cpm(counts, total = 1e6, prior = 1)
  sparse <- .pgx_log_cpm(
    methods::as(Matrix::Matrix(counts, sparse = TRUE), "generalMatrix"),
    total = 1e6,
    prior = 1
  )

  expect_true(all(is.finite(dense)))
  expect_identical(dense[, "empty"], c(g1 = 0, g2 = 0))
  expect_equal(as.matrix(sparse), dense)
  expect_identical(dimnames(sparse), dimnames(counts))
})

test_that("mixed count-scale policy preserves each processed count-like layer", {
  counts <- matrix(
    c(
      10,
      20,
      30,
      20,
      5,
      15,
      25,
      15,
      8,
      12,
      12,
      8,
      0.2,
      0.8,
      0.4,
      0.6
    ),
    nrow = 8L,
    byrow = TRUE,
    dimnames = list(
      c(
        "log:a",
        "log:b",
        "linear:a",
        "linear:b",
        "counts:a",
        "counts:b",
        "meth:a",
        "meth:b"
      ),
      c("s1", "s2")
    )
  )
  processed <- counts
  processed[1:2, ] <- log2(matrix(c(1, 9, 9, 1), 2L, byrow = TRUE) + 1)
  processed[3:4, ] <- matrix(c(2, 8, 8, 2), 2L, byrow = TRUE)
  processed[5:6, ] <- matrix(c(3, 7, 7, 3), 2L, byrow = TRUE)
  processed[7:8, ] <- matrix(c(0.9, 0.1, 0.1, 0.9), 2L, byrow = TRUE)
  dimnames(processed) <- dimnames(counts)
  layers <- rep(c("log", "linear", "counts", "meth"), each = 2L)
  metadata <- list(
    alignment = list(
      rows = lapply(seq_len(nrow(counts)), as.integer),
      cols = seq_len(ncol(counts))
    ),
    space = c(
      log = "log2",
      linear = "linear",
      counts = "counts",
      meth = "beta"
    ),
    prior = c(log = 1, linear = NA, counts = NA, meth = NA),
    layers = layers,
    options = NULL
  )
  pgx <- list(
    counts = counts,
    X = processed,
    settings = list(preprocess = metadata)
  )

  actual <- .pgx_count_scale_matrix(pgx)
  aligned <- playbase.preprocess::pp.alignCounts(
    counts,
    metadata$alignment,
    X = processed
  )
  expected <- aligned
  for (layer in c("log", "linear", "counts")) {
    rows <- which(layers == layer)
    expected[rows, ] <- playbase.preprocess::pp.countScaleMatrix(
      X = processed[rows, , drop = FALSE],
      counts = counts,
      alignment = list(
        rows = metadata$alignment$rows[rows],
        cols = metadata$alignment$cols
      ),
      space = if (layer == "linear") "counts" else metadata$space[[layer]],
      prior = metadata$prior[[layer]]
    )
  }

  expect_identical(dimnames(actual), dimnames(processed))
  expect_equal(actual, expected)
  expect_false(isTRUE(all.equal(actual[1:6, ], aligned[1:6, ])))
  expect_identical(actual[7:8, ], aligned[7:8, ])
})
