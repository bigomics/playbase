## PGX lifecycle around the direct preprocessing result.
##
## Source counts remain literal while X and explicit alignment carry analysis
## filtering. Tests use annotation-free settings to isolate that contract.

make_pgx_fixture <- function() {
  set.seed(21)
  counts <- matrix(
    stats::rpois(80 * 8, 50),
    nrow = 80,
    dimnames = list(paste0("G", 1:80), paste0("S", 1:8))
  )
  counts[1:5, ] <- NA
  samples <- data.frame(
    group = rep(c("a", "b"), each = 4),
    row.names = colnames(counts)
  )
  contrasts <- matrix(
    samples$group,
    ncol = 1L,
    dimnames = list(colnames(counts), "b_vs_a")
  )
  list(counts = counts, samples = samples, contrasts = contrasts)
}

create_test_pgx <- function(
  fixture,
  options = list(),
  filter.genes = FALSE,
  prune.samples = FALSE,
  batch.correct.method = "no_batch_correct"
) {
  defaults <- list(
    normalize = FALSE,
    filter_missing = TRUE,
    filter_threshold = 3,
    remove_outliers = FALSE
  )
  suppressMessages(playbase::pgx.createPGX(
    counts = fixture$counts,
    samples = fixture$samples,
    contrasts = fixture$contrasts,
    organism = "No organism",
    datatype = "RNA-seq",
    preprocess = utils::modifyList(defaults, options),
    batch.correct.method = batch.correct.method,
    add.gmt = FALSE,
    filter.genes = filter.genes,
    only.known = FALSE,
    only.proteincoding = FALSE,
    convert.hugo = FALSE,
    prune.samples = prune.samples
  ))
}

test_that("preprocess NULL preserves the legacy log2 prior through the leaf", {
  counts <- matrix(
    c(2, 4, 8, 16, 3, 6, 12, 24),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(c("G1", "G2"), paste0("S", 1:4))
  )
  samples <- data.frame(
    group = c("a", "a", "b", "b"),
    row.names = colnames(counts)
  )
  contrasts <- matrix(
    samples$group,
    ncol = 1L,
    dimnames = list(colnames(counts), "b_vs_a")
  )

  pgx <- suppressMessages(playbase::pgx.createPGX(
    counts = counts,
    samples = samples,
    contrasts = contrasts,
    organism = "No organism",
    datatype = "RNA-seq",
    preprocess = NULL,
    norm_method = "raw",
    filter.genes = FALSE,
    only.known = FALSE,
    only.proteincoding = FALSE,
    convert.hugo = FALSE,
    add.gmt = FALSE
  ))

  expect_identical(pgx$counts, counts)
  expect_identical(pgx$X, log2(counts + min(counts)))
  expect_identical(pgx$settings$preprocess$prior, min(counts))
  expect_false(pgx$settings$options$normalize)
  expect_false(pgx$settings$options$filter_missing)
  expect_false(pgx$settings$options$impute)
  expect_false(pgx$settings$options$remove_outliers)
  expect_identical(
    pgx$settings$options$batch.correct.method,
    "no_batch_correct"
  )
})

test_that("preprocess NULL rejects counts without a finite positive value", {
  counts <- matrix(
    0,
    nrow = 2L,
    ncol = 4L,
    dimnames = list(c("G1", "G2"), paste0("S", 1:4))
  )
  samples <- data.frame(
    group = c("a", "a", "b", "b"),
    row.names = colnames(counts)
  )
  contrasts <- matrix(
    samples$group,
    ncol = 1L,
    dimnames = list(colnames(counts), "b_vs_a")
  )

  expect_error(
    playbase::pgx.createPGX(
      counts = counts,
      samples = samples,
      contrasts = contrasts,
      organism = "No organism",
      datatype = "RNA-seq",
      preprocess = NULL,
      add.gmt = FALSE
    ),
    "counts must contain a finite positive value to derive the legacy prior",
    fixed = TRUE
  )
})

test_that("back-transformed log2 uploads run the complete preprocessing path", {
  skip_if_not_installed("limma")
  set.seed(42)
  raw.counts <- matrix(
    stats::rpois(30 * 8, 40),
    nrow = 30L,
    dimnames = list(paste0("G", 1:30), paste0("S", 1:8))
  )
  raw.counts[1L, ] <- 0
  uploaded.log2 <- log2(raw.counts + 1)
  counts <- 2^uploaded.log2
  log.prior <- min(counts)
  counts <- counts - log.prior
  samples <- data.frame(
    group = rep(c("a", "b"), each = 4L),
    batch = rep(c("one", "two"), 4L),
    row.names = colnames(counts)
  )
  contrasts <- matrix(
    samples$group,
    ncol = 1L,
    dimnames = list(colnames(counts), "b_vs_a")
  )
  options <- list(
    normalize = TRUE,
    norm_method = "CPM",
    filter_missing = FALSE,
    impute = FALSE,
    remove_outliers = FALSE,
    batch = samples[, "batch", drop = FALSE],
    target = stats::setNames(samples$group, rownames(samples))
  )
  create <- function(method) {
    suppressMessages(playbase::pgx.createPGX(
      counts = counts,
      samples = samples,
      contrasts = contrasts,
      organism = "No organism",
      datatype = "RNA-seq",
      preprocess = options,
      batch.correct.method = method,
      settings = list(input.transform = list(
        source_space = "log2",
        log_prior = log.prior,
        operation = "2^X-min"
      )),
      filter.genes = FALSE,
      only.known = FALSE,
      only.proteincoding = FALSE,
      convert.hugo = FALSE,
      add.gmt = FALSE
    ))
  }

  uncorrected <- create("no_batch_correct")
  corrected <- create("limma")

  expect_equal(counts, raw.counts, tolerance = 1e-12)
  expect_equal(corrected$counts, raw.counts, tolerance = 1e-12)
  expect_true(corrected$settings$options$normalize)
  expect_identical(
    corrected$settings$options$batch.correct.method,
    "limma"
  )
  expect_false(identical(corrected$X, uncorrected$X))
  expect_identical(
    corrected$settings$input.transform,
    list(
      source_space = "log2",
      log_prior = log.prior,
      operation = "2^X-min"
    )
  )
})

test_that("caller-supplied X keeps identity preprocessing metadata", {
  fixture <- make_pgx_fixture()
  X <- log2(fixture$counts + 1)
  pgx <- suppressMessages(playbase::pgx.createPGX(
    counts = fixture$counts,
    X = X,
    samples = fixture$samples,
    contrasts = fixture$contrasts,
    organism = "No organism",
    datatype = "RNA-seq",
    preprocess = list(unknown_option = TRUE),
    filter.genes = FALSE,
    only.known = FALSE,
    only.proteincoding = FALSE,
    convert.hugo = FALSE,
    add.gmt = FALSE
  ))

  expect_identical(pgx$X, X)
  expect_identical(
    pgx$settings$preprocess$alignment$rows,
    lapply(seq_len(nrow(X)), as.integer)
  )
  expect_identical(
    pgx$settings$preprocess$alignment$cols,
    seq_len(ncol(X))
  )
  expect_null(pgx$settings$preprocess$options)
})

test_that("single-cell construction remains independent of the leaf", {
  local_mocked_bindings(
    pgx.createSingleCellPGX = function(...) list(route = "single-cell"),
    .package = "playbase"
  )
  counts <- matrix(
    1:4,
    nrow = 2L,
    dimnames = list(c("G1", "G2"), c("S1", "S2"))
  )
  samples <- data.frame(group = c("a", "b"), row.names = colnames(counts))

  pgx <- suppressMessages(playbase::pgx.createPGX(
    counts = counts,
    samples = samples,
    contrasts = NULL,
    organism = "Human",
    datatype = "scRNA-seq",
    preprocess = list(unknown_option = TRUE)
  ))

  expect_identical(pgx, list(route = "single-cell"))
})

test_that("pgx.createPGX preserves counts and persists final metadata", {
  fixture <- make_pgx_fixture()
  pgx <- create_test_pgx(fixture)

  expect_identical(pgx$counts, fixture$counts)
  expect_lt(nrow(pgx$X), nrow(pgx$counts))
  expect_named(
    pgx$settings$preprocess,
    c("alignment", "space", "prior", "layers", "options")
  )
  expect_identical(
    playbase.preprocess::pp.alignCounts(
      pgx$counts,
      pgx$settings$preprocess$alignment,
      X = pgx$X
    ) |>
      dimnames(),
    dimnames(pgx$X)
  )
})

test_that("count-scale composition follows processed X and explicit alignment", {
  fixture <- make_pgx_fixture()
  pgx <- create_test_pgx(fixture)
  count_scale <- .pgx_count_scale_matrix(pgx)

  expect_identical(dim(count_scale), dim(pgx$X))
  expect_identical(dimnames(count_scale), dimnames(pgx$X))
  aligned <- playbase.preprocess::pp.alignCounts(
    pgx$counts,
    pgx$settings$preprocess$alignment,
    X = pgx$X
  )
  expect_identical(is.na(count_scale), is.na(aligned))
})

test_that("lifecycle row filters preserve source counts", {
  fixture <- make_pgx_fixture()
  fixture$counts[10, ] <- 0
  pgx <- create_test_pgx(
    fixture,
    options = list(filter_missing = FALSE),
    filter.genes = TRUE
  )

  expect_identical(pgx$counts, fixture$counts)
  expect_false("G10" %in% rownames(pgx$X))
  expect_identical(
    dim(playbase.preprocess::pp.alignCounts(
      pgx$counts,
      pgx$settings$preprocess$alignment,
      X = pgx$X
    )),
    dim(pgx$X)
  )
})

test_that("createPGX applies selected batch correction inside preprocessing", {
  skip_if_not_installed("sva")
  skip_if_not_installed("limma")
  fixture <- make_pgx_fixture()
  fixture$counts[10L, ] <- 0
  batch <- data.frame(
    site = rep(c("one", "two"), 4L),
    row.names = colnames(fixture$counts)
  )
  target <- stats::setNames(fixture$samples$group, rownames(fixture$samples))
  options <- list(
    batch = batch,
    target = target,
    normalize = FALSE,
    filter_missing = TRUE,
    remove_outliers = FALSE
  )
  uncorrected <- create_test_pgx(fixture, options, filter.genes = FALSE)

  for (method in c("ComBat", "limma")) {
    batch.args <- if (identical(method, "limma")) {
      list(use_covariates = TRUE)
    } else {
      list()
    }
    corrected <- create_test_pgx(
      fixture,
      utils::modifyList(options, list(batch_args = batch.args)),
      filter.genes = TRUE,
      batch.correct.method = method
    )
    expected.target <- target[colnames(uncorrected$X)]
    expected.batch <- batch[colnames(uncorrected$X), , drop = FALSE]
    expected.full <- do.call(
      playbase.preprocess::pp.batchCorrect,
      c(
        list(
          X = uncorrected$X,
          layers = NULL,
          target = expected.target,
          batch = expected.batch,
          method = method
        ),
        batch.args
      )
    )
    expected <- expected.full[rownames(corrected$X), , drop = FALSE]
    expected.options <- list(
      batch.correct.method = method,
      batch = expected.batch,
      target = expected.target,
      batch_args = batch.args
    )

    expect_identical(corrected$X, expected)
    expect_identical(corrected$counts, fixture$counts)
    expect_identical(corrected$settings$batch.correct.method, method)
    expect_identical(
      corrected$settings$options[names(expected.options)],
      expected.options
    )
    expect_identical(
      corrected$settings$preprocess$options[names(expected.options)],
      expected.options
    )
  }
})

test_that("sex and cycle inference sees the batch-corrected matrix", {
  skip_if_not_installed("limma")
  fixture <- make_pgx_fixture()
  batch <- data.frame(
    site = rep(c("one", "two"), 4L),
    row.names = colnames(fixture$counts)
  )
  target <- stats::setNames(fixture$samples$group, rownames(fixture$samples))
  observed.X <- NULL
  local_mocked_bindings(
    compute_cellcycle_gender = function(pgx) {
      observed.X <<- pgx$X
      pgx
    },
    .package = "playbase"
  )

  corrected <- create_test_pgx(
    fixture,
    options = list(
      normalize = FALSE,
      batch = batch,
      target = target
    ),
    batch.correct.method = "limma"
  )

  expect_identical(observed.X, corrected$X)
  expect_identical(corrected$settings$batch.correct.method, "limma")
})

test_that("multi-omics createPGX applies ComBat by layer before PGX filters", {
  skip_if_not_installed("sva")
  fixture <- make_pgx_fixture()
  rownames(fixture$counts) <- c(
    paste0("gx:g", 1:40),
    paste0("px:p", 1:40)
  )
  fixture$counts[10L, ] <- 0
  batch <- data.frame(
    site = rep(c("one", "two"), 4L),
    row.names = colnames(fixture$counts)
  )
  target <- stats::setNames(fixture$samples$group, rownames(fixture$samples))
  options <- list(
    normalize = FALSE,
    filter_missing = TRUE,
    remove_outliers = FALSE,
    target = target,
    batch = batch
  )
  create <- function(method, preprocessing = options, filter.genes = TRUE) {
    suppressMessages(playbase::pgx.createPGX(
      counts = fixture$counts,
      samples = fixture$samples,
      contrasts = fixture$contrasts,
      organism = "No organism",
      datatype = "multi-omics",
      preprocess = preprocessing,
      batch.correct.method = method,
      filter.genes = filter.genes,
      only.known = FALSE,
      only.proteincoding = FALSE,
      add.gmt = FALSE
    ))
  }
  uncorrected <- create("no_batch_correct", filter.genes = FALSE)
  corrected <- create("ComBat")
  expected.full <- playbase.preprocess::pp.batchCorrect(
    uncorrected$X,
    layers = sub(":.*", "", rownames(uncorrected$X)),
    target = target,
    batch = batch,
    method = "ComBat"
  )
  expected <- expected.full[rownames(corrected$X), , drop = FALSE]
  whole <- playbase.preprocess::pp.batchCorrect(
    uncorrected$X,
    layers = NULL,
    target = target,
    batch = batch,
    method = "ComBat"
  )

  expect_identical(corrected$X, expected)
  expected.options <- list(
    batch.correct.method = c(gx = "ComBat", px = "ComBat"),
    batch = batch,
    target = target,
    batch_args = list()
  )
  expect_identical(
    corrected$settings$options[names(expected.options)],
    expected.options
  )
  expect_identical(
    corrected$settings$preprocess$options[names(expected.options)],
    expected.options
  )
  expect_false(identical(corrected$X, whole[rownames(corrected$X), ]))
  expect_false("gx:g10" %in% rownames(corrected$X))
  expect_identical(corrected$counts, fixture$counts)
  expect_error(
    create(
      "no_batch_correct",
      utils::modifyList(
        options,
        list(
          batch_method = "limma"
        )
      )
    ),
    "select batch correction only with batch.correct.method"
  )
})

test_that("createPGX batches the final sample axis in one leaf call", {
  skip_if_not_installed("limma")
  fixture <- make_pgx_fixture()
  rownames(fixture$counts) <- c(
    paste0("gx:g", 1:40),
    paste0("px:p", 1:40)
  )
  fixture$contrasts[8L, 1L] <- NA
  fixture$samples$site <- rep(c("one", "two"), 4L)
  batch <- data.frame(
    site = fixture$samples$site,
    row.names = colnames(fixture$counts)
  )
  target <- stats::setNames(fixture$samples$group, rownames(fixture$samples))
  options <- list(
    normalize = FALSE,
    filter_missing = TRUE,
    remove_outliers = FALSE,
    target = target,
    batch_args = list(use_covariates = TRUE)
  )
  create <- function(preprocessing, method) {
    suppressMessages(playbase::pgx.createPGX(
      counts = fixture$counts,
      samples = fixture$samples,
      contrasts = fixture$contrasts,
      organism = "No organism",
      datatype = "multi-omics",
      preprocess = preprocessing,
      batch.correct.method = method,
      batch.pars = "site",
      prune.samples = TRUE,
      filter.genes = FALSE,
      only.known = FALSE,
      only.proteincoding = FALSE,
      add.gmt = FALSE
    ))
  }
  uncorrected <- create(options, "no_batch_correct")
  original.preprocess <- playbase.preprocess::pgx.preprocess
  leaf.calls <- list()
  local_mocked_bindings(
    pgx.preprocess = function(counts, ..., options = list()) {
      leaf.calls[[length(leaf.calls) + 1L]] <<- list(
        samples = colnames(counts),
        method = options$batch.correct.method,
        batch = options$batch,
        target = options$target
      )
      original.preprocess(counts = counts, ..., options = options)
    },
    .package = "playbase.preprocess"
  )

  corrected <- create(options, "limma")
  keep <- colnames(corrected$X)
  expected <- playbase.preprocess::pp.batchCorrect(
    uncorrected$X,
    layers = sub(":.*", "", rownames(uncorrected$X)),
    target = target[keep],
    batch = batch[keep, , drop = FALSE],
    method = "limma",
    use_covariates = TRUE
  )
  expected.options <- list(
    batch.correct.method = c(gx = "limma", px = "limma"),
    batch = batch[keep, , drop = FALSE],
    target = target[keep],
    batch_args = list(use_covariates = TRUE)
  )

  expect_length(leaf.calls, 1L)
  expect_identical(leaf.calls[[1L]]$samples, keep)
  expect_identical(leaf.calls[[1L]]$method, "limma")
  expect_identical(leaf.calls[[1L]]$batch, batch[keep, , drop = FALSE])
  expect_identical(leaf.calls[[1L]]$target, target[keep])
  expect_identical(corrected$X, expected)
  expect_identical(corrected$counts, fixture$counts)
  expect_identical(
    rownames(corrected$settings$preprocess$options$batch),
    keep
  )
  expect_identical(
    names(corrected$settings$preprocess$options$target),
    keep
  )
  expect_length(
    corrected$settings$preprocess$options$target,
    length(keep)
  )
  expect_identical(
    corrected$settings$options[names(expected.options)],
    expected.options
  )
  expect_identical(
    corrected$settings$preprocess$options[names(expected.options)],
    expected.options
  )
})

test_that("lifecycle sample pruning preserves source counts", {
  fixture <- make_pgx_fixture()
  fixture$contrasts[8, 1] <- NA
  pgx <- create_test_pgx(fixture, prune.samples = TRUE)

  expect_identical(pgx$counts, fixture$counts)
  expect_identical(colnames(pgx$X), colnames(fixture$counts)[1:7])
  expect_identical(pgx$settings$preprocess$alignment$cols, 1:7)
})

make_duplicate_pgx <- function(average.duplicated) {
  counts <- matrix(
    c(10, 20, 30, 40, 20, 30, 40, 50, 5, 10, 15, 20),
    nrow = 3L,
    byrow = TRUE,
    dimnames = list(c("dup", "dup", "other"), paste0("S", 1:4))
  )
  X <- log2(counts + 1)
  samples <- data.frame(
    group = c("a", "a", "b", "b"),
    row.names = colnames(counts)
  )
  contrasts <- matrix(
    samples$group,
    ncol = 1L,
    dimnames = list(colnames(counts), "b_vs_a")
  )
  pgx <- suppressMessages(playbase::pgx.createPGX(
    counts = counts,
    X = X,
    samples = samples,
    contrasts = contrasts,
    organism = "No organism",
    datatype = "RNA-seq",
    preprocess = NULL,
    average.duplicated = average.duplicated,
    add.gmt = FALSE,
    filter.genes = FALSE,
    only.known = FALSE,
    only.proteincoding = FALSE,
    convert.hugo = FALSE,
    prune.samples = FALSE
  ))
  list(counts = counts, X = X, pgx = pgx)
}

test_that("caller-supplied X averages duplicates without changing counts", {
  fixture <- make_duplicate_pgx(TRUE)
  expected <- playbase.preprocess::pp.deduplicate(
    fixture$X,
    method = "average",
    space = "log2"
  )

  expect_identical(
    serialize(fixture$pgx$counts, NULL),
    serialize(fixture$counts, NULL)
  )
  expect_identical(fixture$pgx$X, expected$X)
  expect_identical(
    fixture$pgx$settings$preprocess$alignment$rows,
    list(c(1L, 2L), 3L)
  )
})

test_that("preprocess NULL averages duplicates through the leaf", {
  fixture <- make_duplicate_pgx(TRUE)
  pgx <- suppressMessages(playbase::pgx.createPGX(
    counts = fixture$counts,
    samples = fixture$pgx$samples,
    contrasts = fixture$pgx$contrasts,
    organism = "No organism",
    datatype = "RNA-seq",
    preprocess = NULL,
    norm_method = "CPM",
    average.duplicated = TRUE,
    add.gmt = FALSE,
    filter.genes = FALSE,
    only.known = FALSE,
    only.proteincoding = FALSE,
    convert.hugo = FALSE,
    prune.samples = FALSE
  ))
  expected <- playbase.preprocess::pp.deduplicate(
    log2(fixture$counts + 1),
    method = "average",
    space = "log2"
  )

  expect_identical(pgx$counts, fixture$counts)
  expect_identical(pgx$X, expected$X)
  expect_identical(pgx$settings$options$dedup, "average")
  expect_identical(
    pgx$settings$preprocess$alignment$rows,
    list(c(1L, 2L), 3L)
  )
})

test_that("caller-supplied X uniquifies duplicates without changing counts", {
  fixture <- make_duplicate_pgx(FALSE)
  expected <- playbase.preprocess::pp.deduplicate(
    fixture$X,
    method = "unique",
    space = "log2"
  )

  expect_identical(
    serialize(fixture$pgx$counts, NULL),
    serialize(fixture$counts, NULL)
  )
  expect_identical(fixture$pgx$X, expected$X)
  expect_identical(
    fixture$pgx$settings$preprocess$alignment$rows,
    lapply(1:3, as.integer)
  )
})

test_that("max-feature bridge subsets X and final metadata together", {
  fixture <- make_pgx_fixture()
  pgx <- create_test_pgx(fixture)
  before <- pgx$counts
  computed <- suppressMessages(suppressWarnings(playbase::pgx.computePGX(
    pgx,
    max.genes = 20,
    gx.methods = "trend.limma",
    gset.methods = character(),
    extra.methods = character(),
    do.cluster = FALSE,
    do.clustergenes = FALSE,
    do.clustergenesets = FALSE
  )))

  expect_identical(computed$counts, before)
  expect_identical(nrow(computed$X), 20L)
  expect_length(computed$settings$preprocess$alignment$rows, 20L)
})
