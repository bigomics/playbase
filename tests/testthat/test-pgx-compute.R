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
  prune.samples = FALSE
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
    add.gmt = FALSE,
    filter.genes = filter.genes,
    only.known = FALSE,
    only.proteincoding = FALSE,
    convert.hugo = FALSE,
    prune.samples = prune.samples
  ))
}

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
