#' Test for compute_cellcycle_gender
test_that("compute_cellcycle_gender adds cell cycle and gender data", {
  # Create mock data
  pgx <- list(
    samples = playbase::SAMPLES,
    counts = as.matrix(playbase::COUNTS),
    genes = playbase::GENES_TABLE,
    organism = "Human"
  )
  pgx$counts <- pgx$counts[!duplicated(rownames(pgx$counts)), , drop = FALSE]
  pgx$genes <- pgx$genes[!pgx$genes$symbol == "", , drop = FALSE]
  pgx$counts <- pgx$counts[rownames(pgx$counts) %in% pgx$genes$symbol, , drop = FALSE]

  result <- playbase::compute_cellcycle_gender(pgx, pgx$counts)

  # Check cell cycle and gender added to sample df
  expect_true(".cell_cycle" %in% colnames(result$samples))
  expect_true(".gender" %in% colnames(result$samples))
})

test_that("compute_cellcycle_gender adds cell cycle and gender data", {
  # Create mock data
  pgx <- list(
    samples = playbase::SAMPLES,
    counts = as.matrix(playbase::COUNTS),
    genes = playbase::GENES_TABLE,
    organism = "Human"
  )
  pgx$counts <- pgx$counts[!duplicated(rownames(pgx$counts)), , drop = FALSE]
  pgx$genes <- pgx$genes[!pgx$genes$symbol == "", , drop = FALSE]
  pgx$counts <- pgx$counts[rownames(pgx$counts) %in% pgx$genes$symbol, , drop = FALSE]
  pgx$X <- pgx$counts
  pgx$settings <- list(
    preprocess = .pgx_identity_preprocess_metadata(
      pgx$counts,
      pgx$X,
      space = "counts",
      prior = NA_real_
    )
  )

  result <- playbase::compute_cellcycle_gender(pgx)
  expect_cell_cycle_stages <- c(
    "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1",
    "S", "S", "S", "S", "S", "S", "S", "S"
  )
  # Check cell cycle and gender added to sample df
  expect_equal(result$samples$.cell_cycle, expect_cell_cycle_stages)
})

test_that("inferred phenotypes follow retained analysis sample names", {
  source_counts <- matrix(
    c(10, 1, 1, 10, 1, 10, 10, 1),
    nrow = 2L,
    dimnames = list(c("XIST", "DDX3Y"), paste0("s", 1:4))
  )
  samples <- data.frame(
    group = rep("a", 4L),
    row.names = colnames(source_counts)
  )
  pgx <- list(
    counts = source_counts,
    samples = samples,
    genes = data.frame(
      symbol = rownames(source_counts),
      row.names = rownames(source_counts)
    ),
    organism = "Human"
  )
  analysis_counts <- source_counts[, c("s3", "s1", "s4"), drop = FALSE]
  testthat::local_mocked_bindings(
    pgx.inferCellCyclePhase = function(counts) c("S", "G1", "G2M"),
    .package = "playbase"
  )

  result <- suppressMessages(playbase::compute_cellcycle_gender(
    pgx,
    rna.counts = analysis_counts
  ))

  expect_identical(result$counts, source_counts)
  expect_identical(rownames(result$samples), rownames(samples))
  expect_identical(result$samples$.cell_cycle, c("G1", NA, "S", "G2M"))
  expect_identical(result$samples$.gender, c("F", NA, "M", "F"))
})

test_that("inferred phenotype assignment rejects ambiguous sample mappings", {
  samples <- data.frame(group = rep("a", 4L), row.names = paste0("s", 1:4))
  result <- .pgx_assign_inferred_samples(
    samples,
    ".gender",
    c(s3 = "M", s1 = "F", s4 = "F"),
    c("s1", "s3", "s4")
  )

  expect_identical(result$.gender, c("F", NA, "M", "F"))
  expect_error(
    .pgx_assign_inferred_samples(
      samples,
      ".gender",
      c("F", "M"),
      c("s1", "s1")
    ),
    "processed samples must match"
  )
  expect_error(
    .pgx_assign_inferred_samples(samples, ".gender", "F", "unknown"),
    "processed samples must match"
  )
  expect_error(
    .pgx_assign_inferred_samples(samples, ".gender", "F", c("s1", "s3")),
    "prediction length must match"
  )
})

test_that("createPGX infers gender after outlier removal without narrowing source data", {
  counts <- as.matrix(playbase::COUNTS)
  samples <- playbase::SAMPLES
  samples$.cell_cycle <- rep("G1", nrow(samples))
  testthat::local_mocked_bindings(
    getProbeAnnotation = function(..., probes) {
      data.frame(symbol = probes, row.names = probes)
    },
    pgx.inferGender = function(X, gene_name = NULL) {
      rep(c("F", "M"), length.out = ncol(X))
    },
    .package = "playbase"
  )

  pgx <- suppressMessages(playbase::pgx.createPGX(
    counts = counts,
    samples = samples,
    contrasts = playbase::CONTRASTS,
    organism = "Human",
    datatype = "RNA-seq",
    preprocess = list(
      norm_method = "CPM",
      remove_outliers = TRUE,
      outlier_threshold = 3,
      impute = FALSE
    ),
    add.gmt = FALSE,
    only.known = FALSE,
    only.proteincoding = FALSE
  ))

  removed <- setdiff(colnames(counts), colnames(pgx$X))
  expect_length(removed, 1L)
  expect_identical(pgx$counts, counts)
  expect_identical(rownames(pgx$samples), rownames(samples))
  expect_true(".gender" %in% names(pgx$samples))
  expect_true(is.na(pgx$samples[removed, ".gender"]))
  expect_identical(
    pgx$samples[colnames(pgx$X), ".gender"],
    rep(c("F", "M"), length.out = ncol(pgx$X))
  )
  expect_identical(
    dimnames(playbase.preprocess::pp.alignCounts(
      pgx$counts,
      pgx$settings$preprocess$alignment,
      X = pgx$X
    )),
    dimnames(pgx$X)
  )

  computed <- withCallingHandlers(
    suppressMessages(playbase::pgx.computePGX(
      pgx,
      gx.methods = "trend.limma",
      gset.methods = character(),
      do.cluster = FALSE,
      do.clustergenes = FALSE,
      do.clustergenesets = FALSE,
      extra.methods = character(),
      ai_features = NULL
    )),
    warning = function(w) {
      if (
        grepl(
          "very small variances detected",
          conditionMessage(w),
          fixed = TRUE
        )
      ) {
        invokeRestart("muffleWarning")
      }
    }
  )
  expect_identical(computed$counts, counts)
  expect_identical(colnames(computed$X), colnames(pgx$X))
})

test_that("ComBat createPGX and compute retain source data after outlier inference", {
  skip_if_not_installed("sva")
  counts <- as.matrix(playbase::COUNTS)
  samples <- playbase::SAMPLES
  samples$.cell_cycle <- rep("G1", nrow(samples))
  samples$site <- rep(c("one", "two"), length.out = nrow(samples))
  testthat::local_mocked_bindings(
    getProbeAnnotation = function(..., probes) {
      data.frame(symbol = probes, row.names = probes)
    },
    pgx.inferGender = function(X, gene_name = NULL) {
      rep(c("F", "M"), length.out = ncol(X))
    },
    .package = "playbase"
  )

  pgx <- suppressMessages(playbase::pgx.createPGX(
    counts = counts,
    samples = samples,
    contrasts = playbase::CONTRASTS,
    organism = "Human",
    datatype = "RNA-seq",
    preprocess = list(
      norm_method = "CPM",
      remove_outliers = TRUE,
      outlier_threshold = 3,
      impute = FALSE
    ),
    batch.correct.method = "ComBat",
    batch.pars = "site",
    add.gmt = FALSE,
    only.known = FALSE,
    only.proteincoding = FALSE
  ))

  removed <- setdiff(colnames(counts), colnames(pgx$X))
  expect_length(removed, 1L)
  expect_identical(pgx$counts, counts)
  expect_identical(rownames(pgx$samples), rownames(samples))
  expect_identical(pgx$settings$batch.correct.method, "ComBat")
  expect_true(is.na(pgx$samples[removed, ".gender"]))

  computed <- withCallingHandlers(
    suppressMessages(playbase::pgx.computePGX(
      pgx,
      gx.methods = "trend.limma",
      gset.methods = character(),
      do.cluster = FALSE,
      do.clustergenes = FALSE,
      do.clustergenesets = FALSE,
      extra.methods = character(),
      ai_features = NULL
    )),
    warning = function(w) {
      if (
        grepl(
          "very small variances detected",
          conditionMessage(w),
          fixed = TRUE
        )
      ) {
        invokeRestart("muffleWarning")
      }
    }
  )
  expect_identical(computed$counts, counts)
  expect_identical(colnames(computed$X), colnames(pgx$X))
})
