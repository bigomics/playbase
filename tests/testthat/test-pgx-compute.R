#' Test for pgx.createFromFiles
#'
#'

# pgx <- playbase::pgx.createPGX(
#     samples = playbase::SAMPLES,
#     counts = playbase::COUNTS,
#     contrasts = playbase::CONTRASTS[,1:2],
#     organism = "Human"
#   )

#' Test for pgx.createPGX
test_that("pgx.createPGX produce all pgx slots", {
  skip("high memory usage, maybe use mini-example?")
  # Call example data
  # pgx_data <- playbase::get_mini_example_data()

  # Create expected outputs
  expected_tests <- c(
    "name", "organism", "version", "date", "creator", "datatype",
    "description", "samples", "counts", "contrasts", "X",
    "total_counts", "counts_multiplier", "genes", "all_genes",
    "probe_type", "filtered", "tsne2d", "tsne3d", "cluster", "cluster.genes"
  )
  total_counts <- apply(playbase::COUNTS, 2, sum)

  gene_table <- data.frame(
    symbol = c("A1BG", "AGAP2", "ANXA4", "ARPC1A", "BATF", "C19orf53"),
    gene_title = c(
      "alpha-1-B glycoprotein", "ArfGAP with GTPase domain, ankyrin repeat and PH domain 2",
      "annexin A4", "actin related protein 2/3 complex subunit 1A", "basic leucine zipper ATF-like transcription factor",
      "chromosome 19 open reading frame 53"
    ),
    gene_biotype = rep("protein_coding", 6),
    chr = c("19", "12", "2", "7", "14", "19"),
    pos = c(58345178, 57723761, 69644425, 99325898, 75522455, 13774456),
    tx_len = c(2134, 5388, 905, 1582, 617, 897),
    map = c("q13.43", "q14.1", "p13.3", "q22.1", "q24.3", "p13.13"),
    source = c(
      "Source:HGNC Symbol;Acc:HGNC:5", "Source:HGNC Symbol;Acc:HGNC:16921", "Source:HGNC Symbol;Acc:HGNC:542",
      "Source:HGNC Symbol;Acc:HGNC:703", "Source:HGNC Symbol;Acc:HGNC:958", "Source:HGNC Symbol;Acc:HGNC:24991"
    ),
    gene_name = c("A1BG", "AGAP2", "ANXA4", "ARPC1A", "BATF", "C19orf53"),
    feature = c("A1BG", "AGAP2", "ANXA4", "ARPC1A", "BATF", "C19orf53")
  )
  rownames(gene_table) <- gene_table$symbol

  # Check output
  ## Check all te test present
  expect_true(all(names(pgx) == expected_tests))

  ## Check the total counts in each sample
  expect_equal(pgx$total_counts, total_counts)

  ## Check multiplier
  expect_equal(pgx$counts_multiplier, 1)

  ## Check contrast
  expect_equal(dim(pgx$contrasts), c(18, 2))


  ## Check that the gene info is generated correctly
  # TODO expect_equal(pgx$genes[c(1, 10, 20, 30, 40, 50), , drop = FALSE], gene_table)

  ## Check cluster.genes
  expect_equal(dim(pgx$cluster$pos$pca2d), c(ncol(playbase::COUNTS), 2))
  expect_equal(dim(pgx$cluster$pos$tsne3d), c(ncol(playbase::COUNTS), 3))
  expect_equal(dim(pgx$cluster$pos$umap2d), c(ncol(playbase::COUNTS), 2))
  expect_equal(dim(pgx$cluster$pos$umap3d), c(ncol(playbase::COUNTS), 3))
})


#' Test for pgx.computePGX
test_that("pgx.computePGX runs without errors", {
  skip("high memory usage, maybe use mini-example?")
  # Run function
  suppressWarnings(pgx_comp <- playbase::pgx.computePGX(pgx))

  # Expected outputs
  expected_slots <- c(
    "name", "organism", "version", "date", "creator", "datatype", "description", "samples",
    "counts", "contrasts", "X", "total_counts", "counts_multiplier", "genes",
    "all_genes", "probe_type", "tsne2d", "tsne3d", "cluster", "cluster.genes",
    "model.parameters", "filtered", "timings", "gx.meta", "gset.meta", "gsetX",
    "GMT", "cluster.gsets", "meta.go"
  )
  # Check output
  expect_equal(names(pgx_comp), expected_slots)
})


#' D-24 / D-42: what pgx.computePGX() is allowed to do to pgx$counts, which is
#' nothing. The `max.genes` shrink takes rows out of `X`; the re-cut that used
#' to follow it (`gg <- intersect(rownames(pgx$counts), rownames(pgx$X))`) is
#' gone. It mattered because upload_server.R:1093-1095 seeds Reanalyse from the
#' COMPUTED object, so a `counts` narrowed by compute made every recompute start
#' from the previous recompute's output -- the ratchet playbase-lh8 measured.
#' The second round below is that Reanalyse, and it is the assertion with teeth.
test_that("pgx.computePGX shrinks X and leaves counts at the upload's shape", {
  counts <- as.matrix(playbase::COUNTS)
  opts <- list(
    datatype = "RNA-seq", norm_method = "CPM",
    remove_outliers = TRUE, outlier_threshold = 2, impute = FALSE
  )
  create <- function(counts, samples, contrasts) {
    suppressMessages(playbase::pgx.createPGX(
      counts = counts, samples = samples, contrasts = contrasts,
      organism = "Human", datatype = "RNA-seq", preprocess = opts,
      add.gmt = FALSE, convert.hugo = FALSE,
      filter.genes = FALSE, only.known = FALSE, only.proteincoding = FALSE
    ))
  }
  compute <- function(pgx) {
    suppressMessages(suppressWarnings(playbase::pgx.computePGX(
      pgx, max.genes = 500, gx.methods = "trend.limma", gset.methods = c(),
      extra.methods = c(), do.cluster = FALSE, do.clustergenes = FALSE,
      do.clustergenesets = FALSE
    )))
  }

  pgx <- compute(create(counts, playbase::SAMPLES, playbase::CONTRASTS))

  ## non-vacuity: both removals really fired, so the shapes really do split
  expect_lt(ncol(pgx$X), ncol(counts))
  expect_lt(nrow(pgx$X), nrow(counts))

  ## counts is the upload, whole, including dimnames and storage mode
  expect_identical(pgx$counts, counts)

  ## the shrink is what capped X, and its audit trail names the rows it dropped
  expect_identical(nrow(pgx$X), 500L)
  dropped <- strsplit(pgx$filtered[["low.variance"]], ";")[[1]]
  expect_identical(length(dropped), nrow(counts) - 500L)
  expect_identical(intersect(dropped, rownames(pgx$X)), character(0))

  ## X's samples and features are a subset of the ones counts still carries
  expect_true(all(colnames(pgx$X) %in% colnames(pgx$counts)))
  expect_true(all(rownames(pgx$X) %in% rownames(pgx$counts)))

  ## Reanalyse: the computed object seeds the next upload. Round two must land
  ## on round one, not inside it.
  re <- compute(create(pgx$counts, pgx$samples, playbase::CONTRASTS))
  expect_identical(re$counts, counts)
  expect_identical(dim(re$X), dim(pgx$X))
  expect_identical(rownames(re$X), rownames(pgx$X))
})

test_that("pgx.countScaleMatrix is the one back-transform rule", {
  set.seed(1)
  counts <- matrix(rpois(200, 50) + 1, 20, 10,
    dimnames = list(paste0("g", 1:20), paste0("s", 1:10))
  )

  ## No provenance record: pgx.ranWithCorrection() says "we do not know", which
  ## is not "it ran", so the upload is returned untouched. This is the answer
  ## every object playbase holds today gets (D-37).
  legacy <- list(counts = counts, X = log2(1 + counts))
  expect_true(is.na(playbase.preprocess::pgx.ranWithCorrection(legacy)))
  expect_identical(playbase::pgx.countScaleMatrix(legacy), counts)

  ## A real record, no correction in its history: still the upload.
  pgx <- playbase.preprocess::pgx.transform(
    list(counts = counts, X = NULL),
    to = "log2", prior = 1
  )
  expect_false(playbase.preprocess::pgx.ranWithCorrection(pgx))
  expect_identical(playbase::pgx.countScaleMatrix(pgx), counts)

  ## The same record with a batchCorrect step: the reconstruction, not the
  ## upload. No verb in either package writes this step today (D-13/D-37), so
  ## it is built by hand -- the branch is live code, not a branch that cannot
  ## be reached.
  pp <- pgx$settings$preprocessing
  corrected <- pgx
  corrected$settings$preprocessing <- playbase.preprocess:::new_pgx_preprocessing(
    space = pp$space, layers = pp$layers, prior = pp$prior,
    invertible = pp$invertible,
    history = c(pp$history, list(list(verb = "batchCorrect", method = "ComBat"))),
    sealed_at = pp$sealed_at, engine = pp$engine
  )
  expect_true(playbase.preprocess::pgx.ranWithCorrection(corrected))
  expect_identical(
    playbase::pgx.countScaleMatrix(corrected),
    playbase.preprocess::pgx.recomputeCounts(corrected)
  )

  ## And it refuses rather than guessing when the record cannot be replayed
  ## against the counts it is handed.
  broken <- corrected
  broken$counts <- counts[1:5, , drop = FALSE]
  expect_error(playbase::pgx.countScaleMatrix(broken), "pgx.alignXtoCounts")
})
