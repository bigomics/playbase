## Canonical cross-package preprocessing boundary.
##
## Playbase re-exports only the ten matrix families and owns one pgx metadata
## wrapper. Removed compatibility names must not reappear in its namespace.

CANONICAL_PREPROCESS <- c(
  "pgx.preprocess",
  "pp.normalize",
  "pp.impute",
  "pp.filterFeatures",
  "pp.removeOutliers",
  "pp.deduplicate",
  "pp.batchCorrect",
  "pp.convertSpace",
  "pp.toCountScale",
  "pp.countScaleMatrix",
  "pp.alignCounts"
)

test_that("playbase exposes the canonical preprocessing surface", {
  exported <- getNamespaceExports("playbase")
  expect_identical(
    setdiff(CANONICAL_PREPROCESS, exported),
    character()
  )
  expect_length(CANONICAL_PREPROCESS, 11L)
})

test_that("the ten pp functions are direct leaf reexports", {
  leaf <- asNamespace("playbase.preprocess")
  playbase <- asNamespace("playbase")
  for (name in setdiff(CANONICAL_PREPROCESS, "pgx.preprocess")) {
    expect_identical(get(name, playbase), get(name, leaf), info = name)
  }
})

test_that("pgx.preprocess delegates with the seven-field result", {
  counts <- matrix(
    c(1, 2, 3, 4, 5, 6),
    nrow = 2,
    dimnames = list(c("g1", "g2"), c("s1", "s2", "s3"))
  )
  samples <- data.frame(
    group = c("a", "a", "b"),
    row.names = colnames(counts)
  )
  contrasts <- matrix(
    samples$group,
    ncol = 1L,
    dimnames = list(colnames(counts), "b_vs_a")
  )
  result <- pgx.preprocess(
    counts,
    samples,
    contrasts,
    options = list(normalize = FALSE)
  )

  expect_named(
    result,
    c("counts", "X", "annot", "prior", "space", "alignment", "options")
  )
  expect_identical(result$counts, counts)
  expect_identical(
    playbase.preprocess::pp.alignCounts(
      result$counts,
      result$alignment,
      X = result$X
    ),
    result$counts
  )
})

test_that("removed preprocessing compatibility names stay absent", {
  removed <- c(
    "getPrior",
    "logCPM",
    "normalizeExpression",
    "normalizeMultiOmics",
    "imputeMissing",
    "imputeMissing.mox",
    "svdImpute2",
    "detectOutlierSamples",
    "plotOutlierScores",
    "betaToM",
    "mToBeta",
    "pgx.countNormalization",
    "pgx.countScaleMatrix",
    "pgx.recomputeCounts",
    "pgx.alignXtoCounts",
    "pgx.ranWithCorrection",
    "pgx.removeLowVariance",
    "limmaCorrect",
    "combatCorrect",
    "ruvCorrect",
    "svaCorrect",
    "nnmCorrect"
  )
  present <- removed[vapply(
    removed,
    exists,
    logical(1),
    envir = asNamespace("playbase"),
    inherits = FALSE
  )]
  expect_identical(present, character())
})
