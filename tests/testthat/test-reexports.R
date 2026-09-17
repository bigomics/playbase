## Completed preprocessing namespace boundary.
##
## Playbase must consume the leaf through explicit namespace calls. It neither
## re-exports pp.* families nor defines a second preprocessing orchestrator.

RETIRED_PLAYBASE_PREPROCESS_EXPORTS <- c(
  "pgx.preprocess",
  "pp.alignCounts",
  "pp.batchCorrect",
  "pp.convertSpace",
  "pp.countScaleMatrix",
  "pp.deduplicate",
  "pp.filterFeatures",
  "pp.impute",
  "pp.inferLayers",
  "pp.alignLayers",
  "pp.normalize",
  "pp.removeOutliers",
  "pp.groupsFromContrasts",
  "pp.batchFromSamples",
  "pp.scaleCounts",
  "pp.toCountScale",
  "imputeMedian",
  "counts.mergeDuplicateFeatures",
  "counts.autoScaling",
  "pgx.filterZeroCounts",
  "pgx.filterLowExpressed"
)

test_that("playbase does not export preprocessing implementations", {
  exported <- getNamespaceExports("playbase")
  leaf_exports <- getNamespaceExports("playbase.preprocess")
  forbidden <- union(RETIRED_PLAYBASE_PREPROCESS_EXPORTS, leaf_exports)
  present <- intersect(forbidden, exported)
  expect_identical(
    present,
    character()
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
