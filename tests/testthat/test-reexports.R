## Cross-package regression for the preprocessing extraction.
##
## The preprocessing layer lives in playbase.preprocess and is re-exported
## from playbase by R/reexport-preprocess.R. These tests pin the promise the
## port was built on: every public name still resolves from playbase::, and
## pgx.preprocess() still accepts the historical (samples, contrasts)
## signature even though the child package now takes a `groups` vector.
##
## If one of these fails, an omicsplayground call site is broken.

## The complete public surface moved to playbase.preprocess. 25 names:
## 24 bare re-exports plus pgx.preprocess(), which is a real wrapper.
REEXPORTED <- c(
  "betaToM", "detectOutlierSamples", "getPrior", "imputeMissing",
  "imputeMissing.mox", "is.multiomics", "log1s", "logCPM", "mToBeta",
  "maxMedianNormalization", "maxSumNormalization", "mofa.get_prefix",
  "nmfImpute", "nmfImpute2", "normalizeExpression", "normalizeMethylation",
  "normalizeMultiOmics", "normalizeRLE", "normalizeTMM", "perseusImpute",
  "pgx.countNormalization", "pgx.preprocess", "plotOutlierScores",
  "referenceNormalization", "svdImpute2"
)

test_that("all 25 preprocessing names are exported from playbase", {
  exported <- getNamespaceExports("playbase")
  missing <- setdiff(REEXPORTED, exported)
  ## named so a failure says which one went missing
  expect_identical(missing, character(0))
  expect_length(REEXPORTED, 25L)
})

test_that("each re-exported name resolves to a function via playbase::", {
  for (nm in REEXPORTED) {
    obj <- get(nm, envir = asNamespace("playbase"))
    expect_true(is.function(obj), info = nm)
  }
})

test_that("re-exports are the child's objects, not stale playbase copies", {
  skip_if_not_installed("playbase.preprocess")
  child <- asNamespace("playbase.preprocess")
  ## pgx.preprocess is deliberately NOT identical -- it is a wrapper.
  bare <- setdiff(REEXPORTED, "pgx.preprocess")
  for (nm in bare) {
    expect_identical(
      get(nm, envir = asNamespace("playbase")),
      get(nm, envir = child),
      info = nm
    )
  }
})

test_that("pgx.preprocess keeps the historical signature", {
  ## omicsplayground and pgx.createPGX() both call positionally.
  expect_identical(
    names(formals(playbase::pgx.preprocess)),
    c("counts", "samples", "contrasts", "annot", "options")
  )
  ## the child underneath takes `groups` instead
  skip_if_not_installed("playbase.preprocess")
  expect_true("groups" %in% names(formals(playbase.preprocess::pgx.preprocess)))
})

## ---------------------------------------------------------------- NULL cases
##
## pgx.createPGX() validates only counts/samples/organism, so `contrasts` is
## legitimately NULL. Before the port that reached
## contrasts.convertToLabelMatrix(NULL, samples), which stops(). The wrapper
## must degrade to groups = NULL instead of erroring.

mk_counts <- function(n = 80, m = 6, seed = 11) {
  set.seed(seed)
  counts <- matrix(stats::rpois(n * m, lambda = 50), n, m)
  dimnames(counts) <- list(paste0("g", seq_len(n)), paste0("s", seq_len(m)))
  counts
}

mk_samples <- function(counts) {
  data.frame(
    group = rep(c("a", "b"), length.out = ncol(counts)),
    row.names = colnames(counts)
  )
}

test_that("wrapper tolerates NULL samples and NULL contrasts", {
  counts <- mk_counts()
  samples <- mk_samples(counts)

  cases <- list(
    "both NULL"        = list(samples = NULL,    contrasts = NULL),
    "samples only"     = list(samples = samples, contrasts = NULL),
    "neither supplied" = NULL
  )

  for (nm in names(cases)) {
    a <- cases[[nm]]
    res <- if (is.null(a)) {
      playbase::pgx.preprocess(counts)
    } else {
      playbase::pgx.preprocess(counts, samples = a$samples, contrasts = a$contrasts)
    }
    expect_type(res, "list")
    expect_named(res, c("counts", "X", "annot", "prior"), ignore.order = TRUE, info = nm)
    expect_identical(dimnames(res$counts), dimnames(res$X), info = nm)
    expect_true(all(is.finite(res$X)), info = nm)
  }
})

test_that("NULL contrasts still filters -- it does not silently no-op or wipe", {
  counts <- mk_counts()
  counts[1:10, ] <- NA               ## 10 all-NA features
  counts[11:15, 1:5] <- NA           ## 5 features with a single observation

  ## f >= 1 selects the per-group count rule with a hardcoded cutoff of 3.
  ## With no design, all samples form one group. The all-NA rows must go and
  ## the fully observed ones must stay: a no-op would keep all 80, a wipe
  ## would keep none.
  res <- playbase::pgx.preprocess(
    counts, samples = NULL, contrasts = NULL,
    options = list(filter_missing = TRUE, filter_threshold = 3)
  )
  expect_gt(nrow(res$X), 0)
  expect_lt(nrow(res$X), nrow(counts))
  expect_false(any(paste0("g", 1:10) %in% rownames(res$X)))
  expect_true(all(paste0("g", 16:80) %in% rownames(res$X)))
})
