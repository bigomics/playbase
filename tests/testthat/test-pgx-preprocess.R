## Parity harness for pgx.preprocess().
##
## `ref_normalize()` below is an INDEPENDENT, line-for-line copy of the Omics
## Playground upload module reactives it replaces:
##   omicsplayground/components/board.upload/R/upload_module_normalization.R
##   imputedX (L49-131) -> normalizedX (L134-162) -> cleanX (L165-191)
## If pgx.preprocess() drifts from that logic, these tests fail. When the app is
## later switched to call pgx.preprocess(), this locks its behavior unchanged.

## ---- golden reference: the app's former pipeline, de-reactived ----
ref_normalize <- function(counts, samples, contrasts, annot = NULL, opt = list()) {
  counts <- as.matrix(counts)

  ## imputedX
  counts[which(is.nan(counts))] <- NA
  counts[which(is.infinite(counts))] <- NA
  if (isTRUE(opt$is_npx)) counts <- 2**counts
  if (any(counts < 0, na.rm = TRUE)) counts <- pmax(counts, 0)
  if (isTRUE(opt$zero_as_na)) counts[which(counts == 0)] <- NA

  is.mox <- playbase::is.multiomics(rownames(counts))
  if (is.mox) {
    X <- counts
    dtypes <- unique(sub(":.*", "", rownames(X)))
    for (i in seq_along(dtypes)) {
      ii <- grep(paste0("^", dtypes[i], ":"), rownames(counts))
      prior <- 1
      if (dtypes[i] != "gx") prior <- playbase::getPrior(counts[ii, ])
      X[ii, ] <- log2(counts[ii, ] + prior)
    }
  } else if (opt$datatype == "methylomics") {
    X <- playbase::mToBeta(counts)
    prior <- 0
  } else {
    prior0 <- playbase::getPrior(counts)
    m <- opt$norm_method
    prior <- ifelse(grepl("CPM|TMM", m), 1, prior0)
    X <- log2(counts + prior)
  }

  if (sum(is.na(X)) > 0 && isTRUE(opt$filter_missing)) {
    f <- opt$filter_threshold
    sample.contrasts <- playbase::contrasts.convertToLabelMatrix(contrasts, samples)
    grp <- apply(sample.contrasts, 1, paste, collapse = "_")
    if (f >= 1) {
      grp.sum <- tapply(seq_len(ncol(counts)), grp, function(i) rowSums(!is.na(counts[, i, drop = FALSE])))
      maxsum <- apply(do.call(cbind, grp.sum), 1, max, na.rm = TRUE)
      sel <- (maxsum >= 3)
    } else if (f < 0) {
      grp.avg <- tapply(seq_len(ncol(counts)), grp, function(i) rowMeans(!is.na(counts[, i, drop = FALSE])))
      maxavg <- apply(do.call(cbind, grp.avg), 1, max, na.rm = TRUE)
      sel <- (maxavg >= 0.5)
    } else {
      sel <- (rowMeans(is.na(X)) <= f)
    }
    X <- X[which(sel), , drop = FALSE]
    counts <- counts[which(sel), , drop = FALSE]
    if (!is.null(annot)) annot <- annot[which(sel), , drop = FALSE]
  }

  if (any(is.na(X)) && isTRUE(opt$impute)) {
    if (is.mox) {
      X <- playbase::imputeMissing.mox(X, method = opt$impute_method)
    } else {
      X <- playbase::imputeMissing(X, method = opt$impute_method)
    }
  }

  imputedX <- X ## what the wizard's imputedX() reactive shows

  ## normalizedX
  if (isTRUE(opt$normalize)) {
    if (opt$datatype == "multi-omics") {
      X <- playbase::normalizeMultiOmics(X)
    } else if (opt$datatype == "methylomics") {
      nX <- try(playbase::normalizeMethylation(X, opt$norm_method, opt$meth_type), silent = TRUE)
      if (!inherits(nX, "try-error") && !is.null(nX)) X <- nX
    } else {
      X <- playbase::normalizeExpression(X, method = opt$norm_method, ref = opt$ref_gene, prior = prior)
    }
  }

  normalizedX <- X ## what the wizard's normalizedX() reactive shows

  ## cleanX
  kk <- intersect(rownames(X), rownames(counts))
  X <- X[kk, , drop = FALSE]
  counts <- counts[kk, , drop = FALSE]
  if (isTRUE(opt$remove_outliers)) {
    if (sum(is.na(X)) > 0) {
      if (is.mox) X <- playbase::imputeMissing.mox(X, method = "SVD2") else X <- playbase::imputeMissing(X, method = "SVD2")
    }
    res <- playbase::detectOutlierSamples(X, plot = FALSE)
    is.outlier <- (res$z.outlier > opt$outlier_threshold)
    if (any(is.outlier) && !all(is.outlier)) {
      X <- X[, which(!is.outlier), drop = FALSE]
      counts <- counts[, colnames(X), drop = FALSE]
    }
  }
  list(
    counts = counts, X = X, annot = annot, prior = prior,
    imputedX = imputedX, normalizedX = normalizedX
  )
}

defaults <- list(
  datatype = "RNA-seq", is_npx = FALSE, zero_as_na = FALSE, normalize = TRUE,
  norm_method = "CPM", ref_gene = NULL, filter_missing = FALSE, filter_threshold = 3,
  impute = FALSE, impute_method = "SVD2", remove_outliers = FALSE, outlier_threshold = 3,
  meth_type = NULL
)

## Assert pgx.preprocess() == app reference for a given option set. If computing
## the reference errors in this environment (missing optional dep), skip.
expect_parity <- function(counts, samples, contrasts, opt) {
  o <- utils::modifyList(defaults, opt)
  ref <- tryCatch(ref_normalize(counts, samples, contrasts, opt = o), error = function(e) e)
  if (inherits(ref, "error")) testthat::skip(paste("reference unavailable:", conditionMessage(ref)))
  got <- pgx.preprocess(counts, samples, contrasts, options = o)
  testthat::expect_equal(got$X, ref$X, tolerance = 1e-8)
  ## `counts` is no longer cut down to `X` (D-24/D-39): it comes back at the
  ## shape it went in, whatever the removals did. The app reference still cuts
  ## it, so parity on counts is now "the rows and samples the reference kept",
  ## located through the alignment rather than by name.
  a <- playbase.preprocess::pgx.alignXtoCounts(list(counts = got$counts, X = got$X))
  testthat::expect_equal(got$counts[a$rows, a$cols, drop = FALSE], ref$counts,
    tolerance = 1e-8)
  ## The tripwire that replaces "counts and X stay row/col-aligned": every row
  ## and sample of X is still one counts has, and the alignment says which.
  testthat::expect_identical(rownames(got$counts)[a$rows], rownames(got$X))
  testthat::expect_identical(colnames(got$counts)[a$cols], colnames(got$X))
}

get_fixture <- function() {
  if (!all(sapply(c("COUNTS", "SAMPLES", "CONTRASTS"), exists, where = asNamespace("playbase")))) {
    testthat::skip("playbase COUNTS/SAMPLES/CONTRASTS not available")
  }
  list(
    counts = playbase::COUNTS,
    samples = playbase::SAMPLES,
    contrasts = playbase::CONTRASTS[, 1:2, drop = FALSE]
  )
}

test_that("pgx.preprocess reproduces the app pipeline across normalization methods", {
  fx <- get_fixture()
  ## keep in sync with normalizeExpression()'s method list
  for (m in c("CPM", "TMM", "quantile", "CPM+quantile", "maxMedian", "maxSum")) {
    expect_parity(fx$counts, fx$samples, fx$contrasts, list(norm_method = m))
  }
  ## 'reference' needs a valid ref feature; both paths agree once one is given.
  expect_parity(fx$counts, fx$samples, fx$contrasts,
    list(norm_method = "reference", ref_gene = rownames(fx$counts)[1]))
})

test_that("pgx.preprocess errors on reference normalization without a ref gene", {
  fx <- get_fixture()
  expect_error(
    pgx.preprocess(fx$counts, fx$samples, fx$contrasts,
      options = list(norm_method = "reference")),
    "ref_gene"
  )
})

test_that("pgx.preprocess reproduces the app pipeline with normalization skipped", {
  fx <- get_fixture()
  expect_parity(fx$counts, fx$samples, fx$contrasts, list(normalize = FALSE))
})

test_that("pgx.preprocess reproduces zero-as-NA + missingness filter + imputation", {
  fx <- get_fixture()
  expect_parity(fx$counts, fx$samples, fx$contrasts, list(zero_as_na = TRUE))
  expect_parity(fx$counts, fx$samples, fx$contrasts, list(zero_as_na = TRUE, filter_missing = TRUE, filter_threshold = 3))
  expect_parity(fx$counts, fx$samples, fx$contrasts, list(zero_as_na = TRUE, impute = TRUE, impute_method = "SVD2"))
})

test_that("pgx.preprocess reproduces outlier-sample removal", {
  fx <- get_fixture()
  expect_parity(fx$counts, fx$samples, fx$contrasts, list(remove_outliers = TRUE, outlier_threshold = 3))
})

## ---- the preview stages ----
## The upload wizard previews the pipeline at three points -- after imputation,
## after normalization, after outlier removal -- and computes each with its own
## copy of the code, so what the user approves is an approximation of what then
## runs. Each point is the same chain stopped early, which pgx.preprocess()
## expresses by switching off the steps that come after it; nothing else has to
## exist for the wizard to show the real computation. This asserts that, stage
## by stage, against the same golden reference the final matrix is checked
## against.
expect_preview_parity <- function(counts, samples, contrasts, opt) {
  o <- utils::modifyList(defaults, opt)
  ref <- tryCatch(ref_normalize(counts, samples, contrasts, opt = o), error = function(e) e)
  if (inherits(ref, "error")) testthat::skip(paste("reference unavailable:", conditionMessage(ref)))
  stage <- function(...) {
    pgx.preprocess(counts, samples, contrasts,
      options = utils::modifyList(o, list(...)))
  }
  imputed <- stage(normalize = FALSE, remove_outliers = FALSE)
  normalized <- stage(remove_outliers = FALSE)
  clean <- stage()
  testthat::expect_equal(imputed$X, ref$imputedX, tolerance = 1e-8)
  testthat::expect_equal(normalized$X, ref$normalizedX, tolerance = 1e-8)
  testthat::expect_equal(clean$X, ref$X, tolerance = 1e-8)
  ## The prior the preview panels label their axes with is the same one.
  testthat::expect_equal(imputed$prior, ref$prior)
  ## Not vacuous: the three stages must actually differ, or this would pass on
  ## a pipeline that did nothing.
  testthat::expect_false(isTRUE(all.equal(imputed$X, normalized$X)))
  testthat::expect_lt(ncol(clean$X), ncol(normalized$X))
}

test_that("each wizard preview stage is the real chain stopped early", {
  fx <- get_fixture()
  ## Every removal fires: the filter cuts features, the threshold of 1 cuts
  ## samples, and normalization changes every value in between.
  expect_preview_parity(fx$counts, fx$samples, fx$contrasts, list(
    zero_as_na = TRUE, filter_missing = TRUE, filter_threshold = 3,
    impute = TRUE, impute_method = "SVD2",
    remove_outliers = TRUE, outlier_threshold = 1
  ))
})

test_that("pgx.preprocess default output is a plausible log-expression matrix", {
  fx <- get_fixture()
  res <- pgx.preprocess(fx$counts, fx$samples, fx$contrasts)
  expect_equal(dim(res$X), dim(fx$counts))
  expect_true(all(is.finite(res$X)))
  expect_true(max(res$X) < 40) ## log2 scale, not raw counts
})

## ---- degenerate designs: no app reference exists, so assert directly ----
## The app always had a contrast design to reduce; pgx.createPGX() does not, so
## a missing or unusable design is a path the golden reference cannot cover.
## What the grouping itself does to the missingness filter is covered by
## playbase.preprocess; what matters here is that the shim reduces a design the
## same way the app did, and degrades to "one single group" rather than erroring
## when it cannot.

## 8 samples, non-NA counts per feature: r1=8, r2=3, r3=4, r4=2, r5=0, r6=3.
mk_missing <- function() {
  counts <- matrix(100, 6, 8,
    dimnames = list(paste0("r", 1:6), paste0("s", 1:8)))
  counts["r2", 4:8] <- NA
  counts["r3", c(3, 4, 7, 8)] <- NA
  counts["r4", 3:8] <- NA
  counts["r5", ] <- NA
  counts["r6", c(3, 4, 6, 7, 8)] <- NA
  counts
}

filter_opts <- function(threshold) {
  list(filter_missing = TRUE, filter_threshold = threshold, normalize = FALSE)
}

## Two arms over the 8 samples, expressed the way playbase callers hold a
## design: a samples table plus a sample-wise label matrix.
mk_design <- function() {
  ids <- paste0("s", 1:8)
  groups <- rep(c("a", "b"), each = 4)
  list(
    groups = groups,
    samples = data.frame(group = groups, row.names = ids),
    contrasts = matrix(groups, ncol = 1, dimnames = list(ids, "b_vs_a"))
  )
}

## Assert the shim reduces (samples, contrasts) to `groups` under every filter
## rule: >=1 group count, <0 group ratio, and the overall NA ratio in between.
expect_reduces_to <- function(samples, contrasts, groups) {
  counts <- mk_missing()
  for (f in c(3, -1, 0.2)) {
    got <- pgx.preprocess(counts, samples, contrasts, options = filter_opts(f))
    ref <- playbase.preprocess::pgx.preprocess(counts, groups, options = filter_opts(f))
    testthat::expect_equal(got$X, ref$X)
  }
}

test_that("pgx.preprocess reduces a contrast design to sample groups", {
  d <- mk_design()
  expect_reduces_to(d$samples, d$contrasts, d$groups)
  ## and that design really is stricter than none: r3/r6 clear 3 non-NA samples
  ## overall, but not within either arm
  expect_equal(
    rownames(pgx.preprocess(mk_missing(), d$samples, d$contrasts, options = filter_opts(3))$X),
    c("r1", "r2")
  )
})

test_that("pgx.preprocess treats a missing design as one single group", {
  d <- mk_design()
  expect_reduces_to(NULL, NULL, NULL)
  expect_reduces_to(d$samples, NULL, NULL)
  suppressMessages(expect_reduces_to(NULL, d$contrasts, NULL))
  expect_equal(
    rownames(pgx.preprocess(mk_missing(), options = filter_opts(3))$X),
    c("r1", "r2", "r3", "r6")
  )
})

test_that("pgx.preprocess degrades on an unusable contrast matrix", {
  d <- mk_design()
  ## contrasts.convertToLabelMatrix() stop()s on a matrix without column names
  bad <- matrix(d$groups, ncol = 1, dimnames = list(rownames(d$samples), NULL))
  expect_message(
    got <- pgx.preprocess(mk_missing(), d$samples, bad, options = filter_opts(3)),
    "cannot derive sample groups"
  )
  expect_equal(got$X, pgx.preprocess(mk_missing(), options = filter_opts(3))$X)
})
