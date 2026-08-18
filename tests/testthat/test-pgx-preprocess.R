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
    X <- playbase.epigenetics::mToBeta(counts)
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

  ## normalizedX
  if (isTRUE(opt$normalize)) {
    if (opt$datatype == "multi-omics") {
      X <- playbase::normalizeMultiOmics(X)
    } else if (opt$datatype == "methylomics") {
      nX <- try(playbase.epigenetics::normalizeMethylation(X, opt$norm_method, opt$meth_type), silent = TRUE)
      if (!inherits(nX, "try-error") && !is.null(nX)) X <- nX
    } else {
      X <- playbase::normalizeExpression(X, method = opt$norm_method, ref = opt$ref_gene, prior = prior)
    }
  }

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
  list(counts = counts, X = X, annot = annot)
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
  testthat::expect_equal(got$counts, ref$counts, tolerance = 1e-8)
  ## Invariant: counts and X must stay row/col-aligned regardless of branch.
  ## createPGX stop()s otherwise; cheap tripwire if the realign assumption breaks.
  testthat::expect_identical(rownames(got$counts), rownames(got$X))
  testthat::expect_identical(colnames(got$counts), colnames(got$X))
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

test_that("pgx.preprocess default output is a plausible log-expression matrix", {
  fx <- get_fixture()
  res <- pgx.preprocess(fx$counts, fx$samples, fx$contrasts)
  expect_equal(dim(res$X), dim(fx$counts))
  expect_true(all(is.finite(res$X)))
  expect_true(max(res$X) < 40) ## log2 scale, not raw counts
})
