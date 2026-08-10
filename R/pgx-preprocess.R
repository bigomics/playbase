## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.

#' @title Preprocess uploaded counts into a normalized expression matrix
#'
#' @description
#' Single source of truth for the data preprocessing that turns a raw uploaded
#' counts matrix into the normalized log-expression matrix `X` used downstream.
#' This reproduces, step for step, the pipeline that previously lived in the
#' Omics Playground Shiny upload module
#' (`upload_module_normalization.R`: `imputedX` -> `normalizedX` -> `cleanX`),
#' so that a bare `pgx.createPGX()`/compute-endpoint call yields the same `X`
#' as the interactive app given the same settings.
#'
#' The app keeps its normalization panels for *previewing* how a setting affects
#' the data, but the data that actually enters compute is produced here.
#'
#' @param counts Raw counts matrix (features x samples). Not modified in place.
#' @param samples Samples/phenotype data.frame. Only used for group-wise NA filtering.
#' @param contrasts Contrast matrix. Only used for group-wise NA filtering.
#' @param annot Optional annotation table aligned to `counts` rows; subset alongside NA filtering.
#' @param options Named list of preprocessing settings (the former `input$...` values):
#'   \describe{
#'     \item{datatype}{"RNA-seq" (default), "methylomics" or "multi-omics"; selects the log/normalize branch.}
#'     \item{is_npx}{Data is on Olink NPX / NULISA NPQ log2 scale; back-transform with 2^x. Default FALSE.}
#'     \item{zero_as_na}{Treat zeros as missing (NA). Default FALSE.}
#'     \item{normalize}{Apply the normalization step. Default TRUE.}
#'     \item{norm_method}{Normalization method passed to `normalizeExpression` (CPM, TMM, quantile, ...). Default "CPM".}
#'     \item{ref_gene}{Reference feature(s) for `norm_method = "reference"`. Default NULL.}
#'     \item{filter_missing}{Filter features by missingness. Default FALSE.}
#'     \item{filter_threshold}{Missingness threshold: >=1 group count (>=3), <0 group ratio (>=0.5), else max NA ratio. Default 3.}
#'     \item{impute}{Impute remaining missing values in X. Default FALSE.}
#'     \item{impute_method}{Imputation method passed to `imputeMissing`. Default "SVD2".}
#'     \item{remove_outliers}{Drop outlier samples. Default FALSE.}
#'     \item{outlier_threshold}{z-score cutoff for `detectOutlierSamples`. Default 3.}
#'     \item{outlier_methods}{z-score methods for `detectOutlierSamples`: any of
#'       "z.correlation", "z.distance", "z.features", "z.isoforest". Default the
#'       first three; "z.isoforest" is opt-in as it fits an isolation forest.}
#'     \item{meth_type}{Methylation array type for `normalizeMethylation`. Default NULL.}
#'   }
#'
#' @return A list with `counts` (aligned to X rows/cols), `X` (normalized log-expression),
#'   `annot` (subset alongside X, or NULL) and `prior` (log2 prior used).
#'
#' @export
pgx.preprocess <- function(counts,
                           samples = NULL,
                           contrasts = NULL,
                           annot = NULL,
                           options = list()) {
  opt <- utils::modifyList(
    list(
      datatype = "RNA-seq",
      is_npx = FALSE,
      zero_as_na = FALSE,
      normalize = TRUE,
      norm_method = "CPM",
      ref_gene = NULL,
      filter_missing = FALSE,
      filter_threshold = 3,
      impute = FALSE,
      impute_method = "SVD2",
      remove_outliers = FALSE,
      outlier_threshold = 3,
      ## NB: not NULL. detectOutlierSamples() reads NULL as "all methods",
      ## which would switch on the isoforest behind the caller's back.
      outlier_methods = c("z.correlation", "z.distance", "z.features"),
      meth_type = NULL
    ),
    options
  )

  counts <- as.matrix(counts)
  info("[pgx.preprocess] input: ", nrow(counts), " features x ", ncol(counts),
    " samples (datatype=", opt$datatype, ")")

  ## -------------------------------------------------------------------
  ## imputedX: NA handling + log2(counts + prior)
  ## (upload_module_normalization.R: imputedX reactive)
  ## -------------------------------------------------------------------
  counts[which(is.nan(counts))] <- NA
  counts[which(is.infinite(counts))] <- NA

  ## Olink NPX / NULISA NPQ are on a log2 scale; back-transform to linear.
  if (isTRUE(opt$is_npx)) counts <- 2**counts

  if (any(counts < 0, na.rm = TRUE)) counts <- pmax(counts, 0)

  if (isTRUE(opt$zero_as_na)) counts[which(counts == 0)] <- NA

  is.mox <- playbase::is.multiomics(rownames(counts))
  if (is.mox) {
    X <- counts
    dtypes <- unique(sub(":.*", "", rownames(X)))
    for (dt in dtypes) {
      ii <- grep(paste0("^", dt, ":"), rownames(counts))
      prior <- if (dt != "gx") playbase::getPrior(counts[ii, ]) else 1
      X[ii, ] <- log2(counts[ii, ] + prior)
    }
  } else if (opt$datatype == "methylomics") {
    X <- playbase::mToBeta(counts)
    prior <- 0
  } else {
    prior0 <- playbase::getPrior(counts)
    prior <- ifelse(grepl("CPM|TMM", opt$norm_method), 1, prior0)
    X <- log2(counts + prior)
  }
  info("[pgx.preprocess] NA/Inf handling + log2 transform done (mox=", is.mox,
    "); ", sum(is.na(X)), " NA(s) remain")

  ## Filter features by missingness.
  if (sum(is.na(X)) > 0 && isTRUE(opt$filter_missing)) {
    f <- opt$filter_threshold
    sample.contrasts <- playbase::contrasts.convertToLabelMatrix(contrasts, samples)
    grp <- apply(sample.contrasts, 1, paste, collapse = "_")
    if (f >= 1) {
      grp.sum <- tapply(seq_len(ncol(counts)), grp, function(i) {
        rowSums(!is.na(counts[, i, drop = FALSE]))
      })
      maxsum <- apply(do.call(cbind, grp.sum), 1, max, na.rm = TRUE)
      sel <- (maxsum >= 3)
    } else if (f < 0) {
      grp.avg <- tapply(seq_len(ncol(counts)), grp, function(i) {
        rowMeans(!is.na(counts[, i, drop = FALSE]))
      })
      maxavg <- apply(do.call(cbind, grp.avg), 1, max, na.rm = TRUE)
      sel <- (maxavg >= 0.5)
    } else {
      sel <- (rowMeans(is.na(X)) <= f)
    }
    info("[pgx.preprocess] filter_missing (threshold=", f, "): ",
      sum(!sel), " of ", length(sel), " feature(s) removed")
    X <- X[which(sel), , drop = FALSE]
    counts <- counts[which(sel), , drop = FALSE]
    if (!is.null(annot)) annot <- annot[which(sel), , drop = FALSE]
  }

  ## Impute remaining missing values. Never runs for methylomics.
  if (any(is.na(X)) && isTRUE(opt$impute)) {
    if (is.mox) {
      X <- playbase::imputeMissing.mox(X, method = opt$impute_method)
    } else {
      X <- playbase::imputeMissing(X, method = opt$impute_method)
    }
    info("[pgx.preprocess] imputed missing values (method=", opt$impute_method, ")")
  }

  ## -------------------------------------------------------------------
  ## normalizedX: column-wise normalization
  ## (upload_module_normalization.R: normalizedX reactive)
  ## -------------------------------------------------------------------
  if (isTRUE(opt$normalize)) {
    if (opt$datatype == "multi-omics") {
      X <- playbase::normalizeMultiOmics(X)
    } else if (opt$datatype == "methylomics") {
      nX <- try(playbase::normalizeMethylation(X, opt$norm_method, opt$meth_type), silent = TRUE)
      if (!inherits(nX, "try-error") && !is.null(nX)) X <- nX
    } else {
      if (identical(opt$norm_method, "reference") &&
        (is.null(opt$ref_gene) || !any(nzchar(opt$ref_gene)))) {
        stop("[pgx.preprocess] norm_method = 'reference' requires options$ref_gene to be set")
      }
      X <- playbase::normalizeExpression(X, method = opt$norm_method, ref = opt$ref_gene, prior = prior)
    }
    info("[pgx.preprocess] normalization done (datatype=", opt$datatype,
      ", method=", opt$norm_method, ")")
  }

  ## -------------------------------------------------------------------
  ## cleanX: align and remove outlier samples
  ## (upload_module_normalization.R: cleanX reactive)
  ## -------------------------------------------------------------------
  ## Realign counts to X: normalization may drop and/or reorder rows, so shrink
  ## counts to the features X kept and keep the two row-aligned (createPGX later
  ## asserts rownames(counts) == rownames(X) element-wise). X is always a strict
  ## row-subset of counts (built from it; normalization never adds/renames rows).
  if (!identical(rownames(X), rownames(counts))) {
    ndrop <- sum(!(rownames(counts) %in% rownames(X)))
    if (anyDuplicated(rownames(counts))) {
      ## Duplicate feature names (e.g. RNA-seq gene symbols): column-wise
      ## normalization preserves row order, so a %in% subset keeps counts
      ## aligned to X and, unlike name-indexing, keeps every duplicate row.
      ## ponytail: assumes no normalizer both reorders AND duplicates names
      ## -- true today (dup names only in RNA-seq, which never reorders).
      counts <- counts[rownames(counts) %in% rownames(X), , drop = FALSE]
    } else {
      ## Unique names: index by name so counts follows any reorder AND drop
      ## from normalization (methylation BMIQ reorders to its probe manifest).
      counts <- counts[rownames(X), , drop = FALSE]
    }
    info("[pgx.preprocess] realigned counts to X: ", ndrop,
      " feature(s) dropped by normalization")
  }
  if (isTRUE(opt$remove_outliers)) {
    ## NB: when X has NAs the app returns the SVD2-imputed X (not the original).
    if (sum(is.na(X)) > 0) {
      if (is.mox) {
        X <- playbase::imputeMissing.mox(X, method = "SVD2")
      } else {
        X <- playbase::imputeMissing(X, method = "SVD2")
      }
    }
    res <- playbase::detectOutlierSamples(X, methods = opt$outlier_methods, plot = FALSE)
    is.outlier <- (res$z.outlier > opt$outlier_threshold)
    if (any(is.outlier) && !all(is.outlier)) {
      X <- X[, which(!is.outlier), drop = FALSE]
      counts <- counts[, colnames(X), drop = FALSE]
      info("[pgx.preprocess] removed ", sum(is.outlier), " outlier sample(s)")
    }
  }

  info("[pgx.preprocess] output: ", nrow(X), " features x ", ncol(X), " samples")
  list(counts = counts, X = X, annot = annot, prior = prior)
}
