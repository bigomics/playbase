## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
## It owns the pgx-to-matrix preprocessing boundary; numerical stages stay in playbase.preprocess.

#' @importFrom playbase.preprocess pp.alignCounts pp.batchCorrect
#' @importFrom playbase.preprocess pp.convertSpace pp.countScaleMatrix
#' @importFrom playbase.preprocess pp.deduplicate pp.filterFeatures pp.impute
#' @importFrom playbase.preprocess pp.normalize pp.removeOutliers pp.toCountScale
NULL

#' @export
playbase.preprocess::pp.alignCounts

#' @export
playbase.preprocess::pp.batchCorrect

#' @export
playbase.preprocess::pp.convertSpace

#' @export
playbase.preprocess::pp.countScaleMatrix

#' @export
playbase.preprocess::pp.deduplicate

#' @export
playbase.preprocess::pp.filterFeatures

#' @export
playbase.preprocess::pp.impute

#' @export
playbase.preprocess::pp.normalize

#' @export
playbase.preprocess::pp.removeOutliers

#' @export
playbase.preprocess::pp.toCountScale

# Reduces a contrast design to one grouping label per source sample.
# Missing or unusable designs mean that no grouping is available.
# The leaf receives only a plain vector in source-column order.
groupsFromContrasts <- function(samples, contrasts) {
  if (is.null(contrasts)) {
    return(NULL)
  }
  tryCatch(
    apply(
      contrasts.convertToLabelMatrix(contrasts, samples),
      1L,
      paste,
      collapse = "_"
    ),
    error = function(e) {
      message(
        "[pgx.preprocess] cannot derive sample groups from contrasts: ",
        conditionMessage(e)
      )
      NULL
    }
  )
}

# Derives explicit batch columns from playbase sample metadata.
# Only columns named as technical batch variables are selected automatically.
# An absent match stays NULL so the leaf can reject an under-specified method.
batchFromSamples <- function(samples) {
  if (is.null(samples) || !ncol(samples)) {
    return(NULL)
  }
  fields <- grep(
    "batch|replicate|donor|clone|plate|run|lane",
    colnames(samples),
    ignore.case = TRUE,
    value = TRUE
  )
  if (!length(fields)) {
    return(NULL)
  }
  samples[, fields, drop = FALSE]
}

# Returns explicit row-layer labels when every feature has a prefix.
# Partial or absent prefixes are treated as a single matrix.
# This helper supports downstream family calls, not pipeline inference.
.pgx_preprocess_layers <- function(X) {
  names <- rownames(X)
  if (
    is.null(names) || !length(names) || !all(grepl("^[A-Za-z0-9]+:", names))
  ) {
    return(NULL)
  }
  sub(":.*", "", names)
}

# Derives final row layers from source-layer labels and explicit alignment.
# Every merged row must come from one source layer.
# The result follows processed X row order.
.pgx_aligned_layers <- function(options, alignment) {
  layers <- options$layers
  if (is.null(layers)) {
    return(NULL)
  }
  vapply(
    alignment$rows,
    function(group) {
      value <- unique(layers[group])
      if (length(value) != 1L) {
        stop(
          "[pgx.preprocess] aligned row spans multiple layers",
          call. = FALSE
        )
      }
      value
    },
    character(1)
  )
}

#' Preprocess uploaded counts through the matrix pipeline
#'
#' Converts playbase sample and contrast metadata to a plain group vector, then
#' delegates once to [playbase.preprocess::pgx.preprocess()]. All numerical
#' stage order and method dispatch remain in the leaf package.
#'
#' @param counts Numeric source matrix with features in rows and samples in
#'   columns. The returned source matrix is byte-identical.
#' @param samples Optional sample annotation with rows in source-column order.
#' @param contrasts Optional contrast design used to derive missingness groups.
#' @param annot Optional source-row-aligned annotation matrix or data frame.
#' @param options Named canonical preprocessing options. Plain `target` and
#'   `batch` values, when needed, belong in this list.
#'
#' @return The seven-field leaf result: `counts`, `X`, `annot`, `prior`,
#'   `space`, `alignment`, and fully resolved `options`.
#' @export
pgx.preprocess <- function(
  counts,
  samples = NULL,
  contrasts = NULL,
  annot = NULL,
  options = list()
) {
  groups <- groupsFromContrasts(samples, contrasts)
  if (isTRUE(options$batch_correct)) {
    if (is.null(options$target)) {
      options["target"] <- list(groups)
    }
    batch_method <- options$batch_method
    if (is.null(batch_method)) {
      batch_method <- "limma"
    }
    supervised_methods <- c("ComBat", "limma")
    if (is.null(options$batch) && any(batch_method %in% supervised_methods)) {
      options["batch"] <- list(batchFromSamples(samples))
    }
  }
  playbase.preprocess::pgx.preprocess(
    counts = counts,
    groups = groups,
    annot = annot,
    options = options
  )
}
