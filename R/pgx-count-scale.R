## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
## It owns PGX consumer policy for count-scale matrices and single-cell CPM.
## Generic space conversion remains owned by playbase.preprocess.

# Computes legacy configurable CPM for downstream single-cell utilities.
# The canonical normalization family fixes its pipeline target at one million.
# Non-pipeline callers that require another total keep this arithmetic private.
.pgx_log_cpm <- function(counts, total = 1e6, prior = 1, log = TRUE) {
  if (is.null(total)) {
    total <- min(mean(Matrix::colSums(counts, na.rm = TRUE)), 1e6)
  }
  if (inherits(counts, "dgCMatrix")) {
    out <- counts
    out@x[is.na(out@x)] <- 0
    totals <- Matrix::colSums(out)
    scales <- ifelse(is.finite(totals) & totals != 0, total / totals, 0)
    out@x <- out@x * rep.int(scales, diff(out@p))
    if (log) {
      out@x <- log2(out@x + prior)
    }
    return(out)
  }
  totals <- Matrix::colSums(counts, na.rm = TRUE)
  invalid <- !is.finite(totals) | totals == 0
  divisor <- totals
  divisor[invalid] <- 1
  out <- sweep(as.matrix(counts), 2L, divisor, "/") * total
  out[, invalid] <- 0
  if (log) {
    out <- log2(out + prior)
  }
  dimnames(out) <- dimnames(counts)
  out
}

# Resolves one scalar or named value for every represented layer.
# Named vectors must cover each layer exactly by name.
# The result follows first appearance order in layers.
.pgx_resolve_layer_values <- function(value, layers, field) {
  layer_names <- unique(layers)
  if (length(value) == 1L) {
    return(stats::setNames(rep(value, length(layer_names)), layer_names))
  }
  if (is.null(names(value)) || anyDuplicated(names(value))) {
    stop(
      "[count scale] ",
      field,
      " must be scalar or uniquely named by layer",
      call. = FALSE
    )
  }
  index <- match(layer_names, names(value))
  if (anyNA(index)) {
    stop(
      "[count scale] ",
      field,
      " does not cover every layer",
      call. = FALSE
    )
  }
  stats::setNames(value[index], layer_names)
}

# Builds the processed PGX matrix on the count scale for model consumers.
# Count-compatible layers retain processed values through leaf conversion.
# Beta and M-value layers retain the established aligned-source policy.
.pgx_count_scale_matrix <- function(pgx) {
  metadata <- pgx$settings$preprocess
  if (is.null(metadata)) {
    stop(
      "[count scale] pgx has no explicit preprocessing metadata",
      call. = FALSE
    )
  }
  aligned <- playbase.preprocess::pp.alignCounts(
    pgx$counts,
    metadata$alignment,
    X = pgx$X
  )
  count_spaces <- c("counts", "linear", "log2")
  non_count_spaces <- c("beta", "mvalue")
  if (is.null(metadata$layers)) {
    if (length(metadata$space) != 1L) {
      stop("[count scale] layered spaces require layer metadata", call. = FALSE)
    }
    if (!metadata$space %in% c(count_spaces, non_count_spaces)) {
      stop("[count scale] unknown processed space", call. = FALSE)
    }
    if (metadata$space %in% non_count_spaces) {
      return(aligned)
    }
    space <- if (metadata$space == "linear") "counts" else metadata$space
    return(playbase.preprocess::pp.countScaleMatrix(
      X = pgx$X,
      counts = pgx$counts,
      alignment = metadata$alignment,
      space = space,
      prior = metadata$prior
    ))
  }
  if (length(metadata$layers) != nrow(pgx$X)) {
    stop("[count scale] layers must follow X rows", call. = FALSE)
  }
  space <- .pgx_resolve_layer_values(metadata$space, metadata$layers, "space")
  prior <- .pgx_resolve_layer_values(metadata$prior, metadata$layers, "prior")
  out <- aligned
  for (layer in unique(metadata$layers)) {
    if (!space[[layer]] %in% c(count_spaces, non_count_spaces)) {
      stop("[count scale] unknown processed space", call. = FALSE)
    }
    if (space[[layer]] %in% non_count_spaces) {
      next
    }
    rows <- which(metadata$layers == layer)
    layer_space <- if (identical(space[[layer]], "linear")) {
      "counts"
    } else {
      space[[layer]]
    }
    out[rows, ] <- playbase.preprocess::pp.countScaleMatrix(
      X = pgx$X[rows, , drop = FALSE],
      counts = pgx$counts,
      alignment = list(
        rows = metadata$alignment$rows[rows],
        cols = metadata$alignment$cols
      ),
      space = layer_space,
      prior = prior[[layer]]
    )
  }
  dimnames(out) <- dimnames(pgx$X)
  out
}
