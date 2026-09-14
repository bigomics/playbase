## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
## It owns private downstream composition around the canonical matrix API.

# Returns one explicit preprocessing metadata record for a leaf result.
# The record contains only final state needed by later matrix operations.
# It is not a history, replay log, or compatibility object.
.pgx_preprocess_metadata <- function(result) {
  list(
    alignment = result$alignment,
    space = result$space,
    prior = result$prior,
    layers = .pgx_aligned_layers(result$options, result$alignment),
    options = result$options
  )
}

# Extracts the text before the first feature-name colon.
# Unprefixed values are returned unchanged for downstream grouping.
# This is identifier parsing, not preprocessing dispatch.
.pgx_feature_prefix <- function(names) {
  sub(":.*", "", names)
}

# Creates explicit identity metadata for a caller-supplied processed matrix.
# The matrices must already have identical dimensions and dimnames.
# No relationship is inferred from partial name overlap.
.pgx_identity_preprocess_metadata <- function(counts, X, space, prior) {
  if (
    !identical(dim(counts), dim(X)) || !identical(dimnames(counts), dimnames(X))
  ) {
    stop(
      "[pgx.createPGX] X without preprocessing must be exactly aligned to counts",
      call. = FALSE
    )
  }
  list(
    alignment = list(
      rows = lapply(seq_len(nrow(counts)), as.integer),
      cols = seq_len(ncol(counts))
    ),
    space = space,
    prior = prior,
    layers = .pgx_preprocess_layers(X),
    options = NULL
  )
}

# Subsets final row metadata with positions in the current processed matrix.
# Source-row groups remain anchored to the pristine counts matrix.
# Layer labels follow the same positional subset.
.pgx_subset_preprocess_rows <- function(metadata, keep_rows) {
  metadata$alignment$rows <- metadata$alignment$rows[keep_rows]
  if (!is.null(metadata$layers)) {
    metadata$layers <- metadata$layers[keep_rows]
  }
  metadata
}

# Subsets final column metadata with positions in the current processed matrix.
# Source-column indices remain anchored to the pristine counts matrix.
# Sample metadata is never matched or guessed here.
.pgx_subset_preprocess_cols <- function(metadata, keep_cols) {
  metadata$alignment$cols <- metadata$alignment$cols[keep_cols]
  metadata
}

# Returns one representative source position for every processed row.
# The first source row supplies names for an averaged duplicate group.
# Numerical alignment continues to retain every group member.
.pgx_first_source_rows <- function(metadata) {
  vapply(metadata$alignment$rows, `[[`, integer(1), 1L)
}

# Locates processed rows wholly contained in a retained source-row set.
# Averaged duplicate groups are kept only when every source member survives.
# The returned index follows current processed row order.
.pgx_rows_for_source_rows <- function(metadata, keep_rows) {
  which(vapply(
    metadata$alignment$rows,
    function(group) all(group %in% keep_rows),
    logical(1)
  ))
}

# Composes current row groups into the source-row alignment.
# Each output group retains all source positions of its current rows.
# Layer labels must remain unambiguous across every group.
.pgx_group_preprocess_rows <- function(metadata, row_groups) {
  metadata$alignment$rows <- lapply(row_groups, function(group) {
    unique(as.integer(unlist(
      metadata$alignment$rows[group],
      use.names = FALSE
    )))
  })
  if (!is.null(metadata$layers)) {
    metadata$layers <- vapply(
      row_groups,
      function(group) {
        layer <- unique(metadata$layers[group])
        if (length(layer) != 1L) {
          stop(
            "[preprocess metadata] row group spans multiple layers",
            call. = FALSE
          )
        }
        layer
      },
      character(1)
    )
  }
  metadata
}

# Composes an analysis-row collapse into the source-row alignment.
# Member names must identify current X rows without ambiguity.
# A collapsed row retains the union of every contributing source group.
.pgx_collapse_preprocess_rows <- function(metadata, current_names, members) {
  if (anyDuplicated(current_names)) {
    stop(
      "[preprocess metadata] current row names must be unique",
      call. = FALSE
    )
  }
  member_names <- strsplit(as.character(members), ";", fixed = TRUE)
  current_groups <- lapply(member_names, match, table = current_names)
  if (anyNA(unlist(current_groups, use.names = FALSE))) {
    stop(
      "[preprocess metadata] collapsed rows contain unknown members",
      call. = FALSE
    )
  }
  .pgx_group_preprocess_rows(metadata, current_groups)
}

# Imputes a downstream matrix with the canonical SVD2 family operation.
# Fully prefixed multi-omics rows are processed independently by layer.
# This is analysis preparation and does not mutate pgx preprocessing state.
.pgx_impute_svd2 <- function(X) {
  playbase.preprocess::pp.impute(
    X,
    layers = .pgx_preprocess_layers(X),
    method = "SVD2"
  )
}

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

# Builds the processed pgx matrix on the count scale for model consumers.
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
    # Linear is already count-like; playbase owns this downstream space claim.
    space <- if (metadata$space == "linear") {
      "counts"
    } else {
      metadata$space
    }
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
  space <- .pgx_resolve_layer_values(
    metadata$space,
    metadata$layers,
    "space"
  )
  prior <- .pgx_resolve_layer_values(
    metadata$prior,
    metadata$layers,
    "prior"
  )
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
