## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
## It owns private PGX preprocessing policy, alignment, and lifecycle metadata.
## Numerical preprocessing remains owned by playbase.preprocess.

# Resolves one complete leaf-options list for PGX construction.
# NULL selects the explicit legacy log2 conversion without optional stages.
# The top-level Playbase batch selector remains authoritative.
.pgx_preprocess_options <- function(
  preprocess,
  counts,
  norm_method,
  average.duplicated,
  batch.correct.method
) {
  # Build the legacy-compatible default path
  if (is.null(preprocess)) {
    positive.counts <- counts[is.finite(counts) & counts > 0]
    if (!length(positive.counts)) {
      stop(
        "[pgx.createPGX] counts must contain a finite positive value to derive the legacy prior",
        call. = FALSE
      )
    }
    prior <- ifelse(
      grepl("CPM|TMM|TPM", norm_method),
      1,
      min(positive.counts)
    )
    message("[pgx.createPGX] creating X as log2(counts+p) with p = ", prior)
    preprocess <- list(
      input_space = "counts",
      output_space = "log2",
      is_npx = FALSE,
      zero_as_na = FALSE,
      filter_missing = FALSE,
      impute = FALSE,
      normalize = FALSE,
      normalize_args = list(prior = prior),
      remove_outliers = FALSE,
      dedup = if (average.duplicated) "average" else "unique",
      max_features = NULL
    )
  }

  # Enforce the single Playbase batch selector
  nested.batch.options <- intersect(
    c("batch_correct", "batch_method", "batch.correct.method"),
    names(preprocess)
  )
  if (length(nested.batch.options)) {
    stop(
      "[pgx.createPGX] select batch correction only with batch.correct.method",
      call. = FALSE
    )
  }
  if (is.null(preprocess$dedup)) {
    preprocess$dedup <- if (average.duplicated) "average" else "unique"
  }
  batch.methods <- c(
    "no_batch_correct",
    "ComBat",
    "limma",
    "RUV",
    "SVA",
    "NPM"
  )
  if (
    !is.character(batch.correct.method) ||
      length(batch.correct.method) != 1L ||
      is.na(batch.correct.method) ||
      !batch.correct.method %in% batch.methods
  ) {
    stop(
      "[pgx.createPGX] batch.correct.method must select one supported method",
      call. = FALSE
    )
  }
  preprocess$batch.correct.method <- batch.correct.method
  preprocess
}

# Resolves Playbase batch policy to leaf-ready sample metadata.
# Metadata is restricted to the final analysis sample axis before fitting.
# Explicit batch and target values remain authoritative over autodetection.
.pgx_resolve_batch_metadata <- function(
  preprocess,
  counts,
  samples,
  contrasts,
  analysis.samples,
  batch.pars
) {
  if (identical(preprocess$batch.correct.method, "no_batch_correct")) {
    preprocess[["batch"]] <- NULL
    preprocess$target <- NULL
    preprocess$batch_args <- list()
    return(preprocess)
  }

  analysis.index <- match(analysis.samples, colnames(counts))
  analysis.samples.data <- samples[analysis.samples, , drop = FALSE]
  analysis.contrasts <- contrasts[analysis.samples, , drop = FALSE]

  # Resolve or align batch metadata
  if (is.null(preprocess[["batch"]])) {
    if (!length(batch.pars)) {
      batch.pars <- "<autodetect>"
    }
    if (any(grepl("<autodetect>", batch.pars))) {
      detection.X <- log2(pmax(counts[, analysis.index, drop = FALSE], 0) + 1)
      batch.model <- get_model_parameters(
        detection.X,
        analysis.samples.data,
        pheno = NULL,
        analysis.contrasts
      )
      batch.pars <- batch.model$batch.pars
    }
    if (any(grepl("<none>", batch.pars))) {
      batch.pars <- character()
    }
    batch.pars <- intersect(batch.pars, colnames(analysis.samples.data))
    if (length(batch.pars)) {
      preprocess[["batch"]] <- analysis.samples.data[,
        batch.pars,
        drop = FALSE
      ]
    }
  } else if (
    is.matrix(preprocess[["batch"]]) ||
      is.data.frame(preprocess[["batch"]])
  ) {
    if (
      !is.null(rownames(preprocess[["batch"]])) &&
        all(analysis.samples %in% rownames(preprocess[["batch"]]))
    ) {
      preprocess[["batch"]] <- preprocess[["batch"]][
        analysis.samples,
        ,
        drop = FALSE
      ]
    } else {
      preprocess[["batch"]] <- preprocess[["batch"]][
        analysis.index,
        ,
        drop = FALSE
      ]
    }
  } else if (
    !is.null(names(preprocess[["batch"]])) &&
      all(analysis.samples %in% names(preprocess[["batch"]]))
  ) {
    preprocess[["batch"]] <- preprocess[["batch"]][analysis.samples]
  } else {
    preprocess[["batch"]] <- preprocess[["batch"]][analysis.index]
  }

  # Align optional biological targets
  if (!is.null(preprocess$target)) {
    if (
      !is.null(names(preprocess$target)) &&
        all(analysis.samples %in% names(preprocess$target))
    ) {
      preprocess$target <- preprocess$target[analysis.samples]
    } else {
      preprocess$target <- preprocess$target[analysis.index]
    }
  }
  preprocess
}

# Returns one explicit preprocessing metadata record for a leaf result.
# The record contains only final state needed by later matrix operations.
# It is not a history, replay log, or compatibility object.
.pgx_preprocess_metadata <- function(result) {
  list(
    alignment = result$alignment,
    space = result$space,
    prior = result$prior,
    layers = playbase.preprocess::pp.alignLayers(
      result$options$layers,
      result$alignment
    ),
    options = result$options
  )
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
    layers = playbase.preprocess::pp.inferLayers(X),
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
