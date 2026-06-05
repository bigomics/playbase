##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##


##---------------------------------------------------------------------
##-------------------- AI REPORT COMMON DATA --------------------------
##---------------------------------------------------------------------

.ai_report_get <- function(pgx,
                           what = c("label", "info", "sample_metadata",
                                    "design_columns", "contrast_names",
                                    "contrast_groups", "contrast_context",
                                    "contrast_matrix_block"),
                           ...) {
  what <- match.arg(what)
  args <- list(...)

  switch(what,
    label = .ai_report_get_label(
      pgx,
      fallback = args$fallback %||% "omics experiment",
      override = args$override
    ),
    info = .ai_report_get_info(
      pgx,
      n_features = args$n_features,
      n_samples = args$n_samples,
      n_contrasts = args$n_contrasts,
      fallback = args$fallback %||% "omics experiment",
      override = args$override
    ),
    sample_metadata = .ai_report_get_sample_metadata(
      pgx,
      max_levels = args$max_levels %||% 8L,
      include_missing = args$include_missing %||% TRUE
    ),
    design_columns = .ai_report_get_design_columns(
      samples = args$samples %||% pgx$samples,
      contrasts = args$contrasts %||% pgx$contrasts
    ),
    contrast_names = .ai_report_get_contrast_names(
      pgx,
      slice = args$slice
    ),
    contrast_groups = .ai_report_get_contrast_groups(
      pgx,
      contrast = args$contrast
    ),
    contrast_context = .ai_report_get_contrast_context(
      pgx,
      contrast = args$contrast,
      max_levels = args$max_levels %||% 8L
    ),
    contrast_matrix_block = .ai_report_get_contrast_matrix_block(pgx)
  )
}

.ai_report_get_label <- function(pgx, fallback = "omics experiment",
                                 override = NULL) {
  override %||% pgx$name %||% pgx$description %||% fallback
}

.ai_report_get_info <- function(pgx, n_features = NULL, n_samples = NULL,
                                n_contrasts = NULL,
                                fallback = "omics experiment",
                                override = NULL) {
  if (is.null(n_samples)) {
    n_samples <- tryCatch(nrow(pgx$samples), error = function(e) NA_integer_)
  }
  if (is.null(n_features)) {
    n_features <- tryCatch(nrow(pgx$X), error = function(e) NA_integer_)
  }
  if (is.null(n_contrasts)) {
    n_contrasts <- length(.ai_report_get_contrast_names(pgx))
  }

  list(
    experiment = .ai_report_get_label(pgx, fallback = fallback,
                                      override = override),
    description = pgx$description %||% "(not supplied)",
    organism = pgx$organism %||% "unknown",
    datatype = paste(pgx$datatype %||% "unknown", collapse = ", "),
    n_samples = as.character(n_samples),
    n_features = as.character(n_features),
    n_contrasts = as.character(n_contrasts)
  )
}

.ai_report_get_sample_metadata <- function(pgx, max_levels = 8L,
                                           include_missing = TRUE) {
  samples <- pgx$samples %||% data.frame()
  rows <- lapply(colnames(samples), function(variable) {
    x <- samples[[variable]]
    if (is.numeric(x)) {
      present <- x[is.finite(x)]
      summary <- if (length(present)) {
        paste0(
          omicsai::omicsai_format_num(min(present), 3L), " to ",
          omicsai::omicsai_format_num(max(present), 3L),
          if (isTRUE(include_missing) && sum(is.na(x))) {
            paste0("; ", sum(is.na(x)), " missing")
          } else {
            ""
          }
        )
      } else {
        "no finite values"
      }
    } else {
      tab <- sort(table(x, useNA = "ifany"), decreasing = TRUE)
      values <- paste0(names(tab), "(", as.integer(tab), ")")
      summary <- paste(head(values, max_levels), collapse = ", ")
      if (length(values) > max_levels) summary <- paste0(summary, ", ...")
    }
    data.frame(variable = variable, summary = summary, stringsAsFactors = FALSE)
  })
  if (!length(rows)) {
    return(data.frame(variable = character(), summary = character()))
  }
  do.call(rbind, rows)
}

.ai_report_get_design_columns <- function(samples, contrasts) {
  sample_cols <- colnames(samples %||% data.frame())
  contrast_cols <- colnames(contrasts %||% data.frame())
  contrast_vars <- unique(sub(":.*", "", contrast_cols))
  nuisance <- sample_cols[grepl(
    "(^sample$|sample_id|Subject_ID|pt_id|cat_id|^rep$|^Plate$|^Position$|^layer$|batch|run|qc|file|id$)",
    sample_cols,
    ignore.case = TRUE
  )]
  exact_primary <- intersect(sample_cols, contrast_vars)
  biological <- setdiff(sample_cols, nuisance)
  list(
    exact_primary = exact_primary,
    biological = biological,
    nuisance = nuisance
  )
}

.ai_report_get_contrast_names <- function(pgx, slice = NULL) {
  slice_names <- names(slice$meta %||% list())
  if (length(slice_names)) return(slice_names)

  contrast_cols <- colnames(pgx$contrasts %||% data.frame())
  if (length(contrast_cols)) return(contrast_cols)

  cm <- pgx$model.parameters$contr.matrix
  colnames(cm %||% matrix(nrow = 0L, ncol = 0L))
}

.ai_report_get_contrast_groups <- function(pgx, contrast) {
  if (is.null(contrast) || is.null(pgx$contrasts) ||
      !contrast %in% colnames(pgx$contrasts)) {
    return("-")
  }
  values <- pgx$contrasts[, contrast]
  values <- values[!is.na(values) & nzchar(values)]
  tab <- table(values)
  if (!length(tab)) return("-")
  paste0(names(tab), " (n=", as.integer(tab), ")", collapse = " vs ")
}

.ai_report_get_contrast_context <- function(pgx, contrast, max_levels = 8L) {
  if (is.null(contrast) || is.null(pgx$contrasts) ||
      !contrast %in% colnames(pgx$contrasts)) {
    return("-")
  }

  labels <- pgx$contrasts[, contrast]
  keep <- !is.na(labels) & nzchar(labels)
  labels <- labels[keep]
  samples_all <- pgx$samples %||% data.frame()
  samples <- samples_all[keep, , drop = FALSE]

  variables <- colnames(samples_all)
  variables <- variables[!vapply(samples_all[, variables, drop = FALSE],
                                 is.numeric, logical(1))]
  variables <- variables[vapply(
    samples_all[, variables, drop = FALSE],
    function(x) length(unique(x[!is.na(x)])) <= max_levels,
    logical(1)
  )]
  if (!length(variables)) return("-")

  groups <- split(seq_along(labels), labels)
  parts <- vapply(names(groups), function(group) {
    idx <- groups[[group]]
    values <- vapply(variables, function(variable) {
      x <- as.character(samples[[variable]][idx])
      x[is.na(x) | !nzchar(x)] <- "<missing>"
      tab <- sort(table(x), decreasing = TRUE)
      paste0(variable, "=", paste0(names(tab), "(", as.integer(tab), ")",
                                    collapse = ", "))
    }, character(1))
    paste0(group, " (n=", length(idx), "): ", paste(values, collapse = "; "))
  }, character(1))
  paste(parts, collapse = " | ")
}

.ai_report_get_contrast_matrix_block <- function(pgx) {
  cm <- pgx$model.parameters$contr.matrix
  if (is.null(cm)) return("(no contrasts available)")

  group_names <- rownames(cm)
  lines <- vapply(colnames(cm), function(cn) {
    pos <- group_names[cm[, cn] > 0]
    neg <- group_names[cm[, cn] < 0]
    sprintf("- %s: %s vs %s", cn,
            paste(pos, collapse = "+"), paste(neg, collapse = "+"))
  }, character(1))
  paste(lines, collapse = "\n")
}
