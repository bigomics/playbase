##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

# =============================================================================
# AI report common data helpers
# =============================================================================

.ai_report_get <- function(pgx,
                           what = c("label", "info", "sample_metadata",
                                    "design_columns", "contrast_names",
                                    "contrast_groups", "contrast_context",
                                    "contrast_matrix_block",
                                    "experiment_info"),
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
    contrast_matrix_block = .ai_report_get_contrast_matrix_block(pgx),
    experiment_info = .ai_report_get_experiment_info(
      pgx,
      slice = args$slice,
      n_features = args$n_features,
      n_samples = args$n_samples,
      n_contrasts = args$n_contrasts,
      fallback = args$fallback %||% "omics experiment",
      override = args$override,
      include_missing = args$include_missing %||% TRUE,
      max_levels = args$max_levels %||% 8L,
      contrasts_block = args$contrasts_block
    )
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
  if (is.null(cm)) {
    contrasts <- .ai_report_get_contrast_names(pgx)
    if (!length(contrasts)) return("(no contrasts available)")
    return(paste(sprintf("- %s", contrasts), collapse = "\n"))
  }

  group_names <- rownames(cm)
  lines <- vapply(colnames(cm), function(cn) {
    pos <- group_names[cm[, cn] > 0]
    neg <- group_names[cm[, cn] < 0]
    sprintf("- %s: %s vs %s", cn,
            paste(pos, collapse = "+"), paste(neg, collapse = "+"))
  }, character(1))
  paste(lines, collapse = "\n")
}

.ai_report_get_experiment_info <- function(pgx, slice = NULL,
                                           n_features = NULL,
                                           n_samples = NULL,
                                           n_contrasts = NULL,
                                           fallback = "omics experiment",
                                           override = NULL,
                                           include_missing = TRUE,
                                           max_levels = 8L,
                                           contrasts_block = NULL) {
  contrasts <- .ai_report_get_contrast_names(pgx, slice = slice)
  info <- .ai_report_get_info(
    pgx,
    n_features = n_features,
    n_samples = n_samples,
    n_contrasts = n_contrasts %||% length(contrasts),
    fallback = fallback,
    override = override
  )
  sample_metadata <- .ai_report_get_sample_metadata(
    pgx,
    max_levels = max_levels,
    include_missing = include_missing
  )
  sample_metadata_table <- if (nrow(sample_metadata)) {
    paste(omicsai::omicsai_format_mdtable(sample_metadata), collapse = "\n")
  } else {
    "(no sample metadata available)"
  }

  template <- omicsai::omicsai_load_template(
    .ai_report_prompt_path("common", "experiment_info.md")
  )
  omicsai::omicsai_substitute_template(template, c(info, list(
    sample_metadata_table = sample_metadata_table,
    contrasts_block = contrasts_block %||% .ai_report_get_contrast_matrix_block(pgx)
  )))
}

#' Collect deterministic methods provenance from a pgx object.
#'
#' Reads what the pipeline actually did for a report slice, so the Methods block
#' can name the real engines, normalization, and imputation instead of a
#' hardcoded transcriptomics list. The engines are the column names of the
#' slice's first meta result (\code{slice$meta[[1]]$q}); normalization is
#' \code{pgx$norm_method}; imputation is the method recorded on
#' \code{pgx$model.parameters$impute_method} (see \code{compute_testGenes}).
#'
#' @param pgx PGX object.
#' @param slice Module slice (e.g. \code{pgx$gx.meta} or \code{pgx$gset.meta}).
#' @return Named list with \code{engines} (character vector), \code{norm}
#'   (string), and \code{impute} (string).
collect_methods <- function(pgx, slice) {
  meta <- slice$meta %||% list()
  engines <- character(0)
  if (length(meta)) {
    qmat <- meta[[1]]$q
    if (!is.null(qmat)) engines <- colnames(qmat)
  }
  list(
    engines = engines,
    norm = pgx$norm_method %||% "not recorded",
    impute = pgx$model.parameters$impute_method %||% "none (no imputation applied)"
  )
}

# Assemble the parameter block shared by every module's methods template:
# experiment identity, counts, the engines that ran (label + engine-filtered
# bibliography), provenance (norm/impute), and datatype vocabulary. norm/impute
# are harmless extras for templates that omit them (e.g. pathways).
.methods_params <- function(pgx, slice, module) {
  prov <- collect_methods(pgx, slice)
  eng <- omicsai::omicsai_methods_engines(prov$engines, module)
  c(list(
    experiment = .ai_report_get(pgx, "label"),
    date = format(Sys.Date(), "%Y-%m-%d"),
    n_contrasts = length(names(slice$meta %||% list())),
    n_samples = tryCatch(nrow(pgx$samples), error = function(e) NA_integer_),
    n_features = tryCatch(nrow(pgx$X), error = function(e) NA_integer_),
    engines = eng$label,
    norm = prov$norm,
    impute = prov$impute,
    bibliography = eng$bibliography
  ), omicsai::omicsai_datatype_vocab(pgx$datatype))
}

#' Build a deterministic methods appendix for a report.
#'
#' Loads the report-specific methods template, substitutes the provided
#' parameters, and optionally appends footer lines.
#'
#' @param report Report namespace under inst/prompts/.
#' @param template_name Template filename inside the report namespace.
#' @param params Named list of placeholder values.
#' @param footer Optional character vector appended after the template.
#' @return Rendered methods appendix text.
build_report_methods <- function(report, template_name = NULL, params = list(),
                                 footer = NULL) {
  if (!requireNamespace("omicsai", quietly = TRUE)) {
    stop("omicsai package required for AI report generation", call. = FALSE)
  }
  template_name <- template_name %||% paste0(report, "_methods.md")
  template <- omicsai::omicsai_load_template(
    .ai_report_prompt_path(report, template_name)
  )
  text <- omicsai::omicsai_substitute_template(template, params)
  if (!is.null(footer) && length(footer) > 0) {
    text <- omicsai::collapse_lines(text, footer, sep = "\n\n")
  }
  text
}

# =============================================================================
# Report text helpers
# =============================================================================
# The functions above assemble the *input* fed into a report prompt; these
# operate on an already-generated report's *output* text (e.g.
# pgx$ai[[slot]]$report), for callers that want to pull out one section
# rather than show the whole thing (e.g. a dataset summary popup, a board
# preview) instead of duplicating markdown-heading parsing at each call site.

#' Extract one section from an AI report's markdown text.
#'
#' Pulls the text under a given ATX heading (e.g. "## Discussion", any
#' level 1-6), up to (but not including) the next heading at the same or a
#' higher level, or the end of the text. Matching is by heading text alone
#' (case-insensitive, `#`s and surrounding `**bold**`/`*italic*` markup
#' ignored) -- the level isn't part of the match, only used to find where
#' the section ends. Line-based, so a `#`-led line inside a fenced code
#' block is (incorrectly) treated as a heading -- not expected in an
#' LLM-generated prose report, but worth knowing if this is ever pointed at
#' arbitrary markdown.
#'
#' @param report Character scalar (or vector, joined with `"\n"`): the
#'   report's markdown text, e.g. `pgx$ai$combined$report`.
#' @param heading Heading text to match, e.g. `"Discussion"`. A leading
#'   `#`/`##`/etc. is stripped if present, so `"## Discussion"` also works.
#' @param include_heading Include the heading line itself in the returned
#'   text. Default `FALSE` (just the section body).
#' @return Character scalar with the section's text (possibly `""` for an
#'   empty section), or `NULL` if `report` is empty or `heading` isn't
#'   found.
#' @export
ai_report_get_section <- function(report, heading, include_heading = FALSE) {
  if (is.null(report) || !length(report) || !nzchar(report[[1L]])) return(NULL)
  heading <- trimws(sub("^#{1,6}[ \t]*", "", heading %||% ""))
  if (!nzchar(heading)) return(NULL)

  txt <- paste(report, collapse = "\n")
  txt <- gsub("\r\n?", "\n", txt)
  lines <- strsplit(txt, "\n", fixed = TRUE)[[1]]

  is_heading <- function(line) grepl("^#{1,6}[ \t]+", line)
  heading_level <- function(line) nchar(regmatches(line, regexpr("^#+", line)))
  heading_text <- function(line) {
    trimws(gsub("\\*+", "", sub("^#{1,6}[ \t]+", "", line)))
  }

  start_idx <- NA_integer_
  start_level <- NA_integer_
  for (i in seq_along(lines)) {
    if (is_heading(lines[i]) && tolower(heading_text(lines[i])) == tolower(heading)) {
      start_idx <- i
      start_level <- heading_level(lines[i])
      break
    }
  }
  if (is.na(start_idx)) return(NULL)

  end_idx <- length(lines)
  if (start_idx < length(lines)) {
    for (i in seq(start_idx + 1L, length(lines))) {
      if (is_heading(lines[i]) && heading_level(lines[i]) <= start_level) {
        end_idx <- i - 1L
        break
      }
    }
  }

  body_start <- if (isTRUE(include_heading)) start_idx else start_idx + 1L
  if (body_start > end_idx) return("")
  trimws(paste(lines[body_start:end_idx], collapse = "\n"))
}
