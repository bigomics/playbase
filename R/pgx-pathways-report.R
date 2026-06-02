##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

# =============================================================================
# Pathway-enrichment AI-report extraction
# =============================================================================
# Owns deterministic PGX -> prompt-data extraction for one integrated pathway
# report across all contrasts. LLM orchestration lives in ai-report.R.


# -----------------------------------------------------------------------------
# Narrow shaping helpers
# -----------------------------------------------------------------------------

.pathways_is_term <- function(term) {
  grepl(
    paste0(
      "^(GOBP|GOMF|GOCC|GO_BP|GO_MF|GO_CC|PATHWAY|PATHWAY_REACTOME|",
      "PATHWAY_WIKI|PATHWAY_BIOPLANET|PATHWAY_NATURE|PATHWAY_PANTHER|",
      "H|HALLMARK):"
    ),
    term
  )
}

.pathways_clean_term <- function(term) {
  term <- sub("^[^:]+:", "", term)
  term <- sub("_Homo sapiens_", " [", term)
  term <- sub("_Mus musculus_", " [", term)
  term <- sub("(R-HSA-[0-9]+|WP[0-9]+)$", "\\1]", term)
  open_bracket <- grepl("\\[[^]]+$", term)
  term[open_bracket] <- paste0(term[open_bracket], "]")
  term
}

.pathways_modeled_context <- function(pgx, contrast) {
  labels <- pgx$contrasts[, contrast]
  keep <- !is.na(labels) & nzchar(labels)
  labels <- labels[keep]
  samples <- pgx$samples[keep, , drop = FALSE]

  variables <- colnames(pgx$samples)
  variables <- variables[!vapply(pgx$samples, is.numeric, logical(1))]
  variables <- variables[vapply(
    pgx$samples[, variables, drop = FALSE],
    function(x) length(unique(x[!is.na(x)])) <= 8L,
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


# -----------------------------------------------------------------------------
# Output-table extractors
# -----------------------------------------------------------------------------

.pathways_sample_metadata <- function(pgx) {
  rows <- lapply(colnames(pgx$samples), function(variable) {
    x <- pgx$samples[[variable]]
    if (is.numeric(x)) {
      present <- x[is.finite(x)]
      summary <- if (length(present)) {
        paste0(
          omicsai::omicsai_format_num(min(present), 3L), " to ",
          omicsai::omicsai_format_num(max(present), 3L)
        )
      } else {
        "no finite values"
      }
    } else {
      tab <- sort(table(x, useNA = "ifany"), decreasing = TRUE)
      values <- paste0(names(tab), "(", as.integer(tab), ")")
      summary <- paste(head(values, 8L), collapse = ", ")
      if (length(values) > 8L) summary <- paste0(summary, ", ...")
    }
    data.frame(variable = variable, summary = summary, stringsAsFactors = FALSE)
  })
  if (!length(rows)) {
    return(data.frame(variable = character(), summary = character()))
  }
  do.call(rbind, rows)
}

.pathways_contrast_landscape <- function(pgx, fx, q, q_threshold = 0.05) {
  rows <- lapply(colnames(fx), function(contrast) {
    effect <- fx[, contrast]
    qvalue <- q[, contrast]
    significant <- is.finite(effect) & is.finite(qvalue) & qvalue < q_threshold
    data.frame(
      contrast = contrast,
      modeled_sample_context = .pathways_modeled_context(pgx, contrast),
      significant_pathways = sum(significant),
      elevated = sum(significant & effect > 0),
      suppressed = sum(significant & effect < 0),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

.pathways_cross_contrast <- function(fx, q, n = 35L, q_threshold = 0.05) {
  sig_count <- rowSums(q < q_threshold, na.rm = TRUE)
  max_effect <- apply(abs(fx), 1L, max, na.rm = TRUE)
  keep <- head(order(sig_count, max_effect, decreasing = TRUE), n)
  values <- apply(
    fx[keep, , drop = FALSE],
    c(1L, 2L),
    omicsai::omicsai_format_num,
    digits = 3L
  )
  data.frame(
    pathway = .pathways_clean_term(rownames(fx)[keep]),
    significant_contrasts = sig_count[keep],
    values,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

.pathways_contrast_evidence <- function(pgx, slice, fx, q, n = 12L,
                                        q_threshold = 0.05) {
  blocks <- vapply(colnames(fx), function(contrast) {
    mx <- slice$meta[[contrast]]
    terms <- rownames(fx)
    effect <- fx[, contrast]
    qvalue <- q[, contrast]
    significant <- is.finite(effect) & is.finite(qvalue) & qvalue < q_threshold
    candidates <- which(if (any(significant)) {
      significant
    } else {
      is.finite(effect) & is.finite(qvalue)
    })
    candidates <- candidates[order(qvalue[candidates], -abs(effect[candidates]))]
    up <- head(candidates[effect[candidates] > 0], n)
    down <- head(candidates[effect[candidates] < 0], n)

    format_evidence <- function(index) {
      if (!length(index)) return("(none)")
      qmat <- as.matrix(mx$q)
      qmat <- qmat[match(terms[index], rownames(mx)), , drop = FALSE]
      table <- data.frame(
        pathway = .pathways_clean_term(terms[index]),
        collection = sub(":.*$", "", terms[index]),
        logFC = omicsai::omicsai_format_num(effect[index], 3L),
        q = omicsai::omicsai_format_pvalue(qvalue[index]),
        significant_methods = paste0(rowSums(qmat < q_threshold, na.rm = TRUE),
                                     "/", ncol(qmat)),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      paste(omicsai::omicsai_format_mdtable(table), collapse = "\n")
    }

    paste0(
      "## ", contrast, "\n\n",
      "Modeled sample context: ", .pathways_modeled_context(pgx, contrast),
      "\n\n",
      if (any(significant)) {
        paste0("Terms below are significant at q < ", q_threshold, ".")
      } else {
        paste0("LOW SIGNAL: no pathway term reaches q < ", q_threshold,
               "; directional trends only.")
      },
      "\n\n### Elevated pathway evidence\n\n", format_evidence(up),
      "\n\n### Suppressed pathway evidence\n\n", format_evidence(down)
    )
  }, character(1))
  paste(blocks, collapse = "\n\n---\n\n")
}


# -----------------------------------------------------------------------------
# Report data and prompt assembly
# -----------------------------------------------------------------------------

#' Build deterministic evidence tables for an integrated pathway AI report.
#'
#' @param slice PGX gset.meta slot.
#' @param pgx Full PGX object.
#' @param ntop Maximum pathway terms per direction and contrast.
#' @param cross_n Maximum rows in the cross-contrast footprint.
#' @param q_threshold Significance threshold for pathway meta q-values.
#' @return List with rendered `text` and structured `data`.
pathways_build_report_tables <- function(slice, pgx, ntop = 12L,
                                         cross_n = 35L, q_threshold = 0.05) {
  if (!requireNamespace("omicsai", quietly = TRUE)) {
    stop("omicsai package required for AI report generation", call. = FALSE)
  }
  contrasts <- names(slice$meta %||% list())
  if (!length(contrasts)) {
    stop("Pathway report requires at least one contrast", call. = FALSE)
  }

  terms <- rownames(slice$meta[[1L]])
  keep <- .pathways_is_term(terms)
  if (!any(keep)) {
    stop("Pathway report requires pathway-oriented terms", call. = FALSE)
  }
  fx <- sapply(slice$meta, function(mx) mx$meta.fx[keep])
  q <- sapply(slice$meta, function(mx) mx$meta.q[keep])
  if (is.null(dim(fx))) {
    fx <- matrix(fx, ncol = 1L)
    q <- matrix(q, ncol = 1L)
  }
  rownames(fx) <- rownames(q) <- terms[keep]
  colnames(fx) <- colnames(q) <- contrasts

  sample_metadata <- .pathways_sample_metadata(pgx)
  contrast_landscape <- .pathways_contrast_landscape(
    pgx, fx, q, q_threshold = q_threshold
  )
  cross_contrast <- .pathways_cross_contrast(
    fx, q, n = cross_n, q_threshold = q_threshold
  )
  contrast_evidence <- .pathways_contrast_evidence(
    pgx, slice, fx, q, n = ntop, q_threshold = q_threshold
  )

  params <- list(
    experiment = pgx$name %||% "(unnamed)",
    description = pgx$description %||% "(not supplied)",
    organism = pgx$organism %||% "unknown",
    datatype = paste(pgx$datatype %||% "unknown", collapse = ", "),
    n_samples = as.character(nrow(pgx$samples)),
    n_contrasts = as.character(length(contrasts)),
    sample_metadata_table = paste(
      omicsai::omicsai_format_mdtable(sample_metadata), collapse = "\n"
    ),
    contrast_landscape_table = paste(
      omicsai::omicsai_format_mdtable(contrast_landscape), collapse = "\n"
    ),
    cross_contrast_table = paste(
      omicsai::omicsai_format_mdtable(cross_contrast), collapse = "\n"
    ),
    contrast_evidence = contrast_evidence
  )

  template <- omicsai::omicsai_load_template(
    omicsai::omicsai_prompt_path("pathways/pathways_report_data.md")
  )
  list(
    text = omicsai::omicsai_substitute_template(template, params),
    data = list(
      sample_metadata = sample_metadata,
      contrast_landscape = contrast_landscape,
      cross_contrast = cross_contrast
    )
  )
}

#' Assemble the integrated pathway-enrichment report prompt.
#'
#' @param slice PGX gset.meta slot.
#' @param pgx Full PGX object.
#' @param ai Resolved AI-report options.
#' @return List with `system` and `board` prompt strings.
#' @keywords internal
pathways_assemble_prompt <- function(slice, pgx, ai) {
  ntop <- min(as.integer(ai$ntop), 12L)
  data_block <- pathways_build_report_tables(slice, pgx, ntop = ntop)$text
  prompt <- omicsai::report_prompt(
    role = omicsai::frag("system_base"),
    task = omicsai::frag("text/report"),
    species = omicsai::omicsai_species_prompt(pgx$organism),
    context = omicsai::frag("pathways/pathways_interpretation"),
    board_rules = omicsai::frag("pathways/pathways_report_rules"),
    data = data_block
  )
  omicsai::build_prompt(prompt)
}
