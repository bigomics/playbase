##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

# =============================================================================
# Pathway-enrichment AI report
# =============================================================================

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

# Per-geneset meta effect for one contrast, with a fallback direction.
# meta.fx (the average member-gene fold-change) is NaN for sets whose member
# genes are absent from the matrix -- pervasive in multi-omics (layer-prefixed
# rows miss the bare-symbol GMT) and sparse proteomics -- which would otherwise
# make the report drop genuinely significant sets. Where meta.fx is undefined
# but per-method effects (mx$fc) exist, fall back to their mean so the set keeps
# a direction and is reported.
.pathways_meta_fx <- function(mx, keep) {
  eff <- mx$meta.fx[keep]
  bad <- !is.finite(eff)
  if (any(bad) && !is.null(mx$fc)) {
    fc <- as.matrix(mx$fc)[keep, , drop = FALSE]
    fallback <- rowMeans(fc, na.rm = TRUE)
    eff[bad] <- fallback[bad]
  }
  eff
}

.pathways_contrast_landscape <- function(pgx, fx, q, q_threshold = 0.05) {
  rows <- lapply(colnames(fx), function(contrast) {
    effect <- fx[, contrast]
    qvalue <- q[, contrast]
    significant <- is.finite(effect) & is.finite(qvalue) & qvalue < q_threshold
    data.frame(
      contrast = contrast,
      modeled_sample_context = .ai_report_get(pgx, "contrast_context",
                                              contrast = contrast),
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
    candidates <- which(is.finite(effect) & is.finite(qvalue))
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
      "Modeled sample context: ", .ai_report_get(pgx, "contrast_context",
                                                 contrast = contrast),
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

pathways_build_report_tables <- function(slice, pgx, ntop = 50L,
                                         cross_n = 50L, q_threshold = 0.05) {
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
  fx <- sapply(slice$meta, .pathways_meta_fx, keep = keep)  # meta.fx + fallback
  q <- sapply(slice$meta, function(mx) mx$meta.q[keep])
  if (is.null(dim(fx))) {
    fx <- matrix(fx, ncol = 1L)
    q <- matrix(q, ncol = 1L)
  }
  rownames(fx) <- rownames(q) <- terms[keep]
  colnames(fx) <- colnames(q) <- contrasts

  sample_metadata <- .ai_report_get(pgx, "sample_metadata",
                                    include_missing = FALSE)
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
    experiment_info = .ai_report_get(
      pgx, "experiment_info",
      slice = slice,
      n_contrasts = length(contrasts),
      fallback = "(unnamed)",
      include_missing = FALSE
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
    .ai_report_prompt_path("pathways", "pathways_report_data.md")
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

pathways_build_methods <- function(slice, pgx) {
  build_report_methods("pathways", "pathways_methods.md",
                       params = .methods_params(pgx, slice, "pathways"))
}

ai.pathways.create_report <- function(pgx, slice, ai) {
  contrasts <- names(slice$meta %||% list())
  if (!length(contrasts)) {
    message("[ai.pathways.create_report] no contrasts in pgx$gset.meta$meta -- skipping")
    return(NULL)
  }
  if (!requireNamespace("omicsai", quietly = TRUE)) {
    stop("omicsai package required for AI report generation", call. = FALSE)
  }

  ntop <- min(as.integer(ai$ntop), 50L)
  data_block <- pathways_build_report_tables(slice, pgx, ntop = ntop)$text
  bp <- .ai_report_build_prompt(pgx, "pathways", data_block)
  out <- .ai_report_run_prompt(bp, ai)
  out$report <- paste(out$report, pathways_build_methods(slice, pgx), sep = "\n\n")
  out
}
