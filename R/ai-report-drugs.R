##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

## ============================================================================
## Drug Connectivity Report Utilities
## ----------------------------------------------------------------------------
## Lightweight extractors and tier classifiers consumed by the drug
## connectivity AI report pipeline.
##
## The functions in this file are intentionally small and pure: they
## operate on pre-computed pgx$drugs[[method]] matrices and return
## data frames or lists, with no Shiny or LLM dependencies.
## ============================================================================

## Internal: split pipe / semicolon / comma separated annotation tokens.
## L1000 / CMap annotations sometimes ship mixed-encoding bytes; iconv with
## sub = "" silently drops invalid sequences so strsplit's regex engine
## can run without "input string is invalid" / "unable to translate" warnings.
.drugs_parseTokens <- function(x) {
  x <- as.character(x %||% "")
  x[is.na(x)] <- ""
  x <- iconv(x, from = "UTF-8", to = "UTF-8", sub = "")
  x[is.na(x)] <- ""
  lapply(x, function(s) trimws(strsplit(s, split = "[\\|;,]")[[1]]))
}


#' @title NES direction semantics for a drug connectivity analysis type
#'
#' @param analysis_type Character. Name of the analysis backing the drug
#'   enrichment result, e.g. \code{"L1000/activity"}, \code{"L1000/gene"},
#'   \code{"CTRPv2/sensitivity"}, \code{"GDSC/sensitivity"}.
#'
#' @return A list with two character scalars:
#'   \describe{
#'     \item{\code{analysis_type}}{the input value (verbatim, or \code{"unknown"}).}
#'     \item{\code{analysis_type_description}}{a paragraph describing the
#'       data source and the direction semantics of NES for that backend.}
#'   }
#'
#' @description Returns a deterministic, human-readable description of a
#'   drug connectivity analysis type, including the meaning of positive vs.
#'   negative NES. L1000-based analyses interpret negative NES as opposition
#'   (reversal candidate); pharmacological sensitivity analyses (CTRPv2,
#'   GDSC) interpret positive NES as predicted drug vulnerability.
#'
#' @details The returned description is intended for inclusion in AI report
#'   prompts so the language model interprets NES sign correctly for the
#'   selected backend. Matching is case-insensitive.
#'
drugs.analysisInfo <- function(analysis_type) {
  key <- tolower(analysis_type %||% "")

  desc <- if (key %in% c("l1000/activity", "l1000_activitys_n20d1011")) {
    paste(
      "L1000 Transcriptional Activity Scores (TAS) connectivity.",
      "Each drug signature is weighted by strength and reproducibility across cell lines.",
      "Negative NES = drug opposes the experimental state (reversal candidate);",
      "positive NES = drug mimics it (shared transcriptional programme)."
    )
  } else if (key == "l1000/gene") {
    paste(
      "L1000 gene-level log-fold-change connectivity (978 landmark genes, per-cell-line).",
      "More granular than TAS but noisier; best used to corroborate TAS findings.",
      "Negative NES = drug opposes the experimental state;",
      "positive NES = drug mimics it."
    )
  } else if (key == "ctrpv2/sensitivity") {
    paste(
      "CTRPv2 pharmacological sensitivity connectivity (AUC dose-response,",
      "~481 compounds, ~860 cancer cell lines).",
      "NES reflects co-variation of drug sensitivity with the experimental signature.",
      "Positive NES = the experimental state predicts drug vulnerability (sensitivity);",
      "negative NES = resistance or insensitivity."
    )
  } else if (key == "gdsc/sensitivity") {
    paste(
      "GDSC pharmacogenomic sensitivity connectivity (IC50-based,",
      "~367 clinical anti-cancer drugs, ~987 cell lines).",
      "Strongest clinical relevance due to approved and late-stage compounds.",
      "Positive NES = the experimental state predicts drug sensitivity;",
      "negative NES = predicted resistance."
    )
  } else {
    paste(
      "Drug connectivity analysis based on the selected pre-computed signature resource.",
      "Interpret results according to the selected analysis type metadata."
    )
  }

  list(
    analysis_type = analysis_type %||% "unknown",
    analysis_type_description = desc
  )
}


#' @title Build a drug enrichment (DSEA) table for one contrast/method
#'
#' @param pgx A pgx object holding pre-computed drug enrichment matrices in
#'   \code{pgx$drugs[[method]]} (X = NES, P = pvalue, Q = qvalue, annot).
#' @param method Character. Name of the entry in \code{pgx$drugs} to read.
#' @param contrast Character. Contrast name (column of \code{pgx$drugs[[method]]$X}).
#' @param annot Optional data.frame of drug annotations with at least the
#'   columns \code{moa} and \code{target} and row names matching drug names.
#'   If \code{NULL}, the function falls back to
#'   \code{pgx$drugs[[method]]$annot}.
#' @param only_annotated Logical. If \code{TRUE}, drop drugs that have empty
#'   MOA and target after the join.
#'
#' @return A data.frame ordered by descending \code{|NES|}, with columns
#'   \code{drug}, \code{NES}, \code{pval}, \code{padj}, \code{moa}, \code{target}.
#'   Returns \code{NULL} if the inputs are missing or the contrast does not
#'   exist.
#'
#' @description Extracts a single-contrast drug enrichment table from
#'   pre-computed \code{pgx$drugs[[method]]} matrices, joining drug-level
#'   MOA and target annotations when available. This is the input
#'   expected by \code{\link{drugs.moaEnrichment}}.
#'
#' @details Compared to \code{\link{pgx.getDrugConnectivityTable}}, this
#'   function accepts an explicit \code{annot} argument so callers can
#'   supply annotations from an external source when
#'   \code{pgx$drugs[[method]]$annot} is missing. NES/pval/padj are
#'   rounded to four digits; missing NES are coerced to 0 and missing
#'   p/q-values to 1.
#'
drugs.enrichmentTable <- function(pgx, method, contrast,
                                  annot = NULL, only_annotated = FALSE) {
  if (is.null(pgx$drugs) || is.null(method) || !nzchar(method)) return(NULL)
  dr <- pgx$drugs[[method]]
  if (is.null(dr) || is.null(contrast) || !contrast %in% colnames(dr$X)) return(NULL)

  nes <- round(dr$X[, contrast], 4)
  pv <- round(dr$P[, contrast], 4)
  qv <- round(dr$Q[, contrast], 4)
  drug <- rownames(dr$X)

  if (is.null(annot)) annot <- dr$annot

  nes[is.na(nes)] <- 0
  qv[is.na(qv)] <- 1
  pv[is.na(pv)] <- 1

  if (is.null(annot)) {
    dt <- data.frame(
      drug = drug, NES = nes, pval = pv, padj = qv,
      moa = "", target = "",
      stringsAsFactors = FALSE
    )
  } else {
    jj <- match(toupper(drug), toupper(rownames(annot)))
    moa_col <- if ("moa" %in% colnames(annot)) "moa" else NA_character_
    target_col <- if ("target" %in% colnames(annot)) "target" else NA_character_
    aa <- data.frame(
      moa = if (!is.na(moa_col)) annot[jj, moa_col] else "",
      target = if (!is.na(target_col)) annot[jj, target_col] else "",
      stringsAsFactors = FALSE
    )
    dt <- data.frame(
      drug = drug, NES = nes, pval = pv, padj = qv, aa,
      stringsAsFactors = FALSE
    )
  }

  dt <- dt[order(-abs(dt$NES)), , drop = FALSE]
  rownames(dt) <- dt$drug

  if (isTRUE(only_annotated)) {
    sel <- which((dt$moa %||% "") != "" | (dt$target %||% "") != "")
    dt <- dt[sel, , drop = FALSE]
  }

  dt
}


#' @title On-the-fly MOA / target enrichment via fgsea
#'
#' @param dsea_table A data.frame produced by \code{\link{drugs.enrichmentTable}},
#'   with at least the columns \code{NES}, \code{moa}, \code{target} and
#'   drug names as row names.
#' @param field One of \code{"moa"} or \code{"target"}. Selects the
#'   annotation column to tokenise into pseudo gene sets.
#' @param nperm Integer. Number of fgsea permutations.
#'
#' @return A data.frame of fgsea results (\code{pathway}, \code{NES},
#'   \code{pval}, \code{padj}, \code{size}, ...) ordered by descending
#'   \code{|NES|}. Returns \code{NULL} if the input is empty or no
#'   tokens are found.
#'
#' @description Computes MOA-level or target-level enrichment by treating
#'   each annotation token as a pseudo gene set over the drug-ranked NES
#'   vector. Unlike \code{\link{pgx.getDrugMOATable}}, which requires
#'   \code{pgx$drugs[[method]]$moa} to be pre-populated, this function
#'   recomputes the enrichment on demand from the drug-level table.
#'
#' @details Tokens are extracted from the selected annotation column by
#'   splitting on \code{|}, \code{;} or \code{,}. Empty / \code{NA} /
#'   \code{"N/A"} tokens are dropped. Each unique token defines a gene set
#'   (the set of drugs annotated with that token), and fgsea ranks them
#'   against the NES vector. Useful for AI report generation where the
#'   report context only has the drug-level table and needs MOA-level signal
#'   on the fly.
#'
drugs.moaEnrichment <- function(dsea_table, field = c("moa", "target"),
                                nperm = 5000) {
  field <- match.arg(field)
  if (is.null(dsea_table) || !is.data.frame(dsea_table) || nrow(dsea_table) == 0) {
    return(NULL)
  }

  terms.list <- .drugs_parseTokens(dsea_table[[field]])
  names(terms.list) <- rownames(dsea_table)
  terms <- setdiff(unique(unlist(terms.list)), c(NA, "", " ", "NA", "N/A"))
  if (length(terms) == 0) return(NULL)

  gmt <- lapply(terms, function(g) {
    names(which(sapply(terms.list, function(t) g %in% t)))
  })
  names(gmt) <- terms

  rnk <- dsea_table$NES
  names(rnk) <- rownames(dsea_table)

  res <- suppressWarnings(tryCatch(
    fgsea::fgsea(gmt, rnk, nperm = nperm),
    error = function(e) NULL
  ))
  if (is.null(res) || !is.data.frame(res) || nrow(res) == 0) return(NULL)

  res[order(-abs(res$NES)), , drop = FALSE]
}


#' @title Classify a result by significance support bucket
#'
#' @param padj Numeric. Adjusted p-value (q-value) of one result.
#' @param pval Numeric. Raw p-value of the same result. Used as a
#'   fallback when \code{padj} is \code{NA} or above 0.05.
#'
#' @return A character scalar: \code{"significant"} when \code{padj < 0.05},
#'   \code{"nominal"} when only \code{pval < 0.05}, otherwise
#'   \code{"unsupported"}.
#'
#' @description Tier classifier used to gate evidence buckets in AI report
#'   generation: significant results drive headline claims, nominal
#'   results are mentioned with hedging, unsupported results are
#'   suppressed.
#'
drugs.supportBucket <- function(padj, pval = NA_real_) {
  if (!is.na(padj) && padj < 0.05) return("significant")
  if (!is.na(pval) && pval < 0.05) return("nominal")
  "unsupported"
}

# =============================================================================
# AI-report builders
# =============================================================================

.drugs_prompt_path <- function(name) {
  .ai_report_prompt_path("drugs", name)
}

dc_percent <- function(x, digits = 1L) {
  sprintf(paste0("%.", digits, "f%%"), 100 * x)
}

## -----------------------------------------------------------------------------
## L1000 annotation fallback - reads the bundled CSV when pgx$drugs[[method]]$annot
## is missing and FILESX is available.
## -----------------------------------------------------------------------------
dc_load_annot_fallback <- function() {
  tryCatch(
    {
      a <- read.csv(
        file.path(FILESX, "cmap/L1000_repurposing_drugs.txt"),
        sep = "\t", comment.char = "#"
      )
      if (!is.null(a) && "pert_iname" %in% colnames(a)) {
        rownames(a) <- a$pert_iname
      }
      a
    },
    error = function(e) NULL
  )
}

## -----------------------------------------------------------------------------
## Direction subset + tier-rank ordering of an MOA / target enrichment table.
## Adds `support` and `support_rank` columns for report context summaries.
## -----------------------------------------------------------------------------
dc_enrichment_direction_subset <- function(dt, direction = c("opposing", "mimicking")) {
  direction <- match.arg(direction)
  if (is.null(dt) || !is.data.frame(dt) || nrow(dt) == 0) return(dt[0, , drop = FALSE])
  sel <- if (direction == "opposing") dt$NES < 0 else dt$NES > 0
  out <- dt[sel, , drop = FALSE]
  if (nrow(out) == 0) return(out)
  out$support <- vapply(seq_len(nrow(out)), function(i) {
    drugs.supportBucket(out$padj[i], out$pval[i])
  }, character(1))
  out$support_rank <- match(out$support, c("significant", "nominal", "unsupported"))
  out[order(out$support_rank, out$padj, out$pval, -abs(out$NES), out$pathway), , drop = FALSE]
}

## -----------------------------------------------------------------------------
## Context builder - extract deterministic data for one contrast/method context.
## -----------------------------------------------------------------------------
extract_drugconnectivity_context_data <- function(pgx, contrast, method,
                                                  dsea_table = NULL,
                                                  moa_class = NULL,
                                                  moa_target = NULL,
                                                  only_annotated = FALSE,
                                                  n_top = 15L) {
  summary_terms <- function(dt, direction, n = 3L) {
    dt <- dc_enrichment_direction_subset(dt, direction = direction)
    if (is.null(dt) || nrow(dt) == 0) {
      empty <- if (is.data.frame(dt)) dt[0, , drop = FALSE] else data.frame()
      return(list(overall = empty, supported = empty,
                  significant = empty, nominal = empty, unsupported = empty))
    }
    significant <- dt[dt$support == "significant", , drop = FALSE]
    nominal <- dt[dt$support == "nominal", , drop = FALSE]
    unsupported <- dt[dt$support == "unsupported", , drop = FALSE]
    supported <- if (nrow(significant)) significant else if (nrow(nominal)) nominal else dt[0, , drop = FALSE]
    list(
      overall = head(dt, n),
      supported = head(supported, n),
      significant = head(significant, n),
      nominal = head(nominal, n),
      unsupported = head(unsupported, n)
    )
  }

  targets_for_moa <- function(moa_terms, direction, n = 3L) {
    if (is.null(moa_target) || !is.data.frame(moa_target) || !nrow(moa_target) ||
        is.null(moa_terms) || !length(moa_terms)) {
      return(data.frame())
    }
    moa_tokens <- .drugs_parseTokens(dsea_table$moa)
    target_tokens <- .drugs_parseTokens(dsea_table$target)
    keep <- vapply(moa_tokens, function(tokens) any(tokens %in% moa_terms), logical(1))
    if (!any(keep)) return(moa_target[0, , drop = FALSE])
    targets <- setdiff(unique(unlist(target_tokens[keep])), c("", " ", NA, "NA", "N/A"))
    if (!length(targets)) return(moa_target[0, , drop = FALSE])
    dt <- dc_enrichment_direction_subset(moa_target, direction = direction)
    if (is.null(dt) || !nrow(dt)) return(dt)
    dt <- dt[dt$pathway %in% targets, , drop = FALSE]
    if (!nrow(dt)) return(dt)
    head(dt[order(dt$support_rank, dt$padj, dt$pval, -abs(dt$NES), dt$pathway), , drop = FALSE], n)
  }

  exemplars_for_moa <- function(moa_terms, direction, n = 3L) {
    if (is.null(moa_terms) || !length(moa_terms)) return(data.frame())
    dt <- dsea_table[!is.na(dsea_table$NES), , drop = FALSE]
    dt <- if (direction == "opposing") dt[dt$NES < 0, , drop = FALSE] else dt[dt$NES > 0, , drop = FALSE]
    if (!nrow(dt)) return(dt)
    keep <- vapply(.drugs_parseTokens(dt$moa), function(tokens) any(tokens %in% moa_terms), logical(1))
    dt <- dt[keep & !(is.na(dt$moa) | dt$moa == ""), , drop = FALSE]
    if (!nrow(dt)) return(dt)
    dt$support_rank <- vapply(seq_len(nrow(dt)), function(i) {
      match(drugs.supportBucket(dt$padj[i], dt$pval[i]),
            c("significant", "nominal", "unsupported"))
    }, integer(1))
    head(dt[order(dt$support_rank, dt$padj, -abs(dt$NES), dt$drug), , drop = FALSE], n)
  }

  if (is.null(dsea_table)) {
    fallback_annot <- if (is.null(pgx$drugs[[method]]$annot)) dc_load_annot_fallback() else NULL
    dsea_table <- drugs.enrichmentTable(
      pgx = pgx, method = method, contrast = contrast,
      annot = fallback_annot, only_annotated = only_annotated
    )
  }

  if (is.null(dsea_table) || !is.data.frame(dsea_table) || nrow(dsea_table) == 0) {
    return(NULL)
  }

  if (is.null(moa_class))  moa_class  <- drugs.moaEnrichment(dsea_table, "moa")
  if (is.null(moa_target)) moa_target <- drugs.moaEnrichment(dsea_table, "target")

  experiment <- .ai_report_get(pgx, "label")
  has_annot <- (dsea_table$moa != "" | dsea_table$target != "")
  has_annot[is.na(has_annot)] <- FALSE

  top_opposing  <- head(dsea_table[order(dsea_table$NES), , drop = FALSE], n_top)
  top_mimicking <- head(dsea_table[order(-dsea_table$NES), , drop = FALSE], n_top)
  top_abs       <- head(dsea_table[order(-abs(dsea_table$NES)), , drop = FALSE], n_top)
  moa_summary <- list(
    opposing = summary_terms(moa_class, "opposing", n = 3L),
    mimicking = summary_terms(moa_class, "mimicking", n = 3L)
  )
  target_summary <- list(
    opposing = list(
      supported = targets_for_moa(moa_summary$opposing$supported$pathway %||% character(0), "opposing"),
      overall = summary_terms(moa_target, "opposing", n = 3L)
    ),
    mimicking = list(
      supported = targets_for_moa(moa_summary$mimicking$supported$pathway %||% character(0), "mimicking"),
      overall = summary_terms(moa_target, "mimicking", n = 3L)
    )
  )
  exemplars <- list(
    opposing = exemplars_for_moa(moa_summary$opposing$supported$pathway %||% character(0), "opposing"),
    mimicking = exemplars_for_moa(moa_summary$mimicking$supported$pathway %||% character(0), "mimicking")
  )
  frac_annotated <- mean(has_annot)
  annotation_confidence <- if (is.na(frac_annotated)) {
    "unknown"
  } else if (frac_annotated < 0.20) {
    "very low"
  } else if (frac_annotated < 0.30) {
    "low"
  } else if (frac_annotated < 0.60) {
    "moderate"
  } else {
    "high"
  }

  caveats <- c(
    "Connectivity reflects L1000 cell-line perturbation signatures, not direct clinical efficacy.",
    "Dose, cell context, and perturbation time in L1000 may differ from study biology."
  )

  list(
    context_id    = paste(contrast, method, sep = "::"),
    contrast      = contrast,
    method        = method,
    experiment    = experiment,
    n_drugs_tested = nrow(dsea_table),
    n_drugs_sig   = sum(dsea_table$padj < 0.05, na.rm = TRUE),
    n_pos         = sum(dsea_table$NES > 0, na.rm = TRUE),
    n_neg         = sum(dsea_table$NES < 0, na.rm = TRUE),
    max_abs_NES   = max(abs(dsea_table$NES), na.rm = TRUE),
    top_opposing  = top_opposing,
    top_mimicking = top_mimicking,
    top_abs       = top_abs,
    moa_class     = moa_class,
    moa_target    = moa_target,
    moa_summary   = moa_summary,
    target_summary = target_summary,
    exemplars     = exemplars,
    reliability = list(
      has_annotations = any(has_annot),
      frac_annotated  = frac_annotated,
      annotation_confidence = annotation_confidence,
      n_sig_moa_classes = if (is.data.frame(moa_class)) sum(moa_class$padj < 0.05, na.rm = TRUE) else 0L,
      n_sig_targets     = if (is.data.frame(moa_target)) sum(moa_target$padj < 0.05, na.rm = TRUE) else 0L
    ),
    caveats = caveats
  )
}

## -----------------------------------------------------------------------------
## Markdown formatters - render context fields into compact tables and prose.
## -----------------------------------------------------------------------------
dc_compact_table <- function(dt, n = 10L, type = c("drug", "moa"), label = "Pathway") {
  type <- match.arg(type)
  empty <- if (type == "drug") "No data available." else "No significant enrichment available."
  if (is.null(dt) || !is.data.frame(dt) || nrow(dt) == 0) return(empty)
  dt <- head(dt, n)
  fmt <- if (type == "drug") {
    data.frame(
      Drug = dt$drug,
      NES = dt$NES,
      `q-value` = dt$padj,
      MOA = ifelse(is.na(dt$moa) | dt$moa == "", "-", dt$moa),
      Target = ifelse(is.na(dt$target) | dt$target == "", "-", dt$target),
      check.names = FALSE, stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      Name = dt$pathway,
      NES = dt$NES,
      `q-value` = dt$padj,
      Size = dt$size,
      check.names = FALSE, stringsAsFactors = FALSE
    )
  }
  if (type == "moa") colnames(fmt)[1] <- label
  paste(omicsai::omicsai_format_mdtable(
    fmt,
    formatters = list(
      NES = function(x) omicsai::omicsai_format_num(x, 3),
      `q-value` = omicsai::omicsai_format_pvalue
    )
  ), collapse = "\n")
}

dc_format_entries <- function(dt, n = 3L, type = c("summary", "exemplar"),
                              include_support = TRUE) {
  type <- match.arg(type)
  if (is.null(dt) || !is.data.frame(dt) || nrow(dt) == 0) return("none")
  dt <- head(dt, n)
  entries <- vapply(seq_len(nrow(dt)), function(i) {
    name <- if (type == "summary") dt$pathway[i] else dt$drug[i]
    parts <- c(name)
    if (type == "exemplar") parts <- c(parts, dt$moa[i])
    parts <- c(parts, paste0("NES=", omicsai::omicsai_format_num(dt$NES[i], 2)))
    if (!is.na(dt$padj[i])) {
      parts <- c(parts, paste0("q=", omicsai::omicsai_format_pvalue(dt$padj[i])))
    }
    if (type == "summary" && isTRUE(include_support) && "support" %in% colnames(dt)) {
      parts <- c(parts, dt$support[i])
    }
    open <- if (type == "summary") " (" else " ["
    close <- if (type == "summary") ")" else "]"
    paste0(parts[1], open, paste(parts[-1], collapse = ", "), close)
  }, character(1))
  paste(entries, collapse = "; ")
}

## Returns a named list of pre-formatted *values* (no labels, no prose).
## The labels live in the inner template (drugconnectivity_contrast_data.md).
dc_evidence_summary_block <- function(ctx) {
  moa_summary    <- ctx$moa_summary
  target_summary <- ctx$target_summary

  list(
    supported_opposing_moa          = dc_format_entries(moa_summary$opposing$supported, n = 3L),
    supported_mimicking_moa         = dc_format_entries(moa_summary$mimicking$supported, n = 3L),
    corroborating_opposing_targets  = dc_format_entries(target_summary$opposing$supported, n = 3L),
    corroborating_mimicking_targets = dc_format_entries(target_summary$mimicking$supported, n = 3L),
    preferred_opposing_exemplars    = dc_format_entries(ctx$exemplars$opposing, n = 3L, type = "exemplar"),
    preferred_mimicking_exemplars   = dc_format_entries(ctx$exemplars$mimicking, n = 3L, type = "exemplar")
  )
}

## -----------------------------------------------------------------------------
## Render a single per-contrast block from the shared inner template
## (drugconnectivity_contrast_data.md). Used by both Summary mode (one block)
## and Report mode (one block per contrast inside the multi-contrast scaffold).
## -----------------------------------------------------------------------------
dc_render_contrast_block <- function(ctx, tier_label, ntop = 10L) {
  if (is.null(ctx)) return("")
  rel <- ctx$reliability
  ev  <- dc_evidence_summary_block(ctx)

  template <- omicsai::omicsai_load_template(
    .drugs_prompt_path("drugconnectivity_contrast_data.md")
  )

  omicsai::omicsai_substitute_template(template, list(
    contrast                        = ctx$contrast,
    tier                            = tier_label,
    n_drugs_tested                  = as.character(ctx$n_drugs_tested),
    n_drugs_sig                     = as.character(ctx$n_drugs_sig),
    n_neg                           = as.character(ctx$n_neg),
    n_pos                           = as.character(ctx$n_pos),
    max_abs_nes                     = omicsai::omicsai_format_num(ctx$max_abs_NES, 2),
    frac_annotated_pct              = dc_percent(rel$frac_annotated),
    n_sig_moa_classes               = as.character(rel$n_sig_moa_classes),
    n_sig_targets                   = as.character(rel$n_sig_targets),
    annotation_confidence           = rel$annotation_confidence,
    supported_opposing_moa          = ev$supported_opposing_moa,
    supported_mimicking_moa         = ev$supported_mimicking_moa,
    corroborating_opposing_targets  = ev$corroborating_opposing_targets,
    corroborating_mimicking_targets = ev$corroborating_mimicking_targets,
    preferred_opposing_exemplars    = ev$preferred_opposing_exemplars,
    preferred_mimicking_exemplars   = ev$preferred_mimicking_exemplars,
    top_opposing_table              = dc_compact_table(ctx$top_opposing, n = ntop),
    top_mimicking_table             = dc_compact_table(ctx$top_mimicking, n = ntop),
    moa_class_table                 = dc_compact_table(ctx$moa_class, n = ntop, type = "moa", label = "MOA Class"),
    moa_target_table                = dc_compact_table(ctx$moa_target, n = ntop, type = "moa", label = "Target")
  ))
}

## -----------------------------------------------------------------------------
## Cross-contrast MOA convergence matrix.
## Rows = top MOA classes (ranked by max |NES| across contrasts).
## Columns = contrasts.  Cells = NES with significance marker (** q<0.05, * p<0.05).
## -----------------------------------------------------------------------------
dc_moa_convergence_matrix <- function(contexts, n_moa = 12L) {
  if (length(contexts) == 0) return(NULL)

  moa_entries <- list()
  for (ct in names(contexts)) {
    mc <- contexts[[ct]]$moa_class
    if (is.null(mc) || nrow(mc) == 0) next
    for (i in seq_len(nrow(mc))) {
      moa <- mc$pathway[i]
      if (!moa %in% names(moa_entries)) moa_entries[[moa]] <- list()
      moa_entries[[moa]][[ct]] <- list(
        NES  = mc$NES[i],
        padj = mc$padj[i],
        pval = mc$pval[i]
      )
    }
  }
  if (length(moa_entries) == 0) return(NULL)

  scores <- vapply(moa_entries, function(entries) {
    sig_n <- sum(vapply(entries, function(e) !is.na(e$padj) && e$padj < 0.05, logical(1)))
    nom_n <- sum(vapply(entries, function(e) {
      !is.na(e$padj) && e$padj >= 0.05 && !is.na(e$pval) && e$pval < 0.05
    }, logical(1)))
    max_abs <- max(abs(vapply(entries, function(e) e$NES %||% 0, numeric(1))), na.rm = TRUE)
    100 * sig_n + 10 * nom_n + max_abs
  }, numeric(1))
  top_moa <- names(sort(scores, decreasing = TRUE))[seq_len(min(n_moa, length(scores)))]
  contrasts <- names(contexts)

  fmt_cell <- function(entry) {
    if (is.null(entry) || is.null(entry$NES) || is.na(entry$NES)) return("-")
    sig <- if (!is.na(entry$padj) && entry$padj < 0.05) "**"
           else if (!is.na(entry$pval) && entry$pval < 0.05) "*"
           else ""
    paste0(omicsai::omicsai_format_num(entry$NES, 2), sig)
  }

  col_names  <- c("MOA Class", contrasts)
  header_row <- paste0("| ", paste(col_names, collapse = " | "), " |")
  sep_row    <- paste0("| ", paste(rep("---", length(col_names)), collapse = " | "), " |")
  data_rows  <- vapply(top_moa, function(moa) {
    cells <- vapply(contrasts, function(ct) fmt_cell(moa_entries[[moa]][[ct]]), character(1))
    paste0("| ", paste(c(moa, cells), collapse = " | "), " |")
  }, character(1))

  paste(c(header_row, sep_row, data_rows), collapse = "\n")
}

## -----------------------------------------------------------------------------
## Orchestrators - assemble template-ready params for Summary and Report modes.
## -----------------------------------------------------------------------------
drugconnectivity_build_summary_params <- function(pgx, contrast, method,
                                                  only_annotated = FALSE,
                                                  ntop = 12L) {
  ctx <- extract_drugconnectivity_context_data(
    pgx = pgx, contrast = contrast, method = method,
    only_annotated = only_annotated, n_top = ntop
  )
  if (is.null(ctx)) return(NULL)

  ranking <- drugconnectivity_rank_contexts(list(default = ctx))
  tier_label <- if (nrow(ranking) > 0) ranking$tier[1] else "unknown"

  at <- drugs.analysisInfo(method)

  list(
    experiment                = ctx$experiment,
    analysis_type             = at$analysis_type,
    analysis_type_description = at$analysis_type_description,
    contrast_detail           = dc_render_contrast_block(ctx, tier_label, ntop = ntop)
  )
}

drugconnectivity_rank_contexts <- function(context_data) {
  if (length(context_data) == 0) return(data.frame())

  rows <- lapply(context_data, function(ctx) {
    rel <- ctx$reliability
    score <-
      0.35 * min(ctx$n_drugs_sig / 25, 1) +
      0.25 * min(ctx$max_abs_NES / 2.5, 1) +
      0.20 * min(rel$n_sig_moa_classes / 10, 1) +
      0.20 * min(rel$n_sig_targets / 10, 1)

    tier <- if (rel$frac_annotated < 0.20 || ctx$n_drugs_sig == 0) {
      "data-limited"
    } else if (score >= 0.70) {
      "strong"
    } else if (score >= 0.45) {
      "moderate"
    } else {
      "weak"
    }

    data.frame(
      contrast = ctx$contrast,
      method = ctx$method,
      score = score,
      tier = tier,
      significant_drugs = ctx$n_drugs_sig,
      max_abs_NES = ctx$max_abs_NES,
      sig_moa = rel$n_sig_moa_classes,
      sig_targets = rel$n_sig_targets,
      annotated_frac = rel$frac_annotated,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  out[order(-out$score), , drop = FALSE]
}

drugconnectivity_build_report_tables <- function(pgx, method,
                                                 only_annotated = FALSE,
                                                 max_contexts = 8L,
                                                 ntop = 10L) {
  dr <- pgx$drugs[[method]]
  if (is.null(dr)) {
    return(list(text = "No pre-computed drug data available.", data = list(), ranking = data.frame()))
  }

  contrasts <- if (!is.null(dr$X)) sort(setdiff(colnames(dr$X), grep("^IA:", colnames(dr$X), value = TRUE))) else character(0)
  if (length(contrasts) == 0) {
    return(list(text = "No contrasts available.", data = list(), ranking = data.frame()))
  }

  contexts <- list()
  for (ct in contrasts) {
    ctx <- extract_drugconnectivity_context_data(
      pgx = pgx, contrast = ct, method = method,
      only_annotated = only_annotated, n_top = ntop
    )
    if (!is.null(ctx)) contexts[[ct]] <- ctx
  }

  if (length(contexts) == 0) {
    return(list(text = "No reportable contexts.", data = list(), ranking = data.frame()))
  }

  ranking <- drugconnectivity_rank_contexts(contexts)
  keep <- head(ranking$contrast, max_contexts)
  contexts <- contexts[keep]

  rank_table <- paste(omicsai::omicsai_format_mdtable(
    head(ranking, max_contexts),
    formatters = list(
      score        = function(x) omicsai::omicsai_format_num(x, 2),
      max_abs_NES  = function(x) omicsai::omicsai_format_num(x, 2),
      annotated_frac = dc_percent
    )
  ), collapse = "\n")

  moa_matrix <- dc_moa_convergence_matrix(contexts, n_moa = 12L) %||% "(no MOA convergence data)"

  at <- drugs.analysisInfo(method)

  contrast_blocks <- vapply(keep, function(ct) {
    tier_row <- ranking[ranking$contrast == ct, , drop = FALSE]
    tier_label <- if (nrow(tier_row) > 0) tier_row$tier[1] else "unknown"
    dc_render_contrast_block(contexts[[ct]], tier_label, ntop = ntop)
  }, character(1))

  header_tmpl <- omicsai::omicsai_load_template(
    .drugs_prompt_path("drugconnectivity_report_data.md")
  )

  text <- omicsai::omicsai_substitute_template(header_tmpl, list(
    experiment               = .ai_report_get(pgx, "label"),
    analysis_type            = at$analysis_type,
    analysis_type_description = at$analysis_type_description,
    rank_table               = rank_table,
    moa_matrix               = moa_matrix,
    contrast_detail          = paste(contrast_blocks, collapse = "\n")
  ))

  list(text = text, data = contexts, ranking = ranking)
}

drugconnectivity_build_methods <- function(pgx, method) {
  at <- drugs.analysisInfo(method)
  params <- list(
    experiment = .ai_report_get(pgx, "label"),
    analysis_type = at$analysis_type,
    analysis_type_description = at$analysis_type_description,
    date = format(Sys.Date(), "%Y-%m-%d")
  )
  build_report_methods("drugs", "drugconnectivity_methods.md", params = params)
}


# =============================================================================
# Section 3 - Prompt assembly (per drug DB / method)
# =============================================================================

#' Assemble drug-connectivity static-report prompt for a single drug DB.
#' @keywords internal
drugs_assemble_prompt <- function(pgx, method, ai) {
  if (!requireNamespace("omicsai", quietly = TRUE)) {
    stop("omicsai package required for AI report generation", call. = FALSE)
  }
  tables <- drugconnectivity_build_report_tables(pgx, method,
                                                 ntop = ai$ntop %||% 10L)
  .ai_report_build_prompt(
    pgx, "drugs", tables$text,
    context_vars = list(experiment = .ai_report_get(pgx, "label"))
  )
}


# =============================================================================
# Section 4 - Entry point
# =============================================================================

#' Generate a drug-connectivity AI report.
#'
#' Loops over the requested drug databases (`ai$drug_dbs`, NULL = all
#' `names(slice)`), calls `drugconnectivity_build_report_tables` per DB,
#' issues one LLM call per DB, and returns one report per DB.
#'
#' @param pgx full pgx object.
#' @param slice pgx$drugs (named list keyed by drug DB).
#' @param ai resolved `ai` list.
#' @keywords internal
.drugs_safe_db_key <- function(db) gsub("[^A-Za-z0-9]+", "_", db)

#' @return Named list keyed by sanitized DB name; each entry is
#'   `list(report, prompt)`. The orchestrator expands each entry into a
#'   top-level `pgx$ai$drugs_<safe_db>` slot - mirrors the wgcna/wgcna_mox
#'   top-level-per-variant pattern. Returns NULL when drugs slot empty.
ai.drugs.create_report <- function(pgx, slice, ai) {
  if (is.null(slice) || length(slice) == 0L) return(NULL)
  if (!requireNamespace("omicsai", quietly = TRUE)) {
    stop("omicsai package required for AI report generation", call. = FALSE)
  }
  dbs <- ai$drug_dbs
  if (is.null(dbs) || length(dbs) == 0L) dbs <- names(slice)
  dbs <- intersect(dbs, names(slice))
  if (length(dbs) == 0L) return(NULL)

  out <- list()
  for (db in dbs) {
    bp  <- drugs_assemble_prompt(pgx, db, ai)
    res <- .ai_report_run_prompt(bp, ai)
    res$report <- paste(res$report, drugconnectivity_build_methods(pgx, db),
                        sep = "\n\n")
    key <- .drugs_safe_db_key(db)
    out[[key]] <- res
  }
  structure(out, class = c("ai_report_multi", "list"))
}
