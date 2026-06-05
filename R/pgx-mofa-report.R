##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

# =============================================================================
# MOFA AI-report extraction module
# =============================================================================
# Owns deterministic PGX -> prompt-data extraction for MOFA reports.
# LLM orchestration is routed through ai-report.R and mofa.create_report().

.mofa_mdtable <- function(df) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return("(none)")
  paste(omicsai::omicsai_format_mdtable(df), collapse = "\n")
}

.mofa_factor_choices <- function(mofa) {
  if (is.null(mofa$W) || is.null(colnames(mofa$W))) character(0) else colnames(mofa$W)
}

#' Classify MOFA sample metadata columns for report interpretation.
#'
#' @param samples PGX sample metadata table.
#' @param contrasts PGX contrasts table.
#' @return List with `exact_primary`, `biological`, and `nuisance` columns.
mofa_classify_design_columns <- function(samples, contrasts) {
  .ai_report_get(list(samples = samples, contrasts = contrasts),
                 "design_columns")
}

.mofa_matrix <- function(mofa, candidates, factors_as_rows = TRUE) {
  factors <- .mofa_factor_choices(mofa)
  for (slot in candidates) {
    x <- mofa[[slot]]
    if (is.null(x) || !(is.matrix(x) || is.data.frame(x))) next
    x <- as.matrix(x)
    if (factors_as_rows) {
      if (all(factors %in% rownames(x))) return(x)
      if (all(factors %in% colnames(x))) return(t(x))
    } else {
      if (all(factors %in% colnames(x))) return(x)
      if (all(factors %in% rownames(x))) return(t(x))
    }
  }
  matrix(numeric(0), nrow = 0L, ncol = 0L)
}

.mofa_trait_rows <- function(trait_names, columns) {
  unique(unlist(lapply(columns, function(col) {
    grep(paste0("^", gsub("([\\W])", "\\\\\\1", col), "="), trait_names, value = TRUE)
  })))
}

#' Rank MOFA factors for reporting.
#'
#' @param factor_summary MOFA factor summary table.
#' @return Character vector of factor names.
mofa_rank_factors <- function(factor_summary) {
  if (is.null(factor_summary) || nrow(factor_summary) == 0) return(character(0))
  score <- 2.0 * abs(factor_summary$trait_r) +
    0.8 * factor_summary$total_variance +
    0.15 * pmin(factor_summary$n_sig_pathways, 20) / 20
  factor_summary$factor[order(score, decreasing = TRUE)]
}

.mofa_report_data <- function(mofa, pgx, n_factors = 6L, ntop = 8L) {
  factors <- .mofa_factor_choices(mofa)
  if (!length(factors)) {
    empty <- data.frame()
    return(list(
      factor_summary = empty, view_variance = empty,
      trait_correlations = empty, factor_features = empty,
      factor_pathways = empty
    ))
  }

  v <- .mofa_matrix(mofa, c("variance", "V"))
  z <- .mofa_matrix(mofa, "Z", factors_as_rows = FALSE)
  design <- mofa_classify_design_columns(pgx$samples, pgx$contrasts)
  trait_names <- rownames(z)
  biological_traits <- .mofa_trait_rows(trait_names, design$biological)
  nuisance_traits <- .mofa_trait_rows(trait_names, design$nuisance)
  if (!length(biological_traits)) biological_traits <- setdiff(trait_names, nuisance_traits)

  summary <- lapply(factors, function(f) {
    vf <- if (f %in% rownames(v)) v[f, , drop = TRUE] else numeric(0)
    if (length(vf) && is.null(names(vf))) names(vf) <- colnames(v)
    zf <- if (f %in% colnames(z)) z[, f, drop = TRUE] else numeric(0)
    zb <- zf[intersect(names(zf), biological_traits)]
    if (!length(zb)) zb <- zf
    top_trait <- "none"
    top_r <- NA_real_
    if (length(zb)) {
      i <- which.max(abs(zb))
      top_trait <- names(zb)[i]
      top_r <- as.numeric(zb[i])
    }
    g <- if (!is.null(mofa$gsea$table) && f %in% names(mofa$gsea$table)) {
      mofa$gsea$table[[f]]
    } else NULL
    if (is.data.frame(g) && all(c("NES", "padj") %in% names(g))) {
      padj <- suppressWarnings(as.numeric(g$padj))
      nes <- suppressWarnings(as.numeric(g$NES))
      sig <- !is.na(padj) & padj < 0.05
      n_sig <- sum(sig)
      max_abs_nes <- if (any(sig)) max(abs(nes[sig]), na.rm = TRUE) else 0
    } else {
      n_sig <- 0L
      max_abs_nes <- 0
    }
    data.frame(
      factor = f,
      total_variance = sum(vf, na.rm = TRUE),
      lead_view = if (length(vf)) names(vf)[which.max(vf)] else "unknown",
      view_class = if (length(vf)) omicsai::omicsai_mofa_view_class(vf) else "view-specific",
      biological_trait = top_trait,
      trait_r = top_r,
      trait_label = omicsai::omicsai_verbalize_r(top_r),
      n_sig_pathways = n_sig,
      max_abs_nes = max_abs_nes,
      max_abs_weight = if (is.null(mofa$W) || !f %in% colnames(mofa$W)) {
        0
      } else {
        max(abs(as.numeric(mofa$W[, f])), na.rm = TRUE)
      },
      stringsAsFactors = FALSE
    )
  })
  summary <- do.call(rbind, summary)
  summary <- summary[match(mofa_rank_factors(summary), summary$factor), , drop = FALSE]
  summary <- head(summary, as.integer(n_factors))
  factors <- summary$factor %||% character(0)

  keep <- intersect(factors, rownames(v))
  view_variance <- if (length(keep)) {
    do.call(rbind, lapply(keep, function(f) data.frame(
      factor = f,
      view = colnames(v),
      variance = as.numeric(v[f, ]),
      variance_pct = paste0(omicsai::omicsai_format_num(100 * as.numeric(v[f, ]), 1L), "%"),
      variance_label = omicsai::omicsai_verbalize_mofa_variance(as.numeric(v[f, ])),
      stringsAsFactors = FALSE
    )))
  } else data.frame()

  rows <- list()
  for (f in intersect(factors, colnames(z))) {
    vals <- z[, f]
    ord <- names(head(sort(abs(vals), decreasing = TRUE), ntop))
    keep <- unique(c(intersect(ord, biological_traits), intersect(ord, nuisance_traits), ord))
    keep <- head(keep, ntop)
    rows[[f]] <- data.frame(
      factor = f,
      trait = keep,
      trait_class = ifelse(keep %in% nuisance_traits, "nuisance/design", "biological"),
      r = as.numeric(vals[keep]),
      label = omicsai::omicsai_verbalize_r(as.numeric(vals[keep])),
      stringsAsFactors = FALSE
    )
  }
  trait_correlations <- if (!length(rows)) data.frame() else do.call(rbind, rows)

  rows <- list()
  if (!is.null(mofa$W) && length(factors)) {
    for (f in intersect(factors, colnames(mofa$W))) {
      w <- as.numeric(mofa$W[, f])
      names(w) <- rownames(mofa$W)
      pos <- names(head(sort(w, decreasing = TRUE), ntop))
      neg <- names(head(sort(w, decreasing = FALSE), ntop))
      feats <- c(pos, neg)
      symbols <- sub("^[^:]+:", "", feats)
      if (!is.null(pgx$genes)) {
        symbols <- tryCatch(
          playbase::probe2symbol(feats, pgx$genes, "symbol", fill_na = TRUE),
          error = function(e) symbols
        )
        symbols <- ifelse(is.na(symbols) | symbols == "", sub("^[^:]+:", "", feats), symbols)
      }
      rows[[f]] <- data.frame(
        factor = f,
        side = rep(c("positive", "negative"), c(length(pos), length(neg))),
        view = sub(":.*", "", feats),
        feature = feats,
        symbol = symbols,
        weight = as.numeric(w[feats]),
        contribution = omicsai::omicsai_verbalize_mofa_weight(as.numeric(w[feats])),
        stringsAsFactors = FALSE
      )
    }
  }
  factor_features <- if (!length(rows)) data.frame() else do.call(rbind, rows)

  rows <- list()
  for (f in factors) {
    g <- if (!is.null(mofa$gsea$table) && f %in% names(mofa$gsea$table)) {
      g <- mofa$gsea$table[[f]]
    } else NULL
    if (!is.data.frame(g) || !all(c("pathway", "NES", "padj") %in% names(g))) next
    g$NES <- suppressWarnings(as.numeric(g$NES))
    g$padj <- suppressWarnings(as.numeric(g$padj))
    g <- g[!is.na(g$padj) & g$padj < 0.05, , drop = FALSE]
    if (!nrow(g)) next
    g <- head(g[order(abs(g$NES), decreasing = TRUE), , drop = FALSE], ntop)
    rows[[f]] <- data.frame(
      factor = f,
      pathway = omicsai::omicsai_clean_pathway_label(g$pathway),
      NES = g$NES,
      direction_strength = omicsai::omicsai_verbalize_mofa_nes(g$NES),
      significance = omicsai::omicsai_verbalize_q(g$padj),
      stringsAsFactors = FALSE
    )
  }
  factor_pathways <- if (!length(rows)) data.frame() else do.call(rbind, rows)

  list(
    factor_summary = summary,
    view_variance = view_variance,
    trait_correlations = trait_correlations,
    factor_features = factor_features,
    factor_pathways = factor_pathways
  )
}

.mofa_report_tables <- function(data) {
  summary <- data$factor_summary
  if (nrow(summary)) {
    summary <- summary[, c("factor", "total_variance", "lead_view", "view_class",
                           "biological_trait", "trait_r", "trait_label", "n_sig_pathways")]
    summary$total_variance <- paste0(
      omicsai::omicsai_format_num(100 * summary$total_variance, 1L), "%"
    )
    summary$trait_r <- omicsai::omicsai_format_num(summary$trait_r, 2L)
  }

  traits <- data$trait_correlations
  if (nrow(traits)) {
    traits <- traits[, c("factor", "trait", "trait_class", "r", "label")]
    traits$r <- omicsai::omicsai_format_num(traits$r, 2L)
  }

  features <- data$factor_features
  if (nrow(features)) {
    features <- features[, c("factor", "side", "view", "symbol", "weight", "contribution")]
    features$weight <- omicsai::omicsai_format_num(features$weight, 2L)
  }

  pathways <- data$factor_pathways
  if (nrow(pathways)) {
    pathways <- pathways[, c("factor", "pathway", "NES", "direction_strength", "significance")]
    pathways$NES <- omicsai::omicsai_format_num(pathways$NES, 2L)
  }

  list(
    factor_summary_table = if (nrow(summary)) .mofa_mdtable(summary) else "(no factors available)",
    view_variance_table = if (nrow(data$view_variance)) {
      .mofa_mdtable(data$view_variance[, c("factor", "view", "variance_pct", "variance_label")])
    } else "(variance not available)",
    trait_correlations_table = if (nrow(traits)) .mofa_mdtable(traits) else "(trait correlations not available)",
    factor_features_table = if (nrow(features)) .mofa_mdtable(features) else "(feature loadings not available)",
    factor_pathways_table = if (nrow(pathways)) .mofa_mdtable(pathways) else "(no significant pathway enrichments)"
  )
}

#' Build deterministic evidence tables for a MOFA AI report.
#'
#' @param mofa MOFA result object.
#' @param pgx Full PGX object.
#' @param n_factors Number of factors selected for report detail.
#' @param ntop Number of feature/pathway/trait rows per factor.
#' @return List with rendered `text` and raw `data` tables.
mofa_build_report_tables <- function(mofa, pgx, n_factors = 6L, ntop = 8L) {
  if (!requireNamespace("omicsai", quietly = TRUE)) {
    stop("omicsai package required for AI report generation", call. = FALSE)
  }
  data <- .mofa_report_data(mofa, pgx, n_factors = n_factors, ntop = ntop)
  tables <- .mofa_report_tables(data)

  tmpl <- omicsai::omicsai_load_template(.ai_report_prompt_path("mofa", "mofa_report_data.md"))
  text <- omicsai::omicsai_substitute_template(tmpl, c(list(
    experiment_info = .ai_report_get(
      pgx, "experiment_info",
      override = mofa$experiment,
      n_samples = if (!is.null(mofa$F)) nrow(mofa$F) else NULL,
      n_features = if (!is.null(mofa$W)) nrow(mofa$W) else NULL
    )
  ), tables))

  list(text = text, data = data)
}

#' Build deterministic methods section for MOFA report.
mofa_build_methods <- function(mofa, pgx) {
  template <- omicsai::omicsai_load_template(
    .ai_report_prompt_path("mofa", "mofa_methods.md")
  )
  report_date <- format(Sys.Date(), "%Y-%m-%d")
  params <- list(
    experiment = .ai_report_get(pgx, "label", override = mofa$experiment),
    n_factors = length(.mofa_factor_choices(mofa)),
    n_samples = if (!is.null(mofa$F)) nrow(mofa$F) else NA_integer_,
    n_features = if (!is.null(mofa$W)) nrow(mofa$W) else NA_integer_,
    date = report_date
  )
  omicsai::collapse_lines(
    omicsai::omicsai_substitute_template(template, params),
    sprintf("_This report was generated with OmicsPlayground (BigOmics, %s)._", report_date),
    "_Note: AI-generated interpretation may contain inaccuracies and must be independently verified._",
    sep = "\n\n"
  )
}

#' Assemble the MOFA report prompt.
#' @keywords internal
mofa_assemble_prompt <- function(slice, pgx, ai) {
  if (!requireNamespace("omicsai", quietly = TRUE)) {
    stop("omicsai package required for AI report generation", call. = FALSE)
  }
  n_factors <- ai$n_factors %||% 6L
  ntop <- min(as.integer(ai$ntop %||% 8L), 12L)
  tables <- mofa_build_report_tables(slice, pgx, n_factors = n_factors, ntop = ntop)
  .ai_report_build_prompt(
    pgx, "mofa", tables$text,
    context_vars = list(
      experiment = .ai_report_get(pgx, "label", override = slice$experiment)
    )
  )
}

#' Generate a MOFA AI report.
#'
#' @param slice PGX `mofa` slot.
#' @param pgx Full PGX object.
#' @param ai Resolved AI-report options.
#' @return List with `report` and `prompt`, or NULL when MOFA is absent.
#' @export
mofa.create_report <- function(slice, pgx, ai) {
  if (!requireNamespace("omicsai", quietly = TRUE)) {
    stop("omicsai package required for AI report generation", call. = FALSE)
  }
  if (is.null(slice)) return(NULL)
  bp <- mofa_assemble_prompt(slice, pgx, ai)
  .ai_report_run_prompt(bp, ai)
}
