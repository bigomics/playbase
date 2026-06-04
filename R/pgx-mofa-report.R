##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

# =============================================================================
# MOFA AI-report extraction module
# =============================================================================
# Owns deterministic PGX -> prompt-data extraction for MOFA reports.
# LLM orchestration is routed through ai-report.R and mofa.create_report().

.mofa_prompt_path <- function(name) {
  if (!requireNamespace("omicsai", quietly = TRUE)) {
    stop("omicsai package required for AI report templates", call. = FALSE)
  }
  omicsai::omicsai_prompt_path(paste0("mofa/", name))
}

.mofa_mdtable <- function(df) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return("(none)")
  paste(omicsai::omicsai_format_mdtable(df), collapse = "\n")
}

.mofa_factor_choices <- function(mofa) {
  if (is.null(mofa$W) || is.null(colnames(mofa$W))) character(0) else colnames(mofa$W)
}

.mofa_feature_view <- function(x) sub(":.*", "", x)

.mofa_feature_name <- function(x) sub("^[^:]+:", "", x)

.mofa_experiment_label <- function(pgx, experiment = NULL) {
  experiment %||% pgx$name %||% pgx$description %||% "omics experiment"
}

.mofa_variance_matrix <- function(mofa) {
  if (!is.null(mofa$variance) &&
      (is.matrix(mofa$variance) || is.data.frame(mofa$variance))) {
    v <- as.matrix(mofa$variance)
    if (all(.mofa_factor_choices(mofa) %in% rownames(v))) return(v)
    if (all(.mofa_factor_choices(mofa) %in% colnames(v))) return(t(v))
  }
  if (!is.null(mofa$V) && (is.matrix(mofa$V) || is.data.frame(mofa$V))) {
    v <- as.matrix(mofa$V)
    if (all(.mofa_factor_choices(mofa) %in% colnames(v))) return(t(v))
    if (all(.mofa_factor_choices(mofa) %in% rownames(v))) return(v)
  }
  matrix(numeric(0), nrow = 0L, ncol = 0L)
}

#' Classify MOFA sample metadata columns for report interpretation.
#'
#' @param samples PGX sample metadata table.
#' @param contrasts PGX contrasts table.
#' @return List with `exact_primary`, `biological`, and `nuisance` columns.
mofa_classify_design_columns <- function(samples, contrasts) {
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

.mofa_trait_rows_for_columns <- function(trait_names, columns) {
  unique(unlist(lapply(columns, function(col) {
    grep(paste0("^", gsub("([\\W])", "\\\\\\1", col), "="), trait_names, value = TRUE)
  })))
}

.mofa_trait_matrix <- function(mofa) {
  if (!is.null(mofa$Z) && (is.matrix(mofa$Z) || is.data.frame(mofa$Z))) {
    z <- as.matrix(mofa$Z)
    if (all(.mofa_factor_choices(mofa) %in% colnames(z))) return(z)
    if (all(.mofa_factor_choices(mofa) %in% rownames(z))) return(t(z))
  }
  matrix(numeric(0), nrow = 0L, ncol = 0L)
}

.mofa_factor_pathway_metrics <- function(mofa, factor, q = 0.05) {
  g <- NULL
  if (!is.null(mofa$gsea$table) && factor %in% names(mofa$gsea$table)) {
    g <- mofa$gsea$table[[factor]]
  }
  if (!is.data.frame(g) || !all(c("pathway", "NES", "padj") %in% names(g))) {
    return(list(n_sig = 0L, max_abs_nes = 0))
  }
  padj <- suppressWarnings(as.numeric(g$padj))
  nes <- suppressWarnings(as.numeric(g$NES))
  sig <- !is.na(padj) & padj < q
  list(
    n_sig = sum(sig),
    max_abs_nes = if (any(sig)) max(abs(nes[sig]), na.rm = TRUE) else 0
  )
}

.mofa_max_abs_weight <- function(mofa, factor) {
  if (is.null(mofa$W) || !factor %in% colnames(mofa$W)) return(0)
  max(abs(as.numeric(mofa$W[, factor])), na.rm = TRUE)
}

#' Rank MOFA factors for reporting.
#'
#' @param factor_summary Output from `mofa_extract_factor_summary()`.
#' @return Character vector of factor names.
mofa_rank_factors <- function(factor_summary) {
  if (is.null(factor_summary) || nrow(factor_summary) == 0) return(character(0))
  score <- 2.0 * abs(factor_summary$trait_r) +
    0.8 * factor_summary$total_variance +
    0.15 * pmin(factor_summary$n_sig_pathways, 20) / 20
  factor_summary$factor[order(score, decreasing = TRUE)]
}

#' Extract MOFA experiment design metadata.
#'
#' @param mofa MOFA result object.
#' @param pgx Full PGX object.
#' @return Data frame with one metadata row per field.
mofa_extract_experiment_design <- function(mofa, pgx) {
  design <- mofa_classify_design_columns(pgx$samples, pgx$contrasts)
  views <- if (!is.null(mofa$ww) && is.list(mofa$ww)) names(mofa$ww) else colnames(.mofa_variance_matrix(mofa))
  data.frame(
    Field = c(
      "Experiment", "Organism", "Datatype", "Samples", "Features", "Views",
      "Factors", "Sample metadata", "Contrast names",
      "Exact primary design candidates", "Biological metadata candidates",
      "Nuisance/design candidates"
    ),
    Value = c(
      .mofa_experiment_label(pgx, mofa$experiment),
      pgx$organism %||% "unknown",
      pgx$datatype %||% "unknown",
      as.character(if (!is.null(mofa$F)) nrow(mofa$F) else nrow(pgx$samples)),
      as.character(if (!is.null(mofa$W)) nrow(mofa$W) else nrow(pgx$X)),
      if (length(views)) paste(views, collapse = ", ") else "unknown",
      as.character(length(.mofa_factor_choices(mofa))),
      paste(colnames(pgx$samples %||% data.frame()), collapse = ", "),
      paste(colnames(pgx$contrasts %||% data.frame()), collapse = ", "),
      paste(design$exact_primary, collapse = ", "),
      paste(design$biological, collapse = ", "),
      paste(design$nuisance, collapse = ", ")
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

#' Extract ranked MOFA factor summary.
#'
#' @param mofa MOFA result object.
#' @param pgx Full PGX object.
#' @param n_factors Maximum factors to return.
#' @return Data frame with factor-level summary rows.
mofa_extract_factor_summary <- function(mofa, pgx, n_factors = 6L) {
  factors <- .mofa_factor_choices(mofa)
  if (!length(factors)) return(data.frame())
  v <- .mofa_variance_matrix(mofa)
  z <- .mofa_trait_matrix(mofa)
  design <- mofa_classify_design_columns(pgx$samples, pgx$contrasts)
  trait_names <- rownames(z)
  biological_traits <- .mofa_trait_rows_for_columns(trait_names, design$biological)
  nuisance_traits <- .mofa_trait_rows_for_columns(trait_names, design$nuisance)
  if (!length(biological_traits)) biological_traits <- setdiff(trait_names, nuisance_traits)

  rows <- lapply(factors, function(f) {
    vf <- if (f %in% rownames(v)) v[f, , drop = TRUE] else numeric(0)
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
    pm <- .mofa_factor_pathway_metrics(mofa, f)
    data.frame(
      factor = f,
      total_variance = sum(vf, na.rm = TRUE),
      lead_view = if (length(vf)) names(vf)[which.max(vf)] else "unknown",
      view_class = if (length(vf)) omicsai::omicsai_mofa_view_class(vf) else "view-specific",
      biological_trait = top_trait,
      trait_r = top_r,
      trait_label = omicsai::omicsai_verbalize_r(top_r),
      n_sig_pathways = pm$n_sig,
      max_abs_nes = pm$max_abs_nes,
      max_abs_weight = .mofa_max_abs_weight(mofa, f),
      stringsAsFactors = FALSE
    )
  })
  df <- do.call(rbind, rows)
  ranked <- mofa_rank_factors(df)
  df <- df[match(ranked, df$factor), , drop = FALSE]
  head(df, as.integer(n_factors))
}

#' Extract MOFA per-view variance table.
#'
#' @param mofa MOFA result object.
#' @param factors Character vector of factors to include.
#' @return Long data frame of factor, view, variance, and verbal label.
mofa_extract_view_variance <- function(mofa, factors) {
  v <- .mofa_variance_matrix(mofa)
  keep <- intersect(factors, rownames(v))
  if (!length(keep)) return(data.frame())
  rows <- list()
  for (f in keep) {
    rows[[f]] <- data.frame(
      factor = f,
      view = colnames(v),
      variance = as.numeric(v[f, ]),
      variance_pct = paste0(omicsai::omicsai_format_num(100 * as.numeric(v[f, ]), 1L), "%"),
      variance_label = omicsai::omicsai_verbalize_mofa_variance(as.numeric(v[f, ])),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

#' Extract selected MOFA factor-trait correlations.
#'
#' @param mofa MOFA result object.
#' @param pgx Full PGX object.
#' @param factors Character vector of factors to include.
#' @param top_n Maximum trait rows per factor.
#' @return Data frame with selected factor-trait correlations.
mofa_extract_trait_correlations <- function(mofa, pgx, factors, top_n = 8L) {
  z <- .mofa_trait_matrix(mofa)
  if (!nrow(z) || !ncol(z)) return(data.frame())
  design <- mofa_classify_design_columns(pgx$samples, pgx$contrasts)
  biological_traits <- .mofa_trait_rows_for_columns(rownames(z), design$biological)
  nuisance_traits <- .mofa_trait_rows_for_columns(rownames(z), design$nuisance)
  rows <- list()
  for (f in intersect(factors, colnames(z))) {
    vals <- z[, f]
    ord <- names(head(sort(abs(vals), decreasing = TRUE), top_n))
    keep <- unique(c(intersect(ord, biological_traits), intersect(ord, nuisance_traits), ord))
    keep <- head(keep, top_n)
    rows[[f]] <- data.frame(
      factor = f,
      trait = keep,
      trait_class = ifelse(keep %in% nuisance_traits, "nuisance/design", "biological"),
      r = as.numeric(vals[keep]),
      label = omicsai::omicsai_verbalize_r(as.numeric(vals[keep])),
      stringsAsFactors = FALSE
    )
  }
  if (!length(rows)) data.frame() else do.call(rbind, rows)
}

#' Extract top positive and negative MOFA loading features.
#'
#' @param mofa MOFA result object.
#' @param factors Character vector of factors to include.
#' @param ntop Number of positive and negative features per factor.
#' @param annot Optional feature annotation table.
#' @return Data frame with factor-feature loading evidence.
mofa_extract_factor_features <- function(mofa, factors, ntop = 8L, annot = NULL) {
  if (is.null(mofa$W) || !length(factors)) return(data.frame())
  rows <- list()
  for (f in intersect(factors, colnames(mofa$W))) {
    w <- as.numeric(mofa$W[, f])
    names(w) <- rownames(mofa$W)
    pos <- names(head(sort(w, decreasing = TRUE), ntop))
    neg <- names(head(sort(w, decreasing = FALSE), ntop))
    feats <- c(pos, neg)
    symbols <- .mofa_feature_name(feats)
    if (!is.null(annot)) {
      symbols <- tryCatch(
        playbase::probe2symbol(feats, annot, "symbol", fill_na = TRUE),
        error = function(e) symbols
      )
      symbols <- ifelse(is.na(symbols) | symbols == "", .mofa_feature_name(feats), symbols)
    }
    rows[[f]] <- data.frame(
      factor = f,
      side = rep(c("positive", "negative"), c(length(pos), length(neg))),
      view = .mofa_feature_view(feats),
      feature = feats,
      symbol = symbols,
      weight = as.numeric(w[feats]),
      contribution = omicsai::omicsai_verbalize_mofa_weight(as.numeric(w[feats])),
      stringsAsFactors = FALSE
    )
  }
  if (!length(rows)) data.frame() else do.call(rbind, rows)
}

#' Extract top MOFA factor pathway enrichments.
#'
#' @param mofa MOFA result object.
#' @param factors Character vector of factors to include.
#' @param ntop Number of pathways per factor.
#' @param q Adjusted p-value threshold.
#' @return Data frame with factor pathway evidence.
mofa_extract_factor_pathways <- function(mofa, factors, ntop = 8L, q = 0.05) {
  rows <- list()
  for (f in factors) {
    g <- NULL
    if (!is.null(mofa$gsea$table) && f %in% names(mofa$gsea$table)) {
      g <- mofa$gsea$table[[f]]
    }
    if (!is.data.frame(g) || !all(c("pathway", "NES", "padj") %in% names(g))) next
    g$NES <- suppressWarnings(as.numeric(g$NES))
    g$padj <- suppressWarnings(as.numeric(g$padj))
    g <- g[!is.na(g$padj) & g$padj < q, , drop = FALSE]
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
  if (!length(rows)) data.frame() else do.call(rbind, rows)
}

.mofa_format_factor_summary <- function(df) {
  if (!nrow(df)) return("(no factors available)")
  x <- df[, c("factor", "total_variance", "lead_view", "view_class",
              "biological_trait", "trait_r", "trait_label", "n_sig_pathways")]
  x$total_variance <- paste0(omicsai::omicsai_format_num(100 * x$total_variance, 1L), "%")
  x$trait_r <- omicsai::omicsai_format_num(x$trait_r, 2L)
  .mofa_mdtable(x)
}

.mofa_format_view_variance <- function(df) {
  if (!nrow(df)) return("(variance not available)")
  .mofa_mdtable(df[, c("factor", "view", "variance_pct", "variance_label")])
}

.mofa_format_trait_correlations <- function(df) {
  if (!nrow(df)) return("(trait correlations not available)")
  x <- df[, c("factor", "trait", "trait_class", "r", "label")]
  x$r <- omicsai::omicsai_format_num(x$r, 2L)
  .mofa_mdtable(x)
}

.mofa_format_features <- function(df) {
  if (!nrow(df)) return("(feature loadings not available)")
  x <- df[, c("factor", "side", "view", "symbol", "weight", "contribution")]
  x$weight <- omicsai::omicsai_format_num(x$weight, 2L)
  .mofa_mdtable(x)
}

.mofa_format_pathways <- function(df) {
  if (!nrow(df)) return("(no significant pathway enrichments)")
  x <- df[, c("factor", "pathway", "NES", "direction_strength", "significance")]
  x$NES <- omicsai::omicsai_format_num(x$NES, 2L)
  .mofa_mdtable(x)
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
  design <- mofa_extract_experiment_design(mofa, pgx)
  summary <- mofa_extract_factor_summary(mofa, pgx, n_factors = n_factors)
  factors <- summary$factor %||% character(0)
  variance <- mofa_extract_view_variance(mofa, factors)
  traits <- mofa_extract_trait_correlations(mofa, pgx, factors, top_n = ntop)
  features <- mofa_extract_factor_features(mofa, factors, ntop = ntop, annot = pgx$genes)
  pathways <- mofa_extract_factor_pathways(mofa, factors, ntop = ntop)
  contrasts <- colnames(pgx$contrasts %||% data.frame())
  contrasts_block <- if (length(contrasts)) paste(sprintf("- %s", contrasts), collapse = "\n") else "(no contrasts available)"

  tmpl <- omicsai::omicsai_load_prompt_template("mofa/mofa_report_data.md")
  text <- omicsai::omicsai_substitute_template(tmpl, list(
    experiment_design_table = .mofa_mdtable(design),
    contrasts_block = contrasts_block,
    factor_summary_table = .mofa_format_factor_summary(summary),
    view_variance_table = .mofa_format_view_variance(variance),
    trait_correlations_table = .mofa_format_trait_correlations(traits),
    factor_features_table = .mofa_format_features(features),
    factor_pathways_table = .mofa_format_pathways(pathways)
  ))

  list(
    text = text,
    data = list(
      experiment_design = design,
      factor_summary = summary,
      view_variance = variance,
      trait_correlations = traits,
      factor_features = features,
      factor_pathways = pathways
    )
  )
}

#' Build deterministic methods section for MOFA report.
mofa_build_methods <- function(mofa, pgx) {
  template <- omicsai::omicsai_load_prompt_template("mofa/mofa_methods.md")
  report_date <- format(Sys.Date(), "%Y-%m-%d")
  params <- list(
    experiment = .mofa_experiment_label(pgx, mofa$experiment),
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
  experiment_label <- .mofa_experiment_label(pgx, slice$experiment)
  rp <- omicsai::report_prompt(
    role = omicsai::frag("system_base"),
    task = omicsai::frag("text/report"),
    species = omicsai::omicsai_species_prompt(pgx$organism),
    context = omicsai::frag("mofa/mofa_interpretation",
                            list(experiment = experiment_label)),
    board_rules = omicsai::frag("mofa/mofa_report_rules"),
    data = tables$text
  )
  omicsai::build_prompt(rp)
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
  cfg <- omicsai::omicsai_config(model = ai$llm_model,
                                 system_prompt = bp$system)
  res <- omicsai::omicsai_gen_text(bp$board, config = cfg)
  list(
    report = res$text,
    prompt = paste0("# SYSTEM\n\n", bp$system,
                    "\n\n---\n\n# BOARD\n\n", bp$board)
  )
}
