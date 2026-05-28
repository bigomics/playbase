##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

# =============================================================================
# DE (differential expression) AI-report module
# =============================================================================
# Phase-2 static path: per-contrast LLM call, concatenated into a single
# pgx$ai$de$report. Source data = pgx$gx.meta.
#
# Layout (mirrors pgx-wgcna-report.R):
#   1. Extractors        — none (single PGX accessor in section 2)
#   2. Builders          — expression_build_ai_params (cp'd verbatim from
#                          v03 components/board.expression/R/
#                          expression_ai_summary.R)
#   3. Prompt assembly   — de_assemble_prompt: load v03 context + report
#                          fragments, substitute {{placeholders}}, return
#                          (system, user) pair (no report_prompt — the v03
#                          DE templates are 2-fragment with direct
#                          substitution, simpler than wgcna's 4-fragment
#                          report_prompt class)
#   4. Entry point       — de.create_report(slice, pgx, ai)
# =============================================================================

# -----------------------------------------------------------------------------
# Section 2 — Builder (cp from v03 expression_build_ai_params, verbatim)
# -----------------------------------------------------------------------------

#' Build AI prompt parameters for expression summary.
#'
#' Extracts parameters from differential expression results to populate the
#' prompt template. Verbatim port of v03 expression_build_ai_params.
#'
#' @param pgx PGX object.
#' @param contrast Character; contrast name (e.g., "groupA_vs_groupB").
#' @param ntop Integer; number of top DE genes to include (default 20).
#' @return Named list with template parameters:
#'   `contrast`, `phenotype`, `top_genes`, `summary_stats`, `experiment`.
#' @export
expression_build_ai_params <- function(pgx, contrast, ntop = 20) {

  # Extract experiment description
  experiment <- ""
  if (!is.null(pgx$name) && nzchar(pgx$name)) {
    experiment <- pgx$name
  } else if (!is.null(pgx$description) && nzchar(pgx$description)) {
    experiment <- pgx$description
  }

  # Extract phenotype from contrast name (readable form)
  phenotype <- gsub("_", " ", contrast)

  # Get gene-level meta results for this contrast
  mx <- pgx$gx.meta$meta[[contrast]]

  top_genes <- ""
  summary_stats <- ""

  if (!is.null(mx)) {
    # Compute meta fold change and meta q-value
    meta_fx <- mx$meta.fx
    meta_q <- mx$meta.q
    names(meta_fx) <- rownames(mx)
    names(meta_q) <- rownames(mx)

    # Remove NA / Inf values
    valid <- !is.na(meta_fx) & !is.infinite(meta_fx) & !is.na(meta_q)
    meta_fx <- meta_fx[valid]
    meta_q <- meta_q[valid]

    # Get gene symbols
    symbols <- rep(NA, length(meta_fx))
    if (!is.null(pgx$genes) && "symbol" %in% colnames(pgx$genes)) {
      symbols <- pgx$genes[names(meta_fx), "symbol"]
      symbols <- ifelse(is.na(symbols) | symbols == "", names(meta_fx), symbols)
    } else {
      symbols <- names(meta_fx)
    }

    # Build data frame and sort by absolute fold change
    df <- data.frame(
      probe = names(meta_fx),
      symbol = symbols,
      logFC = meta_fx,
      meta.q = meta_q,
      stringsAsFactors = FALSE
    )
    df <- df[order(-abs(df$logFC)), , drop = FALSE]

    # Deduplicate by symbol, keeping highest absolute FC
    df <- df[!duplicated(df$symbol), , drop = FALSE]

    # Summary statistics (before subsetting to ntop)
    n_total <- nrow(df)
    n_up <- sum(df$logFC > 0 & df$meta.q < 0.05, na.rm = TRUE)
    n_down <- sum(df$logFC < 0 & df$meta.q < 0.05, na.rm = TRUE)
    n_sig <- n_up + n_down

    summary_stats <- paste0(
      "- **Total genes tested:** ", n_total, "\n",
      "- **Significant genes (q < 0.05):** ", n_sig, "\n",
      "  - Upregulated: ", n_up, "\n",
      "  - Downregulated: ", n_down, "\n",
      "- **Top ", min(ntop, nrow(df)), " genes shown** (ranked by |logFC|)"
    )

    # Take top ntop genes
    df <- head(df, ntop)

    if (nrow(df) > 0) {
      gene_table <- paste0(
        "| Gene | logFC | q-value |\n",
        "|------|-------|---------|\n",
        paste(
          sprintf(
            "| %s | %s | %s |",
            df$symbol,
            omicsai::omicsai_format_num(df$logFC, 3),
            omicsai::omicsai_format_num(df$meta.q, 4)
          ),
          collapse = "\n"
        )
      )

      n_up_top <- sum(df$logFC > 0, na.rm = TRUE)
      n_down_top <- sum(df$logFC < 0, na.rm = TRUE)

      top_genes <- paste0(
        "**Top Differentially Expressed Genes (by |logFC|):**\n\n",
        gene_table, "\n\n",
        "**Direction:** ", n_up_top, " upregulated, ", n_down_top,
        " downregulated (of top ", nrow(df), " shown)"
      )
    }
  }

  list(
    contrast = contrast,
    phenotype = phenotype,
    top_genes = top_genes,
    summary_stats = summary_stats,
    experiment = experiment
  )
}


# -----------------------------------------------------------------------------
# Section 3 — Prompt assembly (per contrast)
# -----------------------------------------------------------------------------

#' Assemble the DE static-report prompt for a single contrast.
#'
#' Returns `list(system, user)` where `system` carries the shared role +
#' species context and `user` carries the v03 context + report templates
#' with {{placeholders}} substituted from `expression_build_ai_params`.
#' Phase 3 will feed the same `(system, user)` pair to an
#' `omicsagentovi::Agent` constructor.
#'
#' @keywords internal
de_assemble_prompt <- function(pgx, contrast, ai) {
  if (!requireNamespace("omicsai", quietly = TRUE)) {
    stop("omicsai package required for AI report generation", call. = FALSE)
  }
  params <- expression_build_ai_params(pgx, contrast, ntop = ai$ntop)

  system_prompt <- paste(
    omicsai::omicsai_instructions("system_base"),
    omicsai::omicsai_species_prompt(pgx$organism),
    sep = "\n\n"
  )

  ctx_tmpl  <- omicsai::omicsai_load_template(
    omicsai::omicsai_prompt_path("de/context.md"))
  task_tmpl <- omicsai::omicsai_load_template(
    omicsai::omicsai_prompt_path("de/report.md"))

  user_msg <- paste(
    omicsai::omicsai_substitute_template(ctx_tmpl, params),
    omicsai::omicsai_substitute_template(task_tmpl, params),
    sep = "\n\n"
  )

  list(system = system_prompt, user = user_msg)
}


# -----------------------------------------------------------------------------
# Section 4 — Entry point (called by pgx.update_reports orchestrator)
# -----------------------------------------------------------------------------

#' Generate a DE AI report (phase-2 static path).
#'
#' Loops over every contrast in `slice$meta` (= `pgx$gx.meta$meta`), issues
#' one `omicsai_gen_text` call per contrast with the assembled prompt, and
#' concatenates the responses into a single markdown report. The cost lever
#' is `ai$ntop` (table size per contrast) and the number of contrasts in
#' the PGX.
#'
#' @param slice PGX gx.meta slot (i.e. `pgx$gx.meta`).
#' @param pgx full pgx object.
#' @param ai resolved `ai` list.
#' @return `list(report = <markdown>, prompt = <markdown>)` or NULL when
#'   there are no contrasts in the slice.
#' @export
de.create_report <- function(slice, pgx, ai) {
  if (!requireNamespace("omicsai", quietly = TRUE)) {
    stop("omicsai package required for AI report generation", call. = FALSE)
  }
  contrasts <- names(slice$meta %||% list())
  if (length(contrasts) == 0L) {
    message("[de.create_report] no contrasts in pgx$gx.meta$meta — skipping")
    return(NULL)
  }

  sections <- character(0)
  prompts  <- character(0)
  for (ctr in contrasts) {
    bp  <- de_assemble_prompt(pgx, ctr, ai)
    cfg <- omicsai::omicsai_config(model = ai$llm_model,
                                   system_prompt = bp$system)
    res <- omicsai::omicsai_gen_text(bp$user, config = cfg)
    sections <- c(sections, paste0("## ", ctr, "\n\n", res$text))
    prompts  <- c(prompts,  paste0("## ", ctr, "\n\n",
                                   "# SYSTEM\n\n", bp$system,
                                   "\n\n---\n\n# USER\n\n", bp$user))
  }

  list(
    report = paste(sections, collapse = "\n\n---\n\n"),
    prompt = paste(prompts,  collapse = "\n\n---\n\n")
  )
}
