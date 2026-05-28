##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

# =============================================================================
# Pathways AI-report module
# =============================================================================
# Phase-2 static path: per-contrast LLM call, concatenated into a single
# pgx$ai$pathways$report. Source data = pgx$gset.meta.
#
# Layout mirrors pgx-de-report.R. The v03 board (board.pathway) carries
# pathway-database filtering (GO / Reactome / WikiPathways) at runtime; at
# pgx-compute time we report on the unified geneset enrichment table (all
# pathway DBs together), which is what `pgx$gset.meta$meta[[ctr]]` exposes.
# =============================================================================

# -----------------------------------------------------------------------------
# Section 1 — Extractor (pgx$gset.meta -> pathway_table data.frame)
# -----------------------------------------------------------------------------

#' Coerce `pgx$gset.meta$meta[[contrast]]` into the pathway_table shape
#' that `pathway_build_ai_params` expects.
#'
#' v03 builds this in Shiny via getReactomeTable / getWikiPathwayTable etc;
#' at pgx-compute we build it directly. Returns NULL when the contrast has
#' no enrichment data.
#'
#' @keywords internal
.pathways_build_table <- function(pgx, contrast) {
  mx <- pgx$gset.meta$meta[[contrast]]
  if (is.null(mx)) return(NULL)
  meta_fx <- mx$meta.fx
  meta_q  <- mx$meta.q
  if (is.null(meta_fx) || length(meta_fx) == 0L) return(NULL)
  data.frame(
    pathway = rownames(mx),
    logFC   = meta_fx,
    meta.q  = meta_q,
    stringsAsFactors = FALSE
  )
}


# -----------------------------------------------------------------------------
# Section 2 — Builder (cp from v03 pathway_build_ai_params, verbatim)
# -----------------------------------------------------------------------------

#' Build AI prompt parameters for pathway enrichment summary.
#'
#' Verbatim port of v03 `pathway_build_ai_params`. Sorts by absolute logFC
#' and formats the top N pathways as a markdown table.
#'
#' @param pgx PGX object.
#' @param contrast Character; selected contrast name.
#' @param pathway_table Data frame with columns `pathway`, `logFC`, `meta.q`.
#' @param pathway_type Character; pathway database type (default "all
#'   enriched pathways" at pgx-compute time, since we don't split by DB).
#' @param ntop Integer; number of top pathways to include (default 20).
#' @return Named list: `contrast`, `pathway_type`, `genesets`, `experiment`.
#' @export
pathway_build_ai_params <- function(pgx,
                                    contrast,
                                    pathway_table,
                                    pathway_type = "all enriched pathways",
                                    ntop = 20) {

  # Build genesets section from pathway table
  genesets <- ""
  if (is.data.frame(pathway_table) && nrow(pathway_table) > 0) {
    # Determine column names based on pathway type
    # GO tables have: id, term, score, logFC, meta.q
    # Reactome/WikiPathway tables have: pathway (or pathway.id + pathway), logFC, meta.q
    has_term <- "term" %in% colnames(pathway_table)
    has_pathway <- "pathway" %in% colnames(pathway_table)
    has_score <- "score" %in% colnames(pathway_table)

    # Get pathway name column
    if (has_term) {
      pathway_names <- pathway_table$term
    } else if (has_pathway) {
      pathway_names <- pathway_table$pathway
    } else {
      pathway_names <- rownames(pathway_table)
    }

    # Clean pathway names (remove HTML links if present)
    pathway_names <- gsub("<[^>]+>", "", pathway_names)
    pathway_names <- trimws(pathway_names)

    # Get logFC values
    logfc <- if ("logFC" %in% colnames(pathway_table)) {
      pathway_table$logFC
    } else {
      rep(NA_real_, nrow(pathway_table))
    }

    # Get q-values
    qval <- if ("meta.q" %in% colnames(pathway_table)) {
      pathway_table$meta.q
    } else {
      rep(NA_real_, nrow(pathway_table))
    }

    # Sort by absolute logFC (most activated/repressed first)
    ord <- order(-abs(logfc))
    pathway_names <- pathway_names[ord]
    logfc <- logfc[ord]
    qval <- qval[ord]

    # Take top N
    n <- min(ntop, length(pathway_names))
    pathway_names <- pathway_names[seq_len(n)]
    logfc <- logfc[seq_len(n)]
    qval <- qval[seq_len(n)]

    # Build markdown table
    if (has_score) {
      score <- pathway_table$score[ord][seq_len(n)]
      gset_table <- paste0(
        "| Pathway | Score | logFC | Q-value |\n",
        "|---------|-------|-------|---------|\n",
        paste(
          sprintf(
            "| %s | %s | %s | %s |",
            pathway_names,
            omicsai::omicsai_format_num(score, 2),
            omicsai::omicsai_format_num(logfc, 2),
            omicsai::omicsai_format_num(qval, 4)
          ),
          collapse = "\n"
        )
      )
    } else {
      gset_table <- paste0(
        "| Pathway | logFC | Q-value |\n",
        "|---------|-------|---------|\n",
        paste(
          sprintf(
            "| %s | %s | %s |",
            pathway_names,
            omicsai::omicsai_format_num(logfc, 2),
            omicsai::omicsai_format_num(qval, 4)
          ),
          collapse = "\n"
        )
      )
    }

    # Summary statistics
    n_up <- sum(logfc > 0, na.rm = TRUE)
    n_down <- sum(logfc < 0, na.rm = TRUE)
    n_sig <- sum(qval < 0.05, na.rm = TRUE)

    genesets <- paste0(
      "**Top Enriched Pathways (by effect size):**\n\n",
      gset_table, "\n\n",
      "**Direction:** ", n_up, " upregulated, ", n_down, " downregulated\n",
      "**Significant (q < 0.05):** ", n_sig, " of ", n, " shown"
    )
  }

  experiment <- pgx$description %||% pgx$name %||% ""

  list(
    contrast = contrast,
    pathway_type = pathway_type,
    genesets = genesets,
    experiment = experiment
  )
}


# -----------------------------------------------------------------------------
# Section 3 — Prompt assembly (per contrast)
# -----------------------------------------------------------------------------

#' Assemble the pathways static-report prompt for a single contrast.
#' @keywords internal
pathways_assemble_prompt <- function(pgx, contrast, ai) {
  if (!requireNamespace("omicsai", quietly = TRUE)) {
    stop("omicsai package required for AI report generation", call. = FALSE)
  }
  tbl    <- .pathways_build_table(pgx, contrast)
  params <- pathway_build_ai_params(pgx, contrast, pathway_table = tbl,
                                    ntop = ai$ntop)

  system_prompt <- paste(
    omicsai::omicsai_instructions("system_base"),
    omicsai::omicsai_species_prompt(pgx$organism),
    sep = "\n\n"
  )

  ctx_tmpl  <- omicsai::omicsai_load_template(
    omicsai::omicsai_prompt_path("pathways/context.md"))
  task_tmpl <- omicsai::omicsai_load_template(
    omicsai::omicsai_prompt_path("pathways/report.md"))

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

#' Generate a pathways AI report (phase-2 static path).
#'
#' Mirrors `de.create_report`: per-contrast LLM call, concatenated.
#'
#' @param slice PGX gset.meta slot (i.e. `pgx$gset.meta`).
#' @param pgx full pgx object.
#' @param ai resolved `ai` list.
#' @return `list(report, prompt)` or NULL when no contrasts.
#' @export
pathways.create_report <- function(slice, pgx, ai) {
  if (!requireNamespace("omicsai", quietly = TRUE)) {
    stop("omicsai package required for AI report generation", call. = FALSE)
  }
  contrasts <- names(slice$meta %||% list())
  if (length(contrasts) == 0L) {
    message("[pathways.create_report] no contrasts in pgx$gset.meta$meta — skipping")
    return(NULL)
  }

  sections <- character(0)
  prompts  <- character(0)
  for (ctr in contrasts) {
    bp  <- pathways_assemble_prompt(pgx, ctr, ai)
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
