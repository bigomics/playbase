##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

# =============================================================================
# AI report prompt helpers
# =============================================================================

.ai_report_prompt_path <- function(module, name) {
  if (!requireNamespace("omicsai", quietly = TRUE)) {
    stop("omicsai package required for AI report templates", call. = FALSE)
  }
  omicsai::omicsai_prompt_path(file.path(module, name))
}

.ai_report_prompt_fragments <- function(module) {
  switch(module,
    de = list(
      context = "de/de_interpretation",
      rules = "de/de_report_rules"
    ),
    pathways = list(
      context = "pathways/pathways_interpretation",
      rules = "pathways/pathways_report_rules"
    ),
    mofa = list(
      context = "mofa/mofa_interpretation",
      rules = "mofa/mofa_report_rules"
    ),
    wgcna = list(
      context = "wgcna/wgcna_interpretation",
      rules = "wgcna/wgcna_report_rules"
    ),
    wgcna_mox = list(
      context = "wgcna/wgcna_interpretation",
      rules = "wgcna/wgcna_report_rules"
    ),
    drugs = list(
      context = "drugs/drugconnectivity_interpretation",
      rules = "drugs/drugconnectivity_report_rules"
    ),
    combined = list(
      context = "combined/combined_interpretation",
      rules = "combined/combined_report_rules"
    ),
    stop("Unknown AI report module: ", module, call. = FALSE)
  )
}

.ai_report_build_prompt <- function(pgx, module, data, context_vars = list()) {
  if (!requireNamespace("omicsai", quietly = TRUE)) {
    stop("omicsai package required for AI report generation", call. = FALSE)
  }
  fragments <- .ai_report_prompt_fragments(module)
  context <- if (length(context_vars)) {
    omicsai::frag(fragments$context, context_vars)
  } else {
    omicsai::frag(fragments$context)
  }
  prompt <- omicsai::report_prompt(
    role = omicsai::frag("system_base"),
    task = omicsai::frag("text/report"),
    species = omicsai::omicsai_species_prompt(pgx$organism),
    context = context,
    board_rules = omicsai::frag(fragments$rules),
    data = data
  )
  omicsai::build_prompt(prompt)
}

.ai_report_run_prompt <- function(bp, ai) {
  if (!requireNamespace("omicsai", quietly = TRUE)) {
    stop("omicsai package required for AI report generation", call. = FALSE)
  }
  cfg <- omicsai::omicsai_config(model = ai$llm_model,
                                 system_prompt = bp$system)
  res <- omicsai::omicsai_gen_text(bp$board, config = cfg)
  list(
    report = .ai_report_normheadings(res$text),
    prompt = paste0("# SYSTEM\n\n", bp$system,
                    "\n\n---\n\n# BOARD\n\n", bp$board)
  )
}
