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
  # Datatype-aware vocabulary: resolved once and injected into every prompt
  # fragment so modality nouns ({{feature}}, {{quantity}}, {{analysis}}, ...)
  # are filled from pgx$datatype rather than hardcoded transcriptomics wording.
  vocab <- omicsai::omicsai_datatype_vocab(pgx$datatype)
  context_params <- utils::modifyList(vocab, context_vars)
  prompt <- omicsai::report_prompt(
    role = omicsai::frag("system_base", vocab),
    task = omicsai::frag("text/report", vocab),
    species = omicsai::omicsai_species_prompt(pgx$organism),
    context = omicsai::frag(fragments$context, context_params),
    board_rules = omicsai::frag(fragments$rules, vocab),
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
  ## Preserve omicsai's usage list verbatim (a plain list from
  ## .omicsai_extract_usage(), not the S7 object) so downstream telemetry reads
  ## its native field names (total_tokens, input_tokens_fresh/cached,
  ## output_tokens, reasoning_tokens, cache_hit_rate); attach the model used so
  ## a recorded report event can attribute cost to it.
  usage <- res$metadata$usage
  if (!is.null(usage)) usage$model <- ai$llm_model
  list(
    report = .ai_report_normheadings(res$text),
    prompt = paste0("# SYSTEM\n\n", bp$system,
                    "\n\n---\n\n# BOARD\n\n", bp$board),
    usage  = usage
  )
}
