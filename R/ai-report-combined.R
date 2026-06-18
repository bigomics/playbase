##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

# =============================================================================
# Combined AI report
# =============================================================================

combined_report_slots <- function(ai_slot) {
  if (is.null(ai_slot) || !is.list(ai_slot)) return(character(0))

  slots <- setdiff(names(ai_slot), c("meta", "combined"))
  slots <- slots[vapply(slots, function(slot) {
    x <- ai_slot[[slot]]
    is.list(x) && is.character(x$report) && nzchar(x$report[1])
  }, logical(1))]
  preferred <- c("de", "pathways", "wgcna", "wgcna_mox", "mofa")
  drug_slots <- sort(grep("^drugs_", slots, value = TRUE))
  ordered <- c(intersect(preferred, slots), drug_slots)
  c(ordered, setdiff(slots, ordered))
}

.combined_report_source_info <- function(slot) {
  if (identical(slot, "de")) {
    return(c(title = "Differential expression",
             role = "differential_expression"))
  }
  if (identical(slot, "pathways")) {
    return(c(title = "Pathway enrichment",
             role = "pathway_enrichment"))
  }
  if (identical(slot, "wgcna") || identical(slot, "wgcna_mox")) {
    return(c(title = "Co-expression network",
             role = "coexpression_network"))
  }
  if (identical(slot, "mofa")) {
    return(c(title = "Latent-factor structure",
             role = "latent_factor_model"))
  }
  if (grepl("^drugs_", slot)) {
    return(c(title = "Drug connectivity",
             role = "drug_connectivity"))
  }
  c(title = slot, role = "other")
}

combined_build_report_data <- function(ai_slot, pgx) {
  if (!requireNamespace("omicsai", quietly = TRUE)) {
    stop("omicsai package required for AI report generation", call. = FALSE)
  }

  slots <- combined_report_slots(ai_slot)
  if (!length(slots)) {
    stop("Combined report requires at least one pgx$ai report slot",
         call. = FALSE)
  }

  module_reports <- vapply(slots, function(slot) {
    info <- .combined_report_source_info(slot)
    paste0(
      "## Source report: ", info[["title"]], "\n\n",
      "SOURCE MODULE: ", slot, "\n",
      "SOURCE ROLE: ", info[["role"]], "\n\n",
      ai_slot[[slot]]$report
    )
  }, character(1))

  template <- omicsai::omicsai_load_template(
    .ai_report_prompt_path("combined", "combined_report_data.md")
  )
  text <- omicsai::omicsai_substitute_template(template, list(
    experiment_info = .ai_report_get(
      pgx,
      "experiment_info",
      n_contrasts = length(.ai_report_get_contrast_names(pgx)),
      fallback = "(unnamed)"
    ),
    module_reports = paste(module_reports, collapse = "\n\n---\n\n")
  ))
  list(text = text, slots = slots)
}

ai.combined.create_report <- function(pgx, slice, ai) {
  slots <- combined_report_slots(slice)
  if (!length(slots)) return(NULL)
  if (!requireNamespace("omicsai", quietly = TRUE)) {
    stop("omicsai package required for AI report generation", call. = FALSE)
  }

  data_block <- combined_build_report_data(slice, pgx)$text
  bp <- .ai_report_build_prompt(pgx, "combined", data_block)
  out <- .ai_report_run_prompt(bp, ai)
  # Normalize model slips like "**Discussion**" or "## Discussion".
  out$report <- gsub(
    "(?m)^\\s*(?:##\\s+)?(?:\\*\\*)?(Discussion|Conclusion)(?:\\*\\*)?\\s*$",
    "## **\\1**",
    out$report,
    perl = TRUE
  )

  method_headings <- c(
    "Differential expression",
    "Pathway enrichment",
    "Co-expression network",
    "Latent-factor structure",
    "Drug connectivity"
  )
  allowed_headings <- unique(vapply(slots, function(slot) {
    .combined_report_source_info(slot)[["title"]]
  }, character(1)))
  for (heading in setdiff(method_headings, allowed_headings)) {
    out$report <- gsub(
      paste0("(?ms)^###\\s+", heading, "\\s*\\n.*?(?=^###\\s+|^##\\s+|\\z)"),
      "",
      out$report,
      perl = TRUE
    )
  }
  out
}
