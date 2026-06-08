##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

.ai_copilot_section <- function(title, value, ntop = 20) {
  if (is.null(value)) return(NULL)
  if (is.data.frame(value) || is.matrix(value)) {
    value <- utils::capture.output(print(utils::head(value, ntop)))
  } else if (is.list(value)) {
    value <- utils::capture.output(str(value, max.level = 2, list.len = ntop))
  }
  value <- paste(as.character(value), collapse = "\n")
  if (!nzchar(value)) return(NULL)
  paste0("## ", title, "\n\n", value, "\n")
}

#' Build compact text context for the legacy copilot.
#'
#' @export
ai.build_report_context <- function(pgx, sections = NULL, collate = TRUE,
                                    ntop = 20, ...) {
  all_sections <- c(
    "description",
    "dataset_info",
    "compute_settings",
    "differential_expression",
    "geneset_enrichment",
    "drug_similarity",
    "pcsf_report",
    "wgcna_report",
    "mofa_report"
  )
  if (is.null(sections)) sections <- all_sections
  sections <- intersect(sections, all_sections)

  out <- list()
  if ("description" %in% sections) {
    out$description <- .ai_copilot_section(
      "Description",
      c(
        name = pgx$name %||% "",
        title = pgx$description %||% "",
        datatype = pgx$datatype %||% "",
        organism = pgx$organism %||% ""
      ),
      ntop = ntop
    )
  }
  if ("dataset_info" %in% sections) {
    out$dataset_info <- .ai_copilot_section(
      "Dataset Info",
      list(
        samples = if (!is.null(pgx$samples)) utils::head(pgx$samples, ntop) else NULL,
        contrasts = if (!is.null(pgx$contrasts)) utils::head(pgx$contrasts, ntop) else NULL
      ),
      ntop = ntop
    )
  }
  if ("compute_settings" %in% sections) {
    out$compute_settings <- .ai_copilot_section("Compute Settings", pgx$settings, ntop)
  }
  if ("differential_expression" %in% sections) {
    out$differential_expression <- .ai_copilot_section(
      "Differential Expression",
      pgx$gx.meta,
      ntop
    )
  }
  if ("geneset_enrichment" %in% sections) {
    out$geneset_enrichment <- .ai_copilot_section("Geneset Enrichment", pgx$gset.meta, ntop)
  }
  if ("drug_similarity" %in% sections) {
    out$drug_similarity <- .ai_copilot_section("Drug Similarity", pgx$drugs, ntop)
  }
  if ("pcsf_report" %in% sections) {
    out$pcsf_report <- .ai_copilot_section("PCSF", pgx$pcsf, ntop)
  }
  if ("wgcna_report" %in% sections) {
    out$wgcna_report <- .ai_copilot_section("WGCNA", pgx$wgcna, ntop)
  }
  if ("mofa_report" %in% sections) {
    out$mofa_report <- .ai_copilot_section("MOFA", pgx$mofa, ntop)
  }

  out <- out[!vapply(out, is.null, logical(1))]
  if (isTRUE(collate)) return(paste(unlist(out), collapse = "\n"))
  out
}
