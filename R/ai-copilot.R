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

#' Build compact text context for Copilot.
#'
#' @export
ai.build_report_context <- function(pgx, sections = NULL, collate = TRUE,
                                    ntop = 20, max_chars = Inf, ...) {
  all.sections <- c(
    "description", "dataset_info", "compute_settings",
    "differential_expression", "geneset_enrichment",
    "drug_similarity", "pcsf_report", "wgcna_report", "mofa_report"
  )
  if (is.null(sections)) sections <- all.sections
  sections <- intersect(sections, all.sections)

  ai <- pgx$ai
  if (is.null(ai) || !is.list(ai)) ai <- list()

  out <- list()

  if ("description" %in% sections) {
    out$description <- .ai_copilot_section(
      "Description",
      c(
        name = if (is.null(pgx$name)) "" else pgx$name,
        title = if (is.null(pgx$description)) "" else pgx$description,
        datatype = if (is.null(pgx$datatype)) "" else pgx$datatype,
        organism = if (is.null(pgx$organism)) "" else pgx$organism
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

  if ("differential_expression" %in% sections &&
      is.list(ai$de) && is.character(ai$de$report) &&
      length(ai$de$report) && nzchar(ai$de$report[[1]])) {
    out$differential_expression <- paste0(
      "## Differential Expression\n\n", ai$de$report[[1]], "\n"
    )
  }

  if ("geneset_enrichment" %in% sections &&
      is.list(ai$pathways) && is.character(ai$pathways$report) &&
      length(ai$pathways$report) && nzchar(ai$pathways$report[[1]])) {
    out$geneset_enrichment <- paste0(
      "## Geneset Enrichment\n\n", ai$pathways$report[[1]], "\n"
    )
  }

  if ("wgcna_report" %in% sections &&
      is.list(ai$wgcna) && is.character(ai$wgcna$report) &&
      length(ai$wgcna$report) && nzchar(ai$wgcna$report[[1]])) {
    out$wgcna_report <- paste0("## WGCNA\n\n", ai$wgcna$report[[1]], "\n")
  }

  if ("mofa_report" %in% sections &&
      is.list(ai$mofa) && is.character(ai$mofa$report) &&
      length(ai$mofa$report) && nzchar(ai$mofa$report[[1]])) {
    out$mofa_report <- paste0("## MOFA\n\n", ai$mofa$report[[1]], "\n")
  }

  if ("drug_similarity" %in% sections) {
    drug.slots <- grep("^drugs_", names(ai), value = TRUE)
    drug.reports <- list()
    for (slot in drug.slots) {
      rpt <- ai[[slot]]
      if (!is.list(rpt) || !is.character(rpt$report) ||
          !length(rpt$report) || !nzchar(rpt$report[[1]])) next

      label <- gsub("_", " ", sub("^drugs_", "", slot))
      dbs <- tryCatch(names(pgx$drugs), error = function(e) character(0))
      if (length(dbs)) {
        safe <- paste0("drugs_", gsub("[^A-Za-z0-9]+", "_", dbs))
        idx <- match(slot, safe)
        if (!is.na(idx)) label <- dbs[[idx]]
      }
      if (grepl("^L1000([_/ -]|$)", label, ignore.case = TRUE)) {
        if (grepl("activ", label, ignore.case = TRUE)) label <- "L1000 Activity"
        if (grepl("gene", label, ignore.case = TRUE)) label <- "L1000 Gene"
      }

      drug.reports[[slot]] <- paste0("## ", label, "\n\n", rpt$report[[1]], "\n")
    }
    if (length(drug.reports)) {
      out$drug_similarity <- paste(unlist(drug.reports), collapse = "\n")
    }
  }

  out <- out[!vapply(out, is.null, logical(1))]

  if (isTRUE(collate)) {
    out <- paste(unlist(out), collapse = "\n")
    if (is.finite(max_chars) && nchar(out) > max_chars) {
      out <- paste0(substr(out, 1, max_chars), "\n\n[Report context truncated]")
    }
    return(out)
  }

  if (is.finite(max_chars)) {
    out <- lapply(out, function(x) {
      if (nchar(x) > max_chars) {
        paste0(substr(x, 1, max_chars), "\n\n[Report context truncated]")
      } else {
        x
      }
    })
  }

  out
}
