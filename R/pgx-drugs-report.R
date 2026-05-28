##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

## ============================================================================
## Drug Connectivity Report Utilities
## ----------------------------------------------------------------------------
## Lightweight extractors and tier classifiers consumed by the drug
## connectivity AI report (omicsplayground board.drugconnectivity) and
## by any future board or pipeline that needs to summarise drug
## enrichment results.
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
#' @export
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
#' @export
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
#'   board only has the drug-level table and needs MOA-level signal on the
#'   fly.
#'
#' @export
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
#' @export
drugs.supportBucket <- function(padj, pval = NA_real_) {
  if (!is.na(padj) && padj < 0.05) return("significant")
  if (!is.na(pval) && pval < 0.05) return("nominal")
  "unsupported"
}
