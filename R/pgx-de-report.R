##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

# =============================================================================
# Differential-expression AI-report extraction
# =============================================================================
# Owns deterministic PGX -> prompt-data extraction for the integrated DE report.
# LLM orchestration lives in ai-report.R alongside the other report dispatch.


# -----------------------------------------------------------------------------
# Metadata extraction
# -----------------------------------------------------------------------------

.de_annotate_genes <- function(pgx, features) {
  annot <- pgx$genes
  out <- data.frame(feature = features, symbol = features, title = "",
                    stringsAsFactors = FALSE)
  if (is.null(annot) || !nrow(annot)) return(out)

  idx <- match(features, rownames(annot))
  if ("symbol" %in% colnames(annot)) {
    symbols <- as.character(annot[idx, "symbol"])
    keep <- !is.na(symbols) & nzchar(symbols)
    out$symbol[keep] <- symbols[keep]
  }
  title_col <- intersect(c("gene_title", "gene_name", "description"),
                         colnames(annot))[1]
  if (!is.na(title_col)) {
    titles <- as.character(annot[idx, title_col])
    titles[is.na(titles)] <- ""
    out$title <- titles
  }
  out
}


# -----------------------------------------------------------------------------
# DE evidence extraction
# -----------------------------------------------------------------------------

.de_gene_table <- function(pgx, mx, direction, n = 12L) {
  keep <- is.finite(mx$meta.fx) & is.finite(mx$meta.q)
  keep <- keep & if (direction == "up") mx$meta.fx > 0 else mx$meta.fx < 0
  significant <- keep & mx$meta.q < 0.05
  candidates <- which(if (any(significant)) significant else keep)
  if (!length(candidates)) return(data.frame())

  ord <- candidates[order(mx$meta.q[candidates], -abs(mx$meta.fx[candidates]),
                          na.last = TRUE)]
  ord <- head(ord, n)
  info <- .de_annotate_genes(pgx, rownames(mx)[ord])

  qmat <- as.matrix(mx$q)
  fcmat <- as.matrix(mx$fc)
  support <- rowSums(qmat[ord, , drop = FALSE] < 0.05, na.rm = TRUE)
  sign_consistent <- rowMeans(
    sign(fcmat[ord, , drop = FALSE]) == sign(mx$meta.fx[ord]),
    na.rm = TRUE
  )

  data.frame(
    gene = info$symbol,
    annotation = info$title,
    logFC = omicsai::omicsai_format_num(mx$meta.fx[ord], 3L),
    q = omicsai::omicsai_format_pvalue(mx$meta.q[ord]),
    avg_first = omicsai::omicsai_format_num(mx$avg.0[ord], 3L),
    avg_second = omicsai::omicsai_format_num(mx$avg.1[ord], 3L),
    significant_methods = paste0(support, "/", ncol(qmat)),
    sign_agreement = paste0(round(100 * sign_consistent), "%"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

.de_contrast_summary <- function(pgx, slice) {
  rows <- lapply(names(slice$meta), function(contrast) {
    mx <- slice$meta[[contrast]]
    good <- is.finite(mx$meta.fx) & is.finite(mx$meta.q)
    significant <- good & mx$meta.q < 0.05
    data.frame(
      contrast = contrast,
      modeled_groups = .ai_report_get(pgx, "contrast_groups",
                                      contrast = contrast),
      tested = sum(good),
      significant = sum(significant),
      up = sum(significant & mx$meta.fx > 0),
      down = sum(significant & mx$meta.fx < 0),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

.de_cross_contrast_table <- function(pgx, slice, n = 35L) {
  contrasts <- names(slice$meta)
  features <- rownames(slice$meta[[1]])
  fx <- sapply(slice$meta, function(mx) mx$meta.fx)
  q <- sapply(slice$meta, function(mx) mx$meta.q)

  if (is.null(dim(fx))) {
    fx <- matrix(fx, ncol = 1L, dimnames = list(features, contrasts))
    q <- matrix(q, ncol = 1L, dimnames = list(features, contrasts))
  } else {
    rownames(fx) <- features
    rownames(q) <- features
  }

  sig_count <- rowSums(q < 0.05, na.rm = TRUE)
  score <- sig_count * 10 + apply(abs(fx), 1L, max, na.rm = TRUE)
  keep <- head(order(score, decreasing = TRUE), n)
  info <- .de_annotate_genes(pgx, rownames(fx)[keep])
  formatted <- apply(fx[keep, , drop = FALSE], c(1L, 2L),
                     omicsai::omicsai_format_num, digits = 3L)

  data.frame(
    gene = info$symbol,
    significant_contrasts = sig_count[keep],
    formatted,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

.de_contrast_detail <- function(pgx, slice, ntop) {
  details <- vapply(names(slice$meta), function(contrast) {
    mx <- slice$meta[[contrast]]
    up <- .de_gene_table(pgx, mx, "up", n = ntop)
    down <- .de_gene_table(pgx, mx, "down", n = ntop)
    up_table <- if (nrow(up)) {
      paste(omicsai::omicsai_format_mdtable(up), collapse = "\n")
    } else {
      "(none)"
    }
    down_table <- if (nrow(down)) {
      paste(omicsai::omicsai_format_mdtable(down), collapse = "\n")
    } else {
      "(none)"
    }
    paste0(
      "## ", contrast, "\n\n",
      "Modeled groups: ", .ai_report_get(pgx, "contrast_groups",
                                         contrast = contrast), "\n\n",
      "### Reliable upregulated markers\n\n",
      up_table, "\n\n",
      "### Reliable downregulated markers\n\n",
      down_table
    )
  }, character(1))
  paste(details, collapse = "\n\n---\n\n")
}


# -----------------------------------------------------------------------------
# Report data and prompt assembly
# -----------------------------------------------------------------------------

#' Build deterministic evidence tables for an integrated DE AI report.
#'
#' @param slice PGX gx.meta slot.
#' @param pgx Full PGX object.
#' @param ntop Maximum marker genes per direction and contrast.
#' @param cross_n Maximum rows in the cross-contrast footprint.
#' @return List with rendered `text` and structured `data`.
de_build_report_tables <- function(slice, pgx, ntop = 12L, cross_n = 35L) {
  if (!requireNamespace("omicsai", quietly = TRUE)) {
    stop("omicsai package required for AI report generation", call. = FALSE)
  }
  contrasts <- names(slice$meta %||% list())
  if (!length(contrasts)) {
    stop("DE report requires at least one contrast", call. = FALSE)
  }

  sample_metadata <- .ai_report_get(pgx, "sample_metadata")
  contrast_summary <- .de_contrast_summary(pgx, slice)
  cross_contrast <- .de_cross_contrast_table(pgx, slice, n = cross_n)
  contrast_detail <- .de_contrast_detail(pgx, slice, ntop = ntop)
  info <- .ai_report_get(pgx, "info",
                         n_contrasts = length(contrasts),
                         fallback = "(unnamed)")

  params <- c(info, list(
    sample_metadata_table = paste(
      omicsai::omicsai_format_mdtable(sample_metadata), collapse = "\n"
    ),
    contrast_summary_table = paste(
      omicsai::omicsai_format_mdtable(contrast_summary), collapse = "\n"
    ),
    cross_contrast_table = paste(
      omicsai::omicsai_format_mdtable(cross_contrast), collapse = "\n"
    ),
    contrast_detail = contrast_detail
  ))

  template <- omicsai::omicsai_load_template(
    .ai_report_prompt_path("de", "de_report_data.md")
  )
  list(
    text = omicsai::omicsai_substitute_template(template, params),
    data = list(
      sample_metadata = sample_metadata,
      contrast_summary = contrast_summary,
      cross_contrast = cross_contrast
    )
  )
}

#' Assemble the integrated DE report prompt.
#'
#' @param slice PGX gx.meta slot.
#' @param pgx Full PGX object.
#' @param ai Resolved AI-report options.
#' @return List with `system` and `board` prompt strings.
#' @keywords internal
de_assemble_prompt <- function(slice, pgx, ai) {
  ntop <- min(as.integer(ai$ntop), 12L)
  data_block <- de_build_report_tables(slice, pgx, ntop = ntop)$text
  .ai_report_build_prompt(pgx, "de", data_block)
}
