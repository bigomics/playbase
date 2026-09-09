##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

# =============================================================================
# Multi-omics WGCNA AI report
# =============================================================================
#
# The multi-omics slot (`pgx$wgcna_mox`) has a different shape from the
# single-omics slot (`pgx$wgcna`), so it gets its own report builders here
# rather than being squeezed through ai-report-wgcna.R:
#
#   pgx$wgcna                          pgx$wgcna_mox
#   -----------------------------      ------------------------------------
#   $datExpr / $datTraits / $net       (absent at top level)
#   $stats / $modTraits                (absent at top level)
#   $me.genes  (MEcolor keys)          $me.genes  (all layers, layer-prefixed
#                                       keys: GXblue, PXblue, ...)
#   $gsea      (per module)            $layers$<ln>$gsea (per layer x module);
#                                       the top-level $gsea covers the first
#                                       layer only, so never read it here
#   -                                  $layers$<ln>  one full single-omics
#                                       WGCNA object per omics type
#   -                                  $class = "multiomics", $datanames
#   $annot     (feature-keyed)         $annot (keyed "<ln>:<feature>", while
#                                       layer features are unprefixed)
#
# Consequences, all handled below:
#   - every per-module primitive runs against the OWNING LAYER, never the
#     top-level object (wgcna.ensureStats / wgcna.getGeneStats / eigengenes
#     all need $net + $datExpr + $datTraits);
#   - annotation must be re-keyed per layer (`.moxwgcna_layer_annot`);
#   - module selection is balanced across layers so a single dominant omics
#     type cannot take every slot in the report;
#   - the cross-module section reports CROSS-LAYER eigengene pairs first --
#     that is the signal that only exists in the multi-omics analysis.
#
# The per-module rendering, hub-gene and enrichment leaf renderers are shared
# with ai-report-wgcna.R (single source of truth); only the layer plumbing and
# the cross-omics sections are new.


# -----------------------------------------------------------------------------
# Layer plumbing
# -----------------------------------------------------------------------------

#' The layer list of a multi-omics WGCNA object.
#'
#' Degrades to a single unnamed-omics layer when handed a single-omics object,
#' so every builder below can assume a `list(<ln> = <wgcna>)` shape.
.moxwgcna_layers <- function(mox) {
  if (!is.null(mox$layers) && length(mox$layers) > 0) return(mox$layers)
  list(gx = mox)
}

#' Normalise a feature identifier for annotation lookup.
#'
#' Layer features are stored unprefixed and proteomics groups arrive brace
#' wrapped (`{I1JLC8.16}`), while `mox$annot` is keyed `"<layer>:<feature>"`.
#' Both sides are reduced to the bare identifier before matching.
.moxwgcna_norm_id <- function(x) {
  x <- gsub("^\\{|\\}$", "", as.character(x))
  sub("^[A-Za-z]{2,4}:", "", x)
}

#' Is this the unassigned (grey) module?
#'
#' Not `.wgcna_is_grey()`: that matches any name ending in "grey", which in a
#' layer-prefixed network wrongly swallows the real WGCNA colour modules
#' `GXdarkgrey` / `PXdarkgrey`. Only a bare `grey`/`gray` after the ME- or
#' layer prefix (or WGCNA's reserved integer label 0) is the unassigned bin.
.moxwgcna_is_grey <- function(module) {
  bare <- tolower(sub("^(ME|[A-Z]{2})", "", as.character(module)))
  label <- suppressWarnings(as.integer(sub("^ME", "", module)))
  bare %in% c("grey", "gray") | (!is.na(label) & label == 0L)
}

#' Format a data.frame as markdown table rows, without the header.
#'
#' `wgcna_report_data.md` ships its own static header for the modules summary
#' table, so only the body rows go into the placeholder.
.moxwgcna_mdrows <- function(df) {
  lines <- omicsai::omicsai_format_mdtable(df)
  if (length(lines) > 2) lines <- lines[-(1:2)]
  paste(lines, collapse = "\n")
}

#' Re-key the shared annotation table to one layer's feature vocabulary.
#'
#' Returns a table whose rownames are exactly `features` (in that order), so
#' `probe2symbol()` / `resolve_functions()` / `rename_by2()` resolve against
#' the identifiers the layer actually uses. Features with no annotation row
#' become an all-NA row, which the callers already treat as "unknown".
.moxwgcna_layer_annot <- function(annot, features) {
  if (is.null(annot) || length(features) == 0) return(annot)
  idx <- match(.moxwgcna_norm_id(features), .moxwgcna_norm_id(rownames(annot)))
  if (all(is.na(idx))) return(annot)
  out <- annot[idx, , drop = FALSE]
  rownames(out) <- features
  if ("feature" %in% colnames(out)) out$feature <- features
  out
}

#' Prepare every layer for extraction.
#'
#' Fills stats/labels lazily (older PGX files ship without them), aliases
#' `$gse` onto `$gsea` -- the shared extractors read `$gse` while the
#' multi-omics compute path writes `$gsea` -- and injects a layer-keyed
#' annotation table.
.moxwgcna_prepare_layers <- function(mox, pgx) {
  layers <- .moxwgcna_layers(mox)
  annot <- mox$annot %||% pgx$genes
  out <- list()
  for (ln in names(layers)) {
    layer <- tryCatch(wgcna.ensureStats(layers[[ln]]),
                      error = function(e) layers[[ln]])
    if (is.null(layer$gse)) layer$gse <- layer$gsea
    features <- unique(c(colnames(layer$datExpr), unlist(layer$me.genes)))
    layer$annot <- .moxwgcna_layer_annot(annot, features) %||% layer$annot
    out[[ln]] <- layer
  }
  out
}

#' Map each module to the layer code that owns it.
#' @return Named character vector: module -> layer name (NA when unmatched).
.moxwgcna_module_layers <- function(layers, modules) {
  owner <- rep(NA_character_, length(modules))
  names(owner) <- modules
  for (ln in names(layers)) {
    hit <- modules %in% names(layers[[ln]]$me.genes)
    owner[hit & is.na(owner)] <- ln
  }
  owner
}

#' Read a compute parameter across layers.
#'
#' Collapses to a single value when all layers agree, otherwise renders
#' `"gx=20, px=14"` so the methods block stays honest about per-layer settings.
.moxwgcna_param <- function(layers, field, default = "NA") {
  vals <- vapply(names(layers), function(ln) {
    v <- layers[[ln]][[field]]
    if (is.null(v) || length(v) == 0) return(NA_character_)
    as.character(v[[1]])
  }, character(1))
  vals <- vals[!is.na(vals)]
  if (length(vals) == 0) return(as.character(default))
  u <- unique(unname(vals))
  if (length(u) == 1) return(u)
  paste(sprintf("%s=%s", names(vals), vals), collapse = ", ")
}

#' Sample grouping used for the qualitative eigengene profile.
#'
#' `extract_module_data()` reads `pgx$samples$group`; multi-omics PGX objects
#' frequently have no such column, so fall back to the first categorical
#' sample annotation with more than one and at most eight levels.
.moxwgcna_sample_groups <- function(pgx) {
  samples <- pgx$samples
  if (is.null(samples) || !NCOL(samples)) return(NULL)
  if (!is.null(samples$group)) return(samples$group)
  for (nm in colnames(samples)) {
    x <- samples[[nm]]
    if (is.numeric(x)) next
    lv <- unique(x[!is.na(x)])
    if (length(lv) > 1 && length(lv) <= 8) return(as.character(x))
  }
  NULL
}


# -----------------------------------------------------------------------------
# Module selection and extraction
# -----------------------------------------------------------------------------

#' Select the modules that get a detail block, balanced across omics layers.
#'
#' `wgcna.getTopModules()` already ranks per layer and concatenates, so one
#' dominant layer can fill every slot. Rather than truncating that list, the
#' per-layer rankings are interleaved round-robin: each layer contributes its
#' strongest module first, then its second, and so on, until `n_modules` is
#' reached. Grey is never selected.
.moxwgcna_select_modules <- function(mox, layers, n_modules = 8L) {
  selected <- tryCatch(
    wgcna.getTopModules(mox, min_modules = max(5L, n_modules)),
    error = function(e) character(0)
  )
  if (length(selected) == 0) {
    selected <- unlist(lapply(layers, function(l) names(l$me.genes)),
                       use.names = FALSE)
  }
  selected <- unique(selected[!.moxwgcna_is_grey(selected)])
  if (length(selected) == 0) return(character(0))

  owner <- .moxwgcna_module_layers(layers, selected)
  per_layer <- split(selected, factor(owner, levels = names(layers)))
  per_layer <- Filter(length, per_layer)
  if (length(per_layer) == 0) return(utils::head(selected, n_modules))

  out <- character(0)
  depth <- 1L
  while (length(out) < n_modules && depth <= max(lengths(per_layer))) {
    for (ln in names(per_layer)) {
      if (length(out) >= n_modules) break
      if (depth <= length(per_layer[[ln]])) {
        out <- c(out, per_layer[[ln]][depth])
      }
    }
    depth <- depth + 1L
  }
  out
}

#' Per-module structured data for a multi-omics object.
#'
#' Delegates to the single-omics `.compute_module_data()` per layer (so the
#' hub-gene / enrichment / trait extraction has exactly one implementation)
#' and tags each entry with its layer code and readable omics label.
#'
#' @return Named list keyed by module, in `modules` order.
.moxwgcna_compute_module_data <- function(mox, layers, pgx, modules,
                                          ntop_enrichment = 20L,
                                          ntop_genes = 50L) {
  if (length(modules) == 0) return(list())

  ## Give the shared extractor a usable grouping for the eigengene profile.
  pgx_l <- pgx
  groups <- .moxwgcna_sample_groups(pgx)
  if (!is.null(groups) && !is.null(pgx_l$samples)) pgx_l$samples$group <- groups

  owner <- .moxwgcna_module_layers(layers, modules)
  module_data <- list()
  for (ln in unique(stats::na.omit(owner))) {
    mods <- modules[!is.na(owner) & owner == ln]
    md <- tryCatch(
      .compute_module_data(layers[[ln]], pgx_l, mods,
                           ntop_enrichment = ntop_enrichment,
                           ntop_genes = ntop_genes),
      error = function(e) {
        info(paste0("[moxwgcna] layer '", ln, "' extraction failed: ",
                    conditionMessage(e)))
        list()
      }
    )
    for (m in names(md)) {
      md[[m]]$layer <- ln
      md[[m]]$omics <- .wgcna_layer_omics_label(ln)
      module_data[[m]] <- md[[m]]
    }
  }
  module_data[intersect(modules, names(module_data))]
}


# -----------------------------------------------------------------------------
# Section builders - one per template slot
# -----------------------------------------------------------------------------

#' Per-layer scope table appended to the `## Experiment` block.
#'
#' Tells the model that modules live in separate per-layer networks and that
#' module names are layer-prefixed, which the cross-omics reading depends on.
.moxwgcna_scope_block <- function(mox, layers) {
  rows <- lapply(names(layers), function(ln) {
    layer <- layers[[ln]]
    mods <- names(layer$me.genes)
    grey <- mods[.moxwgcna_is_grey(mods)]
    n_used <- length(unlist(layer$me.genes))
    if (n_used == 0) {
      n_used <- tryCatch(ncol(layer$datExpr), error = function(e) NA_integer_)
    }
    data.frame(
      Layer            = ln,
      Omics            = .wgcna_layer_omics_label(ln),
      `Features used`  = as.character(n_used),
      Modules          = as.character(length(setdiff(mods, grey))),
      `Grey features`  = as.character(length(unlist(layer$me.genes[grey]))),
      Power            = .moxwgcna_param(layers[ln], "power"),
      stringsAsFactors = FALSE, check.names = FALSE
    )
  })
  table_md <- paste(omicsai::omicsai_format_mdtable(do.call(rbind, rows)),
                    collapse = "\n")
  paste(
    "## Multi-omics WGCNA layers",
    paste0("A separate co-expression network was built per omics layer of the ",
           "same samples; modules are therefore layer-local and their names ",
           "are layer-prefixed and globally unique (`GXblue` belongs to the ",
           "`gx` layer, `PXblue` to `px`). Layers are related to each other ",
           "only through eigengene correlation, reported below as cross-omics ",
           "module coordination."),
    table_md,
    sep = "\n\n"
  )
}

#' Build the `## Overview` placeholder values for a multi-omics object.
.moxwgcna_data_overview <- function(mox, layers, pgx) {
  info <- .ai_report_get(pgx, "info", override = mox$experiment)
  n_features_used <- length(unlist(lapply(layers, function(l) unlist(l$me.genes))))
  if (n_features_used == 0) {
    n_features_used <- sum(vapply(layers, function(l) {
      tryCatch(ncol(l$datExpr), error = function(e) 0L)
    }, integer(1)))
  }
  per_layer <- paste(vapply(names(layers), function(ln) {
    n <- length(unlist(layers[[ln]]$me.genes))
    if (n == 0) n <- tryCatch(ncol(layers[[ln]]$datExpr), error = function(e) 0L)
    sprintf("%s %s", n, .wgcna_layer_omics_label(ln))
  }, character(1)), collapse = "; ")

  list(
    experiment       = info$experiment,
    organism         = info$organism,
    n_samples        = info$n_samples,
    n_features_total = info$n_features,
    n_features_used  = sprintf("%d across %d omics layers (%s)",
                               n_features_used, length(layers), per_layer),
    power            = .moxwgcna_param(layers, "power"),
    min_mod_size     = .moxwgcna_param(layers, "minModSize",
                                       mox$settings$minmodsize %||% "10"),
    merge_cut_height = .moxwgcna_param(layers, "mergeCutHeight",
                                       mox$settings$mergeCutHeight %||% "0.15")
  )
}

#' Build the `## Modules summary` table block.
#'
#' Same columns as the single-omics report; the owning omics layer rides along
#' in the Module cell (`GXblue (transcriptomics)`) so the shared template's
#' fixed header stays valid. One grey row per layer is appended.
.moxwgcna_modsummary <- function(mox, layers, module_data, lead_module) {
  module_order <- names(module_data)
  if (length(module_order) == 0) return(list(table = "", footnote = ""))

  fmt_trait_col <- function(trait, r, is_lead_pos = FALSE) {
    if (is.na(trait) || !nzchar(trait)) return("-")
    verbal <- omicsai::omicsai_verbalize_r(r)
    if (is_lead_pos && !is.na(r)) {
      sprintf("%s (%s, r = %+.2f)", trait, verbal, r)
    } else {
      sprintf("%s (%s)", trait, verbal)
    }
  }

  ov_rows <- lapply(module_order, function(m) {
    md <- module_data[[m]]
    is_lead <- identical(m, lead_module)
    data.frame(
      Module = sprintf("%s (%s)", m, md$omics %||% "unknown layer"),
      Genes  = as.character(md$size),
      `Top correlated trait`      = fmt_trait_col(md$top_pos_trait, md$top_pos_r, is_lead),
      `Top anti-correlated trait` = fmt_trait_col(md$top_neg_trait, md$top_neg_r),
      `Enrichment hits`           = as.character(md$n_sig),
      stringsAsFactors = FALSE, check.names = FALSE
    )
  })

  n_total_non_grey <- 0L
  for (ln in names(layers)) {
    mods <- names(layers[[ln]]$me.genes)
    grey <- mods[.moxwgcna_is_grey(mods)]
    n_total_non_grey <- n_total_non_grey + length(setdiff(mods, grey))
    grey_size <- length(unlist(layers[[ln]]$me.genes[grey]))
    if (grey_size > 0) {
      ov_rows[[length(ov_rows) + 1]] <- data.frame(
        Module = sprintf("%s (%s, unassigned)", paste(grey, collapse = "/"),
                         .wgcna_layer_omics_label(ln)),
        Genes = as.character(grey_size),
        `Top correlated trait` = "-",
        `Top anti-correlated trait` = "-",
        `Enrichment hits` = "-",
        stringsAsFactors = FALSE, check.names = FALSE
      )
    }
  }

  list(
    table = .moxwgcna_mdrows(do.call(rbind, ov_rows)),
    footnote = sprintf("Showing %d of %d modules across %d omics layers.",
                       length(module_order), n_total_non_grey, length(layers))
  )
}

#' Build the `## Module-module eigengene correlations` block.
#'
#' Cross-layer pairs are listed first and labelled as such: they are the only
#' relation the multi-omics analysis adds over running WGCNA per layer.
#' Verbal labels only; raw r is stripped at the rendering layer.
.moxwgcna_eigen_cor <- function(mox, modules, min_abs = 0.7) {
  ME <- tryCatch(.wgcna_mox_eigengene_matrix(mox), error = function(e) NULL)
  if (is.null(ME)) return("(no strong correlations detected)")
  me_cols <- intersect(modules, colnames(ME))
  if (length(me_cols) < 2) return("(no strong correlations detected)")

  layer_of <- .moxwgcna_module_layers(.moxwgcna_layers(mox), me_cols)
  me_cor <- stats::cor(ME[, me_cols, drop = FALSE], use = "pairwise")

  cross <- character(0)
  within <- character(0)
  for (i in seq_len(length(me_cols) - 1)) {
    for (j in (i + 1):length(me_cols)) {
      r <- me_cor[i, j]
      if (is.na(r) || abs(r) < min_abs) next
      a <- me_cols[i]
      b <- me_cols[j]
      same <- identical(unname(layer_of[a]), unname(layer_of[b]))
      line <- if (same) {
        sprintf("%s <-> %s: %s", a, b, omicsai::omicsai_verbalize_r(r))
      } else {
        sprintf("%s (%s) <-> %s (%s): %s", a,
                .wgcna_layer_omics_label(unname(layer_of[a])), b,
                .wgcna_layer_omics_label(unname(layer_of[b])),
                omicsai::omicsai_verbalize_r(r))
      }
      if (same) within <- c(within, line) else cross <- c(cross, line)
    }
  }

  if (length(cross) == 0 && length(within) == 0) {
    return("(no strong correlations detected)")
  }
  blocks <- character(0)
  blocks <- c(blocks, paste0(
    "Cross-omics module pairs (eigengenes of modules in DIFFERENT omics layers):\n",
    if (length(cross)) paste(cross, collapse = "\n")
    else "(none reached the threshold)"
  ))
  if (length(within) > 0) {
    blocks <- c(blocks, paste0(
      "Within-layer module pairs:\n", paste(within, collapse = "\n")))
  }
  paste(blocks, collapse = "\n\n")
}

#' Build per-module detail blocks.
#'
#' Each block is the shared single-omics module block, preceded by its layer
#' header and followed by that module's cross-omics coordination table.
.moxwgcna_module_detail <- function(mox, layers, module_data,
                                    families_text = NULL,
                                    ntop_enrichment = 20L,
                                    ntop_genes = 10L) {
  module_order <- names(module_data)
  if (length(module_order) == 0) return("")

  blocks <- vapply(module_order, function(mod) {
    md <- module_data[[mod]]
    ln <- md$layer
    layer <- if (!is.null(ln) && !is.null(layers[[ln]])) layers[[ln]] else mox
    header <- sprintf("Layer: %s (%s)", md$omics %||% "unknown", ln %||% "?")
    body <- .render_module_block(md, mod, layer,
                                 families_text   = families_text,
                                 ntop_enrichment = ntop_enrichment,
                                 ntop_genes      = ntop_genes)
    cross <- tryCatch(.wgcna_cross_omics_section(mox, mod),
                      error = function(e) "")
    parts <- c(header, body, cross)
    paste(parts[nzchar(parts)], collapse = "\n\n")
  }, character(1))

  paste(blocks, collapse = "\n\n")
}

#' Per-module gene-family enrichment across layers.
#'
#' `pgx$families` is keyed by the prefixed feature id (`gx:GLYMA_...`) while
#' layer modules hold bare ids, so ids are re-prefixed before intersecting.
.moxwgcna_families_text <- function(layers, pgx, module_data) {
  if (is.null(pgx$families)) return(NULL)
  fam_names <- names(pgx$families)
  fam_names <- fam_names[fam_names != "<all>"]
  fam_sizes <- lengths(pgx$families[fam_names])
  fam_names <- fam_names[fam_sizes >= 5 & fam_sizes <= 500]
  if (length(fam_names) == 0) return(NULL)

  out <- list()
  for (mod in names(module_data)) {
    ln <- module_data[[mod]]$layer
    if (is.null(ln) || is.null(layers[[ln]])) next
    mod_genes <- layers[[ln]]$me.genes[[mod]]
    if (is.null(mod_genes) || length(mod_genes) == 0) next
    mod_genes <- unique(c(mod_genes, paste0(ln, ":", mod_genes)))

    overlaps <- vapply(fam_names, function(fn) {
      length(intersect(mod_genes, pgx$families[[fn]]))
    }, integer(1))

    sig_fam <- overlaps[overlaps >= 3]
    if (length(sig_fam) == 0) next
    sig_fam <- utils::head(sort(sig_fam, decreasing = TRUE), 5)

    out[[mod]] <- vapply(names(sig_fam), function(fn) {
      sprintf("%s (%d of %d)", fn, sig_fam[fn], length(pgx$families[[fn]]))
    }, character(1))
  }
  out
}


# -----------------------------------------------------------------------------
# Orchestrator
# -----------------------------------------------------------------------------

#' Build structured report tables from multi-omics WGCNA results.
#'
#' Renders the report data block: a single markdown document substituted into
#' `prompts/wgcna/wgcna_report_data.md`, plus a structured `data` list for any
#' downstream callers that want the raw values.
#'
#' @param mox Multi-omics WGCNA slot (`pgx$wgcna_mox`).
#' @param pgx Full pgx object.
#' @param n_modules Total number of modules to detail, balanced across layers.
#' @param ntop_enrichment Enrichment terms per module block.
#' @param ntop_genes Hub genes retained per module before rendering.
#' @param include_module_cors Include the eigengene-correlation section.
#' @param include_families Include the gene-family section.
#' @param include_contrasts Include the contrast matrix in the experiment block.
#' @return list(text = character, data = list)
moxwgcna_build_report_tables <- function(mox, pgx,
                                         n_modules = 8L,
                                         ntop_enrichment = 20L,
                                         ntop_genes = 50L,
                                         include_module_cors = TRUE,
                                         include_families = TRUE,
                                         include_contrasts = TRUE) {
  layers <- .moxwgcna_prepare_layers(mox, pgx)
  modules <- .moxwgcna_select_modules(mox, layers, n_modules = n_modules)

  module_data <- .moxwgcna_compute_module_data(
    mox, layers, pgx, modules,
    ntop_enrichment = ntop_enrichment, ntop_genes = ntop_genes)

  ## Order by significant-enrichment count, then by strongest absolute trait
  ## correlation. The secondary key matters here: enrichment on a non-model
  ## organism is often flat across every module, which would otherwise leave
  ## the ordering (and the lead module) arbitrary.
  if (length(module_data) > 0) {
    n_sig <- vapply(module_data, function(x) as.integer(x$n_sig %||% 0L), integer(1))
    max_r <- vapply(module_data, function(x) {
      r <- c(x$top_r, x$top_pos_r, x$top_neg_r)
      r <- abs(r[!is.na(r)])
      if (length(r) == 0) 0 else max(r)
    }, numeric(1))
    module_data <- module_data[order(-n_sig, -max_r)]
  }
  module_order <- names(module_data)
  lead_module <- if (length(module_order) > 0) module_order[1] else NA_character_

  families_text <- if (include_families) {
    .moxwgcna_families_text(layers, pgx, module_data)
  } else NULL

  module_cors <- if (include_module_cors) {
    .moxwgcna_eigen_cor(mox, module_order)
  } else "(module-module correlations omitted)"

  overview_params <- .moxwgcna_data_overview(mox, layers, pgx)
  modsum <- .moxwgcna_modsummary(mox, layers, module_data, lead_module)
  per_module <- .moxwgcna_module_detail(mox, layers, module_data,
                                        families_text = families_text,
                                        ntop_enrichment = ntop_enrichment,
                                        ntop_genes = ntop_genes)

  experiment_info <- if (include_contrasts) {
    .ai_report_get(pgx, "experiment_info", override = mox$experiment)
  } else {
    .ai_report_get(pgx, "experiment_info", override = mox$experiment,
                   n_contrasts = 0L, contrasts_block = "(contrasts omitted)")
  }
  experiment_info <- paste(experiment_info,
                           .moxwgcna_scope_block(mox, layers), sep = "\n\n")

  tmpl <- omicsai::omicsai_load_template(
    .wgcna_prompt_path("wgcna_report_data.md")
  )

  text <- omicsai::omicsai_substitute_template(tmpl, c(
    overview_params,
    list(
      experiment_info            = experiment_info,
      modules_summary_table      = modsum$table,
      modules_summary_footnote   = modsum$footnote,
      module_module_correlations = module_cors,
      module_detail              = per_module
    )
  ))

  list(
    text = text,
    data = list(
      experiment       = overview_params$experiment,
      organism         = overview_params$organism,
      n_samples        = overview_params$n_samples,
      n_features_total = overview_params$n_features_total,
      n_features_used  = overview_params$n_features_used,
      layers           = names(layers),
      modules          = module_data
    )
  )
}


# -----------------------------------------------------------------------------
# Methods section (deterministic appendix)
# -----------------------------------------------------------------------------

#' Build the deterministic methods section for a multi-omics WGCNA report.
#'
#' Parameters that differ between layers are rendered as `gx=20, px=14`, and a
#' per-layer table is appended as a footer so nothing is silently averaged.
moxwgcna_build_methods <- function(mox, pgx) {
  layers <- .moxwgcna_layers(mox)

  n_wgcna <- length(unlist(lapply(layers, function(l) unlist(l$me.genes))))
  if (n_wgcna == 0) {
    n_wgcna <- sum(vapply(layers, function(l) {
      tryCatch(ncol(l$datExpr), error = function(e) 0L)
    }, integer(1)))
  }

  feature_type <- "features"
  if (!is.null(pgx$datatype)) {
    if (grepl("multi", pgx$datatype, ignore.case = TRUE)) {
      feature_type <- "features (genes, proteins and other omics measurements)"
    } else if (grepl("prot", pgx$datatype, ignore.case = TRUE)) {
      feature_type <- "proteins"
    } else if (grepl("rna|transcript|gene", pgx$datatype, ignore.case = TRUE)) {
      feature_type <- "genes"
    }
  }

  n_samples <- tryCatch(nrow(layers[[1]]$datExpr), error = function(e) {
    tryCatch(nrow(pgx$samples), error = function(e2) NA_integer_)
  })

  all_mods <- unlist(lapply(layers, function(l) names(l$me.genes)), use.names = FALSE)
  grey_mods <- all_mods[.moxwgcna_is_grey(all_mods)]
  n_modules <- length(setdiff(all_mods, grey_mods))
  grey_size <- sum(vapply(layers, function(l) {
    g <- names(l$me.genes)[.moxwgcna_is_grey(names(l$me.genes))]
    length(unlist(l$me.genes[g]))
  }, integer(1)))

  n_genesets <- tryCatch({
    gse <- layers[[1]]$gse %||% layers[[1]]$gsea
    nrow(gse[[setdiff(names(gse), grey_mods)[1]]])
  }, error = function(e) "NA")

  params <- list(
    n_features_wgcna = as.character(n_wgcna),
    feature_type = feature_type,
    n_samples = as.character(n_samples),
    network_type = .moxwgcna_param(layers, "networktype",
                                   mox$settings$networktype %||% "signed"),
    power = .moxwgcna_param(layers, "power"),
    min_mod_size = .moxwgcna_param(layers, "minModSize",
                                   mox$settings$minmodsize %||% "10"),
    merge_cut_height = .moxwgcna_param(layers, "mergeCutHeight",
                                       mox$settings$mergeCutHeight %||% "0.15"),
    min_kme = .moxwgcna_param(layers, "minKME", mox$settings$minKME %||% "0.3"),
    n_modules = as.character(n_modules),
    grey_size = as.character(grey_size),
    n_genesets_tested = as.character(n_genesets)
  )

  layer_rows <- do.call(rbind, lapply(names(layers), function(ln) {
    l <- layers[[ln]]
    mods <- names(l$me.genes)
    data.frame(
      Layer   = ln,
      Omics   = .wgcna_layer_omics_label(ln),
      Features = as.character(length(unlist(l$me.genes))),
      Modules = as.character(length(setdiff(mods, mods[.moxwgcna_is_grey(mods)]))),
      Power   = .moxwgcna_param(layers[ln], "power"),
      `Min module size` = .moxwgcna_param(layers[ln], "minModSize"),
      stringsAsFactors = FALSE, check.names = FALSE
    )
  }))
  footer <- paste(
    "### **Multi-omics layers**",
    paste0("WGCNA was computed independently per omics layer; modules were ",
           "related across layers by correlating their eigengenes."),
    paste(omicsai::omicsai_format_mdtable(layer_rows), collapse = "\n"),
    sep = "\n\n"
  )

  build_report_methods("wgcna", "wgcna_methods.md", params = params,
                       footer = footer)
}


# -----------------------------------------------------------------------------
# Prompt assembly
# -----------------------------------------------------------------------------

#' Assemble the multi-omics WGCNA static-report prompt from omicsai fragments.
#'
#' @param slice the multi-omics WGCNA result slot (`pgx$wgcna_mox`).
#' @param pgx full pgx object (used for organism + annotation).
#' @param ai resolved `ai` list (already passed through `.ai_resolve_defaults`).
#' @return `list(system = <character>, board = <character>)`
#' @keywords internal
moxwgcna_assemble_prompt <- function(slice, pgx, ai) {
  if (!requireNamespace("omicsai", quietly = TRUE)) {
    stop("omicsai package required for AI report generation", call. = FALSE)
  }
  data_block <- moxwgcna_build_report_tables(slice, pgx)$text
  .ai_report_build_prompt(pgx, "wgcna_mox", data_block)
}


# -----------------------------------------------------------------------------
# Entry point
# -----------------------------------------------------------------------------

#' Generate a multi-omics WGCNA AI report.
#'
#' Objects without a `$layers` list are not multi-omics (an older or
#' mis-populated slot); those are handed to the single-omics report so the
#' caller still gets a report rather than an error.
#'
#' @param pgx full pgx object.
#' @param slice WGCNA multi-omics result slot (`pgx$wgcna_mox`).
#' @param ai resolved `ai` list.
#' @return `list(report = <markdown>, prompt = <markdown>)`.
ai.wgcna_mox.create_report <- function(pgx, slice, ai) {
  if (!requireNamespace("omicsai", quietly = TRUE)) {
    stop("omicsai package required for AI report generation", call. = FALSE)
  }
  if (is.null(slice$layers) || length(slice$layers) == 0) {
    return(ai.wgcna.create_report(pgx, slice, ai))
  }
  bp  <- moxwgcna_assemble_prompt(slice, pgx, ai)
  out <- .ai_report_run_prompt(bp, ai)

  ## Deterministic methods appendix, appended after generation so it cannot
  ## influence the prompt.
  out$report <- paste(out$report, moxwgcna_build_methods(slice, pgx),
                      sep = "\n\n")
  out
}
