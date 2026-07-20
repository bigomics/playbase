##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

# =============================================================================
# Durable per-module WGCNA AI summaries
# =============================================================================
#
# These are distinct from the multi-module WGCNA *report* (ai-report-wgcna.R):
# one short summary per co-expression module, precomputed at compute time and
# stored so the board can render them instantly instead of calling the LLM live.
#
# Storage (sibling of the report slots, keyed by module name):
#   pgx$ai$wgcna$extras$<module>
#   pgx$ai$wgcna_mox$extras$<module>       (module names are layer-prefixed,
#                                           e.g. GXblue/MIblue/PXblue, so a flat
#                                           map is collision-free across layers)
# Each entry: list(summary, prompt, usage, model, created_at, edited, edited_at)
#
# The data block is built from the same per-module primitives the report uses
# (extract_module_data / .compute_module_data / .render_module_block), so there
# is a single source of truth and the board keeps no parallel extraction code.


# -----------------------------------------------------------------------------
# Prompt assembly (per module)
# -----------------------------------------------------------------------------

# Is this the WGCNA grey module (unassigned features)? Vectorised. Handles both
# label conventions: colour names ("grey"/"MEgrey"/"gray") and integer labels,
# where WGCNA reserves label 0 (ME0) for the grey/unassigned module -- the same
# convention playbase uses elsewhere via WGCNA::labels2colors().
.wgcna_is_grey <- function(module) {
  is_colour_grey <- grepl("gr[ae]y$", tolower(module))
  suffix <- suppressWarnings(as.integer(sub("^ME", "", module)))
  is_label0_grey <- !is.na(suffix) & suffix == 0L
  is_colour_grey | is_label0_grey
}

# Ad-hoc, omics-adapted note for the grey module. The grey module is not a true
# co-expression module, so we do not run the biological-interpretation data
# block on it; we inject a short neutral note instead. Substitution goes through
# the omicsai backend (never hand-rolled string assembly), and modality nouns
# come from omicsai_datatype_vocab() so the wording adapts to transcriptomics /
# proteomics / metabolomics / methylomics.
.wgcna_grey_data_block <- function(wgcna, module, vocab) {
  n_genes <- length(wgcna$me.genes[[module]])
  # Name the module by its colour (MEgrey), matching how non-grey blocks label
  # ME_color, even when the stored key is an integer label such as ME0.
  me_color <- display_module(module, wgcna)
  omicsai::omicsai_substitute_template(
    paste(
      "Module {{ME_color}} is the WGCNA grey module. In weighted",
      "gene co-expression network analysis the grey module is not a genuine",
      "co-expression module: it is the container for the {{n_genes}}",
      "{{features}} that could not be confidently assigned to any module and",
      "therefore has no coherent eigengene, trait association, or enrichment",
      "signature. Summarise it in two to three sentences as unassigned",
      "background {{features}}, and explicitly avoid inferring any biological",
      "function, pathway, or trait relevance for this module.",
      sep = " "
    ),
    list(
      ME_color = me_color,
      n_genes = as.character(n_genes),
      features = vocab$features %||% "features"
    )
  )
}

# Build the system+board prompt for one module's summary. Mirrors the live board
# card (summary_prompt with an interpretation context) but takes an
# already-rendered/lazy data block so grey and non-grey share one path. `context`
# selects the WGCNA-flavour interpretation fragment (standard / consensus /
# preservation / multiomics), so every flavour routes through this one assembler.
.wgcna_summary_build_prompt <- function(data_block, vocab, species,
                                        context = "wgcna/wgcna_interpretation") {
  prompt <- omicsai::summary_prompt(
    role    = omicsai::frag("system_base", vocab),
    task    = NULL,
    species = species,
    context = omicsai::frag(context, vocab),
    data    = data_block
  )
  omicsai::build_prompt(prompt)
}

# Shape one LLM result into a durable extras entry.
.wgcna_summary_entry <- function(res, ai) {
  list(
    summary    = res$report,
    prompt     = res$prompt,
    usage      = res$usage,
    model      = ai$llm_model,
    created_at = as.numeric(Sys.time()),
    edited     = FALSE,
    edited_at  = ""
  )
}

# Build the module data block for one module: the grey note for grey modules,
# otherwise the rendered per-module template. `md` is the precomputed
# .compute_module_data() element (recomputed for a single module when absent).
.wgcna_module_data_block <- function(wgcna, pgx, module, vocab, md = NULL) {
  if (.wgcna_is_grey(module)) {
    return(.wgcna_grey_data_block(wgcna, module, vocab))
  }
  if (is.null(md)) md <- .compute_module_data(wgcna, pgx, module)[[module]]
  .render_module_block(md, module, wgcna)
}

# Generate one module's summary entry. `md` is the precomputed
# .compute_module_data() element for this module (NULL for grey, which does not
# use a data block). `wgcna` must already have annotation + stats filled.
.wgcna_one_module_summary <- function(pgx, wgcna, module, ai, vocab, species, md = NULL) {
  data_block <- .wgcna_module_data_block(wgcna, pgx, module, vocab, md = md)
  bp <- .wgcna_summary_build_prompt(data_block, vocab, species)
  res <- .ai_report_run_prompt(bp, ai)
  .wgcna_summary_entry(res, ai)
}


# -----------------------------------------------------------------------------
# Slot-level generation (one single-omics WGCNA object -> named extras list)
# -----------------------------------------------------------------------------

# Summarise every module of a single-omics WGCNA object. `existing` carries
# forward already-stored entries; unless ai$force is TRUE, present modules are
# kept as-is (this also preserves user-edited summaries). Grey is included.
.wgcna_summarize_modules <- function(pgx, wgcna, ai, existing = NULL, annot = NULL) {
  wgcna <- wgcna.ensureStats(wgcna)
  # Multi-omics layers carry no per-layer annotation; inject the parent object's
  # annotation (or pgx$genes) so symbol/function resolution works per layer.
  if (is.null(wgcna$annot)) wgcna$annot <- annot %||% pgx$genes

  modules <- names(wgcna$me.genes)
  if (length(modules) == 0) return(existing %||% list())

  out <- existing %||% list()
  vocab   <- omicsai::omicsai_datatype_vocab(pgx$datatype)
  species <- omicsai::omicsai_species_prompt(pgx$organism)

  todo <- if (isTRUE(ai$force)) modules else modules[!modules %in% names(out)]
  if (length(todo) == 0) return(out)

  # Compute per-module data once for all non-grey todo modules (shares top/M).
  nongrey <- todo[!.wgcna_is_grey(todo)]
  md_all <- if (length(nongrey) > 0) {
    .compute_module_data(wgcna, pgx, nongrey)
  } else {
    list()
  }

  for (mod in todo) {
    out[[mod]] <- tryCatch(
      .wgcna_one_module_summary(pgx, wgcna, mod, ai, vocab, species, md = md_all[[mod]]),
      error = function(e) {
        msg <- paste0("[pgx.update_wgcna_summaries] module '", mod,
                      "' failed: ", conditionMessage(e))
        if (identical(ai$on_error, "abort")) stop(msg, call. = FALSE)
        if (identical(ai$on_error, "warn")) warning(msg, call. = FALSE)
        else info(msg)
        out[[mod]]  # keep any prior entry on failure
      }
    )
  }
  out
}

# -----------------------------------------------------------------------------
# Public entry points
# -----------------------------------------------------------------------------

#' Build the system + board prompt for one WGCNA module summary.
#'
#' Exposes the exact prompt the durable generator uses, so the board's live and
#' "Regenerate" paths assemble their prompt from playbase (single source of
#' truth) instead of re-deriving module data. The returned parts feed straight
#' into `omicsai::omicsai_config(system_prompt = )` and the AI text card.
#'
#' `board_type` dispatches to a WGCNA-flavour-specific data block and
#' interpretation context: `"standard"` (single-omics), `"multiomics"` (one
#' layer's module plus its cross-omics coordination), `"consensus"` (trait links
#' preserved across groups), and `"preservation"` (whether a reference module
#' reproduces in the test group). Each non-standard flavour has a different
#' object shape and dedicated extractors.
#'
#' @param pgx PGX object (organism, datatype, annotation fallback).
#' @param wgcna A WGCNA object: a standard slot for `"standard"`, a full
#'   multi-omics object (`pgx$wgcna_mox`) for `"multiomics"`, a consensus object
#'   for `"consensus"`, a preservation object for `"preservation"`.
#' @param module Module name such as `"MEblue"` / `"GXblue"`.
#' @param annot Optional annotation table; used when `wgcna$annot` is absent
#'   (standard only).
#' @param board_type One of `"standard"`, `"multiomics"`, `"consensus"`,
#'   `"preservation"`.
#' @return `list(system = <character>, board = <character>)`.
#' @export
wgcna.module_summary_prompt <- function(pgx, wgcna, module, annot = NULL,
                                        board_type = c("standard", "multiomics",
                                                       "consensus",
                                                       "preservation")) {
  board_type <- match.arg(board_type)
  if (identical(board_type, "multiomics")) {
    return(.wgcna_multiomics_prompt(pgx, wgcna, module))
  }
  if (identical(board_type, "consensus")) {
    return(.wgcna_consensus_prompt(pgx, wgcna, module))
  }
  if (identical(board_type, "preservation")) {
    return(.wgcna_preservation_prompt(pgx, wgcna, module))
  }
  wgcna <- wgcna.ensureStats(wgcna)
  if (is.null(wgcna$annot)) wgcna$annot <- annot %||% pgx$genes
  vocab   <- omicsai::omicsai_datatype_vocab(pgx$datatype)
  species <- omicsai::omicsai_species_prompt(pgx$organism)
  data_block <- .wgcna_module_data_block(wgcna, pgx, module, vocab)
  .wgcna_summary_build_prompt(data_block, vocab, species)
}

#' Find the multi-omics WGCNA layer that owns a (layer-prefixed) module.
#'
#' Module names in a `wgcna_mox` object are layer-prefixed and globally unique
#' (e.g. `GXblue` in the `gx` layer). Returns the layer sub-object that contains
#' the module, or `NULL` if none does.
#'
#' @param mox A multi-omics WGCNA object with a `$layers` list.
#' @param module Layer-prefixed module name.
#' @return The single-omics layer object, or `NULL`.
#' @export
wgcna.mox_layer_for_module <- function(mox, module) {
  if (is.null(mox$layers)) return(NULL)
  for (ln in names(mox$layers)) {
    if (module %in% names(mox$layers[[ln]]$me.genes)) return(mox$layers[[ln]])
  }
  NULL
}

#' Precompute durable per-module AI summaries for a PGX object.
#'
#' Iterates the standard (`pgx$wgcna`) and multi-omics (`pgx$wgcna_mox`, per
#' layer x module) WGCNA slots and stores one summary per module under
#' `pgx$ai$<slot>$extras$<module>`, alongside the existing report slots.
#' Consensus and Preservation are intentionally excluded: their objects are
#' computed live in the board from user input and are handled ephemerally.
#'
#' Existing entries are preserved unless `ai$force` is TRUE (which regenerates
#' all), so a compute-time refresh never clobbers user-edited summaries.
#'
#' @param pgx PGX object.
#' @param ai Report-generation options; `NULL` returns `pgx` unchanged.
#' @return PGX object with `pgx$ai$wgcna$extras` / `pgx$ai$wgcna_mox$extras`
#'   populated where those slots exist.
#' @export
pgx.update_wgcna_summaries <- function(pgx, ai = NULL) {
  if (is.null(ai)) return(pgx)
  ai <- .ai_resolve_defaults(ai)
  if (is.null(pgx$ai)) pgx$ai <- list()

  if (!is.null(pgx$wgcna)) {
    info("[pgx.update_wgcna_summaries] summarising standard WGCNA modules...")
    pgx$ai$wgcna$extras <- .wgcna_summarize_modules(
      pgx, pgx$wgcna, ai,
      existing = pgx$ai$wgcna$extras
    )
  }

  if (!is.null(pgx$wgcna_mox)) {
    info("[pgx.update_wgcna_summaries] summarising multi-omics WGCNA modules...")
    pgx$ai$wgcna_mox$extras <- .wgcna_mox_summarize(
      pgx, pgx$wgcna_mox, ai,
      existing = pgx$ai$wgcna_mox$extras
    )
  }

  pgx
}
