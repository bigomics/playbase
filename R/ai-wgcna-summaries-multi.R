##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

# =============================================================================
# Consensus / Preservation WGCNA per-module summary extraction
# =============================================================================
#
# Consensus and Preservation WGCNA objects have a different shape from the
# standard single-omics object: modules live in $net$colors / $colors, traits in
# $modTraits + per-group $stats, and preservation in $Zsummary / $medianRank.
# The standard extract_module_data() reads $me.genes (absent here) and would
# report every module as empty. These builders extract the fields that make each
# board distinctive: consensus = which trait links hold across groups;
# preservation = whether a reference module's network reproduces in the test set.


# -- shared helpers -----------------------------------------------------------

# Render a named vector of trait correlations as a two-column markdown table
# (trait, verbal strength), keeping only |r| >= min_abs, top ntop.
.wgcna_cor_table <- function(cors, min_abs = 0.3, ntop = 6L, empty = "None.") {
  cors <- cors[!is.na(cors)]
  cors <- cors[abs(cors) >= min_abs]
  if (length(cors) == 0) return(empty)
  cors <- cors[order(-abs(cors))]
  cors <- utils::head(cors, ntop)
  df <- data.frame(
    trait = names(cors),
    strength = omicsai::omicsai_verbalize_r(unname(cors)),
    stringsAsFactors = FALSE
  )
  paste(omicsai::omicsai_format_mdtable(df), collapse = "\n")
}

# Resolve feature IDs to a symbol + short-function two-column hub table.
.wgcna_hub_symbol_table <- function(features, pgx, empty = "None.") {
  if (length(features) == 0) return(empty)
  symbols <- probe2symbol(features, pgx$genes, "symbol", fill_na = TRUE)
  funcs <- resolve_functions(features, pgx$genes)
  paste(omicsai::omicsai_format_mdtable(data.frame(
    gene = symbols,
    known_function = funcs,
    stringsAsFactors = FALSE
  )), collapse = "\n")
}

.wgcna_zsummary_verbal <- function(z) {
  if (is.na(z)) return("preservation not computed")
  if (z >= 10) return("strong preservation")
  if (z >= 2) return("weak to moderate preservation")
  "not preserved"
}


# -- consensus ----------------------------------------------------------------

# Traits whose per-group correlations disagree: opposite sign with either side
# meaningful, or strong in one group and near-zero in the other. Returned as a
# markdown table (trait, one verbal column per group) or a placeholder string.
.wgcna_consensus_divergent_table <- function(gcor, grp) {
  if (length(gcor) < 2) return("None.")
  a <- gcor[[1]]
  b <- gcor[[2]]
  traits <- intersect(names(a), names(b))
  rows <- list()
  for (tr in traits) {
    ha <- a[[tr]]
    hb <- b[[tr]]
    if (is.na(ha) || is.na(hb)) next
    hi <- max(abs(ha), abs(hb))
    lo <- min(abs(ha), abs(hb))
    divergent <- (sign(ha) != sign(hb) && hi >= 0.3) || (hi >= 0.3 && lo < 0.15)
    if (divergent) {
      rows[[tr]] <- data.frame(
        trait = tr,
        g1 = omicsai::omicsai_verbalize_r(ha),
        g2 = omicsai::omicsai_verbalize_r(hb),
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) == 0) return("None.")
  df <- do.call(rbind, rows)
  colnames(df) <- c("trait", grp[1], grp[2])
  paste(omicsai::omicsai_format_mdtable(df), collapse = "\n")
}

# Hub table for a consensus module, using cross-group gene stats. Genes with a
# consistent (consensus == "Y") membership are preferred and ranked by the
# weaker of the two per-group memberships (so both groups must agree).
.wgcna_consensus_hub <- function(cons, module, trait, pgx, ntop = 12L) {
  if (is.na(trait)) return(list(table = "None.", n_consensus = 0L))
  gs <- tryCatch(
    wgcna.getConsensusGeneStats(cons, stats = cons$stats, trait = trait, module = module),
    error = function(e) NULL
  )
  full <- gs$full
  if (is.null(full) || nrow(full) == 0) return(list(table = "None.", n_consensus = 0L))
  n_consensus <- sum(full$consensus == "Y", na.rm = TRUE)
  mmcols <- grep("^moduleMembership", colnames(full), value = TRUE)
  if (length(mmcols) >= 2) {
    min_mm <- apply(abs(full[, mmcols, drop = FALSE]), 1, min, na.rm = TRUE)
    full <- full[order(full$consensus != "Y", -min_mm), , drop = FALSE]
  }
  list(
    table = .wgcna_hub_symbol_table(utils::head(full$feature, ntop), pgx),
    n_consensus = n_consensus
  )
}

# Build the placeholder values for the consensus module data template.
.wgcna_consensus_module_params <- function(cons, module, pgx) {
  colors <- paste0("ME", cons$net$colors)
  mod_feats <- names(cons$net$colors)[colors == module]
  n_genes <- length(mod_feats)

  cons_cor <- if (module %in% rownames(cons$modTraits)) cons$modTraits[module, ] else numeric(0)
  grp <- names(cons$stats)
  gcor <- lapply(cons$stats, function(s) {
    m <- s$moduleTraitCor
    if (!is.null(m) && module %in% rownames(m)) m[module, ] else numeric(0)
  })

  top_trait <- if (length(cons_cor) && any(!is.na(cons_cor))) {
    names(cons_cor)[which.max(abs(cons_cor))]
  } else {
    NA_character_
  }
  hub <- .wgcna_consensus_hub(cons, module, top_trait, pgx)

  list(
    ME_color = module,
    n_genes = as.character(n_genes),
    groups = paste(grp, collapse = ", "),
    consensus_traits_table = .wgcna_cor_table(
      cons_cor, min_abs = 0.3, ntop = 6L,
      empty = "No trait is consistently associated across all groups."
    ),
    divergent_traits_table = .wgcna_consensus_divergent_table(gcor, grp),
    n_consensus_genes = as.character(hub$n_consensus),
    hub_genes_table = hub$table,
    enrichment_summary = .wgcna_multi_enrichment(cons, module)
  )
}


# -- preservation -------------------------------------------------------------

# Hub table for a preservation module, from the reference-group module
# membership matrix ($stats$moduleMembership), restricted to the module members.
.wgcna_preservation_hub <- function(pres, module, pgx, ntop = 12L) {
  mm <- pres$stats$moduleMembership
  if (is.null(mm)) return("None.")
  col <- if (module %in% colnames(mm)) module else sub("^ME", "", module)
  if (!col %in% colnames(mm)) return("None.")
  color <- sub("^ME", "", module)
  feats <- if (!is.null(pres$colors) && "Consensus" %in% colnames(pres$colors)) {
    rownames(pres$colors)[pres$colors[, "Consensus"] == color]
  } else {
    rownames(mm)
  }
  feats <- intersect(feats, rownames(mm))
  if (length(feats) == 0) feats <- rownames(mm)
  feats <- feats[order(-abs(mm[feats, col]))]
  .wgcna_hub_symbol_table(utils::head(feats, ntop), pgx)
}

# Build the placeholder values for the preservation module data template.
.wgcna_preservation_module_params <- function(pres, module, pgx) {
  ref <- names(pres$datExpr)[1]
  test <- names(pres$datExpr)[2]
  color <- sub("^ME", "", module)

  n_genes <- if (!is.null(pres$moduleSize) && module %in% names(pres$moduleSize)) {
    pres$moduleSize[[module]]
  } else if (!is.null(pres$colors) && "Consensus" %in% colnames(pres$colors)) {
    sum(pres$colors[, "Consensus"] == color, na.rm = TRUE)
  } else {
    NA_integer_
  }

  z <- if (!is.null(pres$Zsummary) && module %in% rownames(pres$Zsummary)) {
    pres$Zsummary[module, 1]
  } else {
    NA_real_
  }
  mr <- if (!is.null(pres$medianRank) && module %in% rownames(pres$medianRank)) {
    pres$medianRank[module, 1]
  } else {
    NA_real_
  }
  n_modules <- if (!is.null(pres$medianRank)) nrow(pres$medianRank) else NA_integer_

  mt <- if (!is.null(pres$stats$moduleTraitCor) &&
            module %in% rownames(pres$stats$moduleTraitCor)) {
    pres$stats$moduleTraitCor[module, ]
  } else {
    numeric(0)
  }

  list(
    ME_color = module,
    n_genes = as.character(n_genes),
    reference = ref,
    test = test,
    zsummary = if (is.na(z)) "n/a" else omicsai::omicsai_format_num(z, 1L),
    zsummary_verbal = .wgcna_zsummary_verbal(z),
    median_rank = if (is.na(mr)) "n/a" else as.character(mr),
    n_modules = as.character(n_modules),
    module_trait_table = .wgcna_cor_table(
      mt, min_abs = 0.3, ntop = 6L,
      empty = "No trait clears |r| >= 0.3 in the reference group."
    ),
    hub_genes_table = .wgcna_preservation_hub(pres, module, pgx),
    enrichment_summary = .wgcna_multi_enrichment(pres, module)
  )
}


# -- enrichment (best-effort, shared) -----------------------------------------

# Consensus/Preservation enrichment is board-computed but stored per shape and
# not always present; surface it when available, otherwise a neutral note.
.wgcna_multi_enrichment <- function(obj, module) {
  gse <- obj$gse
  df <- NULL
  if (is.data.frame(gse) && "module" %in% colnames(gse)) {
    df <- gse[!is.na(gse$module) & gse$module == module, , drop = FALSE]
  } else if (is.list(gse) && !is.null(gse[[module]])) {
    df <- gse[[module]]
  }
  if (!is.data.frame(df) || nrow(df) == 0 || !"geneset" %in% colnames(df)) {
    return("not available for this analysis.")
  }
  if ("q.value" %in% colnames(df)) {
    df <- df[!is.na(df$q.value) & df$q.value < 0.05, , drop = FALSE]
    df <- df[order(df$q.value), , drop = FALSE]
  }
  if (nrow(df) == 0) return("no gene set reached q < 0.05.")
  paste(utils::head(df$geneset, 8), collapse = "; ")
}


# -- prompt builders ----------------------------------------------------------

.wgcna_consensus_prompt <- function(pgx, cons, module) {
  vocab <- omicsai::omicsai_datatype_vocab(pgx$datatype)
  species <- omicsai::omicsai_species_prompt(pgx$organism)
  data_block <- omicsai::frag(
    "wgcna/wgcna_consensus_module_data",
    .wgcna_consensus_module_params(cons, module, pgx)
  )
  .wgcna_summary_build_prompt(data_block, vocab, species,
                              context = "wgcna/wgcna_consensus_interpretation")
}

.wgcna_preservation_prompt <- function(pgx, pres, module) {
  vocab <- omicsai::omicsai_datatype_vocab(pgx$datatype)
  species <- omicsai::omicsai_species_prompt(pgx$organism)
  data_block <- omicsai::frag(
    "wgcna/wgcna_preservation_module_data",
    .wgcna_preservation_module_params(pres, module, pgx)
  )
  .wgcna_summary_build_prompt(data_block, vocab, species,
                              context = "wgcna/wgcna_preservation_interpretation")
}


# -- multi-omics (per layer x module, with cross-omics coordination) ----------

# Map a 2-letter omics layer code to a readable label.
.wgcna_layer_omics_label <- function(code) {
  code <- tolower(code)
  map <- c(gx = "transcriptomics", rna = "transcriptomics", tx = "transcriptomics",
           px = "proteomics", pr = "proteomics", mi = "miRNA", mr = "miRNA",
           me = "methylation", cn = "copy number", mb = "metabolomics",
           mx = "metabolomics", gs = "gene sets", gt = "gene sets",
           ph = "phenotype", cl = "clinical")
  out <- unname(map[code])
  out[is.na(out)] <- code[is.na(out)]
  out
}

# The layer name (e.g. "gx") that owns a given module.
.wgcna_module_layer_code <- function(mox, module) {
  for (ln in names(mox$layers)) {
    if (module %in% names(mox$layers[[ln]]$me.genes)) return(ln)
  }
  NA_character_
}

# Concatenated, sample-aligned eigengene matrix across all layers. Columns are
# the (already layer-prefixed) module names; rows are the shared samples.
.wgcna_mox_eigengene_matrix <- function(mox) {
  mats <- lapply(names(mox$layers), function(ln) {
    layer <- mox$layers[[ln]]
    mes <- if (!is.null(layer$datME)) layer$datME else layer$net$MEs
    if (is.null(mes)) return(NULL)
    pfx <- toupper(ln)
    cn <- colnames(mes)
    if (!all(grepl(paste0("^", pfx), cn) | grepl(paste0("^ME", pfx), cn))) {
      colnames(mes) <- paste0(pfx, sub("^ME", "", cn))
    }
    mes
  })
  mats <- Filter(Negate(is.null), mats)
  if (length(mats) == 0) return(NULL)
  common <- Reduce(intersect, lapply(mats, rownames))
  if (length(common) < 3) return(NULL)
  do.call(cbind, lapply(mats, function(m) m[common, , drop = FALSE]))
}

# Text section listing the modules in OTHER omics layers whose eigengene
# co-varies with `module` (the multi-omics coordination signal). Phenotype and
# gene-set layers are excluded from the partner set.
.wgcna_cross_omics_section <- function(mox, module, ntop = 6L, min_abs = 0.3) {
  ME <- .wgcna_mox_eigengene_matrix(mox)
  if (is.null(ME) || !module %in% colnames(ME)) return("")
  own <- .wgcna_module_layer_code(mox, module)
  drop_layers <- toupper(grep("^(gs|gt|ph|pheno|gset)$", names(mox$layers),
                              value = TRUE, ignore.case = TRUE))
  cmat <- stats::cor(ME, use = "pairwise")
  r <- cmat[module, ]
  r <- r[names(r) != module]
  pfx <- toupper(substr(names(r), 1, 2))
  keep <- pfx != toupper(own) & !pfx %in% drop_layers
  r <- r[keep]
  r <- r[!is.na(r) & abs(r) >= min_abs]
  if (length(r) == 0) {
    return("Cross-omics module coordination: no module in another omics layer correlates at |r| >= 0.3.")
  }
  r <- r[order(-abs(r))]
  r <- utils::head(r, ntop)
  df <- data.frame(
    module = names(r),
    omics = .wgcna_layer_omics_label(substr(names(r), 1, 2)),
    coordination = omicsai::omicsai_verbalize_r(unname(r)),
    stringsAsFactors = FALSE
  )
  paste0(
    "Cross-omics module coordination (eigengene correlation with modules in other omics layers):\n",
    paste(omicsai::omicsai_format_mdtable(df), collapse = "\n")
  )
}

# Data block for one multi-omics module: the standard single-layer block (grey
# aware) plus a layer header and the cross-omics coordination section.
.wgcna_multiomics_data_block <- function(mox, module, pgx, vocab) {
  layer <- wgcna.mox_layer_for_module(mox, module)
  if (is.null(layer)) {
    return(sprintf("Module %s was not found in any omics layer.", module))
  }
  layer <- wgcna.ensureStats(layer)
  if (is.null(layer$annot)) layer$annot <- mox$annot %||% pgx$genes
  code <- .wgcna_module_layer_code(mox, module)
  header <- sprintf("Layer: %s (%s)", .wgcna_layer_omics_label(code), code)
  own_block <- .wgcna_module_data_block(layer, pgx, module, vocab)
  cross <- if (.wgcna_is_grey(module)) "" else .wgcna_cross_omics_section(mox, module)
  paste(c(header, own_block, cross)[nzchar(c(header, own_block, cross))],
        collapse = "\n\n")
}

.wgcna_multiomics_prompt <- function(pgx, mox, module) {
  vocab <- omicsai::omicsai_datatype_vocab(pgx$datatype)
  species <- omicsai::omicsai_species_prompt(pgx$organism)
  data_block <- .wgcna_multiomics_data_block(mox, module, pgx, vocab)
  .wgcna_summary_build_prompt(data_block, vocab, species,
                              context = "wgcna/wgcna_multiomics_interpretation")
}

# Compute-time durable summaries for a multi-omics object: one entry per module
# (all layers), each cross-omics aware. Preserves existing/edited entries unless
# ai$force. Keyed by the layer-prefixed, globally-unique module name.
.wgcna_mox_summarize <- function(pgx, mox, ai, existing = NULL) {
  modules <- names(mox$me.genes)
  if (length(modules) == 0) {
    modules <- unique(unlist(lapply(mox$layers, function(l) names(l$me.genes))))
  }
  out <- existing %||% list()
  todo <- if (isTRUE(ai$force)) modules else modules[!modules %in% names(out)]
  for (mod in todo) {
    out[[mod]] <- tryCatch(
      .wgcna_summary_entry(
        .ai_report_run_prompt(.wgcna_multiomics_prompt(pgx, mox, mod), ai), ai
      ),
      error = function(e) {
        msg <- paste0("[pgx.update_wgcna_summaries] mox module '", mod,
                      "' failed: ", conditionMessage(e))
        if (identical(ai$on_error, "abort")) stop(msg, call. = FALSE)
        if (identical(ai$on_error, "warn")) warning(msg, call. = FALSE)
        else info(msg)
        out[[mod]]
      }
    )
  }
  out
}
