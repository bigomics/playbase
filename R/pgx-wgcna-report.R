## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.

# =============================================================================
# WGCNA AI-report extraction module
# =============================================================================
# Helpers consumed by the omicsplayground board.wgcna AI-report layer plus
# multi-omics / consensus variants. This module owns all "describe the
# computed WGCNA result" logic; pgx-wgcna.R retains only the compute side
# (wgcna.compute, wgcna.computeGeneStats, plot helpers, etc).
#
# Migration history (epic playbase-fad):
#   - Original symbols lived in pgx-wgcna.R lines 5599-5989 alongside the
#     compute code. Moved here verbatim in the consolidation step.
#   - `wgcna.selectTopModules` from the initial duplication step has been
#     absorbed into the new `wgcna.getTopModules` (with `min_modules` floor).

# -----------------------------------------------------------------------------
# Module-trait correlation
# -----------------------------------------------------------------------------

#' Module-trait correlation matrix
#'
#' Returns the cached `wgcna$modTraits` matrix when present, otherwise computes
#' the Pearson correlation between module eigengenes (`wgcna$net$MEs`) and the
#' trait matrix (`wgcna$datTraits`). NA cells are filled with zero.
#'
#' @param wgcna WGCNA result object.
#' @return Numeric matrix (modules × traits).
#' @export
wgcna.get_modTraits <- function(wgcna) {
  if(!is.null(wgcna$modTraits)) {
    M <- wgcna$modTraits
  } else {
    M <- cor( wgcna$net$MEs, wgcna$datTraits, use="pairwise")
  }
  M[is.na(M)] <- 0
  return(M)
}


# -----------------------------------------------------------------------------
# Top genes / sets / phenotypes (single-omics, multi-omics, consensus)
# -----------------------------------------------------------------------------

#' Top genes, sets, and phenotypes per module
#'
#' Dispatch is shape-based: consensus objects (with `$layers` and a list-valued
#' `$datExpr`) go to `wgcna.getConsensusTopGenesAndSets()`; multi-omics objects
#' (with `$layers` but no `$datExpr`) go to `wgcna.getMultiTopGenesAndSets()`;
#' everything else is treated as single-omics.
#'
#' @param wgcna  WGCNA result object.
#' @param annot  Annotation table; optional but recommended.
#' @param module Optional module subset.
#' @param ntop   Integer; top N per module.
#' @param psig   Numeric; p-value cutoff for moduleMembership significance.
#' @param level  "gene" or "geneset".
#' @param rename Annotation column to use when renaming features.
#' @return list(sets, genes, pheno, neg.pheno).
#' @export
wgcna.getTopGenesAndSets <- function(wgcna, annot=NULL, module=NULL, ntop=40,
                                     psig = 0.05, level="gene", rename="symbol") {

  is.consensus <- "layers" %in% names(wgcna) && class(wgcna$datExpr) == "list"
  is.multi <- "layers" %in% names(wgcna) && is.null(wgcna$datExpr)

  if(is.consensus) {
    top <- wgcna.getConsensusTopGenesAndSets(wgcna, annot=annot,
      module=module,  ntop=ntop, rename=rename)
    return(top)
  }

  if(is.multi) {
    top <- wgcna.getMultiTopGenesAndSets(
      wgcna$layers, annot=annot, module=module, psig=psig, ntop=ntop,
      level=level, rename=rename)
    return(top)
  }

  stats <- NULL
  if(!"stats" %in% names(wgcna)) {
    stats <- wgcna.computeGeneStats(wgcna$net, wgcna$datExpr, wgcna$datTraits,
      wgcna$svTOM)
  } else {
    stats <- wgcna$stats
  }
  if(!any(c("gse","gsea") %in% names(wgcna))) {
    warning("object has no enrichment results (gsea)")
  }

  ## get top genes by centrality-weighted-meanFC2
  mm <- stats$moduleMembership
  mm.sig <- 1*(stats$MMPvalue <= psig)
  ff <- sqrt(rowMeans(stats$foldChange**2, na.rm=TRUE))
  mm <- mm * mm.sig * ff
  if(!is.null(annot)) {
    annot$gene_title <- paste0(annot$gene_title," (",annot$symbol,")")
    mm <- rename_by2(mm, annot, new_id=rename)
  }
  gg <- rownames(mm)
  mm <- as.list(data.frame(mm))
  if(!is.null(module)) mm <- mm[which(names(mm) %in% module)]
  for(i in 1:length(mm)) names(mm[[i]]) <- gg
  mm <- lapply(mm, function(x) x[x!=0])
  topgenes <- lapply(mm, function(x) names(head(sort(-x),ntop)))

  ## top genesets
  topsets <- NULL
  if(any(c("gse","gsea") %in% names(wgcna))) {
    if(!is.null(wgcna$gsea)) ee <- wgcna$gsea
    if(!is.null(wgcna$gse)) ee <- wgcna$gse
    if(!is.null(module)) ee <- ee[which(names(ee) %in% module)]
    topsets <- lapply(ee,function(x) head(rownames(x),ntop))
  }

  ## top correlated phenotypes
  M <- wgcna.get_modTraits(wgcna)
  top.pheno <- apply(M, 1, function(x) names(which(x > 0.8*max(x, na.rm=TRUE))))
  top.negpheno <- apply(M, 1, function(x) names(which(x < 0.8*min(x, na.rm=TRUE))))

  if(level=="geneset") {
    topsets <- topgenes
    topgenes <- NULL
  }

  list(sets = topsets, genes = topgenes, pheno = top.pheno, neg.pheno = top.negpheno)
}

#' Top genes/sets/phenotypes for a multi-omics WGCNA object (layered).
#' @export
wgcna.getMultiTopGenesAndSets <- function(multi_wgcna, annot=NULL, module=NULL,
                                          psig=0.05, ntop=40, level=NULL,
                                          rename="symbol") {

  if("layers" %in% names(multi_wgcna)) {
    multi_wgcna <- multi_wgcna$layers
  }

  ## set level
  nw <- length(multi_wgcna)
  if(!is.null(level)) {
    level <- head(rep(level, nw),nw)
  } else {
    level <- c("gene","geneset")[1 + 1*grepl("^gs|^gset|geneset",names(multi_wgcna))]
  }
  names(level) <- names(multi_wgcna)

  toplist <- list()
  k <- names(multi_wgcna)[1]
  for (k in names(multi_wgcna)) {
    topk <- wgcna.getTopGenesAndSets(
      multi_wgcna[[k]],  module=module,  annot=annot,
      ntop=ntop, psig=psig, level=level[[k]], rename=rename)
    if(!is.null(module)) {
      topk <- lapply( topk, function(s) s[which(names(s) %in% module)] )
    }
    toplist[[k]] <- topk
  }

  top <- list()
  top$genes <- lapply(toplist, function(t) t[['genes']])
  names(top$genes) <- NULL
  top$genes <- unlist(top$genes, recursive = FALSE)

  top$sets <- lapply(toplist, function(t) t[["sets"]])
  names(top$sets) <- NULL
  top$sets <- unlist(top$sets, recursive = FALSE)

  top$pheno <- lapply(toplist, function(t) t[["pheno"]])
  names(top$pheno) <- NULL
  top$pheno <- unlist(top$pheno, recursive = FALSE)

  top$neg.pheno <- lapply(toplist, function(t) t[["neg.pheno"]])
  names(top$neg.pheno) <- NULL
  top$neg.pheno <- unlist(top$neg.pheno, recursive = FALSE)

  return(top)
}


#' Top genes/sets/phenotypes for a consensus WGCNA object.
#' @export
wgcna.getConsensusTopGenesAndSets <- function(cons, annot=NULL, module=NULL, ntop=40,
                                              level=c("gene","geneset")[1],
                                              rename="symbol" ) {
  if(!"stats" %in% names(cons)) stop("object has no stats")
  if(!any(c("gse","gsea") %in% names(cons))) {
    warning("object has no enrichment results (gsea)")
  }

  if(!is.null(annot)) {
    annot$gene_title <- paste0(annot$gene_title," (",annot$symbol,")")
  }

  ## get top genes (highest kME)
  topgenesx <- list()
  for(i in 1:length(cons$stats)) {
    mm <- cons$stats[[i]]$moduleMembership
    if(!is.null(annot)) {
      mm <- rename_by2(mm, annot, rename)
    }
    gg <- rownames(mm)
    mm <- as.list(data.frame(mm))
    if (!is.null(module)) mm <- mm[module]
    sel.topgenes <- lapply(mm, function(x) head(order(-x), 3 * ntop))
    topgenesx[[i]] <- lapply(sel.topgenes, function(i) gg[i])
  }

  ## intersect topgenes across all datatypes
  topgenes <- topgenesx[[1]]
  k <- 2
  for (k in 2:length(topgenesx)) {
    topgenes <- mapply(intersect, topgenes, topgenesx[[k]], SIMPLIFY = FALSE)
  }
  topgenes <- lapply(topgenes, head, ntop)

  if(!is.null(module)) {
    sel <- intersect(names(topgenes),module)
    topgenes <- topgenes[sel]
  }

  ## top genesets (as symbol!)
  topsets <- NULL
  if(any(c("gse","gsea") %in% names(cons))) {
    if(!is.null(cons$gsea)) ee <- cons$gsea
    if(!is.null(cons$gse)) ee <- cons$gse
    ee <- ee[match(names(topgenes),names(ee))]
    names(ee) <- names(topgenes)
    topsets <- lapply(ee,function(x) head(rownames(x),ntop))
  }

  ## module traits
  M <- lapply(cons$net$multiMEs, function(x) as.matrix(x$data))
  Y <- lapply(M, function(m) cons$datTraits[rownames(m),])
  R <- mapply( function(x,y) abs(cor(x,y,use="pairwise")), M, Y, SIMPLIFY=FALSE)
  R <- Reduce('+', R)
  top.pheno <- apply(R, 1, function(x) names(which(x > 0.8*max(x,na.rm=TRUE))),
    simplify = FALSE)
  top.negpheno <- apply(R, 1, function(x) names(which(x < 0.8*min(x,na.rm=TRUE))),
    simplify = FALSE)

  if (level == "geneset") {
    topsets <- topgenes
    topgenes <- NULL
  }

  list(sets = topsets, genes = topgenes, pheno = top.pheno, neg.pheno = top.negpheno)
}


# -----------------------------------------------------------------------------
# Module selection (with minimum-count floor)
# -----------------------------------------------------------------------------

#' Select top modules with a minimum-count floor
#'
#' A module qualifies if (a) it has at least one trait p-value below `psig`, or
#' (b) its row in the module-trait correlation matrix has at least one column
#' meeting `topratio * max(R)` (clipped at `minrho`).
#'
#' If the resulting selection has fewer than `min_modules` non-grey modules,
#' the selection is padded with the next-strongest non-grey modules ranked by
#' `rowMeans(abs(R))`. If the input has `<= min_modules` non-grey modules in
#' total, all non-grey modules are returned unfiltered.
#'
#' Replaces the older `wgcna.getTopModules()` variant whose `kx` argument was
#' deprecated; the new signature drops `kx` entirely and adds `min_modules`.
#'
#' @param wgcna       WGCNA result object (single-omics or layered).
#' @param min_modules Integer; minimum number of modules to return.
#' @param topratio    Numeric.
#' @param psig        Numeric.
#' @param minrho      Numeric.
#' @param rm.grey     Logical; drop grey from the selection.
#' @param multi       Logical or NULL; force multi-layer dispatch.
#' @return Character vector of module names.
#' @export
wgcna.getTopModules <- function(wgcna, min_modules = 5L,
                                topratio = 0.85, psig = 0.05,
                                minrho = 0.1, rm.grey = TRUE,
                                multi = NULL) {

  ## --- Resolve the candidate non-grey module pool from me.genes -------------
  ## Multi-layer objects may not carry me.genes at the top level; skip the
  ## bypass/floor in that case and let the per-layer selection below drive.
  raw <- names(wgcna$me.genes)
  non_grey <- character(0)
  if (length(raw) > 0) {
    suf <- suppressWarnings(as.integer(sub("^ME", "", raw)))
    display <- raw
    if (length(suf) > 0 && all(!is.na(suf))) {
      display <- paste0("ME", WGCNA::labels2colors(suf))
    }
    is_grey <- raw %in% c("MEgrey", "grey") | display %in% c("MEgrey", "grey")
    non_grey <- raw[!is_grey]

    ## Bypass: too few non-grey modules to filter meaningfully — return all.
    if (length(non_grey) <= min_modules) {
      return(non_grey)
    }
  }

  ## --- Primary selection (significance + topratio gating, per-layer) --------
  if (is.null(topratio)) topratio <- 0.85
  if (is.null(multi) && !is.null(wgcna$layers)) multi <- TRUE
  if (is.null(multi)) multi <- FALSE
  if (!multi) {
    ww <- list(gx = wgcna)
  } else if (!is.null(wgcna$layers)) {
    ww <- wgcna$layers
  } else {
    ww <- wgcna
  }

  R <- list()
  P <- list()
  for (i in seq_along(ww)) {
    me <- ww[[i]]$net$MEs
    dt <- ww[[i]]$datTraits
    R1 <- cor(me, dt, use = "pairwise")
    ndim <- colSums(!is.na(dt))
    P1 <- sapply(seq_len(ncol(dt)), function(j) cor.pvalue(R1[, j], ndim[j]))
    colnames(P1) <- colnames(dt)
    R[[i]] <- R1
    P[[i]] <- P1
  }

  selected <- c()
  for (i in seq_along(R)) {
    idx1 <- which(rowSums(P[[i]] <= psig) > 0)
    rmax <- topratio * pmax(apply(R[[i]], 2, max, na.rm = TRUE), 0)
    rmax <- pmax(rmax, minrho)
    idx2 <- which(colSums(t(R[[i]]) >= rmax) > 0)
    idx <- setdiff(unique(c(idx1, idx2)), 0)
    selected <- c(selected, rownames(R[[i]])[idx])
  }

  if (rm.grey) {
    sel.grey <- grepl("[A-Z]{2}grey$", selected)
    selected <- selected[!sel.grey]
  }
  selected <- as.character(selected)

  ## --- Floor padding --------------------------------------------------------
  if (length(selected) >= min_modules) {
    return(selected)
  }

  me <- wgcna$net$MEs
  dt <- wgcna$datTraits
  if (!is.data.frame(me)) me <- tryCatch(as.data.frame(me), error = function(e) NULL)
  if (!is.data.frame(dt)) dt <- tryCatch(as.data.frame(dt), error = function(e) NULL)
  if (is.null(me) || is.null(dt) || ncol(me) == 0 || ncol(dt) == 0) {
    message("[wgcna.getTopModules] cannot pad to floor — module-trait ",
            "matrix unavailable. Returning ", length(selected), " modules.")
    return(selected)
  }

  Rpad <- cor(me, dt, use = "pairwise")
  mx <- rowMeans(abs(Rpad), na.rm = TRUE)
  mx <- mx[intersect(names(mx), non_grey)]
  mx <- sort(mx, decreasing = TRUE)

  need <- min_modules - length(selected)
  candidates <- setdiff(names(mx), selected)
  c(selected, head(candidates, need))
}


# -----------------------------------------------------------------------------
# Ensure stats are populated (lazy fill for older PGX files)
# -----------------------------------------------------------------------------

#' Ensure a WGCNA object has the fields needed for downstream extraction
#'
#' Older PGX files sometimes ship without `wgcna$stats` or with empty
#' `wgcna$net$labels`, which breaks `wgcna.getGeneStats()`. This helper
#' lazily populates both. Idempotent.
#'
#' @param wgcna WGCNA result object.
#' @return The (possibly mutated) `wgcna` object.
#' @export
wgcna.ensureStats <- function(wgcna) {
  if (!is.null(wgcna$layers) && is.null(wgcna$datExpr)) {
    return(wgcna)
  }

  needs_labels <- is.null(wgcna$net$labels) || length(wgcna$net$labels) == 0
  needs_stats  <- is.null(wgcna$stats) || length(wgcna$stats) == 0

  if (needs_labels && !is.null(wgcna$net$colors)) {
    ## wgcna.getGeneStats() expects ME-prefixed labels regardless of whether
    ## the colors are integer-shape (example-data) or color-name shape.
    raw <- as.character(wgcna$net$colors)
    needs_prefix <- !grepl("^ME", raw)
    raw[needs_prefix] <- paste0("ME", raw[needs_prefix])
    wgcna$net$labels <- raw
  }

  if (needs_stats) {
    if (is.null(wgcna$net) || is.null(wgcna$datExpr) || is.null(wgcna$datTraits)) {
      message("[wgcna.ensureStats] cannot compute stats — net/datExpr/datTraits missing")
      return(wgcna)
    }
    tom <- if (!is.null(wgcna$TOM)) wgcna$TOM
           else if (!is.null(wgcna$svTOM)) wgcna$svTOM
           else if (!is.null(wgcna$wTOM)) wgcna$wTOM
           else NULL
    s <- tryCatch(
      wgcna.computeGeneStats(
        net       = wgcna$net,
        datExpr   = wgcna$datExpr,
        datTraits = wgcna$datTraits,
        TOM       = tom
      ),
      error = function(e) {
        message("[wgcna.ensureStats] wgcna.computeGeneStats failed: ",
                conditionMessage(e))
        NULL
      }
    )
    if (!is.null(s)) wgcna$stats <- s
  }

  wgcna
}

.wgcna_mod_traits <- function(wgcna) {
  if (!is.null(wgcna$layers) && is.null(wgcna$datExpr)) {
    mats <- lapply(wgcna$layers, function(layer) {
      tryCatch(playbase:::wgcna.get_modTraits(layer), error = function(e) NULL)
    })
    mats <- Filter(function(x) is.matrix(x) || is.data.frame(x), mats)
    if (!length(mats)) return(NULL)
    cols <- Reduce(union, lapply(mats, colnames))
    mats <- lapply(mats, function(x) {
      x <- as.matrix(x)
      out <- matrix(NA_real_, nrow = nrow(x), ncol = length(cols),
                    dimnames = list(rownames(x), cols))
      out[, colnames(x)] <- x
      out
    })
    return(do.call(rbind, mats))
  }
  tryCatch(playbase:::wgcna.get_modTraits(wgcna), error = function(e) NULL)
}

# =============================================================================
# AI-report builders (round 2: ported from v03 components/board.wgcna/R/
# ai.report/ai_report_data.R via cp; surgical edits noted below).
#   - BOARD_PROMPTS_DIR replaced by .wgcna_prompt_path()
#   - top-level MODULE_SUMMARY load deferred to lazy accessor
# All other symbols verbatim from v03.
# =============================================================================

.wgcna_prompt_path <- function(name) {
  .ai_report_prompt_path("wgcna", name)
}

.wgcna_module_summary_template <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      cache <<- omicsai::omicsai_load_template(.wgcna_prompt_path("wgcna_module_data.md"))
    }
    cache
  }
})


# -----------------------------------------------------------------------------
# Extraction primitives
# -----------------------------------------------------------------------------

#' Look up function descriptions for features (NULL-safe)
resolve_functions <- function(features, annot, max_chars = 60L) {
  funcs <- rep("", length(features))
  if (is.null(annot)) return(funcs)
  func_col <- intersect(c("gene_title", "gene_name", "description"), colnames(annot))
  if (length(func_col) == 0) return(funcs)
  idx <- match(features, rownames(annot))
  valid <- !is.na(idx)
  funcs[valid] <- as.character(annot[idx[valid], func_col[1]])
  substr(funcs, 1, max_chars)
}

#' Render module identifiers in canonical MEcolor form.
#'
#' Datasets keyed by integer suffix (`ME0`/`ME1`/...) get remapped via
#' `WGCNA::labels2colors`; datasets already in MEcolor form pass through.
#' Identifiers not present in `wgcna$me.genes` are returned unchanged.
display_module <- function(x, wgcna) {
  me_names <- names(wgcna$me.genes)
  suf <- suppressWarnings(as.integer(sub("^ME", "", me_names)))
  if (length(suf) == 0 || any(is.na(suf))) return(x)
  map <- setNames(paste0("ME", WGCNA::labels2colors(suf)), me_names)
  out <- map[x]
  out[is.na(out)] <- x[is.na(out)]
  unname(out)
}


#' Extract all data for a single WGCNA module.
#'
#' Returns a list with raw values used by both the summary mode and the
#' multi-module report mode. The caller is expected to cache `top` and `M`
#' across modules.
#'
#' Handles two `wgcna$gse` shapes: per-module list (mox-brca, Bruker_Generoso)
#' and flat data.frame with a `module` column (example-data and similar).
extract_module_data <- function(wgcna, module, pgx,
                                top = NULL, M = NULL,
                                max_traits = 3L) {
  annot <- wgcna$annot

  if (is.null(top)) {
    top <- tryCatch(
      playbase::wgcna.getTopGenesAndSets(wgcna, annot = annot, ntop = 40,
                                         level = "gene", rename = "gene_title"),
      error = function(e) list(pheno = list(), genes = list(), sets = list())
    )
  }

  phenotypes <- "None"
  if (module %in% names(top$pheno)) {
    phenotypes <- paste(top$pheno[[module]], collapse = ", ")
  }
  fallback_genes <- top$genes[[module]]
  if (is.null(fallback_genes)) fallback_genes <- character(0)
  fallback_sets <- top$sets[[module]]
  if (is.null(fallback_sets)) fallback_sets <- character(0)

  size <- length(wgcna$me.genes[[module]])

  if (is.null(M)) M <- .wgcna_mod_traits(wgcna)

  top_trait <- ""
  top_r <- NA_real_
  top_pos_trait <- NA_character_
  top_pos_r <- NA_real_
  top_neg_trait <- NA_character_
  top_neg_r <- NA_real_
  if (!is.null(M) && module %in% rownames(M)) {
    cors <- M[module, ]
    top_idx <- which.max(abs(cors))
    if (length(top_idx) > 0) {
      top_trait <- names(cors)[top_idx]
      top_r <- cors[top_idx]
    }
    ## Dual-trait split: best positive + best negative trait, each gated by
    ## |r| >= 0.5 so direction cannot be inferred from a single absolute-value
    ## pick (the silent contradiction we hit in v01).
    pos_idx <- which.max(cors)
    if (length(pos_idx) > 0 && !is.na(cors[pos_idx]) && cors[pos_idx] >= 0.5) {
      top_pos_trait <- names(cors)[pos_idx]
      top_pos_r     <- cors[pos_idx]
    }
    neg_idx <- which.min(cors)
    if (length(neg_idx) > 0 && !is.na(cors[neg_idx]) && cors[neg_idx] <= -0.5) {
      top_neg_trait <- names(cors)[neg_idx]
      top_neg_r     <- cors[neg_idx]
    }
  }

  traits <- character(0)
  if (phenotypes != "None" && nzchar(phenotypes)) {
    traits <- trimws(strsplit(phenotypes, ",\\s*")[[1]])
    traits <- traits[nzchar(traits) & traits != "None"]
    traits <- head(traits, max_traits)
  }

  eigengene_profile <- NULL
  datME <- if (!is.null(wgcna$datME)) wgcna$datME else wgcna$net$MEs
  groups <- if (!is.null(pgx$samples$group)) pgx$samples$group else NULL
  if (!is.null(datME) && module %in% colnames(datME) && !is.null(groups)) {
    eigengene_profile <- tapply(datME[, module], groups, mean)
  }

  ## Enrichment: resolve both gse shapes (per-module list vs flat df).
  n_sig <- 0L
  n_total <- 0L
  gse_raw <- wgcna$gse
  gse <- NULL
  if (is.data.frame(gse_raw) && "module" %in% colnames(gse_raw)) {
    gse <- gse_raw[!is.na(gse_raw$module) & gse_raw$module == module, , drop = FALSE]
  } else if (is.list(gse_raw)) {
    gse <- gse_raw[[module]]
  }
  if (is.data.frame(gse) && nrow(gse) > 0 && "q.value" %in% colnames(gse)) {
    n_sig <- sum(gse$q.value < 0.05, na.rm = TRUE)
    n_total <- nrow(gse)
  }

  trait_stats <- list()
  stat_traits <- if (length(traits) > 0) traits else if (nzchar(top_trait)) top_trait
  for (tr in stat_traits) {
    gs <- tryCatch(
      playbase::wgcna.getGeneStats(wgcna, module = module, trait = tr, plot = FALSE),
      error = function(e) NULL
    )
    if (is.data.frame(gs) && nrow(gs) > 0) {
      trait_stats[[tr]] <- gs
    }
  }

  list(
    size = size,
    phenotypes = phenotypes,
    traits = traits,
    fallback_genes = fallback_genes,
    fallback_sets = fallback_sets,
    eigengene_profile = eigengene_profile,
    top_trait = top_trait,
    top_r = top_r,
    top_pos_trait = top_pos_trait,
    top_pos_r = top_pos_r,
    top_neg_trait = top_neg_trait,
    top_neg_r = top_neg_r,
    n_sig = n_sig,
    n_total = n_total,
    gse = gse,
    trait_stats = trait_stats
  )
}


#' Compute per-module structured data for a set of modules.
#'
#' Wraps extract_module_data() over a list of modules, attaching hub-gene
#' tables (resolved symbols + descriptions) and the trait selected for the
#' single-trait report view.
.compute_module_data <- function(wgcna, pgx, modules,
                                 ntop_enrichment = 20L, ntop_genes = 50L) {
  annot <- if (!is.null(wgcna$annot)) wgcna$annot else pgx$genes

  M <- .wgcna_mod_traits(wgcna)
  top <- tryCatch(
    playbase::wgcna.getTopGenesAndSets(wgcna, annot = annot, ntop = 40,
                                       level = "gene", rename = "gene_title"),
    error = function(e) list(pheno = list(), genes = list(), sets = list())
  )

  module_data <- list()
  for (mod in modules) {
    md <- extract_module_data(wgcna, mod, pgx, top = top, M = M)

    ## Hub genes for the single-trait report view (first phenotype-derived
    ## trait if available, else the absolute-top trait).
    hub_genes_df <- NULL
    trait_for_genes <- if (length(md$traits) > 0) md$traits[1] else md$top_trait

    gs <- if (nzchar(trait_for_genes) && trait_for_genes %in% names(md$trait_stats)) {
      md$trait_stats[[trait_for_genes]]
    } else {
      NULL
    }

    if (is.data.frame(gs) && nrow(gs) > 0) {
      if ("moduleMembership" %in% colnames(gs)) {
        gs <- gs[order(-abs(gs$moduleMembership)), , drop = FALSE]
      }
      top_gs <- head(gs, ntop_genes)

      symbols <- playbase::probe2symbol(top_gs$feature, annot, "symbol", fill_na = TRUE)
      funcs <- resolve_functions(top_gs$feature, annot)

      mm_vals <- if ("moduleMembership" %in% colnames(top_gs)) {
        top_gs$moduleMembership
      } else {
        rep(NA_real_, nrow(top_gs))
      }
      logfc_vals <- if ("foldChange" %in% colnames(top_gs)) {
        top_gs$foldChange
      } else {
        rep(NA_real_, nrow(top_gs))
      }

      hub_genes_df <- data.frame(
        symbol = symbols,
        MM = mm_vals,
        logFC = logfc_vals,
        func = funcs,
        stringsAsFactors = FALSE
      )
    }

    module_data[[mod]] <- list(
      size = md$size,
      eigengene_profile = md$eigengene_profile,
      top_trait = md$top_trait,
      top_r = md$top_r,
      top_pos_trait = md$top_pos_trait,
      top_pos_r = md$top_pos_r,
      top_neg_trait = md$top_neg_trait,
      top_neg_r = md$top_neg_r,
      n_sig = md$n_sig,
      n_total = md$n_total,
      full_gse = md$gse,
      hub_genes = hub_genes_df,
      trait_for_genes = trait_for_genes
    )
  }

  module_data
}


#' Per-module gene-family enrichment (≥3 members in a family of 5–500).
#' Returns a named list of character vectors, keyed by module.
.compute_families_text <- function(wgcna, pgx, modules) {
  if (is.null(pgx$families)) return(NULL)
  fam_names <- names(pgx$families)
  fam_names <- fam_names[fam_names != "<all>"]
  fam_sizes <- lengths(pgx$families[fam_names])
  fam_names <- fam_names[fam_sizes >= 5 & fam_sizes <= 500]
  if (length(fam_names) == 0) return(NULL)

  out <- list()
  for (mod in modules) {
    mod_genes <- wgcna$me.genes[[mod]]
    if (is.null(mod_genes) || length(mod_genes) == 0) next

    overlaps <- vapply(fam_names, function(fn) {
      length(intersect(mod_genes, pgx$families[[fn]]))
    }, integer(1))

    sig_fam <- overlaps[overlaps >= 3]
    if (length(sig_fam) == 0) next
    sig_fam <- head(sort(sig_fam, decreasing = TRUE), 5)

    out[[mod]] <- vapply(names(sig_fam), function(fn) {
      sprintf("%s (%d of %d)", fn, sig_fam[fn], length(pgx$families[[fn]]))
    }, character(1))
  }
  out
}


# -----------------------------------------------------------------------------
# Section builders — one per template slot
# -----------------------------------------------------------------------------

#' Build the `## Overview` placeholder values.
#'
#' Returns a list of named substitutions for the Overview block of the
#' data template. Caller merges these into the substitute() arg list.
data_overview <- function(wgcna, pgx) {
  info <- .ai_report_get(pgx, "info", override = wgcna$experiment)
  n_features_used <- length(unlist(wgcna$me.genes))
  if (is.null(n_features_used) || n_features_used == 0) {
    n_features_used <- tryCatch(ncol(wgcna$datExpr), error = function(e) NA_integer_)
  }

  power <- if (!is.null(wgcna$power)) wgcna$power
           else if (!is.null(wgcna$net$power)) wgcna$net$power
           else "NA"
  min_mod_size <- if (!is.null(wgcna$minModSize)) wgcna$minModSize else "20"
  merge_cut_height <- if (!is.null(wgcna$mergeCutHeight)) wgcna$mergeCutHeight else "0.15"

  list(
    experiment       = info$experiment,
    organism         = info$organism,
    n_samples        = info$n_samples,
    n_features_total = info$n_features,
    n_features_used  = as.character(n_features_used),
    power            = as.character(power),
    min_mod_size     = as.character(min_mod_size),
    merge_cut_height = as.character(merge_cut_height)
  )
}


#' Build the `## Experimental contrasts` block.
#'
#' Returns a markdown bullet list (one line per contrast) or a placeholder
#' string if no contrast matrix is available.
data_contrast <- function(pgx) {
  .ai_report_get(pgx, "contrast_matrix_block")
}


#' Build the `## Modules summary` table block.
#'
#' Dual-trait columns: best positive (>= 0.5) and best negative (<= -0.5).
#' Lead module (first row) also carries `(r = ±0.NN)` next to the verbal
#' label on its top correlated trait. Grey is appended as a final row.
data_modsummary <- function(wgcna, module_data, lead_module) {
  module_order <- names(module_data)
  if (length(module_order) == 0) return(list(table = "", footnote = ""))

  fmt_trait_col <- function(trait, r, is_lead_pos = FALSE) {
    if (is.na(trait) || !nzchar(trait)) return("—")
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
      Module = display_module(m, wgcna),
      Genes  = as.character(md$size),
      `Top correlated trait`      = fmt_trait_col(md$top_pos_trait, md$top_pos_r, is_lead),
      `Top anti-correlated trait` = fmt_trait_col(md$top_neg_trait, md$top_neg_r),
      `Enrichment hits`           = as.character(md$n_sig),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })

  grey_size <- length(wgcna$me.genes[["MEgrey"]])
  if (is.null(grey_size) || grey_size == 0) {
    grey_size <- length(wgcna$me.genes[["grey"]])
  }
  if (!is.null(grey_size) && grey_size > 0) {
    ov_rows[[length(ov_rows) + 1]] <- data.frame(
      Module = "grey", Genes = as.character(grey_size),
      `Top correlated trait` = "—",
      `Top anti-correlated trait` = "—",
      `Enrichment hits` = "—",
      stringsAsFactors = FALSE, check.names = FALSE
    )
  }

  overview_df <- do.call(rbind, ov_rows)
  table_md <- paste(omicsai::omicsai_format_mdtable(overview_df), collapse = "\n")

  ## Footnote: "Showing X of N modules."
  raw <- names(wgcna$me.genes)
  suf <- suppressWarnings(as.integer(sub("^ME", "", raw)))
  display <- raw
  if (length(suf) > 0 && all(!is.na(suf))) {
    display <- paste0("ME", WGCNA::labels2colors(suf))
  }
  is_grey <- raw %in% c("MEgrey", "grey") | display %in% c("MEgrey", "grey")
  n_total_non_grey <- sum(!is_grey)
  footnote <- sprintf("Showing %d of %d modules.",
                      length(module_order), n_total_non_grey)

  list(table = table_md, footnote = footnote)
}


#' Build the `## Module-module eigengene correlations` block.
#'
#' Verbal labels only; raw r stripped at the rendering layer.
data_eigen_cor <- function(wgcna, modules) {
  me_data <- if (!is.null(wgcna$datME)) wgcna$datME else wgcna$net$MEs
  if (is.null(me_data)) return("(no strong correlations detected)")

  me_cols <- intersect(modules, colnames(me_data))
  if (length(me_cols) < 2) return("(no strong correlations detected)")

  me_cor <- cor(me_data[, me_cols])
  pairs <- character(0)
  for (i in seq_len(length(me_cols) - 1)) {
    for (j in (i + 1):length(me_cols)) {
      r <- me_cor[i, j]
      if (!is.na(r) && abs(r) >= 0.7) {
        pairs <- c(pairs, sprintf("%s ↔ %s: %s",
                                  display_module(me_cols[i], wgcna),
                                  display_module(me_cols[j], wgcna),
                                  omicsai::omicsai_verbalize_r(r)))
      }
    }
  }
  if (length(pairs) == 0) return("(no strong correlations detected)")
  paste(pairs, collapse = "\n")
}


#' Build per-module detail blocks (v03 layout).
#'
#' For each module: heading, eigengene profile, dual-trait coordination,
#' hub genes (BEFORE enrichment), enrichment top-N with Hub overlap column,
#' optional gene-family enrichment.
data_module_detail <- function(wgcna, pgx, module_data,
                               families_text = NULL,
                               ntop_enrichment = 20L,
                               ntop_genes = 10L) {
  module_order <- names(module_data)
  if (length(module_order) == 0) return("")

  blocks <- vapply(module_order, function(mod) {
    .render_module_block(module_data[[mod]], mod, wgcna,
                         families_text   = families_text,
                         ntop_enrichment = ntop_enrichment,
                         ntop_genes      = ntop_genes)
  }, character(1))

  paste(blocks, collapse = "\n\n")
}


# -----------------------------------------------------------------------------
# Orchestrator
# -----------------------------------------------------------------------------

#' Build structured report tables from WGCNA results.
#'
#' Renders the v03 data block: a single markdown document substituted into
#' `prompts/wgcna_report_data.md`, plus a structured `data` list for any
#' downstream callers that want the raw values.
#'
#' @return list(text = character, data = list)
wgcna_build_report_tables <- function(wgcna, pgx,
                                      n_modules = 8L,
                                      ntop_enrichment = 20L,
                                      ntop_genes = 50L,
                                      include_module_cors = TRUE,
                                      include_families = TRUE,
                                      include_contrasts = TRUE) {
  ## Lazy-fill labels + stats so older PGX files extract cleanly.
  wgcna <- playbase::wgcna.ensureStats(wgcna)

  modules <- playbase::wgcna.getTopModules(wgcna, min_modules = 5L)
  modules <- head(modules, n_modules)

  module_data <- .compute_module_data(wgcna, pgx, modules,
                                      ntop_enrichment = ntop_enrichment,
                                      ntop_genes = ntop_genes)

  ## Order modules by significant-enrichment count, descending.
  module_order <- names(sort(
    vapply(module_data, function(x) x$n_sig, integer(1)),
    decreasing = TRUE
  ))
  module_data <- module_data[module_order]
  lead_module <- if (length(module_order) > 0) module_order[1] else NA_character_

  families_text <- if (include_families) {
    .compute_families_text(wgcna, pgx, module_order)
  } else NULL

  module_cors <- if (include_module_cors) data_eigen_cor(wgcna, module_order)
                 else "(module-module correlations omitted)"

  overview_params <- data_overview(wgcna, pgx)
  modsum <- data_modsummary(wgcna, module_data, lead_module)
  per_module <- data_module_detail(wgcna, pgx, module_data,
                                   families_text = families_text,
                                   ntop_enrichment = ntop_enrichment,
                                   ntop_genes = ntop_genes)

  tmpl <- omicsai::omicsai_load_template(
    .wgcna_prompt_path("wgcna_report_data.md")
  )

  text <- omicsai::omicsai_substitute_template(tmpl, c(
    overview_params,
    list(
      experiment_info            = if (include_contrasts) {
        .ai_report_get(pgx, "experiment_info", override = wgcna$experiment)
      } else {
        .ai_report_get(pgx, "experiment_info", override = wgcna$experiment,
                       n_contrasts = 0L,
                       contrasts_block = "(contrasts omitted)")
      },
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
      modules          = module_data
    )
  )
}




# -- Leaf renderers: one per template placeholder that needs derivation. -----
# Each returns a single string. No prose, only data → string formatting.

.eigengene_profile_str <- function(md) {
  if (is.null(md$eigengene_profile)) return("—")
  paste(sprintf("%s=%s", names(md$eigengene_profile),
                omicsai::omicsai_format_num(md$eigengene_profile, 2)),
        collapse = ", ")
}

.hub_genes_sorted <- function(md, ntop) {
  if (is.null(md$hub_genes) || nrow(md$hub_genes) == 0) return(NULL)
  hg <- md$hub_genes
  mm  <- if ("MM"    %in% colnames(hg)) hg$MM    else rep(NA_real_, nrow(hg))
  lfc <- if ("logFC" %in% colnames(hg)) hg$logFC else rep(NA_real_, nrow(hg))
  head(hg[order(-abs(mm), -abs(lfc), na.last = TRUE), , drop = FALSE], ntop)
}

.hub_genes_table_str <- function(hg) {
  if (is.null(hg) || nrow(hg) == 0) return("—")
  paste(omicsai::omicsai_format_mdtable(data.frame(
    Gene             = paste0("*", hg$symbol, "*"),
    `Known function` = hg$func,
    stringsAsFactors = FALSE, check.names = FALSE
  )), collapse = "\n")
}

.enrichment_table_str <- function(md, hub_symbols, ntop) {
  full_gse <- md$full_gse
  if (!is.data.frame(full_gse) || nrow(full_gse) == 0) return("—")

  ## Sort by overlap (preferred), then q-value, then score.
  full_gse <- if ("overlap" %in% colnames(full_gse)) {
    full_gse[order(full_gse$overlap, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
  } else if ("overlap_count" %in% colnames(full_gse)) {
    full_gse[order(full_gse$overlap_count, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
  } else if ("q.value" %in% colnames(full_gse)) {
    full_gse[order(full_gse$q.value, na.last = TRUE), , drop = FALSE]
  } else if ("score" %in% colnames(full_gse)) {
    full_gse[order(full_gse$score, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
  } else full_gse

  top_n <- head(full_gse, ntop)
  hub_overlap_col <- if (length(hub_symbols) == 0) {
    rep("—", nrow(top_n))
  } else if ("genes" %in% colnames(top_n)) {
    vapply(top_n$genes, function(g) {
      gset <- strsplit(as.character(g), "|", fixed = TRUE)[[1]]
      sprintf("%d/%d", length(intersect(gset, hub_symbols)), length(hub_symbols))
    }, character(1), USE.NAMES = FALSE)
  } else {
    rep(sprintf("0/%d", length(hub_symbols)), nrow(top_n))
  }

  paste(omicsai::omicsai_format_mdtable(data.frame(
    Rank         = seq_len(nrow(top_n)),
    Geneset      = as.character(top_n$geneset),
    `Hub overlap` = hub_overlap_col,
    stringsAsFactors = FALSE, check.names = FALSE
  )), collapse = "\n")
}


#' Render one module's data dict into a MODULE_SUMMARY block.
#'
#' The `md` argument is one element from `.compute_module_data()`'s output.
#' This function builds only the value dict; MODULE_SUMMARY owns all prose.
.render_module_block <- function(md, mod_name, wgcna,
                                 families_text = NULL,
                                 ntop_enrichment = 20L,
                                 ntop_genes = 10L) {
  hg <- .hub_genes_sorted(md, ntop_genes)
  hub_symbols <- if (!is.null(hg)) as.character(hg$symbol) else character(0)

  fam <- if (!is.null(families_text) && !is.null(families_text[[mod_name]])) {
    paste(families_text[[mod_name]], collapse = "; ")
  } else "—"

  overlaps <- if (!is.null(md$enrichment_overlap) && length(md$enrichment_overlap) > 0) {
    paste(md$enrichment_overlap, collapse = "; ")
  } else "—"

  na_dash <- function(x) if (is.na(x) || !nzchar(x)) "—" else x

  omicsai::omicsai_substitute_template(.wgcna_module_summary_template(), list(
    ME_color                      = display_module(mod_name, wgcna),
    n_genes                       = as.character(md$size),
    tier                          = "—",                ## tier classification dropped
    eigengene_profile_qualitative = .eigengene_profile_str(md),
    top_pos_trait                 = na_dash(md$top_pos_trait),
    top_pos_verbal                = omicsai::omicsai_verbalize_r(md$top_pos_r),
    top_neg_trait                 = na_dash(md$top_neg_trait),
    top_neg_verbal                = omicsai::omicsai_verbalize_r(md$top_neg_r),
    n_sig_terms                   = as.character(md$n_sig),
    n_total_terms                 = as.character(md$n_total),
    enrichment_themes_table       = .enrichment_table_str(md, hub_symbols, ntop_enrichment),
    n_hub                         = as.character(if (is.null(hg)) 0L else nrow(hg)),
    hub_genes_table               = .hub_genes_table_str(hg),
    gene_families_summary         = fam,
    enrichment_overlaps           = overlaps
  ))
}

# =============================================================================
# Section 3 — Prompt assembly (omicsai fragments)
# =============================================================================

#' Assemble the WGCNA static-report prompt from omicsai fragments.
#'
#' Returns a `list(system, board)` ready to feed into `omicsai_gen_text`.
#' Phase 2 only. Phase 3 (agentic / deep) will reuse the same fragments via
#' `omicsagentovi::Agent(..., system_prompt = bp$system)`; see
#' `.active_plans/edge_merge/round2_plan.md` §"Phase 3 outlook".
#'
#' @param slice the WGCNA result slot (`pgx$wgcna` or `pgx$wgcna_mox`). Shape
#'   dispatch (single / multi / consensus) is handled inside the builders.
#' @param pgx full pgx object (used for organism + annotation).
#' @param ai resolved `ai` list (already passed through `.ai_resolve_defaults`).
#' @return `list(system = <character>, board = <character>)`
#' @keywords internal
wgcna_assemble_prompt <- function(slice, pgx, ai) {
  if (!requireNamespace("omicsai", quietly = TRUE)) {
    stop("omicsai package required for AI report generation", call. = FALSE)
  }
  data_block <- wgcna_build_report_tables(slice, pgx)$text
  .ai_report_build_prompt(pgx, "wgcna", data_block)
}


# =============================================================================
# Section 4 — Entry point (called by pgx.update_reports orchestrator)
# =============================================================================

#' Generate a WGCNA AI report (phase-2 static path).
#'
#' Phase 2: assembles a single prompt via `wgcna_assemble_prompt` and issues one
#' `omicsai::omicsai_gen_text` call. Phase 3 (`ai$report_type = "deep"`) will
#' branch here onto an `omicsagentovi::Agent` path; until that lands the
#' orchestrator hard-errors on "deep" so this function only ever sees "normal".
#'
#' Slot dispatch (single-omics vs multi-omics) is delegated to
#' `wgcna_build_report_tables` (shape-based via `wgcna.getTopGenesAndSets`).
#' Callers route both `pgx$wgcna` and `pgx$wgcna_mox` to this single entry
#' point — see `.ai_dispatch_module` in `ai-report.R`.
#'
#' @param slice WGCNA result slot.
#' @param pgx full pgx object.
#' @param ai resolved `ai` list.
#' @return `list(report = <markdown>, prompt = <markdown>)`.
#' @export
wgcna.create_report <- function(slice, pgx, ai) {
  if (!requireNamespace("omicsai", quietly = TRUE)) {
    stop("omicsai package required for AI report generation", call. = FALSE)
  }
  bp  <- wgcna_assemble_prompt(slice, pgx, ai)
  .ai_report_run_prompt(bp, ai)
}
