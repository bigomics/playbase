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


# -----------------------------------------------------------------------------
# Describe modules (legacy LLM helper)
# -----------------------------------------------------------------------------

#' Per-module LLM description helper.
#'
#' Used by `wgcna.create_report()` and other pre-AI-report-tab features. The
#' newer AI-report tab in board.wgcna does not call this; it builds its own
#' structured prompts via the omicsai stack.
#'
#' @export
wgcna.describeModules <- function(wgcna, ntop=50, psig = 0.05,
                                  annot=NULL, multi=FALSE, modules=NULL,
                                  experiment=NULL, verbose=1, model=DEFAULT_LLM,
                                  docstyle = "detailed summary", numpar = 2,
                                  level="gene")  {

  if(is.null(annot)) {
    message("[wgcna.describeModules] WARNING. user annot table is recommended.")
  }

  if(multi) {
    top <- wgcna.getMultiTopGenesAndSets(wgcna, annot=annot, ntop=ntop,
      psig=psig, level=NULL, rename="gene_title")
  } else {
    top <- wgcna.getTopGenesAndSets(wgcna, annot=annot, ntop=ntop,
      psig=psig, level=level, rename="gene_title")
  }

  if(is.null(modules)) {
    modules <- union(names(top$genes), names(top$sets))
  }

  if(is.null(experiment) && !is.null(wgcna$experiment)) experiment <- wgcna$experiment
  if(is.null(experiment)) experiment <- ""

  if(length(modules)==0) {
    info("[wgcna.describeModules] warning: empty module list!")
    return(NULL)
  }

  ## If no LLM is available we do just a manual summary
  model <- setdiff(model, c("", NA))
  if (is.null(model) || length(model) == 0) {
    desc <- list()
    for(m in modules) {
      ss=gg=pp=nn="<none>"

      if(!is.null(top$genes[[m]])) {
        gg <- paste( top$genes[[m]], collapse=', ')
      }
      if(!is.null(top$sets[[m]])) {
        ss <- paste( sub(".*:","",top$sets[[m]]), collapse='; ')
      }
      if(m %in% names(top$pheno)) {
        pp <- paste( top$pheno[[m]], collapse='; ')
      }
      if(m %in% names(top$neg.pheno)) {
        nn <- paste( top$neg.pheno[[m]], collapse='; ')
      }
      d <- ""
      if(!is.null(pp)) d <- paste(d, "**Positively correlated phenotypes**:", pp, "\n\n")
      if(!is.null(nn)) d <- paste(d, "**Negatively correlated phenotypes**:", nn, "\n\n")
      if(!is.null(gg) && gg!="") {
        d <- paste(d, "**Key genes**:", gg, "\n\n")
      }
      if(!is.null(ss) && ss!="") {
        d <- paste(d, "**Top enriched gene sets**:", ss, "\n\n")
      }
      desc[[m]] <- d
    }

    res <- list(
      prompt = NULL,
      questions = NULL,
      answers = desc
    )
    return(res)
  }

  prompt <- paste("Give a",docstyle,"of the main overall biological function of the following top enriched genesets belonging to module <MODULE>. After that, shortly discuss if any of these key genes/proteins/metabolites might be involved in the biological function. No need to mention all, just a few. Discuss the possible relationship with phenotypes <PHENOTYPES> of this experiment about \"<EXPERIMENT>\". Use maximum",numpar,"paragraphs. Use prose, do not use any bullet points or tables. \n\nHere is list of enriched gene sets:\n <GENESETS>\n\n")

  prompt <- paste("These are part of the results of a WGCNA analysis of an experiment about \"<EXPERIMENT>\". Give a", docstyle, "of the main overall biological function of the following top enriched genesets belonging to module <MODULE>. Discuss the possible relationship with positively correlated phenotypes <PHENOTYPES> and, if not obvious, negatively correlated phenotypes <NEGPHENOTYPES>. Use maximum", numpar, "paragraphs. Do not use any bullet points. \n\nHere is list of enriched gene sets: <GENESETS>\n")

  if (verbose > 1) cat(prompt)

  desc <- list()
  questions <- list()
  for (k in modules) {
    if (verbose > 0) message("Describing module ", k)

    ss=gg=pp=nn=""
    if(length(top$sets[[k]])>0) {
      ss <- sub( ".*:","", top$sets[[k]] )
      ss <- paste(ss, collapse=';')
    } else {
      ss <- "[no significant genesets]"
    }

    if(k %in% names(top$pheno)) {
      pp <- paste0("'",top$pheno[[k]],"'")
      pp <- paste( pp, collapse=';')
    }
    if(k %in% names(top$neg.pheno)) {
      nn <- paste0("'",top$neg.pheno[[k]],"'")
      nn <- paste( nn, collapse=';')
    }

    q <- prompt

    if(length(top$genes[[k]])>0) {
      gg <- paste( top$genes[[k]], collapse=';')
      q <- paste(q, "\nHere is the list of key genes/proteins/metabolites, or so-called 'features'. Only use features that are in this list in your answer. Do not mention features not in this list. : <KEYGENES>\n")
    }

    q <- sub("<MODULE>", k, q)
    q <- sub("<PHENOTYPES>", pp, q)
    q <- sub("<NEGPHENOTYPES>", nn, q)
    q <- sub("<EXPERIMENT>", experiment, q)
    q <- sub("<GENESETS>", ss, q)
    q <- sub("<KEYGENES>", gg, q)

    answer <- ""
    for (m in model) {
      if (verbose > 0) message("  ...asking LLM model ", m)
      a <- ai.ask(q, model = m)
      a <- paste0(a, "\n\n[AI generated using ", m, "]\n")
      if (length(model) > 1) a <- paste0("\n-------------------------------\n\n", a)
      answer <- paste0(answer, a)
    }

    desc[[k]] <- answer
    questions[[k]] <- q
  }

  res <- list(
    prompt = prompt,
    questions = questions,
    answers = desc
  )
  return(res)
}
