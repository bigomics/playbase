##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##
##

#' @title WGCNA network construction and module detection
#'
#' @param pgx PGX object containing gene expression data
#' @param minmodsize Minimum module size cutoff
#' @param power Soft thresholding power for network construction
#' @param cutheight Cut height for module dendrogram
#' @param deepsplit Number of splits for module dendrogram
#' @param ngenes Number of genes to use (most variable)
#'
#' @return A list with the following components:
#' \describe{
#'   \item{datExpr}{Transposed expression matrix used for WGCNA.}
#'   \item{datTraits}{Numeric trait matrix derived from sample metadata.}
#'   \item{svTOM}{Singular vectors of the TOM matrix (low-rank approximation).}
#'   \item{net}{WGCNA network object from \code{blockwiseModules}.}
#'   \item{me.genes}{Gene-to-module membership vectors.}
#'   \item{me.colors}{Module color assignments.}
#'   \item{W}{Module eigengene matrix.}
#'   \item{modTraits}{Module-trait correlation results.}
#'   \item{stats}{Module-level statistics and test results.}
#'   \item{power}{Effective soft-thresholding power used.}
#'   \item{minModSize}{Effective minimum module size (may be capped for small datasets).}
#'   \item{mergeCutHeight}{Merge cut height for module merging.}
#'   \item{minKME}{Minimum eigengene connectivity to stay in a module.}
#'   \item{networktype}{Network type used (\code{"signed"}, \code{"unsigned"}, etc.).}
#'   \item{tomtype}{TOM type used.}
#'   \item{clust}{Clustering results from dimensionality reduction on TOM.}
#'   \item{annot}{Gene annotation table from the PGX object.}
#'   \item{experiment}{Experiment description from the PGX object.}
#' }
#'
#' @description Constructs a weighted gene coexpression network and detects
#' modules using WGCNA on a PGX object.
#'
#' @details This function takes a PGX object containing gene expression data.
#' It constructs a weighted gene coexpression network using the WGCNA package.
#' Soft thresholding power is set by \code{power}.
#'
#' Modules are detected by cutting the module dendrogram at \code{cutheight} and
#' with \code{deepsplit} splits. Only the \code{ngenes} most variable genes are used.
#'
#' The output is a list containing the WGCNA network object and module results,
#' including module assignments, colors, and summary statistics.
#'
#' @export
pgx.wgcna <- function(
  pgx,
  minmodsize = 20,
  power = 12,
  cutheight = 0.15,
  deepsplit = 2,
  minKME = 0.3,
  networktype = "signed",
  tomtype = "signed",
  numericlabels = FALSE,
  ngenes = 2000,
  maxBlockSize = 9999,
  gset.filter = "PATHWAY|HALLMARK|^GO|^C[1-9]",
  compute.enrichment = TRUE,
  summary = NULL,
  verbose = 1,
  progress = NULL
) {
  ## minmodsize=20;power=NULL;cutheight=0.15;deepsplit=2;ngenes=4000;networktype="signed";tomtype="signed";numericlabels=FALSE;ngenes=2000;gset.filter=NULL;minKME=0.8;maxBlockSize=5000

  if(!is.null(summary)) message("[pgx.wgcna] warning: summary is deprecated")
  
  samples <- pgx$samples
  contrasts <- pgx$contrasts
  ## no dot pheno
  samples <- samples[, grep("^[.]", colnames(samples), invert = TRUE), drop = FALSE]
  X <- pgx$X

  ## WGCNA does not work very well with scRNAseq due to sparsity.
  ## To bypass the issue, hdWGCNA computes metacells.
  ## Here we compute supercells too.
  is.singlecell <- !is.null(pgx$datatype) && pgx$datatype == "scRNAseq"
  is.toobig <- ncol(X) > 1000
  if (is.singlecell && is.toobig) {
    message("[pgx.wgcna] scRNAseq: >1K cells. Computing supercells.")
    counts <- pmax(2**X - 1, 0)
    ct <- samples[, "celltype"]
    group <- paste0(ct, ":", apply(pgx$contrasts, 1, paste, collapse = "_"))
    if ("batch" %in% colnames(samples)) {
      group <- paste0(group, ":", samples[, "batch"])
    }
    nb <- round(ncol(counts) / 500)
    message("[pgx.wgcna] running SuperCell. nb = ", nb)
    sc <- pgx.supercell(counts, samples, group = group, gamma = nb)
    message("[pgx.wgcna] SuperCell done: ", ncol(counts), " ->", ncol(sc$counts))
    message("[pgx.wgcna] Normalizing supercell matrix (logCPM)")
    X <- as.matrix(logCPM(sc$counts, total = 1e4, prior = 1))
    samples <- sc$meta
    remove(counts, ct, group, nb, sc)
    gc()
  }

  if (!is.null(pgx$datatype) && pgx$datatype == "multi-omics") {
    message("[pgx.wgcna] Performing multi-omics normalization")
    X <- normalizeMultiOmics(X)
  }

  if (!is.null(progress)) progress$set(message = "Calculating WGCNA...", value = 0.2)
  message("[pgx.wgcna] start wgcna.compute...")
  wgcna <- wgcna.compute(
    X = X,
    samples = samples,
    contrasts = contrasts,
    minmodsize = minmodsize, # default: min(20,...)
    power = power, # default: 12 (for signed)
    mergeCutHeight = cutheight, # default: 0.15
    deepsplit = deepsplit, # default: 2
    minKME = minKME, # default: 0.30
    networktype = networktype, # default: unsigned (but signed is better...)
    tomtype = tomtype, # default: signed
    numericlabels = numericlabels,
    maxBlockSize = maxBlockSize,
    ngenes = ngenes,
    verbose = verbose
  )
  message("[pgx.wgcna] finished computing wgcna.compute!")

  ## ---------------------------------------------------
  ## compute dimensionality reductions using TOM matrix
  ## ---------------------------------------------------
  wTOM <- NULL
  if ("TOM" %in% names(wgcna)) wTOM <- wgcna$TOM
  if (is.null(wTOM) && "svTOM" %in% names(wgcna)) {
    wTOM <- tcrossprod(wgcna$svTOM)
  }
  dissTOM <- 1 - wTOM
  rownames(dissTOM) <- colnames(dissTOM) <- colnames(wgcna$datExpr)
  clust <- pgx.clusterBigMatrix(dissTOM, methods = c("umap", "tsne", "pca"), dims = c(2))
  if ("cluster.genes" %in% names(pgx)) {
    posx <- pgx$cluster.genes$pos[["umap2d"]]
    clust[["umap2d"]] <- posx[colnames(wgcna$datExpr), ]
  }
  remove(dissTOM)

  ## ----------------------------------------------------
  ## Do geneset analysis
  ## ----------------------------------------------------
  if (compute.enrichment) {
    if (!is.null(progress)) progress$set(message = "Computing enrichment...", value = 0.4)
    message("computing module enrichment...")

    ## We augment pgx$GMT with original GSETxGENE to get better results
    ## for all modules.
    GMT0 <- getPlaydataGMT()
    GMT0 <- rename_by2(GMT0, pgx$genes, "symbol")
    GMT <- merge_sparse_matrix(pgx$GMT, GMT0)

    wgcna$gsea <- wgcna.computeModuleEnrichment(
      wgcna,
      annot = pgx$genes,
      GMT = GMT,
      methods = c("fisher", "gsetcor", "xcor"),
      ntop = 1000,
      xtop = 100,
      filter = gset.filter
    )
  }

  ## update
  wgcna <- wgcna.init(wgcna, annot=pgx$genes)
  
  ## save setting in object
  settings <- list(
    minmodsize = minmodsize,
    power = power,
    mergeCutHeight = cutheight,
    deepsplit = deepsplit,
    minKME = minKME,
    networktype = "signed",
    tomtype = "signed",
    # ngenes = 2000,
    # maxBlockSize = 9999,
    # gset.filter = "PATHWAY|HALLMARK|^GO|^C[1-9]",
    # compute.enrichment = TRUE,
    NULL
  )

  ## add to results object
  wgcna$clust <- clust
  wgcna$networktype <- networktype
  wgcna$tomtype <- tomtype
  wgcna$annot <- pgx$genes
  wgcna$experiment <- pgx$description
  wgcna$settings <- settings

  return(wgcna)
}


#' Combined GSETxGENE/MSETxMETABOLITE geneset matrix from playdata
#' @seealso WGCNAplus::getPlaydataGMT
getPlaydataGMT <- function() {
  WGCNAplus::getPlaydataGMT()
}

#' Core WGCNA computation (lower-level workhorse called by \code{\link{pgx.wgcna}})
#' @param X Gene expression matrix (genes x samples).
#' @param samples Sample metadata data.frame.
#' @param contrasts Contrast matrix (optional).
#' @param ngenes Number of most-variable genes to retain.
#' @param minmodsize Minimum module size (may be capped for small datasets).
#' @param power Soft-thresholding power.
#' @param mergeCutHeight Cut height for merging similar modules.
#' @param deepsplit Deep-split sensitivity for module detection.
#' @param minKME Minimum eigengene connectivity to retain a gene in a module.
#' @param treeCut Relative or absolute tree cut height.
#' @param treeCutCeiling Ceiling quantile used to rescale a relative \code{treeCut}.
#' @param networktype Network type (\code{"signed"}, \code{"unsigned"}, etc.).
#' @param tomtype TOM type (\code{"signed"}, \code{"unsigned"}, etc.).
#' @param clustMethod Hierarchical clustering method.
#' @param cutMethod Tree cut method (\code{"hybrid"}, \code{"tree"}, \code{"static"}).
#' @param calcMethod TOM calculation method (\code{"fast"}, \code{"adjacency"}, \code{"full"}).
#' @param lowrank Low-rank approximation dimension for fast TOM calculation.
#' @param numericlabels Use numeric module labels instead of colors.
#' @param maxBlockSize Maximum block size for module detection.
#' @param merge.dendro Merge per-block dendrograms into one.
#' @param compute.stats Compute gene-level statistics.
#' @param prefix Two-letter module label prefix.
#' @param sv.tom Number of singular vectors used for the reduced TOM.
#' @param drop.ref Drop reference level when expanding phenotype matrix.
#' @param net Precomputed network object (skips module detection if given).
#' @param is.multiomics Whether \code{X} is a multi-omics matrix (auto-detected if \code{NULL}).
#' @param verbose Verbosity level.
#' @return A list with \code{datExpr}, \code{datTraits}, \code{svTOM}, \code{net},
#'   \code{me.genes}, \code{me.colors}, \code{W}, \code{modTraits}, \code{stats},
#'   \code{power}, \code{minModSize}, \code{mergeCutHeight}, \code{minKME},
#'   \code{networktype}, \code{tomtype}.
#' @seealso WGCNAplus::computeWGCNA
#' @export
wgcna.compute <- function(X,
                          samples,
                          contrasts = NULL,
                          ngenes = 2000,
                          minmodsize = 20,
                          power = 12,
                          mergeCutHeight = 0.15,
                          deepsplit = 2,
                          minKME = 0.3,
                          treeCut = 0.99,
                          treeCutCeiling = 1,
                          networktype = "signed",
                          tomtype = "signed",
                          clustMethod = "average",
                          cutMethod = "hybrid",
                          calcMethod = "fast",
                          lowrank = 40,
                          numericlabels = FALSE,
                          maxBlockSize = 9999,
                          merge.dendro = TRUE,
                          compute.stats = TRUE,
                          prefix = "ME",
                          sv.tom = 40,
                          drop.ref = FALSE,
                          net = NULL,
                          is.multiomics = NULL,
                          verbose = 0) {
  WGCNAplus::computeWGCNA(
    X = X,
    samples = samples,
    contrasts = contrasts,
    ngenes = ngenes,
    minmodsize = minmodsize,
    power = power,
    mergeCutHeight = mergeCutHeight,
    deepsplit = deepsplit,
    minKME = minKME,
    treeCut = treeCut,
    treeCutCeiling = treeCutCeiling,
    networktype = networktype,
    tomtype = tomtype,
    clustMethod = clustMethod,
    cutMethod = cutMethod,
    calcMethod = calcMethod,
    lowrank = lowrank,
    numericlabels = numericlabels,
    maxBlockSize = maxBlockSize,
    merge.dendro = merge.dendro,
    compute.stats = compute.stats,
    prefix = prefix,
    sv.tom = sv.tom,
    drop.ref = drop.ref,
    net = net,
    is.multiomics = is.multiomics,
    verbose = verbose
  )
}

#' Create lasagna graph for wgcna object
#'
#'
wgcna.create_lasagna_model <- function(wgcna, layers = NULL) {
  if (is.null(layers)) {
    if (!is.null(wgcna$layers)) {
      layers <- wgcna$layers
    } else {
      layers <- list(gx = wgcna)
    }
  }

  ## Get eigengene matrices, remove grey modules
  ww <- lapply(layers, function(w) t(w$net$MEs))
  ww <- lapply(ww, function(w) w[!grepl("[A-Z]{2}grey$", rownames(w)), , drop = FALSE])
  ww <- ww[which(sapply(ww, nrow) > 0)]

  datTraits <- layers[[1]]$datTraits
  gdata <- list(X = ww, samples = datTraits)

  ## Create lasagna model
  lasagna <- lasagna.create_model(
    gdata,
    pheno = "expanded",
    ntop = 2000,
    nc = 20,
    add.sink = FALSE,
    intra = TRUE,
    fully_connect = FALSE,
    add.revpheno = TRUE
  )

  ## Multi-condition edge weighting
  graph <- lasagna.multisolve(
    lasagna,
    min_rho = 0.1,
    max_edges = 1000,
    # value.type = "logFC",
    fc.weight = TRUE,
    sp.weight = FALSE,
    prune = TRUE
  )

  lasagna$graph <- graph

  return(lasagna)
}

#' Init/update function
#'
#' @export
wgcna.init <- function(wgcna, annot = NULL,
                       sv.tom = 40, progress = NULL) {
  if(is.null(wgcna)) return(NULL)
  is.multi <- ("layers" %in% names(wgcna))

  if (is.null(wgcna$svTOM) && !is.null(wgcna$TOM)) {
    ## sv.tom <- ceiling(min(sv.tom,dim(datExpr)/2))
    message("[wgcna.init] computing reduced svTOM")
    if (!is.null(progress)) progress$inc(0.1, "computing reduced svTOM...")
    sv.tom <- min(sv.tom, ncol(TOM) - 1)
    sv <- irlba::irlba(wgcna$TOM, nv = sv.tom)
    svTOM <- sv$v %*% diag(sqrt(sv$d))
    rownames(svTOM) <- colnames(datExpr)
    wgcna$svTOM <- svTOM
    wgcna$TOM <- NULL
  }

  if (is.null(wgcna$modTraits) && !is.multi) {
    wgcna$modTraits <- cor(wgcna$net$MEs, wgcna$datTraits)
  }

  if (is.null(wgcna$stats)) {
    is.single <- all(c("net","datExpr","datTraits","svTOM") %in% names(wgcna)) 
    is.mox <- all(c("layers") %in% names(wgcna))
    if(is.single) {
      ## single-omics WGCNA
      message("[wgcna.init] computing missing gene statistics...")
      if (!is.null(progress)) progress$inc(0.1, "computing gene statistics...")
      wgcna$stats <- wgcna.computeGeneStats(
        wgcna$net, wgcna$datExpr, wgcna$datTraits, TOM = wgcna$svTOM)
    } else if(is.mox) {
      if (!is.null(progress)) progress$inc(0.1, "computing gene statistics...")    
      wgcna$stats <- wgcna.computeMultiGeneStats(wgcna)       
    } else {
      message("[wgcna.init] WARNING: could not calculate stats")
    }
  }

  if (is.null(wgcna$graph)) {
    message("[wgcna.init] computing graph...")
    if (!is.null(progress)) progress$inc(0.1, "computing graph...")
    wgcna$lasagna <- wgcna.create_lasagna_model(wgcna)
    wgcna$graph <- wgcna$lasagna$graph
  }

  wgcna
}


#' Core multi-omics WGCNA computation (wraps WGCNAplus::computeWGCNA_multiomics)
#'
#' WGCNAplus's version differs from this one in three ways that are
#' preserved here: (1) it does not auto-load the default playdata geneset
#' matrix when \code{GMT} is \code{NULL} (it silently disables enrichment
#' instead), (2) it does not truncate \code{names(dataX)} to two letters,
#' and (3) it returns a \code{report} field instead of an \code{annot}
#' field. \code{report} generation is disabled here to keep behavior and
#' return shape unchanged for existing callers.
#' @export
wgcna.compute_multiomics <- function(dataX,
                                     samples,
                                     contrasts = NULL,
                                     power = 12,
                                     ngenes = 2000,
                                     datanames = NULL,
                                     clustMethod = "average",
                                     cutMethod = "hybrid",
                                     minmodsize = 10,
                                     minKME = 0.3,
                                     deepsplit = 2,
                                     mergeCutHeight = 0.15,
                                     compute.enrichment = TRUE,
                                     xref = NULL,
                                     annot = NULL,
                                     GMT = NULL,
                                     drop.ref = FALSE,
                                     add.pheno = FALSE,
                                     add.gsets = FALSE,
                                     do.consensus = FALSE,
                                     gset.methods = c("fisher", "gsetcor", "xcor"),
                                     gset.ntop = 1000,
                                     gset.xtop = 100,
                                     experiment = "",
                                     verbose = 1,
                                     progress = NULL) {
  ## default-load the playdata GMT when needed (WGCNAplus does not do this)
  if (!is.null(annot) && !is.null(GMT)) {
    GMT <- rename_by2(GMT, annot, "symbol")
  }
  if (add.gsets || compute.enrichment) {
    GMT0 <- getPlaydataGMT()
    if (!is.null(annot)) GMT0 <- rename_by2(GMT0, annot, "symbol")
    if (!is.null(GMT)) {
      GMT <- merge_sparse_matrix(GMT, GMT0)
    } else {
      GMT <- GMT0
    }
  }

  ## WGCNAplus expects pre-truncated 2-letter datatype names
  names(dataX) <- substring(names(dataX), 1, 2)

  if (!is.null(progress)) {
    progress$set(message = paste("computing WGCNA modules..."), value = 0.33)
  }

  out <- WGCNAplus::computeWGCNA_multiomics(
    dataX = dataX,
    samples = samples,
    contrasts = contrasts,
    power = power,
    ngenes = ngenes,
    datanames = datanames,
    clustMethod = clustMethod,
    cutMethod = cutMethod,
    minmodsize = minmodsize,
    minKME = minKME,
    deepsplit = deepsplit,
    mergeCutHeight = mergeCutHeight,
    compute.enrichment = compute.enrichment,
    xref = xref,
    annot = annot,
    GMT = GMT,
    drop.ref = drop.ref,
    add.pheno = add.pheno,
    add.gsets = add.gsets,
    do.consensus = do.consensus,
    gset.methods = gset.methods,
    gset.ntop = gset.ntop,
    gset.xtop = gset.xtop,
    report = FALSE,
    experiment = experiment,
    verbose = verbose,
    progress = progress
  )

  if (!is.null(progress)) {
    progress$set(message = paste("computing module enrichment..."), value = 0.66)
  }

  ## restore `annot` field (WGCNAplus's version drops it in favor of `report`)
  out$report <- NULL
  out$annot <- annot

  out
}


#' Reimplementation for WGCNA::blockwiseModules(). This returns exact
#' same object as WGCNA::blockwiseModules() but is faster and allows
#' different clustering linkage methods (ward, complete).
#' @seealso WGCNAplus::computeModules
wgcna.computeModules <- function(
  datExpr,
  power = 6,
  networkType = "signed",
  TOM = NULL,
  TOMType = "signed",
  calcMethod = "fast",
  lowrank = 20,
  clustMethod = "average",
  cutMethod = "hybrid", ## hybrid, tree, static
  deepSplit = 2,
  treeCut = 0.99,
  treeCutCeiling = 1,
  minModuleSize = 20,
  minModuleSize2 = minModuleSize,
  mergeCutHeight = 0.15,
  minKMEtoStay = 0.3,
  numericLabels = FALSE, ## numeric or 'color' labels
  maxBlockSize = 9999,
  returnTOM = FALSE,
  verbose = 1
) {
  WGCNAplus::computeModules(
    datExpr = datExpr,
    power = power,
    networkType = networkType,
    TOM = TOM,
    TOMType = TOMType,
    calcMethod = calcMethod,
    lowrank = lowrank,
    clustMethod = clustMethod,
    cutMethod = cutMethod,
    deepSplit = deepSplit,
    treeCut = treeCut,
    treeCutCeiling = treeCutCeiling,
    minModuleSize = minModuleSize,
    minModuleSize2 = minModuleSize2,
    mergeCutHeight = mergeCutHeight,
    minKMEtoStay = minKMEtoStay,
    numericLabels = numericLabels,
    maxBlockSize = maxBlockSize,
    returnTOM = returnTOM,
    verbose = verbose
  )
}

#' Faster implementation of TOM computation using low-rank SVD
#' approximation.
#' @seealso WGCNAplus::fastTOMsimilarity
#' @export
fastTOMsimilarity <- function(A, tomtype = "signed", lowrank = 20) {
  WGCNAplus::fastTOMsimilarity(A, tomtype = tomtype, lowrank = lowrank)
}

#' Merge WGCNA modules whose eigengenes are closer than a cut height
#' @seealso WGCNAplus::mergeCloseModules
wgcna.mergeCloseModules <- function(datExpr, colors, cutHeight, MEs = NULL) {
  WGCNAplus::mergeCloseModules(datExpr = datExpr, colors = colors, cutHeight = cutHeight, MEs = MEs)
}

## ---------------------------------------------------------------------
## Gene statistics
## ---------------------------------------------------------------------

#' @seealso WGCNAplus::cortest
cortest <- function(X, Y) {
  WGCNAplus::cortest(X, Y)
}

#' Compute general feature statistics after WGCNA results.
#' @seealso WGCNAplus::computeGeneStats
wgcna.computeGeneStats <- function(net, datExpr=NULL, datTraits=NULL, TOM=NULL) {
  WGCNAplus::computeGeneStats(net = net, datExpr = datExpr, datTraits = datTraits, TOM = TOM)
}


wgcna.computeMultiGeneStats <- function(multi) {
  stats <- list()
  for(k in names(multi$layers)) {
    stats[[k]] <- wgcna.computeGeneStats(multi$layers[[k]]) 
  }

  for(k in names(stats)) {
    for(i in 3:length(stats[[k]])) {
      if(!is.null(rownames(stats[[k]][[i]]))) {
        rownames(stats[[k]][[i]]) <- paste0(k,":",rownames(stats[[k]][[i]]))
      }
      if(!is.null(names(stats[[k]][[i]]))) {
        names(stats[[k]][[i]]) <- paste0(k,":",names(stats[[k]][[i]]))
      }
    }    
  }

  xstats <- list()
  for(m in c("moduleTraitCor","moduleTraitPvalue","traitSignificance","TSPvalue",
    "foldChange","foldChangePvalue")) {
    xstats[[m]] <- do.call(rbind, lapply(stats, '[[', m))
  }
  for(m in c("moduleMembership","MMPvalue")) {
    ss <- list()
    for(k in names(stats)) {
      mm <- stats[[k]][[m]]
      all.mods <- rownames(xstats$moduleTraitCor)
      mods <- setdiff(all.mods,colnames(mm))
      mm.zero <- matrix(NA, nrow=nrow(mm), ncol=length(mods))
      colnames(mm.zero) <- mods
      rownames(mm.zero) <- rownames(mm)
      mm <- cbind(mm, mm.zero)
      mm <- mm[, all.mods]
      ss[[k]] <- mm
    }
    xstats[[m]] <- do.call(rbind, ss)    
  }  
  names(xstats)
  xstats[["geneCentrality"]] <- Reduce(c,sapply(stats, '[[', "geneCentrality"))
  xstats <- xstats[names(stats[[1]])]
  return(xstats)
}


#' Compute per-gene significance table for a module/trait
#' @seealso WGCNAplus::getGeneStats
#' @export
wgcna.getGeneStats <- function(wgcna, trait, module = NULL, plot = TRUE,
                               stats = NULL, labels = NULL, showlogFC = TRUE,
                               col = NULL, main = NULL) {
  WGCNAplus::getGeneStats(
    wgcna = wgcna, trait = trait, module = module, plot = plot,
    stats = stats, labels = labels, showlogFC = showlogFC,
    col = col, main = main
  )
}

## ----------------------------------------------------
## Perform geneset analysis on modules
## ----------------------------------------------------

#' Merge multiple ME matrices into one. Allow different dimensions.
#'
#' WGCNAplus::mergeME() has no explicit handling for \code{mlist = NULL,
#' me2 = <matrix>}; that edge case is preserved here as an early return.
#' @seealso WGCNAplus::mergeME
wgcna.mergeME <- function(mlist, me2 = NULL, prefix = FALSE) {
  if (is.null(mlist) && !is.null(me2)) {
    return(me2)
  }
  WGCNAplus::mergeME(mlist = mlist, me2 = me2, prefix = prefix)
}

#' Compute enrichment of each WGCNA module using various
#' methods. Handles single-type and multi-omics WGCNA objects.
#' @seealso WGCNAplus::computeModuleEnrichment
#' @export
wgcna.computeModuleEnrichment <- function(wgcna,
                                          GMT,
                                          multi = FALSE,
                                          methods = c("fisher", "gsetcor", "xcor"),
                                          ntop = 200,
                                          xtop = 100,
                                          annot = NULL,
                                          xref = NULL,
                                          min.genes = 3,
                                          min.rho = 0.8,
                                          filter = NULL) {
  ## Exempt dataset-specific sets from the reference-geneset name filter,
  ## which would otherwise leave wgcna$gsea empty. CUSTOM/TEST always: for
  ## custom organisms those are the only applicable sets. METABOLITE only as
  ## a fallback: METABOLITE_PATHWAY already matches the filter, so keeping
  ## METABOLITE_ONTOLOGY/CHEMCLASS unconditionally would drown real pathway
  ## hits. Mirrors the CUSTOM exemption in pgx.compute(). Testing colnames(GMT)
  ## matches WGCNAplus, which greps the same colnames after row-subsetting.
  if (!is.null(filter)) {
    if (!any(grepl(filter, colnames(GMT)))) filter <- paste0(filter, "|^METABOLITE")
    filter <- paste0(filter, "|^CUSTOM|^TEST:")
  }

  ## add.wgcna=FALSE keeps the WGCNAplus return shape identical to this
  ## function's historical contract (a plain gsea list, not the wgcna
  ## object with $gsea attached, which is WGCNAplus's default).
  WGCNAplus::computeModuleEnrichment(
    wgcna = wgcna, GMT = GMT, multi = multi, methods = methods,
    ntop = ntop, xtop = xtop, annot = annot, xref = xref,
    min.genes = min.genes, min.rho = min.rho, filter = filter,
    add.wgcna = FALSE
  )
}


#' @seealso WGCNAplus::run_enrichment_methods
wgcna.run_enrichment_methods <- function(ME, me.genes, GMT, geneX,
                                         methods = c("fisher", "gsetcor", "xcor"),
                                         ntop = 400, xtop = 100, min.genes = 3,
                                         min.rho = 0.8) {
  WGCNAplus::run_enrichment_methods(
    ME = ME, me.genes = me.genes, GMT = GMT, geneX = geneX,
    methods = methods, ntop = ntop, xtop = xtop, min.genes = min.genes,
    min.rho = min.rho
  )
}


#' Merge per-block dendrograms into a single dendrogram
#' @seealso WGCNAplus::merge_block_dendrograms
wgcna.merge_block_dendrograms <- function(net, X, method = 1) {
  WGCNAplus::merge_block_dendrograms(net = net, X = X, method = method)
}

#' Match module colors of a wgcna result to a reference set of colors
#' @seealso WGCNAplus::matchColors
#' @export
wgcna.matchColors <- function(wgcna, refcolors) {
  WGCNAplus::matchColors(wgcna = wgcna, refcolors = refcolors)
}

#' Find genes most correlated with each module eigengene across layers
#' @seealso WGCNAplus::getModuleCrossGenes
#' @export
wgcna.getModuleCrossGenes <- function(wgcna, ref = NULL, ngenes = 100,
                                      multi = TRUE, modules = NULL) {
  WGCNAplus::getModuleCrossGenes(wgcna = wgcna, ref = ref, ngenes = ngenes,
    multi = multi, modules = modules)
}


## =========================================================================
## CONSENSUS WGCNA
## =========================================================================

#' Run consensus WGCNA across a list of expression matrices
#'
#' WGCNAplus::runConsensusWGCNA() differs in ways preserved here: (1) it
#' defaults to \code{summary=TRUE}, which triggers an AI-based module
#' annotation call - disabled here to keep behavior unchanged; (2) it only
#' computes enrichment when \code{GMT} is explicitly supplied, whereas this
#' function auto-loads the default playdata GMT when enrichment is
#' requested; (3) it never calls the \code{progress} callback.
#' @seealso WGCNAplus::runConsensusWGCNA
#' @export
wgcna.runConsensusWGCNA <- function(exprList,
                                    phenoData,
                                    contrasts = NULL,
                                    GMT = NULL,
                                    annot = NULL,
                                    ngenes = 2000,
                                    power = 12,
                                    minModuleSize = 20,
                                    minKME = 0.3,
                                    mergeCutHeight = 0.15,
                                    deepSplit = 2,
                                    maxBlockSize = 9999,
                                    addCombined = FALSE,
                                    calcMethod = "fast",
                                    drop.ref = FALSE,
                                    cons.psig = 0.05,
                                    compute.stats = TRUE,
                                    compute.enrichment = TRUE,
                                    gsea.mingenes = 10,
                                    gsea.ntop = 1000,
                                    gset.methods = c("fisher", "gsetcor", "xcor"),
                                    verbose = 1,
                                    progress = NULL) {
  if (!is.null(progress)) progress$inc(0.1, "Computing layers...")

  if (compute.enrichment) {
    GMT0 <- getPlaydataGMT()
    GMT <- if (!is.null(GMT)) merge_sparse_matrix(GMT, GMT0) else GMT0
  }

  if (!is.null(progress)) progress$inc(0.2, "Computing consensus and enrichment...")

  WGCNAplus::runConsensusWGCNA(
    exprList = exprList,
    phenoData = phenoData,
    contrasts = contrasts,
    GMT = GMT,
    annot = annot,
    ngenes = ngenes,
    power = power,
    minModuleSize = minModuleSize,
    minKME = minKME,
    mergeCutHeight = mergeCutHeight,
    deepSplit = deepSplit,
    maxBlockSize = maxBlockSize,
    addCombined = addCombined,
    calcMethod = calcMethod,
    drop.ref = drop.ref,
    cons.psig = cons.psig,
    compute.stats = compute.stats,
    compute.enrichment = compute.enrichment,
    summary = FALSE,
    gsea.mingenes = gsea.mingenes,
    gsea.ntop = gsea.ntop,
    gset.methods = gset.methods,
    verbose = verbose,
    progress = progress
  )
}


#' @seealso WGCNAplus::createConsensusLayers
#' @export
wgcna.createConsensusLayers <- function(exprList,
                                        samples,
                                        contrasts = NULL,
                                        ngenes = 2000,
                                        power = 12,
                                        minModuleSize = 20,
                                        deepSplit = 2,
                                        mergeCutHeight = 0.15,
                                        minKME = 0.3,
                                        maxBlockSize = 9999,
                                        prefix = NULL,
                                        verbose = 1) {
  WGCNAplus::createConsensusLayers(
    exprList = exprList, samples = samples, contrasts = contrasts,
    ngenes = ngenes, power = power, minModuleSize = minModuleSize,
    deepSplit = deepSplit, mergeCutHeight = mergeCutHeight, minKME = minKME,
    maxBlockSize = maxBlockSize, prefix = prefix, verbose = verbose
  )
}

#' Compute gene statistics with original datExpr but with consensus
#' colors/labels for each layers. A separate function
#' wgcna.getConsensusGeneStats() extracts clean tables from this
#' results object.
#' @seealso WGCNAplus::computeConsensusGeneStats
#' @export
wgcna.computeConsensusGeneStats <- function(cons) {
  WGCNAplus::computeConsensusGeneStats(cons)
}

#' Extract clean consensus gene-statistics tables for a trait/module
#'
#' WGCNAplus::getConsensusGeneStats() has no guard for the case where
#' \code{trait} isn't present in a layer's stats (which makes
#' \code{wgcna.getGeneStats} return \code{NULL} for that layer); this
#' function checks for that upfront and propagates \code{NULL} instead of
#' erroring downstream.
#' @seealso WGCNAplus::getConsensusGeneStats
#' @export
wgcna.getConsensusGeneStats <- function(cons, stats, trait, module = NULL) {
  labels <- paste0("ME", cons$net$colors)
  gstats <- lapply(stats, function(s) {
    wgcna.getGeneStats(
      wgcna = NULL, stats = s, labels = labels, trait = trait,
      plot = FALSE, module = module, col = NULL, main = NULL
    )
  })
  if (any(sapply(gstats, is.null))) {
    return(NULL)
  }
  WGCNAplus::getConsensusGeneStats(cons = cons, stats = stats, trait = trait, module = module)
}


#' Compute consensus matrix from list of matrices. The consensus
#' matrix checks for consistent sign and minimal threshold for each
#' matrix. Optionally filters on consistent p-value.
#' @param ydim original dimension of data
#' @seealso WGCNAplus::computeConsensusMatrix
#' @export
wgcna.computeConsensusMatrix <- function(matlist, ydim, psig = 0.05, consfun = "min") {
  WGCNAplus::computeConsensusMatrix(matlist = matlist, ydim = ydim, psig = psig, consfun = consfun)
}

#' Compute consensus matrix from list of matrices. The consensus
#' matrix checks for consistent sign and minimal threshold for each
#' matrix. Optionally filters on consistent p-value.
#'
#' @export
wgcna.computeDistinctMatrix <- function(matlist, ydim, psig = 0.05, min.diff = 0.3,
                                        consmax = 0) {
  ## create difference module-trait matrix
  pv <- mapply(function(z, n) WGCNA::corPvalueStudent(z, n),
    matlist, ydim,
    SIMPLIFY = FALSE
  )
  matsign <- lapply(matlist, sign)
  Q <- matlist
  i <- 1
  for (i in 1:length(matlist)) {
    ## Any entry not significant is anyway invalid
    notsig <- (pv[[i]] > psig)
    Q[[i]][notsig] <- NA

    ## Any entry with too small difference with others is invalid
    refmat <- Reduce("+", matlist[-i]) / (length(matlist) - 1)
    diff <- matlist[[i]] - refmat
    notdiff <- (abs(diff) < min.diff)
    Q[[i]][notdiff] <- NA

    ## any entry that has consensus is invalid
    cons <- mapply(function(P, S) (P < 0.05) * (S == matsign[[i]]),
      pv[-i], matsign[-i],
      SIMPLIFY = FALSE
    )
    cons <- (Reduce("+", cons) > consmax) ## or function
    Q[[i]][cons] <- NA
  }
  return(Q)
}

#' Compute consensus enrichment by calculating overlapping enriched
#' terms.
#' @seealso WGCNAplus::computeConsensusModuleEnrichment
wgcna.computeConsensusModuleEnrichment <- function(cons,
                                                   GMT,
                                                   annot,
                                                   methods = c("fisher", "gsetcor", "xcor"),
                                                   min.genes = 3,
                                                   ntop = 400) {
  WGCNAplus::computeConsensusModuleEnrichment(
    cons = cons, GMT = GMT, annot = annot, methods = methods,
    min.genes = min.genes, ntop = ntop
  )
}

## =========================================================================
## PRESERVATION WGCNA
## =========================================================================

#' Run WGCNA module preservation analysis across a list of expression matrices
#'
#' WGCNAplus::runPreservationWGCNA() only computes enrichment when GMT is
#' explicitly supplied; this function auto-loads the default playdata GMT
#' (renamed via the annotation's "human_ortholog" column, matching the
#' original cross-species convention) when enrichment is requested.
#' @seealso WGCNAplus::runPreservationWGCNA
#' @export
wgcna.runPreservationWGCNA <- function(exprList,
                                       phenoData,
                                       contrasts = NULL,
                                       power = 12,
                                       reference = 1,
                                       add.merged = FALSE,
                                       ngenes = 2000,
                                       minModuleSize = 20,
                                       deepSplit = 2,
                                       annot = NULL,
                                       compute.stats = TRUE,
                                       compute.enrichment = TRUE,
                                       GMT = NULL,
                                       gset.methods = c("fisher", "gsetcor", "xcor")) {
  if (compute.enrichment) {
    GMT0 <- getPlaydataGMT()
    if (!is.null(annot)) GMT0 <- rename_by2(GMT0, annot, "human_ortholog")
    GMT <- if (!is.null(GMT)) merge_sparse_matrix(GMT, GMT0) else GMT0
  }

  WGCNAplus::runPreservationWGCNA(
    exprList = exprList,
    phenoData = phenoData,
    contrasts = contrasts,
    power = power,
    reference = reference,
    add.merged = add.merged,
    ngenes = ngenes,
    minModuleSize = minModuleSize,
    deepSplit = deepSplit,
    annot = annot,
    compute.stats = compute.stats,
    compute.enrichment = compute.enrichment,
    GMT = GMT,
    gset.methods = gset.methods
  )
}


#' Converts WGCNA labels (numeric or color) to colors.
#' @seealso WGCNAplus::labels2colors
wgcna.labels2colors <- function(colors, ...) {
  WGCNAplus::labels2colors(colors, ...)
}


#' Filter color vector by minimum KME and mergeCutHeight. Set color of
#' features with KME smaller than minKME to grey (or 0) group. Merge
#' similar modules with (module) correlation larger than
#' (1-mergeCutHeight) together.
#'
#' @export
wgcna.filterColors <- function(X, colors, minKME = 0.3, mergeCutHeight = 0.15,
                               minmodsize = 20, ntop = -1) {
  ## minKME=0.3;mergeCutHeight=0.15;minmodsize=20;ntop=-1

  sX <- X + 1e-8 * matrix(rnorm(length(X)), nrow(X), ncol(X))
  sX <- t(scale(t(sX)))

  ## get singular vectors and correct sign
  vv <- tapply(1:nrow(sX), colors, function(i) svd(sX[i, ], nv = 1)$v[, 1])
  mm <- tapply(1:nrow(sX), colors, function(i) colMeans(sX[i, ]))
  vv.sign <- mapply(function(a, b) sign(cor(a, b)), mm, vv)
  vv <- mapply(function(a, b) a * b, vv, vv.sign, SIMPLIFY = FALSE)

  kme <- rep(NA, nrow(X))
  names(kme) <- rownames(X)
  names(colors) <- rownames(X)

  grey.val <- NULL
  is.color <- mean(colors %in% WGCNA::standardColors(435)) > 0.8
  if (is.numeric(colors)) {
    colors <- as.integer(colors)
    grey.val <- 0
  } else {
    colors <- as.character(colors)
    grey.val <- "---"
    if (is.color) grey.val <- "grey"
  }
  names(colors) <- rownames(X)
  new.colors <- colors

  if (minKME > 0) {
    i <- 1
    for (i in 1:length(vv)) {
      ii <- which(colors == names(vv)[i])
      r <- cor(t(X[ii, ]), vv[[i]])[, 1]
      max(r)
      jj <- ii[which(r < minKME)]
      if (length(jj)) {
        new.colors[jj] <- NA
      }
      kme[ii] <- r
    }
    new.colors[is.na(new.colors)] <- grey.val
  }

  ## merge groups
  if (mergeCutHeight > 0) {
    mx <- rowmean(X, new.colors)
    rr <- cor(t(mx))
    diag(rr) <- 0
    merge.idx <- which(rr > (1 - mergeCutHeight), arr.ind = TRUE)
    if (nrow(merge.idx) > 0) {
      i <- 1
      for (i in 1:nrow(merge.idx)) {
        aa <- rownames(rr)[merge.idx[i, ]]
        jj <- which(new.colors %in% aa)
        max.color <- names(which.max(table(new.colors[jj])))
        new.colors[jj] <- max.color
      }
    }
  }

  ## remove small groups
  modsize <- table(new.colors)
  modsize
  if (min(modsize) < minmodsize) {
    small.mod <- names(which(modsize < minmodsize))
    sel <- which(new.colors %in% small.mod)
    new.colors[sel] <- NA
  }

  ## Filter by KME score
  if (ntop > 0) {
    keep <- tapply(names(kme), new.colors, function(i) head(names(sort(-kme[i])), ntop))
    keep <- unlist(keep)
    not.keep <- setdiff(names(kme), keep)
    if (length(not.keep)) new.colors[not.keep] <- NA
  }

  new.colors[which(is.na(new.colors))] <- grey.val
  ## if(!is.numeric(colors)) new.colors <- factor(new.colors)

  return(new.colors)
}

#' Wrapper to hclust from matrix using default WGCNA parameters.
#'
wgcna.tomclust <- function(X, power = 6) {
  A <- WGCNA::adjacency(t(X), power = power, type = "signed")
  TOM <- fastTOMsimilarity(A, tomtype = "signed", lowrank = 40)
  hc <- fastcluster::hclust(as.dist(1 - TOM), method = "average")
  hc
}


#' @seealso WGCNAplus::checkDendroHeights
wgcna.checkDendroHeights <- function(datExpr, n = 200, powers = NULL, maxpower = 20) {
  WGCNAplus::checkDendroHeights(datExpr = datExpr, n = n, powers = powers, maxpower = maxpower)
}


#' Better (?) method to pick soft threshold (aka power).
#' @seealso WGCNAplus::pickSoftThreshold
#' @export
wgcna.pickSoftThreshold <- function(datExpr, sft = NULL, rcut = 0.85,
                                    method = c("sft", "iqr")[1],
                                    nmax = -1, powers = NULL,
                                    verbose = 1) {
  WGCNAplus::pickSoftThreshold(datExpr = datExpr, sft = sft, rcut = rcut,
    method = method, nmax = nmax, powers = powers, verbose = verbose)
}

#' Scale a list of TOM matrices so that the quantiles (default p=0.95)
#' are equal after scaling with respect to the first TOM matrix.
#'
#'
wgcna.scaleTOMs <- function(TOMs, scaleP = 0.95) {
  nGenes <- nrow(TOMs[[1]])
  nSets <- length(TOMs)
  # Sample sufficiently large number of TOM entries
  nSamples <- as.integer(1 / (1 - scaleP) * 1000)
  # Choose the sampled TOM entries
  scaleSample <- sample(nGenes * (nGenes - 1) / 2, size = nSamples)
  TOMScalingSamples <- list()
  # These are TOM values at reference percentile
  scaleQuant <- rep(1, nSets)
  # Scaling powers to equalize reference TOM values
  scalePowers <- rep(1, nSets)
  # Loop over sets
  set <- 1
  for (set in 1:nSets)
  {
    # Select the sampled TOM entries
    tval <- as.dist(TOMs[[set]])[scaleSample]
    # Calculate the 95th percentile
    scaleQuant[set] <- quantile(tval, probs = scaleP, type = 8)
    TOMScalingSamples[[set]] <- tval

    # Scale the TOM
    if (set > 1) {
      scalePowers[set] <- log(scaleQuant[1]) / log(scaleQuant[set])
      TOMs[[set]] <- TOMs[[set]]^scalePowers[set]
    }
  }
  return(TOMs)
}


#' Calculate compound significance scores per gene
#' @param wgcna WGCNA result object with stats.
#' @return Data frame of compound significance scores.
#' @export
wgcna.calculateSignificanceScore <- function(wgcna, collapse = TRUE, sort.by = "score",
                                             digits = 4, annot = NULL, rownames = NULL,
                                             annot.cols = c("feature", "symbol", "gene_title")) {
  Q <- list()
  if (!is.null(wgcna$layers)) {
    ww <- wgcna$layers
  } else {
    ww <- list(gx = wgcna)
  }
  names(ww)
  k=1
  for (k in names(ww)) {
    stats <- ww[[k]]$stats
    if(is.null(stats)) {
      w <- ww[[k]]
      stats <- wgcna.computeGeneStats(w$net, w$datExpr, w$datTraits, TOM=w$svTOM)
    }
    m1 <- stats$moduleMembership    
    t1 <- stats$traitSignificance
    f1 <- stats$foldChange
    c1 <- ww[[k]]$net$labels[rownames(m1)]
    x1 <- stats$moduleMembership[cbind(1:nrow(m1), match(c1, colnames(m1)))]
    rxs <- function(x, k = 2) apply(x**k, 1, max, na.rm = TRUE)^(1 / k)
    # Q1 <- data.frame(c1, rxs(m1,k=1), rxs(t1), rxs(f1))
    Q1 <- data.frame(c1, x1, rxs(t1), rxs(f1))
    colnames(Q1) <- c("module", "MM", "max.TS", "max.FC")
    Q1$score <- apply(Q1[, c(2, 3, 4)], 1, prod)
    if (sort.by %in% colnames(Q1)) Q1 <- Q1[order(-Q1[, sort.by]), ]
    Q1[, 2:ncol(Q1)] <- round(Q1[, 2:ncol(Q1)], digits = digits)
    Q[[k]] <- Q1
  }

  if (is.null(annot)) annot <- wgcna$annot
  if (!is.null(annot.cols) && length(annot.cols) && !is.null(annot)) {
    i <- 1
    for (i in 1:length(Q)) {
      Q1 <- Q[[i]]
      Q1 <- rename_by2(Q1, annot, "feature", na.rm = FALSE)
      rr <- rownames(Q1)
      kk <- match(rr, rownames(annot))
      sel <- intersect(annot.cols, colnames(annot))
      aa <- annot[kk, sel, drop = FALSE]
      ## if feature and symbol are same drop
      if (all(c("feature", "symbol") %in% colnames(aa))) {
        if (mean(aa$symbol == aa$feature, na.rm = TRUE)) {
          aa$symbol <- NULL
        }
      }
      rr <- mofa.strip_prefix(rr)
      Q[[i]] <- data.frame(aa, Q[[i]], row.names = rr)
    }
  }

  if (length(Q) > 1) {
    for (k in 1:length(Q)) {
      rownames(Q[[k]]) <- paste0(
        names(Q)[k], ":",
        rownames(Q[[k]])
      )
    }
  }
  names(Q) <- NULL
  Q <- do.call(rbind, Q)
  if (sort.by %in% colnames(Q)) Q <- Q[order(-Q[, sort.by]), ]

  if (is.null(rownames)) rownames <- !("feature" %in% colnames(Q))
  if (!rownames) {
    rownames(Q) <- NULL
  }

  if (!collapse) {
    ## split by module
    Q <- tapply(1:nrow(Q), Q$module, function(i) Q[i, ])
  }

  return(Q)
}
