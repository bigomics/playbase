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
wgcna.mergeME <- function(mlist, me2 = NULL, prefix = FALSE) {
  if (is.null(mlist) && !is.null(me2)) {
    return(me2)
  }
  if (!is.null(me2) && !inherits(mlist, "list")) {
    mlist <- list(mlist, me2)
  }
  all.samples <- sapply(mlist, rownames, simplify = FALSE)
  all.samples <- unique(unlist(all.samples))
  if (prefix) {
    for (i in 1:length(mlist)) {
      colnames(mlist[[i]]) <- paste0(names(mlist)[i], ":", colnames(mlist[[i]]))
    }
  }
  is.mat <- all(sapply(mlist, inherits, what = "matrix"))
  all.me <- sapply(mlist, colnames, simplify = FALSE)
  all.me <- unique(unlist(all.me))
  M <- matrix(NA, nrow = length(all.samples), ncol = length(all.me))
  M <- as.data.frame(M)
  rownames(M) <- all.samples
  colnames(M) <- all.me
  i <- 1
  for (i in 1:length(mlist)) {
    ii <- match(rownames(mlist[[i]]), rownames(M))
    jj <- match(colnames(mlist[[i]]), colnames(M))
    M[ii, jj] <- mlist[[i]]
  }
  if (is.mat) M <- as.matrix(M)
  return(M)
}

#' Compute enrichment of each WGCNA module using various
#' methods. Handles single-type and multi-omics WGCNA objects.
#'
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
  if (!multi) {
    wgcna <- list(gx = wgcna)
    if (!is.null(annot)) rownames(annot) <- paste0("gx:", rownames(annot))
  }

  if (is.null(GMT)) {
    message("ERROR: must provide GMT")
    return(NULL)
  }

  ## create fake annotation table if no user annotation table is
  ## given.
  if (is.null(annot)) {
    gg <- lapply(wgcna, function(w) colnames(w$datExpr))
    ff <- list()
    for (i in 1:length(gg)) ff[[i]] <- paste0(names(wgcna)[i], ":", gg[[i]])
    gg <- unlist(gg)
    ff <- unlist(ff)
    annot <- data.frame(feature = ff, symbol = gg)
    rownames(annot) <- NULL
  }

  ## make sure GMT features  are in symbols
  symbol.col <- intersect(c("symbol", "gene_name"), colnames(annot))[1]
  GMT <- rename_by2(GMT, annot, symbol.col)

  ## add cross-referencing data??
  xref <- intersect(xref, names(wgcna))
  if (length(xref)) {
    message("[wgcna.computeModuleEnrichment] NOTE. Adding cross-correlation datatype: ", xref)
  }

  gsea <- list()
  dtype <- names(wgcna)[1]
  for (dtype in names(wgcna)) {
    ## collapse features to symbol
    sel <- unique(c(dtype, xref))
    datExpr <- lapply(wgcna[sel], function(w) t(as.matrix(w$datExpr)))
    geneX <- mofa.merge_data2(datExpr, merge.rows = "prefix", merge.cols = "union")
    geneX <- rename_by2(geneX, annot, symbol.col)

    ## check if overlap exists
    bg <- intersect(rownames(geneX), rownames(GMT))
    if (length(bg) == 0) {
      message("[wgcna.computeModuleEnrichment] WARNING. no overlapping genes for ", dtype)
      next()
    }

    G1 <- GMT[bg, , drop = FALSE]
    if (!is.null(filter)) {
      sel <- grep(filter, colnames(G1))
      if (length(sel)) G1 <- G1[, sel, drop = FALSE]
    }
    G1 <- G1[, which(Matrix::colSums(G1 != 0) >= min.genes), drop = FALSE]

    if (nrow(G1) >= 3 && ncol(G1) >= 3) {
      ## get eigengene members. convert to symbols.
      me.genes <- wgcna[[dtype]]$me.genes
      me.genes <- lapply(me.genes, function(g) probe2symbol(g, annot, query = symbol.col))

      ## Get eigengene matrix
      ME <- as.matrix(wgcna[[dtype]]$net$MEs)

      dt.gsea <- wgcna.run_enrichment_methods(
        ME = ME,
        me.genes = me.genes,
        GMT = G1,
        geneX = geneX,
        methods = methods,
        ntop = ntop,
        xtop = xtop,
        min.rho = min.rho
      )

      ## add to results
      gsea <- c(gsea, dt.gsea)
    }
  }

  return(gsea)
}


wgcna.run_enrichment_methods <- function(ME, me.genes, GMT, geneX,
                                         methods = c("fisher", "gsetcor", "xcor"),
                                         ntop = 400, xtop = 100, min.genes = 3,
                                         min.rho = 0.8) {
  rho.list <- list()
  pval.list <- list()

  ## align matrices
  bg <- intersect(rownames(GMT), rownames(geneX))
  GMT <- GMT[bg, ]
  geneX <- geneX[bg, ]

  ## select on minimum genes
  sel <- which(Matrix::colSums(GMT != 0) >= min.genes)
  GMT <- GMT[, sel]

  ## create gsetX
  gsetX <- plaid::plaid(geneX, matG = GMT)
  message("Computing enrichment for ", nrow(gsetX), " genesets")

  ## Add highly cross-correlated genes. limit xtop if geneX is too
  ## small.
  if (xtop > 0) {
    message("Adding high cross-correlated features. min.rho = ", min.rho)
    xtop <- min(xtop, round(nrow(geneX) / 4))
    nbx.genes <- list()
    for (k in colnames(ME)) {
      ss <- intersect(colnames(geneX), rownames(ME))
      if (length(ss) < 3) next()
      gx <- geneX[, ss, drop = FALSE]
      mx <- ME[ss, k, drop = FALSE]
      cx <- cor(t(gx), mx, use = "pairwise")[, 1]
      cx <- cx[!is.na(cx) & abs(cx) > min.rho]
      if (length(cx)) {
        nbx.genes[[k]] <- head(names(sort(-cx)), xtop)
      } else {
        nbx.genes[[k]] <- c()
      }
    }

    ## add to ME genes
    for (k in names(me.genes)) {
      me.genes[[k]] <- unique(c(me.genes[[k]], nbx.genes[[k]]))
    }
  }

  ## Here we correlate geneset score (averageCLR) directly with the
  ## module eigengene (ME). This should select genesets correlated
  ## with the ME.
  if ("gsetcor" %in% methods) {
    message("[wgcna.run_enrichment_methods] calculating single-sample geneset correlation...")
    rc.rho <- matrix(NA, ncol(GMT), ncol(ME))
    rc.pvalue <- matrix(NA, ncol(GMT), ncol(ME))
    dimnames(rc.rho) <- list(colnames(GMT), colnames(ME))
    dimnames(rc.pvalue) <- list(colnames(GMT), colnames(ME))
    jj <- which(rownames(gsetX) %in% colnames(GMT))
    kk <- intersect(colnames(gsetX), rownames(ME)) ## common samples
    tt <- cortest(t(gsetX[jj, kk]), ME[kk, ])
    rho.jj <- tt$rho
    pvalue.jj <- tt$pvalue
    ii <- match(rownames(gsetX)[jj], rownames(rc.rho))
    rc.rho[ii, ] <- rho.jj
    rc.pvalue[ii, ] <- pvalue.jj
    rho.list[["gsetcor"]] <- rc.rho
    pval.list[["gsetcor"]] <- rc.pvalue
  }

  ## Here we correlate the module eigengene (ME) with genes and then
  ## do a gset.rankcor() on the ME correlation.
  if ("xcor" %in% methods) {
    message("[wgcna.run_enrichment_methods] calculating eigengene correlation...")
    ss <- intersect(colnames(geneX), rownames(ME))
    rho <- cor(t(geneX[, ss, drop = FALSE]), ME[ss, , drop = FALSE], use = "pairwise")
    rho[is.na(rho)] <- 0
    rc <- gset.rankcor(rho, GMT, compute.p = TRUE) ## NEEDS CHECK!!!
    rho.list[["xcor"]] <- rc$rho
    pval.list[["xcor"]] <- rc$p.value
  }

  gmt <- mat2gmt(GMT)
  if (1) {
    ## we pre-select to make this faster
    Pmin <- sapply(pval.list, function(P) apply(P, 1, min))
    sel <- head(order(rowMeans(apply(Pmin, 2, rank))), 5 * ntop)
    message("[wgcna.run_enrichment_methods] pre-selecting ", length(sel), " sets for fgsea/Fisher test...")
    sel <- rownames(Pmin)[sel]
    gmt <- gmt[sel]
  }

  ## fGSEA. Compute gene correlation to eigengenes ME, then do
  ## pre-ranked enrichment on rho value.
  if ("fgsea" %in% methods) {
    message("[wgcna.run_enrichment_methods] calculating module fgsea...")
    ss <- intersect(colnames(geneX), rownames(ME))
    xrho <- cor(t(geneX[, ss, drop = FALSE]), ME[ss, , drop = FALSE], use = "pairwise")
    xrho[is.na(xrho)] <- 0
    res <- list()
    i <- 1
    for (i in 1:ncol(xrho)) {
      k <- colnames(xrho)[i]
      res[[k]] <- fgsea::fgsea(gmt, xrho[, i]) ## NEEDS CHECK!!!
    }
    pw <- res[[1]]$pathway
    res <- lapply(res, function(r) r[match(pw, r$pathway), ])
    nes <- sapply(res, function(r) r$NES)
    pval <- sapply(res, function(r) r$pval)
    rownames(nes) <- rownames(pval) <- pw
    colnames(nes) <- colnames(pval) <- names(res)
    rho.list[["fgsea"]] <- nes
    pval.list[["fgsea"]] <- pval
  }

  ## Perform fisher-test on (extended) ME genes. The ME genes might
  ## have been extended with most correlated genes.
  if ("fisher" %in% methods) {
    message("[wgcna.run_enrichment_methods] calculating Fisher tests...")
    rho <- matrix(NA, length(gmt), ncol(ME))
    pval <- matrix(NA, length(gmt), ncol(ME))
    dimnames(rho) <- list(names(gmt), colnames(ME))
    dimnames(pval) <- list(names(gmt), colnames(ME))

    ## perform Fisher test for all modules using the module genes
    i <- 1
    for (i in 1:ncol(rho)) {
      k <- colnames(rho)[i]
      gg <- me.genes[[k]]
      rr <- try(gset.fisher(gg, GMT,
        background = bg, fdr = 1,
        min.genes = -1, verbose = 0, sort.by = "none", no.pass = 1
      ))

      if (!"try-error" %in% class(rr)) {
        rr <- rr[match(rownames(rho), rownames(rr)), ]
        rho[, i] <- rr$odd.ratio
        pval[, i] <- rr$p.value
      }
    }

    ## handle infinite or NA
    rho[is.infinite(rho)] <- 2 * max(rho, na.rm = TRUE) ## Inf odd.ratio
    pval[is.na(pval)] <- 1
    rho[is.na(rho)] <- 0

    rho.list[["fisher"]] <- rho
    pval.list[["fisher"]] <- pval
  }

  lapply(rho.list, dim)

  ## ensure dimensions
  gsets <- Reduce(intersect, lapply(rho.list, rownames))
  modules <- Reduce(intersect, lapply(rho.list, colnames))
  rho.list <- lapply(rho.list, function(x) x[gsets, modules, drop = FALSE])
  pval.list <- lapply(pval.list, function(x) x[gsets, modules, drop = FALSE])

  ## Compute meta rank and pval. Handle NA for failing methods.
  pvalNA <- lapply(pval.list, function(x) {
    x[is.na(x)] <- 0
    x
  })
  ## pvalNA <- lapply(pval.list, function(x) {x[is.na(x)]=1;x})
  meta.p <- Reduce(pmax, pvalNA) ## NEED RETHINK!!!
  meta.q <- apply(meta.p, 2, p.adjust, method = "fdr")

  ## NEED RETHINK: how about negative FC???
  rnk.list <- lapply(rho.list, function(x) apply(x, 2, rank, na.last = "keep") / nrow(x))
  meta.rnk <- Reduce("+", rnk.list) / length(rnk.list)
  rnk.NAZERO <- lapply(rnk.list, function(x) {
    x[is.na(x)] <- 0
    x
  })
  rnk.NSUM <- Reduce("+", lapply(rnk.list, function(x) !is.na(x)))
  meta.rnk <- Reduce("+", rnk.NAZERO) / rnk.NSUM

  ## create dataframe by module
  message("[wgcna.run_enrichment_methods] creating dataframes...")
  gse.list <- list()
  i <- 1
  for (i in 1:ncol(meta.p)) {
    k <- colnames(meta.p)[i]
    pv <- sapply(pval.list, function(x) x[, i])
    colnames(pv) <- paste0("p.", colnames(pv))
    df <- data.frame(
      module = k,
      geneset = rownames(meta.p),
      score = meta.rnk[, i],
      p.value = meta.p[, i],
      q.value = meta.q[, i],
      pv
    )
    df <- df[order(-abs(df$score)), ]
    df <- head(df, ntop)
    gse.list[[k]] <- df
  }

  ## add genes
  ## gse.genes <- list()
  k <- names(gse.list)[1]
  for (k in names(gse.list)) {
    gset <- rownames(gse.list[[k]])
    gg <- me.genes[[k]]
    set.genes <- lapply(gmt[gset], function(s) intersect(s, gg))
    n0 <- sapply(gmt[gset], length)
    n1 <- sapply(set.genes, length)
    ## gse.genes[[k]] <- sort(table(unlist(set.genes)), decreasing = TRUE)
    set.genes <- sapply(set.genes, function(g) paste(sort(g), collapse = "|"))
    gse.list[[k]]$overlap <- paste0(n1, "/", n0)
    gse.list[[k]]$genes <- set.genes
  }

  return(gse.list)
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

#'
#'
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
  ## if(0) {
  ##   power=6;minKME=0.5;cutheight=0.15;deepSplit=2;maxBlockSize=5000;verbose=1;calcMethod="fast";addCombined=0;ngenes=2000;minModuleSize=20;mergeCutHeight=0.15
  ##   gsea.mingenes=20;gset.methods = c("fisher","gsetcor","xcor")
  ## }

  colors <- NULL

  ## Align and reduce matrices if needed
  gg <- Reduce(intersect, lapply(exprList, rownames))
  exprList <- lapply(exprList, function(x) x[gg, , drop = FALSE])
  if (length(gg) > ngenes) {
    sdx <- Reduce("*", lapply(exprList, function(x) matrixStats::rowSds(x)))
    ii <- head(order(-sdx), ngenes)
    exprList <- lapply(exprList, function(x) x[ii, , drop = FALSE])
  }

  if (addCombined) {
    exprList[["Combined"]] <- do.call(cbind, exprList)
  }

  exprsamples <- unlist(lapply(exprList, colnames))
  if (!all(exprsamples %in% rownames(phenoData))) {
    stop("samples mismatch for exprList and phenoData")
  }

  multiExpr <- WGCNA::list2multiData(lapply(exprList, Matrix::t))
  cor <- WGCNA::cor ## needed...

  if (!is.null(power) && length(power) == 1) {
    power <- rep(power, length(multiExpr))
  }

  # module detection procedure
  layers <- list()
  if (!is.null(progress)) progress$inc(0.1, "Computing layers...")
  for (i in 1:length(multiExpr)) {
    k <- names(multiExpr)[i]
    message("[wgcna.runConsensusWGCNA] >>> computing WGCNA for ", k)
    X <- Matrix::t(multiExpr[[i]]$data)
    layers[[k]] <- wgcna.compute(
      X = X,
      samples = phenoData,
      contrasts = contrasts,
      ngenes = ngenes,
      power = power[i],
      minmodsize = minModuleSize,
      calcMethod = calcMethod,
      deepsplit = deepSplit,
      mergeCutHeight = mergeCutHeight,
      numericlabels = FALSE,
      minKME = minKME,
      maxBlockSize = maxBlockSize,
      compute.stats = compute.stats,
      sv.tom = 40,
      verbose = verbose
    )
  }

  # now we run automatic consensus module detection
  message("[wgcna.runConsensusWGCNA] >>> computing CONSENSUS modules...")
  if (!is.null(progress)) progress$inc(0.1, "Computing consensus...")
  consensusPower <- unlist(sapply(layers, function(w) w$net$power))
  if (is.null(consensusPower) && !is.null(power)) {
    consensusPower <- power
  }
  if (is.null(consensusPower)) {
    consensusPower <- rep(12, length(layers))
  }

  sel <- setdiff(names(multiExpr), c("Combined"))
  cons <- WGCNA::blockwiseConsensusModules(
    multiExpr[sel],
    power = as.numeric(consensusPower),
    networkType = "signed",
    TOMType = "signed",
    minModuleSize = as.integer(minModuleSize),
    deepSplit = as.integer(deepSplit),
    mergeCutHeight = as.numeric(mergeCutHeight),
    numericLabels = FALSE,
    minKMEtoStay = as.numeric(minKME),
    maxBlockSize = as.integer(maxBlockSize),
    saveTOMs = FALSE,
    useDiskCache = FALSE,
    verbose = verbose
  )
  cons$power <- consensusPower

  ## create and match colors
  for (i in 1:length(layers)) {
    layers[[i]] <- wgcna.matchColors(layers[[i]], cons$colors)
  }

  layers.colors <- sapply(layers, function(r) r$net$colors)
  colors <- cbind(Consensus = cons$colors, layers.colors)

  ## add labels to dendrogram
  for (i in 1:length(cons$dendrograms)) {
    ii <- which(cons$goodGenes & cons$blocks == i)
    xnames <- names(cons$colors)
    cons$dendrograms[[i]]$labels <- xnames[ii]
  }

  ## merge dendrograms ????
  message("[wgcna.runConsensusWGCNA] merge_block_dendrograms...")
  multiX <- Matrix::t(do.call(rbind, lapply(exprList, function(x) scale(t(x)))))
  merged <- try(wgcna.merge_block_dendrograms(cons, multiX))
  if (!inherits(merged, "try-error")) {
    cons$merged_dendro <- merged
  } else {
    cons$merged_dendro <- NULL
  }

  ## create module-trait matrices for each set
  message("[wgcna.runConsensusWGCNA] >>> computing module-traits matrices...")
  datTraits <- 1 * expandPhenoMatrix(
    phenoData,
    drop.ref = drop.ref,
    keep.numeric = TRUE
  )
  if (!is.null(contrasts)) {
    message("[wgcna.runConsensusWGCNA] adding contrasts to datTraits")
    ctx <- makeContrastsFromLabelMatrix(contrasts)
    ctx <- sign(ctx)
    ctx[ctx == 0] <- NA
    ctx[ctx == -1] <- 0
    datTraits <- cbind(datTraits, ctx)
  }

  zlist <- list()
  k <- 1
  for (k in names(cons$multiME)) {
    M <- (cons$multiME[[k]][[1]])
    Z <- datTraits
    kk <- intersect(rownames(M), rownames(Z))
    zrho <- cor(M[kk, ], Z[kk, ], use = "pairwise")
    zrho[is.na(zrho)] <- 0 ## NEED RETHINK!!
    zlist[[k]] <- zrho
  }

  ## create consensus module-trait matrix
  ydim <- sapply(exprList, ncol)
  consZ <- wgcna.computeConsensusMatrix(zlist, ydim = ydim, psig = cons.psig)
  avgZ <- Reduce("+", zlist) / length(zlist)

  ## add slots
  datExpr <- lapply(exprList, Matrix::t)

  res <- list(
    net = cons,
    layers = layers,
    datExpr = datExpr,
    datTraits = datTraits,
    modTraits = avgZ,
    consModTraits = consZ,
    dendro = cons$merged_dendro,
    colors = colors,
    zlist = zlist,
    ydim = ydim,
    class = "consensus"
  )

  ## run stats
  if (compute.stats) {
    message("[wgcna.runConsensusWGCNA] >>> computing gene statistics...")
    res$stats <- wgcna.computeConsensusGeneStats(res)
  }

  ## run enrichment
  if (compute.enrichment) {
    if (!is.null(progress)) progress$inc(0.2, "Computing enrichment...")
    message("[wgcna.runConsensusWGCNA] >>> computing module enrichment...")
    if (!is.null(GMT)) {
      GMT0 <- getPlaydataGMT()
      if (!is.null(annot)) GMT0 <- rename_by2(GMT0, annot, "symbol")
      GMT <- merge_sparse_matrix(GMT, GMT0)
      remove(GMT0)
    } else {
      GMT <- getPlaydataGMT()
      if (!is.null(annot)) GMT <- rename_by2(GMT, annot, "symbol")
    }
    res$gsea <- wgcna.computeConsensusModuleEnrichment(
      res,
      GMT = GMT,
      method = gset.methods,
      annot = annot,
      min.genes = gsea.mingenes,
      ntop = gsea.ntop
    )
  }

  res
}


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
  if (0) {
    ngenes <- 2000
    power <- 12
    minModuleSize <- 5
    deepSplit <- 2
    mergeCutHeight <- 0.15
    minKME <- 0.3
    maxBlockSize <- 9999
    verbose <- 1
    prefix <- NULL
  }

  if (is.null(prefix)) prefix <- names(exprList)
  nx <- length(exprList)
  prefix <- head(rep(prefix, nx), nx)

  ## reduce
  message("[wgcna.computeConsensusLayers] Aligning matrices...")
  gg <- Reduce(intersect, lapply(exprList, rownames))
  exprList <- lapply(exprList, function(x) x[gg, ])

  if (length(gg) > ngenes) {
    message("[wgcna.computeConsensusLayers] Reducing to ", ngenes, " genes")
    sdx <- Reduce("*", lapply(exprList, function(x) matrixStats::rowSds(x)))
    ii <- head(order(-sdx), ngenes)
    exprList <- lapply(exprList, function(x) x[ii, ])
  }
  multiExpr <- WGCNA::list2multiData(lapply(exprList, Matrix::t))

  ## determine power vector
  if (is.null(power) || any(is.na(power))) power <- "sft"
  if (as.character(power[1]) %in% c("sft", "iqr")) {
    ## Estimate best power
    power <- power[1]
    message("[wgcna.createConsensusLayers] optimal power method = ", power)
    est.power <- rep(NA, length(exprList))
    i <- 1
    for (i in 1:length(exprList)) {
      p <- wgcna.pickSoftThreshold(
        Matrix::t(exprList[[i]]),
        sft = NULL, rcut = 0.85, powers = NULL,
        method = power, nmax = 1000, verbose = 0
      )
      if (length(p) == 0 || is.null(p)) p <- NA
      est.power[i] <- p
    }
    est.power
    power <- ifelse(is.na(est.power), 12, est.power)
  } else {
    power <- as.numeric(power)
  }
  nw <- length(exprList)
  power <- head(rep(power, nw), nw)
  names(power) <- names(exprList)

  message("[wgcna.computeConsensusLayers] Computing consensus modules...")
  cons <- WGCNA::blockwiseConsensusModules(
    multiExpr,
    power = as.numeric(power),
    networkType = "signed",
    TOMType = "signed",
    minModuleSize = as.integer(minModuleSize),
    deepSplit = as.integer(deepSplit),
    mergeCutHeight = as.numeric(mergeCutHeight),
    numericLabels = FALSE,
    minKMEtoStay = as.numeric(minKME),
    maxBlockSize = as.integer(maxBlockSize),
    saveTOMs = FALSE,
    useDiskCache = FALSE,
    verbose = verbose
  )

  ##
  message("[wgcna.computeConsensusLayers] Creating consensus layers...")
  aligned <- list()
  i <- 1
  for (i in 1:length(exprList)) {
    k <- names(exprList)[i]
    sel <- c(
      "colors", "unmergedColors", "goodSamples", "goodGenes",
      "dendrograms", "blockGenes", "blocks"
    )
    net <- cons[sel]
    net$power <- power[i]
    X <- exprList[[i]]
    w <- wgcna.compute(
      X = exprList[[i]],
      samples = samples,
      contrasts = contrasts,
      prefix = prefix[i],
      ngenes = -1,
      net = net,
      calcMethod = "fast",
      sv.tom = 0
    )
    aligned[[k]] <- w
  }

  return(aligned)
}

#' Compute gene statistics with original datExpr but with consensus
#' colors/labels for each layers. A separate function
#' wgcna.getConsensusGeneStats() extracts clean tables from this
#' results object.
#'
#' @export
wgcna.computeConsensusGeneStats <- function(cons) {
  k <- names(cons$layers)[1]
  stats <- list()
  for (k in names(cons$layers)) {
    w <- cons$layers[[k]]
    colors <- cons$net$colors
    wMEs <- cons$net$multiMEs[[k]]$data
    wnet <- list(MEs = wMEs, colors = colors)
    stats[[k]] <- wgcna.computeGeneStats(
      wnet, w$datExpr, w$datTraits,
      TOM = NULL
    )
  }
  return(stats)
}

#'
#'
#' @export
wgcna.getConsensusGeneStats <- function(cons, stats, trait, module = NULL) {
  ## create extended color vector
  labels <- paste0("ME", cons$net$colors)
  gstats <- list()
  for (k in names(stats)) {
    gstats[[k]] <- wgcna.getGeneStats(
      wgcna = NULL,
      stats = stats[[k]],
      labels = labels,
      trait = trait,
      plot = FALSE,
      module = module,
      col = NULL,
      main = NULL
    )
  }

  ## If any layer returned NULL (trait not in stats), propagate NULL
  if (any(sapply(gstats, is.null))) {
    return(NULL)
  }

  ## Align rows
  ff <- gstats[[1]]$feature
  for (k in names(gstats)) {
    ii <- match(ff, gstats[[k]]$feature)
    gstats[[k]] <- gstats[[k]][ii, ]
  }

  ## Compute consensus statistics. Consensus statistics are computed
  ## as geometric mean of score variables, and/or maximum pvalue for
  ## p.value columns.
  xcols <- c(3, 4, 6, 8)
  pcols <- c(10, 5, 7, 9)
  pcols1 <- c(5, 7, 9)
  xcols <- c("score", "moduleMembership", "traitSignificance", "foldChange")
  pcols <- c("scorePvalue", "MMPvalue", "TSPvalue", "foldChangePvalue")
  pcols1 <- pcols[-1]
  for (i in 1:length(gstats)) {
    gstats[[i]][, "scorePvalue"] <- apply(gstats[[i]][, pcols1], 1, max, na.rm = TRUE)
  }
  ## xc <- lapply(gstats, function(x) log(abs(x[,xcols])*(x[,pcols]<0.05)))
  xc <- lapply(gstats, function(x) log(abs(x[, xcols])))
  xc <- exp(Reduce("+", xc) / length(xc))
  xp <- Reduce(pmax, lapply(gstats, function(x) x[, pcols]))
  df3 <- data.frame(gstats[[1]][, 1:2], xc, xp)
  df3 <- df3[, colnames(gstats[[1]])]
  head(df3)

  ## Determine consensus status. Feature is 'C' (concordant) if sign
  ## in all layers are equal and significant. 'D' (discordant) if sign
  ## if not equal in all layers but significant. 'N' is any is
  ## non-significant.
  sign.pos <- Reduce("*", lapply(gstats, function(g) sign(g$score) == 1))
  sign.neg <- Reduce("*", lapply(gstats, function(g) sign(g$score) == -1))
  allsig <- Reduce("*", lapply(gstats, function(g) (g$scorePvalue) < 0.05))
  table(allsig)
  consensus <- c("D", "C")[1 + 1 * (sign.pos | sign.neg)]
  consensus[which(allsig == 0)] <- "N"
  cons.df <- data.frame(df3[, 1:2], consensus, df3[, -c(1, 2)])
  head(cons.df)

  ## This creates the full stats matrix (all subgroups)
  df1 <- gstats[[1]][, c("feature", "module")]
  df2 <- gstats[[1]][, 0]
  cols <- colnames(gstats[[1]])[-c(1:2)]
  for (k in cols) {
    xx <- sapply(gstats, function(g) g[, k])
    df2[[k]] <- I(xx)
  }
  df2 <- do.call(cbind, lapply(df2, unclass))
  newcols <- unlist(lapply(cols, function(k) paste0(k, ".", names(gstats))))
  colnames(df2) <- newcols
  full.df <- data.frame(df1, consensus = cons.df$consensus, df2)

  ## sort??
  ii <- order(-cons.df$score * sign(mean(cons.df$score, na.rm = TRUE)))
  cons.df <- cons.df[ii, ]
  full.df <- full.df[ii, ]

  list(
    consensus = cons.df,
    full = full.df
  )
}


#' Compute consensus matrix from list of matrices. The consensus
#' matrix checks for consistent sign and minimal threshold for each
#' matrix. Optionally filters on consistent p-value.
#'
#' @param ydim original dimension of data
#'
#'
#' @export
wgcna.computeConsensusMatrix <- function(matlist, ydim, psig = 0.05, consfun = "min") {
  if (length(ydim) == 1) ydim <- rep(ydim[1], length(matlist))
  pv <- mapply(function(z, n) {
    WGCNA::corPvalueStudent(z, n)
  }, matlist, ydim, SIMPLIFY = FALSE)
  for (i in 1:length(pv)) pv[[i]][is.na(pv[[i]])] <- 1 ## missing???

  ## create consensus module-trait matrix
  matsign <- list()
  for (i in 1:length(matlist)) {
    matsign[[i]] <- sign(matlist[[i]]) * (pv[[i]] <= psig)
  }
  matsign <- lapply(matsign, function(x) {
    x[is.na(x)] <- 0
    x
  })
  all.pos <- Reduce("*", lapply(matsign, function(z) (z >= 0)))
  all.neg <- Reduce("*", lapply(matsign, function(z) (z <= 0)))
  concordant <- (all.pos | all.neg)

  matlistN <- Reduce("+", lapply(matlist, function(x) !is.na(x)))
  matlist0 <- lapply(matlist, function(x) {
    x[is.na(x)] <- 0
    x
  })

  zsign <- sign(Reduce("+", matsign)) ## mean sign??
  if (consfun == "min") {
    pminFUN <- function(...) pmin(..., na.rm = TRUE)
    consZ <- Reduce(pminFUN, lapply(matlist, abs)) * zsign
  } else if (consfun == "gmean") {
    ## geometric mean
    matlistG <- lapply(matlist, function(x) {
      x <- log(abs(x))
      x[is.na(x)] <- 0
      x
    })
    consZ <- exp(Reduce("+", matlistG) / matlistN)
    consZ <- consZ * zsign
  } else {
    ## mean
    consZ <- Reduce("+", matlist0) / matlistN
  }
  consZ[!concordant] <- NA

  if (psig < 1) {
    ## enforce strong consensus. All layers must be strictly
    ## significant.
    all.sig <- Reduce("*", lapply(pv, function(p) 1 * (p <= psig)))
    consZ[!all.sig] <- NA
  }
  return(consZ)
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
#'
wgcna.computeConsensusModuleEnrichment <- function(cons,
                                                   GMT,
                                                   annot,
                                                   methods = c("fisher", "gsetcor", "xcor"),
                                                   min.genes = 3,
                                                   ntop = 400) {
  if (0) {
    methods <- c("fisher", "gsetcor", "xcor")
    min.genes <- 3
    ntop <- 400
    annot <- NULL
    GMT <- Matrix::t(playdata::GSETxGENE)
  }

  if (is.null(GMT)) {
    message("ERROR: must provide GMT")
    return(NULL)
  }

  gseaX <- list()
  i <- 1
  for (i in 1:length(cons$datExpr)) {
    geneX <- t(cons$datExpr[[i]])
    dim(geneX)

    ## Rename everything to symbols
    if (!is.null(annot)) {
      geneX <- rename_by2(geneX, annot, "symbol")
      GMT <- rename_by2(GMT, annot, "symbol")
    }
    ng <- length(intersect(rownames(geneX), rownames(GMT)))
    if (ng == 0) {
      message("[wgcna.computeConsensusModuleEnrichment] ERROR. No symbol overlap.")
      return(NULL)
    }
    symbols <- intersect(rownames(GMT), rownames(geneX))
    message("[wgcna.computeConsensusModuleEnrichment] number of symbols: ", length(symbols))
    geneX <- geneX[symbols, ]
    GMT <- GMT[symbols, ]

    ## select on minimum gene sets size
    sel <- which(Matrix::colSums(GMT != 0) >= min.genes)
    GMT <- GMT[, sel]

    ## Create extended Eigengene matrix (ME). ME should be nicely
    ## normalized/scaled so we just rbind across datasets
    ME <- cons$net$multiMEs[[i]]$data
    dim(ME)

    ## get genes in modules
    me.genes <- tapply(names(cons$net$colors), cons$net$colors, list)
    names(me.genes) <- paste0("ME", names(me.genes))
    if (!is.null(annot)) {
      me.genes <- lapply(me.genes, function(gg) probe2symbol(gg, annot))
    }
    me.genes <- lapply(me.genes, function(g) intersect(g, symbols))
    rownames(ME)
    colnames(geneX) <- rownames(ME)

    k <- names(cons$datExpr)[i]
    gseaX[[k]] <- wgcna.run_enrichment_methods(
      ME,
      me.genes = me.genes,
      GMT = GMT,
      geneX = geneX,
      methods = methods,
      min.genes = min.genes,
      ntop = ntop
    )
  }

  cons.gsea <- list()
  m <- 1
  for (m in names(gseaX[[1]])) {
    xx <- lapply(gseaX, function(g) g[[m]])
    sel <- Reduce(intersect, lapply(xx, rownames))
    if (length(sel) > 0) {
      if (length(sel) == 1) sel <- c(sel, sel) ## length==1 crashes...
      xx <- lapply(xx, function(x) x[sel, , drop = FALSE])
      xx.score <- sapply(xx, function(x) x[, "score"])
      colnames(xx.score) <- paste0("score.", colnames(xx.score))

      xx.pvalue <- lapply(xx, function(x) x[, grep("^p", colnames(x))])
      xx.pvalue <- do.call(cbind, xx.pvalue)

      m.score <- rowMeans(xx.score, na.rm = TRUE)
      m.pvalue <- apply(sapply(xx, function(x) x[, "p.value"]), 1, max, na.rm = TRUE)
      m.qvalue <- p.adjust(m.pvalue)
      df <- data.frame(
        module = xx[[1]]$module,
        geneset = xx[[1]]$geneset,
        score = m.score,
        xx.score,
        p.value = m.pvalue,
        q.value = m.qvalue,
        overlap = xx[[1]]$overlap,
        genes = xx[[1]]$genes,
        xx.pvalue
      )
      df <- df[order(df$p.value), ]
      # df <- df[!duplicated(df$geneset),,drop=FALSE]
      cons.gsea[[m]] <- df
    }
  }

  return(cons.gsea)
}

## =========================================================================
## PRESERVATION WGCNA
## =========================================================================

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
  if (is.character(reference)) {
    reference <- match(reference, names(exprList))
  }
  if (reference > 0) {
    reference.name <- names(exprList)[reference]
  } else {
    reference.name <- "Consensus"
  }

  ## multiset WGCNA
  pres <- wgcna.runConsensusWGCNA(
    exprList,
    phenoData = phenoData,
    contrasts = contrasts,
    GMT = NULL, ## no enrichment now
    annot = NULL, ## no enrichment now
    ngenes = ngenes,
    power = power,
    minModuleSize = minModuleSize,
    minKME = 0.3,
    mergeCutHeight = 0.15,
    deepSplit = deepSplit,
    maxBlockSize = 9999,
    addCombined = FALSE,
    calcMethod = "fast",
    drop.ref = FALSE,
    compute.stats = FALSE,
    compute.enrichment = FALSE,
    gsea.mingenes = 10,
    gset.methods = gset.methods
  )

  colorList <- lapply(pres$layers, function(w) w$net$colors)
  names(colorList) <- names(pres$layers)
  exprList <- lapply(pres$layers, function(w) t(w$datExpr))

  if (add.merged || reference == 0) {
    message("[wgcna.runPreservationWGCNA] adding merged layer...")
    cX <- lapply(exprList, function(x) x - rowMeans(x))
    merged <- do.call(cbind, cX)
    exprList$Merged <- NULL
    exprList <- c(list(Merged = merged), exprList)
    cons.colors <- pres$net$colors
    colorList <- c(list(Consensus = cons.colors), colorList)
    reference <- reference + 1
  }

  message("[wgcna.runPreservationWGCNA] running WGCNA::modulePreservation...")
  multiExpr <- WGCNA::list2multiData(lapply(exprList, Matrix::t))
  mp <- WGCNA::modulePreservation(
    multiExpr,
    colorList,
    referenceNetworks = reference,
    nPermutations = 10,
    networkType = "signed",
    quickCor = 0,
    verbose = 2,
    indent = 0
  )

  ## Zsummary tables
  mp.tables <- mp$preservation$Z[[1]][-reference]
  Z <- sapply(mp.tables, function(mat) mat[, "Zsummary.pres"])
  rownames(Z) <- rownames(mp.tables[[1]])
  rownames(Z) <- paste0("ME", rownames(Z))
  colnames(Z) <- names(multiExpr)[-reference]

  ## median rank
  mp.tables <- mp$preservation$observed[[1]][-reference]
  M <- sapply(mp.tables, function(mat) mat[, "medianRank.pres"])
  rownames(M) <- rownames(mp.tables[[1]])
  rownames(M) <- paste0("ME", rownames(M))
  colnames(M) <- names(multiExpr)[-reference]

  ## module size
  moduleSize <- mp.tables[[1]][, "moduleSize"]
  names(moduleSize) <- rownames(Z)

  ## module-traits. We need to recompute the MEs (module eigengenes)
  ## using the color coding of the reference set.
  refColors <- colorList[[1]]
  MEx <- lapply(exprList, function(x) {
    WGCNA::moduleEigengenes(t(x), colors = refColors)$eigengenes
  })

  ## Compute module-trait correlation matrices
  Y <- lapply(pres$layers, function(w) w$datTraits)
  names(Y)
  if ("Merged" %in% names(MEx) && !"Merged" %in% names(Y)) {
    kk <- rownames(MEx[["Merged"]])
    Y[["Merged"]] <- pres$datTraits[kk, ]
    Y <- Y[names(MEx)]
  }
  kk <- Reduce(union, lapply(Y, colnames))
  Y <- lapply(Y, function(y) y[, match(kk, colnames(y)), drop = FALSE])
  for (i in 1:length(Y)) colnames(Y[[i]]) <- kk
  R <- mapply(cor, MEx, Y, use = "pairwise", SIMPLIFY = FALSE)
  ## for(i in 1:length(R)) colnames(R[[i]]) <- paste0(names(R)[i],":",colnames(R[[i]]))

  ## gene statistics of reference layer
  if (compute.stats) {
    message("[wgcna.runPreservationWGCNA] computing gene statistics...")
    ref <- reference.name
    wnet <- list(MEs = MEx[[ref]], colors = pres$colors[, ref])
    pres$stats <- wgcna.computeGeneStats(wnet, pres$datExpr[[ref]],
      pres$datTraits,
      TOM = NULL
    )
  }

  ## geneset enrichment of reference layer
  if (compute.enrichment) {
    message("[wgcna.runPreservationWGCNA] computing geneset enrichment...")
    if (!is.null(GMT)) {
      GMT0 <- getPlaydataGMT()
      if (!is.null(annot)) GMT0 <- rename_by2(GMT0, annot, "human_ortholog")
      GMT <- merge_sparse_matrix(GMT, GMT0)
      remove(GMT0)
    } else {
      GMT <- getPlaydataGMT()
      if (!is.null(annot)) GMT <- rename_by2(GMT, annot, "human_ortholog")
    }

    ## we should check here if GMT and X overlap....
    pres$gsea <- wgcna.computeModuleEnrichment(
      pres$layers[[ref]],
      GMT = GMT,
      annot = annot,
      methods = gset.methods,
      ntop = 1000,
      xtop = 100,
      filter = NULL
    )
  }

  pres$modulePreservation <- mp
  pres$Zsummary <- Z
  pres$medianRank <- M
  pres$moduleSize <- moduleSize
  pres$modTraits <- R
  pres$MEs <- MEx

  return(pres)
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
