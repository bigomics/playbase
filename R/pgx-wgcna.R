##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
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
#' @return List containing WGCNA network and module results
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
    ngenes = 4000,
    maxBlockSize = 9999,
    gset.filter = "PATHWAY|HALLMARK|^GO|^C[1-9]",
    verbose = 1,
    progress = NULL
    ) {

  ##minmodsize=10;power=NULL;cutheight=0.15;deepsplit=2;ngenes=4000;networktype="signed";tomtype="signed";numericlabels=FALSE;ngenes=2000;gset.filter=NULL;minKME=0.8;maxBlockSize=5000

  samples <- pgx$samples
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
    group <- paste0(ct, ":", apply(pgx$contrasts, 1, paste, collapse = '_'))
    if ("batch" %in% colnames(samples))
      group <- paste0(group, ":", samples[, "batch"])
    nb <- round(ncol(counts)/500)
    message("[pgx.wgcna] running SuperCell. nb = ", nb)    
    sc <- pgx.supercell(counts, samples, group = group, gamma = nb)
    message("[pgx.wgcna] SuperCell done: ", ncol(counts), " ->", ncol(sc$counts))
    message("[pgx.wgcna] Normalizing supercell matrix (logCPM)")
    X <- as.matrix(logCPM(sc$counts, total = 1e4, prior = 1))
    samples <- sc$meta
    remove(counts, ct, group, nb, sc)
    gc()
  }
  
  if(!is.null(pgx$datatype) && pgx$datatype == "multi-omics") {
    message("[pgx.wgcna] Performing multi-omics ComBat on datatype.")
    X <- playbase::normalizeMultiOmics(X, method = "combat")
  }

  message("[pgx.wgcna] start wgcna.compute...")
  wgcna <- wgcna.compute(
    X = X,
    samples = samples,
    minmodsize = minmodsize, # default: min(20,...)
    power = power, # default: 12 (for signed)
    cutheight = cutheight, # default: 0.15
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
  
  ##---------------------------------------------------
  ## compute dimensionality reductions using TOM matrix
  ##---------------------------------------------------
  wTOM <- NULL
  if("TOM" %in% names(wgcna)) wTOM <- wgcna$TOM
  if(is.null(wTOM) && "svTOM" %in% names(wgcna)) wTOM  <- wgcna$svTOM %*% t(wgcna$svTOM)
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
  message("computing module enrichment...")
  res.gse <- wgcna.compute_enrichment(
    wgcna, pgx,
    method = c("fisher", "gsetcor", "xcor"),
    ntop = 1000,
    xtop = 100,
    filter = gset.filter
  )

  ## add to results object
  wgcna$gse <- res.gse
  wgcna$clust <- clust
  wgcna$networktype <- networktype
  wgcna$tomtype <- tomtype

  return(wgcna)
}




#' @export
wgcna.compute <- function(X,
                          samples,
                          ngenes = 4000,
                          minmodsize = 20,
                          power = 12,
                          #treeCut = NULL,
                          mergeCutHeight = 0.15,                          
                          deepsplit = 2,
                          minKME = 0.3,
                          treeCut = 0.99,
                          treeCutCeiling = 1,
                          networktype = "signed",
                          tomtype = "signed",
                          clustMethod = "average",
                          cutMethod = "hybrid",
                          calcMethod = "full",
                          lowrank = 40,
                          numericlabels = FALSE,
                          maxBlockSize = 8000,
                          merge.dendro = TRUE,
                          compute.stats = TRUE,
                          prefix = "ME",
                          sv.tom = 40,
                          drop.ref = TRUE,
                          verbose = 0
                          ) {
  
  #minmodsize=20;power=12;mergeCutHeight=0.15;deepsplit=2;ngenes=2000;networktype="signed";tomtype="signed";numericlabels=FALSE;prefix="ME";minKME=0.3;maxBlockSize=5000;calcMethod="fast";lowrank=40

  require(WGCNA)

  if(nchar(prefix)!=2) {
    stop("prefix must be two capital letters")
  }

  kk <- intersect(colnames(X),rownames(samples))
  X <- as.matrix(X[,kk])
  samples <- as.data.frame(samples, check.names=FALSE)
  samples <- samples[kk,,drop=FALSE]
  
  nmissing <- sum(is.na(X))
  if (nmissing > 0) {
    message("Found ", nmissing, " missing values in X. Imputing prior to WGCNA.")
    X <- svdImpute2(X)
  }

  X <- X[!duplicated(rownames(X)), ]
  
  ## restrict number of genes
  if (ngenes > 0 && nrow(X) > ngenes) {
    is.multiomics <- all(grepl(":",rownames(X)))
    if(is.multiomics) {
      message("[wgcna.compute] multiomics topSD = ",ngenes)
      X <- mofa.topSD(X, ngenes)    
    } else {
      message("[wgcna.compute] topSD = ",ngenes)
      sdx <- matrixStats::rowSds(X, na.rm = TRUE)
      X <- X[sdx > 0.1 * mean(sdx, na.rm = TRUE), ] ## filter low SD??
      X <- X[order(-matrixStats::rowSds(X, na.rm = TRUE)), ]
      X <- utils::head(X, ngenes)
    }
  }

  message("[wgcna.compute] dim(X) = ", paste(dim(X), collapse = " x "))
  message("[wgcna.compute] dim(samples) = ", paste(dim(samples), collapse = " x "))
  
  datExpr <- t(X)

  ## adapt for small datasets (also done in WGCNA package)
  minmodsize <- min(minmodsize, ncol(datExpr) / 2)

  message("[wgcna.compute] minmodsize = ", minmodsize)
  message("[wgcna.compute] number of features = ", nrow(X))
  message("[wgcna.compute] minKME = ", minKME)
  message("[wgcna.compute] power = ", power)
  message("[wgcna.compute] mergeCutHeight = ", mergeCutHeight)
  message("[wgcna.compute] calcMethod = ", calcMethod)  
  
  ##WGCNA::enableWGCNAThreads()
  cor <- WGCNA::cor ## needed...
  net <- wgcna.computeModules(
    datExpr,
    power = power,
    networkType = networktype,
    TOMType = tomtype,
    calcMethod = calcMethod,  
    lowrank = lowrank,
    clustMethod = clustMethod,
    cutMethod = cutMethod,
    deepSplit = deepsplit,
    minModuleSize = minmodsize,
    mergeCutHeight = mergeCutHeight,
    minKMEtoStay = minKME,
    treeCut = treeCut,
    treeCutCeiling = treeCutCeiling,      
    numericLabels = numericlabels,
    maxBlockSize = maxBlockSize,
    returnTOM = TRUE,
    verbose = verbose )
  cor <- stats::cor
  
  ## Substitue prefix="ME"
  if(prefix!="ME") names(net$MEs) <- sub("^ME", prefix, names(net$MEs))
  net$labels <- paste0(prefix, net$colors)
  
  ## clean up traits matrix
  datTraits <- data.frame(samples, check.names = FALSE)
  isdate <- apply(datTraits, 2, is.Date)
  datTraits <- datTraits[, !isdate, drop = FALSE]
    
  ## Expand multi-class discrete phenotypes into binary vectors
  datTraits <- utils::type.convert(datTraits, as.is = TRUE)
  datTraits <- expandPhenoMatrix(datTraits, keep.numeric=TRUE,
    drop.ref=drop.ref)

  ## list of genes in modules
  me.genes <- tapply(names(net$colors), net$colors, list)
  names(me.genes) <- paste0(prefix, names(me.genes))
  me.genes <- me.genes[names(net$MEs)]

  ## get colors of eigengene modules
  color1 <- wgcna.labels2colors(net$colors)
  me.colors <- color1[!duplicated(color1)]
  me.labels <- net$labels[!duplicated(color1)]
  names(me.colors) <- me.labels
  me.colors <- me.colors[names(net$MEs)]

  ## compute TOM matrix (we need for some plots)
  if(!is.null(net$TOM)) {
    TOM <- net$TOM
    net$TOM <- NULL
  } else {
    message("[wgcna.compute] recomputing TOM matrix...")
    TOM <- WGCNA::TOMsimilarityFromExpr(
      datExpr,
      power = power,
      TOMType = tomtype,
      networkType = networktype,
      verbose = verbose
    )
  }

  ## instead of the huge TOM matrix we save a smaller SVD.
  svTOM <- NULL
  if(sv.tom>0) {
    ## sv.tom <- ceiling(min(sv.tom,dim(datExpr)/2))
    message("[wgcna.compute] reducing TOM. sv.tom = ", sv.tom)
    rownames(TOM) <- colnames(TOM) <- colnames(datExpr)
    sv.tom <- min(sv.tom, ncol(TOM)-1)
    sv <- irlba::irlba(TOM, nv = sv.tom)
    svTOM <- sv$v %*% diag(sqrt(sv$d))
    rownames(svTOM) <- colnames(datExpr)
  }
  
  ## compute module eigenvectors (loading matrix)
  message("[wgcna.compute] computing module eigenvectors...")
  MVs <- list()
  for (clr in unique(net$colors)) {
    ii <- which(net$colors == clr)
    mX <- datExpr[, ii, drop=FALSE]
    mX <- scale(mX) ## NOTE: seems WGCNA is using full scaling
    if(ncol(mX)>1) {
      res <- irlba::irlba(mX, nv = 1, nu = 1)
      sv1 <- res$v[,1] * res$d[1]
    } else {
      sv1 <- 1
    }
    sv <- rep(0, ncol(datExpr))
    sv[ii] <- sv1
    MVs[[paste0(prefix, clr)]] <- sv
  }
  MVs <- as.matrix(do.call(cbind, MVs[names(net$MEs)]))
  rownames(MVs) <- colnames(datExpr)

  ## compute gene statistics
  stats <- NULL
  if(compute.stats) {
    message("[wgcna.compute] computing gene statistics...")
    stats <- wgcna.compute_geneStats(net, datExpr, datTraits, TOM=TOM) 
  }
  
  ## module-traits matrix
  message("[wgcna.compute] computing module-traits matrix...")
  modTraits <- cor(net$MEs, datTraits, use="pairwise")
  
  ## merge dendrograms
  if(merge.dendro) {
    message("[wgcna.compute] merge_block_dendrograms...")  
    merged <- try(wgcna.merge_block_dendrograms(net, TOM))
    if(!inherits(merged,"try-error")) {
      net$merged_dendro <- merged
    } else {
      net$merged_dendro <- NULL
    }
  }
  remove(TOM)  ## big
  
  ## construct results object
  results <- list(
    datExpr = datExpr,
    datTraits = datTraits,
    ## TOM = TOM,  ## this can be BIG!!! generally no need, just for plotting
    svTOM = svTOM,  ## smaller singular vectors
    net = net,
    #power = power,
    me.genes = me.genes,
    me.colors = me.colors,
    W = MVs,
    modTraits = modTraits,
    stats = stats
  )

  return(results)
}

#' @export
wgcna.compute_multiomics <- function(dataX,
                                     samples,
                                     power = 12,
                                     ngenes = 4000,
                                     clustMethod = "average",
                                     cutMethod = "hybrid",
                                     minmodsize = 10,
                                     minKME = 0.3,
                                     deepsplit = 2,
                                     compute.enrichment = TRUE,
                                     xtop = 100,
                                     annot = NULL,
                                     GMT = NULL,
                                     gsetX = NULL,
                                     verbose = 1,
                                     progress = NULL
                                     ) {

  ## preprocessing
  dataX <- lapply(dataX, imputeMissing, method="SVD2")
  names(dataX) <- substring(names(dataX),1,2)
  lapply(dataX, dim)
  dataX <- mofa.topSD(dataX, ngenes)

  if(!is.null(progress)) {
    progress$set(message = paste("computing WGCNA modules..."), value = 0.33)
  }
  
  ## This runs WGCNA on an expression list. 
  wgcna <- list()
  for(dt in names(dataX)) {
    cat("******** computing WGCNA for",dt,"********\n")
    
    minKME <- ifelse(dt=='ph', 0, minKME)
    minmodsize <- ifelse( dt=='ph', 1, minmodsize)      
    
    wgcna[[dt]] <- wgcna.compute(
      X = dataX[[dt]],
      samples = samples, 
      ngenes = -1,
      calcMethod = "fast",
      power = power,
      lowrank = 40,
      clustMethod = clustMethod,
      cutMethod = cutMethod,
      deepsplit = deepsplit,
      minKME = minKME,
      minmodsize = minmodsize,
      mergeCutHeight = 0.15,
      compute.stats = TRUE,
      sv.tom = 40,  
      prefix = toupper(dt),
      verbose = verbose
    )
  }

  ## Compute enrichment
  if(compute.enrichment) {
    message("computing module enrichment...")
    if(!is.null(progress)) {
      progress$set(message = paste("computing module enrichment..."), value = 0.66)
    }

    gse <- wgcna.compute_enrichment(
      wgcna = wgcna,
      multi = TRUE,
      pgx = NULL,
      method = c("fisher","gsetcor","xcor"),
      ntop = 400,
      xtop = xtop,
      annot = annot,      
      GMT = GMT,
      gsetX = gsetX,
      filter = NULL
      #filter = "^PATHWAY|^HALLMARK|^GO|^C[1-9]"  
    )

    ## split up results
    for(k in names(wgcna)) {
      mm <- names(wgcna[[k]]$me.genes)
      wgcna[[k]]$gsea <- gse[mm]
    }
  } ## end if-compute-enrichment

  return(wgcna)
}


#' Reimplementation for WGCNA::blockwiseModules(). This returns exact
#' same object as WGCNA::blockwiseModules() but is faster and allows
#' different clustering linkage methods (ward, complete).
#' 
wgcna.computeModules <- function(
  datExpr,
  power = 6,
  networkType = "signed",
  TOM = NULL,
  TOMType = "signed",
  calcMethod = "full",  
  lowrank = 20,
  clustMethod = "average",
  cutMethod = "hybrid",  ## hybrid, tree, static
  deepSplit = 2,
  treeCut = 0.99,
  treeCutCeiling = 1,  
  minModuleSize = 20,
  mergeCutHeight = 0.15,
  minKMEtoStay = 0.3,
  numericLabels = FALSE, ## numeric or 'color' labels
  maxBlockSize = 9999,
  returnTOM = FALSE,
  verbose = 1 ) {

  #power=6;networkType=TOMType="signed";minModuleSize=20;mergeCutHeight=0.15;minKMEtoStay=0.3;numericLabels=FALSE;clustMethod="average";deepSplit = 2;treeCut = 0.99;treeCutCeiling = 1;
  
  cor <- WGCNA::cor  ## needed...

  if (is.null(power)) {
    ## Estimate best power
    message("[wgcna.compute] estimating optimal power...")
    powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
    powers <- c(powers, seq(from = 25, to = 50, by = 5))
    power <- wgcna.pickSoftThreshold(datExpr, sft=NULL, rcut=0.85,
      plot=FALSE, verbose=0) 
    if (is.na(power)) power <- 6
  }
  
  if(calcMethod == "blockwise") {
    message("[wgcna.compute] computing blockwiseModules...")
    net <- WGCNA::blockwiseModules(
      datExpr,
      power = power,
      networkType = networkType,
      TOMType = TOMType,
      minModuleSize = minModuleSize,
      mergeCutHeight = mergeCutHeight,
      minKMEtoStay = minKMEtoStay,
      numericLabels = numericLabels, ## numeric or 'color' labels
      deepSplit = deepSplit,
      maxBlockSize = maxBlockSize,
      verbose = verbose
    )
    return(net)
  }
  
  clustMethod <- sub("^ward$","ward.D",clustMethod)
  ## define distance matrix
  if(is.null(TOM)) {
    adjacency <- WGCNA::adjacency(datExpr, power = power, type = networkType) 
    adjacency[is.na(adjacency)] <- 0
    
    if(calcMethod == "fast") {
      if(verbose>0) message("Computing TOM matrix using fast method...")
      TOM <- fastTOMsimilarity(adjacency, tomtype=TOMType,
        lowrank=lowrank)
    } else if(calcMethod == "adjacency") {
      if(verbose>0) message("Computing using adjacency as TOM matrix...")    
      TOM <- adjacency
    } else if(calcMethod == "full") {
      if(verbose>0) message("Computing full TOM matrix...")
      ## SLOW!!!
      TOM <- WGCNA::TOMsimilarity(adjacency, TOMType=TOMType, verbose=verbose) 
    } else {
      stop("ERROR: invalid calcMethod parameter:", calcMethod)
    }
    dimnames(TOM) <- dimnames(adjacency)
  }

  ## transform to dissTOM
  dissTOM <- 1 - TOM

  ## clustering
  if(verbose>0) message("Clustering features using ", clustMethod, " linkage")
  ##geneTree <- flashClust::flashClust(as.dist(dissTOM), method=clustMethod)
  geneTree <- stats::hclust(as.dist(dissTOM), method=clustMethod)  

  ## sometimes there is a height error. following is a fix.
  geneTree$height <- round(geneTree$height, 6) 
  
  ## exception1
  if(minModuleSize <= 1 && cutMethod!="static") {
    message("WARNING: minModuleSize==1. Changing to static cutting")
    cutMethod <- "static"
  }

  if(treeCut > 1) cutMethod <- "static"  
  if(treeCut <= 1) {
    ## transform from relative to actual
    qq <- quantile( geneTree$height, probs = c(0.05, treeCutCeiling))
    treeCut <- qq[1] + treeCut * diff(qq)
  }

  ## exception2
  extest <- mad(geneTree$height) / median(geneTree$height)
  if(extest < 0.01 ) {
    cutMethod <- "static"
    treeCut <- 0.9 * min(geneTree$height)
    message("WARNING: very low similarity. Switching to static cut h=", treeCut)
  }
  
  if(cutMethod %in% c("hybrid","tree")) {

    if(cutMethod=="tree") deepSplit <- (deepSplit > 0)

    if(verbose>0) message("Dynamic cutting method: ", cutMethod)
    if(verbose>0) message("treeCut = ", treeCut)
    if(verbose>0) message("deepSplit = ", deepSplit)    

    label <- dynamicTreeCut::cutreeDynamic(
      geneTree,
      method = cutMethod,
      cutHeight = treeCut,
      distM = dissTOM,
      deepSplit = deepSplit,
#      pamRespectsDendro = FALSE,
      minClusterSize = minModuleSize
    )
  } else if( cutMethod == "static" && treeCut <= 1 ) {
    if(verbose>0) message("Static cutting tree at fixed H = ", treeCut)    
    label <- cutree( geneTree, h = treeCut)
  } else if(cutMethod == "static" && treeCut > 1 ) {
    if(verbose>0) message("Static cutting tree with fixed K = ", treeCut)
    label <- cutree( geneTree, k = treeCut)
  } else {
    stop("ERROR: could not determine cutMethod")
  } 
  label <- as.integer(label)
  table(label)
  nmodules <- length(unique(label))
  if(verbose>0) message("Found ", nmodules, " modules")
  
  ## Eigengenes
  if(verbose>0) message("Calculating eigengenes...")
  colors <- WGCNA::labels2colors(label)
  MEs = WGCNA::moduleEigengenes(datExpr, colors = colors)$eigengenes
  table(colors)
  
  ## prune using minKME
  if(minKMEtoStay > 0) {
    if(verbose>0) message("Pruning features using minKME = ",minKMEtoStay)
    ngrey <- sum(colors == 'grey')
    for(k in unique(colors)) {
      ii <- which(colors == k)
      if(length(ii)>1) {
        eg <- MEs[,paste0("ME",k)]
        kme <- cor(datExpr[,ii, drop=FALSE], eg, use="pairwise")[,1]
        kme.sign <- as.numeric(names(which.max(table(sign(kme)))))
        kme <- kme * kme.sign
        ii <- ii[which(kme < minKMEtoStay)]
        if(length(ii)) colors[ii] <- 'grey'
      }
    }
    ngrey <- sum(colors == 'grey') - ngrey
    if(verbose>0) message("Pruned ",ngrey, " low KME features")    
  }
  
  ## Merge similar modules
  unmergedColors <- colors
  if(mergeCutHeight>0 && length(MEs)>1) {
    if(verbose>0) message("Merging similar modules: mergeCutHeight = ",mergeCutHeight)
    #merge <- WGCNA::mergeCloseModules(datExpr, colors, cutHeight=mergeCutHeight, verbose=0)
    merge <- wgcna.mergeCloseModules(datExpr, colors, cutHeight=mergeCutHeight, MEs=MEs)
    unmergedColors <- colors
    colors <- merge$colors
    MEs <- merge$MEs
    if(verbose>0) {
      n0 <- length(unique(unmergedColors))
      n1 <- length(unique(colors))      
      message("Merged ", n0 - n1, " modules")
    }
  } 

  ## filter on minModuleSize
  if(minModuleSize>1) {
    too.small <- names(which(table(colors) < minModuleSize))
    too.small
    if(length(too.small)) {
      if(verbose>0) message("Removing ",length(too.small)," too small modules")
      colors[colors %in% too.small] <- "grey"
    }
  }

  ## Update MEs
  if("grey" %in% colors && !"MEgrey" %in% names(MEs)) {
    MEs <- WGCNA::moduleEigengenes(datExpr, colors = colors)$eigengenes
  }
  MEs <- MEs[sub("^ME","",names(MEs)) %in% colors]
  
  # Rename to numeric
  if(numericLabels) {
    if(verbose>0) message("Renaming to numeric labels")    
    colorOrder <- names(sort(table(colors),decreasing=TRUE))
    colorOrder <- unique(c("grey", colorOrder))
    colors <- match(colors, colorOrder)-1
    unmergedColors <- match(unmergedColors, colorOrder)-1
    mecolor <- sub("^ME","",names(MEs))
    names(MEs) <- paste0("ME",match(mecolor, colorOrder)-1)
  } else {
    # Rename to standard colors, most frequent first
    if(verbose>0) message("Renaming to standard colors")
    colorOrder <- names(sort(table(colors),decreasing=TRUE))
    colorOrder <- unique(c("grey", colorOrder))
    newcolor <- setdiff( WGCNA::standardColors(), "grey")
    n0 <- length(colorOrder)
    n1 <- length(newcolor)
    if(n1 < n0) newcolor <- make.unique(rep(newcolor,ceiling(n0/n1)))
    newcolor <- unique(c("grey", newcolor))
    colors <- newcolor[match(colors, colorOrder)]
    unmergedColors <- newcolor[match(unmergedColors, colorOrder)]
    mecolor <- sub("^ME","",names(MEs))
    names(MEs) <- paste0("ME",newcolor[match(mecolor, colorOrder)])
  } 
  
  names(colors) <- colnames(datExpr)
  names(unmergedColors) <- colnames(datExpr)
  
  net <- list()
  net$colors <- colors
  net$unmergedColors <- unmergedColors
  net$MEs <- MEs
  net$goodSamples <- rep(TRUE, nrow(datExpr))
  net$goodGenes <- rep(TRUE, ncol(datExpr))
  net$dendrograms <- list( geneTree )
  net['TOMFiles'] <- list(NULL)
  net$blockGenes <- list(1:ncol(datExpr))
  net$blocks <- rep(1, ncol(datExpr))
  net$MEsOK <- TRUE
  net$power <- power
  if(returnTOM) net$TOM <- TOM
  return(net)
}

#' Faster implementation of TOM computation using low-rank SVD
#' approximation.
#' 
fastTOMsimilarity <- function(A, tomtype="signed", lowrank=20) {
  #https://stackoverflow.com/questions/56574729
  #
  #Given square symmetric adjacency matrix A, its possible to
  #calculate the TOM matrix W without the use of for loops, which
  #speeds up the process tremendously
  if(!tomtype %in% c("signed","unsigned")) {
    stop("only works for signed and unsigned tomtype")
  }
  
  ## Adjacency matrix A can be approximated with SVD. This can make
  ## TOM calculation much faster.
  diag(A) <- 0
  if(ncol(A) < 2*lowrank) lowrank <- -1
  if(lowrank>0) {
    res <- irlba::irlba(A, nv=lowrank)
    U <- res$u %*% diag(sqrt(res$d))
    L <- U %*% (Matrix::t(U) %*% U) %*% Matrix::t(U)
  } else {
    L <- A %*% A
  }
  k <- Matrix::colSums(A)
  Kmat <- outer(k, k, FUN = function(x, y) pmin(x, y))
  D <- (Kmat + 1 - A)
  W <- (L + A) / D
  diag(W) <- 1
  W <- as.matrix(W)  
  dimnames(W) <- dimnames(A)
  W <- pmax(W, 0)  ## sometimes has negative values...
  return(W)
}


wgcna.mergeCloseModules <- function(datExpr, colors, cutHeight, MEs=NULL) {
  if(is.null(MEs)) {
    MEs = WGCNA::moduleEigengenes(datExpr, colors = colors)$eigengenes
    dim(MEs)
  }
  hc <- hclust(as.dist(1 - cor(MEs)), method="average")
  idx <- cutree(hc, h=cutHeight)
  names(idx) <- sub("^ME","",names(idx))
  table(idx)
  new.colors <- colors
  m=2
  for(m in unique(idx)) {
    cc <- names(which(idx == m))
    cc <- setdiff(cc, "grey") ## never merge grey
    ii <- which(colors %in% cc)
    new.colors[ii] <- cc[1]
  }
  new.MEs = WGCNA::moduleEigengenes(datExpr, colors = new.colors)$eigengenes
  list(
    colors = new.colors,
    MEs = new.MEs
  )
}

##---------------------------------------------------------------------
## Gene statistics
##---------------------------------------------------------------------

wgcna.compute_geneStats <- function(net, datExpr, datTraits, TOM) {
  ## Define numbers of genes and samples
  nGenes <- ncol(datExpr)
  nSamples <- nrow(datExpr)

  ## Recalculate MEs with color labels
  moduleTraitCor <- cor(net$MEs, datTraits, use = "pairwise.complete")
  moduleTraitPvalue <- WGCNA::corPvalueStudent(moduleTraitCor, nSamples)

  ## Module membership correlation (with p-values)
  moduleMembership <- cor(datExpr, net$MEs, use = "pairwise.complete")
  MMPvalue <- WGCNA::corPvalueStudent(as.matrix(moduleMembership), nSamples)

  ## Gene-trait significance (trait correlation) (with p-values)
  traitSignificance <- cor(datExpr, datTraits, use = "pairwise.complete")
  GSPvalue <- WGCNA::corPvalueStudent(as.matrix(traitSignificance), nSamples)

  ## Fold-change
  foldChange <- NULL
  foldChangePvalue <- NULL  
  is.binary <- apply(datTraits, 2, function(a) length(unique(a[!is.na(a)])) == 2)
  if(any(is.binary)) {
    binY <- datTraits[, which(is.binary), drop = FALSE]
    lm <- apply(binY, 2, function(y) {
      gx.limma(t(datExpr), y,
        lfc = 0, fdr = 1,
        sort.by = "none", verbose = 0
      )
    })
    foldChange <- sapply(lm, function(m) m$logFC)
    foldChangePvalue <- sapply(lm, function(m) m$P.Value)
    if(length(lm)==1) {
      foldChange <- cbind(foldChange)
      foldChangePvalue <- cbind(foldChangePvalue)
    }
    rownames(foldChange) <- rownames(lm[[1]])
    rownames(foldChangePvalue) <- rownames(lm[[1]])
    colnames(foldChange) <- colnames(binY)
    colnames(foldChangePvalue) <- colnames(binY)    
  }

  # Continuous traits (not always present)
  contY <- datTraits[, which(!is.binary), drop = FALSE]
  foldChange.cont <- NULL
  foldChangePvalue.cont <- NULL  
  if(NCOL(contY) > 0) {
    contlm <- apply(contY, 2, function(y) {
      rho <- cor(datExpr, y, use="pairwise")[,1]
      #P.Value <- cor.pvalue(rho, length(y))
      P.Value <- WGCNA::corPvalueStudent(rho, n = length(y))
      data.frame(rho, P.Value)
    })
    foldChange.cont <- sapply(contlm, function(m) m$rho)
    foldChangePvalue.cont <- sapply(contlm, function(m) m$P.Value)
    if(length(contlm)==1) {
      foldChange.cont <- cbind(foldChange.cont)
      foldChangePvalue.cont <- cbind(foldChangePvalue.cont)
    }
    rownames(foldChange.cont) <- rownames(contlm[[1]])
    rownames(foldChangePvalue.cont) <- rownames(contlm[[1]])
    colnames(foldChange.cont) <- colnames(contY)
    colnames(foldChangePvalue.cont) <- colnames(contY)    
  }

  # Merge
  foldChange <- cbind(foldChange, foldChange.cont)
  foldChangePvalue <- cbind(foldChangePvalue, foldChangePvalue.cont)

  ## Gene Centrality. Compute centrality of gene in Module subgraph
  ## using TOM matrix.
  if(is.null(dimnames(TOM))) dimnames(TOM) <- list(colnames(datExpr),colnames(datExpr))
  adj <- TOM
  diag(adj) <- NA
  adj[which(abs(adj) < 0.01)] <- 0
  gr <- igraph::graph_from_adjacency_matrix(
    adj, mode = "undirected", weighted = TRUE, diag = FALSE
  )
  geneCentrality <- rep(NA, nrow(adj))
  names(geneCentrality) <- rownames(adj)
  me.genes <- tapply(names(net$colors), net$colors, list)
  gg <- me.genes[[1]]
  for (gg in me.genes) {
    gr1 <- igraph::subgraph(gr, gg)
    ct <- igraph::page_rank(gr1, weights = NULL)$vector
    ct <- ct / mean(ct, na.rm = TRUE)
    geneCentrality[gg] <- ct
  }

  ## propotion variance explained (PVE) per module
  propVarExplained <- WGCNA::propVarExplained(
    datExpr,
    colors = net$colors, MEs = net$MEs
  )

  ## force align. Sometime they are shorter for some reason...
  gg <- rownames(moduleMembership)
  matMatch <- function(m) {
    m <- m[match(gg,rownames(m)),,drop=FALSE]
    rownames(m) <- gg
    return(m)
  }
  ##moduleMembership <- matMatch(moduleMembership)
  MMPvalue <- matMatch(MMPvalue)
  traitSignificance = matMatch(traitSignificance)
  GSPvalue = matMatch(GSPvalue)
  foldChange = matMatch(foldChange)
  foldChangePvalue = matMatch(foldChangePvalue)
  
  stats <- list(
    moduleTraitCor = moduleTraitCor,
    moduleTraitPvalue = moduleTraitPvalue,
    moduleMembership = moduleMembership,
    MMPvalue = MMPvalue,
    traitSignificance = traitSignificance,
    GSPvalue = GSPvalue,
    foldChange = foldChange,
    foldChangePvalue = foldChangePvalue,
    geneCentrality = geneCentrality
  )

  return(stats)
}


## ----------------------------------------------------
## Perform geneset analysis on modules
## ----------------------------------------------------

wgcna.compute_enrichment <- function(wgcna,
                                     pgx,
                                     multi = FALSE,
                                     method = c("fisher","gsetcor","xcor"),
                                     ntop = 200,
                                     xtop = 100,
                                     annot = NULL,
                                     GMT = NULL,
                                     gsetX = NULL,
                                     filter = "^PATHWAY|^HALLMARK|^GO|^C[1-9]") {


  if(!multi) {
    wgcna <- list(gx = wgcna)
  }

  if(!any( c("gx","px") %in% names(wgcna))) {
    message("ERROR: datasets must have gx or px datatype!")
    return(NULL)
  }
  
  ## collapse features to symbol
  k <- intersect(c("gx","px"), names(wgcna))[1]
  geneX <- t(as.matrix(wgcna[[k]]$datExpr))
  symbol.col <- NULL
  if(!is.null(pgx)) {
    GMT <- pgx$GMT
    gsetX <- pgx$gsetX
    annot <- pgx$genes
  }

  ## limit xtop if geneX is too small
  xtop <- min(xtop, round(nrow(geneX)/3))
  
  if(is.null(annot)) {
    gg <- lapply(wgcna, function(w) colnames(w$datExpr))
    gg1 <- lapply(names(wgcna), function(a) paste0(a,":",gg[[a]]))
    gg <- as.character(unlist(gg))
    gg1 <- as.character(unlist(gg1))    
    annot <- data.frame(feature = gg1, symbol = gg)
  }
  
  symbol.col <- intersect(c("symbol","gene_name"),colnames(annot))[1]
  symbol.col
  geneX <- rename_by2(geneX, annot, symbol.col)
  
  if(is.null(GMT)) {
    message("Using playdata GSETxGENE genesets")
    GMT <- Matrix::t(playdata::GSETxGENE)
    if(!is.null(annot)) GMT <- rename_by2( GMT, annot, symbol.col)
  }

  bg <- rownames(geneX)
  bg <- intersect(bg, rownames(GMT))
  if(length(bg)==0) {
    message("FATAL. no overlapping genes")
    return(NULL)
  }
  G1 <- GMT[bg,, drop=FALSE]
  if(!is.null(filter)) {
    sel <- grep(filter,colnames(G1))
    if(length(sel)) G1 <- G1[,sel, drop=FALSE]
  }
  G1 <- G1[, which(Matrix::colSums(G1 != 0) >= 4), drop = FALSE]

  if(is.null(gsetX)) {
    gsetX <- plaid::plaid( geneX, G1)
  }

  ## align dimensions
  ss <- intersect(rownames(gsetX), colnames(G1))
  G1 <- G1[,ss,drop=FALSE]
  gsetX <- gsetX[ss,]
  gmt <- mat2gmt(G1)

  dbg("[wgcna.compute_enrichment] dim.G1 = ",dim(G1))
  dbg("[wgcna.compute_enrichment] len.gmt = ",length(gmt))
  dbg("[wgcna.compute_enrichment] dim.gsetX = ",dim(gsetX))
  
  ME <- lapply(wgcna, function(w) as.matrix(w$net$MEs))
  ME <- do.call(cbind, ME)
  me.genes <- lapply(wgcna, function(w) w$me.genes)
  names(me.genes) <- NULL
  me.genes <- do.call(c, me.genes)
  me.genes <- lapply(me.genes, function(g) probe2symbol(g, annot, query = symbol.col))

  ## compute most closest gx/px genes
  nbx.genes <- list()
  for(k in 1:length(wgcna)) {
    mx <- wgcna[[k]]$net$MEs
    nbx.cor <- cor(t(geneX), mx)
    for(i in colnames(nbx.cor)) {
      if(is.null(xtop)) {
        n <- length(wgcna[[k]]$me.genes[[i]])
      } else {
        n <- xtop
      }
      nbx.genes[[i]] <- head(names(sort(-nbx.cor[,i])), n)
    }
  }

  ## add nearest expression neighbors to module genes
  nbx.genes <- nbx.genes[names(me.genes)]
  for(i in names(me.genes)) {
    me.genes[[i]] <- unique(c(me.genes[[i]], nbx.genes[[i]]))
  }
  
  rho.list <- list()
  pval.list <- list()

  ## Here we correlate geneset score (averageCLR) with the module
  ## eigengene (ME). This should select genesets correlated with the
  ## ME.
  if ("gsetcor" %in% method && !is.null(gsetX)) {
    message("[wgcna.compute_enrichment] calculating single-sample geneset correlation...")
    rc.rho <- matrix(NA, ncol(G1), ncol(ME))
    rc.pvalue <- matrix(NA, ncol(G1), ncol(ME))    
    dimnames(rc.rho) <- list( colnames(G1), colnames(ME))
    dimnames(rc.pvalue) <- list( colnames(G1), colnames(ME))    
    jj <- which(rownames(gsetX) %in% colnames(G1))
    kk <- intersect(colnames(gsetX), rownames(ME)) ## common samples
    rho.jj <- cor(t(gsetX[jj,kk]), ME[kk,], use = "pairwise")
    #rc.pvalue <- cor.pvalue(rc.rho, n = length(kk))
    pvalue.jj <- WGCNA::corPvalueStudent(rho.jj, n = length(kk))
    ii <- match(rownames(gsetX)[jj], rownames(rc.rho))
    rc.rho[ii,] <- rho.jj
    rc.pvalue[ii,] <- pvalue.jj    
    rho.list[["gsetcor"]] <- rc.rho
    pval.list[["gsetcor"]] <- rc.pvalue
  }

  ## Here we correlate the module eigengene (ME) with genes and then
  ## do a gset.rankcor() on the ME correlation.
  if ("xcor" %in% method) {
    message("[wgcna.compute_enrichment] calculating Eigengene correlation...")
    mm <- cor(t(geneX), ME, use = "pairwise")
    rc <- gset.rankcor(mm, G1, compute.p = TRUE) ## NEEDS CHECK!!!
    rho.list[["xcor"]] <- rc$rho
    pval.list[["xcor"]] <- rc$p.value
  }

  ## Perform fisher-test on ME genes
  if ("fisher" %in% method) {
    message("[wgcna.compute_enrichment] calculating Fisher tests...")
    
    ## NOTE: should we pre-select to make this faster?? prioritize
    ## upon previous tests?
    sel <- 1:length(gmt)
    if(0) {
      message("[wgcna.compute_enrichment] selecting sets...")
      gsets <- Reduce(intersect, lapply( rho.list, rownames))
      modules <- Reduce(intersect, lapply( rho.list, colnames))  
      rho.list <- lapply( rho.list, function(x) x[gsets,modules])
      pval.list <- lapply( pval.list, function(x) x[gsets,modules])  
      pval.NAZERO <- lapply(pval.list, function(x) {x[is.na(x)]=0;x})
      meta.p <- Reduce(pmax, pval.NAZERO) ## NEED RETHINK!!!
      ## select maximum 5*ntop genesets
      sel <- head(order(apply(meta.p, 1, min, na.rm=TRUE)), 5*ntop)
    }    
    rho <- matrix(NA, ncol(G1), ncol(ME))
    pval <- matrix(NA, ncol(G1), ncol(ME))    
    dimnames(rho) <- list( colnames(G1), colnames(ME))
    dimnames(pval) <- list( colnames(G1), colnames(ME))    

    i <- 1
    for (i in 1:length(me.genes)) {
      gg <- me.genes[[i]]      
      ## perform Fisher test
      rr <- try(gset.fisher(gg, gmt[sel], background = bg, fdr = 1,
        min.genes = 4, verbose = 0))
      if (!"try-error" %in% class(rr)) {
        rr <- rr[match(rownames(rho), rownames(rr)), ]
        rho[,i] <- rr$odd.ratio
        pval[,i] <- rr$p.value
      }
    }

    ## match dimensions with prevous
    rho <- rho[, match(names(me.genes),colnames(rho))]
    pval <- pval[, match(names(me.genes),colnames(pval))]    

    ## handle infinite
    rho[is.infinite(rho)] <- 2 * max(rho, na.rm = TRUE)  ## Inf odd.ratio
    rho.list[["fisher"]] <- rho
    pval.list[["fisher"]] <- pval
  }

  ## ensure dimensions
  gsets <- Reduce(intersect, lapply( rho.list, rownames))
  modules <- Reduce(intersect, lapply( rho.list, colnames))  
  rho.list <- lapply( rho.list, function(x) x[gsets,modules])
  pval.list <- lapply( pval.list, function(x) x[gsets,modules])  
  
  ## Compute meta rank and pval. Handle NA for failing methods.
  pval.NAZERO <- lapply(pval.list, function(x) {x[is.na(x)]=0;x})
  meta.p <- Reduce(pmax, pval.NAZERO) ## NEED RETHINK!!!
  meta.q <- apply(meta.p, 2, p.adjust, method = "fdr")

  ## NEED RETHINK: how about negative FC???
  rnk.list <- lapply(rho.list, function(x) apply(x, 2, rank, na.last='keep') / nrow(x))
  meta.rnk <- Reduce("+", rnk.list) / length(rnk.list)
  rnk.NAZERO  <- lapply(rnk.list, function(x) {x[is.na(x)]=0;x})
  rnk.NSUM <- Reduce('+', lapply(rnk.list, function(x) !is.na(x)))
  meta.rnk <- Reduce('+', rnk.NAZERO) / rnk.NSUM

  ## create dataframe by module
  message("[wgcna.compute_enrichment] creating dataframes...")
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
    df <- df[order(-df$score), ]
    if (!is.null(ntop) && ntop > 0) {
      sel <- c(head(order(-df$score),ntop), tail(order(-df$score),ntop))
      df <- df[unique(sel),]
    }
    gse.list[[k]] <- df
  }

  ## add genes
  gse.genes <- list()
  k <- names(gse.list)[1]
  for (k in names(gse.list)) {
    gset <- rownames(gse.list[[k]])
    gg <- me.genes[[k]]
    set.genes <- lapply(gmt[gset], function(s) intersect(s, gg))
    n0 <- sapply(gmt[gset], length)
    n1 <- sapply(set.genes, length)
    gse.genes[[k]] <- sort(table(unlist(set.genes)), decreasing = TRUE)
    set.genes <- sapply(set.genes, function(g) paste(sort(g), collapse = "|"))
    gse.list[[k]]$overlap <- paste0(n1, "/", n0)
    gse.list[[k]]$genes <- set.genes
  }

  return(gse.list)
}



#'
#'
#' @export
wgcna.getGeneStats <- function(wgcna, trait, module=NULL, plot = TRUE,
                               showallmodules= is.null(module),
                               col = NULL, main = NULL) {

  if(is.null(module)) showallmodules <- TRUE
  
  p1 <- c("moduleMembership", "MMPvalue")
  p2 <- c("traitSignificance", "GSPvalue", "foldChange", "foldChangePvalue")
  p3 <- c("geneCentrality")
  nrow <- ncol(wgcna$datExpr)
  ##df <- wgcna$stats[[p1[1]]][, 0] ## empty data.frame
  df <- data.frame(
    feature = colnames(wgcna$datExpr),
    module = wgcna$net$labels
  )
  rownames(df) <- colnames(wgcna$datExpr)
  
  ## get moduleMembership
  mm.stats <- wgcna$stats[p1]
  labels <- wgcna$net$labels
  idx <- cbind(1:length(labels), match(labels, colnames(mm.stats[[1]])))
  A1 <- sapply( mm.stats, function(x) x[idx])  
  df <- cbind(df, A1)

  
  ## get traitSig columns for trait
  tt.cols <- colnames(wgcna$stats[[p2[1]]])  
  if(is.null(trait)) trait <- tt.cols
  trait <- intersect(trait, tt.cols)

  if (length(trait)) {
    A2 <- sapply(wgcna$stats[p2], function(x) x[, trait])
    df <- cbind(df, A2)
  }

  
  A3 <- wgcna$stats[[p3]]
  df <- cbind(df, centrality = A3)

  ## calculate score
  sel <- c("moduleMembership", "traitSignificance", "foldChange", "centrality")
  sel <- intersect(sel, colnames(df))
  #df1 <- pmax(as.matrix(df[, sel]), 1e-8)
  df1 <- as.matrix(abs(df[, sel]))
  score <- exp(rowMeans(log(1e-8 + df1))) * sign(df[,"foldChange"])

  df <- data.frame(df[,1:2], score = score, df[,-c(1,2)])
  df <- df[order(-df$score), ]
  if (!is.null(module) && !showallmodules) {
    sel <- which(df$module == module)
    df <- df[sel, , drop = FALSE]
  }
  
  if (plot) {
    cols <- c("moduleMembership", "traitSignificance", "foldChange", "centrality")
    cols <- intersect(cols, colnames(df))
    df1 <- df[, cols]
    col1 <- wgcna.labels2colors(wgcna$net$colors[rownames(df1)])
    if (!is.null(col)) col1 <- col
    pairs(df1, col = col1)
    if (is.null(main)) {
      main <- paste("Gene significance for module", module, "and trait", trait)
    }
    title(main, line = 3, cex.main = 1.15)
  }

  df
}

wgcna.merge_block_dendrograms <- function(net, X, method=1) {
  ## This function is fragile: it can give a C stack limit error. In
  ## case that happens you can increase the stack limit by running on
  ## the cmd line: >>> ulimit -s unlimited
  ## 
  hc = net$dendrograms
  
  ## merge block dendrogram
  mx <- list()
  for(b in 1:length(net$dendrograms)) {
    ii <- which(net$goodGenes & net$blocks==b)
    mx[[b]] <- colMeans(X[ii,])
    hc[[b]]$labels <- rownames(X)[ii]
  }

  if(length(hc)==1) {
    return(hc[[1]])
  }

  ## compute parent dendrogram
  M <- do.call(rbind, mx)
  hclust_p = hclust(dist(M),method="average")
  dend_p = as.dendrogram(hclust_p)
  dend.list = lapply(hc, as.dendrogram)
  
  if(method==1) {
    merged = ComplexHeatmap::merge_dendrogram(dend_p, dend.list)
  } else {
    mrg <- hclust_p$merge    
    merged_branch <- list()
    k=1
    for(k in 1:nrow(mrg)) {
      i <- mrg[k,1]
      j <- mrg[k,2]
      if(i<0) d1 <- dend.list[[-i]]
      if(i>0) d1 <- merged_branch[[i]]      
      if(j<0) d2 <- dend.list[[-j]]
      if(j>0) d2 <- merged_branch[[j]]
      merged_branch[[k]] <- merge(d1,d2)  ## actual merge
    }
    merged <- merged_branch[[k]]
  }
  
  merged_hclust <- as.hclust(merged)
  merged_hclust
}


##=========================================================================
## CONSENSUS WGCNA
##=========================================================================

if(0) {
  power = 12
  minKME = 0.8
  cutheight = 0.15
  deepsplit = 2
  maxBlockSize = 5000
  addCombined = FALSE
}


wgcna.netReplaceColors <- function(net, newcolors) {
  oldcolors <- net$colors
  lut <- table(oldcolors, newcolors)
  old2new <- colnames(lut)[max.col(lut)]
  names(old2new) <- rownames(lut)

  ## rename colors
  ss <- names(net$colors)
  if("colors" %in% names(net)) {    
    net$colors <- old2new[net$colors]
    names(net$colors) <- ss
  }
  ## rename unmergedColors
  if("unmergedColors" %in% names(net)) {
    net$unmergedColors <- old2new[net$unmergedColors]
    names(net$unmergedColors) <- ss
  }
  ## rename MEs
  if("MEs" %in% names(net)) {
    klr <- sub("^ME","",colnames(net$MEs))
    colnames(net$MEs) <- paste0("ME",old2new[klr])
  }
  return(net)
}


#'
#'
#' @export
wgcna.runConsensusWGCNA <- function(exprList,
                                    phenoData,
                                    ngenes = 4000,
                                    power = 12,
                                    consensusPower = power,
                                    minModuleSize = 20,
                                    minKME = 0.3, 
                                    mergeCutHeight = 0.15,
                                    deepSplit = 2,
                                    maxBlockSize = 5000,
                                    addCombined = FALSE,
                                    calcMethod = "fast",
                                    drop.ref = TRUE,
                                    verbose = 1
                                    ) {

  #power=6;minKME=0.5;cutheight=0.15;deepSplit=2;maxBlockSize=5000;verbose=1;
  #calcMethod="fast"
  
  if(addCombined) {
    exprList[['Combined']] <- do.call(cbind, exprList)
  }
  names(exprList)

  exprsamples <- unlist(lapply(exprList, colnames))
  if(!all(exprsamples %in% rownames(phenoData))) {
    stop("samples mismatch for exprList and phenoData")
  }
  
  multiExpr = WGCNA::list2multiData(lapply(exprList, Matrix::t))
  cor <- WGCNA::cor  ## needed...

  if(!is.null(power) && length(power)==1) {
    power <- rep(power, length(multiExpr))
  }
  
  # run module detection procedure
  netList <- list()
  k=names(multiExpr)[1]
  for (i in 1:length(multiExpr)) {
    message("[wgcna.runConsensusWGCNA] computing WGCNA for ", k)
    k <- names(multiExpr)[i]
    netList[[k]] <- wgcna.compute(    
      X = Matrix::t(multiExpr[[i]]$data),
      samples = phenoData,
      ngenes = ngenes,
      power = power[i],
      minmodsize = minModuleSize,
      calcMethod = calcMethod,
      deepsplit = deepSplit,
      mergeCutHeight = mergeCutHeight, 
      numericlabels = FALSE,
      minKME = minKME,
      maxBlockSize = maxBlockSize,
      sv.tom = 40,
      verbose = verbose
    )
  }
  
  # now we run automatic consensus module detection
  message("[wgcna.runConsensusWGCNA] computing CONSENSUS modules...")

  if(is.null(consensusPower)) {
    consensusPower <- unlist(sapply(netList, function(w) w$net$power))
    if(!is.null(consensusPower)) {
      message("[wgcna.runConsensusWGCNA] setting to mean consensusPower")
      consensusPower <- mean(consensusPower)
    } else {
      consensusPower <- 6
    }
  }
  message("[wgcna.runConsensusWGCNA] consensusPower = ", consensusPower)
  
  sel <- setdiff(names(multiExpr), c("Combined"))
  cons = WGCNA::blockwiseConsensusModules(
    multiExpr[sel],
    power = consensusPower, 
    networkType = "signed",
    TOMType = "signed",
    minModuleSize = minModuleSize,
    deepSplit = deepSplit,
    mergeCutHeight = mergeCutHeight, 
#    detectCutHeight = 1,    
    numericLabels = FALSE,
    minKMEtoStay = minKME,
    maxBlockSize = maxBlockSize,
    saveTOMs = FALSE,
    verbose = 1
  )
  cons$power = consensusPower
  
  ## create and match colors
  old.colors <- sapply(netList, function(net) net$net$colors)
  c0 <- cons$colors
  matched.colors <- apply(old.colors, 2, function(k) WGCNA::matchLabels(k,c0))
  colors <- cbind(Consensus=c0, matched.colors)

  ## replace also colors in net solutions
  for(i in 1:length(netList)) {
    net <- netList[[i]]$net
    netList[[i]]$net <- wgcna.netReplaceColors(net, matched.colors[,i]) 
  }
  
  ## add labels to dendrogram
  for(i in 1:length(cons$dendrograms)) {
    ii <- which(cons$goodGenes & cons$blocks==i)
    xnames <- names(cons$colors)
    cons$dendrograms[[i]]$labels <- xnames[ii]
  }
  
  ## merge dendrograms
  message("[wgcna.compute] merge_block_dendrograms...")  
  multiX <- Matrix::t(do.call(rbind,lapply(exprList,function(x)scale(t(x)))))
  merged <- try(wgcna.merge_block_dendrograms(cons, multiX))
  if(!inherits(merged,"try-error")) {
    cons$merged_dendro <- merged
  } else {
    cons$merged_dendro <- NULL
  }
  
  ## create module-trait matrices for each set
  message("[wgcna.compute] computing module-traits matrices...")
  datTraits <- expandPhenoMatrix( phenoData, drop.ref=drop.ref )
  zlist <- list()
  k=1
  for(k in names(cons$multiME)) {
    M <- (cons$multiME[[k]][[1]])
    Z <- datTraits
    ##rownames(Z) <- gsub("[-]",".",rownames(Z)) ## :(
    kk <- intersect(rownames(M),rownames(Z))
    zlist[[k]] <- cor(M[kk,], Z[kk,], use="pairwise")
  }
  z.check <- sapply(zlist, function(z) colSums(z != 0, na.rm = TRUE) > 0)
  sel <- names(which(rowMeans(z.check) == 1))
  zlist <- lapply(zlist, function(z) z[, sel, drop = FALSE])

  ## create consensus module-trait matrix
  message("[wgcna.compute] computing consensus module-traits matrix...")
  consZ <- wgcna.computeConsensusMatrix(zlist, psig = 0.05) 
  
  ## add slots
  ydim <- sapply(exprList, ncol) 
  datExpr <- lapply(exprList, Matrix::t)
  
  res <- list(
    net = cons,
    netList = netList,
    datExpr = datExpr,
    datTraits = datTraits,
    modTraits = consZ,
    dendro = cons$merged_dendro,    
    colors = colors,
    zlist = zlist,
    ydim = ydim,
    class = "consensus"
  )

  res
}

#' Compute consensus matrix from list of matrices. The consensus
#' matrix checks for consistent sign and minimal threshold for each
#' matrix. Optionally filters on consistent p-value.
#'
#' 
wgcna.computeConsensusMatrix <- function(matlist, psig = 0.05) {
  ## create consensus module-trait matrix
  all.pos <- Reduce("*", lapply(matlist, function(z) (sign(z) == 1) ))
  all.neg <- Reduce("*", lapply(matlist, function(z) (sign(z) == -1) )) 
  disconcordant <- !(all.pos | all.neg)
  zsign <- sign(Reduce("+", lapply(matlist, sign)))
  consZ <- Reduce(pmin, lapply(matlist, abs)) * zsign
  consZ[disconcordant] <- NA

  if(psig < 1) {
    mdim <- sapply(exprList, ncol)
    pv <- mapply(function(z, n) corPvalueStudent(z, n), matlist, mdim, SIMPLIFY = FALSE)
    all.sig <- Reduce("*", lapply(pv, function(p) 1 * (p < psig)))
    consZ[!all.sig] <- NA
  }
  return(consZ)
}

#' Compute consensus matrix from list of matrices. The consensus
#' matrix checks for consistent sign and minimal threshold for each
#' matrix. Optionally filters on consistent p-value.
#'
#' 
wgcna.computeDistinctMatrix <- function(matlist, psig = 0.05, min.diff=0.3) {
  ## create difference module-trait matrix
  mdim <- sapply(exprList, ncol)
  pv <- mapply(function(z, n) corPvalueStudent(z, n), matlist, mdim, SIMPLIFY = FALSE)
  Q <- matlist
  for(i in 1:length(matlist)) {
    refmat <- Reduce("+", matlist[-i]) / (length(matlist)-1)
    diff <- matlist[[i]] - refmat
    notsig <- ( pv[[i]] > psig | abs(diff) < min.diff )
    Q[[i]][notsig] <- NA
  }
  return(Q)
}


wgcna.plotConsensusOverlapHeatmap <- function(net1, net2,
                                              setLabels=NULL, setpar=TRUE) {
  
  firstColors  <- wgcna.labels2colors(net1$colors)
  secondColors <- wgcna.labels2colors(net2$colors)
  overlap <- overlapTable(firstColors, secondColors)
  names(overlap)

  T1 = overlap$countTable
  T2 = table(firstColors, secondColors)
  T3 = table(net1$colors, net2$colors)

  if(is.null(setLabels)) setLabels <- c("Set1","Set2")
  if(length(setLabels)==1) setLabels <- paste0(setLabels,1:2)

  firstModTotals <- rowSums(overlap$countTable)
  secondModTotals <- colSums(overlap$countTable)
  firstModules <- rownames(overlap$countTable)
  secondModules <- colnames(overlap$countTable)

  # Truncate p values smaller than 10^{-50} to 10^{-50} 
  pTable <- -log10(overlap$pTable)
  pTable[is.infinite(pTable)] = 1.3*max(pTable[is.finite(pTable)]);
  pTable[pTable>50 ] = 50 ;
  
  if(setpar) {
    par(mfrow=c(1,1));
    par(cex = 1.0);
    par(mar=c(10, 12.4, 2.7, 1)+0.3);
  }

  # Use function labeledHeatmap to produce the color-coded table with all the trimmings
  WGCNA::labeledHeatmap(
    Matrix = t(pTable),
    xLabels = paste(" ", firstModules),
    yLabels = paste(" ", secondModules),
    colorLabels = TRUE,
    #xSymbols = paste0(setLabels[1],":", firstModules, " (", firstModTotals,")"),
    #ySymbols = paste0(setLabels[2],":", secondModules, " (", secondModTotals, ")"),
    xSymbols = paste0(firstModules, " (", firstModTotals,")"),
    ySymbols = paste0(secondModules, " (", secondModTotals, ")"),
    textMatrix = t(overlap$countTable),
    colors = WGCNA::blueWhiteRed(100)[50:100],
    main = paste("Correspondence of",setLabels[1],"and ",
      setLabels[2],"modules"),
    cex.text = 1.0, cex.lab = 1.0, setStdMargins = FALSE);   
  mtext(toupper(setLabels[1]), side=1, line=9, cex=1.6)
  mtext(toupper(setLabels[2]), side=2, line=11, cex=1.6)

}


##=========================================================================
## PLOTTING FUNCTIONS
##=========================================================================

#'
#'
#' @export
wgcna.plotTOM <- function(wgcna, justdata = FALSE, block=NULL,
                          legend=TRUE, downsample=NULL) {

  datExpr <- wgcna$datExpr
  wTOM <- NULL  
  if(!is.null(wgcna$TOM)) {
    wTOM <- wgcna$TOM
  }
  ## if SV of TOM is stored, reconstruct TOM
  if(is.null(wTOM) && !is.null(wgcna$svTOM)) {
    wTOM  <- wgcna$svTOM %*% t(wgcna$svTOM)
  }
  if(is.null(wTOM)) {
    message("[wgcna.plotTOM] ERROR. no TOM matrix")
    return(NULL)
  }
  dissTOM <- 1 - wTOM
  rownames(dissTOM) <- colnames(dissTOM) <- colnames(datExpr)
  
  ## clustering wgcnaults
  moduleColors <- NULL
  if(is.null(block) && "merged_dendro" %in% names(wgcna$net)) {
    geneTree <- wgcna$net$merged_dendro
    gg <- geneTree$labels    
    if(length(gg)>0) {
      dissTOM <- dissTOM[gg,gg]
      moduleColors <- wgcna.labels2colors(wgcna$net$colors[gg])
    }
  }

  if(is.null(moduleColors)) {
    if(is.null(block)) block <- 1
    geneTree <- wgcna$net$dendrograms[[block]]
    ii <- which(wgcna$net$blocks == block & wgcna$net$goodGenes==TRUE)
    gg <- names(wgcna$net$color)[ii]
    dissTOM <- dissTOM[gg,gg]
    moduleColors <- wgcna.labels2colors(wgcna$net$colors[gg])
  }

  if (justdata) {
    return(dissTOM)
  }
 
  if(!is.null(downsample) && ncol(dissTOM)>downsample) {
    ii <- seq(1,ncol(dissTOM),length.out=downsample)
    dissTOM <- dissTOM[ii,ii]
    moduleColors <- moduleColors[ii]
    geneTree <- fastcluster::hclust(as.dist(dissTOM),
                                    method = "average")
  }
  
  if(legend) {
    par(oma = c(1,0,0,0), mar=c(0,0,0,0))
    plotly::layout(
      matrix(c(
        0, 0, 5, 0,
        0, 0, 2, 0,
        4, 1, 3, 6
      ), nr = 3, byrow = TRUE),
      widths = c(2.3, 0.5, 10, 3),
      heights = c(2.3, 0.5, 10)
    )
    WGCNA::TOMplot(
      dissTOM^7,
      geneTree,
      moduleColors,
      setLayout = FALSE,
      main = NULL
    )
    
    ## add color legend
    frame()
    legend(
      -0.1, 1,
      fill = wgcna$me.colors,
      legend = names(wgcna$me.colors),
      cex = 1.2, bty = "n", x.intersp = 0.5
    )
  } else {
    WGCNA::TOMplot(
      dissTOM^7,
      geneTree,
      moduleColors,
      setLayout = TRUE,
      main = NULL
    )
  }
  
}

#'
#'
#' @export
wgcna.plotDendroAndColors <- function(wgcna, main=NULL, block=NULL,
                                      extra.colors=NULL,
                                      marAll = c(0.2, 5, 1, 0.2),
                                      setLayout=TRUE) {

  if("net" %in% names(wgcna)) {
    net <- wgcna$net
    if("colors" %in% names(wgcna)) {
      net$netcolors <- net$colors
      net$colors <- wgcna$colors
    }
  } else {
    net <- wgcna
  }
  
  if(length(net$dendrograms)>1) {
    message("warning: this wgcna has multiple blocks")
  }

  if(is.null(block) && "merged_dendro" %in% names(net)) {
    geneTree <- net$merged_dendro
  } else {
    if(is.null(block)) block <- 1
    geneTree <- net$dendrograms[[block]]
  }
  colors <- cbind(wgcna.labels2colors(net$colors))
  if(NCOL(colors)==1) colnames(colors)[1] <- "Module colors"

  gg <- geneTree$labels
  if(is.null(gg) && !is.null(block)) {
    ii <- which(net$blocks == block & net$goodGenes==TRUE)
    gg <- names(net$color)[ii]
  }
  if(is.null(gg) && is.null(block)) {
    ii <- which(net$goodGenes==TRUE)
    gg <- names(net$color)[ii]
  }
  colors <- colors[gg,,drop=FALSE]
  if(!is.null(extra.colors)) {
    jj <- match(gg, rownames(extra.colors))
    colors <- cbind(colors, extra.colors[jj,])
  }
  
  if(is.null(main)) main <- "Cluster dendrogram and module colors"
  ## Plot the dendrogram and the module colors underneath
  WGCNA::plotDendroAndColors(
    dendro = geneTree,
    colors = colors,
    groupLabels = colnames(colors),
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = FALSE,
    guideHang = 0.05,
    marAll = marAll,
    setLayout = setLayout,
    main = main
  )
}

wgcna.plotDendroAndTraitCorrelation <- function(wgcna,
                                                multi = FALSE,
                                                traits = NULL,
                                                main=NULL,
                                                block=NULL,
                                                rm.na = TRUE,
                                                marAll = c(0.2, 5, 1, 0.2),
                                                setLayout=TRUE)
{

  ## if consensus output do this
  is.cons <- ("class" %in% names(wgcna) && wgcna$class == "cons")
  is.cons2 <- (all(c("netList","zlist") %in% names(wgcna)))
  if(is.cons || is.cons2) {
    wgcna.plotDendroAndTraitCorrelation_cons(
      cons = wgcna,
      traits = traits,
      main = main,
      rm.na = rm.na,
      marAll = marAll,
      setLayout = setLayout)
    return()
  }
  
  X <- wgcna$datExpr
  Y <- wgcna$datTraits
  traitSig <- cor(X, Y, use="pairwise")
  if(rm.na) {
    sel <- colMeans(is.na(traitSig)) < 1
    traitSig <- traitSig[, sel, drop=FALSE]
  }
  traitColors <- rho2bluered(traitSig)

  moduleColors <- wgcna$net$colors
  moduleColors <- cbind(moduleColors)
  if(NCOL(moduleColors)==1) colnames(moduleColors) <- "Module"

  colors <- cbind( moduleColors, 0, traitColors)
  
  geneTree <- wgcna$net$dendrograms[[1]]
  colors <- colors[geneTree$labels, ]
  
  if(is.null(main)) main = "Cluster Dendrogram, Modules and Trait Correlation"

  WGCNA::plotDendroAndColors(
    geneTree,
    colors = colors,
    colnames(colors),
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05,
    marAll = marAll,
    main = main
  )

}

#' wgcna.plotDendroAndTraits for Consensus output
#'
#' 
wgcna.plotDendroAndTraitCorrelation_cons <- function(cons,
                                                     traits = NULL,
                                                     main = NULL,
                                                     rm.na = TRUE,
                                                     use.tree = 0,
                                                     marAll = c(0.2, 5, 1, 0.2),
                                                     setLayout=TRUE)
{

  moduleColors <- cons$colors

  traitSig <- list()
  nsets <- length(cons$datExpr)
  i=1
  for(i in 1:nsets) {
    Y <- cons$datTraits
    sel <- colnames(Y)
    if(!is.null(traits)) sel <- intersect(sel, traits)
    X <- cons$datExpr[[i]]
    kk <- intersect(rownames(X), rownames(Y))
    traitSig[[i]] <- cor(X[kk,], Y[kk,sel], use="pairwise")
  }

  if(rm.na) {
    for(i in 1:length(traitSig)) {
      sel <- colMeans(is.na(traitSig[[i]])) < 1
      traitSig[[i]] <- traitSig[[i]][, sel, drop=FALSE]
    }
  }
  
  sets <- names(cons$datExpr)
  ## prepend datatype/set name
  for(i in 1:length(traitSig)) {
    colnames(traitSig[[i]]) <- paste0(sets[i],":",colnames(traitSig[[i]]))
  }
  traitSig2 <- c()
  for(i in 1:length(traitSig)) {
    traitSig2 <- cbind(traitSig2, traitSig[[i]])
    if(i < length(traitSig)) traitSig2 <- cbind(traitSig2, 0)
  }
  traitColors <- rho2bluered(traitSig2, f=0.95)
  ii <- which(colnames(traitColors)=='')
  if(length(ii)) traitColors[,ii] <- "#FFFFFF"
  
  colors <- cbind( moduleColors, 0, traitColors)
  if(use.tree == 0) {
    geneTree <- cons$net$dendrograms[[1]]
  } else {
    geneTree <- cons$netList[[use.tree]]$net$dendrograms[[1]]
  }
  
  if(is.null(main)) main = "Cluster Dendrogram, Modules and Trait Correlation"
  
  WGCNA::plotDendroAndColors(
    geneTree,
    colors = colors,
    colnames(colors),
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05,
    marAll = marAll,
    main = main
  )

}

#' Converts correlation values [-1;1] to blue-white-red colors. Good
#' for creating color labels for labeledHeatmaps that expect colors.
#' 
rho2bluered <- function(R, a=1, f=0.95) {
  BLUERED <- WGCNA::blueWhiteRed(100)
  if(a!=1) R <- sign(R) * abs(R)**a
  if(is.null(ncol(R))) {
    col <- BLUERED[1+round(99*(1+R)/2)]
  } else {
    col <- apply(R, 2, function(x) BLUERED[1+round(99*(1+x)/2)])
    dimnames(col) <- dimnames(R)
  }
  if(f < 1) {
    col <- apply(col, 2, adjustcolor, red.f=f, green.f=f, blue.f=f)
  }
  col
}


#' Converts WGCNA labels (numeric or color) to colors.
#'
wgcna.labels2colors <- function(colors, ...) {
  if (all(is.numeric(colors))) {
    colors <- WGCNA::labels2colors(colors, ...)
    return(colors)
  }
  stdColors <- c("grey", WGCNA::standardColors())
  if (all(colors %in% stdColors)) {
    return(colors)
  }
  icolors <- as.integer(factor(as.character(colors)))
  colors <- WGCNA::standardColors()[icolors]
  return(colors)
}



#' Plot membership correlation vs gene signficance (correlation with
#' trait) to discover biomarkers/driver genes.
#'
#' @export
wgcna.plotMMvsGS <- function(wgcna, module, trait, abs = TRUE, par = TRUE,
                             plotlib = "base") {
  ## module="ME3";trait="activated=act"
  moduleGenes <- wgcna$me.genes[[module]]
  nSamples <- nrow(wgcna$datExpr)

  ## Module membership correlation (with p-values)
  if ("stats" %in% names(wgcna)) {
    moduleMembership <- wgcna$stats$moduleMembership
    MMPvalue <- wgcna$stats$MMPvalue
  } else {
    moduleMembership <- as.data.frame(cor(wgcna$datExpr, wgcna$net$MEs, use = "p"))
    MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(moduleMembership), nSamples))
  }

  ## Gene-trait significance (trait correlation) (with p-values)
  if ("stats" %in% names(wgcna)) {
    traitSignificance <- wgcna$stats$traitSignificance
    GSPvalue <- wgcna$stats$GSPvalue
  } else {
    traitSignificance <- as.data.frame(cor(wgcna$datExpr, wgcna$datTraits, use = "p"))
    GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(traitSignificance), nSamples))
  }

  x <- (moduleMembership[moduleGenes, module])
  y <- (traitSignificance[moduleGenes, trait])
  if (abs == TRUE) {
    x <- abs(x)
    y <- abs(y)
  }
  ##
  px <- MMPvalue[moduleGenes, module]
  py <- GSPvalue[moduleGenes, trait]
  qx <- p.adjust(px, method = "fdr")
  qy <- p.adjust(py, method = "fdr")
  is.sig <- (qx < 0.05 & qy < 0.05)
  sigx <- (qx < 0.05)
  sigy <- (qy < 0.05)
  ii <- which(is.sig)
  qv <- quantile(x[ii], prob = 0.1)[1]
  qh <- quantile(y[ii], prob = 0.1)[1]

  pos <- cbind(x, y)
  rownames(pos) <- moduleGenes
  is.sig1 <- c("notsig", "onesig", "sig")[1 + 1 * sigx + 1 * sigy]
  hi1 <- NULL
  ## hi1 <- head(rownames(pos),10)
  col1 <- c("grey70", "grey20", "red2")

  if (par) par(mfrow = c(1, 1), mar = c(5, 5, 3, 2))
  if (plotlib == "ggplot") {
    pgx.scatterPlotXY.GGPLOT(
      pos,
      var = is.sig1, hilight = hi1, col = col1,
      xlab = paste("Module membership in", module, "module"),
      ylab = paste("Gene significance for trait", trait),
      title = paste("Module membership vs. gene significance\n"),
      cex.title = 0.9,
      girafe = FALSE
    )    
  } else if (plotlib == "girafe") {
    pgx.scatterPlotXY.GGPLOT(
      pos,
      var = is.sig1, hilight = hi1, col = col1,
      xlab = paste("Module membership in", module, "module"),
      ylab = paste("Gene significance for trait", trait),
      title = paste("Module membership vs. gene significance\n"),
      cex.title = 0.7, cex.axis = 0.7,
      girafe = TRUE
    )
  } else {
    ii <- which(is.sig1 == "notsig")
    verboseScatterplot(
      x[-ii], y[-ii],
      xlab = paste("Module membership in", module, "module"),
      ylab = paste("Gene significance for trait", trait),
      main = paste("Module membership vs. gene significance\n"),
      cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = col1[1]
    )
    ii <- which(is.sig1 == "onesig")
    points(x[ii], y[ii], col = col1[2])
    ii <- which(is.sig1 == "sig")
    points(x[ii], y[ii], col = col1[3])
    abline(v = qv, h = qh, col = "darkred")
  }
}

#' @export
wgcna.plotModuleTraitHeatmap <- function(wgcna, setpar = TRUE, cluster = FALSE,
                                         multi = FALSE, main = NULL, justdata = FALSE,
                                         transpose = FALSE, colorlabel = TRUE,
                                         nmax = -1, text = TRUE, pstar = TRUE) {
  
  if(!multi) wgcna <- list(wgcna)
  nSamples <- nrow(wgcna[[1]]$datExpr)
  MEs <- do.call(cbind, lapply(wgcna, function(w) as.matrix(w$net$MEs))) 
  Y <- wgcna[[1]]$datTraits
  
  moduleTraitCor <- cor(MEs, Y, use = "pairwise.complete")
  if(nmax > 0) {
    sel <- head(order(-apply(abs(moduleTraitCor), 1, max, na.rm=TRUE)),nmax)
    moduleTraitCor <- moduleTraitCor[sel,,drop=FALSE]
  }

  if(transpose) {
    moduleTraitCor <- t(moduleTraitCor)
  }

  wgcna.plotLabeledCorrelationHeatmap(
    R = moduleTraitCor,
    nSamples = nSamples,
    setpar = setpar,
    cluster = cluster,
    text = text,
    main = main,
    justdata = justdata,
    colorlabel =  colorlabel,
    pstar = pstar
  ) 
}


#' Plot cluster dendrogram with eigengenes and traits.
#'
#' @export
wgcna.plotEigenGeneClusterDendrogram <- function(wgcna,
                                                 add_traits = TRUE,
                                                 horiz = FALSE,
                                                 setMargins = TRUE,
                                                 method = 'wgcna',
                                                 showlabels = TRUE,
                                                 plot = TRUE,
                                                 multi = FALSE,
                                                 main = NULL) {
  # Matrix with eigengenes and traits
  if(multi) {
    ME <- lapply(wgcna, function(w) as.matrix(w$net$MEs))
    ME <- do.call(cbind, ME)
    Y <- wgcna[[1]]$datTraits
  } else {
    ME <- wgcna$net$MEs
    Y <- wgcna$datTraits
  }

  if (length(add_traits)==1 && is.logical(add_traits) && add_traits==TRUE) {
    ME <- cbind(ME, Y)
  } else if (length(add_traits) > 0 && !is.logical(add_traits)) {
    sel <- intersect(add_traits, colnames(Y))
    if(length(sel)) ME <- cbind(ME, Y[,sel])
  }

  ME <- WGCNA::orderMEs(ME)
  if (NCOL(ME) <= 2) ME <- cbind(ME, ME) ## error if ncol(ME)<=2 !!!!
  if (is.null(main)) main <- "Eigengene cluster dendrogram"

  hc <- NULL
  if(method == 'wgcna') {
    ## plot dendrogram with WGCNA function
    WGCNA::plotEigengeneNetworks(
      ME, main,
      setMargins = setMargins,
      marDendro = c(0, 4, 2, 0),
      plotHeatmaps = FALSE
    )

  } else {
    ## plot dendrogram with hclust function
    if(setMargins && horiz) par(mar=c(4,4,4,8))
    if(setMargins && !horiz) par(mar=c(8,4,4,1))    
    hc <- hclust( as.dist(1 - cor(ME)), method="average")
    if(plot) {
      save.labels <- hc$labels
      if(!showlabels) hc$labels <- rep("",ncol(ME))
      plot( as.dendrogram(hc), horiz = horiz)
      hc$labels <- save.labels
    }
  }

  invisible(hc)
}

#' Plot the adjacency correlation heatmap matrix of eigengenes with or
#' without traits. This can show how traits cluster together with the
#' eigengenes.
#'
#' @export
wgcna.plotEigenGeneAdjacencyHeatmap <- function(wgcna, add_traits = TRUE,
                                                marx = 1, main = NULL,
                                                multi = FALSE, 
                                                colorlabel = TRUE,
                                                text = FALSE,
                                                pstar = TRUE,
                                                setMargins = TRUE,
                                                mar1 = c(5.5, 5, 1.6, 1),
                                                mar2 = c(8, 10, 4, 2),
                                                plotDendro = TRUE,
                                                plotHeatmap = TRUE, 
                                                dendro.horiz = TRUE,
                                                dendro.width = 0.3,
                                                dendro.labels = TRUE,
                                                method = 2,
                                                nmax = -1,
                                                mask.intra = FALSE,
                                                justdata = FALSE) {

  if(!multi) wgcna <- list(wgcna)

  # Matrix with eigengenes and traits
  ME <- lapply(wgcna, function(w) as.matrix(w$net$MEs))
  ME <- do.call(cbind, ME)
  Y <- wgcna[[1]]$datTraits    

  if (length(add_traits)==1 && is.logical(add_traits) && add_traits==TRUE) {
    ME <- cbind(ME, Y)
  } else if (length(add_traits) > 0 && !is.logical(add_traits)) {
    sel <- intersect(add_traits, colnames(Y))
    if(length(sel)) ME <- cbind(ME, Y[,sel])
  }
  
  if (NCOL(ME) <= 2) ME <- cbind(ME, ME) ## error if ncol(ME)<=2 !!!!

  ## eigengene correlation
  R <- cor(ME, use="pairwise")

  if(nmax>0) {
    ii <- head( order(-Matrix::rowMeans(abs(R))), nmax)
    R <- R[ii,ii]
  }
  
  if (justdata) {
    return(R)
  }
  
  # Plot the correlation heatmap matrix (note: this plot will overwrite
  # the dendrogram plot)
  if (is.null(main)) main <- "Eigengene adjacency heatmap"

##  plotDendro=TRUE;plotHeatmap=TRUE;setMargins=TRUE;marx=1
  
  if(method == 1) {
    ME <- as.data.frame(ME)  ## argh.... WGCNA like dataframes...
    ME <- WGCNA::orderMEs(ME)
    names(ME) <- paste0(" ",names(ME))    
    WGCNA::plotEigengeneNetworks(
      ME,
      main,
      marHeatmap = c(8 * marx, 10 * marx, 1, 2),
      marDendro = c(0, 8 * marx, 2, 2),
      plotDendrograms = plotDendro,
      plotHeatmaps = plotHeatmap,
      #excludeGrey = FALSE,
      #greyLabel = "grey", 
      colorLabels = FALSE,
      setMargins = setMargins,
      xLabelsAngle = 45,
      yLabels = colnames(ME),
      xLabels = colnames(ME),
      xSymbols = colnames(ME),
      ySymbols = colnames(ME)
    )
  }

  if(method == 2) {

    if(plotDendro && plotHeatmap) {
      layout.matrix <- matrix( 1:2, nrow = 1, ncol = 2)
      layout(mat = layout.matrix, heights = 1, widths = c(dendro.width, 1))
      if(dendro.horiz && dendro.labels) {
        mar1[4] <- mar2[2] ## copy left margin
      }
    }
    if(plotDendro) par(mar=mar1)

    hc <- wgcna.plotEigenGeneClusterDendrogram(
      wgcna,
      multi = TRUE,
      add_traits = add_traits,
      method = "hclust",
      showlabels = dendro.labels,
      horiz = dendro.horiz,
      plot = plotDendro,
      setMargins = FALSE
    )

    if(plotHeatmap) {
      ii <- hc$labels[hc$order]
      ii <- intersect(ii, rownames(R))
      R1 <- R[rev(ii), ii]
      ##image(R1)
      nsamples <- nrow(Y)
      par(mar=mar2)
      wgcna.plotLabeledCorrelationHeatmap(
        R1,
        #is.dist = TRUE,
        nSamples = nsamples,
        text = text,
        pstar = pstar,
        cluster = FALSE,
        setpar = FALSE,
        main = "Eigengene correlation heatmap"
      )
    }
    
  }


  
}

#' @export
wgcna.plotEigenGeneGraph <- function(wgcna, add_traits = TRUE, main = NULL,
                                     multi=FALSE, vcex=1, labcex=1) {

  ## require(igraph)
  if(multi) {
    ME <- lapply(wgcna, function(w) as.matrix(w$net$MEs))
    ME <- do.call(cbind, ME)
    if (add_traits)  ME <- cbind(ME, wgcna[[1]]$datTraits)    
  } else {
    ME <- wgcna$net$MEs
    if (add_traits)  ME <- cbind(ME, wgcna$datTraits)
  }
  if (NCOL(ME) <= 2) ME <- cbind(ME, ME) ## error if ncol(ME)<=2 !!!!

  sdx <- matrixStats::colSds(as.matrix(ME), na.rm = TRUE)
  if (any(sdx == 0)) ME <- ME + runif(length(ME), 0, 1e-5)
  
  ## Recalculate MEs with color as labels
  clust <- hclust(dist(t(scale(ME))))
  phylo <- ape::as.phylo(clust)
  gr <- igraph::as.igraph(phylo, directed = FALSE)

  is.node <- grepl("Node", igraph::V(gr)$name)
  module.name <- igraph::V(gr)$name
  if(multi) {
    module.size <- lapply(wgcna, function(w) table(w$net$labels))
    names(module.size) <- NULL
    module.size <- unlist(module.size)
    module.colors <- sapply(wgcna, function(w) w$me.colors)
    names(module.colors) <- NULL
    module.colors <- unlist(module.colors)        
  } else {
    module.size <- table(wgcna$net$labels)
    module.colors <- wgcna$me.colors
  }
  module.size <- module.size / mean(module.size)
  igraph::V(gr)$label <- igraph::V(gr)$name
  igraph::V(gr)$label[is.node] <- NA
  igraph::V(gr)$color <- module.colors[module.name]
  igraph::V(gr)$size <- vcex * 18 * (module.size[module.name])**0.4
  igraph::V(gr)$size[is.na(igraph::V(gr)$size)] <- 0

  ## par(mfrow = c(1, 1), mar = c(1, 1, 1, 1) * 0)
  igraph::plot.igraph(
    gr,
    layout = igraph::layout.kamada.kawai,
    vertex.label.cex = 0.85 * labcex,
    edge.width = 3
  )
  if (!is.null(main)) title(main, line = -1.5)
}


#' Plot Multi-dimensional scaling (MDS) of centered data matrix.
#'
#' @export
wgcna.plotMDS <- function(wgcna, main = NULL, scale = FALSE) {
  cc <- wgcna.labels2colors(wgcna$net$color)
  pc <- svd(t(scale(wgcna$datExpr, scale = scale)), nv = 2)$u[, 1:2]
  # pc <- svd(t(scale(wgcna$datExpr)),nv=1)$u[,1:2]
  colnames(pc) <- c("MDS-x", "MDS-y")
  if (is.null(main)) main <- "MDS of features"
  plot(pc, col = cc, main = main)
}

#'
#'
#' @export
wgcna.plotFeatureUMAP <- function(wgcna, nhub = 3, method = "clust",
                                  scale = FALSE, main = NULL,
                                  plotlib = "base") {
  if (method == "clust" && "clust" %in% names(wgcna)) {
    pos <- wgcna$clust[["umap2d"]]
  } else if (method == "umap") {
    cX <- t(scale(wgcna$datExpr, scale = scale)) ## WGCNA uses correlation
    pos <- uwot::umap(cX)
    colnames(pos) <- c("UMAP-x", "UMAP-y")
    rownames(pos) <- colnames(wgcna$datExpr)
    ##  } else if(method=="mds") {
  } else {
    pos <- svd(t(scale(wgcna$datExpr, scale = scale)), nv = 2)$u[, 1:2]
    colnames(pos) <- c("MDS-x", "MDS-y")
    rownames(pos) <- colnames(wgcna$datExpr)
  }

  if (is.null(main)) main <- "Feature UMAP colored by module"

  ## get top hub genes
  mm <- wgcna$stats$moduleMembership
  hubgenes <- apply(mm, 2, function(x) head(names(sort(-x)), 3), simplify = FALSE)
  hubgenes
  sel <- which(names(hubgenes) != "MEgrey")
  hubgenes <- unlist(hubgenes[sel])
  col1 <- wgcna$net$colors
  genes1 <- names(which(col1 != "grey"))
  pgx.scatterPlotXY(
    pos,
    var = col1, col = sort(unique(col1)),
    hilight = genes1, hilight2 = hubgenes,
    cex.lab = 1.2, label.clusters = FALSE,
    title = main,
    plotlib = plotlib
  )
}


#' Plot module significance.
#'
#' @export
wgcna.plotModuleSignificance <- function(wgcna, trait, main = NULL, abs = FALSE) {
  ## cc <- paste0("ME",wgcna$net$color)
  cc <- wgcna.labels2colors(wgcna$net$color)
  if ("stats" %in% names(wgcna)) {
    traitSignificance <- wgcna$stats$traitSignificance
  } else {
    traitSignificance <- as.data.frame(cor(wgcna$datExpr, wgcna$datTraits, use = "p"))
    names(traitSignificance) <- names(wgcna$datTraits)
    rownames(traitSignificance) <- colnames(wgcna$datExpr)
  }
  geneSig <- traitSignificance[, trait]
  if (is.null(main)) main <- paste("Module significance with", trait)
  if (abs) geneSig <- abs(geneSig)
  WGCNA::plotModuleSignificance(
    geneSig,
    colors = cc, main = main, boxplot = FALSE
  )
}

#'
#'
#' @export
wgcna.plotSampleDendroAndColors <- function(wgcna, input.type="wgcna",
                                            what = c("me", "traits", "both")[3],
                                            datTraits = NULL, datExpr = NULL,
                                            clust.expr = TRUE, setLayout = TRUE,
                                            marAll = c(0.2, 7, 1.5, 0.5),
                                            main = NULL, justdata = FALSE) {
  
  if(input.type == 'net') {
    ME0 <- wgcna$MEs
    if(is.null(datExpr)) stop("must supply datExpr")
    if(is.null(datTraits)) stop("must supply datTraits")    
  } else {
    ME0 <- wgcna$net$MEs
    datTraits <- wgcna$datTraits
    datExpr <- wgcna$datExpr
  }

  ME <- ME0[, 0]
  samples <- rownames(ME)
  if (any(what %in% c("me", "both"))) {
    ME <- cbind(ME, ME0)
  }
  if (any(what %in% c("traits", "both"))) {
    ME <- cbind(ME, datTraits[samples,])
  }
  if (NCOL(ME) <= 2) ME <- cbind(ME, ME) ## error if ncol(ME)<=2 !!!!

  sdx <- matrixStats::colSds(as.matrix(ME), na.rm = TRUE)
  ME <- ME[, which(sdx>0), drop=FALSE]
  
  ## Recalculate MEs with color as labels
  if (clust.expr) {
    sampleTree <- hclust(as.dist(1 - cor(t(datExpr))), method = "average")
  } else {
    sampleTree <- hclust(as.dist(1 - cor(t(ME0))), method = "average")
  }
  ii <- sampleTree$order
  jj <- hclust(dist(t(scale(ME))))$order
  colors <- WGCNA::numbers2colors(ME[, jj])

  if (justdata) {
    return(ME)
  }

  if (is.null(main)) {
    if (what == "me") main <- "Sample dendrogram and module heatmap"
    if (what == "traits") main <- "Sample dendrogram and trait heatmap"
    if (what == "both") main <- "Sample dendrogram and module+traits heatmap"
  }

  ## Plot the dendrogram and the module colors underneath
  WGCNA::plotDendroAndColors(
    dendro = sampleTree,
    colors = colors,
    groupLabels = colnames(ME)[jj],
    dendroLabels = rownames(ME),
    hang = 0.03,
    addGuide = FALSE,
    guideHang = 0.05,
    setLayout = setLayout,
    marAll = marAll,
    main = main
  )
}


#' @export
wgcna.plotLabeledCorrelationHeatmap <- function(R, nSamples, setpar = TRUE,
                                                cluster = TRUE, text = TRUE,
                                                main = NULL, justdata = FALSE,
                                                colorlabel = TRUE, pstar = TRUE,
                                                cex.text = 0.7, cex.lab = NULL,
                                                is.dist = FALSE) {

  ## Define numbers of genes and samples
  if (cluster && nrow(R)>1 && ncol(R)>1) {
    R0 <- R
    R0[is.na(R0)] <- 0
    if(is.dist) {
      ii <- hclust(as.dist(abs(R0)))$order
      jj <- ii
    } else {
      ii <- hclust(dist(R0), method="average")$order
      jj <- hclust(dist(t(R0)), method="average")$order
    }
    R <- R[ii, jj]
  } 
  Pvalue <- corPvalueStudent(R, nSamples)

  if (justdata) {
    return(R)
  }

  ## Will display correlations and their p-values
  if (pstar) {
    textPv <- cut(Pvalue,
      breaks = c(-1, 0.001, 0.01, 0.05, 99),
      #labels = c("★★★", "★★", "★", "")
      #labels = c("***", "**", "*", "")
      labels = c("+++", "++", "+", "")
    )
  } else {
    textPv <- paste0("(", signif(Pvalue, 1), ")")
  }

  textMatrix <- NULL
  if(text && pstar)  textMatrix <- paste0(signif(R, 2), "\n", textPv)
  if(text && !pstar)  textMatrix <- paste0(signif(R, 2))
  if(!text && pstar)  textMatrix <- textPv
  if(!text && !pstar) textMatrix <- NULL
  
  if(!is.null(textMatrix)) {
    textMatrix[which(is.na(R))] <- ""
    dim(textMatrix) <- dim(R)
  }

  if(!colorlabel) {
    colnames(R) <- paste0(" ", colnames(R))
    rownames(R) <- paste0(" ", rownames(R))
  }
  
  if (setpar) par(mar = c(8, 8, 3, 3))
  if (is.null(main)) main <- "Correlation heatmap"
  ## Display the correlation values within a heatmap plot
  WGCNA::labeledHeatmap(
    Matrix = R,
    xLabels = colnames(R),
    yLabels = rownames(R),
    xSymbols = colnames(R),
    ySymbols = rownames(R),
    colorLabels = FALSE,
    colors = WGCNA::blueWhiteRed(50),
    textMatrix = textMatrix,
    setStdMargins = FALSE,
    cex.text = cex.text,
    cex.lab = cex.lab,
    zlim = c(-1, 1),
    main = main
  )
}


#' Plot module hub genes as graph in circle layout. This shows the
#' main connectivity structure for the modules. It can be used for
#' (edge) preservation analysis for comparing different group of
#' samples..
#'
#' @export
wgcna.plotModuleHubGenes <- function(wgcna, modules = NULL,
                                     alpha = 0.5, setpar = TRUE) {
  if (is.null(modules)) {
    modules <- colnames(wgcna$stats$moduleMembership)
  }
  modules <- intersect(modules, colnames(wgcna$stats$moduleMembership))
  if (length(modules) == 0) {
    message("ERROR. no valid modules")
    return(NULL)
  }

  if (setpar) {
    nr <- floor(length(modules)**0.5)
    nc <- ceiling(length(modules) / nr)
    nr
    nc
    par(mfrow = c(nr, nc), mar = c(0, 1, 2.5, 1))
  }
  for (k in modules) {
    mm <- wgcna$stats$moduleMembership[, k]
    mm.score <- head(sort(mm, decreasing = TRUE), 30)
    topgenes <- names(mm.score)
    A <- cor(wgcna$datExpr[, topgenes])
    diag(A) <- NA
    A <- (A - min(A, na.rm = TRUE)) / (max(A, na.rm = TRUE) - min(A, na.rm = TRUE))
    gr <- igraph::graph_from_adjacency_matrix(
      A,
      mode = "undirected", weighted = TRUE, diag = FALSE
    )
    norm.mm.score <- (mm.score - min(mm.score)) / (max(mm.score) - min(mm.score))
    clr <- sub("ME", "", k)
    if (!is.na(as.integer(clr))) clr <- as.integer(clr)
    if (clr == "black") clr <- "grey40"
    plot(gr,
      layout = igraph::layout_in_circle,
      edge.width = 6 * igraph::E(gr)$weight**8,
      vertex.size = 5 + 15 * norm.mm.score,
      vertex.color = adjustcolor(clr, alpha.f = alpha),
      vertex.frame.color = clr
    )
    title(paste("Module", k), line = 0.33)
  }
}

#' @title Map module colors to rainbow palette
#'
#' @description
#' Maps the module colors from a WGCNA network to a rainbow palette.
#'
#' @param net A WGCNA network object.
#'
#' @details
#' This function takes a WGCNA network object and maps the module colors
#' to a rainbow palette based on the modules' hierarchical clustering order.
#'
#' It extracts the hierarchical clustering dendrogram order, gets the number
#' of unique module colors, ranks the module color means based on the
#' dendrogram order, and assigns rainbow colors accordingly.
#'
#' This allows easier visualization and interpretation of modules in the
#' standard rainbow palette order.
#'
#' @return
#' A named character vector mapping the original module colors to rainbow colors.
#'
labels2rainbow.BAK <- function(net) {
  hc <- net$dendrograms[[1]]
  nc <- length(unique(net$colors))
  n <- length(net$colors)
  ii <- rep(NA, n)
  ii[net$goodGenes] <- hc$order
  col1 <- WGCNA::labels2colors(net$colors)
  col.rnk <- rank(tapply(1:n, col1[ii], mean))
  new.col <- grDevices::rainbow(nc)[col.rnk]
  names(new.col) <- names(col.rnk)
  new.col["grey"] <- "#AAAAAA"
  new.col <- new.col[col1]
  names(new.col) <- net$colors
  return(new.col)
}


#' Filter color vector by minimum KME and mergeCutHeight. Set color of
#' features with KME smaller than minKME to grey (or 0) group. Merge
#' similar modules with (module) correlation larger than
#' (1-mergeCutHeight) together.
#'
#' @export
wgcna.filterColors <- function(X, colors, minKME=0.3, mergeCutHeight=0.15,
                               minmodsize = 20, ntop=-1 ) {
  ##minKME=0.3;mergeCutHeight=0.15;minmodsize=20;ntop=-1

  sX <- X + 1e-8*matrix(rnorm(length(X)),nrow(X),ncol(X))
  sX <- t(scale(t(sX)))

  ## get singular vectors and correct sign
  vv <- tapply(1:nrow(sX), colors, function(i) svd(sX[i,],nv=1)$v[,1])
  mm <- tapply(1:nrow(sX), colors, function(i) colMeans(sX[i,]))
  vv.sign <- mapply( function(a,b) sign(cor(a,b)), mm, vv)
  vv <- mapply( function(a,b) a*b, vv, vv.sign, SIMPLIFY=FALSE)

  kme <- rep(NA,nrow(X))
  names(kme) <- rownames(X)
  names(colors) <- rownames(X)

  grey.val <- NULL
  is.color <- mean(colors %in% WGCNA::standardColors(435)) > 0.8  
  if(is.numeric(colors)) {
    colors <- as.integer(colors)
    grey.val <- 0
  } else {
    colors <- as.character(colors)
    grey.val <- "---"
    if(is.color) grey.val <- "grey"
  }
  names(colors) <- rownames(X)  
  new.colors <- colors

  if(minKME > 0) {
    i=1
    for(i in 1:length(vv)) {
      ii <- which(colors == names(vv)[i])
      r <- cor(t(X[ii,]), vv[[i]])[,1]
      max(r)
      jj <- ii[which(r < minKME)]
      if(length(jj)) {
        new.colors[jj] <- NA
      }
      kme[ii] <- r
    }  
    new.colors[is.na(new.colors)] <- grey.val
  }

  ## merge groups
  if(mergeCutHeight > 0) {
    mx <- rowmean(X, new.colors)
    rr <- cor(t(mx))
    diag(rr) <- 0
    merge.idx <- which(rr > (1 - mergeCutHeight), arr.ind=TRUE)
    if(nrow(merge.idx)>0) {
      i=1
      for(i in 1:nrow(merge.idx)) {
        aa <- rownames(rr)[merge.idx[i,]]
        jj <- which(new.colors  %in% aa)
        max.color <- names(which.max(table(new.colors[jj])))
        new.colors[jj] <- max.color
      }
    }
  }

  ## remove small groups
  modsize <- table(new.colors)
  modsize
  if( min(modsize) < minmodsize ) {
    small.mod <- names(which(modsize < minmodsize))
    sel <- which( new.colors %in% small.mod)
    new.colors[sel] <- NA
  }

  ## Filter by KME score
  if(ntop>0) {
    keep <- tapply( names(kme), new.colors, function(i) head(names(sort(-kme[i])),ntop) )
    keep <- unlist(keep)
    not.keep <- setdiff(names(kme), keep)
    if(length(not.keep)) new.colors[not.keep] <- NA
  }

  new.colors[which(is.na(new.colors))] <- grey.val
  ##if(!is.numeric(colors)) new.colors <- factor(new.colors)

  return(new.colors)
}

#' Wrapper to hclust from matrix using default WGCNA parameters.
#' 
wgcna.tomclust <- function(X, power=6) {
  A  <- WGCNA::adjacency(t(X), power=power, type = "signed")
  TOM <- fastTOMsimilarity(A, tomtype="signed", lowrank=40)
  hc <- hclust(as.dist(1-TOM), method="average")    
  hc
}


wgcna.checkDendroHeights <- function(datExpr, n=200, powers=NULL, maxpower=20) {
  ii <- 1:ncol(datExpr)
  if(n < ncol(datExpr)) {
    ##ii <- sample(1:ncol(datExpr), n)
    ii <- head(order(-matrixStats::colSds(datExpr)),n)
  }
  tX <- datExpr[,ii]
  ht <- list()
  p=9
  p=24
  if(is.null(powers)) {
    powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
    if(maxpower>20) {
      powers <- c(powers, seq(from = 20, to = maxpower, by = 5))
    }
  }
  for(i in 1:length(powers)) {
    A  <- WGCNA::adjacency(tX, power=powers[i], type = "signed")
    TOM <- fastTOMsimilarity(A, tomtype="signed", lowrank=40)
    hc <- hclust(as.dist(1-TOM), method="average")    
    ht[[i]] <- hc$height
  }
  names(ht) <- paste0("p=",powers)
  S <- sapply(ht, quantile, probs=c(0.25,0.5,0.75))
  iqr <- S[3,] - S[1,]
  optK <- powers[which.max(iqr)]

  list(
    quantiles = S,
    IQR = iqr,
    optK = optK
  )
}

#'
#'
#' @export
wgcna.plotPowerAnalysis <- function(datExpr, networktype = "signed",
                                    cex=1, maxpower = 20, main=NULL,
                                    plots=c("sft.modelfit", "mean.k",
                                      "dendro.IQR"),
                                    RsquaredCut = 0.85, setPar=TRUE) {

  RsquaredCut <- RsquaredCut[1]
  
  ## Choose a set of soft-thresholding powers
  powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
  if(maxpower>20) {
    powers <- c(powers, seq(from = 20, to = maxpower, by = 5))
  }

  ## Call the network topology analysis function
  sft <- WGCNA::pickSoftThreshold(
    datExpr,
    powerVector = powers,
    RsquaredCut = RsquaredCut, 
    networkType = networktype,
    verbose = 0
  )

  ## This is more robust
  optPower <- wgcna.pickSoftThreshold(datExpr=NULL, sft=sft,
    rcut=RsquaredCut, plot = FALSE) 
  
  if(setPar) {
    np <- length(plots)
    nc <- ceiling(sqrt(np))
    par(mfrow = c(nc, nc), mar = c(3.3, 3.5, 1, 1), mgp = c(2, 0.9, 0))
    par(mfrow = c(1, np), mar = c(3.8, 3.8, 1, 1), mgp = c(2.4, 0.95, 0))
  }
  
  ## Plot the results:
  if("sft.modelfit" %in% plots) {
    ## Scale-free topology fit index as a function of the soft-thresholding power
    y <- -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]
    base::plot(
      x = sft$fitIndices[, 1],
      y = y,
      ylim = c(min(y),1),
      type = "n",
      xlab = "Soft threshold (power)",
      ylab = "SFT model fit (signed R^2)",
      main = main
    )
    abline(h = 0, col = "black", lty=3)
    text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
      labels = powers, cex = cex, col = "red"
    )
    ## this line corresponds to using an R^2 cut-off of h
    abline(h = RsquaredCut, col = "red", lty=2)
    legend("bottomright", legend=paste("opt. power =",optPower))
  }
  
  ## Mean connectivity as a function of the soft-thresholding power
  if("mean.k" %in% plots) {
    base::plot(sft$fitIndices[, "Power"], sft$fitIndices[, "mean.k."],
      type = "n",
      xlab = "Soft threshold (power)",
      ylab = "Mean connectivity",
      main = main
    )
    text(sft$fitIndices[,"Power"], sft$fitIndices[, "mean.k."], labels = powers,
      cex = cex, col = "red")
  }

  ht <- NULL
  if("dendro.IQR" %in% plots) {
    ht <- wgcna.checkDendroHeights(datExpr, n=100, powers=powers)
    base::plot(
      sft$fitIndices[, 1], ht$IQR,
      type = "n",
      xlab = "Soft threshold (power)",
      ylab = "Dendrogram height IQR",
      main = main
    )
    text(sft$fitIndices[, 1], ht$IQR, labels = powers,
      cex = cex, col = "red")
  }
}

#' Better (?) method to pick soft threshold (aka power). 
#'
#' 
wgcna.pickSoftThreshold <- function(datExpr, sft=NULL, rcut=0.85, plot=TRUE,
                                    powers = NULL, verbose = 1) {

  if(is.null(powers)) {
    powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
    powers <- c(powers, seq(from = 25, to = 50, by = 5))
  }
  
  if(is.null(sft)) {
    sft <- WGCNA::pickSoftThreshold(
      datExpr,
      powerVector = powers,
      networkType = "signed",
      verbose = verbose
    )
  }

  sqr <- -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]
  sqr0 <- sqr
  pwr <- sft$fitIndices[, 1]
  ## remove initial value that are possible negative
  if(sqr[1] < 0.05) {
    for(i in 1:length(sqr)) sqr[i] <- ifelse(sqr[i] < 0.05, NA, sqr[i])
  }
  ds <- 0.5*median(abs(diff(sqr)),na.rm=TRUE) ##small step
  if(any(diff(sqr) < -ds, na.rm=TRUE)) {
    i <- min(which(diff(sqr) < -ds))+1
    sqr[i:length(sqr)] <- NA
  }
  if(plot) {
    plot(pwr, sqr0, ylim=c(min(sqr0),1), cex=1.5)
    abline(h=0, lty=3)
    points(pwr, sqr, pch=19, cex=1.5)
  }
  
  if( max(sqr, na.rm=TRUE) >= rcut ) argmax <- min(which(sqr >= rcut))
  if( max(sqr, na.rm=TRUE) < rcut ) argmax <- which.max(sqr)
  optPower <- sft$fitIndices[argmax, 1]
  optPower  

  if(verbose>0) {
    message("[wgcna.pickSoftThreshold] sft$powerEstimate = ", sft$powerEstimate)
    message("[wgcna.pickSoftThreshold] optPower = ", optPower)
  }
  optPower
}


#' Scale a list of TOM matrices so that the quantiles (default p=0.95)
#' are equal after scaling with respect to the first TOM matrix.
#'
#' 
wgcna.scaleTOMs <- function(TOMs, scaleP=0.95) {
  nGenes <- nrow(TOMs[[1]])
  nSets <- length(TOMs)
  # Sample sufficiently large number of TOM entries
  nSamples = as.integer(1/(1-scaleP) * 1000);
  # Choose the sampled TOM entries
  scaleSample = sample(nGenes*(nGenes-1)/2, size = nSamples)
  TOMScalingSamples = list();
  # These are TOM values at reference percentile
  scaleQuant = rep(1, nSets)
  # Scaling powers to equalize reference TOM values
  scalePowers = rep(1, nSets)
  # Loop over sets
  set=1
  for (set in 1:nSets)
  {
    # Select the sampled TOM entries
    tval = as.dist(TOMs[[set]])[scaleSample]
    # Calculate the 95th percentile
    scaleQuant[set] = quantile(tval, probs = scaleP, type = 8);
    TOMScalingSamples[[set]] <- tval

    # Scale the TOM
    if (set>1)
    {
      scalePowers[set] = log(scaleQuant[1])/log(scaleQuant[set]);
      TOMs[[set]] = TOMs[[set]]^scalePowers[set];
    }
  }
  return(TOMs)
}

