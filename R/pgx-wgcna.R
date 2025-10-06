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
    ngenes = 2000,
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
  if (is.null(wTOM) && "svTOM" %in% names(wgcna)) wTOM <- wgcna$svTOM %*% t(wgcna$svTOM)
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
  res.gse <- wgcna.computeModuleEnrichment(
    wgcna,
    annot = pgx$genes,
    # GMT = pgx$GMT,
    # gsetX = pgx$gsetX,
    methods = c("fisher", "gsetcor", "xcor"),
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
                          verbose = 0
                          ) {

  if(0) {
    ngenes = 2000;
    minmodsize = 20;
    power = 12;
    mergeCutHeight = 0.15;                          
    deepsplit = 2;
    minKME = 0.3;
    treeCut = 0.99;
    treeCutCeiling = 1;
    networktype = "signed";
    tomtype = "signed";
    clustMethod = "average";
    cutMethod = "hybrid";
    calcMethod = "fast";
    lowrank = 40;
    numericlabels = FALSE;
    maxBlockSize = 9999;
    merge.dendro = TRUE;
    compute.stats = TRUE;
    prefix = "ME";
    sv.tom = 40;
    drop.ref = FALSE;
    net = NULL;
    is.multiomics = NULL;
    verbose = 0
  }
    
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
    if(is.null(is.multiomics)) is.multiomics <- all(grepl(":",rownames(X))) ## not robust!!
    if(is.multiomics) {
      message("[wgcna.compute] topSD = ",ngenes, " (multi-omics)")
      X <- mofa.topSD(X, ngenes)    
    } else {
      message("[wgcna.compute] topSD = ", ngenes, " (single omics)")
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
  if(!is.null(power)) power <- power[1]
  
  message("[wgcna.compute] minmodsize = ", minmodsize)
  message("[wgcna.compute] number of features = ", nrow(X))
  message("[wgcna.compute] minKME = ", minKME)
  message("[wgcna.compute] power = ", power)
  message("[wgcna.compute] mergeCutHeight = ", mergeCutHeight)
  message("[wgcna.compute] calcMethod = ", calcMethod)  
  
  ##WGCNA::enableWGCNAThreads()
  if(is.null(net)) {
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
      verbose = verbose
    )
  }

  if(!"MEs" %in% names(net)) {
    net$MEs = WGCNA::moduleEigengenes(datExpr, colors = net$colors)$eigengenes
  }
  
  ## Substitue prefix="ME"
  if(prefix!="ME") names(net$MEs) <- sub("^ME", prefix, names(net$MEs))
  net$labels <- paste0(prefix, net$colors)
  names(net$labels) <- colnames(datExpr)
  
  ## clean up traits matrix
  datTraits <- data.frame(samples, check.names = FALSE)
  isdate <- apply(datTraits, 2, is.Date)
  datTraits <- datTraits[, !isdate, drop = FALSE]

  ## Expand multi-class discrete phenotypes into binary vectors
  datTraits <- utils::type.convert(datTraits, as.is = TRUE)
  datTraits <- expandPhenoMatrix(datTraits, keep.numeric=TRUE,
    drop.ref=drop.ref)

  if(is.null(datTraits)) {
    message("WARNING:: no valid traits. creating random traits.")
    ##random.trait <- sample(c(0,1), nrow(samples), replace=TRUE)
    random.trait <- head(rep(c(0,1), nrow(samples)), nrow(samples))
    datTraits <- data.frame( random = random.trait )
    rownames(datTraits) <- rownames(samples)
  }
  
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
    stats <- wgcna.computeGeneStats(net, datExpr, datTraits, TOM=TOM) 
  }
  
  ## module-traits matrix
  message("[wgcna.compute] computing module-traits matrix...")
  modTraits <- cor(net$MEs, datTraits, use = "pairwise")

  ## merge dendrograms
  if (merge.dendro) {
    message("[wgcna.compute] merge_block_dendrograms...")
    merged <- try(wgcna.merge_block_dendrograms(net, TOM))
    if (!inherits(merged, "try-error")) {
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
                                     ngenes = 2000,
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
                                     drop.ref = FALSE,
                                     add.pheno = FALSE,
                                     add.gsets = FALSE,
                                     do.consensus = FALSE,
                                     gset.methods = c("fisher","gsetcor","xcor"),
                                     verbose = 1,
                                     progress = NULL
                                     ) {

  if(0) {
    do.consensus = 1
    cutMethod = "hybrid"
    deepsplit = 2
    power = 12
    ngenes = 2000
    minmodsize = 10
    minKME = 0.3
    compute.enrichment = TRUE
    xtop = 100
    annot = pgx$genes
    GMT = pgx$GMT  ##??
    gsetX = pgx$gsetX
    progress = NULL
  }
  
  ## preprocessing
  if(!is.null(annot)) {
    dataX <- lapply(dataX, function(x) rename_by2(x, annot, "symbol"))    
  }

  ## add phenomatrix??
  if(add.gsets) {
    if(is.null(GMT)) {
      GMT <- Matrix::t(playdata::GSETxGENE)
    }
    if(is.null(gsetX)) {
      X <- do.call(rbind, dataX)
      kk <- intersect(rownames(X), rownames(GMT))
      if(length(kk)) {
        gsetX <- plaid::plaid(X[kk,], GMT[kk,])
      }
    }
    if(!is.null(gsetX)) {
      dataX$gs <- gsetX
    }
  }

  ## add phenomatrix??
  if(add.pheno) {
    phenoX <- expandPhenoMatrix(samples, keep.numeric=TRUE, drop.ref=drop.ref)
    dataX$ph <- t(phenoX)
  }

  dataX <- lapply(dataX, imputeMissing, method="SVD2")
  names(dataX) <- substring(names(dataX),1,2)

  lapply(dataX, dim)
  #dataX <- mofa.topSD(dataX, ngenes)
  
  if(!is.null(progress)) {
    progress$set(message = paste("computing WGCNA modules..."), value = 0.33)
  }

  if(is.null(power)) power <- "sft"
  if(as.character(power[1]) %in% c("sft","iqr")) {
    ## Estimate best power
    message("[wgcna.compute] estimating optimal power with method = ", power[1])
    est.power <- rep(NA, length(dataX))    
    i=1
    for(i in 1:length(dataX)) {
      p <- wgcna.pickSoftThreshold(
        Matrix::t(dataX[[i]]), sft=NULL, rcut=0.85, powers = NULL,
        method=power[1], nmax=1000, verbose=1)
      if(length(p)==0 || is.null(p) ) p <- NA
      est.power[i] <- p
    }
    est.power
    power <- ifelse( is.na(est.power), 12, est.power)
  } else {
    power <- as.numeric(power)
  }
  nw <- length(dataX)
  power <- head(rep(power, nw),nw)
  names(power) <- names(dataX)
  
  ## This runs WGCNA on an expression list. 
  wgcna <- list()
  has.gxpx <- all(c("gx","px") %in% names(dataX))
  if( do.consensus && has.gxpx ) {
    cat("----------- computing consensus WGCNA for GX+PX ------------\n")
    gx <- dataX[['gx']]
    px <- dataX[['px']]    
    nn <- mean(rownames(gx) %in% rownames(px))
    if( nn < 0.10) {
      message("ERROR: gx and px features do not overlap")
    } else {
      wgcna <- wgcna.createConsensusLayers(
        dataX[c('gx','px')],
        samples = samples,
        prefix = c('GX','PX'),
        ngenes = ngenes,
        power = power[c('gx','px')],
        minModuleSize = minmodsize,
        deepSplit = deepsplit,
        mergeCutHeight = 0.15,
        minKME = minKME,
        maxBlockSize = 9999,
        verbose = 1
      )
    }
  }

  names(wgcna)
  dtlist <- setdiff(names(dataX), names(wgcna))
  dtlist
  for(dt in dtlist) {
    cat("------------ computing WGCNA for",dt,"-------------\n")
    minKME <- ifelse(dt=='ph', 0, minKME)
    minmodsize <- ifelse( dt=='ph', 1, minmodsize)      
    
    wgcna[[dt]] <- wgcna.compute(
      X = dataX[[dt]],
      samples = samples, 
      ngenes = ngenes,
      calcMethod = "fast",
      power = power[dt],
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
      drop.ref = drop.ref,
      is.multiomics = FALSE,
      verbose = verbose
    )
  }

  wgcna <- wgcna[names(dataX)]
  
  ## Compute enrichment
  if(compute.enrichment) {
    message("computing module enrichment...")
    if(!is.null(progress)) {
      progress$set(message = paste("computing module enrichment..."), value = 0.66)
    }
    
    gse <- wgcna.computeModuleEnrichment(
      wgcna = wgcna,
      multi = TRUE,
      methods = gset.methods,
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
  calcMethod = "fast",  
  lowrank = 20,
  clustMethod = "average",
  cutMethod = "hybrid",  ## hybrid, tree, static
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
  verbose = 1 ) {

  #power=6;networkType=TOMType="signed";minModuleSize=20;mergeCutHeight=0.15;minKMEtoStay=0.3;numericLabels=FALSE;clustMethod="average";deepSplit = 2;treeCut = 0.99;treeCutCeiling = 1;
  
  cor <- WGCNA::cor  ## needed...
  deepSplit <- as.integer(deepSplit)
  lowrank <- as.integer(lowrank)

  if(is.null(power)) power <- "sft"
  auto.power <- power[1] %in% c("sft","iqr")
  if (auto.power) {
    ## Estimate best power
    message("[wgcna.compute] estimating optimal power with method = ", power[1])
    powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
    powers <- c(powers, seq(from = 25, to = 50, by = 5))
    power <- wgcna.pickSoftThreshold(datExpr, sft=NULL, rcut=0.85,
      method=power[1], nmax=2000, verbose=0) 
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

  if(verbose>0) {
    message("Tree cut method: ", cutMethod)
    message("treeCut = ", treeCut)
    message("deepSplit = ", deepSplit)
    message("minModuleSize = ", minModuleSize)
  }
  
  if(cutMethod %in% c("hybrid","tree")) {

    if(cutMethod=="tree") deepSplit <- (deepSplit > 0)
    label <- dynamicTreeCut::cutreeDynamic(
      geneTree,
      method = cutMethod,
      cutHeight = treeCut,
      distM = dissTOM,
      deepSplit = deepSplit,
      minClusterSize = minModuleSize2
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

#' Compute general feature statistics after WGCNA results.
#'
#'
#' 
wgcna.computeGeneStats <- function(net, datExpr, datTraits, TOM) {

  ## align
  kk <- intersect(rownames(datExpr), rownames(datTraits))
  datExpr <- datExpr[kk,]
  datTraits <- datTraits[kk,]

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
  TSPvalue <- WGCNA::corPvalueStudent(as.matrix(traitSignificance), nSamples)

  ## Fold-change
  foldChange <- NULL
  foldChangePvalue <- NULL  
  is.binary <- apply(datTraits, 2, function(a) length(unique(a[!is.na(a)])) == 2)
  is.binary
  binY <- NULL
  if(any(is.binary)) {
    binY <- datTraits[, which(is.binary), drop = FALSE]
    nmin <- apply(binY, 2, function(x) min(table(x)))
    binY <- binY[,which(nmin>=2),drop=FALSE]
  }
  if(!is.null(binY) && NCOL(binY)>0) {
    lm <- list()
    i=1
    for(i in 1:ncol(binY)) {
      y <- 1*binY[,i]
      X <- t(datExpr)
      suppressWarnings(suppressMessages(
        res <- try(gx.limma(X, y, lfc=0, fdr=1, sort.by="none", verbose=0, max.na=1))
      ))
      if(!"try-error" %in% class(res)) {
        k <- colnames(binY)[i]
        lm[[k]] <- res
      }
    }
    lm <- lm[!sapply(lm,is.null)]
    foldChange <- sapply(lm, function(m) m$logFC)
    foldChangePvalue <- sapply(lm, function(m) m$P.Value)
    if(length(lm)==1) {
      foldChange <- cbind(foldChange)
      foldChangePvalue <- cbind(foldChangePvalue)
    }
    rownames(foldChange) <- rownames(lm[[1]])
    rownames(foldChangePvalue) <- rownames(lm[[1]])
  }
  
  # Continuous traits (not always present)
  contY <- datTraits[, which(!is.binary), drop = FALSE]
  foldChange.cont <- NULL
  foldChangePvalue.cont <- NULL  
  dim(contY)
  if(NCOL(contY) > 0) {
    contlm <- apply(contY, 2, function(y) {
      rho <- cor(datExpr, y, use="pairwise")[,1]
      P.Value <- WGCNA::corPvalueStudent(rho, n = length(y))
      data.frame(rho, P.Value)
    })
    contlm <- contlm[!sapply(contlm, is.null)]
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
  geneCentrality <- NULL
  if(!is.null(TOM)) {
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
  }

  ## force align. Sometime they are shorter for some reason...
  gg <- rownames(moduleMembership)
  matMatch <- function(m, gg) {
    if(is.null(m)) return(NULL)
    m <- m[match(gg,rownames(m)),,drop=FALSE]
    rownames(m) <- gg
    return(m)
  }
  
  ##moduleMembership <- matMatch(moduleMembership)
  MMPvalue <- matMatch(MMPvalue, gg)
  traitSignificance = matMatch(traitSignificance, gg)
  TSPvalue = matMatch(TSPvalue, gg)
  foldChange = matMatch(foldChange, gg)
  foldChangePvalue = matMatch(foldChangePvalue, gg)
  
  stats <- list(
    moduleTraitCor = moduleTraitCor,
    moduleTraitPvalue = moduleTraitPvalue,
    moduleMembership = moduleMembership,
    MMPvalue = MMPvalue,
    traitSignificance = traitSignificance,
    TSPvalue = TSPvalue,
    foldChange = foldChange,
    foldChangePvalue = foldChangePvalue,
    geneCentrality = geneCentrality
  )

  return(stats)
}

#'
#'
#' @export
wgcna.getGeneStats <- function(wgcna, trait, module=NULL, plot = TRUE,
                               stats = NULL, labels = NULL,
                               col = NULL, main = NULL) {

  if(!is.null(stats)) {
    features <- rownames(stats[['moduleMembership']])
    if(is.null(stats)) stop("must give stats")
    if(is.null(labels)) stop("must give labels")    
  } else if(!is.null(wgcna)) {
    labels <- wgcna$net$labels
    features <- names(wgcna$net$colors)
    stats <- wgcna$stats
  } else {
    stop("must supply wgcna or trait")
  }

  mean(labels %in% WGCNA::standardColors())
  setdiff(labels, WGCNA::standardColors())

  is.color <- mean(labels %in% c("grey",WGCNA::standardColors())) > 0.9
  if(is.color) {
    prefix <- substring(rownames(stats[['moduleTraitCor']]),1,2)[1]
    labels <- paste0(prefix, labels)
  }
  
  p1 <- c("moduleMembership", "MMPvalue")
  p2 <- c("traitSignificance", "TSPvalue", "foldChange", "foldChangePvalue")
  p3 <- c("geneCentrality")

  df <- data.frame(
    feature = features,
    module = labels
  )
  head(df)
  
  ## get moduleMembership
  mm.stats <- stats[p1]
  mm.label <- colnames(mm.stats[[1]])
  idx <- cbind(1:length(labels), match(labels, mm.label))
  A1 <- sapply( mm.stats, function(x) x[idx])  
  rownames(A1) <- labels
  df <- cbind(df, A1)

  ## get traitSig columns for trait
  tt.cols <- colnames(stats[[p2[1]]])  
  if(is.null(trait)) trait <- tt.cols
  trait <- intersect(trait, tt.cols)
  
  if (length(trait)>1) {
    A2 <- lapply(stats[p2], function(x) x[,trait])
    for(i in 1:length(A2)) colnames(A2[[i]]) <- paste0(names(A2)[i],".",colnames(A2[[i]]))
    A2 <- do.call(cbind, A2)
    df <- cbind(df, A2)
  } else if(length(trait)==1) {
    A2 <- lapply(stats[p2], function(x) x[,trait])
    A2 <- do.call(cbind, A2)
    df <- cbind(df, A2)
  } else {
    message("[wgcna.getGeneStats] ERROR: trait not in stats object")
    return(NULL)
  }
  
  A3 <- stats[[p3]]
  if(!is.null(A3)) {
    df <- cbind(df, centrality = A3)
  }
  
  ## calculate score
  sel <- c("moduleMembership", "traitSignificance", "foldChange", "centrality")
  sel <- intersect(sel, colnames(df))
  #df1 <- pmax(as.matrix(df[, sel]), 1e-8)
  df1 <- as.matrix(abs(df[, sel]))
  score <- exp(rowMeans(log(1e-8 + df1))) * sign(df[,"foldChange"])

  df <- data.frame(df[,1:2], score = score, df[,-c(1,2)])
  score.sign <- sign(median(df$score,na.rm=TRUE))
  df <- df[order(-df$score * score.sign), ]
  
  if (!is.null(module)) {
    ##sel <- which(df$module == module)
    sel <- grep( paste0(module,"$"), df$module)
    df <- df[sel, , drop = FALSE]
  }

  if (plot) {
    cols <- c("moduleMembership", "traitSignificance", "foldChange", "centrality")
    cols <- intersect(cols, colnames(df))
    df1 <- df[, cols]
    col1 <- wgcna.labels2colors(labels[rownames(df1)])
    if (!is.null(col)) col1 <- col
    pairs(df1, col = col1, oma=c(1,1,1,1)*1.8)
    if (is.null(main)) {
      main <- paste("Gene significance for module", module, "and trait", trait)
    }
    title(main, line = 3, cex.main = 1.15)
  }

  rownames(df) <- df$feature
  df
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
  for(k in names(cons$layers)) {
    w <- cons$layers[[k]]
    colors <- cons$net$colors
    wMEs <- cons$net$multiMEs[[k]]$data
    wnet <- list( MEs = wMEs, colors = colors)
    stats[[k]] <- wgcna.computeGeneStats(
      wnet, w$datExpr, w$datTraits, TOM=NULL)
  }
  return(stats)
}

#'
#'
#' @export
wgcna.getConsensusGeneStats <- function(cons, stats, trait, module=NULL) {

  ## create extended color vector
  labels = paste0("ME",cons$net$colors)
  gstats <- list()
  for(k in names(stats)) {
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

  ## Align rows
  ff <- gstats[[1]]$feature
  for(k in names(gstats)) {
    ii <- match(ff, gstats[[k]]$feature)
    gstats[[k]] <- gstats[[k]][ii,]
  }

  ## Compute consensus statistics
  xcols <- c(3,4,6,8)
  pcols <- c(10,5,7,9)
  pcols1 <- c(5,7,9)
  xcols <- c("score","moduleMembership","traitSignificance","foldChange")
  pcols <- c("scorePvalue","MMPvalue","TSPvalue","foldChangePvalue")
  pcols1 <- pcols[-1]
  for(i in 1:length(gstats)) {
    gstats[[i]][,'scorePvalue'] <- apply(gstats[[i]][,pcols1],1,max,na.rm=TRUE)
  }
  ##xc <- lapply(gstats, function(x) log(abs(x[,xcols])*(x[,pcols]<0.05)))
  xc <- lapply(gstats, function(x) log(abs(x[,xcols])))
  xc <- exp(Reduce('+', xc) / length(xc))
  xp <- Reduce(pmax, lapply(gstats, function(x) x[,pcols]))
  df3 <- data.frame( gstats[[1]][,1:2], xc, xp)
  df3 <- df3[,colnames(gstats[[1]])]  
  head(df3)  
  sign.pos <- Reduce('*',lapply(gstats,function(g) sign(g$score) == 1))
  sign.neg <- Reduce('*',lapply(gstats,function(g) sign(g$score) == -1))
  allsig   <- Reduce('*',lapply(gstats,function(g) (g$scorePvalue) < 0.05))    
  table(allsig)  
  consensus <- c("D","C")[ 1 + 1*(sign.pos | sign.neg)]
  consensus[which(allsig==0)] <- 'N'
  cons.df <- data.frame(df3[,1:2], consensus, df3[,-c(1,2)])
  head(cons.df)
  
  ## This creates the full stats matrix (all subgroups)
  df1 <- gstats[[1]][,c("feature","module")]
  df2 <- gstats[[1]][,0]
  cols <- colnames(gstats[[1]])[-c(1:2)]
  for(k in cols ) {
    xx <- sapply(gstats, function(g) g[,k])
    df2[[k]] <- I(xx)
  }  
  df2 <- do.call(cbind, lapply(df2,unclass))
  newcols <- unlist(lapply(cols, function(k) paste0(k,'.',names(gstats))))
  colnames(df2) <- newcols
  full.df <- data.frame(df1, consensus=cons.df$consensus, df2)


  ## sort??
  ii <- order(-cons.df$score * sign(mean(cons.df$score,na.rm=TRUE)))
  cons.df <- cons.df[ii,]
  full.df <- full.df[ii,]
  
  list(
    consensus = cons.df,
    full = full.df
  )
}


## ----------------------------------------------------
## Perform geneset analysis on modules
## ----------------------------------------------------

wgcna.computeModuleEnrichment <- function(wgcna,
                                          multi = FALSE,
                                          methods = c("fisher","gsetcor","xcor"),
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
  
  if(is.null(annot)) {
    gg <- lapply(wgcna, function(w) colnames(w$datExpr))
    gg1 <- lapply(names(wgcna), function(a) paste0(a,":",gg[[a]]))
    gg <- as.character(unlist(gg))
    gg1 <- as.character(unlist(gg1))    
    annot <- data.frame(feature = gg1, symbol = gg)
  }
    
  if(is.null(GMT)) {
    message("Using playdata GSETxGENE genesets")
    GMT <- Matrix::t(playdata::GSETxGENE)
  }

  if(!is.null(annot)) {
    symbol.col <- intersect(c("symbol","gene_name"),colnames(annot))[1]
    geneX <- rename_by2(geneX, annot, symbol.col)
    GMT <- rename_by2( GMT, annot, symbol.col)
  }
  if(length(intersect(rownames(geneX),rownames(GMT)))==0) {
    message("ERROR: no overlap for geneX and GMT features. Please add annotation.")
    return(NULL)
  }
  
  bg <- rownames(geneX)
  bg <- intersect(bg, rownames(GMT))
  if (length(bg) == 0) {
    message("FATAL. no overlapping genes")
    return(NULL)
  }
  G1 <- GMT[bg, , drop = FALSE]
  if (!is.null(filter)) {
    sel <- grep(filter, colnames(G1))
    if (length(sel)) G1 <- G1[, sel, drop = FALSE]
  }
  G1 <- G1[, which(Matrix::colSums(G1 != 0) >= 4), drop = FALSE]

  if(is.null(gsetX)) {
    gsetX <- plaid::plaid( geneX, G1)
  }

  ## align dimensions
  ss <- intersect(rownames(gsetX), colnames(G1))
  G1 <- G1[,ss,drop=FALSE]
  gsetX <- gsetX[ss,]
  
  ME <- lapply(wgcna, function(w) as.matrix(w$net$MEs))
  me.genes <- lapply(wgcna, function(w) w$me.genes)
  names(me.genes) <- NULL
  me.genes <- do.call(c, me.genes)
  me.genes <- lapply(me.genes, function(g) probe2symbol(g, annot, query = symbol.col))

  ## compute most correlated gx/px genes. limit xtop if geneX is too
  ## small
  xtop <- min(xtop, round(nrow(geneX)/4))
  nbx.genes <- list()
  for(i in 1:length(wgcna)) {
    nbx.cor <- cor(t(geneX), ME[[i]])
    for(k in colnames(nbx.cor)) {
      if(is.null(xtop)) {
        n <- length(wgcna[[i]]$me.genes[[k]])
      } else {
        n <- xtop
      }
      nbx.genes[[k]] <- head(names(sort(-nbx.cor[,k])), n)
    }
  }

  ## add nearest expression neighbors to module genes
  nbx.genes <- nbx.genes[names(me.genes)]
  for(k in names(me.genes)) {
    me.genes[[k]] <- unique(c(me.genes[[k]], nbx.genes[[k]]))
  }

  ## make single ME matrix
  MEx <- do.call(cbind, ME)

  gse <- wgcna.run_enrichment_methods(
    ME = MEx,
    me.genes = me.genes,
    G = G1,
    geneX = geneX,
    gsetX = gsetX,
    methods = methods,
    ntop = ntop
  )
  
  return(gse)
}

#'
#'
#' @export
wgcna.getModuleCrossGenes <-  function(wgcna, ref=NULL, ngenes = 100,
                                       multi=TRUE, modules=NULL)
{

  if(!multi) {
    wgcna <- list(gx = wgcna)
    ref <- 'gx'
  }
      
  if(is.null(ref)) ref <- head(intersect(names(wgcna),c("gx","px")),1)
  if(is.null(ref) || !ref %in% names(wgcna)) ref <- names(wgcna)[1]

  W <- wgcna[[ref]]
  geneX <- W$datExpr

  MEx <- sapply(wgcna, function(w) as.matrix(w$net$MEs))
  MEx <- do.call(cbind, MEx)

  if(!is.null(modules)) {
    modules <- intersect(modules, colnames(MEx))
    MEx <- MEx[,modules,drop=FALSE]
  }
  
  nbx.cor <- cor(geneX, MEx)

  nbx.list <- list()
  for(k in colnames(nbx.cor)) {
    ii <- head(order(-nbx.cor[,k]), ngenes)
    rho <- nbx.cor[ii,k]
    gene <- rownames(nbx.cor)[ii]
    me <- W$net$labels[gene]
    nbx.list[[k]] <- data.frame( gene = gene, rho = rho, module = me)
  }

  ##if(length(nbx.list)==1) nbx.list <- nbx.list[[1]]
  return(nbx.list)
}


#'
#'
#' 
wgcna.computeConsensusModuleEnrichment <-
  function(cons,
           GMT,
           annot,
           methods = c("fisher","gsetcor","xcor"),
           min.genes = 10,
           ntop = 400 )
  {
    
    lapply(cons$datExpr, dim)
    geneX <- t(do.call(rbind, cons$datExpr))

    ## batch correct on datasets using ComBat
    batch <- rep(1:length(cons$datExpr), sapply(cons$datExpr, nrow))
    geneX <- sva::ComBat( geneX, batch = batch )
    
    ## for correlation-based enrichment, we used scaledX
    scaledX <- t(do.call(rbind, lapply(cons$datExpr, function(x) scale(x))))
    
    ## Rename everything to symbols
    if(!is.null(annot)) {
      geneX <- rename_by2(geneX, annot, "symbol")
      GMT   <- rename_by2(GMT, annot, "symbol")
    }
    ng <- length(intersect(rownames(geneX), rownames(GMT)))
    if(ng == 0) {
      message("[wgcna.computeConsensusModuleEnrichment] ERROR. No symbol overlap.")
      return(NULL)
    }
    symbols <- intersect(rownames(GMT),rownames(geneX))
    message("[wgcna.computeConsensusModuleEnrichment] number of symbols: ", length(symbols))
    geneX <- geneX[symbols,]
    GMT <- GMT[symbols,]
    
    ## select on minimum gene sets size
    sel <- which( Matrix::colSums(GMT!=0) >= min.genes )
    GMT <- GMT[,sel]
    
    ## Compute single-samples gene set expression
    gsetX <- plaid::plaid(geneX, matG=GMT)
    dim(gsetX)
    
    ## Create extended Eigengene matrix (ME). ME should be nicely
    ## normalized/scaled so we just rbind across datasets
    ME <- lapply( cons$net$multiMEs, function(m) m$data)
    ME <- do.call(rbind, ME)
    dim(ME)
    
    ## get genes in modules
    me.genes <- tapply(names(cons$net$colors), cons$net$colors, list)
    names(me.genes) <- paste0("ME",names(me.genes))
    if(!is.null(annot)) {
      me.genes <- lapply(me.genes, function(gg) probe2symbol(gg, annot))
    }
    me.genes <- lapply(me.genes, function(g) intersect(g, symbols))
    
    colnames(scaledX) <- rownames(ME)
    colnames(gsetX) <- rownames(ME)

    geneX=scaledX    

    gseaX <- wgcna.run_enrichment_methods(
      ME, me.genes = me.genes, GMT=GMT,
      geneX=scaledX, gsetX=gsetX, methods = methods,
      min.genes = min.genes, ntop = ntop)
    
    return(gseaX)
  }


wgcna.run_enrichment_methods <- function(ME, me.genes, GMT, geneX, gsetX,
                                         methods = c("fisher","gsetcor","xcor"),
                                         ntop = 400, min.genes = 3
                                         )
{

  rho.list <- list()
  pval.list <- list()

  ## align matrices
  bg <- intersect(rownames(GMT), rownames(geneX))
  ss <- intersect(colnames(GMT), rownames(gsetX))
  GMT <- GMT[bg,ss]
  geneX <- geneX[bg,]
  gsetX <- gsetX[ss,]

  ## select on minimum genes
  sel <- which( Matrix::colSums(GMT!=0) >= min.genes )
  GMT <- GMT[,ss]
  gsetX <- gsetX[ss,]
  
  message("Computing enrichment for ", length(ss), " genesets")

  ## Here we correlate geneset score (averageCLR) with the module
  ## eigengene (ME). This should select genesets correlated with the
  ## ME.
  if ("gsetcor" %in% methods && !is.null(gsetX)) {
    message("[wgcna.computeModuleEnrichment] calculating single-sample geneset correlation...")
    rc.rho <- matrix(NA, ncol(GMT), ncol(ME))
    rc.pvalue <- matrix(NA, ncol(GMT), ncol(ME))    
    dimnames(rc.rho) <- list( colnames(GMT), colnames(ME))
    dimnames(rc.pvalue) <- list( colnames(GMT), colnames(ME))    
    jj <- which(rownames(gsetX) %in% colnames(GMT))
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
  if ("xcor" %in% methods) {
    message("[wgcna.computeModuleEnrichment] calculating eigengene GBA correlation...")
    gba <- cor(t(geneX), ME, use = "pairwise")
    rc <- gset.rankcor(gba, GMT, compute.p = TRUE) ## NEEDS CHECK!!!
    rho.list[["xcor"]] <- rc$rho
    pval.list[["xcor"]] <- rc$p.value
  }

  gmt <- mat2gmt(GMT)
  if(1) {
    ## we pre-select to make this faster
    Pmin <- sapply(pval.list, function(P) apply(P,1,min))
    sel <- head(order(rowMeans(apply(Pmin, 2, rank))), 5*ntop)
    message("[wgcna.computeModuleEnrichment] pre-selecting ",length(sel)," sets for fgsea/Fisher test...")
    sel <- rownames(Pmin)[sel]
    gmt <- gmt[sel]
  }    
    
  ## fGSEA
  if ("fgsea" %in% methods) {
    message("[wgcna.computeModuleEnrichment] calculating module fgsea...")
    xrho <- cor(t(geneX), ME, use = "pairwise")
    res <- list()
    i=1
    for(i in 1:ncol(xrho)) {
      k <- colnames(xrho)[i]
      res[[k]] <- fgsea::fgsea(gmt, xrho[,i]) ## NEEDS CHECK!!!
    }
    pw <- res[[1]]$pathway
    res <- lapply(res, function(r) r[match(pw,r$pathway),])
    nes <- sapply(res, function(r) r$NES)
    pval <- sapply(res, function(r) r$pval)
    rownames(nes) <- rownames(pval) <- pw
    colnames(nes) <- colnames(pval) <- names(res)
    rho.list[["fgsea"]] <- nes
    pval.list[["fgsea"]] <- pval
  }
  
  ## Perform fisher-test on ME genes
  if ("fisher" %in% methods) {

    message("[wgcna.computeModuleEnrichment] calculating Fisher tests...")    
    rho <- matrix(NA, length(gmt), ncol(ME))
    pval <- matrix(NA, length(gmt), ncol(ME))    
    dimnames(rho) <- list( names(gmt), colnames(ME))
    dimnames(pval) <- list( names(gmt), colnames(ME))    

    ## perform Fisher test for all modules using the module genes
    i=1
    for (i in 1:ncol(rho)) {
      k <- colnames(rho)[i]
      gg <- me.genes[[k]]      
      rr <- try(gset.fisher(gg, gmt, background = bg, fdr = 1,
        min.genes = -1, verbose = 0, sort.by='none', no.pass=1))

      if (!"try-error" %in% class(rr)) {        
        rr <- rr[match(rownames(rho), rownames(rr)), ]
        rho[,i] <- rr$odd.ratio
        pval[,i] <- rr$p.value
      }
    }
    
    ## handle infinite or NA
    rho[is.infinite(rho)] <- 2*max(rho, na.rm=TRUE)  ## Inf odd.ratio
    pval[is.na(pval)] <- 1
    rho[is.na(rho)] <- 0

    rho.list[["fisher"]] <- rho
    pval.list[["fisher"]] <- pval
  }

  lapply(rho.list, dim)
  
  ## ensure dimensions
  gsets <- Reduce(intersect, lapply( rho.list, rownames))
  modules <- Reduce(intersect, lapply( rho.list, colnames))  
  rho.list <- lapply( rho.list, function(x) x[gsets,modules])
  pval.list <- lapply( pval.list, function(x) x[gsets,modules])  

  lapply(rho.list, dim)
  
  ## Compute meta rank and pval. Handle NA for failing methods.
  pvalNA <- lapply(pval.list, function(x) {x[is.na(x)]=0;x})
  ##pvalNA <- lapply(pval.list, function(x) {x[is.na(x)]=1;x})
  meta.p <- Reduce(pmax, pvalNA) ## NEED RETHINK!!!
  meta.q <- apply(meta.p, 2, p.adjust, method = "fdr")

  ## NEED RETHINK: how about negative FC???
  rnk.list <- lapply(rho.list, function(x) apply(x, 2, rank, na.last='keep') / nrow(x))
  meta.rnk <- Reduce("+", rnk.list) / length(rnk.list)
  rnk.NAZERO  <- lapply(rnk.list, function(x) {x[is.na(x)]=0;x})
  rnk.NSUM <- Reduce('+', lapply(rnk.list, function(x) !is.na(x)))
  meta.rnk <- Reduce('+', rnk.NAZERO) / rnk.NSUM

  ## create dataframe by module
  message("[wgcna.computeModuleEnrichment] creating dataframes...")
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


wgcna.merge_block_dendrograms <- function(net, X, method = 1) {
  ## This function is fragile: it can give a C stack limit error. In
  ## case that happens you can increase the stack limit by running on
  ## the cmd line: >>> ulimit -s unlimited
  ##
  hc <- net$dendrograms

  ## merge block dendrogram
  mx <- list()
  for (b in 1:length(net$dendrograms)) {
    ii <- which(net$goodGenes & net$blocks == b)
    mx[[b]] <- colMeans(X[ii, ])
    hc[[b]]$labels <- rownames(X)[ii]
  }

  if (length(hc) == 1) {
    return(hc[[1]])
  }

  ## compute parent dendrogram
  M <- do.call(rbind, mx)
  hclust_p <- hclust(dist(M), method = "average")
  dend_p <- as.dendrogram(hclust_p)
  dend.list <- lapply(hc, as.dendrogram)

  if (method == 1) {
    merged <- ComplexHeatmap::merge_dendrogram(dend_p, dend.list)
  } else {
    mrg <- hclust_p$merge
    merged_branch <- list()
    k <- 1
    for (k in 1:nrow(mrg)) {
      i <- mrg[k, 1]
      j <- mrg[k, 2]
      if (i < 0) d1 <- dend.list[[-i]]
      if (i > 0) d1 <- merged_branch[[i]]
      if (j < 0) d2 <- dend.list[[-j]]
      if (j > 0) d2 <- merged_branch[[j]]
      merged_branch[[k]] <- merge(d1, d2) ## actual merge
    }
    merged <- merged_branch[[k]]
  }

  merged_hclust <- as.hclust(merged)
  merged_hclust
}


##=========================================================================
## CONSENSUS WGCNA
##=========================================================================

#'
#'
#' @export
wgcna.runConsensusWGCNA <- function(exprList,
                                    phenoData,
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
                                    compute.stats = TRUE,
                                    compute.enrichment = TRUE,
                                    gsea.mingenes = 10,
                                    gset.methods = c("fisher","gsetcor","xcor"),
                                    verbose = 1
                                    ) {
  if(0) {
    power=6;minKME=0.5;cutheight=0.15;deepSplit=2;maxBlockSize=5000;verbose=1;calcMethod="fast";addCombined=0;ngenes=2000;minModuleSize=20;mergeCutHeight=0.15
    gsea.mingenes=20;gset.methods = c("fisher","gsetcor","xcor")
  }

  colors <- NULL
  
  ## Align and reduce matrices if needed
  gg <- Reduce(intersect, lapply(exprList,rownames))
  exprList <- lapply(exprList, function(x) x[gg,])
  if( length(gg) > ngenes ) {
    sdx <- Reduce('*',lapply(exprList, function(x) matrixStats::rowSds(x)))
    ii <- head(order(-sdx), ngenes)
    exprList <- lapply(exprList, function(x) x[ii,])
  }
  
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
  layers <- list()
  k=names(multiExpr)[1]
  for (i in 1:length(multiExpr)) {
    k <- names(multiExpr)[i]
    message("[wgcna.runConsensusWGCNA] >>> computing WGCNA for ", k)
    X = Matrix::t(multiExpr[[i]]$data)
    layers[[k]] <- wgcna.compute(    
      X = X,
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
      compute.stats = compute.stats,
      sv.tom = 40,
      verbose = verbose
    )
  }
  
  # now we run automatic consensus module detection
  message("[wgcna.runConsensusWGCNA] >>> computing CONSENSUS modules...")
  consensusPower <- unlist(sapply(layers, function(w) w$net$power))
  if(is.null(consensusPower) && !is.null(power)) {
    consensusPower <- power
  }
  if(is.null(consensusPower)) {
    consensusPower <- rep(12, length(layers))
  }
  
  sel <- setdiff(names(multiExpr), c("Combined"))
  cons = WGCNA::blockwiseConsensusModules(
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
    verbose = verbose
  )
  cons$power = consensusPower
  
  ## create and match colors
  i=1
  for(i in 1:length(layers)) {
    layers[[i]] <- wgcna.matchColors(layers[[i]], cons$colors) 
  }

  #
  layers.colors <- sapply( layers, function(r) r$net$colors )
  colors <- cbind( Consensus = cons$colors, layers.colors)
    
  ## add labels to dendrogram
  for(i in 1:length(cons$dendrograms)) {
    ii <- which(cons$goodGenes & cons$blocks==i)
    xnames <- names(cons$colors)
    cons$dendrograms[[i]]$labels <- xnames[ii]
  }
  
  ## merge dendrograms ????
  message("[wgcna.runConsensusWGCNA] merge_block_dendrograms...")  
  multiX <- Matrix::t(do.call(rbind,lapply(exprList,function(x)scale(t(x)))))
  merged <- try(wgcna.merge_block_dendrograms(cons, multiX))
  if(!inherits(merged,"try-error")) {
    cons$merged_dendro <- merged
  } else {
    cons$merged_dendro <- NULL
  }
  
  ## create module-trait matrices for each set
  message("[wgcna.runConsensusWGCNA] >>> computing module-traits matrices...")
  datTraits <- 1*expandPhenoMatrix(
    phenoData,
    drop.ref = drop.ref,
    keep.numeric = TRUE 
  )
  zlist <- list()
  k=1
  for(k in names(cons$multiME)) {
    M <- (cons$multiME[[k]][[1]])
    Z <- datTraits
    kk <- intersect(rownames(M),rownames(Z))
    zrho <- cor(M[kk,], Z[kk,], use="pairwise")
    zrho[is.na(zrho)] <- 0  ## NEED RETHINK!!
    zlist[[k]] <- zrho
  }

  ## create consensus module-trait matrix
  ydim <- sapply(exprList, ncol) 
  consZ <- wgcna.computeConsensusMatrix(zlist, ydim=ydim, psig=0.05) 
  
  ## add slots
  datExpr <- lapply(exprList, Matrix::t)
  
  res <- list(
    net = cons,
    layers = layers,
    datExpr = datExpr,
    datTraits = datTraits,
    modTraits = consZ,
    dendro = cons$merged_dendro,    
    colors = colors,
    zlist = zlist,
    ydim = ydim,
    class = "consensus"
  )

  ## run stats
  message("[wgcna.runConsensusWGCNA] >>> computing gene statistics...")  
  res$stats <- wgcna.computeConsensusGeneStats(res)

  ## run enrichment
  if(compute.enrichment) {
    message("[wgcna.runConsensusWGCNA] >>> computing module enrichment...")  
    if(is.null(GMT)) GMT <- Matrix::t(playdata::GSETxGENE)
    res$gsea <- wgcna.computeConsensusModuleEnrichment(
      res,
      GMT = GMT,
      method = gset.methods,
      annot = annot,
      min.genes = gsea.mingenes
    ) 
  }
  
  res
}


wgcna.matchColors <- function(wgcna, refcolors) {
 
  oldcolors <- wgcna$net$colors
  newcolors <- WGCNA::matchLabels(oldcolors, refcolors)
  lut <- table(oldcolors, newcolors)
  old2new <- colnames(lut)[max.col(lut)]
  names(old2new) <- rownames(lut)
  prefix <- substring(names(wgcna$me.colors),1,2)[1]
  old2newME <- paste0(prefix,old2new)
  names(old2newME) <- paste0(prefix,names(old2new))
  old2new <- c(old2new, old2newME)

  newcol <- function(x) {
    array( old2new[x], dimnames = list(names(x)))
  }
  
  ## rename everything in net object
  wgcna$net$colors <- newcol(wgcna$net$colors)
  if("labels" %in% names(wgcna$net)) {
    wgcna$net$labels <- newcol(wgcna$net$labels)
  }

  ## rename unmergedColors
  if("unmergedColors" %in% names(wgcna$net)) {
    wgcna$net$unmergedColors <- newcol(wgcna$net$unmergedColors)
  }
  names(wgcna$net$MEs) <- newcol(names(wgcna$net$MEs))
  
  ## rename everything in wgcna object
  names(wgcna$me.genes) <- newcol(names(wgcna$me.genes))
  wgcna$me.colors <- newcol(wgcna$me.colors)
  names(wgcna$me.colors) <- newcol(names(wgcna$me.colors))
  colnames(wgcna$W) <- newcol(colnames(wgcna$W))
  rownames(wgcna$modTraits) <- newcol(rownames(wgcna$modTraits))

  ## rename everything in stats object  
  if("stats" %in% names(wgcna)) {
    rownames(wgcna$stats[['moduleTraitCor']]) <- newcol(rownames(wgcna$stats[['moduleTraitCor']]))
    rownames(wgcna$stats[['moduleTraitPvalue']]) <- newcol(rownames(wgcna$stats[['moduleTraitPvalue']]))
    colnames(wgcna$stats[['moduleMembership']]) <- newcol(colnames(wgcna$stats[['moduleMembership']]))
    colnames(wgcna$stats[['MMPvalue']]) <- newcol(colnames(wgcna$stats[['MMPvalue']]))
  }
  return(wgcna)
}

wgcna.createConsensusLayers <- function(exprList,
                                        samples,
                                        ngenes = 2000,
                                        power = 12,
                                        minModuleSize = 5,
                                        deepSplit = 2,
                                        mergeCutHeight = 0.15,
                                        minKME = 0.3,
                                        maxBlockSize = 9999,
                                        prefix = NULL,
                                        verbose = 1
                                        ) {

  if(0) {
    ngenes = 2000
    power = 12
    minModuleSize = 5
    deepSplit = 2
    mergeCutHeight = 0.15
    minKME = 0.3
    maxBlockSize = 9999
    verbose = 1
    prefix=NULL
  }

  if(is.null(prefix)) prefix <- names(exprList)
  nx <- length(exprList)
  prefix <- head(rep(prefix, nx),nx)
    
  ## reduce
  message("[wgcna.computeConsensusLayers] Aligning matrices...")
  gg <- Reduce(intersect, lapply(exprList,rownames))
  exprList <- lapply(exprList, function(x) x[gg,])

  if( length(gg) > ngenes ) {
    message("[wgcna.computeConsensusLayers] Reducing to ",ngenes," genes")    
    sdx <- Reduce('*',lapply(exprList, function(x) matrixStats::rowSds(x)))
    ii <- head(order(-sdx), ngenes)
    exprList <- lapply(exprList, function(x) x[ii,])
  }
  multiExpr = WGCNA::list2multiData(lapply(exprList, Matrix::t))

  ## determine power vector
  if(is.null(power)) power <- "sft"
  if(as.character(power[1]) %in% c("sft","iqr")) {
    ## Estimate best power
    power <- power[1]
    message("[wgcna.createConsensusLayers] optimal power method = ", power)
    est.power <- rep(NA, length(exprList))    
    i=1
    for(i in 1:length(exprList)) {
      p <- wgcna.pickSoftThreshold(
        Matrix::t(exprList[[i]]), sft=NULL, rcut=0.85, powers = NULL,
        method=power, nmax=1000, verbose=0)
      if(length(p)==0 || is.null(p) ) p <- NA
      est.power[i] <- p
    }
    est.power
    power <- ifelse( is.na(est.power), 12, est.power)
  } else {
    power <- as.numeric(power)
  }
  nw <- length(exprList)
  power <- head(rep(power, nw),nw)
  names(power) <- names(exprList)
  
  message("[wgcna.computeConsensusLayers] Computing consensus modules...")
  cons = WGCNA::blockwiseConsensusModules(
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
    verbose = verbose
  )
  
  ##
  message("[wgcna.computeConsensusLayers] Creating consensus layers...")
  aligned <- list()
  i=1
  for(i in 1:length(exprList)) {
    k <- names(exprList)[i]
    sel <- c("colors","unmergedColors","goodSamples","goodGenes",
      "dendrograms","blockGenes","blocks")
    net <- cons[sel]
    net$power <- power[i]
    X = exprList[[i]]
    w <- wgcna.compute(
      X = exprList[[i]],
      samples,
      prefix = prefix[i],
      ngenes = -1,
      net = net,
      calcMethod = "fast",
      sv.tom=0
    )
    aligned[[k]] <- w
  }

  return(aligned)
}


#' Compute consensus matrix from list of matrices. The consensus
#' matrix checks for consistent sign and minimal threshold for each
#' matrix. Optionally filters on consistent p-value.
#'
#' @param ydim original dimension of data
#'
#' 
#' @export
wgcna.computeConsensusMatrix <- function(matlist, ydim, psig = 0.05, consfun="min") {

  ## create consensus module-trait matrix
  matsign <- lapply(matlist, sign)
  matsign <- lapply(matsign, function(x) {x[is.na(x)]=0; x})
  all.pos <- Reduce("*", lapply(matsign, function(z) (z >= 0) ))
  all.neg <- Reduce("*", lapply(matsign, function(z) (z <= 0) )) 
  concordant <- (all.pos | all.neg)

  matlistN <- Reduce("+", lapply(matlist, function(x) !is.na(x)))
  matlist0 <- lapply(matlist, function(x) {x[is.na(x)]=0; x})
  
  zsign <- sign(Reduce("+", matsign)) ## mean sign??
  if(consfun=="min") { 
    pminFUN <- function(...) pmin(..., na.rm=TRUE)  
    consZ <- Reduce(pminFUN, lapply(matlist, abs)) * zsign
  } else if(consfun=="gmean") {
    ## geometric mean
    matlistG <- lapply(matlist, function(x) {x=log(abs(x));x[is.na(x)]=0;x})    
    consZ <- exp(Reduce('+', matlistG) / matlistN)
    consZ <- consZ * zsign
  } else {
    ## mean
    consZ <- Reduce('+', matlist0) / matlistN
  }
  consZ[!concordant] <- NA

  if(psig < 1) {
    if(length(ydim) == 1) ydim <- rep(ydim[1], length(matlist))
    pv <- mapply(function(z, n)
      WGCNA::corPvalueStudent(z, n), matlist, ydim, SIMPLIFY = FALSE)
    for(i in 1:length(pv)) pv[[i]][is.na(pv[[i]])] <- 0 ## missing???
    all.sig <- Reduce("*", lapply(pv, function(p) 1 * (p < psig)))
    consZ[!all.sig] <- NA
  }
  return(consZ)
}

#' Compute consensus matrix from list of matrices. The consensus
#' matrix checks for consistent sign and minimal threshold for each
#' matrix. Optionally filters on consistent p-value.
#'
#' @export
wgcna.computeDistinctMatrix <- function(matlist, ydim, psig = 0.05, min.diff=0.3,
                                        consmax = 0) {
  ## create difference module-trait matrix
  pv <- mapply(function(z, n) corPvalueStudent(z, n),
    matlist, ydim, SIMPLIFY = FALSE)
  matsign <- lapply( matlist, sign )
  Q <- matlist
  i=1
  for(i in 1:length(matlist)) {

    ## Any entry not significant is anyway invalid
    notsig <- (pv[[i]] > psig)
    Q[[i]][notsig] <- NA

    ## Any entry with too small difference with others is invalid
    refmat <- Reduce("+", matlist[-i]) / (length(matlist)-1)
    diff <- matlist[[i]] - refmat
    notdiff <- ( abs(diff) < min.diff )
    Q[[i]][notdiff] <- NA

    ## any entry that has consensus is invalid
    cons <- mapply(function(P, S) (P<0.05)*(S==matsign[[i]]),
      pv[-i], matsign[-i], SIMPLIFY=FALSE)
    cons <- (Reduce('+', cons) > consmax) ## or function
    Q[[i]][cons] <- NA
  }
  return(Q)
}


#' @export
wgcna.plotConsensusOverlapHeatmap <- function(net1, net2, 
                                              setLabels=NULL,
                                              lab.line = c(8,8),
                                              plotDendro = FALSE,
                                              setpar=TRUE) {

  if(is.null(setLabels))
    setLabels <- c("Set1","Set2")
  if(length(setLabels)==1) setLabels <- paste0(setLabels,1:2)
  
  if(plotDendro) {

    layout.matrix <- matrix( c(1,2,5,3,4,5), nrow = 3, ncol = 2)
    layout(mat = layout.matrix, heights= c(0.8,0.2,2.5), widths = c(1,1))

    WGCNA::plotDendroAndColors(
      dendro = net1$dendrograms[[1]],
      colors = net1$colors,
      dendroLabels = FALSE,
      hang = 0.03,
      addGuide = FALSE,
      guideHang = 0.05,
      #marAll = marAll,
      setLayout = FALSE,
      main = setLabels[1]
    )

    WGCNA::plotDendroAndColors(
      dendro = net2$dendrograms[[1]],
      colors = net2$colors,
      dendroLabels = FALSE,
      hang = 0.03,
      addGuide = FALSE,
      guideHang = 0.05,
      #marAll = marAll,
      setLayout = FALSE,
      main = setLabels[2]
    )
  }

  
  firstColors  <- wgcna.labels2colors(net1$colors)
  secondColors <- wgcna.labels2colors(net2$colors)
  overlap <- overlapTable(firstColors, secondColors)
  names(overlap)

  T1 = overlap$countTable
  T2 = table(firstColors, secondColors)
  T3 = table(net1$colors, net2$colors)

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
  mtext(toupper(setLabels[1]), side=1, line=lab.line[1], cex=1.1)
  mtext(toupper(setLabels[2]), side=2, line=lab.line[2], cex=1.1)

}

##=========================================================================
## PRESERVATION WGCNA
##=========================================================================

#' @export
wgcna.runPreservationWGCNA <- function(exprList, phenoData,
                                       reference = 1,
                                       add.merged=FALSE,
                                       ngenes = 2000,
                                       annot = NULL,
                                       compute.stats = TRUE,
                                       compute.enrichment = TRUE,
                                       gset.methods = c("fisher","gsetcor","xcor")
                                       ) {

  if(is.character(reference)) {
    reference <- match(reference, names(exprList))
  }
  if(reference > 0) {
    reference.name <- names(exprList)[reference]
  } else {
    reference.name <- "Consensus"
  }
  
  ## multiset WGCNA
  pres <- wgcna.runConsensusWGCNA(
    exprList,
    phenoData,
    GMT = NULL,
    annot = NULL,
    ngenes = ngenes,
    power = 12,
    minModuleSize = 20,
    minKME = 0.3, 
    mergeCutHeight = 0.15,
    deepSplit = 2,
    maxBlockSize = 9999,
    addCombined = FALSE,
    calcMethod = "fast",
    drop.ref = FALSE,
    compute.stats = TRUE,
    compute.enrichment = FALSE,
    gsea.mingenes = 10,
    gset.methods = gset.methods
  )
  
  colorList <- lapply(pres$layers, function(w) w$net$colors)
  names(colorList) <- names(pres$layers)
  exprList <- lapply( pres$layers, function(w) t(w$datExpr))
  
  if(add.merged || reference==0) {
    message("[wgcna.runPreservationWGCNA] adding merged layer...")
    cX <- lapply(exprList, function(x) x - rowMeans(x))
    merged <- do.call( cbind, cX)
    exprList$Merged <- NULL
    exprList <- c( list( Merged = merged), exprList)
    cons.colors <- pres$net$colors
    colorList <- c( list( Consensus = cons.colors ), colorList)
    reference <- reference + 1
  }
  
  message("[wgcna.runPreservationWGCNA] running WGCNA::modulePreservation...")
  multiExpr = WGCNA::list2multiData(lapply(exprList,Matrix::t))
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
  Z <- sapply(mp.tables, function(mat) mat[,"Zsummary.pres"])
  rownames(Z) <- rownames(mp.tables[[1]])
  rownames(Z) <- paste0("ME",rownames(Z))
  colnames(Z) <- names(multiExpr)[-reference]

  ## module size
  moduleSize <- mp.tables[[1]][,"moduleSize"]
  names(moduleSize) <- rownames(Z)
  
  ## median rank
  mp.tables <- mp$preservation$observed[[1]][-reference]
  M <- sapply(mp.tables, function(mat) mat[,"medianRank.pres"])
  rownames(M) <- rownames(mp.tables[[1]])
  rownames(M) <- paste0("ME",rownames(M))
  colnames(M) <- names(multiExpr)[-reference]

  ## module-traits. We need to recompute the MEs (module eigengenes)
  ## using the color coding of the reference set.
  refColors <- colorList[[1]]
  MEx = lapply(exprList, function(x)
    WGCNA::moduleEigengenes(t(x), colors = refColors)$eigengenes
  )
  names(MEx)

  ## Compute module-trait correlation matrices
  Y <- lapply(pres$layers, function(w) w$datTraits)
  names(Y)
  if("Merged" %in% names(MEx) && !"Merged" %in% names(Y) ) {
    kk <- rownames(MEx[['Merged']])
    Y[['Merged']] <- pres$datTraits[kk,]
    Y <-Y[names(MEx)]
  }
  kk <- Reduce( union, lapply(Y, colnames))
  Y <- lapply(Y, function(y) y[,match(kk,colnames(y))])
  for(i in 1:length(Y)) colnames(Y[[i]]) <- kk
  R <- mapply(cor, MEx, Y, use="pairwise", SIMPLIFY=FALSE)
  for(i in 1:length(R)) {
    colnames(R[[i]]) <- paste0(names(R)[i],":",colnames(R[[i]]))
  }
  
  ## gene statistics of reference layer
  if(compute.stats) {
    message("[wgcna.runPreservationWGCNA] computing gene statistics...")
    ref = reference.name
    wnet <- list( MEs = MEx[[ref]], colors = pres$colors[,ref])
    pres$stats <- wgcna.computeGeneStats(wnet, pres$datExpr[[ref]],
      pres$datTraits, TOM=NULL)
  }
  
  ## geneset enrichment of reference layer
  if(compute.enrichment) {
    message("[wgcna.runPreservationWGCNA] computing geneset enrichment...")
    pres$gsea <- wgcna.computeModuleEnrichment(
      pres$layers[[ref]],
      GMT = NULL,
      gsetX = NULL,
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

#' @export
wgcna.plotPreservationSummaries <- function(pres, setpar=TRUE) {

  # Create a simple bar plot of Zsummary:
  Z <- pres$Zsummary
  ntest <- ncol(Z)

  if(setpar) {
    par(mfrow=c(3, ntest), mar=c(5,5,4,1))
  }
  xylist <- list(
    c("moduleSize", "Zsummary.pres"),
    c("moduleSize", "medianRank.pres"),
    c("Zsummary.pres", "medianRank.pres")  
  )

  for(xy in xylist) {
    for(k in colnames(Z)) {
      X <- data.frame(
        color = substring(names(pres$moduleSize),3,99),
        moduleSize = pres$moduleSize,
        Zsummary.pres = pres$Zsummary[,k],
        medianRank.pres = pres$medianRank[,k]
      )
      xvar <- xy[1]
      yvar <- xy[2]
      ylim <- c(0,max(X[,yvar]))
      if(yvar == "medianRank.pres") ylim <- rev(ylim)
      plot(
        X[,xvar],
        X[,yvar],
        pch = 21,
        cex=2,
        bg = X$color,
        ylim = ylim,
        xlab = xvar,
        ylab = yvar
      )
      title(yvar, cex.main=1.4, line=2.2)  
      sub <- paste( k, "vs.", "reference")
      title(sub, cex.main=1, line=0.9)  
      if(yvar=="Zsummary.pres") abline(h = c(2,10), lty = 2)
    }
  }
}

#' @export
wgcna.plotPreservationModuleTraits <- function(pres, subplots=FALSE,
                                               setpar=TRUE, rm.na=FALSE) {  

  if(setpar) {
    par(mfrow=c(2,2), mar=c(14,12,4,2))
    if(subplots) {
      par(mfrow=c(3,3), mar=c(14,12,4,2))
    }
  }
  
  ##--------------------------------------
  ## Zsummary heatmap
  ##--------------------------------------

  Z <- pres$Zsummary
  WGCNA::labeledHeatmap(
    Matrix = Z,
    xLabels = colnames(Z),
    yLabels = rownames(Z),
    ySymbols = rownames(Z),
    colors = tail(WGCNA::blueWhiteRed(100),50),
    colorLabels = TRUE,
    setStdMargins = FALSE
  )
  title("Module preservation (Zsummary)", line=1.2, cex.main=1.2)
  
  ##--------------------------------------
  ## Merged Module-Trait
  ##--------------------------------------

  ## merge module-traits
  R <- pres$modTraits
  sel <- unique(unlist(lapply(R, function(x) which(colMeans(is.na(x))<1))))
  R <- lapply(R, function(r) r[,sel,drop=FALSE])
  R <- do.call(cbind, R)
  
  ii <- intersect(rownames(Z),rownames(R))
  R <- R[ii,,drop=FALSE]
  Z <- Z[ii,,drop=FALSE]  
  ##R[is.na(R)] <- 0

  ## collapse by consensus
  group <- sub(":.*","",colnames(R))
  cR <- tapply( 1:ncol(R), group, function(i) R[,i]) 
  for(i in 1:length(cR)) colnames(cR[[i]]) <- sub(".*:","",colnames(cR[[i]]))
  cR <- lapply(cR, function(x) x[,sort(colnames(cR[[1]]))])
    
  validcol <- function(R) which( colMeans(is.na(R)) < 1 &
                                   matrixStats::colSds(R,na.rm=TRUE) > 0.01 )
  ii <- 1:ncol(R)
  if(rm.na) ii <- validcol(R)
  
  WGCNA::labeledHeatmap(
    Matrix = R[,ii,drop=FALSE],
    xLabels = colnames(R)[ii],
    yLabels = rownames(R),
    ySymbols = rownames(R),
    colors = WGCNA::blueWhiteRed(100),
    colorLabels = TRUE,
    setStdMargins = FALSE
  )
  title("Module-Trait correlation", line=1.2, cex.main=1.2)

  if(subplots) {
    for(i in 1:length(cR)) {
      ## set color scale
      R1 <- cR[[i]]
      clim <- max(abs(R1),na.rm=TRUE)
      cval <- seq( -clim, clim, length.out=201)
      ii <- which( cval >= min(R1,na.rm=TRUE) & cval <= max(R1,na.rm=TRUE) )
      col1 <- WGCNA::blueWhiteRed(201)[ii]
      jj <- 1:ncol(R1)
      if(rm.na) jj <- validcol(R1)
      WGCNA::labeledHeatmap(
        Matrix = R1[,jj,drop=FALSE],
        xLabels = colnames(R1)[jj],
        yLabels = rownames(R1),
        ySymbols = rownames(R1),
        colors = col1,
        colorLabels = TRUE,
        setStdMargins = FALSE
      )
      title(names(cR)[i], line=1.2, cex.main=1.2)
    }   
  }

  ##--------------------------------------
  ## Consensus Module-Trait
  ##--------------------------------------
  ydim <- sapply(pres$layers, function(w) nrow(w$datTraits))
  consZ <- wgcna.computeConsensusMatrix(cR, ydim, psig=1, consfun='gmean')
  consZ <- consZ[rownames(cR[[1]]), colnames(cR[[1]])]

  clim <- max(abs(consZ),na.rm=TRUE)
  cval <- seq( -clim, clim, length.out=201)
  ii <- which( cval >= min(consZ,na.rm=TRUE) & cval <= max(consZ,na.rm=TRUE) )
  col2 <- WGCNA::blueWhiteRed(201)[ii]

  jj <- 1:ncol(consZ)
  if(rm.na) jj <- validcol(consZ)

  WGCNA::labeledHeatmap(
    Matrix = consZ[,jj,drop=FALSE],
    xLabels = colnames(consZ)[jj],
    yLabels = rownames(consZ),
    ySymbols = rownames(consZ),
    colors = col2,
    colorLabels = TRUE,
    setStdMargins = FALSE
  )
  title("Consensus Module-Traits", line=1.2, cex.main=1.2)

  ##--------------------------------------
  ## preservation-weighted COnsensus Module-Trait
  ##--------------------------------------

  wz <- rowMeans(Z**2, na.rm=TRUE)
  wz <- wz / max(wz)
  consW <- consZ * wz

  clim <- max(abs(consW),na.rm=TRUE)
  cval <- seq( -clim, clim, length.out=201)
  ii <- which( cval >= min(consW,na.rm=TRUE) & cval <= max(consW,na.rm=TRUE) )
  col2 <- WGCNA::blueWhiteRed(201)[ii]

  jj <- 1:ncol(consW)
  if(rm.na) jj <- validcol(consW)
  
  WGCNA::labeledHeatmap(
    Matrix = consW[,jj,drop=FALSE],
    xLabels = colnames(consW)[jj],
    yLabels = rownames(consW),
    ySymbols = rownames(consW),
    colors = col2,
    colorLabels = TRUE,
    setStdMargins = FALSE
  )
  title("Preservation-weighted Consensus\nModule-Traits", line=1, cex.main=1.2)
}



##=========================================================================
## PLOTTING FUNCTIONS
##=========================================================================


#' @export  
wgcna.plotTopModules <- function(wgcna, trait, nmax=16, setpar=TRUE) {

  MEx <- wgcna$net$MEs
  Y <- wgcna$datTrait
  kk <- intersect(rownames(MEx), rownames(Y))
  MEx <- MEx[kk,]
  Y <- Y[kk,,drop=FALSE]
  
  rho <- cor(MEx, Y, use="pairwise")
  sel <- order(-abs(rho[,trait]))
  sel <- head(sel, nmax)
  n <- length(sel)
  nr <- ceiling(sqrt(n))
  nc <- ceiling(n / nr)

  yclass <- sapply(as.data.frame(Y), class)
  is.binary <- apply(Y, 2, function(x) all(x %in% c(TRUE,FALSE,0,1,NA)))
  yclass[which(is.binary)] <- "factor"
  
  if(setpar==1) par(mfrow=c(nr,nc), mgp=c(2.6,0.85,0), mar=c(4,4,2.5,1))
  if(setpar==2) par(mfrow=c(nc,nr), mgp=c(2.6,0.85,0), mar=c(4,4,2.5,1))  

  i=1
  for(i in head(sel, nmax)) {
    x <- Y[, trait]  
    y <- MEx[,i]
    label <- colnames(MEx)[i]
    col <- substring(label, 3, 99)
    col1 <- adjustcolor(col, alpha.f=0.5)

    if( yclass[trait] == "factor" ) {
      boxplot(y ~ x, main = label, col = col1,
        xlab = trait, ylab = "ME score")
      points(1 + x + 0.04*rnorm(length(x)), y,
        pch=21, bg=col1, lwd=0.5)
    }

    if( yclass[trait] == "numeric" ) {
      plot( x, y, main = label,
        pch=21, cex = 1.1, col=1, bg = col1, lwd=0.25,
        xlab = trait, ylab = "ME score")
      abline(h=0, lty=2, lwd=0.5)
      r <- cor(x,y,use="pairwise")
      if(abs(r) > 0.3) {
        abline( lm(y ~ x), col=1, lwd=0.6 )
        legend("bottomright", legend=paste("r=",round(r,3)))
      }
    }

  }

}


#' Plot top modules most correlated with trait for multi expression
#' data.
#'
#' @export  
wgcna.plotTopModules_multi <- function(multi, trait, nmax=16,
                                       plotlib = "base", setpar=TRUE)
{

  ## we 'just' concatenate the ME matrices
  MEx <- do.call(rbind, lapply(multi$MEs, as.matrix))
  me.samples <- lapply(multi$MEs, rownames)
  Y <- multi$datTrait
  batch <- max.col(sapply(me.samples, function(s) rownames(Y) %in% s))
  batch <- names(multi$MEs)[batch]
  names(batch) <- rownames(Y)

  kk <- intersect(rownames(MEx), rownames(Y))
  MEx <- MEx[kk,]
  Y <- Y[kk,]
  batch <- batch[kk]
  nbatch <- length(multi$MEs)
  
  rho <- cor(MEx, Y, use="pairwise")
  sel <- order(-abs(rho[,trait]))
  sel <- head(sel, nmax)
  n <- length(sel)
  nr <- ceiling(sqrt(n))
  nc <- ceiling(n / nr)

  if(setpar==1) par(mfrow=c(nr,nc), mgp=c(2.6,0.85,0), mar=c(4,4,2.5,1))
  if(setpar==2) par(mfrow=c(nc,nr), mgp=c(2.6,0.85,0), mar=c(4,4,2.5,1))  

  yclass <- sapply(as.data.frame(Y), class)
  is.binary <- apply(Y, 2, function(x) all(x %in% c(TRUE,FALSE,0,1,NA)))
  yclass[which(is.binary)] <- "logical"
  yclass
  
  i=sel[1]
  for(i in head(sel, nmax)) {

    x <- Y[, trait]  
    y <- MEx[,i]
    label <- colnames(MEx)[i]
    col <- substring(label, 3, 99)
    col1 <- adjustcolor(col, alpha.f=0.66)
    col2 <- col1

    yclass[trait]
    if(yclass[trait] %in% c("factor","logical","binary")) {

      if(yclass[i] %in% c("logical","binary")) {
        x <- (x == 1)
      }
      df <- data.frame(x=x, y=y, group=factor(batch))
    
      if(plotlib=="base") {
        par(mgp=c(2.4,0.9,0))
        
        aa <- seq(0.45, 0.1, length.out=nbatch)
        col1 <- sapply(aa, function(a) adjustcolor(col, alpha.f=a))
        
        nb <- nbatch
        atx <- c(1,2, 4,5)
        atx <- c( seq(1,2, length.out=nb), seq(4,5, length.out=nb) )
        atmid <- c(2, 4)
        
        boxplot(
          #df$y ~ df$x + df$group,
          df$y ~ df$group + df$x,        
          #df$y ~ df$x,        
          at = atx,
          cols = col1,
          main = label,
          col = col1,
          boxwex=0.8,
          xaxt = 'n',
          xlab = trait,
          ylab = "ME score"
        )
        mtext( c("FALSE","TRUE"), side=1, line=0.75, cex=0.75, at=atmid)
        bb <- names(multi$MEs)
        legend("topright", legend=bb, fill=col1, cex=0.9,
          y.intersp=0.85)        
      }
      
      # Boxplot with default spaces
      if(plotlib=="ggplot") {
        ggplot2::ggplot(df,                          
          ggplot2::aes(x = x, y = y, fill = group)) +
          ggplot2::geom_boxplot() +
          ggplot2::xlab(trait) +
          ggplot2::ylab("ME score")       
      }
    }  ## end of if factor

    if(yclass[trait] %in% c("numeric","integer")) {
      
      df <- data.frame(x=x, y=y, group=factor(batch))
    
      if(plotlib=="base") {
        par(mgp=c(2.4,0.9,0))
        
        aa <- seq(0.55, 0.15, length.out=nbatch)
        col1 <- sapply(aa, function(a) adjustcolor(col, alpha.f=a))
        names(col1) <- sort(unique(batch))
        
        nb <- nbatch
        colx <- col1[as.integer(factor(batch))]
        
        plot(
          df$x , df$y,
          pch = 21,
          cex = 1,
          lwd = 0.4,
          bg = colx,
          main = label,
          xlab = trait,
          ylab = "ME score"
        )
        bb <- names(multi$MEs)
        legend("topright", legend=bb, fill=col1, cex=0.9,
          y.intersp=0.85)

        ## add regression lines
        b=df$group[1]
        col2 <- adjustcolor(col1, red.f=0.5, green.f=0.5, blue.f=0.5)
        names(col2) <- names(col1)
        for(b in unique(df$group)) {
          ii <- which(df$group == b)
          abline( lm(df$y[ii] ~ df$x[ii]), lwd=1, lty=1, col=col2[1] )
        }
        
      }
      
      # Boxplot with default spaces
      if(plotlib=="ggplot") {
        ggplot2::ggplot(df,                          
          ggplot2::aes(x = x, y = y, fill = group)) +
          ggplot2::geom_point() +
          ggplot2::xlab(trait) +
          ggplot2::ylab("ME score")       
      }
    }  ## end of if continuous

  }

    

  
}


#'
#'
#' @export
wgcna.plotTOM <- function(wgcna, justdata = FALSE, block = NULL,
                          legend = TRUE, downsample = NULL) {
  datExpr <- wgcna$datExpr
  wTOM <- NULL
  if (!is.null(wgcna$TOM)) {
    wTOM <- wgcna$TOM
  }
  ## if SV of TOM is stored, reconstruct TOM
  if (is.null(wTOM) && !is.null(wgcna$svTOM)) {
    wTOM <- wgcna$svTOM %*% t(wgcna$svTOM)
  }
  if (is.null(wTOM)) {
    message("[wgcna.plotTOM] ERROR. no TOM matrix")
    return(NULL)
  }
  dissTOM <- 1 - wTOM
  rownames(dissTOM) <- colnames(dissTOM) <- colnames(datExpr)

  ## clustering wgcnaults
  moduleColors <- NULL
  if (is.null(block) && "merged_dendro" %in% names(wgcna$net)) {
    geneTree <- wgcna$net$merged_dendro
    gg <- geneTree$labels
    if (length(gg) > 0) {
      dissTOM <- dissTOM[gg, gg]
      moduleColors <- wgcna.labels2colors(wgcna$net$colors[gg])
    }
  }

  if (is.null(moduleColors)) {
    if (is.null(block)) block <- 1
    geneTree <- wgcna$net$dendrograms[[block]]
    ii <- which(wgcna$net$blocks == block & wgcna$net$goodGenes == TRUE)
    gg <- names(wgcna$net$color)[ii]
    dissTOM <- dissTOM[gg, gg]
    moduleColors <- wgcna.labels2colors(wgcna$net$colors[gg])
  }

  if (justdata) {
    return(dissTOM)
  }

  if (!is.null(downsample) && ncol(dissTOM) > downsample) {
    ii <- seq(1, ncol(dissTOM), length.out = downsample)
    dissTOM <- dissTOM[ii, ii]
    moduleColors <- moduleColors[ii]
    geneTree <- fastcluster::hclust(as.dist(dissTOM),
      method = "average"
    )
  }

  if (legend) {
    par(oma = c(1, 0, 0, 0), mar = c(0, 0, 0, 0))
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
                                      marAll = c(0.2, 8, 2, 0.2),
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
  
  if(is.null(main)) main <- "Gene dendrogram and module colors"
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

#'
#'
#' @export
wgcna.plotDendroAndTraitCorrelation <- function(wgcna,
                                                traits = NULL,
                                                show.traits = TRUE,
                                                main=NULL,
                                                block=NULL,
                                                rm.na = TRUE,
                                                use.tree = 0,
                                                marAll = c(0.2, 8, 2, 0.2),
                                                setLayout=TRUE,
                                                ... )
{

  ## if consensus output do this
  is.cons <- ("class" %in% names(wgcna) && wgcna$class == "cons")
  is.cons2 <- (all(c("layers","zlist") %in% names(wgcna)))
  if(is.cons || is.cons2) {
    wgcna.plotDendroAndTraitCorrelation_cons(
      cons = wgcna,
      traits = traits,
      main = main,
      rm.na = rm.na,
      show.traits = show.traits,
      marAll = marAll,
      use.tree = use.tree,
      setLayout = setLayout,
      ... )
    return()
  }

  moduleColors <- cbind(wgcna$net$colors)
  if(NCOL(moduleColors)==1) colnames(moduleColors) <- "Module"
  colors <- moduleColors

  if(show.traits) {
    X <- wgcna$datExpr
    Y <- wgcna$datTraits
    traitSig <- cor(X, Y, use="pairwise")
    if(rm.na) {
      sel <- colMeans(is.na(traitSig)) < 1
      traitSig <- traitSig[, sel, drop=FALSE]
    }
    traitColors <- rho2bluered(traitSig)
    colors <- cbind(moduleColors, 0, traitColors)
  }
  
  geneTree <- wgcna$net$dendrograms[[1]]
  geneTree$labels <- names(wgcna$net$colors)
  colors <- colors[geneTree$labels, ]
  
  if(is.null(main)) main = "Gene Dendrogram, Modules and Trait Correlation"

  WGCNA::plotDendroAndColors(
    geneTree,
    colors = colors,
    colnames(colors),
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05,
    marAll = marAll,
    main = main,
    setLayout = setLayout,
    ...
  )

}

#' wgcna.plotDendroAndTraits for Consensus output
#'
#' 
wgcna.plotDendroAndTraitCorrelation_cons <- function(cons,
                                                     show.traits = TRUE,
                                                     traits = NULL,
                                                     main = NULL,
                                                     rm.na = TRUE,
                                                     use.tree = 0,
                                                     marAll = c(0.2, 8, 2, 0.2),
                                                     setLayout=TRUE,
                                                     ... )
{

  if(0) {
    show.traits = TRUE
    traits = NULL
    main = NULL
    rm.na = TRUE
    use.tree = 0
    marAll = c(0.2, 8, 2, 0.2)
    setLayout=TRUE
  }

  ## quick hack to use wgcna.plotDendroAndTraitCorrelation_multi()
  multi <- c( list(Consensus = cons), cons$layers )
  use.tree0 <- use.tree
  if(use.tree %in% 0:99) use.tree <- as.integer(use.tree)
  if(is.character(use.tree)) {
    use.tree <- match(use.tree, names(multi))
  } else {
    use.tree <- as.integer(use.tree) + 1
  }
  if(is.na(use.tree)) {
    message("ERROR: invalid class(use.tree) = ", class(use.tree0))
    message("ERROR: invalid use.tree = ", use.tree0)    
    return(NULL)
  }

  wgcna.plotDendroAndTraitCorrelation_multi(
    multi,
    show.traits = show.traits,
    traits = traits,
    main = main,
    rm.na = rm.na,
    use.tree = use.tree + 1,
    marAll = marAll,
    setLayout = setLayout,
    ...
  )
  
}


#' @export
wgcna.plotDendroAndTraitCorrelation_multi <- function(multi,
                                                      show.traits = TRUE,
                                                      traits = NULL,
                                                      main = NULL,
                                                      rm.na = TRUE,
                                                      use.tree = 1,
                                                      marAll = c(0.2, 8, 2, 0.2),
                                                      setLayout=TRUE,
                                                      ... )
{

  ## module colors  
  colors <- sapply( multi, function(m) m$net$colors )
  
  if(show.traits) {

    traitSig <- list()
    nsets <- length(multi)
    i=1
    for(k in names(multi)) {
      if( k == "Consensus") next
      w <- multi[[k]]
      Y <- w$datTraits
      sel <- colnames(Y)
      if(!is.null(traits)) sel <- intersect(sel, traits)
      X <- w$datExpr
      kk <- intersect(rownames(X), rownames(Y))
      traitSig[[k]] <- cor(X[kk,], Y[kk,sel], use="pairwise")
    }
    
    if(rm.na) {
      for(i in 1:length(traitSig)) {
        sel <- colMeans(is.na(traitSig[[i]])) < 1
        traitSig[[i]] <- traitSig[[i]][, sel, drop=FALSE]
      }
    }
    
    ## prepend datatype/set name
    for(k in names(traitSig)) {
      colnames(traitSig[[k]]) <- paste0(k,":",colnames(traitSig[[k]]))
    }

    traitSig2 <- c()
    for(i in 1:length(traitSig)) {
      traitSig2 <- cbind(traitSig2, traitSig[[i]])
      if(i < length(traitSig)) traitSig2 <- cbind(traitSig2, 0)
    }

    traitColors <- rho2bluered(traitSig2, f=0.95)
    ii <- which(colnames(traitColors)=='')
    if(length(ii)) traitColors[,ii] <- "#FFFFFF"
    if(is.null(colors)) {
      colors <- traitColors
    } else {
      colors <- cbind( colors, 0, traitColors)
    }
    
  }

  message("using tree of layer: ", names(multi)[use.tree] )
  geneTree <- multi[[use.tree]]$net$dendrograms[[1]]
  
  if(is.null(main)) main = "Gene Dendrogram, Modules and Trait Correlation"
  
  WGCNA::plotDendroAndColors(
    geneTree,
    colors = colors,
    colnames(colors),
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05,
    marAll = marAll,
    main = main,
    ...
  )

}


#' Converts correlation values [-1;1] to blue-white-red colors. Good
#' for creating color labels for labeledHeatmaps that expect colors.
#' NOTE: use WGCNA::numbers2colors???
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
    TSPvalue <- wgcna$stats$TSPvalue
  } else {
    traitSignificance <- as.data.frame(cor(wgcna$datExpr, wgcna$datTraits, use = "p"))
    TSPvalue <- as.data.frame(corPvalueStudent(as.matrix(traitSignificance), nSamples))
  }

  x <- (moduleMembership[moduleGenes, module])
  y <- (traitSignificance[moduleGenes, trait])
  if (abs == TRUE) {
    x <- abs(x)
    y <- abs(y)
  }
  ##
  px <- MMPvalue[moduleGenes, module]
  py <- TSPvalue[moduleGenes, trait]
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
                                         nmax = -1, tmax = -1,
                                         text = TRUE, pstar = TRUE) {
  
  if(!multi) wgcna <- list(wgcna)
  nSamples <- nrow(wgcna[[1]]$datExpr)
  MEs <- do.call(cbind, lapply(wgcna, function(w) as.matrix(w$net$MEs))) 
  Y <- wgcna[[1]]$datTraits
  
  moduleTraitCor <- cor(MEs, Y, use = "pairwise.complete")
  if(nmax > 0) {
    sel <- head(order(-apply(abs(moduleTraitCor), 1, max, na.rm=TRUE)),nmax)
    moduleTraitCor <- moduleTraitCor[sel,,drop=FALSE]
  }
  if(tmax > 0) {
    sel <- head(order(-apply(abs(moduleTraitCor), 2, max, na.rm=TRUE)),tmax)
    moduleTraitCor <- moduleTraitCor[,sel,drop=FALSE]
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
wgcna.plotEigenGeneClusterDendrogram <- function(wgcna = NULL,
                                                 ME = NULL,
                                                 add_traits = TRUE,
                                                 horiz = FALSE,
                                                 setMargins = TRUE,
                                                 method = 'wgcna',
                                                 showlabels = TRUE,
                                                 plot = TRUE,
                                                 multi = FALSE,
                                                 main = NULL) {
  # Matrix with eigengenes and traits
  if(is.null(wgcna) && is.null(ME)) {
    stop("ERROR: wgcna or ME must be given")
  }
  
  if(is.null(ME)) {

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
  }

  ME <- WGCNA::orderMEs(ME)
  if (NCOL(ME) <= 2) ME <- cbind(ME, ME) ## error if ncol(ME)<=2 !!!!
  if (is.null(main)) main <- "Eigengene Dendrogram"

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
      plot( as.dendrogram(hc), horiz = horiz, main = main)
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
wgcna.plotEigenGeneAdjacencyHeatmap <- function(wgcna,
                                                add_traits = TRUE,
                                                traits = NULL,
                                                add_me = TRUE,
                                                marx = 1, main = NULL,
                                                multi = FALSE, 
                                                phenotype = NULL,
                                                colorlabel = TRUE,
                                                text = FALSE,
                                                pstar = TRUE,
                                                setMargins = TRUE,
                                                mar1 = c(5.5, 5, 1.6, 1),
                                                mar2 = c(8, 10, 4, 2),
                                                cex.lab = 0.8,
                                                cex.text = 0.7,
                                                plotDendro = TRUE,
                                                plotHeatmap = TRUE, 
                                                dendro.horiz = TRUE,
                                                dendro.width = 0.3,
                                                dendro.labels = TRUE,
                                                nmax = -1,
                                                fixclust = FALSE,
                                                mask.intra = FALSE,
                                                justdata = FALSE) {

  if(0) {
    add_traits = TRUE;
    traits = NULL;
    marx = 1; main = NULL;
    multi = FALSE; 
    phenotype = NULL;
    colorlabel = TRUE;
    text = FALSE;
    pstar = TRUE;
    setMargins = TRUE;
    mar1 = c(5.5, 5, 1.6, 1);
    mar2 = c(8, 10, 4, 2);
    cex.lab = 0.8;
    cex.text = 0.7;
    plotDendro = TRUE;
    plotHeatmap = TRUE; 
    dendro.horiz = TRUE;
    dendro.width = 0.3;
    dendro.labels = TRUE;
    nmax = -1;
    fixclust = FALSE;
    mask.intra = FALSE;
    justdata = FALSE
    add_me = TRUE
  }

  
  if(!multi) wgcna <- list(gx=wgcna)

  # Matrix with eigengenes and traits
  ME <- c()
  if(add_me) {
    ME <- lapply(wgcna, function(w) as.matrix(w$net$MEs))
    ME <- do.call(cbind, ME)
  }
  Y <- wgcna[[1]]$datTraits    
  
  if (add_traits) {
    sel <- colnames(Y)
    if(!is.null(traits)) {
      sel <- intersect(traits, sel)
    }
    ME <- cbind(ME, Y[,sel,drop=FALSE])
  } 
  
  if (!add_traits && !is.null(phenotype)) {
    ME <- cbind(ME, Y[,phenotype,drop=FALSE])
  } 
  
  if (NCOL(ME) <= 2) ME <- cbind(ME, ME) ## error if ncol(ME)<=2 !!!!

  ## eigengene correlation
  R <- cor(ME, use="pairwise")
  R0 <- R

  ## If phenotype is given we condition the heatmap using the
  ## correlation to the phenotype.
  if(!is.null(phenotype)) {
    ## proper sign in case of inhibitor layer (like miRNA)
    layersign <- rep(1, length(wgcna))
    names(layersign) <- names(wgcna)
    layersign[grep("^mi",names(wgcna),ignore.case=TRUE)] <- -1
    ff <- list()
    for(k in names(wgcna)) {
      rho <- cor(ME, Y[,phenotype], use="pairwise")[,1]
      ##rho <- wgcna[[k]]$modTraits[,phenotype]      
      ff[[k]] <- layersign[k] * rho
    }
    names(ff) <- NULL
    ff <- unlist(ff)
    ff <- 0.5 * (1 + ff)  ## signed...
    ff <- ff[match(rownames(R),names(ff))]
    names(ff) <- rownames(R)
    ff[is.na(ff)] <- 1  ## really??? NEED RETHINK
    ww <- outer(ff, ff)            
    ##ww[is.na(ww)] <- 0
    ww <- ww / max(ww,na.rm=TRUE)
    R <- R * ww
  }
  
  if(nmax>0) {
    if(!is.null(phenotype)) {
      rho <- cor(ME, Y[,phenotype], use="pairwise")[,1]
      ii <- head(order(-abs(rho)), nmax)
    } else {
      ii <- head( order(-Matrix::rowMeans(R**2)), nmax)
    }
    R <- R[ii,ii]
  }
  
  if (justdata) {
    return(R)
  }
  
  # Plot the correlation heatmap matrix (note: this plot will overwrite
  # the dendrogram plot)
  if (is.null(main)) main <- "Eigengene Adjacency Heatmap"
  
  if(plotDendro && plotHeatmap) {
    layout.matrix <- matrix( 1:2, nrow = 1, ncol = 2)
    layout(mat = layout.matrix, heights = 1, widths = c(dendro.width, 1))
    if(dendro.horiz && dendro.labels) {
      mar1[4] <- mar2[2] ## copy left margin
    }
  }
  if(plotDendro) par(mar=mar1)

  #fixclust=FALSE
  if(fixclust) {
    ii <- rownames(R)
    hc <- hclust(as.dist(1 - R0[ii,ii]), method="average")    
  } else {
    hc <- hclust(as.dist(1 - R), method="average")    
  }
  if(plotDendro) {
    plot( as.dendrogram(hc), horiz = TRUE, ylab="Eigengene dendrogram")
  }
  
  if(plotHeatmap) {
    ii <- hc$labels[hc$order]
    ii <- intersect(ii, rownames(R))
    R1 <- R[rev(ii), ii]
    nsamples <- nrow(Y)
    par(mar=mar2)
    wgcna.plotLabeledCorrelationHeatmap(
      R1,
      nSamples = nsamples,
      text = text,
      pstar = pstar,
      colorlabel = colorlabel,
      cluster = FALSE,
      setpar = FALSE,
      main = main,
      cex.lab = cex.lab,
      cex.text = cex.text
    )
  }  
}

#' @export
wgcna.plotMultiEigengeneCorrelation <- function(wgcna, addtraits = TRUE,
                                                phenotype=NULL, nmax=-1, main = NULL,
                                                showvalues = FALSE, showsig = TRUE,
                                                cex.text = 0.7, cex.lab = 0.8,
                                                fixcluster = TRUE, setpar = TRUE) {

  ## Show inter-correlation of modules
  me <- lapply(wgcna, function(w) w$net$MEs)        
  if(length(me) == 1) {
    me <- list(me[[1]], me[[1]])
  }
  
  comb <- combn(length(me),2)
  ncomb <- ncol(comb)
  nsamples <- nrow(wgcna[[1]]$datExpr)
  Y <- wgcna[[1]]$datTraits
  
  ## for miRNA we need to flip sign
  msign <- c(1,-1)[1 + 1*(names(wgcna) %in% c("mi","mir"))]

  if(setpar) {
    nc <- ceiling(sqrt(ncomb))
    nr <- ceiling(ncomb / nc)
    par(mfrow=c(nr,nc), mar=c(8,10,3,1))
  }

  k=1
  for(k in 1:ncol(comb)) {
    i <- comb[1,k]
    j <- comb[2,k]
    M1 <- me[[i]]
    M2 <- me[[j]]

    if(addtraits) {
      M1 <- cbind(M1, Y)
      M2 <- cbind(M2, Y)
    }
    if(FALSE && !addtraits && !is.null(phenotype)) {
      y <- Y[,phenotype,drop=FALSE]
      M1 <- cbind(M1, y)
      M2 <- cbind(M2, y)
    }


    R1 <- cor(M1, M2, use="pairwise.complete")

    if(nmax > 0) {
      ii <- head(order(-apply(abs(R1), 1, max)), nmax)
      jj <- head(order(-apply(abs(R1), 2, max)), nmax)          
      R1 <- R1[ii,jj]
    }
    
    ## cluster unweighted matrix
    ii <- hclust(dist(R1), method="average")$order
    jj <- hclust(dist(t(R1)), method="average")$order
    R1 <- R1[ii, jj]
    
    ## This conditions the correlation on phenotype. Important.
    do.condition <- !is.null(phenotype)
    if(do.condition) {
      y <- Y[,phenotype]
      w1 <- cor( M1[,rownames(R1)], y, use="pairwise")[,1]
      w2 <- cor( M2[,colnames(R1)], y, use="pairwise")[,1]      
      if(msign[i]!=0) w1 <- msign[i] * w1
      if(msign[j]!=0) w2 <- msign[j] * w2            
      w1 <- pmax(w1,0)
      w2 <- pmax(w2,0)            
      ww <- outer(w1, w2)            
      ww <- ww / max(ww, na.rm=TRUE)
      R1 <- R1 * ww
    }
    
    main <- paste(names(me)[i],"vs.",names(me)[j])
    if(do.condition) main <- paste(main, "(conditioned)")
    
    playbase::wgcna.plotLabeledCorrelationHeatmap(
      R1,
      nsamples,
      text = showvalues,
      pstar = showsig,
      is.dist = FALSE,
      cluster = !fixcluster,
      setpar = FALSE,
      main = main,
      cex.text = cex.text,
      cex.lab = cex.lab
    )
    
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

  sdx <- matrixStats::colSds(as.matrix(ME*1), na.rm = TRUE)
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
wgcna.plotConsensusSampleDendroAndColors <- function(cons, i,
                                                     what = c("both","me","traits")[1],
                                                     clust.expr = TRUE,
                                                     setLayout = TRUE,
                                                     marAll = c(0.2, 7, 1.5, 0.5),
                                                     colorHeightMax = 0.6,
                                                     main = NULL)
{

  wgcna.plotSampleDendroAndColors(
    wgcna = cons$layers[[i]],
    main = main,
    datExpr = cons$datExpr[[i]],
    datTraits = cons$datTraits,
    datME = cons$net$multiME[[i]]$data,
    what = what,
    marAll = marAll,
    clust.expr = clust.expr,
    setLayout = setLayout,
    colorHeightMax = colorHeightMax
  ) 
  
}

#'
#'
#' @export
wgcna.plotSampleDendroAndColors <- function(wgcna, input.type="wgcna",
                                            what = c("me", "traits", "both")[3],
                                            datTraits = NULL, datExpr = NULL, datME = NULL, 
                                            clust.expr = TRUE, setLayout = TRUE,
                                            marAll = c(0.2, 7, 1.5, 0.5),
                                            colorHeightMax = 0.6,
                                            main = NULL, justdata = FALSE) {
  
  if(input.type == 'net') {
    ME0 <- wgcna$MEs
    if(is.null(datExpr)) stop("must supply datExpr")
    if(is.null(datTraits)) stop("must supply datTraits")    
  } else {
    ME0 <- wgcna$net$MEs
    datTraits <- 1*wgcna$datTraits
    datExpr <- wgcna$datExpr
  }
  
  if(!is.null(datME)) {
    ME0 <- datME
  }
  
  ME <- ME0[, 0]
  samples <- rownames(ME)
  if (any(what %in% c("me", "both"))) {
    ME <- cbind(ME, ME0)
  }
  
  if (any(what %in% c("traits", "both"))) {
    ME <- cbind(ME, datTraits[samples,,drop=FALSE])
  }
  
  if (NCOL(ME) <= 2) ME <- cbind(ME, ME) ## error if ncol(ME)<=2 !!!!  
  sdx <- matrixStats::colSds(as.matrix(ME*1), na.rm = TRUE)
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
    main = main,
    colorHeightMax = colorHeightMax
  )
}


#' @export
wgcna.plotLabeledCorrelationHeatmap <- function(R, nSamples, 
                                                cluster = TRUE, text = TRUE,
                                                main = NULL, justdata = FALSE,
                                                colorlabel = TRUE, pstar = TRUE,
                                                cex.text = 0.7, cex.lab = NULL,
                                                setpar = TRUE, is.dist = FALSE) {

  ## Define numbers of genes and samples
  if (cluster && nrow(R)>1 && ncol(R)>1) {
    R0 <- R
    R0[is.na(R0)] <- 0
    is.sym <- nrow(R) == ncol(R) && all(rownames(R)==colnames(R))
    if(is.dist) {
      ii <- hclust(as.dist(abs(R0)))$order
      jj <- ii
    } else if(is.sym) {
      ii <- hclust(dist(R0), method="average")$order
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

  ## set colorscale. make sure 0 is white if non-symmetric
  col1 <- WGCNA::blueWhiteRed(201)
  rlim <- max(abs(R),na.rm=TRUE)
  if(rlim>0) {
    rval <- seq( -rlim, rlim, length.out=201)
    ii <- which( rval >= min(R,na.rm=TRUE) & rval <= max(R,na.rm=TRUE) )
    col1 <- WGCNA::blueWhiteRed(201)[ii]
  }
  
  ## Display the correlation values within a heatmap plot
  WGCNA::labeledHeatmap(
    Matrix = R,
    #xLabels = paste0(1:ncol(R),":",colnames(R)),
    xLabels = colnames(R),
    #xLabels = paste0(" ",colnames(R)),
    yLabels = rownames(R),
    xSymbols = colnames(R),
    ySymbols = rownames(R),
    colorLabels = TRUE,
    colors = col1,
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

#' @export
wgcna.plotGeneNetwork <- function(wgcna, genes, col=NULL,
                                  edge.alpha = 0.3,
                                  rgamma = 4,
                                  edge.width = 6,
                                  alpha = 0.5,
                                  min.rho = 0.5,
                                  setpar = TRUE) {

  A <- cor(wgcna$datExpr[, genes])
  A <- A * (abs(A) > min.rho)

  gr <- igraph::graph_from_adjacency_matrix(
    A, mode = "undirected", weighted = TRUE, diag = FALSE
  )
  vcex <- matrixStats::colVars(wgcna$datExpr[, genes])
  vcex <- vcex / max(abs(vcex))
  if(is.null(col)) {
    col <- wgcna$net$color[genes]
  }
  col <- sub("black","grey40",col)
  ecol <- c("darkred","darkgreen")[ 1 + 1*(igraph::E(gr)$weight > 0)]
  table(ecol)  
  ecol <- adjustcolor(ecol, alpha.f = edge.alpha)
  ewt <- abs(igraph::E(gr)$weight)
  ewt <- (ewt / max(abs(ewt)))**rgamma  
  plot(gr,
    layout = igraph::layout_in_circle,
    edge.width = edge.width * ewt,
    edge.color = ecol,
    vertex.size = 5 + 20 * vcex,
    vertex.color = adjustcolor(col, alpha.f = alpha),
      vertex.frame.color = col
  )
}

#' @export
wgcna.plotModuleHeatmap <- function(wgcna, module, genes=NULL,
                                    rgamma = 4, min.rho = 0,
                                    cex = 0.8, nmax=-1,
                                    type = c("expression","correlation")[1],
                                    heatmap.mar = c(7,7),
                                    diag = TRUE,
                                    main = NULL) {

  if(!is.null(module) && is.null(genes) ) {
    genes <- wgcna$me.genes[[module]]
  }
  if(is.null(genes) || length(genes)==0) {
    stop("must specify genes or module")
  }

  if(nmax > 0) {
    sdx <- matrixStats::colSds(wgcna$datExpr[,genes])
    genes <- head( genes[order(-sdx)], nmax)
  }

  if(type == "expression") {
    X <- t(wgcna$datExpr[,genes])
    annot <- wgcna$datTraits
    gx.heatmap(X, nmax=nmax,
      ##col.annot = annot,
      key=FALSE, keysize=0.5, mar=heatmap.mar)

  }

  if(type == "correlation") {
    R <- cor(wgcna$datExpr[,genes])
    R <- sign(R) * abs(R)**rgamma
    ii <- hclust(as.dist(1-R), method="average")$order
    R <- R[ii,ii]
    R[abs(R) < min.rho] <- NA
    if(!diag) diag(R) <- NA
    image(R)
    mtext( rownames(R), side=4, at=seq(0,1,1/(nrow(R)-1)),
      las=2, adj = 0, cex=cex, line=0.5 )
  }
  
  if(is.null(main) && !is.null(module)) main <- module
  if(is.null(main)) main <- "Module Heatmap"
  if(!is.null(main) && main!="") {
    title( main, line=2, cex.main=1.3)  
  }
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
  iqr <- (S[3,] - S[1,])
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
                                    cex=1, maxpower = 20, nmax = 2000,
                                    plots=c("sft.modelfit", "mean.k",
                                      "dendro.IQR"), main=NULL, 
                                    RsquaredCut = 0.85, setPar=TRUE) {

  RsquaredCut <- RsquaredCut[1]
  
  ## Choose a set of soft-thresholding powers
  powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
  if(maxpower>20) {
    powers <- c(powers, seq(from = 20, to = maxpower, by = 5))
  }

  ## subsample for speed
  if(ncol(datExpr) > nmax && nmax > 0) {
    ii <- sample(1:ncol(datExpr),nmax)
    datExpr <- datExpr[,ii]
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
    ##if(legend) legend("bottomright", legend=paste("opt. power =",optPower))
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
    ht <- wgcna.checkDendroHeights(datExpr, n=200, powers=powers)
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

#'
#'
#' @export
wgcna.plotPowerAnalysis_multi <- function(exprList,
                                          cex=1, maxpower = 20,
                                          nmax = 2000,
                                          networktype = "signed",
                                          plots=c("sft.modelfit", "mean.k",
                                            "dendro.IQR"),
                                          main=NULL,
                                          cex.legend = 1,
                                          RsquaredCut = 0.85,
                                          setPar=TRUE) {
  if(0) {
    networktype = "signed"
    cex=1; maxpower = 20; nmax = 2000;
    plots=c("sft.modelfit", "mean.k","dendro.IQR");
    main=NULL; 
    RsquaredCut = 0.85
  }
  
  RsquaredCut <- RsquaredCut[1]
  
  ## Choose a set of soft-thresholding powers
  powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
  if(maxpower>20) {
    powers <- c(powers, seq(from = 20, to = maxpower, by = 5))
  }

  ## process each list
  sft <- list()
  for(i in 1:length(exprList)) {
    datExpr <- Matrix::t(exprList[[i]])

    ## subsample for speed
    if( nmax > 0 && nrow(exprList[[i]]) > nmax ) {
      ii <- sample(1:ncol(datExpr),nmax)
      datExpr <- datExpr[,ii]
    }
  
    ## Call the network topology analysis function
    k <- names(exprList)[i]
    sft[[k]] <- WGCNA::pickSoftThreshold(
      datExpr,
      powerVector = powers,
      RsquaredCut = RsquaredCut, 
      networkType = networktype,
      verbose = 0
    )
  }
  
  if(setPar) {
    np <- length(plots)
    par(mfrow = c(1, np), mar = c(3.8, 4.5, 3, 1), mgp = c(2.6, 0.95, 0))
  }
  
  ## Plot the results:
  if("sft.modelfit" %in% plots) {
    ## Scale-free topology fit index as a function of the soft-thresholding power
    Y <- c()
    for(i in 1:length(sft)) {
      y1 <- -sign(sft[[i]]$fitIndices[,"slope"]) * sft[[i]]$fitIndices[,"SFT.R.sq"]
      Y <- cbind(Y, y1)
    }
    colnames(Y) <- names(sft)
    x <- sft[[1]]$fitIndices[, "Power"]
    Y <- pmax(Y,0)
    matplot(
      x = x,
      y = Y,
      ylim = c(0,1),      
      type = "l",
      col = 2:99,
      lty = 1,
      lwd = 0.6,
      xlab = "Soft threshold (power)",
      ylab = "SFT model fit (signed R^2)",
      main = main
    )
    #abline(h = 0, col = "black", lty=1)
    for(i in 1:ncol(Y)) {      
      text(powers, Y[,i], labels = "█", cex=cex, col="white")      
      text(powers, Y[,i], labels = powers, cex=cex, col = 1+i)
    }
    ## this line corresponds to using an R^2 cut-off of h
    abline(h = RsquaredCut, col = "grey10", lty=2)
    legend("bottomright", legend=colnames(Y), fill=2:10,
      cex=cex.legend, y.intersp=0.9)
    title("SFT model fit")
  }
  
  ## Mean connectivity as a function of the soft-thresholding power
  if("mean.k" %in% plots) {
    Y <- sapply(sft, function(s) s$fitIndices[, "mean.k."])
    matplot(
      powers,
      Y,
      type = "l",
      col = 2:99,
      lty = 1,
      lwd = 0.6,
      xlab = "Soft threshold (power)",
      ylab = "Mean connectivity",
      main = main
    )
    for(i in 1:ncol(Y)) {
      text(powers, Y[,i], labels = "█", cex=cex, col="white")            
      text(powers, Y[,i], labels = powers, cex=cex, col = 1+i)
    }
    legend("topright", legend=colnames(Y), fill=2:10,
      cex=cex.legend, y.intersp=0.9)
    title("Mean connectivity")    
  }

  ht <- NULL
  if("dendro.IQR" %in% plots) {
    ht <- list()
    for(i in 1:length(exprList)) {
      ht[[i]] <- wgcna.checkDendroHeights(
        Matrix::t(exprList[[i]]), n=200, powers=powers)
    }
    Y <- sapply(ht, function(h) h$IQR)
    matplot(
      powers,
      Y,
      type = "l",
      col = 2:99,
      lty = 1,
      lwd = 0.6,
      xlab = "Soft threshold (power)",
      ylab = "Dendrogram height IQR",
      main = main
    )
    for(i in 1:ncol(Y)) {
      text(powers, Y[,i], labels = "█", cex=cex, col="white")            
      text(powers, Y[,i], labels = powers, cex = cex, col = 1+i)
    }
    legend("bottomright", legend=names(exprList), fill=2:10,
      cex=cex.legend, y.intersp=0.9)
    title("Dendrogram IQR")    
  }  

}



#' Better (?) method to pick soft threshold (aka power). 
#'
#' 
wgcna.pickSoftThreshold <- function(datExpr, sft=NULL, rcut=0.85,
                                    method = c("sft","iqr")[1],
                                    nmax = -1,  powers = NULL,
                                    verbose = 1) {

  if(is.null(powers)) {
    powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
    # powers <- c(powers, seq(from = 25, to = 50, by = 5))
  }

  ## subsample for speed
  if(ncol(datExpr) > nmax && nmax > 0) {
    ii <- sample(1:ncol(datExpr),nmax)
    datExpr <- datExpr[,ii]
  }
  
  if(is.null(sft)) {
    sft <- WGCNA::pickSoftThreshold(
      datExpr,
      powerVector = powers,
      networkType = "signed",
      verbose = verbose
    )
  }

  optPower <- NULL
  if(method == "sft") {
    ## Pick power according to scale-free (SFT) parameter
    sqr <- -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]
    if(max(sqr, na.rm=TRUE) >= rcut) {
      optPower <- min(powers[which(sqr >= rcut)])
    } else {
      ## remove initial value that are possible negative
      if(sqr[1] < 0.05) {
        for(i in 1:length(sqr)) sqr[i] <- ifelse(sqr[i] < 0.05, NA, sqr[i])
      }
      ds <- 0.5*median(abs(diff(sqr)),na.rm=TRUE) ##small step
      if(any(diff(sqr) < -ds, na.rm=TRUE)) {
        i <- min(which(diff(sqr) < -ds))+1
        sqr[i:length(sqr)] <- NA
      }
      optPower <- powers[which.max(sqr)]
    }
  } else if(method == "iqr") {
    ## Pick power according to IQR
    ht <- wgcna.checkDendroHeights(datExpr, n=200, powers = powers)
    ##base::plot( powers, ht$IQR)
    optPower <- powers[ which.max(ht$IQR) ]
  } else {
    stop("[wgcna.pickSoftThreshold] invalid method = ", method)
  }
  
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

