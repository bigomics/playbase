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
    power = 6,
    cutheight = 0.15,
    deepsplit = 4,
    networktype = "signed",
    tomtype = "signed",                          
    ngenes = 1000) {
  ## if we dont source WGCNA, blockwiseModules does not work..
  require(WGCNA)

  # minmodsize=30;power=6;cutheight=0.25;deepsplit=2;ngenes=1000
  if(is.null(networktype)) networktype <- "signed"
  if(is.null(tomtype)) tomtype <- "signed"  
  
  res <- wgcna.compute(
    X = pgx$X,
    samples = pgx$samples,
    minmodsize = minmodsize,   # default: min(20,...)
    power = power,             # default: 6
    cutheight = cutheight,     # default: 0.15
    deepsplit = deepsplit,     # default: 2
    networktype = networktype, # default: unsigned 
    tomtype = tomtype,         # default: signed                  
    ngenes = ngenes)

  ##---------------------------------------------------
  ## compute clustering based on TOM matrix
  ##---------------------------------------------------
  dissTOM <- 1 - res$TOM
  rownames(dissTOM) <- colnames(dissTOM) <- colnames(res$datExpr)
  clust <- pgx.clusterBigMatrix(dissTOM, methods = c("umap", "tsne", "pca"), dims = c(2))
  if ("cluster.genes" %in% names(pgx)) {
    posx <- pgx$cluster.genes$pos[["umap2d"]]
    clust[["umap2d"]] <- posx[colnames(res$datExpr), ]
  }

  ##----------------------------------------------------
  ## Do quick geneset analysis 
  ##----------------------------------------------------
  gmt <- getGSETS_playbase(pattern = "HALLMARK|GOBP|^C[1-9]|GO_BP")
  gse <- NULL
  bg <- rownames(pgx$X)
  bg <- probe2symbol(bg, pgx$genes, query = "human_ortholog")
  bg <- toupper(bg[!is.na(bg)])
  ## bg <- toupper(rownames(pgx$X))

  gmt1 <- lapply(gmt, function(m) intersect(m, bg))
  gmt1.size <- sapply(gmt1, length)
  gmt1 <- gmt1[gmt1.size >= 10]

  ## Perform fisher-test (fastest method)
  i <- 1
  me.genes <- res$me.genes
  for (i in 1:length(me.genes)) {
    ## gg <- toupper(me.genes[[i]])
    gg <- probe2symbol(me.genes[[i]], pgx$genes, query = "human_ortholog")
    gg <- toupper(gg)
    rr <- try(gset.fisher(gg, gmt1, background = bg, fdr = 1, min.genes = 10))
    if (!"try-error" %in% class(rr)) {
      rr <- cbind(module = names(me.genes)[i], geneset = rownames(rr), rr)
      rr <- rr[order(rr$p.value), , drop = FALSE]
      if (is.null(gse)) {
        gse <- rr
      } else {
        gse <- rbind(gse, rr)
      }
    }
  }
  rownames(gse) <- NULL

  ## construct results object
  return(
    list(
      datExpr = res$datExpr,
      datTraits = res$datTraits,
      net = res$net,
      gse = gse,
      clust = clust,
      power = power,
      networktype = networktype,
      tomtype = tomtype,
      me.genes = res$me.genes,
      me.colors = res$me.colors
    )
  )
}


#' @export
wgcna.compute <- function(X, samples,
                          minmodsize = 20,
                          power = 6,
                          cutheight = 0.15,
                          deepsplit = 4,
                          networktype = "signed",
                          tomtype = "signed",                          
                          ngenes = 1000) {

  ##minmodsize=20;power=6;cutheight=0.15;deepsplit=2;ngenes=1000;networktype="signed";tomtype = "signed"

  nmissing <- sum(is.na(X))
  if (nmissing > 0) {
    message("Found ", nmissing, " missing values in X. Imputing prior to WGCNA.")
    ##X <- X[complete.cases(X), , drop = FALSE]
    X <- svdImpute2(X)
  }
  X <- X[!duplicated(rownames(X)), ]

  if( nrow(X) > ngenes ) {
    sdx <- matrixStats::rowSds(X, na.rm = TRUE)
    X <- X[sdx > 0.1 * mean(sdx, na.rm = TRUE), ] ## filter low SD??
    X <- X[order(-matrixStats::rowSds(X, na.rm = TRUE)), ]
    X <- utils::head(X, ngenes)
  }

  datExpr <- t(X)

  ## adapt for small datasets
  minmodsize = min(minmodsize, ncol(datExpr)/2 )
  minmodsize
  message("[wgcna.compute] minmodsize = ", minmodsize)
  
  WGCNA::enableWGCNAThreads()
  cor <- WGCNA::cor
  ## we iterate over minmodsize, until we get less than 10 clusters
  net <- WGCNA::blockwiseModules(
    datExpr,
    power = power,
    networkType = networktype,
    TOMType = tomtype,
    minModuleSize = minmodsize,
    # reassignThreshold = 0,
    mergeCutHeight = cutheight,
    numericLabels = TRUE,
    deepSplit = deepsplit,
    verbose = 0
  )
  table(net$colors)
  cor <- stats::cor
  
  ## clean up traits matrix
  datTraits <- samples
  ## no dates please...
  isdate <- apply(datTraits, 2, is.Date)
  datTraits <- datTraits[, !isdate, drop = FALSE]

  ## Expand multi-class discrete phenotypes into binary vectors
  tr.class <- sapply(utils::type.convert(datTraits, as.is = TRUE), class)
  sel1 <- which(tr.class %in% c("factor", "character"))
  sel2 <- which(tr.class %in% c("integer", "numeric"))

  tr1 <- datTraits[, 0]
  if (length(sel1)) {
    tr1 <- expandPhenoMatrix(datTraits[, sel1, drop = FALSE], drop.ref = FALSE)
    if (is.null(tr1)) tr1 <- datTraits[, 0]
  }
  ## keeping numeric phenotypes
  tr2 <- datTraits[, sel2, drop = FALSE]
  datTraits <- cbind(tr1, tr2)

  ## get colors of eigengene modules
  me.genes <- tapply(names(net$colors), net$colors, list)
  names(me.genes) <- paste0("ME", names(me.genes))
  color1 <- labels2rainbow(net)
  me.colors <- color1[!duplicated(color1)]
  names(me.colors) <- paste0("ME", names(me.colors))
  me.colors <- me.colors[names(me.genes)]

  ## compute clustering based on TOM matrix
  TOM <- WGCNA::TOMsimilarityFromExpr(
    datExpr,
    power = power,
    TOMType = tomtype,
    networkType = networktype    
  )
  rownames(TOM) <- colnames(TOM) <- colnames(datExpr)

  ## compute module eigenvectors (loading matrix)
  MVs <- list()
  for(clr in unique(net$colors)) {
    ii <- which(net$colors == clr)
    mX <- datExpr[,ii]
    mX <- scale(mX)  ## NOTE: seems WGCNA is using full scaling
    ##mX <- scale(mX, scale=FALSE)  ## better????  
    sv1 <- irlba::irlba(mX, nv=1, nu=1)
    ##plot(sv1$u, net$MEs[[paste0("ME",clr)]])
    ##plot(rowMeans(mX), net$MEs[[paste0("ME",clr)]])
    sv <- rep(0, ncol(datExpr))
    sv[ii] <- sv1$v[,1] * sv1$d[1]
    MVs[[paste0("ME",clr)]] <- sv
  }
  MVs <- as.matrix(do.call(cbind, MVs[names(net$MEs)]))
  rownames(MVs) <- colnames(datExpr)

  ## create loading matrix W
#  W <- sapply( me.genes, function(gg) 1*(colnames(datExpr) %in% gg))
#  rownames(W) <- colnames(datExpr)
#  W <- W * colMeans(datExpr)
  
  ## construct results object
  results <- list(
      datExpr = datExpr,
      datTraits = datTraits,
      TOM = TOM,
      net = net,
      power = power,
      me.genes = me.genes,
      me.colors = me.colors,
      W = MVs
  )

  return(results)
}


#'
#'
#' @export
wgcna.plotTOM <- function(results, power=NULL, networktype="signed",
                          tomtype="signed", nSelect=1000) {

  datExpr <- results$datExpr
  MEs <- results$net$MEs
  moduleColors <- labels2rainbow(results$net)
  if(is.null(power) && !is.null(results$power)) {
    power <- results$power
  }
  if(is.null(power)) power <- 6
  
  ## Calculate topological overlap anew: this could be done
  ## more efficiently by saving the TOM calculated during
  ## module detection, but let us do it again here.
  TOM <- TOMsimilarityFromExpr(
    datExpr,
    power = power,
    networkType = networktype, 
    TOMType = tomtype,
    verbose = 0
  )
  dissTOM <- 1 - TOM
  rownames(dissTOM) <- colnames(dissTOM) <- colnames(datExpr)
  
  #nSelect <- 999999
  #nSelect <- 800
  ## For reproducibility, we set the random seed
  set.seed(10)
  select <- head(1:ncol(dissTOM), nSelect)
  selectTOM <- dissTOM[select, select]
  ## There’s no simple way of restricting a clustering tree
  ## to a subset of genes, so we must re-cluster.
  #
  selectTree <- hclust(as.dist(selectTOM), method = "average")
  selectColors <- moduleColors[select]

  ## Taking the dissimilarity to a power, say 10, makes the plot
  ## more informative by effectively changing the color palette;
  ## setting the diagonal to NA also improves the clarity of the
  ## plot
  plotDiss <- selectTOM^7
  diag(plotDiss) <- NA
  myheatcol <- gplots::colorpanel(250, "red", "orange", "lemonchiffon")
  myheatcol <- gplots::colorpanel(250, "lemonchiffon", "orange", "red")
  
  par(oma = c(2, 0, 0, 0))
  plotly::layout(
    matrix(c(
      0, 0, 5, 0,
      0, 0, 2, 0,
      4, 1, 3, 6
    ), nr = 3, byrow = TRUE),
    widths = c(2.3, 0.5, 10, 1.8),
    heights = c(2.3, 0.5, 10)
  )
  
  WGCNA::TOMplot(
    plotDiss,
    selectTree,
    selectColors,
    col = myheatcol,
    setLayout = FALSE,
    main = NULL
  )
  
  ## add color legend
  frame()
  me.names <- colnames(MEs)
  me.nr <- as.integer(sub("ME", "", me.names))
  ii <- order(me.nr)
  label.colors <- labels2rainbow(results$net)
  me.colors <- label.colors[!duplicated(names(label.colors))]
  me.colors <- me.colors[as.character(me.nr)]
  
  legend(
    -0.1, 1,
    legend = me.names[ii], fill = me.colors[ii],
    cex = 1.2, bty = "n", x.intersp = 0.5
  )

}


#'
#'
#' @export
wgcna.plotDendroAndColors <- function(results, power=NULL, networktype="signed",
                                      tomtype="signed", nSelect=1000) {

  datExpr <- results$datExpr
  MEs <- results$net$MEs
  moduleColors <- labels2rainbow(results$net)
  if(is.null(power) && !is.null(results$power)) {
    power <- results$power
  }
  if(is.null(power)) power <- 6
  
  ## Calculate topological overlap anew: this could be done
  ## more efficiently by saving the TOM calculated during
  ## module detection, but let us do it again here.
  TOM <- TOMsimilarityFromExpr(
    datExpr,
    power = as.numeric(power),
    networkType = networktype, 
    TOMType = tomtype,
    verbose = 0
  )
  dissTOM <- 1 - TOM
  rownames(dissTOM) <- colnames(dissTOM) <- colnames(datExpr)
  
  ##nSelect <- 999999
  ##nSelect <- 800
  ## For reproducibility, we set the random seed
  set.seed(10)
  select <- head(1:ncol(dissTOM), nSelect)
  selectTOM <- dissTOM[select, select]
  ## There’s no simple way of restricting a clustering tree
  ## to a subset of genes, so we must re-cluster.
  #
  selectTree <- hclust(as.dist(selectTOM), method = "average")
  selectColors <- moduleColors[select]

  ## Convert labels to colors for plotting
  ## Plot the dendrogram and the module colors underneath
  plotDendroAndColors(
    dendro = selectTree,
    colors = selectColors,
    dendroLabels = FALSE, hang = 0.03,
    addGuide = FALSE, guideHang = 0.05,
    marAll = c(0.2, 5, 0.4, 0.2),
    main = NULL
  )
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
#' @export
labels2rainbow <- function(net) {
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


#' @export
wgcna.plotModuleTraitHeatmap <- function(results) {
  
  ## Define numbers of genes and samples
  nGenes = ncol(results$datExpr);
  nSamples = nrow(results$datExpr);
  ## Recalculate MEs with color labels
  moduleColors <- results$net$colors
  MEs0 = moduleEigengenes(results$datExpr, moduleColors)$eigengenes
  MEs = orderMEs(MEs0)
  moduleTraitCor = cor(MEs, results$datTraits, use = "pairwise.complete");
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
  
  sizeGrWindow(10,6)
  ## Will display correlations and their p-values
  textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                      signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
  par(mar = c(6, 8.5, 3, 3));
  ## Display the correlation values within a heatmap plot
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(results$datTraits),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))

}
