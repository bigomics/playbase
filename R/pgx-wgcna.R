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
    deepsplit = 2,
    networktype = "signed",
    tomtype = "signed",
    numericlabels = FALSE,
    ngenes = 1000) {

  ##minmodsize=10;power=NULL;cutheight=0.25;deepsplit=2;ngenes=1000;networktype="signed";tomtype="signed";numericlabels=FALSE

  samples <- pgx$samples
  samples <- samples[, grep("^[.]",colnames(samples)),drop=FALSE]  ## no dot pheno
  
  res <- wgcna.compute(
    X = pgx$X,
    samples = samples,
    minmodsize = minmodsize,   # default: min(20,...)
    power = power,             # default: 6
    cutheight = cutheight,     # default: 0.15
    deepsplit = deepsplit,     # default: 2
    networktype = networktype, # default: unsigned (but signed is better...)
    tomtype = tomtype,         # default: signed                  
    numericlabels = numericlabels,
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
  gse <- wgcna.compute_enrichment(res, pgx, method="fisher") 
  
  ## add to results object
  res$gse <- gse
  res$clust <- clust
  res$networktype <- networktype
  res$tomtype <- tomtype
  
  return(res)
}


#' @export
wgcna.compute <- function(X, samples,
                          minmodsize = 20,
                          power = 6,
                          cutheight = 0.15,
                          deepsplit = 2,
                          networktype = "signed",
                          tomtype = "signed",                          
                          ngenes = 1000,
                          numericlabels = TRUE,
                          prefix = "ME"
                          ) {

  #minmodsize=20;power=6;cutheight=0.15;deepsplit=2;ngenes=1000;networktype="signed";tomtype="signed";numericlabels=FALSE;prefix="ME"

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

  if(is.null(power)) {
    ## Estimate best power
    powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
    sft <- WGCNA::pickSoftThreshold(
      datExpr,
      powerVector = powers,
      networkType = networktype,
      verbose = 0
    )
    power <- sft$powerEstimate
    if(is.na(power)) power <- 20
    message("[wgcna.compute] estimated power = ", power)    
  } else {
    message("[wgcna.compute] power = ", power)    
  }
  
  ## adapt for small datasets (also done in WGCNA package)
  minmodsize = min( minmodsize, ncol(datExpr)/2 )
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
    numericLabels = numericlabels, ## numeric or 'color' labels
    deepSplit = deepsplit,
    verbose = 0
  )
  cor <- stats::cor

##  prefix="ME"
  table(net$colors)
  names(net$MEs) <- sub("^ME",prefix,names(net$MEs))
  net$labels <- paste0(prefix,net$colors)

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
  
  ## compute clustering based on TOM matrix
  TOM <- WGCNA::TOMsimilarityFromExpr(
    datExpr,
    power = power,
    TOMType = tomtype,
    networkType = networktype,
    verbose = 0
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
    sv <- rep(0, ncol(datExpr))
    sv[ii] <- sv1$v[,1] * sv1$d[1]
    MVs[[paste0(prefix,clr)]] <- sv
  }
  MVs <- as.matrix(do.call(cbind, MVs[names(net$MEs)]))
  rownames(MVs) <- colnames(datExpr)

  ## compute gene statistics
  stats <- wgcna.compute_geneStats(net, datExpr, datTraits, TOM=TOM) 
  
  ## construct results object
  results <- list(
      datExpr = datExpr,
      datTraits = datTraits,
      TOM = TOM,
      net = net,
      power = power,
      me.genes = me.genes,
      me.colors = me.colors,
      W = MVs,
      stats = stats
  )

  return(results)
}

wgcna.compute_geneStats <- function(net, datExpr, datTraits, TOM) {
  
  ## Define numbers of genes and samples
  nGenes = ncol(datExpr);
  nSamples = nrow(datExpr);  

  ## Recalculate MEs with color labels
  moduleTraitCor = cor(net$MEs, datTraits, use = "pairwise.complete");
  moduleTraitPvalue = WGCNA::corPvalueStudent(moduleTraitCor, nSamples);

  ## Module membership correlation (with p-values)
  moduleMembership = cor(datExpr, net$MEs, use = "p");
  MMPvalue = WGCNA::corPvalueStudent(as.matrix(moduleMembership), nSamples);
  
  ## Gene-trait significance (trait correlation) (with p-values)
  traitSignificance = cor(datExpr, datTraits, use = "p");
  GSPvalue = WGCNA::corPvalueStudent(as.matrix(traitSignificance), nSamples);

  ## Fold-change
  is.binary <- apply(datTraits,2,function(a) length(unique(a[!is.na(a)]))==2)
  Y <- datTraits[,which(is.binary)]
  lm <- lapply(Y, function(y) gx.limma( t(datExpr), y, lfc=0, fdr=1,
                                               sort.by='none', verbose=0))
  foldChange <- sapply(lm, function(m) m$logFC)
  foldChangePvalue <- sapply(lm, function(m) m$P.Value)  
  rownames(foldChange) <- rownames(lm[[1]])
  rownames(foldChangePvalue) <- rownames(lm[[1]])  

  ## Gene Centrality. Compute centrality of gene in Module subgraph
  ## using TOM matrix.
  adj <- TOM
  diag(adj) <- NA
  adj[which(abs(adj) < 0.01)] <- 0
  gr <- igraph::graph_from_adjacency_matrix(
    adj, mode="undirected", weighted=TRUE, diag=FALSE )
  geneCentrality <- rep(NA, nrow(adj))
  names(geneCentrality) <- rownames(adj)
  me.genes <- tapply( names(net$colors), net$colors, list)
  gg <- me.genes[[1]]
  for(gg in me.genes) {
    gr1 <- igraph::subgraph(gr, gg)
    ct <- igraph::page_rank(gr1, weights=NULL)$vector
    ct <- ct / mean(ct,na.rm=TRUE)
    geneCentrality[gg] <- ct
  }

  ## propotion variance explained (PVE) per module
  propVarExplained <- WGCNA::propVarExplained(
    datExpr, colors=net$colors, MEs=net$MEs)
  
  ##
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

wgcna.compute_enrichment <- function(res, pgx, method="fisher", ntop=1000) {

  ##----------------------------------------------------
  ## Do quick geneset analysis 
  ##----------------------------------------------------
  symbol.col <- intersect(c("symbol","gene_name"),colnames(pgx$genes))[1]
  bg <- probe2symbol(rownames(pgx$X), pgx$genes, query = symbol.col)  
  bg <- intersect(bg, rownames(pgx$GMT))
  G1 <- pgx$GMT[bg,]
  
  sel <- grep("PATHWAY|HALLMARK|^GO|^C[1-9]",colnames(pgx$GMT))
  G1 <- pgx$GMT[,sel]
  G1 <- G1[, which(Matrix::colSums(G1!=0) >= 4)]

  W <- as.matrix(res$W)
  W <- rename_by(W, pgx$genes, symbol.col)  

  X <- t(res$datExpr)
  X <- rename_by(X, pgx$genes, symbol.col)  
  
  if(method == "gsetcor") {
    gsetX <- pgx$gsetX[colnames(G1),]
    dim(gsetX)
    ME <- as.matrix(res$net$MEs)
    rc <- gset.rankcor(ME, gsetX, compute.p=TRUE)  ## NEEDS CHECK!!!
    head(rc$rho)
    gse <- reshape2::melt(rc$rho)
    gs.p <- reshape2::melt(rc$p.value)
    gs.q <- reshape2::melt(rc$q.value)
    colnames(gse) <- c("geneset","module","rho")
    gse$p.value <- gs.p$value
    gse$q.value <- gs.q$value
    gse$score <- gse$rho * -log10(gse$p.value)
    head(gse)
  }
  
  if(method == "rankcor") {
    
    mm <- cor( t(X), res$net$MEs)
    rc <- gset.rankcor(mm, G1, compute.p=TRUE)  ## NEEDS CHECK!!!
    head(rc$rho)
    gse  <- reshape2::melt(rc$rho)
    gs.p <- reshape2::melt(rc$p.value)
    gs.q <- reshape2::melt(rc$q.value)
    colnames(gse) <- c("geneset","module","rho")
    gse$p.value <- gs.p$value
    gse$q.value <- gs.q$value
    gse$score <- gse$rho * -log10(gse$p.value)
    head(gse)
  }

  if(method == "fisher") {
    ## Perform fisher-test (fastest method)
    gmt1 <- mat2gmt(G1)
    i=1
    gse <- NULL
    me.genes <- res$me.genes
    for (i in 1:length(me.genes)) {
      ## gg <- toupper(me.genes[[i]])
      gg <- probe2symbol(me.genes[[i]], pgx$genes, query = symbol.col)
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
      rownames(gse) <- NULL
    }
    ## handle infinite
    gse$odd.ratio[is.infinite(gse$odd.ratio)] <- 99
    gse$score <- gse$odd.ratio * -log10(gse$p.value)
  }
  
  dbg("[pgx.wgcna] dim.gse = ", dim(gse))
  gse <- gse[order(-gse$score),]
  gse.list <- tapply( 1:nrow(gse), gse$module, function(ii) gse[ii,] )  
  if(!is.null(ntop) && ntop>0) gse.list <- lapply(gse.list, head, n=ntop)
  return(gse.list)
}


#'
#'
#' @export
wgcna.getGeneStats <- function(res, module, trait, plot=TRUE,
                               showallmodules=TRUE, col=NULL,
                               main=NULL) {
  # module="MEblue";trait="activated=act"
  p1 <- c("moduleMembership","MMPvalue")
  p2 <- c("traitSignificance","GSPvalue","foldChange","foldChangePvalue")
  p3 <- c("geneCentrality")
  nrow <- ncol(res$datExpr)
  df <- res$stats[[p1[1]]][,0]  ## empty data.frame
  rownames(df) <- colnames(res$datExpr)
  if(!is.null(module)) {
    A1 <- sapply( res$stats[p1], function(x) x[,module])
    df <- cbind(df, A1)
  }
  if(!is.null(trait)) {
    A2 <- sapply( res$stats[p2], function(x) x[,trait])
    df <- cbind(df, A2)
  }
  A3 <- res$stats[[p3]]
  df <- cbind(df, centrality=A3)

  sel <- c("moduleMembership","traitSignificance","foldChange","centrality")
  sel <- intersect(sel, colnames(df))
  df1 <- pmax( as.matrix(df[,sel]), 1e-8)
  score <- exp(rowMeans(log(df1)))

  labels <- res$net$labels
  df <- data.frame( module = labels, score=score, df )
  df <- df[order(-df$score),]
  if(!is.null(module) && !showallmodules) {
    sel <- which( df$module == module )
    df <- df[sel,,drop=FALSE]
  }
  
  if(plot) {
    cols <- c("moduleMembership","traitSignificance","foldChange","centrality")
    cols <- intersect(cols, colnames(df))
    df1 <- df[,cols]
    col1 <- wgcna.labels2colors(res$net$colors[rownames(df1)])
    if(!is.null(col)) col1 <- col
    pairs(df1, col=col1)
    if(is.null(main)) {
      main <- paste("Gene significance for module",module,"and trait",trait)
    }
    title(main, line=3, cex.main=1.15)
  }

  df
}


#'
#'
#' @export
wgcna.plotTOM <- function(res, justdata=FALSE) {

  datExpr <- res$datExpr
  MEs <- res$net$MEs
  moduleColors <- wgcna.labels2colors(res$net$colors)  

  ## Topological overlap dissimilarity matrix
  dissTOM <- 1 - res$TOM
  rownames(dissTOM) <- colnames(dissTOM) <- colnames(datExpr)

  ## clustering results
  geneTree <- res$net$dendrograms[[1]]
  
  ## Taking the dissimilarity to a power, say 10, makes the plot
  ## more informative by effectively changing the color palette;
  ## setting the diagonal to NA also improves the clarity of the
  ## plot
  plotDiss <- dissTOM^7
  diag(plotDiss) <- NA
  myheatcol <- gplots::colorpanel(250, "red", "orange", "lemonchiffon")
  myheatcol <- gplots::colorpanel(250, "lemonchiffon", "orange", "red")

  if(justdata) return(plotDiss)
  
  par(oma = c(2, 0, 0, 0))
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
    plotDiss,
    geneTree,
    moduleColors,
    #col = myheatcol,
    setLayout = FALSE,
    main = NULL
  )
  
  ## add color legend
  frame()
  legend(
    -0.1, 1,
    fill = res$me.colors,
    legend = names(res$me.colors),
    cex = 1.2, bty = "n", x.intersp = 0.5
  )

}


#'
#'
#' @export
wgcna.plotDendroAndColors <- function(res, main=NULL, unmerged=FALSE) {

  colors <- wgcna.labels2colors(res$net$colors)
  groupLabels <- "Module colors"
  geneTree = res$net$dendrograms[[1]]
  if(unmerged) {
    colors <- cbind(colors,
      wgcna.labels2colors(res$net$unmergedColors))
      groupLabels <- c( "Merged colors", "Unmerged colors")
  }
  
  if(is.null(main)) main <- "Gene dendrogram and module colors"
  ## Plot the dendrogram and the module colors underneath
  plotDendroAndColors(
    dendro = geneTree,
    colors = colors,
    groupLabels = groupLabels,
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = FALSE, guideHang = 0.05,
    marAll = c(0.2, 5, 1, 0.2),
    main = main
  )

}


#' @export
wgcna.labels2colors <- function(colors, ...) {
  if(all(is.numeric(colors))) {
    colors <- WGCNA::labels2colors(colors, ...)
    return(colors)
  }
  stdColors <- c("grey",WGCNA::standardColors())
  if(all(colors %in% stdColors)) {
    return(colors)
  }
  icolors <- as.integer(factor(as.character(colors)))
  colors <- WGCNA::standardColors()[icolors]
  return(colors)
}


#' @export
wgcna.plotModuleTraitHeatmap <- function(res, setpar=TRUE, cluster=FALSE,
                                         main = NULL, justdata=FALSE,
                                         pstar = TRUE) {
  
  ## Define numbers of genes and samples
  nGenes = ncol(res$datExpr);
  nSamples = nrow(res$datExpr);
  if("stats" %in% names(res)) {
    moduleTraitCor = res$stats$moduleTraitCor
    moduleTraitPvalue <- res$stats$moduleTraitPvalue
  } else {
    moduleTraitCor = cor(res$net$MEs, res$datTraits, use = "pairwise.complete");
    moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
  }
  
  if(cluster) {
    ii <- hclust(dist(moduleTraitCor))$order
    jj <- hclust(dist(t(moduleTraitCor)))$order
  } else {
    ii <- rev(order(rownames(moduleTraitCor)))
    jj <- order(colnames(moduleTraitCor))  
  }
  moduleTraitCor <- moduleTraitCor[ii,jj]
  moduleTraitPvalue <- moduleTraitPvalue[ii,jj]    

  if(justdata) return(moduleTraitCor)
  
  ## Will display correlations and their p-values
  if(pstar) {
    textPv = cut(moduleTraitPvalue, breaks=c(-1,0.001,0.01,0.05,99),
                 labels = c("★★★","★★","★",""))
  } else {
    textPv = paste0("(",signif(moduleTraitPvalue, 1), ")")
  }
  
  textMatrix = paste0(signif(moduleTraitCor, 2), "\n", textPv)
  dim(textMatrix) = dim(moduleTraitCor)
  if(setpar) par(mar = c(8, 8, 3, 3));
  if(is.null(main)) main <- "Module-trait relationships"
  ## Display the correlation values within a heatmap plot
  WGCNA::labeledHeatmap(Matrix = moduleTraitCor,
                        xLabels = colnames(moduleTraitCor),
                        yLabels = rownames(moduleTraitCor),
                        ySymbols = rownames(moduleTraitCor),
                        colorLabels = FALSE,
                        colors = blueWhiteRed(50),
                        textMatrix = textMatrix,
                        setStdMargins = FALSE,
                        cex.text = 0.7,
                        zlim = c(-1,1),
                        main = main)
  
}

#' Plot membership correlation vs gene signficance (correlation with
#' trait) to discover biomarkers/driver genes.
#' 
#' @export
wgcna.plotMMvsGS <- function(res, module, trait, abs=TRUE, par=TRUE,
                             plotlib = "base") {
  ##module="ME3";trait="activated=act"
  moduleGenes = res$me.genes[[module]]
  nSamples = nrow(res$datExpr)
  
  ## Module membership correlation (with p-values)
  if("stats" %in% names(res)) {
    moduleMembership <- res$stats$moduleMembership
    MMPvalue <- res$stats$MMPvalue
  } else {
    moduleMembership = as.data.frame(cor(res$datExpr, res$net$MEs, use = "p"));
    MMPvalue = as.data.frame(corPvalueStudent(as.matrix(moduleMembership), nSamples));
  }

  ## Gene-trait significance (trait correlation) (with p-values)
  if("stats" %in% names(res)) {
    traitSignificance <- res$stats$traitSignificance
    GSPvalue <- res$stats$GSPvalue
  } else {
    traitSignificance = as.data.frame(cor(res$datExpr, res$datTraits, use = "p"));
    GSPvalue = as.data.frame(corPvalueStudent(as.matrix(traitSignificance), nSamples));
  }
    
  x <- (moduleMembership[moduleGenes, module])
  y <- (traitSignificance[moduleGenes, trait])
  if(abs==TRUE) {
    x <- abs(x)
    y <- abs(y)
  }
  ## 
  px <- MMPvalue[moduleGenes, module]
  py <- GSPvalue[moduleGenes, trait]
  qx <- p.adjust(px, method="fdr")
  qy <- p.adjust(py, method="fdr")
  is.sig <- ( qx < 0.05 & qy < 0.05 )
  sigx <- (qx < 0.05)
  sigy <- (qy < 0.05)  
  ii <- which(is.sig)
  qv <- quantile(x[ii],prob=0.1)[1]
  qh <- quantile(y[ii],prob=0.1)[1]

  pos <- cbind(x,y)
  rownames(pos) <- moduleGenes
  is.sig1 <- c("notsig","onesig","sig")[1 + 1*sigx + 1*sigy]
  hi1 <- NULL
  ##hi1 <- head(rownames(pos),10)
  col1 <- c("grey70","grey20","red2")

  if(par) par(mfrow = c(1,1), mar=c(5,5,3,2));
  if(plotlib == "ggplot") {
    pgx.scatterPlotXY.GGPLOT(
      pos, var=is.sig1, hilight=hi1, col=col1,
      xlab = paste("Module membership in", module, "module"),
      ylab = paste("Gene significance for trait",trait),
      title = paste("Module membership vs. gene significance\n"),
      cex.title = 0.9, 
      girafe=FALSE)
  } else if(plotlib == "girafe") {
    pgx.scatterPlotXY.GGPLOT(
      pos, var=is.sig1, hilight=hi1, col=col1,
      xlab = paste("Module membership in", module, "module"),
      ylab = paste("Gene significance for trait",trait),
      title = paste("Module membership vs. gene significance\n"),
      cex.title = 0.7, cex.axis = 0.7,
      girafe=TRUE)
  } else {
    ii <- which( is.sig1 == "notsig")
    verboseScatterplot(
      x[-ii], y[-ii],
      xlab = paste("Module membership in", module, "module"),
      ylab = paste("Gene significance for trait",trait),
      main = paste("Module membership vs. gene significance\n"),
      cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = col1[1])
    ii <- which( is.sig1 == "onesig")
    points( x[ii], y[ii], col=col1[2] )
    ii <- which( is.sig1 == "sig")
    points( x[ii], y[ii], col=col1[3] )    
    abline(v=qv, h=qh, col="darkred")  
  }  
}


#' Plot cluster dendrogram with eigengenes and traits.
#' 
#' @export
wgcna.plotEigenGeneClusterDendrogram <- function(res, add_traits=TRUE,
                                                 main = NULL) {
  # Matrix with eigengenes and traits
  MET <- res$net$MEs
  if(add_traits) {
    MET <- cbind(MET, res$datTraits)
  }
  MET = WGCNA::orderMEs(MET)
  if (NCOL(MET) <= 2) MET <- cbind(MET, MET) ## error if ncol(MET)<=2 !!!!  
  if(is.null(main)) main <- "Eigengene cluster dendrogram"
  WGCNA::plotEigengeneNetworks(
    MET, main,
    marDendro = c(0,4,2,0),
    plotHeatmaps = FALSE)

}

#' Plot the adjacency correlation heatmap matrix of eigengenes with or
#' without traits. This can show how traits cluster together with the
#' eigengenes.
#' 
#' @export
wgcna.plotEigenGeneAdjacencyHeatmap <- function(res, add_traits=TRUE,
                                                marx=1, main=NULL,
                                                justdata=FALSE) {
  # Matrix with eigengenes and traits
  MET <- res$net$MEs
  if(add_traits) {
    MET <- cbind(MET, res$datTraits)
  }
  colnames(MET) <- paste0(" ",colnames(MET))
  MET = WGCNA::orderMEs(MET)
  if (NCOL(MET) <= 2) MET <- cbind(MET, MET) ## error if ncol(MET)<=2 !!!!

  if(justdata) {
    R <- (1 + cor(MET))/2
    return(R)
  }
  
  # Plot the correlation heatmap matrix (note: this plot will overwrite
  # the dendrogram plot)
  if(is.null(main)) main <- "Eigengene adjacency heatmap"
  WGCNA::plotEigengeneNetworks(
    MET, main,
    marHeatmap = c(8*marx,10*marx,2,2),
    plotDendrograms = FALSE,
    colorLabels = TRUE,
    xLabelsAngle = 45)
  
}

#' @export
wgcna.plotEigenGeneGraph <- function(res, add_traits=TRUE, main=NULL) {
  ##require(igraph)
  net <- res$net
  MET <- net$MEs
  if(add_traits) {
    MET <- cbind(MET, res$datTraits)
  }
  if (NCOL(MET) <= 2) MET <- cbind(MET, MET) ## error if ncol(MET)<=2 !!!!
  
  ## Recalculate MEs with color as labels
  clust <- hclust(dist(t(scale(MET))))
  clust
  phylo <- ape::as.phylo(clust)
  gr <- igraph::as.igraph(phylo, directed = FALSE)
  
  is.node <- grepl("Node", igraph::V(gr)$name)
  module.name <- igraph::V(gr)$name
  module.size <- table(res$net$labels)
  module.size <- module.size / mean(module.size)
  module.size

  igraph::V(gr)$label <- igraph::V(gr)$name
  igraph::V(gr)$label[is.node] <- NA
  igraph::V(gr)$color <- res$me.colors[module.name]
  igraph::V(gr)$size <- 20 * (module.size[module.name])**0.4
  igraph::V(gr)$size[is.na(igraph::V(gr)$size)] <- 0
  
  ##par(mfrow = c(1, 1), mar = c(1, 1, 1, 1) * 0)
  igraph::plot.igraph(
    gr,
    layout = igraph::layout.kamada.kawai,
    vertex.label.cex = 0.8,
    edge.width = 3
  )
  if(!is.null(main)) title(main, line=-1.5)  
}


#' Plot Multi-dimensional scaling (MDS) of centered data matrix.
#'
#' @export
wgcna.plotMDS <- function(res, main=NULL, scale=FALSE) {
  cc <- wgcna.labels2colors(res$net$color)
  pc <- svd(t(scale(res$datExpr,scale=scale)),nv=2)$u[,1:2]
  #pc <- svd(t(scale(res$datExpr)),nv=1)$u[,1:2]
  colnames(pc) <- c("MDS-x","MDS-y")
  if(is.null(main)) main <- "MDS of features"
  plot(pc, col=cc, main=main)

}

#' Plot Multi-dimensional scaling (MDS) of centered data matrix.
#'
#' @export
wgcna.plotModuleSignificance <- function(res, trait, main=NULL, abs=FALSE) {
  ##cc <- paste0("ME",res$net$color)
  cc <- wgcna.labels2colors(res$net$color)
  if("stats" %in% names(res)) {
    traitSignificance <- res$stats$traitSignificance
  } else {
    traitSignificance = as.data.frame(cor(res$datExpr, res$datTraits, use = "p"))
    names(traitSignificance) = names(res$datTraits)
    rownames(traitSignificance) <- colnames(res$datExpr)
  }
  geneSig <- traitSignificance[,trait]
  if(is.null(main)) main <- paste("Module significance with",trait)
  if(abs) geneSig <- abs(geneSig) 
  WGCNA::plotModuleSignificance(
    geneSig, colors=cc, main = main, boxplot=FALSE)
}

#' @export
wgcna.plotSampleDendroAndColors <- function(res,
                                            what=c("me","traits","both")[3],
                                            main=NULL, justdata=FALSE) {

  MET0 <- res$net$MEs
  MET <- MET0[,0]
  if(any(what %in% c("me","both"))) {
     MET <- cbind(MET, MET0)
  }
  if(any(what %in% c("traits","both"))) {
    MET <- cbind(MET, res$datTraits)
  }
  if (NCOL(MET) <= 2) MET <- cbind(MET, MET) ## error if ncol(MET)<=2 !!!!  
  
  ## Recalculate MEs with color as labels
  sampleTree <- hclust(dist(scale(MET0)))  
  ii <- sampleTree$order
  jj <- hclust(dist(t(scale(MET))))$order  
  colors <- WGCNA::numbers2colors(MET[,jj])

  if(justdata) {
    return(MET)
  }
  
  if(is.null(main)) {
    if(what=="me") main <- "Sample dendrogram and module heatmap"
    if(what=="traits") main <- "Sample dendrogram and trait heatmap"
    if(what=="both") main <- "Sample dendrogram and module+traits heatmap"    
  }

  ## Plot the dendrogram and the module colors underneath
  plotDendroAndColors(
    dendro = sampleTree,
    colors = colors,
    groupLabels = colnames(MET)[jj],
    dendroLabels = rownames(MET),
    hang = 0.03,
    addGuide = FALSE,
    guideHang = 0.05,
    marAll = c(0.2, 7, 1.5, 0.5),
    main = main
  )
}


#' @export
wgcna.plotLabeledCorrelationHeatmap <- function(R, nSamples, setpar=TRUE,
                                                cluster=FALSE,
                                                main = NULL, justdata=FALSE,
                                                pstar = TRUE) {
  
  ## Define numbers of genes and samples
  Pvalue = corPvalueStudent(R, nSamples);
  if(cluster) {
    ii <- hclust(dist(R))$order
    jj <- hclust(dist(t(R)))$order
  } else {
    ii <- rev(order(rownames(R)))
    jj <- order(colnames(R))  
  }
  R <- R[ii,jj]
  Pvalue <- Pvalue[ii,jj]    

  if(justdata) return(R)
  
  ## Will display correlations and their p-values
  if(pstar) {
    textPv = cut(Pvalue, breaks=c(-1,0.001,0.01,0.05,99),
                 labels = c("★★★","★★","★",""))
  } else {
    textPv = paste0("(",signif(Pvalue, 1), ")")
  }
  
  textMatrix = paste0(signif(R, 2), "\n", textPv)
  dim(textMatrix) = dim(R)
  if(setpar) par(mar = c(8, 8, 3, 3));
  if(is.null(main)) main <- "Correlation heatmap"
  ## Display the correlation values within a heatmap plot
  WGCNA::labeledHeatmap(Matrix = R,
                        xLabels = colnames(R),
                        yLabels = rownames(R),
                        xSymbols = colnames(R),
                        ySymbols = rownames(R),                        
                        colorLabels = FALSE,
                        colors = blueWhiteRed(50),
                        textMatrix = textMatrix,
                        setStdMargins = FALSE,
                        cex.text = 0.7,
                        zlim = c(-1,1),
                        main = main)

}



#'
#' @export
notWGCNA.compute <- function(X, Y, samples, ntop, nx, ny,
                             xlabel='X', ylabel='Y') {

  X <- head(X[order(-matrixStats::rowSds(X)),],ntop)
  Y <- head(Y[order(-matrixStats::rowSds(Y)),],ntop)
  X <- t(scale(t(X)))
  Y <- t(scale(t(Y)))

  hx <- hclust(dist(X))
  hy <- hclust(dist(Y))
  idx <- cutree(hx, nx)
  idy <- cutree(hy, ny)
  idx <- paste0(xlabel,".",idx)
  idy <- paste0(ylabel,".",idy)  
  
  mX <- rowmean( X, idx )
  mY <- rowmean( Y, idy )

  Z <- expandPhenoMatrix(samples, drop.ref=FALSE)
  P <- cor(t(mX), Z)
  Q <- cor(t(mY), Z)
  R <- cor(t(mX), t(mY))
    
  xgenes <- tapply(1:nrow(X), idx, function(i) rownames(X)[i])
  ygenes <- tapply(1:nrow(Y), idy, function(i) rownames(Y)[i])
  
  res <- list(
    mX = mX,
    mY = mY,    
    Z = Z,
    P = P,
    Q = Q,    
    R = R,
    xgenes = xgenes,
    ygenes = ygenes
  )
  return(res)
}

notWGCNA.condition_on_phenotype <- function(res, pheno, k=2) {
  p.wt <- ((1 + res$P[,pheno])/2)**k
  q.wt <- ((1 + res$Q[,pheno])/2)**k
  outer(p.wt, q.wt)**k 
}
