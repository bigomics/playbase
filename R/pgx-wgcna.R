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

  ##minmodsize = 30;power = 6;cutheight = 0.25;deepsplit = 2;ngenes = 1000;

  nmissing <- sum(is.na(X))
  if (nmissing > 0) {
    message("Found ", nmissing, " missing values in X. Removing prior to WGCNA.")
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
  
  ## construct results object
  results <- list(
      datExpr = datExpr,
      datTraits = datTraits,
      TOM = TOM,
      net = net,
      power = power,
      me.genes = me.genes,
      me.colors = me.colors
  )

  return(results)
}


#'
#'
#' @export
plotTOMfromResults <- function(results, power=NULL, networktype="signed",
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
plotDendroFromResults <- function(results, power=NULL, networktype="signed",
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


#' WMFCNA: weighted multi-omics feature correlation
#' analysis... :). Essentially a WGCNA redo but redone-simpler-faster.
#'
#' 
wmfcna.compute <- function(X, y=NULL, k=20, power=1, vex=1,
                           ewidth = 1, min.cor = 0.7) {

  if(!is.null(y)) {
    Y <- expandPhenoMatrix(data.frame( PHENO = y),drop.ref=FALSE)
    X <- rbind(X, t(Y))
  }

  ## cluster features
  TOM <- cor(t(X),use="pairwise")
  if(power!=1) TOM <- abs(TOM)**power * sign(TOM)
  hc <- fastcluster::hclust(as.dist(1 - TOM), method="complete")

  ## name modules. sort by highest variance contribution
  me.idx <- cutree(hc, k=k)
  cX <- X - rowMeans(X, na.rm=TRUE)
  varx <- rowsum( rowSums(cX**2), me.idx)[,1]
  me.idx <- rank(-varx)[me.idx]  ## sort by variance
  if(k>=10) me.idx <- sprintf("%02d", me.idx)
  me <- paste0("Module",me.idx)
  names(me) <- rownames(TOM)
  table(me)
  
  ##
  tapply(names(me), me, c)
  me.size <- table(me)
  me.size
  table(me.size)
  
  ## cluster-averaged correlation
  datSD <- matrixStats::rowSds(X,na.rm=TRUE)
  stdExpr <- (X - rowMeans(X,na.rm=TRUE)) / datSD
  meExpr <- playbase::rowmean(stdExpr, me)
  meCor <- cor(t(meExpr))
  
  ## Module-averaged correlation
  F <- t(meExpr)  ## scores
  W <- matrix(0, nrow(X), ncol(F)) ## loadings
  idx <- cbind(1:nrow(W), match(me, colnames(F)))
  W[idx] <- 1
  W <- W * datSD
  rownames(W) <- rownames(X)
  colnames(W) <- colnames(F) 

  ## cluster-reduced graph
  create_graph <- function(gx, y, min.cor) {
    rho <- cor(t(gx))      
    gr <- igraph::graph_from_adjacency_matrix(
      rho, weighted=TRUE, diag=FALSE, mode = "undirected")
    
    ew <- igraph::E(gr)$weight
    ew <- (ew - min(ew,na.rm=TRUE)) /
      (max(ew,na.rm=TRUE)-min(ew,na.rm=TRUE))
    igraph::E(gr)$weight <- ew
    gr <- igraph::delete_edges(gr, edges = which(ew < min.cor))
    
    if(!is.null(y)) {
      y <- factor(as.character(y))
      if(length(unique(y[!is.na(y)]))==2) {
        val <- cor(t(gx), as.integer(y))[,1]
      } else {
        res <- gx.limmaF(gx, y, fdr=1, lfc=0)
        val <- res$logFC
        names(val) <- rownames(res)
      }
      ##val <- tapply(val, me, mean)
      val <- val[igraph::V(gr)$name]
      igraph::V(gr)$score <- val
      igraph::V(gr)$color <- colorscale(val, gamma=2)
    }
    igraph::E(gr)$width <- 3 * ewidth * abs(igraph::E(gr)$weight)**2
    ##V(gr)$size <-  5 * vsize*(me.size)**0.5
    ##gr <- subgraph(gr, vids = which(igraph::degree(gr)>1))
    ##gr <- mst(gr, 1/E(gr)$weight)
    gr
  }

  ## member labels
  labels <- tapply(names(me), me, function(s) paste(s,collapse=" "))
  labels <- sapply(sapply(labels, strwrap, 60), paste, collapse="\n")
    
  ## create factor graph
  meGraph <- create_graph(meExpr, y, min.cor )
  vname <- igraph::V(meGraph)$name
  vsize <- me.size[vname]
  vsize[is.na(vsize)] <- median(vsize,na.rm=TRUE)
  vsize <- vsize / mean(vsize)
  igraph::V(meGraph)$size <-  10 * vex * vsize**0.5

  labels <- ifelse( is.na(labels[vname]), vname, labels[vname])
  igraph::V(meGraph)$label <- labels
  igraph::V(meGraph)$title <- labels  

  ## color phenotype nodes
  has.pheno <- grepl("PHENO", labels)
  igraph::V(meGraph)$shape <- ifelse(has.pheno, "star", "dot")
  igraph::V(meGraph)$size[which(has.pheno)] <- 0.9 * max(igraph::V(meGraph)$size)
  
  subgraphs <- list()
  table(me)
  k=me[1]
  modules <- unique(me)
  for(k in modules) {
    sel <- which(me == k)
    if(length(sel)>1) {
      gx <- X[sel,,drop=FALSE]
      gr <- create_graph(gx, y, min.cor )
      is.pheno <- grepl("PHENO", igraph::V(gr)$name)
      igraph::V(gr)$shape <- ifelse(is.pheno,"star","dot")
      igraph::V(gr)$color <- ifelse(is.pheno,"orange",igraph::V(gr)$color)
      val <- abs(igraph::V(gr)$score)
      igraph::V(gr)$size <- 10 * vex * (val / max(val))
      igraph::V(gr)$size[which(is.pheno)] <- 0.9 * max(igraph::V(gr)$size)
      subgraphs[[k]] <- gr
    } else {
      subgraphs[[k]] <- NULL
    }
  }

  members <- tapply(names(me), me, c)
  
  return(
    list(
      X = X,
      Y = Y,
      #F = F,
      W = W,
      ## samples = samples,
      TOM = TOM,
      meCor = meCor,
      meExpr = meExpr,      
      cluster = me,
      members = members,
      graph = meGraph,
      subgraphs = subgraphs 
    )
  )

  
}
