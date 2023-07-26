##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##
##


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
  ii <- hc$order
  col1 <- labels2colors(net$colors)
  col.rnk <- rank(tapply(1:n, col1[ii], mean))
  new.col <- rainbow(nc)[col.rnk]
  names(new.col) <- names(col.rnk)
  new.col["grey"] <- "#AAAAAA"
  new.col <- new.col[col1]
  names(new.col) <- net$colors
  return(new.col)
}


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
    minmodsize = 30,
    power = 6,
    cutheight = 0.25,
    deepsplit = 2,
    ngenes = 1000) {
  X <- as.matrix(pgx$X)
  X <- X[order(-apply(X, 1, sd, na.rm = TRUE)), ]
  X <- X[!duplicated(rownames(X)), ]

  datExpr <- t(head(X, ngenes))
  require(WGCNA) # fun fact: if we dont source WGCNA, blockwiseModules does not work
  net <- WGCNA::blockwiseModules(
    datExpr,
    power = power,
    TOMType = "unsigned", minModuleSize = minmodsize,
    reassignThreshold = 0, mergeCutHeight = cutheight,
    numericLabels = TRUE, pamRespectsDendro = FALSE,
    deepSplit = deepsplit,
    ## saveTOMs = TRUE, saveTOMFileBase = "WCNA.tom",
    verbose = 3
  )

  ## clean up traits matrix
  datTraits <- pgx$samples
  ## no dates please...
  isdate <- apply(datTraits, 2, is.Date)
  datTraits <- datTraits[, !isdate, drop = FALSE]

  ## Expand multi-class discrete phenotypes into binary vectors
  tr.class <- sapply(type.convert(datTraits), class)
  sel1 <- which(tr.class %in% c("factor", "character"))
  sel2 <- which(tr.class %in% c("integer", "numeric"))

  tr1 <- datTraits[, 0]
  if (length(sel1)) {
    tr1 <- expandPhenoMatrix(datTraits[, sel1, drop = FALSE], drop.ref = FALSE)
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

  X1 <- t(datExpr)
  X1 <- t(scale(datExpr))
  dissTOM <- 1 - WGCNA::TOMsimilarityFromExpr(datExpr, power = power)
  rownames(dissTOM) <- colnames(dissTOM) <- colnames(datExpr)
  clust <- pgx.clusterBigMatrix(dissTOM, methods = c("umap", "tsne", "pca"), dims = c(2))
  if ("cluster.genes" %in% names(pgx)) {
    clust[["umap2d"]] <- pgx$cluster.genes$pos[["umap2d"]][colnames(datExpr), ]
  }

  gmt <- getGSETS_playbase(pattern = "HALLMARK|GOBP|^C[1-9]")
  gse <- NULL
  bg <- toupper(rownames(pgx$X))
  i <- 1
  for (i in 1:length(me.genes)) {
    gg <- toupper(me.genes[[i]])
    rr <- gset.fisher(gg, gmt, background = bg, fdr = 1)
    rr <- cbind(
      module = names(me.genes)[i],
      geneset = rownames(rr), rr
    )
    rr <- rr[order(rr$p.value), , drop = FALSE]
    if (i == 1) gse <- rr
    if (i > 1) gse <- rbind(gse, rr)
  }
  rownames(gse) <- NULL

  ## construct results object
  return(
    list(
      datExpr = datExpr,
      datTraits = datTraits,
      net = net,
      gse = gse,
      clust = clust,
      me.genes = me.genes,
      me.colors = me.colors
    )
  )
}
