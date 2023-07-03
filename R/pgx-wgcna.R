##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##
##

#' @export
labels2rainbow <- function(net) {
  hc <- net$dendrograms[[1]]
  nc <- length(unique(net$colors))
  n <- length(net$colors)
  ii <- hc$order
  col1 <- labels2colors(net$colors)
  col.rnk <- rank(tapply(1:n, col1[ii], mean))
  new.col <- rainbow(nc)[col.rnk]
  #
  names(new.col) <- names(col.rnk)
  new.col["grey"] <- "#AAAAAA"
  new.col
  new.col <- new.col[col1]
  names(new.col) <- net$colors
  new.col
}

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
  names(net)
  table(net$colors)

  ## clean up traits matrix
  datTraits <- pgx$samples
  ## no dates please...
  isdate <- apply(datTraits, 2, playbase::is.Date)
  datTraits <- datTraits[, !isdate, drop = FALSE]

  ## Expand multi-class discrete phenotypes into binary vectors
  #
  tr.class <- sapply(type.convert(datTraits), class)
  sel1 <- which(tr.class %in% c("factor", "character"))
  sel2 <- which(tr.class %in% c("integer", "numeric"))

  tr1 <- datTraits[, 0]
  if (length(sel1)) {
    tr1 <- playbase::expandPhenoMatrix(datTraits[, sel1, drop = FALSE], drop.ref = FALSE)
  }
  ## keeping numeric phenotypes
  tr2 <- datTraits[, sel2, drop = FALSE]
  datTraits <- cbind(tr1, tr2)

  ## get colors of eigengene modules
  me.genes <- tapply(names(net$colors), net$colors, list)
  names(me.genes) <- paste0("ME", names(me.genes))
  color1 <- playbase::labels2rainbow(net)
  me.colors <- color1[!duplicated(color1)]
  names(me.colors) <- paste0("ME", names(me.colors))
  me.colors <- me.colors[names(me.genes)]
  

  
  

  X1 <- t(datExpr)
  X1 <- t(scale(datExpr))
  #
  #
  dissTOM <- 1 - WGCNA::TOMsimilarityFromExpr(datExpr, power = power)
  rownames(dissTOM) <- colnames(dissTOM) <- colnames(datExpr)
  clust <- playbase::pgx.clusterBigMatrix(dissTOM, methods = c("umap", "tsne", "pca"), dims = c(2))
  #
  #
  names(clust)
  if ("cluster.genes" %in% names(pgx)) {
    clust[["umap2d"]] <- pgx$cluster.genes$pos[["umap2d"]][colnames(datExpr), ]
  }
  

  
  

  gmt <- getGSETS_playbase(pattern = "HALLMARK|GOBP|^C[1-9]")
  gse <- NULL
  #
  bg <- toupper(rownames(pgx$X))
  i <- 1
  for (i in 1:length(me.genes)) {
    gg <- toupper(me.genes[[i]])
    rr <- playbase::gset.fisher(gg, gmt, background = bg, fdr = 1)
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
