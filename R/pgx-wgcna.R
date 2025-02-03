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
pgx.wgcna <- function(pgx,
                      minmodsize = 30,
                      power = 6,
                      cutheight = 0.25,
                      deepsplit = 2,
                      ngenes = 1000) {

  ## if we dont source WGCNA, blockwiseModules does not work..
  require(WGCNA)
  
  ##-------------------------------
  ## WGCNA does not work very well with scRNAseq due to sparsity.
  ## To bypass the issue, hdWGCNA compute metacells to bypass the problem.
  ## Here i compute metacells too using Supercell

  X <- NULL
  samples <- NULL
  if (pgx$datatype %in% c("scRNAseq", "scRNA-seq")) {
    message("[playbase::pgx.wgcna] WGCNA: scRNAseq. Computing supercells with SuperCell.")
    counts <- pgx$counts
    group <- pgx$samples[, "celltype"]
    q10 <- quantile(table(group), probs=0.25)
    if (ncol(counts) > 2000) {
      nb <- round( ncol(counts) / 2000 )
    } else {
      d <- round(ncol(counts)/8, 1)
      nb <- round( ncol(counts) / d) ## temporary...
    }
    message("[pgx.createSingleCellPGX]=======================================")
    message("[pgx.createSingleCellPGX] running SuperCell. nb = ", nb)    
    sc <- pgx.supercell(counts, pgx$samples, group = group, gamma = nb)
    message("[pgx.createSingleCellPGX] SuperCell done: ", ncol(counts), " ->", ncol(sc$counts))
    message("[pgx.createSingleCellPGX]=======================================")
    message("[pgx.createSingleCellPGX] Normalizing supercell matrix (logCPM)")
    X <- playbase::logCPM(sc$counts, total = 1e6, prior = 1)
    X <- as.matrix(X)
    samples <- sc$meta
    remove(counts, group, sc); gc()
  }
  
  if (is.null(X)) X <- as.matrix(pgx$X)
  if (is.null(samples)) samples <- pgx$samples
  
  ## get topSD matrix
  nmissing <- sum(is.na(X))
  if (nmissing > 0) {
    message("Found ", nmissing, " missing values in X. Removing prior to WGCNA.")
    X <- X[complete.cases(X), , drop = FALSE]
  }
  sdx <- matrixStats::rowSds(X, na.rm = TRUE)
  X <- X[sdx > 0.1 * mean(sdx, na.rm = TRUE), , drop = FALSE] ## filter low SD
  X <- X[order(-matrixStats::rowSds(X, na.rm = TRUE)), , drop = FALSE]
  X <- X[!duplicated(rownames(X)), , drop = FALSE]
  X <- utils::head(X, ngenes)
  datExpr <- t(X)

  ## adapt for small datasets
  minmodsize <- 30
  minmodsize <- min(minmodsize, nrow(X) / 10)
  minmodsize

  ## we iterate over smaller power, until we get some clusters
  message("[playbase::pgx.wgcna] Computing WGCNA blockwiseModules...")
  ncolors <- 1
  i <- 1
  while (ncolors == 1 && i < 100) {
    p <- 1 + (power - 1) / i
    net <- WGCNA::blockwiseModules(
      datExpr,
      power = p,
      TOMType = "unsigned",
      minModuleSize = minmodsize,
      reassignThreshold = 0,
      mergeCutHeight = cutheight,
      numericLabels = TRUE,
      pamRespectsDendro = FALSE,
      deepSplit = deepsplit,
      ## saveTOMs = TRUE, saveTOMFileBase = "WCNA.tom",
      verbose = 3
    )
    ncolors <- length(table(net$colors))
    i <- i + 1
  }
  table(net$colors)

  ## clean up traits matrix
  ## datTraits <- pgx$samples
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
  message("[playbase::pgx.wgcna] Computing clustering based on TOM matrix...")
  X1 <- t(datExpr)
  X1 <- t(scale(datExpr))
  dissTOM <- 1 - WGCNA::TOMsimilarityFromExpr(datExpr, power = power)
  rownames(dissTOM) <- colnames(dissTOM) <- colnames(datExpr)
  saveRDS(dissTOM, "~/Desktop/MNT/dissTOM.RDS")
  clust <- pgx.clusterBigMatrix(dissTOM, methods = c("umap", "tsne", "pca"), dims = c(2))
  if ("cluster.genes" %in% names(pgx)) {
    clust[["pca2d"]] <- pgx$cluster.genes$pos[["pca2d"]][colnames(datExpr), ]
    clust[["umap2d"]] <- pgx$cluster.genes$pos[["umap2d"]][colnames(datExpr), ]
    clust[["tsne2d"]] <- pgx$cluster.genes$pos[["tsne2d"]][colnames(datExpr), ]
  }

  ## Do quick geneset analysis using fisher-test (fastest method)
  message("[playbase::pgx.wgcna] WGCNA gene set analysis using Fisher...")
  gmt <- getGSETS_playbase(pattern = "HALLMARK|GOBP|^C[1-9]|GO_BP")
  gse <- NULL
  bg <- rownames(pgx$X)
  bg <- probe2symbol(bg, pgx$genes, query = "human_ortholog")
  bg <- toupper(bg[!is.na(bg)])
  ## bg <- toupper(rownames(pgx$X))

  gmt1 <- lapply(gmt, function(m) intersect(m, bg))
  gmt1.size <- sapply(gmt1, length)
  gmt1 <- gmt1[gmt1.size >= 10]

  i <- 1
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
  message("[playbase::pgx.wgcna] WGCNA completed. Returning object.")
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
