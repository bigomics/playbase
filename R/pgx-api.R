##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


#' Get model groups from PGX object
#'
#' @param pgx A PGX object
#'
#' @return A character vector of model group labels
#'
#' @description Retrieves the model group labels defined in a PGX object.
#'
#' @details This function extracts the model group labels from the design matrix
#' or group labels in a PGX object. It first checks if a model formula is defined in
#' pgx$model.parameters$design. If so, it takes the column names corresponding to the
#' model terms. If not, it simply uses the group labels defined in pgx$model.parameters$exp.matrix.
#'
#' @examples
#' \dontrun{
#' mypgx <- pgx.create(counts, samples, contrast)
#' groups <- pgx.getModelGroups(mypgx)
#' }
#'
#' @export
pgx.getModelGroups <- function(pgx) {
  model <- pgx$model.parameters$design
  if (is.null(model)) {
    group <- rownames(pgx$model.parameters$exp.matrix)
  } else {
    group <- colnames(model)[max.col(model)]
  }
  group
}

#' @title Get meta-analysis results matrix
#'
#' @description
#' Returns a matrix of meta-analysis results for a PGX object.
#'
#' @param pgx A PGX object
#' @param methods Character vector of meta-analysis methods to include.
#' Default is meta.
#' @param level Data level to extract results for ("gene", "geneset").
#'
#' @details
#' This function extracts a matrix of meta-analysis results from a PGX object.
#' It returns fold-changes and adjusted p-values from the specified
#' meta-analysis methods.
#'
#' The results can be extracted at the gene, exon or junction level, depending
#' on the \code{level} argument. By default it returns results for all
#' meta-analysis methods found in the PGX object.
#'
#' A subset of methods can be specified with the \code{methods} argument.
#'
#' @return
#' A matrix with rows corresponding to features at the specified level,
#' and columns for fold-changes and adjusted p-values from each meta-analysis
#' method.
#'
#' @examples
#' \dontrun{
#' mypgx <- pgx.create(counts, samples)
#' meta.res <- pgx.getMetaMatrix(mypgx)
#' head(meta.res)
#' }
#' @export
pgx.getMetaMatrix <- function(pgx, methods = "meta", level = "gene") {
  fc0 <- NULL
  qv0 <- NULL
  if (level == "gene") {
    all.methods <- colnames(unclass(pgx$gx.meta$meta[[1]]$fc))
    all.methods
    if (is.null(methods)) methods <- all.methods
    if (any(methods %in% all.methods)) {
      methods <- intersect(methods, all.methods)
      fc0 <- sapply(pgx$gx.meta$meta, function(x) {
        rowMeans(unclass(x$fc)[, methods, drop = FALSE], na.rm = TRUE)
      })
      qv0 <- sapply(pgx$gx.meta$meta, function(x) {
        apply(unclass(x$q)[, methods, drop = FALSE], 1, max)
      })
      rownames(fc0) <- rownames(qv0) <- rownames(pgx$gx.meta$meta[[1]])
    } else if (methods[1] == "meta") {
      fc0 <- sapply(pgx$gx.meta$meta, function(x) x$meta.fx)
      qv0 <- sapply(pgx$gx.meta$meta, function(x) x$meta.q)
      rownames(fc0) <- rownames(qv0) <- rownames(pgx$gx.meta$meta[[1]])
    } else {
      cat("WARNING:: pgx.getMetaFoldChangeMatrix: unknown method")
      return(NULL)
    }
  }
  if (level == "geneset") {
    all.methods <- colnames(unclass(pgx$gset.meta$meta[[1]]$fc))
    if (is.null(methods)) methods <- all.methods
    if (any(methods %in% all.methods)) {
      fc0 <- sapply(pgx$gset.meta$meta, function(x) {
        rowMeans(unclass(x$fc)[, methods, drop = FALSE], na.rm = TRUE)
      })
      qv0 <- sapply(pgx$gset.meta$meta, function(x) {
        apply(unclass(x$q)[, methods, drop = FALSE], 1, max)
      })
      rownames(fc0) <- rownames(qv0) <- rownames(pgx$gset.meta$meta[[1]])
    } else if (methods[1] == "meta") {
      fc0 <- sapply(pgx$gset.meta$meta, function(x) x$meta.fx)
      qv0 <- sapply(pgx$gset.meta$meta, function(x) x$meta.q)
      rownames(fc0) <- rownames(qv0) <- rownames(pgx$gset.meta$meta[[1]])
    } else {
      cat("WARNING:: pgx.getMetaFoldChangeMatrix: unknown method")
      return(NULL)
    }
  }
  res <- list(fc = fc0, qv = qv0)
  return(res)
}


#' @describeIn pgx.getMetaMatrix Get the Meta-Fold ChangeMatrix from PGX
#' @export
pgx.getMetaFoldChangeMatrix <- function(pgx, what = "meta", level = "gene") {
  pgx.getMetaMatrix(pgx, methods = what, level = level)
}

#' @describeIn pgx.getMetaMatrix Get the contrast from PGX
#' @export
pgx.getContrasts <- function(pgx) {
  names(pgx$gx.meta$meta)
}

#' @describeIn pgx.getMetaMatrix Get the contrast matrix from PGX
#' @export
pgx.getContrastMatrix <- function(pgx) {
  ct.matrix <- pgx$contrasts
  if (is.null(ct.matrix)) {
    ct.matrix <- pgx$model.parameters$exp.matrix
    ct.matrix <- contrastAsLabels(ct.matrix)
  }
  ct.matrix
}

#' @describeIn pgx.getMetaMatrix get the top genesets from PGX
#' @export
pgx.getTopGeneSets <- function(pgx, n = 10, ng = 100, dir = 0, sym = FALSE, filt = NULL) {
  ## Gets top marker genes for all comparisons
  ##
  ##
  F <- pgx.getMarkerGenes(pgx, n = ng, dir = dir, sym = sym)
  G <- sapply(pgx$gset.meta$meta, function(m) m$meta.fx)
  rownames(G) <- rownames(pgx$gset.meta$meta[[1]])
  if (!is.null(filt)) {
    G <- G[grep(filt, rownames(G)), , drop = FALSE]
  }
  if (dir == 0) G <- abs(G)
  if (dir < 0) G <- -1 * (G)
  topgs <- apply(G, 2, function(x) list(names(Matrix::head(sort(-x), n))))
  topgs <- unlist(topgs, recursive = FALSE)
  gs.genes <- lapply(topgs, function(gs) {
    lapply(apply(pgx$GMT[, gs] != 0, 2, which), names)
  })
  top.genes <- gs.genes
  for (i in 1:length(gs.genes)) {
    gg <- lapply(gs.genes[[i]], function(gs) intersect(F[[i]], gs))
    top.genes[[i]] <- gg
  }
  out <- list(gsets = topgs, genes = gs.genes, top.genes = top.genes)

  return(out)
}


#' @title Get Marker Genes
#'
#' @param pgx A PGX object containing results
#' @param n The number of top marker genes to return for each comparison. Default 10.
#' @param dir The direction to sort by. 0 for absolute value, -1 for negative, 1 for positive. Default 0.
#' @param sym Whether to return equal numbers of up and down genes. Default FALSE.
#' @param filt Filter to apply to gene names before selecting markers. Default NULL.
#'
#' @details
#' This function extracts the top n marker genes (based on absolute t-statistic or fold change)
#' for each comparison stored in the PGX object. It returns a list with one element per comparison,
#' containing the top marker genes.
#'
#' The number and direction of sorting can be controlled with n and dir arguments.
#' Setting sym=TRUE returns an equal number of up and down genes.
#' A filter can be applied to gene names before selecting markers.
#'
#' @return A named list with top marker genes for each comparison
#'
#' @describeIn pgx.getMetaMatrix Get Marker Genes
#' @export
pgx.getMarkerGenes <- function(pgx, n = 10, dir = 0, sym = FALSE, filt = NULL) {
  ## Gets top marker genes for all comparisons
  ##
  ##
  F <- sapply(pgx$gx.meta$meta, function(m) m$meta.fx)
  #
  rownames(F) <- rownames(pgx$gx.meta$meta[[1]])
  if (!is.null(filt)) {
    F <- F[grep(filt, rownames(F)), , drop = FALSE]
  }
  if (sym == TRUE) {
    n2 <- round(n / 2)
    markers <- apply(F, 2, function(x) {
      list(c(names(Matrix::head(sort(x), n2)), rev(names(Matrix::tail(sort(x), n2)))))
    })
  } else {
    if (dir == 0) F <- abs(F)
    if (dir < 0) F <- -1 * (F)
    markers <- apply(F, 2, function(x) list(names(Matrix::head(sort(-x), n))))
  }
  markers <- unlist(markers, recursive = FALSE)
  return(markers)
}


#' @title Get Families
#'
#' @param pgx A PGX object containing gene expression data
#' @param nmin The minimum number of genes for a family to be returned
#' @param extended Whether to return extended families (TRUE) or only core families (FALSE)
#'
#' @return A character vector of gene family names
#'
#' @description Retrieves the gene families represented in a PGX object that have at least nmin members.
#'
#' @details Searches for gene families in the PGX gene sets that match either the core family names (FAMILY*) or
#' extended family regex (e.g. <FAMILY.ABC.*). Returns the names of families with at least nmin members
#' present in the PGX gene expression matrix. Useful for checking which gene families are represented in the data.
#'
#' @describeIn pgx.getMetaMatrix Get Families from PGX
#' @export
pgx.getFamilies <- function(pgx, nmin = 10, extended = FALSE) {
  if (extended) {
    fam.pattern <- "^[<].*|^FAMILY|^COMPARTMENT|^CUSTOM"
  } else {
    fam.pattern <- "^[<].*|^FAMILY|^CUSTOM"
  }


  xgenes <- pgx$genes[, "symbol"]
  xgenes <- unique(xgenes)
  gsets <- playbase::getGSETS_playbase(pattern = fam.pattern)
  jj <- which(sapply(gsets, function(x) sum(x %in% xgenes)) >= nmin)
  return(sort(names(gsets)[jj])) ## sort alphabetically
}
