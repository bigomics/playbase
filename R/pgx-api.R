##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##



#' Title
#'
#' @param pgx value
#'
#' @return
#' @export
#'
#' @examples
pgx.getModelGroups <- function(pgx) {
  model <- pgx$model.parameters$design
  if (is.null(model)) {
    group <- rownames(pgx$model.parameters$exp.matrix)
  } else {
    group <- colnames(model)[max.col(model)]
  }
  group
}

#' Title
#'
#' @param pgx value
#' @param methods value
#' @param level value
#'
#' @return
#' @export
#'
#' @examples
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
      }) ## maxQ
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


#' Title
#'
#' @param pgx value
#' @param what value
#' @param level value
#'
#' @return
#' @export
#'
#' @examples
pgx.getMetaFoldChangeMatrix <- function(pgx, what = "meta", level = "gene") {
  pgx.getMetaMatrix(pgx, methods = what, level = level)
}


#' Title
#'
#' @param pgx value
#'
#' @return
#' @export
#'
#' @examples
pgx.getContrasts <- function(pgx) {
  names(pgx$gx.meta$meta)
}



#' Title
#'
#' @param pgx value
#' @param n value
#' @param ng value
#' @param dir value
#' @param sym value
#' @param filt value
#'
#' @return
#' @export
#'
#' @examples
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
  list(gsets = topgs, genes = gs.genes, top.genes = top.genes)
}


#' Title
#'
#' @param pgx value
#' @param n value
#' @param dir value
#' @param sym value
#' @param filt value
#'
#' @return
#' @export
#'
#' @examples
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
  markers
}


#' Title
#'
#' @param pgx value
#' @param nmin value
#' @param extended value
#'
#' @return
#' @export
#'
#' @examples
pgx.getFamilies <- function(pgx, nmin = 10, extended = FALSE) {
  if (extended) {

    fam.pattern <- "^[<].*|^FAMILY|^COMPARTMENT|^CUSTOM"
  } else {
    fam.pattern <- "^[<].*|^FAMILY|^CUSTOM"
  }
  xgenes <- toupper(rownames(pgx$X))
  xgenes <- toupper(pgx$genes$gene_name)
  gsets <- getGSETS_playbase(pattern = fam.pattern)
  jj <- which(sapply(gsets, function(x) sum(x %in% xgenes)) >= nmin)
  sort(names(gsets)[jj]) ## sort alphabetically
}
