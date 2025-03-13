##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

## ========================================================================
## ================ rank correlation based geneset testing ================
## ========================================================================


#' Calculate gene set correlation
#'
#' Computes correlation between a gene rank vector and gene sets.
#'
#' @param rnk Numeric vector of gene ranks.
#' @param gset List of gene sets. Each element is a character vector of gene names.
#' @param compute.p Logical indicating whether to compute p-values.
#'
#' @return Matrix of correlation values between rnk and each gene set.
#'
#' @details This function calls gset.rankcor() with use.rank=FALSE to calculate
#' Pearson correlation between the gene rank vector and membership in each gene set.
#'
#' P-values can be optionally computed to assess correlation significance.
#'
#' @examples
#' \dontrun{
#' librart(playbase)
#' ranks <- sample(1:10000, 1000, replace = TRUE)
#' names(ranks) <- replicate(1000, paste(sample(LETTERS, 4, replace = TRUE), collapse = ""))
#' genesets <- list(
#'   GS1 = sample(names(ranks), 100),
#'   GS2 = sample(names(ranks), 50)
#' )
#'
#' gset.cor(ranks, genesets)
#' }
#' @export
gset.cor <- function(rnk, gset, compute.p = FALSE) {
  gset.rankcor(rnk = rnk, gset = gset, compute.p = compute.p, use.rank = FALSE)
}

#' Calculate gene set rank correlation
#'
#' Compute rank correlation between a gene rank vector/matrix and gene sets
#'
#' @param rnk Numeric vector or matrix of gene ranks, with genes as row names
#' @param gset Numeric matrix of gene sets, with genes as row/column names
#' @param compute.p Logical indicating whether to compute p-values
#' @param use.rank Logical indicating whether to rank transform rnk before correlation
#'
#' @return Named list with components:
#' \itemize{
#'  \item rho - Matrix of correlation coefficients between rnk and gset
#'  \item p.value - Matrix of p-values for correlation (if compute.p = TRUE)
#'  \item q.value - Matrix of FDR adjusted p-values (if compute.p = TRUE)
#' }
#'
#' @details This function calculates sparse rank correlation between rnk and each
#' column of gset using \code{qlcMatrix::corSparse()}. It handles missing values in
#' rnk by computing column-wise correlations.
#'
#' P-values are computed from statistical distribution
#'
#' @examples
#' \dontrun{
#' librart(playbase)
#' ranks <- sample(1:10000, 1000, replace = TRUE)
#' names(ranks) <- replicate(1000, paste(sample(LETTERS, 4, replace = TRUE), collapse = ""))
#' genesets <- matrix(rnorm(1000 * 20), ncol = 20)
#' rownames(genesets) <- names(ranks)
#'
#' gset.rankcor(ranks, genesets, compute.p = TRUE)
#' }
#' @export
gset.rankcor <- function(rnk, gset, compute.p = FALSE, use.rank = TRUE) {
  if (ncol(gset) == 0 || NCOL(rnk) == 0) {
    if (ncol(gset) == 0) message("ERROR. gset has zero columns")
    if (NCOL(rnk) == 0) message("ERROR. rnk has zero columns")
    return(NULL)
  }

  #  if (!any(class(gset) %in% c("Matrix", "dgCMatrix", "lgCMatrix", "matrix", "array"))) {
  #    stop("gset must be a matrix")
  #  }
  if (!inherits(gset, "Matrix")) stop("gset must be a matrix")

  is.vec <- (NCOL(rnk) == 1 && !any(class(rnk) %in% c("matrix", "Matrix")))
  if (is.vec && is.null(names(rnk))) stop("rank vector must be named")
  if (!is.vec && is.null(rownames(rnk))) stop("rank matrix must have rownames")
  if (is.vec) rnk <- matrix(rnk, ncol = 1, dimnames = list(names(rnk), "rnk"))
  n1 <- sum(rownames(rnk) %in% colnames(gset), na.rm = TRUE)
  n2 <- sum(rownames(rnk) %in% rownames(gset), na.rm = TRUE)
  if (n1 > n2) gset <- Matrix::t(gset)

  gg <- intersect(rownames(gset), rownames(rnk))
  rnk1 <- rnk[gg, , drop = FALSE]
  gset <- gset[gg, , drop = FALSE]

  if (use.rank) {
    ##  rnk1 <- apply(rnk1, 2, base::rank, na.last = "keep", ties.method="random")
    rnk1 <- t(matrixStats::colRanks(rnk1, na.last = "keep", ties.method = "random"))
  }

  ## two cases: (1) in case no missing values, just use corSparse on
  ## whole matrix. (2) in case the rnk matrix has missing values, we
  ## must proceed 1-column at time and do reduced corSparse on
  ## intersection of genes.
  if (sum(is.na(rnk1)) == 0) {
    rho1 <- qlcMatrix::corSparse(gset, rnk1)
  } else {
    message("rank matrix has missing values: computing column-wise reduced rankcor")
    rankcorSparse.vec <- function(X, y) {
      y <- y[!is.na(y)]
      gg <- intersect(rownames(X), names(y))
      y <- rank(y, na.last = "keep")
      qlcMatrix::corSparse(X[gg, , drop = FALSE], cbind(y[gg]))
    }
    rho1 <- lapply(1:ncol(rnk1), function(i) rankcorSparse.vec(gset, rnk1[, i]))
    rho1 <- do.call(cbind, rho1)
  }

  rownames(rho1) <- colnames(gset)
  colnames(rho1) <- colnames(rnk1)
  rho1[is.nan(rho1)] <- NA ## ??

  ## compute p-value
  .cor.pvalue <- function(x, n) 2 * stats::pnorm(-abs(x / ((1 - x**2) / (n - 2))**0.5))
  if (compute.p) {
    pv <- apply(rho1, 2, function(x) .cor.pvalue(x, n = nrow(rnk1)))
    pv[is.nan(pv)] <- NA ## ??
    qv <- apply(pv, 2, stats::p.adjust, method = "fdr")
    df <- list(rho = rho1, p.value = pv, q.value = qv)
  } else {
    df <- list(rho = rho1, p.value = NA, q.value = NA)
  }
  df
}


#' Compute geneset expression as the average log-ration of genes in
#' the geneset. Requires log-expression matrix X and (sparse) geneset
#' matrix matG.
#'
#' @export
gset.averageCLR <- function(X, matG, center = TRUE, use.rank = FALSE) {
  if (NCOL(X) == 1) X <- cbind(X)
  if (center && NCOL(X) > 1) X <- X - rowMeans(X, na.rm = TRUE)
  if (use.rank) X <- colSignedRanks(X) ## playbase
  gg <- intersect(rownames(X), rownames(matG))
  if (length(gg) == 0) {
    message("[gset.averageCLR] WARNING. features of gene matrix X do not overlap with geneset matrix.")
  }
  X <- X[gg, , drop = FALSE]
  matG <- matG[gg, , drop = FALSE]
  ngenes <- Matrix::colSums(matG != 0)
  gsetX <- (Matrix::t(matG != 0) %*% X) / pmax(ngenes, 1e-3)
  as.matrix(gsetX)
}

#' Fast version of average CLR for sparse matrices.
#'
#' @export
avgCLR_sparse <- function(X, matG) {
  sumG <- 1e-8 + Matrix::colSums(matG != 0, na.rm = TRUE)
  nG <- Matrix::colScale(1 * (matG != 0), 1 / sumG)
  mx <- Matrix::rowMeans(X, na.rm = TRUE)
  Matrix::t(nG) %*% X - Matrix::colSums(nG * mx)
}

#' Fast version of average ranked-CLR for sparse matrices.
#'
#' @export
avgRCLR_sparse <- function(X, matG) {
  sumG <- 1e-8 + Matrix::colSums(matG != 0, na.rm = TRUE)
  nG <- Matrix::colScale(1 * (matG != 0), 1 / sumG)
  ## ????
  mx <- Matrix::rowMeans(X, na.rm = TRUE)
  Matrix::t(nG) %*% X - Matrix::colSums(nG * mx)
}



#' Single sample genesets expression scores. See Foroutan 2018. This
#' seems to be a nice and efficient algorithm.
#'
#' @export
gset.singscore <- function(X, geneSets, center = FALSE, return.score = TRUE) {
  gslist <- lapply(geneSets, GSEABase::GeneSet)
  for (i in 1:length(gslist)) GSEABase::setName(gslist[[i]]) <- names(geneSets)[i]
  gscolc <- GSEABase::GeneSetCollection(gslist)
  if (center) {
    X <- X - rowMeans(X, na.rm = TRUE)
  }
  ranked <- singscore::rankGenes(X)
  sing <- singscore::multiScore(ranked, upSetColc = gscolc)
  ## sing has 'Scores' and 'Dispersions'
  if (return.score) {
    return(sing$Scores)
  }
  return(sing)
}

#' Single sample genesets expression scores. See Foroutan 2018. This
#' seems to be a nice and efficient algorithm.
#'
#' @export
gset.gsva <- function(X, geneSets, method = "gsva") {
  if (!method %in% c("gsva", "ssgsea", "plage", "zscore")) {
    stop("ERROR. unknown method", method)
  }
  gsvapar <- switch(method,
    gsva = GSVA::gsvaParam(X, geneSets),
    plage = GSVA::plageParam(X, geneSets),
    zscore = GSVA::zscoreParam(X, geneSets),
    ssgsea = GSVA::ssgseaParam(X, geneSets)
  )
  es <- GSVA::gsva(gsvapar, verbose = FALSE)
  return(es)
}
