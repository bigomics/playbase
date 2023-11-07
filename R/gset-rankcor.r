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
  if (!any(class(gset) %in% c("Matrix", "dgCMatrix", "matrix", "array"))) {
    stop("gset must be a matrix")
  }
  is.vec <- (NCOL(rnk) == 1 && !class(rnk) %in% c("matrix", "Matrix"))
  if (is.vec && is.null(names(rnk))) stop("rank vector must be named")
  if (!is.vec && is.null(rownames(rnk))) stop("rank matrix must have rownames")
  if (is.vec) rnk <- matrix(rnk, ncol = 1, dimnames = list(names(rnk), "rnk"))
  n1 <- sum(rownames(rnk) %in% colnames(gset), na.rm = TRUE)
  n2 <- sum(rownames(rnk) %in% rownames(gset), na.rm = TRUE)
  if (n1 > n2) gset <- Matrix::t(gset)

  gg <- intersect(rownames(gset), rownames(rnk))
  rnk1 <- rnk[gg, , drop = FALSE]
  gset <- gset[gg, ]

  if (use.rank) {
    rnk1 <- apply(rnk1, 2, base::rank, na.last = "keep")
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
      qlcMatrix::corSparse(X[gg, ], cbind(y[gg]))
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
