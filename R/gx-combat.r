##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

########################################################################
##
## Functions for batch correction
##
########################################################################

#' Nearest neighbor matching batch correction
#'
#' Correct for batch effects in a gene expression matrix using nearest neighbor matching.
#'
#' @param X Numeric matrix of gene expression values (genes in rows, samples in columns).
#' @param y Factor vector indicating batch for each sample.
#' @param use.design Logical for whether to include batch design matrix in limma correction.
#' @param dist.method Distance metric to use for matching ('cor' or 'euclidean').
#' @param center.x Logical for whether to center gene expression by row means.
#' @param center.m Logical for whether to center expression by batch means.
#' @param sdtop Number of top variable genes to use for correlation.
#' @param replace Logical for whether to replace missing pairs.
#'
#' @return List containing:
#' \itemize{
#'   \item X - Batch corrected gene expression matrix
#'   \item pairings - Matrix of sample pairings used for correction
#' }
#'
#' @details This function performs batch correction using the following steps:
#' \enumerate{
#'   \item Compute distance matrix between all samples
#'   \item Find nearest neighbor matches between batches
#'   \item Construct full paired dataset with matches
#'   \item Apply limma batch correction to paired data
#'   \item Average paired samples back to original samples
#' }
#'
#' It allows replacing missing pairs and using batch design matrices.
#'
#' @seealso
#' \code{\link[limma]{removeBatchEffect}} for the batch correction method used.
#'
#' @examples
#' # TODO
#'
#' @export
gx.nnmcorrect <- function(X, y, use.design = TRUE, dist.method = "cor",
                          center.x = TRUE, center.m = TRUE, sdtop = 1000,
                          replace = FALSE) {
  ## Nearest-neighbour matching for batch correction. This
  ## implementation creates a fully paired dataset with nearest
  ## matching neighbours when pairs are missing.

  ## compute distance matrix for NNM-pairing

  y1 <- paste0("y=", y)
  dX <- X

  if (center.x) {
    dX <- dX - rowMeans(dX, na.rm = TRUE)
  }
  if (center.m) {
    ## center per condition group (takes out batch differences)
    mX <- tapply(1:ncol(dX), y1, function(i) rowMeans(dX[, i, drop = FALSE]))
    mX <- do.call(cbind, mX)
    dX <- dX - mX[, y1]
  }

  if (dist.method == "cor") {
    message("[gx.nnmcorrect] computing correlation matrix D...")
    sdx <- apply(dX, 1, stats::sd)
    ii <- Matrix::head(order(-sdx), sdtop)

    D <- 1 - crossprod(scale(dX[ii, ])) / (length(ii) - 1) ## faster
  } else {
    message("[gx.nnmcorrect] computing distance matrix D...\n")
    D <- as.matrix(stats::dist(t(dX)))
  }
  remove(dX)
  D[is.na(D)] <- 0 ## might have NA

  ## find neighbours
  message("[gx.nnmcorrect] finding nearest neighbours...")
  B <- t(apply(D, 1, function(r) tapply(r, y1, function(s) names(which.min(s)))))
  rownames(B) <- colnames(X)
  Matrix::head(B)

  ## imputing full paired data set
  kk <- match(as.vector(B), rownames(B))
  full.y <- y1[kk]
  full.pairs <- rep(rownames(B), ncol(B))
  full.X <- X[, kk]
  dim(full.X)

  ## remove pairing effect
  message("[gx.nnmcorrect] remove pairing effect...")
  if (use.design) {
    design <- stats::model.matrix(~full.y)
    full.X <- limma::removeBatchEffect(full.X,
      batch = full.pairs,
      design = design
    )
  } else {
    full.X <- limma::removeBatchEffect(full.X, batch = full.pairs)
  }

  ## now contract to original samples
  message("[gx.nnmcorrect] matching result...")
  full.idx <- rownames(B)[kk]
  cX <- do.call(cbind, tapply(
    1:ncol(full.X), full.idx,
    function(i) rowMeans(full.X[, i, drop = FALSE])
  ))
  cX <- cX[, colnames(X)]

  res <- list(X = cX, pairings = B)
  return(res)
}

#' Nearest neighbor matching batch correction
#'
#' Correct for batch effects in a gene expression matrix using nearest neighbor matching.
#'
#' @param x Numeric matrix of gene expression values (genes in rows, samples in columns).
#' @param y Factor vector indicating batch for each sample.
#' @param k Number of nearest neighbors to use (default 3).
#'
#' @return Batch corrected gene expression matrix.
#'
#' @details This function performs batch correction by matching each sample to its
#' k-nearest neighbors from other batches based on expression profile correlation.
#' The batch effect for a sample is estimated as the mean difference between it and its
#' matched neighbors. This difference is subtracted from the sample's expression profile.
#'
#' @seealso
#' \code{\link[limma]{removeBatchEffect}} for an alternative batch correction method
#'
#' @examples
#' \dontrun{
#' x <- matrix(rnorm(100 * 30), 100, 30) # random expression matrix
#' y <- gl(3, 10) # 3 batches of 10 samples each
#' xcorr <- gx.nnmcorrect.SAVE(x, y)
#' }
#' @export
gx.nnmcorrect.SAVE <- function(x, y, k = 3) {
  ## -----------------------------------------------------
  ## nearest-neighbour matching for batch correction
  ## -----------------------------------------------------
  xcor <- stats::cor(x)
  diag(xcor) <- 0
  nx <- x
  j <- 1
  for (j in 1:ncol(x)) {
    nj <- which(y != y[j])
    nn <- intersect(order(-xcor[j, ]), nj)
    nn <- Matrix::head(nn, k)
    nx[, j] <- x[, j] - rowMeans(x[, nn, drop = FALSE])
  }
  nx <- nx + rowMeans(x)
  return(nx)
}


#' Quantile normalize gene expression matrix
#'
#' This function performs quantile normalization on a gene expression matrix.
#'
#' @param X Gene expression matrix with genes in rows and samples in columns.
#'
#' @return The quantile normalized gene expression matrix.
#'
#' @details
#' Quantile normalization is a method to make the distribution of gene expression values the same across all samples.
#' It transforms the values so that the quantiles of each column match the quantiles of a reference distribution.
#'
#' This function first checks if the input matrix contains log2-transformed values by looking at the range of values.
#' If max < 40 and min > 0, it assumes log2 transformed input and applies quantile normalization on the 2^X matrix.
#' Otherwise, it directly applies quantile normalization on X.
#'
#' The normalize.qspline function from the affy package is used internally to perform quantile normalization.
#'
#' @seealso
#' \code{\link[affy]{normalize.qspline}} for the quantile normalization method used.
#'
#' @examples
#' \dontrun{
#' X <- matrix(rnorm(100 * 10), 100, 10)
#' X_norm <- gx.qnormalize(X)
#' }
#' @export

# The following function is not used anywhere in the package or omicsplayground
gx.qnormalize <- function(X) {
  ## -----------------------------------------------------
  ## quantile normalization
  ## -----------------------------------------------------

  if (max(X, na.rm = TRUE) < 40 && min(X, na.rm = TRUE) > 0) {
    # The function normalize.qspline likely comes as product of affy package but better confirm
    tmp <- normalize.qspline(2**X, na.rm = TRUE, verbose = FALSE)
    tmp <- log2(tmp)
  } else {
    tmp <- normalize.qspline(X, na.rm = TRUE, verbose = FALSE)
  }
  rownames(tmp) <- rownames(X)
  colnames(tmp) <- colnames(X)
  X <- tmp
}
