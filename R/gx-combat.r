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
#' @param dist.method Distance metric to use for matching ('cor' or 'euclidean').
#' @param center.x Logical for whether to center gene expression by row means.
#' @param center.m Logical for whether to center expression by batch means.
#' @param sdtop Number of top variable genes to use for correlation.
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
#'
#' @seealso
#' \code{\link[limma]{removeBatchEffect}} for the batch correction method used.
#'
#' @examples
#' # TODO
#'
#' @export
gx.nnmcorrect <- function(X, y, dist.method = "cor",
                          center.x = TRUE, center.m = TRUE, sdtop = 1000) {
  ## Nearest-neighbour matching for batch correction. This
  ## implementation creates a fully paired dataset with nearest
  ## matching neighbours when pairs are missing.

  ## use.design = TRUE;dist.method = "cor";center.x = TRUE;center.m = TRUE; sdtop = 1000

  ## compute distance matrix for NNM-pairing
  y1 <- paste0("y=", y)
  dX <- X

  ## reduce for speed
  sdx <- apply(dX, 1, stats::sd)
  ii <- Matrix::head(order(-sdx), sdtop)
  dX <- dX[ii, ]

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
    ## D <- 1 - crossprod(scale(dX)) / (nrow(dX) - 1) ## faster
    D <- 1 - cor(dX)
  } else {
    message("[gx.nnmcorrect] computing distance matrix D...\n")
    D <- as.matrix(stats::dist(t(dX)))
  }
  ## remove(dX)
  D[is.na(D)] <- 0 ## might have NA

  ## find neighbours
  message("[gx.nnmcorrect] finding nearest neighbours...")
  B <- t(apply(D, 1, function(r) tapply(r, y1, function(s) names(which.min(s)))))
  rownames(B) <- colnames(X)
  Matrix::head(B)

  ## ensure sample is always present in own group
  idx <- cbind(1:nrow(B), match(y1, colnames(B)))
  B[idx] <- rownames(B)

  ## imputing full paired data set
  kk <- match(as.vector(B), rownames(B))
  full.y <- y1[kk]
  full.pairs <- rep(rownames(B), ncol(B))
  full.X <- X[, kk]
  dim(full.X)

  ## remove pairing effect
  message("[gx.nnmcorrect] correcting for pairing effects...")
  use.batch <- TRUE
  if (use.batch) {
    design <- stats::model.matrix(~full.y)
    full.X <- limma::removeBatchEffect(
      full.X,
      batch = full.pairs,
      design = design
    )
  } else {
    V <- model.matrix(~full.pairs)
    design <- stats::model.matrix(~full.y)
    full.X <- limma::removeBatchEffect(
      full.X,
      covariates = V,
      design = design
    )
  }

  ## now contract to original samples
  message("[gx.nnmcorrect] matching result...")
  full.idx <- rownames(B)[kk]
  cX <- do.call(cbind, tapply(
    1:ncol(full.X), full.idx,
    function(i) rowMeans(full.X[, i, drop = FALSE])
  ))
  cX <- cX[, colnames(X)]

  ## retain original row means
  cX <- cX - rowMeans(cX, na.rm = TRUE) + rowMeans(X, na.rm = TRUE)

  res <- list(X = cX, pairings = B)
  return(res)
}


#' @export
gx.nnmcorrect2 <- function(X, y, r = 0.35, use.design = TRUE, dist.method = "cor",
                           center.x = TRUE, center.m = TRUE, mode = "", sdtop = 1000) {
  ## use.design = TRUE;dist.method = "cor";center.x = TRUE;center.m = TRUE; sdtop = 1000;r=0.35

  ## compute distance matrix for NNM-pairing
  y1 <- paste0("y=", y)
  dX <- X

  ## reduce for speed
  sdx <- apply(dX, 1, stats::sd)
  ii <- Matrix::head(order(-sdx), sdtop)
  dX <- dX[ii, ]

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
    message("[gx.nnmcorrect2] computing correlation matrix D...")
    ## D <- 1 - crossprod(scale(dX)) / (nrow(dX) - 1) ## faster
    D <- 1 - cor(dX, use = "pairwise")
  } else {
    message("[gx.nnmcorrect2] computing distance matrix D...\n")
    D <- as.matrix(stats::dist(t(dX)))
  }
  ## remove(dX)
  D[is.na(D)] <- 0 ## might have NA

  ## find neighbours
  message("[gx.nnmcorrect2] finding nearest neighbours...")
  B <- t(apply(D, 1, function(x) tapply(x, y1, function(s) names(which.min(s)))))
  rownames(B) <- colnames(X)

  ## ensure sample is always present in own group
  idx <- cbind(1:nrow(B), match(y1, colnames(B)))
  B[idx] <- rownames(B)

  ## create batch design matrix manually
  idx <- apply(B, 2, function(x) match(x, rownames(B)))
  jj <- as.vector(t(idx))
  ii <- as.vector(mapply(rep, 1:nrow(idx), ncol(idx)))
  P <- Matrix::sparseMatrix(
    i = ii, j = jj, x = rep(1, length(ii)),
    dims = c(nrow(B), nrow(B))
  )

  ## correct for pairing effect
  message("[gx.nnmcorrect2] correcting for pairing effects...")
  P1 <- P
  if (mode == "sym") P1 <- P + Matrix::t(P) ## make symmetric
  if (mode == "tr") P1 <- Matrix::t(P) ## transposed
  if (r < 1) {
    k <- round(min(r * dim(X), dim(X) - 1)) ## critical
    k <- max(k, 1)
    if (r > 0.2) {
      sv <- svd(P1, nu = k, nv = 0)
    } else {
      sv <- irlba::irlba(P1, nu = k, nv = 0)
    }
    P1 <- sv$u
  }

  design <- stats::model.matrix(~y1)
  cX <- limma::removeBatchEffect(X, covariates = scale(P1), design = design)

  ## retain original row means
  cX <- cX - rowMeans(cX, na.rm = TRUE) + rowMeans(X, na.rm = TRUE)

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


