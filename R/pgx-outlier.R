##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

#' @export
detectOutlierSamples <- function(X,
                                 methods = c("z.correlation", "z.distance",
                                   "z.features", "z.isoforest")[1:3],
                                 col = "grey70", plot = TRUE, par = TRUE) {

  all.methods <- c("z.correlation", "z.distance", "z.features", "z.isoforest")
  if (is.null(methods)) methods <- all.methods
  ## errors on a typo. without this an unmatched name silently yields a
  ## zero-column Z and NaN z-scores, which read as "no outliers".
  methods <- match.arg(methods, all.methods, several.ok = TRUE)

  ## correlation and distance
  X <- head(X[order(-matrixStats::rowSds(X, na.rm = TRUE)), ], 1000)
  X <- X - median(X, na.rm = TRUE)
  #corX <- HiClimR::fastCor(X, optBLAS = TRUE)
  corX <- cor(X, use="pairwise.complete")  
  distX <- as.matrix(dist(t(X)))

  z1=z2=z3=z4=NULL

  ## robust z-score. mad() is 0 as soon as over half the values tie (e.g. a
  ## sample uploaded twice), which used to give 0/0 = NaN for every sample and
  ## crash the caller. Fall back to the mean absolute deviation, and to zeros
  ## when every value really is identical (nothing deviates, so nothing is an
  ## outlier). Unchanged whenever mad() > 0.
  zscore <- function(x) {
    dx <- abs(x - median(x, na.rm = TRUE))
    s <- mad(x, na.rm = TRUE)
    if (s == 0) s <- mean(dx, na.rm = TRUE) / 0.7979 ## consistency constant
    if (s == 0) return(dx * 0)
    dx / s
  }

  ## z-score based on correlation
  cor.median <- apply(abs(corX), 1, median, na.rm = TRUE)
  x1 <- (cor.median - mean(cor.median, na.rm = TRUE))
  z1 <- zscore(x1)

  ## z-score based on euclidean distance
  dist.max <- apply(distX, 1, max, na.rm = TRUE)
  dist.q10 <- apply(distX, 1, quantile, probs = 0.1, na.rm = TRUE)
  dist.r <- dist.q10 / dist.max
  z2 <- zscore(dist.r)

  ## gene-wise z-score
  xz <- abs(X - rowMeans(X, na.rm = TRUE)) / matrixStats::rowSds(X, na.rm = TRUE)
  xz <- colMeans(xz, na.rm = TRUE)
  z3 <- zscore(xz)

  ## isoforest z-score. only on request: fitting 10k trees is expensive
  if ("z.isoforest" %in% methods) {
    z4 <- outlier.isoforest_zscore(X, ndim=2, ntrees=10000)
  }

  ## NULL columns are dropped by cbind(), so unrequested methods vanish here
  Z <- cbind(z.correlation = z1, z.distance = z2, z.features = z3,
    z.isoforest = z4)
  Z <- Z[, which(colnames(Z) %in% methods), drop = FALSE]

  zz <- rowMeans(Z, na.rm = TRUE)
  z0 <- 0.1 * mean(Z, na.rm = TRUE)
  zz2 <- exp(rowMeans(log(Z + z0), na.rm = TRUE)) - z0

  res <- list(z.outlier = zz, z.outlier2 = zz2, Z = Z)
  if (plot) plotOutlierScores(res, par = par, col=col)
  return(res)
}

#' @export
plotOutlierScores <- function(res.outliers, z.threshold = c(3, 6, 9),
                              col = "grey70", par = TRUE) {
  if (par) {
    ## restore the caller's layout, we are a guest in their device
    opar <- graphics::par(mfrow = c(2, 3), mar = c(8, 4, 2, 2))
    on.exit(graphics::par(opar))
  }
  Z <- res.outliers$Z
  zz <- res.outliers$z.outlier
  zz2 <- res.outliers$z.outlier2
  barplot2 <- function(x, ...) {
    barplot(x, col = col, las = 3,
      ## finite only: an all-NA/Inf column gave "need finite 'ylim' values".
      ## max(10, numeric(0)) is 10, so a fully non-finite Z still plots.
      ylim = c(0, max(10, Z[is.finite(Z)])),
      ylab = "z-score", ...)
    abline(h = z.threshold, lty = 3, col = "red")
  }
  barplot2(zz, main = "z.outlier (mean)")
  barplot2(zz2, main = "z.outlier (geom.mean)")
  for (i in seq_len(ncol(Z))) {
    z1 <- Z[, i]
    barplot2(z1, main = colnames(Z)[i])
  }
}

#' 
outlier.isoforest_zscore <- function(X, ndim=2, ntrees=10000) {
  cX <- X - rowMeans(X, na.rm=TRUE)
  cX <- cX[complete.cases(cX),,drop=FALSE]
  ## ponytail: plain svd. callers cap cX at 1000 rows so a truncated
  ## solver buys nothing, and irlba(nv=3) aborted below 4 samples.
  V <- svd(cX, nu = 0, nv = min(3, ncol(cX)))$v
  ndim <- min(ndim, ncol(V))
  model <- isotree::isolation.forest(V, ndim=ndim, ntrees=ntrees, nthreads=8)
  scores <- predict(model, V)
  scores <- abs(scores - mean(scores, na.rm=TRUE))
  scores <- scores / sd(scores, na.rm=TRUE)  
  scores 
}
