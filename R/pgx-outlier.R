##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

#' @export
detectOutlierSamples <- function(X, plot = TRUE, par=NULL) {
  ## correlation and distance
  ## X <- playbase::logCPM(playbase::COUNTS)
  ## X <- safe.logCPM(2**X) ## quick normalization
  X <- head(X[order(-matrixStats::rowSds(X, na.rm = TRUE)), ], 4000)
  X <- X - median(X, na.rm = TRUE)
  corX <- cor(X, use = "pairwise")
  distX <- as.matrix(dist(t(X)))

  ## z-score based on correlation
  cor.min <- apply(abs(corX), 1, min, na.rm = TRUE)
  cor.max <- apply(abs(corX), 1, max, na.rm = TRUE)
  cor.median <- apply(abs(corX), 1, median, na.rm = TRUE)
  cor.q10 <- apply(abs(corX), 1, quantile, probs = 0.1, na.rm = TRUE)
  x1 <- (cor.median - mean(cor.median, na.rm = TRUE))
  z1 <- abs(x1 - median(x1, na.rm = TRUE)) / mad(x1, na.rm = TRUE)
  z1

  ## z-score based on euclidean distance
  dist.max <- apply(distX, 1, max, na.rm = TRUE)
  dist.median <- apply(distX, 1, median, na.rm = TRUE)
  dist.q90 <- apply(distX, 1, quantile, probs = 0.9, na.rm = TRUE)
  dist.q10 <- apply(distX, 1, quantile, probs = 0.1, na.rm = TRUE)
  dist.r <- dist.q10 / dist.max
  z2 <- abs(dist.r - median(dist.r, na.rm = TRUE)) / mad(dist.r, na.rm = TRUE)
  z2

  ## gene-wise z-score
  xz <- abs(X - rowMeans(X, na.rm = TRUE)) / matrixStats::rowSds(X, na.rm = TRUE)
  xz <- colMeans(xz, na.rm = TRUE)
  z3 <- abs(xz - median(xz, na.rm = TRUE)) / mad(xz, na.rm = TRUE)

  Z <- cbind(z1, z2, z3)
  colnames(Z) <- c("z.correlation", "z.distance", "z.features")
  zz <- rowMeans(Z, na.rm = TRUE)
  z0 <- 0.1 * mean(Z, na.rm = TRUE)
  zz2 <- exp(rowMeans(log(Z + z0), na.rm = TRUE)) - z0

  res <- list(z.outlier = zz, z.outlier2 = zz2, Z = Z)

  if (plot) {
    plotOutlierScores(res, par = par)
  }

  res
}

#' @export
plotOutlierScores <- function(res.outliers, z.threshold = c(3, 6, 9), par = TRUE) {
  if (par) par(mfrow = c(2, 3), mar=c(8,4,2,2))
  Z <- res.outliers$Z
  zz <- res.outliers$z.outlier
  zz2 <- res.outliers$z.outlier2
  barplot2 <- function(x, ...) {
    barplot(x, ylim = c(0, max(10, max(Z))), ylab = "z-score", ...)
    abline(h = z.threshold, lty = 3, col = "red")
  }
  barplot2(zz, main = "z.outlier (mean)", las = 3)
  barplot2(zz2, main = "z.outlier (geom.mean)", las = 3)
  for (i in 1:ncol(Z)) {
    z1 <- Z[, i]
    barplot2(z1, main = colnames(Z)[i], las = 3)
  }
}
