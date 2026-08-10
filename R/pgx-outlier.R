##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

#' @export
detectOutlierSamples <- function(X,
                                 methods = c("z.correlation", "z.distance",
                                   "z.features", "z.isoforest")[1:3],
                                 col = "grey70", plot = TRUE, par = TRUE) {

  if(is.null(methods)) {
    methods = c("z.correlation", "z.distance", "z.features", "z.isoforest")
  }

  ## correlation and distance
  X <- head(X[order(-matrixStats::rowSds(X, na.rm = TRUE)), ], 1000)
  X <- X - median(X, na.rm = TRUE)
  #corX <- HiClimR::fastCor(X, optBLAS = TRUE)
  corX <- cor(X, use="pairwise.complete")  
  distX <- as.matrix(dist(t(X)))

  z1=z2=z3=z4=NULL
  
  ## z-score based on correlation
  cor.median <- apply(abs(corX), 1, median, na.rm = TRUE)
  x1 <- (cor.median - mean(cor.median, na.rm = TRUE))
  z1 <- abs(x1 - median(x1, na.rm = TRUE)) / mad(x1, na.rm = TRUE)

  ## z-score based on euclidean distance
  dist.max <- apply(distX, 1, max, na.rm = TRUE)
  dist.q10 <- apply(distX, 1, quantile, probs = 0.1, na.rm = TRUE)
  dist.r <- dist.q10 / dist.max
  z2 <- abs(dist.r - median(dist.r, na.rm = TRUE)) / mad(dist.r, na.rm = TRUE)

  ## gene-wise z-score
  xz <- abs(X - rowMeans(X, na.rm = TRUE)) / matrixStats::rowSds(X, na.rm = TRUE)
  xz <- colMeans(xz, na.rm = TRUE)
  z3 <- abs(xz - median(xz, na.rm = TRUE)) / mad(xz, na.rm = TRUE)

  ## isoforest z-score. only on request: fitting 10k trees is expensive
  if ("z.isoforest" %in% methods) {
    z4 <- outlier.isoforest_zscore(X, scale=TRUE, ndim=2, ntrees=10000)
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
  if (par) par(mfrow = c(2, 3), mar = c(8, 4, 2, 2))
  Z <- res.outliers$Z
  zz <- res.outliers$z.outlier
  zz2 <- res.outliers$z.outlier2
  barplot2 <- function(x, ...) {
    barplot(x, col = col, las = 3,
      ylim = c(0, max(10, max(Z))),
      ylab = "z-score", ...)
    abline(h = z.threshold, lty = 3, col = "red")
  }
  barplot2(zz, main = "z.outlier (mean)")
  barplot2(zz2, main = "z.outlier (geom.mean)")
  for (i in 1:ncol(Z)) {
    z1 <- Z[, i]
    barplot2(z1, main = colnames(Z)[i])
  }
}

#' 
outlier.isoforest_zscore <- function(X, scale=TRUE, ndim=2, ntrees=10000) {
  cX <- X - rowMeans(X, na.rm=TRUE)
  cX <- cX[complete.cases(cX),,drop=FALSE]
  res <- irlba::irlba(cX, nv=3)
  V <- res$v
  ndim <- min(ndim, min(dim(cX)))
  model <- isotree::isolation.forest(V, ndim=ndim, ntrees=ntrees, nthreads=8)
  scores <- predict(model, V)
  scores <- abs(scores - mean(scores, na.rm=TRUE))
  scores <- scores / sd(scores, na.rm=TRUE)  
  scores 
}
