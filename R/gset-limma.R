##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2025 BigOmics Analytics SA. All rights reserved.
##

## ========================================================================
## =============== single-sample based enrichment tests ===================
## ========================================================================

plaid.ttest <- function(X, G, y, ref) {
  if (length(unique(y[!is.na(y)])) > 2) message("[plaid.test] warning: more than 2 classes")
  if (length(unique(y[!is.na(y)])) < 2) stop("[plaid.test] warning: less than 2 classes")
  sel <- which(!is.na(y))
  y1 <- y[sel]
  X1 <- X[, sel, drop = FALSE]
  gsetX <- plaid::plaid(X1, G, nsmooth = 0)
  res <- matrixTests::row_t_welch(
    gsetX[, which(y1 != ref), drop = FALSE],
    gsetX[, which(y1 == ref), drop = FALSE]
  )
  res
}

plaid.limma <- function(X, G, y, ref) {
  if (length(unique(y[!is.na(y)])) > 2) message("[plaid.test] warning: more than 2 classes")
  if (length(unique(y[!is.na(y)])) < 2) stop("[plaid.test] warning: less than 2 classes")
  sel <- which(!is.na(y))
  y1 <- y[sel]
  X1 <- X[, sel, drop = FALSE]
  gsetX <- plaid::plaid(X1, G, nsmooth = 0)
  res <- playbase::gx.limma(gsetX, y, fdr = 1, lfc = 0, ref = ref)
  res
}
