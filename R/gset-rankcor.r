##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

## ========================================================================
## ================ rank correlation based geneset testing ================
## ========================================================================


#' Title
#'
#' @param rnk value
#' @param gset value
#' @param compute.p value
#'
#' @return
#' @export
#'
#' @examples
gset.cor <- function(rnk, gset, compute.p = FALSE) {
  gset.rankcor(rnk = rnk, gset = gset, compute.p = compute.p, use.rank = FALSE)
}


#' Title
#'
#' @param rnk value
#' @param gset value
#' @param compute.p value
#' @param use.rank value
#'
#' @return
#' @export
#'
#' @examples
gset.rankcor <- function(rnk, gset, compute.p = FALSE, use.rank = TRUE) {
  if (!any(class(gset) %in% c("Matrix", "dgCMatrix", "matrix", "array"))) {
    stop("gset must be a matrix")

  }
  is.vec <- (NCOL(rnk) == 1 && !class(rnk) %in% c("matrix", "array", "Matrix"))
  if (is.vec && is.null(names(rnk))) stop("rank vector must be named")
  if (!is.vec && is.null(rownames(rnk))) stop("rank matrix must have rownames")
  if (is.vec) rnk <- matrix(rnk, ncol = 1, dimnames = list(names(rnk), "rnk"))
  n1 <- sum(rownames(rnk) %in% colnames(gset), na.rm = TRUE)
  n2 <- sum(rownames(rnk) %in% rownames(gset), na.rm = TRUE)
  if (n1 > n2) gset <- t(gset)

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
    rho1 <- parallel::mclapply(1:ncol(rnk1), function(i) rankcorSparse.vec(gset, rnk1[, i]))
    rho1 <- do.call(cbind, rho1)
  }
  rownames(rho1) <- colnames(gset)
  colnames(rho1) <- colnames(rnk1)
  rho1[is.nan(rho1)] <- NA ## ??


  ## compute p-value by permutation
  pv <- qv <- NULL
  if (compute.p) {
    message("computing permutation-based p-values...")
    ## Permutation of the GSET matrix
    idx <- Matrix::which(gset != 0, arr.ind = TRUE)
    slist <- list()
    if (ncol(gset) < 1000) {
      nk <- ceiling(1000 / ncol(gset))
      i <- 1
      for (i in 1:nk) {
        ii <- sample(idx[, 1])
        jj <- sample(idx[, 2])
        slist[[i]] <- Matrix::sparseMatrix(ii, jj,
          x = rep(1, nrow(idx)),
          dims = dim(gset), dimnames = dimnames(gset)
        )
      }
      S <- do.call(cbind, slist)
    }

    rho2 <- gset.rankcor(rnk1, S, compute.p = FALSE, use.rank = use.rank)$rho ## warning:: do not recurse!!
    t.abs.rho2 <- t(abs(rho2))
    pv <- t(apply(rho1, 1, function(x) rowMeans(t.abs.rho2 > abs(x), na.rm = TRUE)))
    pv[is.nan(pv)] <- NA ## ??
    qv <- apply(pv, 2, p.adjust, method = "fdr")
    df <- list(rho = rho1, p.value = pv, q.value = qv)
  } else {
    df <- list(rho = rho1, p.value = NA, q.value = NA)
  }
  df
}
