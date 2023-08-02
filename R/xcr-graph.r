##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

## ===================================================================================
## ============================== GRAPH METHODS ======================================
## ===================================================================================




## Calculate edge values from X
#' @export
calc.edge.similarity <- function(ee, X, nk = 4000, mc.cores = 1) {
  if (mc.cores > 1) {
    return(calc.edge.similarityMC(ee, X, nk = nk, mc.cores = mc.cores))
  } else {
    return(calc.edge.similarityKFOLD(ee, X, nk = nk))
  }
}

#' @export
calc.edge.similarityMC <- function(ee, X, nk = 5000, mc.cores = 4) {
  kfold <- (nrow(ee) %/% nk) + 1
  kfold
  idx <- list()
  i <- 1
  for (i in 1:kfold) {
    j0 <- (i - 1) * nk + 1
    j1 <- min(nrow(ee), (i - 1) * nk + nk)
    if (j0 > nrow(ee)) next
    idx[[i]] <- j0:j1
  }
  cx <- rep(NA, nrow(ee))

  compute.cx <- function(ii) {
    ## NEED RETHINK!!! can be optimized bit more
    cx0 <- rep(NA, length(ii))
    which.x <- which(!is.na(X[ee[ii, 1], 1]) & !is.na(X[ee[ii, 2], 1]))
    jj <- ii[which.x]
    x1 <- X[ee[jj, 1], , drop = FALSE]
    x2 <- X[ee[jj, 2], , drop = FALSE]
    n1 <- Matrix::rowSums(x1**2 * (!is.na(x2)), na.rm = TRUE)**0.5
    n2 <- Matrix::rowSums(x2**2 * (!is.na(x1)), na.rm = TRUE)**0.5
    n1 <- (1e-20 + n1)
    n2 <- (1e-20 + n2)
    xx <- x1 * x2
    ## sparseMatrix does not propagate NA properly!!
    cx1 <- Matrix::rowSums(xx, na.rm = TRUE) / n1 / n2
    cx0[which.x] <- cx1
    cx0
  }
  res <- parallel::mclapply(idx, compute.cx, mc.cores = mc.cores)
  cx <- as.numeric(unlist(res))
  return(cx)
}
## Calculate edge values from X
#' @export
calc.edge.similarityKFOLD <- function(ee, X, nk = 4000) {
  calc.similarity <- function(ee, X) {
    cx0 <- rep(NA, nrow(ee))
    jj <- which(!is.na(X[ee[, 1], 1]) & !is.na(X[ee[, 2], 1])) ## check NA
    x1 <- X[ee[jj, 1], , drop = FALSE]
    x2 <- X[ee[jj, 2], , drop = FALSE]
    n1 <- Matrix::rowSums(x1**2 * (!is.na(x2)), na.rm = TRUE)**0.5
    n2 <- Matrix::rowSums(x2**2 * (!is.na(x1)), na.rm = TRUE)**0.5
    n1 <- (1e-20 + n1)
    n2 <- (1e-20 + n2)
    xx <- x1 * x2
    ## sparseMatrix does not propagate NA properly!!
    cx1 <- Matrix::rowSums(xx, na.rm = TRUE) / n1 / n2
    cx0[jj] <- cx1
    return(cx0)
  }
  kfold <- (nrow(ee) %/% nk) + 1
  kfold
  cx <- rep(NA, nrow(ee))
  i <- 1
  for (i in 1:kfold) {
    j0 <- (i - 1) * nk + 1
    j1 <- min(nrow(ee), (i - 1) * nk + nk)
    if (j0 > nrow(ee)) next
    jj <- j0:j1
    ee1 <- ee[jj, ]
    cx0 <- calc.similarity(ee1, X)
    cx[jj] <- cx0
  }
  return(cx)
}




#' @export
graph.cut_crossings <- function(g, idx, max.wt = 9999) {
  ## Cut graph given indices of membership
  ee <- igraph::get.edgelist(g)
  dim(ee)
  jj <- which(idx[ee[, 1]] != idx[ee[, 2]] & igraph::E(g)$weight < max.wt)
  length(jj)
  if (length(jj) > 0) g <- igraph::delete.edges(g, igraph::E(g)[jj])
  return(g)
}

#' @export
itercluster_louvain <- function(g, n = 3) {
  i <- 1
  idx <- rep(1, length(igraph::V(g)))
  K <- c()
  for (i in 1:n) {
    k <- max(idx)
    newidx <- idx
    for (i in 1:k) {
      ii <- which(idx == i)
      g1 <- igraph::induced_subgraph(g, ii)
      newidx[ii] <- paste(i, ":", cluster_louvain(g1)$membership)
    }
    levels <- names(sort(table(newidx), decreasing = TRUE))
    idx <- as.integer(factor(newidx, levels = levels))
    K <- cbind(K, idx)
  }
  rownames(K) <- igraph::V(g)$name
  colnames(K) <- NULL
  table(idx)
  return(K)
}

# ;n=3;k=10





#' @export
hclust_graph <- function(g, k = NULL, mc.cores = 2) {
  ## Hierarchical clustering of graph using iterative Louvain
  ## clustering on different levels. If k=NULL iterates until
  ## convergences.
  ##

  idx <- rep(1, length(igraph::V(g)))
  K <- c()
  maxiter <- 100
  if (!is.null(k)) maxiter <- k
  iter <- 1
  ok <- 1
  idx.len <- -1
  while (iter <= maxiter && ok) {
    old.len <- idx.len
    newidx0 <- newidx <- idx
    i <- idx[1]
    if (mc.cores > 1 && length(unique(idx)) > 1) {
      idx.list <- tapply(1:length(idx), idx, list)
      mc.cores
      system.time(newidx0 <- parallel::mclapply(idx.list, function(ii) {
        subg <- igraph::induced_subgraph(g, ii)
        subi <- igraph::cluster_louvain(subg)$membership
        return(subi)
      }, mc.cores = mc.cores))
      newidx0 <- lapply(1:length(newidx0), function(i) paste0(i, "-", newidx0[[i]]))
      newidx0 <- as.vector(unlist(newidx0))
      newidx <- rep(NA, length(idx))
      newidx[as.vector(unlist(idx.list))] <- newidx0
    } else {
      for (i in unique(idx)) {
        ii <- which(idx == i)
        subg <- igraph::induced_subgraph(g, ii)
        subi <- igraph::cluster_louvain(subg)$membership
        newidx[ii] <- paste(i, subi, sep = "-")
      }
    }
    vv <- names(sort(table(newidx), decreasing = TRUE))
    idx <- as.integer(factor(newidx, levels = vv))
    K <- cbind(K, idx)
    idx.len <- length(table(idx))
    ok <- (idx.len > old.len)
    iter <- iter + 1
  }
  rownames(K) <- igraph::V(g)$name
  if (!ok && is.null(k)) K <- K[, 1:(ncol(K) - 1)]
  dim(K)

  colnames(K) <- NULL
  return(K)
}









## ===================================================================================
## ============================== END OF FILE ========================================
## ===================================================================================
