##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

## ===================================================================================
## ============================== GRAPH METHODS ======================================
## ===================================================================================


#' @title Calculate edge values from X
#' Calculate edge similarity
#'
#' @title Calculate edge similarity
#' 
#' @param ee Edge matrix 
#' @param X Expression matrix
#' @param nk Number of edges per chunk. Default 4000. 
#' @param mc.cores Number of cores for parallel processing. Default 1.
#'
#' @return Matrix of edge similarity values
#' 
#' @description Calculates similarity values for a graph edge matrix based on gene expression data.
#'
#' @details This function takes an edge matrix \code{ee} and expression matrix \code{X}, and calculates similarity values for each edge.
#' It splits the edges into \code{nk} sized chunks and processes them in parallel using \code{mc.cores} cores.
#' 
#' For each edge, it extracts the linked genes, gets their expression profiles from \code{X}, and calculates correlation or distance.
#' The output is a matrix with updated edge weights based on gene expression similarity.
#'
#' Parallel processing improves performance for large graphs.
#' 
#' @examples 
#' \dontrun{
#' ee <- matrix(c(1,2, 2,3), ncol=2) # example edge matrix
#' X <- matrix(rnorm(100*50), 100, 50) # expression matrix 
#' w <- calc.edge.similarity(ee, X)
#' }
#' @export
calc.edge.similarity <- function(ee, X, nk = 4000, mc.cores = 1) {
  if (mc.cores > 1) {
    return(calc.edge.similarityMC(ee, X, nk = nk, mc.cores = mc.cores))
  } else {
    return(calc.edge.similarityKFOLD(ee, X, nk = nk))
  }
}


#' @describeIn calc.edge.similarity Calculates the similarity between edges in a graph in parallel using multiple cores
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
#' @describeIn calc.edge.similarity Calculates the similarity between edges in a graph in parallel using Kfold algorithm
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


#' Cut graph crossings
#'
#' @title Cut graph crossings 
#'
#' @param g An igraph graph object
#' @param idx Vector of node membership indices
#' @param max.wt Maximum edge weight to cut. Default 9999.
#'
#' @return An igraph graph with cross-cluster edges removed  
#'
#' @description Cuts inter-cluster edges in a graph based on node membership
#'
#' @details This function takes an igraph graph object \code{g} and a vector of node membership indices \code{idx}. 
#' It identifies edges that connect nodes of different clusters based on \code{idx}.
#' Any edges with weight less than \code{max.wt} that link across clusters are removed from the graph.
#' The resulting pruned graph is returned.
#' 
#' By default only edges with very high weight (9999) are kept. This effectively cuts crossings between clusters.
#' The \code{max.wt} parameter can be lowered to allow more inter-cluster edges.
#'
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


#' Louvain Clustering for Graphs
#'
#' @title Louvain Clustering for Graphs
#'
#' @param g An igraph graph object to cluster
#' @param n The number of Louvain iterations to perform. Default is 3.
#' 
#' @return The updated graph object with a "louvain" vertex attribute containing cluster assignments.
#'
#' @description Performs Louvain clustering on an igraph graph object to detect communities.
#'
#' @details This function implements the Louvain algorithm for community detection on a graph.
#' It takes an igraph graph object \code{g} and performs \code{n} iterations of Louvain clustering.
#' In each iteration it groups nodes into communities that maximize modularity.
#' 
#' The algorithm optimizes modularity in a greedy fashion by first assigning each node to its own community. 
#' It then goes through nodes repeatedly to evaluate moving them to neighboring communities. If a move increases modularity it is accepted.
#' This local optimization is applied iteratively to hierarchically build communities.
#'
#' The number of iterations \code{n} controls the granularity of the detected communities.
#' More iterations lead to more fine-grained communities.
#' 
#' The graph object \code{g} is updated in-place by adding a "louvain" vertex attribute containing the cluster assignment of each node after \code{n} iterations.
#' The updated \code{g} is returned by the function.
#'
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


#' Hierarchical Clustering of Graph
#'
#' @title Hierarchical Clustering of Graph 
#'
#' @param g An igraph graph object to cluster
#' @param k The number of hierarchical levels. If NULL iterates until convergence.
#' @param mc.cores Number of cores for parallel processing. Default 2.
#'
#' @return A matrix with hierarchical clustering assignments for nodes.
#' 
#' @description Performs hierarchical clustering on a graph using iterative Louvain clustering.
#'
#' @details This function takes an igraph graph object \code{g} and performs hierarchical clustering to detect communities.
#' It uses iterative Louvain clustering, optimizing modularity at each level of the hierarchy.
#' 
#' At each iteration, it runs Louvain clustering on the communities from the previous level.
#' This splits up the communities into smaller sub-communities in a hierarchical fashion.
#' 
#' The number of levels \code{k} can be specified, otherwise it iterates until convergence.
#' Parallel processing with \code{mc.cores} is used to speed up the computations.
#' 
#' The algorithm returns a matrix containing the hierarchical clustering assignments for each node.
#' The columns represent the clustering at each level of the hierarchy.
#'
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

  colnames(K) <- NULL
  return(K)
}

## ===================================================================================
## ============================== END OF FILE ========================================
## ===================================================================================