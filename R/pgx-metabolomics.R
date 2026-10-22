##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##


#' Extend metabolite sets in M by adding neighbour-metabolites that
#' are connected (neighbours) according to given PPI. Only
#' metabolite-metabolite neighbours are considered.
#'
#' @export
extend_metabolite_sets <- function(M, ppi, add = TRUE, postfix = "(extended)",
                                   maxcost = 0.33) {
  ## get metabolite-metabolite edges from GRAPHITE
  ## ppi <- .pb_libx("graphite_ppi")
  ppi[, 1] <- ifelse(grepl("CHEBI", ppi[, 1]), ppi[, 1], paste0("SYMBOL:", ppi[, 1]))
  ppi[, 2] <- ifelse(grepl("CHEBI", ppi[, 2]), ppi[, 2], paste0("SYMBOL:", ppi[, 2]))
  sel <- which(grepl("CHEBI", ppi[, 1]) & grepl("CHEBI", ppi[, 2]) & ppi[, 3] <= maxcost)
  gr <- igraph::graph_from_edgelist(as.matrix(ppi[sel, 1:2]), directed = FALSE)
  MMI <- as.matrix(gr)
  table(colnames(M) %in% rownames(MMI))
  table(rownames(MMI) %in% colnames(M))

  ## build extended sparse matrix
  idx <- Matrix::which(M != 0, arr.ind = TRUE)
  x <- M[idx]
  features <- unique(c(colnames(M), rownames(MMI)))
  dims <- c(nrow(M), length(features))
  dimnames <- list(rownames(M), features)
  M2 <- Matrix::sparseMatrix(idx[, 1], idx[, 2],
    x = x,
    dims = dims, dimnames = dimnames
  )

  ## MMI neigborhood matrix for metabolites. match with extended matrix M2
  mmi0 <- rbind("na" = 0, cbind("na" = 0, MMI))
  ii <- match(colnames(M2), rownames(mmi0))
  ii[is.na(ii)] <- 1
  B <- mmi0[ii, ii]
  diag(B) <- 1
  colnames(B) <- rownames(B) <- colnames(M2)

  ## propagate neighbor
  extM <- M2 %*% B
  extM <- 1 * (extM != 0)
  extM

  ## row merge extended gene sets with original
  rownames(extM) <- paste(rownames(extM), postfix)
  if (add == TRUE) {
    extM <- Matrix::t(merge_sparse_matrix(Matrix::t(M), Matrix::t(extM)))
  }
  extM <- extM[order(rownames(extM)), ]
  return(extM)
}


#' Extend metabolite sets in M by adding neighbouring-metabolites that
#' are connected (neighbours) with members in metabolite set according
#' to given PPI. Any metabolite-metabolite or protein-metabolite
#' neighbours are considered.
#'
#' @export
extend_metabolite_sets2 <- function(M, ppi, add = TRUE, postfix = "(extended)", maxcost = 0.33) {
  ## get metabolite-metabolite edges from GRAPHITE
  ppi[, 1] <- ifelse(grepl("CHEBI", ppi[, 1]), ppi[, 1], paste0("SYMBOL:", ppi[, 1]))
  ppi[, 2] <- ifelse(grepl("CHEBI", ppi[, 2]), ppi[, 2], paste0("SYMBOL:", ppi[, 2]))
  sel <- which((grepl("CHEBI", ppi[, 1]) | grepl("CHEBI", ppi[, 2])) & ppi[, 3] <= maxcost)
  gr <- igraph::graph_from_edgelist(as.matrix(ppi[sel, 1:2]), directed = FALSE)
  PMI <- as.matrix(gr)

  ## build extended sparse matrix
  idx <- Matrix::which(M != 0, arr.ind = TRUE)
  x <- M[idx]
  features <- unique(c(colnames(M), rownames(PMI)))
  dims <- c(nrow(M), length(features))
  dimnames <- list(rownames(M), features)
  M2 <- Matrix::sparseMatrix(idx[, 1], idx[, 2],
    x = x,
    dims = dims, dimnames = dimnames
  )

  ## PMI neigborhood matrix. match with extended matrix M2
  mmi0 <- rbind("na" = 0, cbind("na" = 0, PMI))
  ii <- match(colnames(M2), rownames(mmi0))
  ii[is.na(ii)] <- 1 ## no matches will map to row/col 1 (all zero)
  B <- mmi0[ii, ii] ## B is now aligned to M2
  diag(B) <- 1
  colnames(B) <- rownames(B) <- colnames(M2)

  ## propagate all nodes to neighbors
  extM <- M2 %*% B
  extM <- 1 * (extM != 0)
  dim(extM)

  ## do not allow gene-gene neigbors. Add extended metabolites to
  ## original matrix M2
  ii <- grepl("SYMBOL", colnames(extM))
  extM[, ii] <- 0
  extM <- 1 * ((M2 + extM) != 0)

  ## row merge extended gene sets with original
  rownames(extM) <- paste(rownames(extM), postfix)
  if (add == TRUE) {
    extM <- Matrix::t(merge_sparse_matrix(Matrix::t(M2), Matrix::t(extM)))
  }
  extM <- extM[order(rownames(extM)), ]

  return(extM)
}
