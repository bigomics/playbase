##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


#' Create annotation for drug combinations
#'
#' @title Create annotation for drug combinations
#'
#' @param combo Character vector of drug combinations (separated by '+')
#' @param annot0 Data frame with annotation for individual drugs
#'
#' @return Data frame with annotation for the drug combinations
#'
#' @description Creates annotation for drug combinations by combining annotation
#'   from individual drugs.
#'
#' @details This function takes a vector of drug combinations and an annotation
#'   data frame for individual drugs. It splits the combinations into individual
#'   drugs, looks up their annotation, and combines it into annotation for the combos.
#'
#'   For example, if combo is 'DrugA + DrugB', it will look up mode of action (moa)
#'   and targets for DrugA and DrugB, and combine them into moa and targets for
#'   the combination. NA values are removed.
#'
#' @export
pgx.createComboDrugAnnot <- function(combo, annot0) {
  drugs <- strsplit(combo, split = "[+]")
  drug1 <- sapply(drugs, "[", 1)
  drug2 <- sapply(drugs, "[", 2)
  j1 <- match(drug1, rownames(annot0))
  j2 <- match(drug2, rownames(annot0))
  cmoa <- paste(annot0[j1, "moa"], "+", annot0[j2, "moa"])
  ctarget <- paste(annot0[j1, "target"], "+", annot0[j2, "target"])
  cmoa <- gsub("NA [+] NA", "", cmoa)
  ctarget <- gsub("NA [+] NA", "", ctarget)
  annot1 <- data.frame(drug = combo, moa = cmoa, target = ctarget)
  rownames(annot1) <- combo
  return(annot1)
}


#' @title Compute drug enrichment from expression data
#'
#' @param X Numeric gene expression matrix
#' @param drugs Character vector of drug names
#' @param nmin Minimum number of genes per drug set
#' @param methods Methods for enrichment (rank, gsea)
#'
#' @return List with enrichment results for each drug
#'
#' @description Computes enrichment of drug gene sets in expression data using
#' rank correlation or GSEA.
#'
#' @details This function takes a gene expression matrix and a set of drugs. It extracts
#' gene sets for each drug from MSigDB. Enrichment is computed by:
#'
#' 1) Rank correlation of drug set ranks with experiment ranks.
#' 2) GSEA using the drug sets as gene sets.
#'
#' Drugs with fewer than nmin genes are filtered out. The output is a list
#' containing the enrichment results for each drug.
#'
#' @export
pgx.computeDrugEnrichment <- function(obj, X, xdrugs, drug_info = NULL,
                                      methods = c("GSEA", "cor"),
                                      nmin = 15, nprune = 250, contrast = NULL) {
  ## 'obj'   : can be ngs object or fold-change matrix
  ## X       : drugs profiles (may have multiple for one drug)
  ## xdrugs  : drug associated with profile

  if(is.null(X)) {
    X <- playdata::L1000_ACTIVITYS_N20D1011
    dim(X)
  }
  
  if ("gx.meta" %in% names(obj)) {
    F <- pgx.getMetaMatrix(obj)$fc
    ## check if multi-omics
    is.multiomics <- any(grepl("\\[gx\\]|\\[mrna\\]", rownames(F)))
    is.multiomics
    if (is.multiomics) {
      jj <- grep("\\[gx\\]|\\[mrna\\]", rownames(F))
      F <- F[jj, , drop = FALSE]
    }
    rownames(F) <- toupper(sub(".*:|.*\\]", "", rownames(F)))

    F <- F[order(-rowMeans(F**2)), , drop = FALSE]
    F <- F[!duplicated(rownames(F)), , drop = FALSE]
  } else {
    ## it is a matrix
    F <- obj
  }

  if (is.null(contrast)) {
    contrast <- colnames(F)
  }
  contrast <- intersect(contrast, colnames(F))
  contrast
  F <- F[, contrast, drop = FALSE]

  ## create drug meta sets
  meta.gmt <- tapply(colnames(X), xdrugs, list)
  meta.gmt <- meta.gmt[which(sapply(meta.gmt, length) >= nmin)]

  if (length(meta.gmt) == 0) {
    message("WARNING::: pgx.computeDrugEnrichment : no valid genesets!!")
    return(NULL)
  }

  ## first level (rank) correlation
  message("Calculating first level rank correlation ...")
  gg <- intersect(rownames(X), rownames(F))

  if (length(gg) < 20) {
    message("WARNING::: pgx.computeDrugEnrichment : not enough common genes!!")
    return(NULL)
  }
  if (any(class(X) == "dgCMatrix")) {
    ## gene set enrichment by rank correlation
    fx <- apply(F[gg, , drop = FALSE], 2, rank)
    R1 <- qlcMatrix::corSparse(X[gg, ], fx)
  } else {
    rnk1 <- apply(X[gg, , drop = FALSE], 2, rank, na.last = "keep")
    rnk2 <- apply(F[gg, , drop = FALSE], 2, rank, na.last = "keep")
    system.time(R1 <- stats::cor(rnk1, rnk2, use = "pairwise"))
  }
  R1 <- as.matrix(R1)
  R1[is.nan(R1)] <- 0
  R1[is.infinite(R1)] <- 0
  R1 <- R1 + 1e-8 * matrix(stats::rnorm(length(R1)), nrow(R1), ncol(R1))
  colnames(R1) <- colnames(F)
  rownames(R1) <- colnames(X)

  ## experiment to drug
  results <- list()
  if ("cor" %in% methods) {
    message("Calculating drug enrichment using rank correlation ...")
    D <- Matrix::sparse.model.matrix(~ 0 + xdrugs)
    colnames(D) <- sub("^xdrugs", "", colnames(D))
    rownames(D) <- colnames(X) ## not necessary..
    rho2 <- qlcMatrix::corSparse(D, R1)
    rownames(rho2) <- colnames(D)
    colnames(rho2) <- colnames(R1)
    rho2 <- rho2[order(-rowMeans(rho2**2)), , drop = FALSE]
    .cor.pvalue <- function(x, n) 2 * stats::pnorm(-abs(x / ((1 - x**2) / (n - 2))**0.5))
    P <- apply(rho2, 2, .cor.pvalue, n = nrow(D))
    Q <- apply(P, 2, stats::p.adjust, method = "fdr")
    results[["cor"]] <- list(X = rho2, Q = Q, P = P)
  }

  if ("GSEA" %in% methods) {
    message("Calculating drug enrichment using GSEA ...")
    res0 <- list()
    i <- 1
    for (i in 1:ncol(R1)) {
      suppressWarnings(res0[[i]] <- fgsea::fgseaSimple(meta.gmt, stats = R1[, i], nperm = 10000))
    }
    names(res0) <- colnames(R1)

    mNES <- sapply(res0, function(x) x$NES)
    mQ <- sapply(res0, function(x) x$padj)
    mP <- sapply(res0, function(x) x$pval)
    if (length(res0) == 1) {
      mNES <- cbind(mNES)
      mP <- cbind(mP)
      mQ <- cbind(mQ)
    }

    pw <- res0[[1]]$pathway
    rownames(mNES) <- rownames(mQ) <- rownames(mP) <- pw
    colnames(mNES) <- colnames(mQ) <- colnames(mP) <- colnames(F)
    msize <- res0[[1]]$size
    results[["GSEA"]] <- list(X = mNES, Q = mQ, P = mP, size = msize)
  }

  ## level2 and level3
  if(!is.null(drug_info)) {





  }

  
  ## this takes only the top matching drugs for each comparison to
  ## reduce the size of the matrices
  nprune
  if (nprune > 0) {
    message("[pgx.computeDrugEnrichment] pruning : nprune = ", nprune)
    k <- 1
    for (k in 1:length(results)) {
      res <- results[[k]]
      ## reduce solution set with top-N of each comparison??

      rx <- apply(abs(res$X), 2, rank)
      rownames(rx) <- rownames(res$X)
      mtop <- names(utils::head(sort(rowMeans(rx), decreasing = TRUE), nprune))
      top.idx <- unique(unlist(meta.gmt[mtop]))

      results[[k]]$X <- res$X[mtop, , drop = FALSE]
      results[[k]]$P <- res$P[mtop, , drop = FALSE]
      results[[k]]$Q <- res$Q[mtop, , drop = FALSE]
      results[[k]]$size <- res$size[mtop]
    }
  }

  ## reduce large stats object
  sel.drugs <- unique(unlist(sapply(results, function(res) rownames(res$X))))
  sel <- which(xdrugs %in% sel.drugs)
  results$stats <- R1[sel, , drop = FALSE]

  message("[pgx.computeDrugEnrichment] done!")

  return(results)
}
