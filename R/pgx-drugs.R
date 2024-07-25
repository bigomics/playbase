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

  if (is.null(X)) {
    X <- playdata::L1000_ACTIVITYS_N20D1011
    dim(X)
  }

  if ("gx.meta" %in% names(obj)) {
    FC <- pgx.getMetaMatrix(obj)$fc
    ## check if multi-omics
    is.multiomics <- any(grepl("\\[gx\\]|\\[mrna\\]", rownames(FC)))
    if (is.multiomics) {
      jj <- grep("\\[gx\\]|\\[mrna\\]", rownames(FC))
      FC <- FC[jj, , drop = FALSE]
    }
    rownames(FC) <- toupper(sub(".*:|.*\\]", "", rownames(FC)))

    FC <- FC[order(-rowMeans(FC**2, na.rm = TRUE)), , drop = FALSE]
    FC <- FC[!duplicated(rownames(FC)), , drop = FALSE]
  } else {
    ## it is a matrix
    FC <- obj
  }

  if (is.null(contrast)) {
    contrast <- colnames(FC)
  }
  contrast <- intersect(contrast, colnames(FC))

  FC <- FC[, contrast, drop = FALSE]

  if (!obj$organism %in% c("Human", "human")) {
    human_genes <- ifelse(!is.na(obj$genes$human_ortholog),
      obj$genes$human_ortholog,
      obj$genes$symbol
    )
    rownames(FC) <- human_genes
  }

  ## create drug meta sets
  meta.gmt <- tapply(colnames(X), xdrugs, list)
  meta.gmt <- meta.gmt[which(sapply(meta.gmt, length) >= nmin)]

  if (length(meta.gmt) == 0) {
    message("WARNING::: pgx.computeDrugEnrichment : no valid genesets!!")
    return(NULL)
  }

  ## first level (rank) correlation
  message("Calculating first level rank correlation ...")
  jj <- match(rownames(FC), rownames(obj$genes)) ## Enable enrichment for proteomics
  rownames(FC) <- obj$genes$symbol[jj]
  gg <- intersect(rownames(X), rownames(FC))

  if (length(gg) < 20) {
    message("WARNING::: pgx.computeDrugEnrichment : not enough common genes!!")
    return(NULL)
  }
  if (any(class(X) == "dgCMatrix")) {
    ## gene set enrichment by rank correlation
    fx <- apply(FC[gg, , drop = FALSE], 2, rank)
    R1 <- qlcMatrix::corSparse(X[gg, ], fx)
  } else {
    rnk1 <- t(matrixStats::colRanks(X[gg, , drop = FALSE], ties.method = "average", na.last = "keep"))
    rnk2 <- t(matrixStats::colRanks(FC[gg, , drop = FALSE], ties.method = "average", na.last = "keep"))
    system.time(R1 <- stats::cor(rnk1, rnk2, use = "pairwise"))
  }
  R1 <- as.matrix(R1)
  R1[is.nan(R1)] <- 0
  R1[is.infinite(R1)] <- 0
  R1 <- R1 + 1e-8 * matrix(stats::rnorm(length(R1)), nrow(R1), ncol(R1))
  colnames(R1) <- colnames(FC)
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
    rho2 <- rho2[order(-rowMeans(rho2**2, na.rm = TRUE)), , drop = FALSE]
    .cor.pvalue <- function(x, n) 2 * stats::pnorm(-abs(x / ((1 - x**2) / (n - 2))**0.5))
    P <- apply(rho2, 2, .cor.pvalue, n = nrow(D))
    Q <- apply(P, 2, stats::p.adjust, method = "fdr")
    results[["cor"]] <- list(X = rho2, Q = Q, P = P)
  }

  if ("GSEA" %in% methods) {
    bpparam <- BiocParallel::MulticoreParam(1)
    message("Calculating drug enrichment using GSEA ...")
    res0 <- list()
    i <- 1
    for (i in 1:ncol(R1)) {
      suppressWarnings(res0[[i]] <- fgsea::fgseaSimple(meta.gmt, stats = R1[, i], nperm = 10000, BPPARAM = bpparam))
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
    colnames(mNES) <- colnames(mQ) <- colnames(mP) <- colnames(FC)
    msize <- res0[[1]]$size
    results[["GSEA"]] <- list(X = mNES, Q = mQ, P = mP, size = msize)
  }

  ## this takes only the top matching drugs for each comparison to
  ## reduce the size of the matrices
  if (nprune > 0) {
    message("[pgx.computeDrugEnrichment] pruning : nprune = ", nprune)
    k <- 1
    for (k in 1:length(results)) {
      res <- results[[k]]
      ## reduce solution set with top-N of each comparison??

      rx <- apply(abs(res$X), 2, rank)
      rownames(rx) <- rownames(res$X)
      mtop <- names(utils::head(sort(rowMeans(rx, na.rm = TRUE), decreasing = TRUE), nprune))
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


#' @export
pgx.plotDrugConnectivity <- function(pgx, contrast,
                                     db = "L1000/activity",
                                     moatype = c("class", "gene", "drug")[1],
                                     ntop = 10) {
  if (!"drugs" %in% names(pgx)) {
    stop("pgx does not have drug enrichment results")
  }

  if (!db %in% names(pgx$drugs)) {
    stop("pgx$drugs does not have database db = ", db)
  }

  ## common getData-esque function for drug connectivity plots / tables
  getActiveDSEA <- function(pgx, contrast, db) {
    dr <- pgx$drugs[[db]]
    nes <- round(dr$X[, contrast], 4)
    pv <- round(dr$P[, contrast], 4)
    qv <- round(dr$Q[, contrast], 4)
    drug <- rownames(dr$X)
    stats <- dr$stats[, contrast]
    annot <- dr$annot
    nes[is.na(nes)] <- 0
    qv[is.na(qv)] <- 1
    pv[is.na(pv)] <- 1

    ## !!!SHOULD MAYBE BE DONE IN PREPROCESSING???
    if (is.null(annot)) {
      stop("[getActiveDSEA] WARNING:: missing drug annotation in PGX file!")
    }

    ## compile results matrix
    jj <- match(toupper(drug), toupper(rownames(annot)))
    annot <- annot[jj, c("moa", "target")]
    dt <- data.frame(drug = drug, NES = nes, pval = pv, padj = qv, annot)
    dt <- dt[order(-dt$NES), ]

    filter_empty <- FALSE
    if (filter_empty) {
      sel <- which(dt$moa != "" | dt$target != "")
      dt <- dt[sel, , drop = FALSE]
    }
    dsea <- list(table = dt, clust = dr$clust, stats = stats)
    return(dsea)
  }

  getMOA.target <- function(dsea) {
    ## meta-GSEA on molecular targets
    dt <- dsea$table
    targets.list <- lapply(
      enc2utf8(as.character(dt$target)),
      function(s) trimws(strsplit(s, split = "[\\|;,]")[[1]])
    )
    names(targets.list) <- rownames(dt)
    targets <- setdiff(unique(unlist(targets.list)), c(NA, "", " "))
    gmt <- lapply(targets, function(g) {
      names(which(sapply(targets.list, function(t) (g %in% t))))
    })
    names(gmt) <- targets

    rnk <- dt$NES
    names(rnk) <- rownames(dt)
    suppressWarnings(
      moa.target <- fgsea::fgsea(gmt, rnk, nperm = 20000)
    )
    moa.target <- moa.target[order(-abs(moa.target$NES)), ]
    return(moa.target)
  }

  getMOA.class <- function(dsea) {
    ## meta-GSEA on MOA terms
    dt <- dsea$table
    moa.list <- lapply(
      enc2utf8(as.character(dt$moa)),
      function(s) trimws(strsplit(s, split = "[\\|;,]")[[1]])
    )
    names(moa.list) <- rownames(dt)
    moa <- setdiff(unlist(moa.list), c("", NA, " "))
    gmt <- lapply(moa, function(g) names(which(sapply(moa.list, function(t) (g %in% t)))))
    names(gmt) <- moa
    rnk <- dt$NES
    names(rnk) <- rownames(dt)
    suppressWarnings(
      moa.class <- fgsea::fgsea(gmt, rnk)
    )
    moa.class <- moa.class[order(-abs(moa.class$NES)), ]
    return(moa.class)
  }

  plotTopEnriched <- function(res, ntop) {
    res$score <- res$NES

    qweight <- FALSE
    if (qweight) {
      res$score <- res$NES * (1 - res$padj) * (1 - 1 / res$size**1)
      yaxistitle <- "score (qNES)"
    }
    jj <- unique(c(head(order(-res$score), ntop), tail(order(-res$score), ntop)))
    moa.top <- res$score[jj]
    names(moa.top) <- res$pathway[jj]

    df <- data.frame(
      x = factor(names(moa.top), levels = names(moa.top)),
      y = as.numeric(moa.top)
    )

    barplot(
      height = df$y,
      names.arg = df$x,
      ylab = "score (NES)",
      xlab = "",
      las = 3,
      ylim = c(-1.1, 1.1) * max(abs(as.numeric(moa.top)))
    )
  }

  plotTopDrugs <- function(db, ntop = 15) {
    dr <- pgx$drugs[[db]]
    sel <- 1:nrow(dr$annot)
    dd <- rownames(dr$X)
    ##  sel <- grepl("[a-z]{4}", dd) & !is.na(dr$annot[dd, "moa"])
    sel <- grepl("[a-z]{4}", dd)
    dx <- sort(dr$X[sel, 1], decreasing = TRUE)
    dx.top <- c(head(dx, ntop), tail(dx, ntop))
    barplot(dx.top, las = 3, horiz = FALSE, ylab = "similarity (NES)")
  }

  ## this should move to pgx.computeDrugEnrichment...
  dsea <- getActiveDSEA(pgx, contrast, db)
  if (moatype == "gene") {
    res <- getMOA.target(dsea)
    plotTopEnriched(res, ntop = ntop)
  } else if (moatype == "class") {
    res <- getMOA.class(dsea)
    plotTopEnriched(res, ntop = ntop)
  } else if (moatype == "drug") {
    plotTopDrugs(db, ntop = ntop)
  }
}
