##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

## ========================================================================
## ======================== Fisher test based =============================
## ========================================================================

#' Perform Fisher's exact test on gene sets
#'
#' This function performs Fisher's exact test on two sets of genes, `genes.up` and
#' `genes.dn`, within a given set of gene sets (`genesets`). It returns a data frame
#' containing the results of the test, including the sign of the fold change (positive
#' or negative) and relevant statistics such as p-values, q-values, and odds ratios.
#'
#' @param genes.up A character vector containing the names of genes in the "up" set.
#' @param genes.dn A character vector containing the names of genes in the "down" set.
#' @param genesets A list of gene sets, where each element is a character vector
#'                 representing a gene set.
#' @param background A character vector containing the names of genes in the background set.
#'                   Defaults to `NULL`, which means all genes are considered.
#' @param fdr The false discovery rate (FDR) threshold for multiple testing adjustment.
#'            Defaults to 0.05.
#' @param mc A logical value indicating whether to perform multiple testing adjustment
#'           using Monte Carlo simulation. Defaults to `TRUE`.
#' @param sort.by The statistic used to sort the results. Defaults to "zratio".
#' @param nmin The minimum number of genes required in a gene set for the test to be performed.
#' @param verbose A numeric value indicating the level of verbosity. Defaults to 1.
#' @param min.genes The minimum number of genes in a gene set to be considered. Defaults to 15.
#' @param max.genes The maximum number of genes in a gene set to be considered. Defaults to 500.
#' @param method The method used for computing p-values. Defaults to "fast.fisher".
#' @param check.background A logical value indicating whether to check the presence of genes in
#'                         the background set. Defaults to `TRUE`.
#' @param common.genes A logical value indicating whether to use only genes common to both the
#'                     input gene sets and the background set. Defaults to `TRUE`.
#'
#' @export
#'
#' @return A data frame containing the results of the Fisher's exact test. The data frame
#'         includes columns such as "sign" (fold change sign), "odd.ratio" (odds ratio),
#'         "p.value" (p-value), "q.value" (adjusted p-value), and others.
#'
gset.fisher2 <- function(genes.up, genes.dn, genesets, background = NULL,
                         fdr = 0.05, mc = TRUE, sort.by = "zratio", nmin = 3, verbose = 1,
                         min.genes = 15, max.genes = 500, method = "fast.fisher",
                         check.background = TRUE, common.genes = TRUE) {
  ft.up <- gset.fisher(
    genes = genes.up, genesets = genesets, background = background,
    fdr = 1, mc = mc, sort.by = sort.by, nmin = nmin, verbose = verbose,
    min.genes = min.genes, max.genes = max.genes, method = method,
    check.background = check.background, common.genes = common.genes
  )
  ft.dn <- gset.fisher(
    genes = genes.dn, genesets = genesets, background = background,
    fdr = 1, mc = mc, sort.by = sort.by, nmin = nmin, verbose = verbose,
    min.genes = min.genes, max.genes = max.genes, method = method,
    check.background = check.background, common.genes = common.genes
  )
  ft.up <- ft.up[rownames(ft.dn), ]
  ft.sign <- c(-1, 1)[1 + 1 * (ft.up$p.value < ft.dn$p.value)]
  ft1 <- cbind(sign = ft.sign, ft.up)[which(ft.sign > 0), , drop = FALSE]
  ft2 <- cbind(sign = ft.sign, ft.dn)[which(ft.sign < 0), , drop = FALSE]
  ft.res <- rbind(ft1, ft2)
  ft.res$sign <- ft.res$sign * ft.res$odd.ratio
  ft.res <- ft.res[which(ft.res$q.value <= fdr), , drop = FALSE]
  return(ft.res)
}


#' Perform Fisher's exact test on gene sets
#'
#' This function performs Fisher's exact test on a set of genes within a given set of gene sets.
#' It returns a data frame containing the results of the test, including p-values, q-values,
#' odds ratios, and gene set overlaps.
#'
#' @param genes A character vector containing the names of genes.
#' @param genesets A list of gene sets, where each element is a character vector
#'                 representing a gene set.
#' @param background A character vector containing the names of genes in the background set.
#'                   Defaults to `NULL`, which means all genes are considered.
#' @param fdr The false discovery rate (FDR) threshold for multiple testing adjustment.
#'            Defaults to 0.05.
#' @param mc A logical value indicating whether to perform multiple testing adjustment
#'           using Monte Carlo simulation. Defaults to `TRUE`.
#' @param sort.by The statistic used to sort the results. Defaults to "zratio".
#' @param nmin The minimum number of genes required in a gene set for the test to be performed.
#' @param min.genes The minimum number of genes in a gene set to be considered. Defaults to 15.
#' @param max.genes The maximum number of genes in a gene set to be considered. Defaults to 500.
#' @param method The method used for computing p-values. Defaults to "fast.fisher".
#' @param check.background A logical value indicating whether to check the presence of genes in
#'                         the background set. Defaults to `TRUE`.
#' @param common.genes A logical value indicating whether to use only genes common to both the
#'                     input gene set and the background set. Defaults to `TRUE`.
#' @param verbose A numeric value indicating the level of verbosity. Defaults to 1.
#'
#' @export
#'
#' @return A data frame containing the results of the Fisher's exact test. The data frame
#'         includes columns such as "p.value" (p-value), "q.value" (adjusted p-value),
#'         "odd.ratio" (odds ratio), "overlap" (gene set overlap), and optionally "genes"
#'         (common genes).
#'
gset.fisher <- function(genes, genesets, background = NULL,
                        fdr = 0.05, mc = TRUE, sort.by = "zratio", nmin = 3,
                        min.genes = 15, max.genes = 500, method = "fast.fisher",
                        check.background = TRUE, common.genes = TRUE,
                        no.pass=NA, verbose = 1) {
  if (is.null(background)) {
    background <- unique(unlist(genesets))
    if (verbose > 0) {
      cat("setting background to ", length(background), "genes covered\n")
    }
  }

  if (check.background) {
    ## restrict on background
    genes <- intersect(genes, background)
    genesets <- lapply(genesets, function(s) intersect(s, background))
  }

  ## select
  if (!is.null(min.genes) && min.genes > 0) {
    genesets.len <- sapply(genesets, length)
    genesets <- genesets[order(-genesets.len)]
    if (sum(duplicated(names(genesets))) > 0) {
      if (verbose > 0) cat("warning: duplicated gene set names. taking largest.\n")
      genesets <- genesets[which(!duplicated(names(genesets)))]
    }
    genesets.len <- sapply(genesets, length)
    genesets <- genesets[genesets.len >= min.genes & genesets.len <= max.genes]
    if (verbose > 0) {
      cat(
        "testing", length(genesets), "genesets with", length(genes),
        "genes (background", length(background), "genes)\n"
      )
    }
    length(genesets)
    if (length(genesets) == 0) {
      cat("warning: no gene sets passed size filter\n")
      rr <- data.frame(
        p.value = NA, q.value = NA, ratio0 = NA, ratio1 = NA,
        zratio = NA, n.size = NA, n.overlap = NA, genes = NA
      )
      rownames(rr) <- NULL
      return(rr[0, ])
    }
  }

  ## odd ratio
  ## see http://jura.wi.mit.edu/bio/education/hot_topics/enrichment/Gene_list_enrichment_Mar10.pdf
  n.size <- sapply(genesets, length)
  bg0 <- setdiff(background, genes)
  nbackground0 <- length(background)
  nbackground1 <- length(bg0)

  ## this can become slow
  a <- unlist(lapply(genesets, function(x) sum(x %in% genes)))
  b <- (n.size - a)
  c <- unlist(lapply(genesets, function(x) sum(!(genes %in% x))))
  d <- (nbackground1 - b)
  odd.ratio <- (a / c) / (b / d) ## note: not exactly same as from fishertest

  ## intersection genes (slow..)
  commongenes <- NULL
  if (common.genes) {
    commongenes <- unlist(lapply(genesets, function(x) paste(sort(intersect(genes, x)), collapse = "|")))
  }

  ## compute fisher-test (should be one-sided?)
  test.fisher <- function(gs) {
    a0 <- table(background %in% gs, background %in% genes)
    if (NCOL(a0) == 1 || colSums(a0)[2] == 0) {
      return(NA)
    }
    stats::fisher.test(a0, alternative = "greater")$p.value
  }
  test.chisq <- function(gs) {
    a0 <- table(background %in% gs, background %in% genes)
    if (NCOL(a0) == 1 || colSums(a0)[2] == 0) {
      return(NA)
    }
    stats::chisq.test(a0)$p.value
  }
  pv <- rep(NA, length(genesets))
  names(pv) <- names(genesets)
  if (method == "fast.fisher") {
    ## this is really-really-really fast...
    pv <- rep(NA, length(a))
    ii <- 1:length(a)
    ii <- which((a + c) > 0)
    d1 <- d + 1 * (d == 0) ## hack to avoid crash...
    b1 <- b + 1 * (b == 0) ## hack to avoid crash...
    pv1 <- try(
      corpora::fisher.pval(a[ii], (a + b1)[ii], c[ii], (c + d1)[ii], alternative = "greater"),
      silent = TRUE
    )
    if (class(pv1) != "try-error") {
      pv[ii] <- pv1
    } else {
      message("[playbase::gset.fisher] fast.fisher failed. Testing with standard fisher.")
      method <- "fisher"
    }
  } else if (method == "fisher") {
    if (mc) {
      pv <- unlist(lapply(genesets, test.fisher))
    } else {
      i <- 1
      for (i in 1:length(genesets)) {
        pv[i] <- test.fisher(genesets[[i]])
      }
    }
  } else if (method == "chisq") {
    if (mc) {
      pv <- unlist(lapply(genesets, test.chisq))
    } else {
      for (i in 1:length(genesets)) {
        pv[i] <- test.chisq(genesets[[i]])
      }
    }
  } else {
    stop("unknown method")
  }

  ## replace NA values
  if(any(is.na(pv))) {
    pv[is.na(pv)] <- no.pass
  }
  
  ## compute q-value
  qv <- rep(NA, length(pv))
  qv <- stats::p.adjust(pv, method = "fdr")

  ## results
  v1 <- as.character(paste0(a, "/", n.size))
  rr <- data.frame(p.value = pv, q.value = qv, odd.ratio = odd.ratio, overlap = v1)
  if (!is.null(commongenes)) {
    rr <- cbind(rr, genes = commongenes)
  }
  rownames(rr) <- names(genesets)

  ## sort
  if (nrow(rr) > 0) {
    ## filter
    jj <- which(rr$q.value <= fdr & n.size >= nmin)
    rr <- rr[jj, ]
    ## sort
    if (sort.by %in% c("pvalue", "p.value", "p")) {
      rr <- rr[order(rr$p.value), ]
    } else {
      rr <- rr[order(rr$odd.ratio, decreasing = TRUE), ]
    }
  }
  dim(rr)
  rr
}
