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
#' @param report.genes A logical value indicating whether to report genes in genesets. Defaults to `FALSE`.
#'
#' @export
#'
#' @return A data frame containing the results of the Fisher's exact test. The data frame
#'         includes columns such as "sign" (fold change sign), "odd.ratio" (odds ratio),
#'         "p.value" (p-value), "q.value" (adjusted p-value), and others.
#'
gset.fisher2 <- function(genes.up, genes.dn, genesets, background = NULL,
                         fdr = 0.05, mc = TRUE, sort.by = "p.value", nmin = 3, verbose = 1,
                         min.genes = 15, max.genes = 500, method = "fast.fisher",
                         check.background = TRUE, report.genes = FALSE) {

  ft.up <- gset.fisher(
    genes = genes.up, genesets = genesets, background = background,
    fdr = 1, mc = mc, sort.by = sort.by, nmin = nmin, verbose = verbose,
    min.genes = min.genes, max.genes = max.genes, method = method,
    check.background = check.background, report.genes = report.genes
  )
  ft.dn <- gset.fisher(
    genes = genes.dn, genesets = genesets, background = background,
    fdr = 1, mc = mc, sort.by = sort.by, nmin = nmin, verbose = verbose,
    min.genes = min.genes, max.genes = max.genes, method = method,
    check.background = check.background, report.genes = report.genes
  )
  ft.up <- ft.up[rownames(ft.dn), ]
  ft.sign <- c(-1, 1)[1 + 1 * (ft.up$p.value < ft.dn$p.value)]
  ft1 <- cbind(sign = ft.sign, ft.up)[which(ft.sign > 0), , drop = FALSE]
  ft2 <- cbind(sign = ft.sign, ft.dn)[which(ft.sign < 0), , drop = FALSE]
  ft.res <- rbind(ft1, ft2)
  ft.res$sign <- ft.res$sign * ft.res$odd.ratio
  ft.res <- ft.res[which(ft.res$q.value <= fdr), , drop = FALSE]
  head(ft.res)
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
#' @param sort.by The statistic used to sort the results. Defaults to "p.value".
#' @param nmin The minimum number of genes required in a gene set for the test to be performed.
#' @param min.genes The minimum number of genes in a gene set to be considered. Defaults to 15.
#' @param max.genes The maximum number of genes in a gene set to be considered. Defaults to 500.
#' @param method The method used for computing p-values. Defaults to "fast.fisher".
#' @param check.background A logical value indicating whether to check the presence of genes in
#'                         the background set. Defaults to `TRUE`.
#' @param report.genes A logical value indicating whether to report genes in gene sets. Defaults to `TRUE`.
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
                        fdr = 0.25, mc = TRUE, sort.by = "p.value", nmin = 3,
                        min.genes = 15, max.genes = 500, method = "fast.fisher",
                        check.background = TRUE, report.genes = FALSE,
                        no.pass=NA, verbose = 1) {

  ## switch according to geneset class
  if(class(genesets) == "list") {
    gsnames <- names(genesets)
    res <- gset.fisherLIST(
      genes = genes, genesets = genesets, background = background,
      fdr = fdr, mc = mc, sort.by = sort.by, nmin = nmin,
      min.genes = min.genes, max.genes = max.genes, method = method,
      check.background = check.background, report.genes = report.genes,
      no.pass = no.pass, verbose = verbose
    )
  } else if(inherits(genesets, "Matrix")) {
    gsnames <- colnames(genesets)
    res <- gset.fastFET(genes, G = genesets, bg = background,
      report.genes = report.genes)
  } else {
    stop("[gset.fisher] FATAL ERROR")
  }
  
  ## filter results
  if (nrow(res) > 0) {
    size.ok <- res$size >= min.genes & res$size <= max.genes
    jj <- which(res$q.value <= fdr & res$N >= nmin & size.ok)
    res <- res[jj, ]
  }

  ## sort results
  if (nrow(res) > 0) {
    if (sort.by %in% colnames(res)) {
      order.sign <- ifelse(sort.by %in% c("p.value","q.value"),+1,-1)
      res <- res[order(order.sign*res[,sort.by]), ]
    } else {
      gsnames <- intersect(gsnames, rownames(res))
      res <- res[gsnames,]
    }
  }
  
  return(res)
}


gset.fisherLIST <- function(genes, genesets, background = NULL,
                        fdr = 0.05, mc = TRUE, sort.by = "p.value", nmin = 3,
                        min.genes = 15, max.genes = 500, method = "fast.fisher",
                        check.background = TRUE, report.genes = FALSE,
                        no.pass=NA, verbose = 1) {

  bgNULL <- FALSE
  if (is.null(background)) {
    message("[gset.fisher] note: it is recommended to specify background")
    background <- unique(unlist(genesets))
    if (verbose > 0) {
      cat("setting background to ", length(background), "genes covered in genesets\n")
    }
    bgNULL <- TRUE
  }

  if (check.background && !bgNULL) {
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
        p.value = NA, q.value = NA, odd.ratio = NA, N = NA, size = NA,
        n.overlap = NA, genes = NA
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

  ## report intersection genes (slow..) like leading edge list.
  gsgenes <- NULL
  if (report.genes) {
    gsgenes <- unlist(lapply(genesets, function(x) paste(sort(intersect(genes, x)), collapse = "|")))
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
  }

  if (method == "fisher") {
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
  } else if (method == "fast.fisher") {
    ## done
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
  rr <- data.frame(p.value = pv, q.value = qv, odd.ratio = odd.ratio,
    N=a, size=n.size, overlap=v1 )

  if (!is.null(gsgenes)) {
    rr <- cbind(rr, genes = gsgenes)
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


#' Calculate fast Fisher exact test.  
#'
#' @param genes Vector of significant genes
#' @param G    Sparse matrix containing gene sets
#' @param bg   Vector of genes as background
#' @param report.genes  Logical to report gene set genes in output
#'
#' @export
gset.fastFET <- function(genes, G, bg, report.genes=FALSE) {
  bgnull <- FALSE
  if(is.null(bg)) {
    message("[gset.fisher] note: it is recommended to specify background")
    bg <- unique(c(genes,rownames(G)))
    bgnull <- TRUE
  }
  if(length(bg)>1 && !bgnull) {
    genes <- intersect(genes, bg)
    G <- G[intersect(bg,rownames(G)),,drop=FALSE]
  }  
  length.bg <- NULL
  if(length(bg)==1 && is.integer(bg)) {
    length.bg <- as.integer(bg)
  } else if(length(bg)>1) {
    length.bg <- length(bg)
  }

  if(is.null(length.bg)) stop("error: invalid background. bg:", head(bg))
  if(length(genes)==0) stop("error: zero genes length")
  if(nrow(G)==0) stop("error: empty gene set matrix G")
  
  genes <- intersect(genes, rownames(G))  
  gsize <- Matrix::colSums(G!=0)
  genes <- intersect(genes, rownames(G))
  a <- Matrix::colSums(G[genes,]!=0)
  b <- length(genes) - a
  c <- gsize - a
  d <- length.bg - (a+b+c) 
  pv <- corpora.fastFET(a,b,c,d)

  names(pv) <- colnames(G)
  odd.ratio <- (a/b)/(c/d)
  qv <- p.adjust(pv, method="fdr")
  overlap <- paste0(a,"/",gsize)

  gsgenes <- NULL
  if(report.genes) {
    gsgenes <- apply( G[genes,], 2, function(x) paste(sort(genes[which(x!=0)]),collapse="|") )
  }

  df <- data.frame(p.value=pv, q.value=qv, odd.ratio=odd.ratio,
    N=a, size=gsize, overlap=overlap )
  
  if(!is.null(gsgenes)) df <- cbind(df, genes=gsgenes)
  return(df)
}



#' Wrapper superfast version of Fisher Exact Test from 'corpora' R
#' package. This is the fastest implementation currently
#' available. Uses phyper inside.
#' 
#'            setAn ¬setA
#'        setB  a     b | a+b
#'       ¬setB  c     d | c+d
#'          ------------|-----
#'             a+c   b+d| a+b+c+d
#' 
corpora.fastFET <- function(a, b, c, d, alternative = c("two.sided", "less", 
    "greater")[3], log.p = FALSE) {
  ## this is really-really-really fast...
  pv <- rep(NA, length(a))
  ii <- 1:length(a)
  ii <- which((a + b) > 0)
  d1 <- d + 1 * (d == 0) ## hack to avoid crash...
  c1 <- c + 1 * (c == 0) ## hack to avoid crash...

  k1 <- a[ii]
  n1 <- (a + c1)[ii]
  k2 <- b[ii]
  n2 <- (b + d1)[ii]
  
  .match.len <-  function (vars, len = NULL, adjust = FALSE, check.numeric = TRUE, 
                           envir = parent.frame()) {
    vecs <- setNames(lapply(vars, get, envir = envir), vars)
    ok <- sapply(vecs, is.numeric)
    if (check.numeric && any(!ok)) 
      stop("argument(s) ", paste(vars[!ok], collapse = ", "), 
        " must be numeric vector(s)")
    if (is.null(len)) 
      len <- max(sapply(vecs, length))
    for (v in vars) {
      if (length(vecs[[v]]) == 1) {
        if (adjust) 
          assign(v, rep(vecs[[v]], len), envir = envir)
      }
      else if (length(vecs[[v]]) != len) {
        stop(sprintf("argument %s should be of length %d or a scalar (%s must have same length)", 
          v, len, paste(vars, collapse = ", ")))
      }
    }
    invisible(len)
  }
  
  alternative <- match.arg(alternative)
  l <- .match.len(c("k1", "n1", "k2", "n2"), adjust = TRUE)
  if (any(k1 < 0) || any(k1 > n1) || any(n1 <= 0)) 
    stop("k1 and n1 must be integers with 0 <= k1 <= n1")
  if (any(k2 < 0) || any(k2 > n2) || any(n2 <= 0)) 
    stop("k2 and n2 must be integers with 0 <= k2 <= n2")
  if (any(k1 + k2 <= 0)) 
    stop("either k1 or k2 must be non-zero")
  k <- k1 + k2
  if (alternative == "two.sided") {
    if (log.p) {
      pval <- pmin(phyper(k1 - 1, n1, n2, k, lower.tail = FALSE, 
        log.p = TRUE), phyper(k1, n1, n2, k, lower.tail = TRUE, 
          log.p = TRUE)) + log(2)
      pval <- pmin(pval, 0)
    }
    else {
      pval <- 2 * pmin(phyper(k1 - 1, n1, n2, k, lower.tail = FALSE), 
        phyper(k1, n1, n2, k, lower.tail = TRUE))
      pval <- pmax(0, pmin(1, pval))
    }
  }
  else if (alternative == "greater") {
    pval <- phyper(k1 - 1, n1, n2, k, lower.tail = FALSE, 
      log.p = log.p)
  }
  else if (alternative == "less") {
    pval <- phyper(k1, n1, n2, k, lower.tail = TRUE, log.p = log.p)
  }
  pval
}


