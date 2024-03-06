##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


method <- "edgeR"
method <- "ttest"
method <- "DESeq2"
pgx.testGenes <- function(pgx, method = "ttest") {
  message("pgx.testGenes: method = ", method)
  Y <- pgx$contrasts
  X <- pgx$X
  if (method %in% c("edgeR", "DESeq2")) {
    X <- pgx$counts[rownames(pgx$X), ]
  }

  colnames(Y)
  yref <- sub(".*_vs_", "", colnames(Y))
  yref

  res <- list()
  i <- 2
  for (i in 1:ncol(Y)) {
    y <- Y[, i]
    ii <- which(!is.na(y))
    res[[i]] <- stats.test(X[, ii], y[ii],
      ref = yref[i],
      add.avg = TRUE, method = method
    )
  }
  names(res) <- colnames(Y)
  length(res)

  ## gather matrices
  F <- do.call(cbind, lapply(res, function(x) x$logFC))
  P <- do.call(cbind, lapply(res, function(x) x$pvalue))
  Q <- do.call(cbind, lapply(res, function(x) x$qvalue))
  rownames(F) <- rownames(X)
  rownames(P) <- rownames(X)
  rownames(Q) <- rownames(X)
  dim(F)

  list(F = F, P = P, Q = Q)
}

stats.test <- function(X, y, ref = NULL, add.avg = TRUE, method = "ttest", ...) {
  switch(method,
    "ttest" = stats.ttest(X, y, ref, add.avg, ...),
    "cortest" = stats.cortest(X, y, ref, add.avg, ...),
    "limma" = stats.limma(X, y, ref, add.avg, ...),
    "edgeR" = stats.edgeR(X, y, ref, add.avg, ...),
    "DESeq2" = stats.DESeq2(X, y, ref, add.avg, ...)
  )
}

stats.runAllMethods <- function(X, counts, y, ref = NULL, add.avg = TRUE, ...) {
  results <- list()
  bench <- microbenchmark::microbenchmark(
    "ttest" = {
      results[["ttest"]] <- stats.ttest(X, y, ref, add.avg, ...)
    },
    "cortest" = {
      results[["cortest"]] <- stats.cortest(X, y, ref, add.avg, ...)
    },
    "limma" = {
      results[["limma"]] <- stats.limma(X, y, ref, add.avg, ...)
    },
    "edgeR.QLF" = {
      results[["edgeR.QLF"]] <- stats.edgeR(counts, y, ref, add.avg, test = "QLF", ...)
    },
    "edgeR.LRT" = {
      results[["edgeR.LRT"]] <- stats.edgeR(counts, y, ref, add.avg, test = "LRT", ...)
    },
    "DESeq2.Wald" = {
      results[["DESeq2.Wald"]] <- stats.DESeq2(counts, y, ref, add.avg, test = "Wald", ...)
    },
    "DESeq2.LRT" = {
      results[["DESeq2.LRT"]] <- stats.DESeq2(counts, y, ref, add.avg, test = "LRT", ...)
    },
    times = 1
  )
  list(results = results, timing = bench)
}

stats.ttest <- function(X, y, ref = NULL, add.avg = TRUE) {
  if (is.null(ref)) ref <- sort(unique(y))[1]
  xx <- X[, which(y == ref)]
  yy <- X[, which(y != ref)]
  res <- matrixTests::row_t_welch(xx, yy, conf.level = NA)
  res <- res[rownames(X), ]
  res$qvalue <- p.adjust(res$pvalue, method = "fdr")
  res$mean <- rowMeans(X, na.rm = TRUE)
  res$mean.diff <- -res$mean.diff ## negative!!
  top <- res[, c("mean.diff", "mean", "pvalue", "qvalue", "mean.x", "mean.y")]
  notref <- paste(setdiff(y, ref), collapse = "|")
  colnames(top) <- c("logFC", "mean", "pvalue", "qvalue", ref, notref)
  top
}

stats.cortest <- function(X, y, ref = NULL, add.avg = TRUE) {
  if (length(unique(y[!is.na(y)])) > 2) stop("only 2 levels allowed")
  if (is.null(ref)) ref <- sort(unique(y[!is.na(y)]))[1]
  res <- data.frame(gene = rownames(X))
  ii <- which(!is.na(y) & !colMeans(is.na(X)) == 1)
  X <- X[, ii, drop = FALSE]
  y <- y[ii]
  dy <- (y != ref) / sum(y != ref, na.rm = TRUE) - (y == ref) / sum(y == ref, na.rm = TRUE)
  res$logFC <- (X %*% dy)[, 1]
  res$mean <- rowMeans(X, na.rm = TRUE)
  res$rho <- cor(t(X), dy, use = "pairwise")[, 1]
  res$pvalue <- cor.pvalue(res$rho, ncol(X))
  res$qvalue <- stats::p.adjust(res$pvalue, method = "fdr")
  if (add.avg) {
    res[[ref]] <- rowMeans(X[, y == ref, drop = FALSE], na.rm = TRUE)
    not.ref <- paste(setdiff(y, ref), collapse = "|")
    res[[not.ref]] <- rowMeans(X[, y != ref, drop = FALSE], na.rm = TRUE)
  }
  rownames(res) <- res$gene
  res$gene <- NULL
  res
}

stats.limma <- function(X, y, ref = NULL, test = c("B", "Treat"), add.avg = TRUE,
                        trend = TRUE) {
  if (length(unique(y[!is.na(y)])) > 2) stop("only 2 levels allowed")
  if (ncol(X) != length(y)) stop("ncol(X) not equal length y")
  test <- test[1]
  X <- X[, !is.na(y), drop = FALSE]
  y <- y[!is.na(y)]
  fy <- factor(as.character(y))
  if (!is.null(ref)) fy <- relevel(fy, ref = as.character(ref))
  ref <- levels(fy)[1]
  coldata <- data.frame(y = fy)
  design <- model.matrix(~y, data = coldata)
  lmfit <- limma::lmFit(X, design)
  efit <- limma::eBayes(lmfit, trend = trend) ## trend is good...
  if (test == "Treat") {
    top <- limma::topTreat(efit, sort.by = "none", coef = 2, number = Inf)
  } else {
    top <- limma::topTable(efit, sort.by = "none", coef = 2, number = Inf)
  }
  top <- top[rownames(X), c("logFC", "AveExpr", "P.Value", "adj.P.Val")]
  colnames(top) <- c("logFC", "mean", "pvalue", "qvalue")
  top$logFC <- -top$logFC ## flipped
  if (add.avg) {
    avg <- tapply(
      colnames(X), fy,
      function(k) rowMeans(X[, k, drop = FALSE])
    )
    avg <- do.call(cbind, avg)
    avg <- avg[, c(ref, setdiff(y, ref))] ## reorder ref first
    ## colnames(avg) <- paste0("AveExpr.",colnames(avg))
    top <- cbind(top, avg)
  }
  top
}

stats.edgeR <- function(counts, y, ref = NULL, add.avg = TRUE, test = "QLF", ...) {
  if (length(unique(y[!is.na(y)])) > 2) stop("only 2 levels allowed")
  dge <- edgeR::DGEList(counts)
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  fy <- factor(as.character(y))
  if (!is.null(ref)) fy <- relevel(fy, ref = as.character(ref))
  ref <- levels(fy)[1]
  df <- data.frame(y = fy)
  design <- model.matrix(~y, data = df)
  dge <- edgeR::estimateDisp(dge, design, robust = TRUE)
  if (test == "LRT") {
    fit <- edgeR::glmFit(dge, design, robust = TRUE)
    ## glm <- edgeR::glmLRT(fit, coef=2)
    glm <- edgeR::glmLRT(fit, coef = 2, ...)
  } else if (test == "Treat") {
    fit <- edgeR::glmFit(dge, design, robust = TRUE)
    glm <- edgeR::glmTreat(fit, lfc = log2(1.2), coef = 2, ...)
  } else if (test == "QLF") {
    fit <- edgeR::glmQLFit(dge, design, robust = TRUE)
    glm <- edgeR::glmQLFTest(fit, ...)
  } else {
    stop("[stats.edgeR] unknown test type", test)
  }
  topgenes <- edgeR::topTags(glm, n = Inf, sort.by = "none")
  top <- topgenes$table
  top <- top[rownames(counts), ]
  top <- top[, c("logFC", "logCPM", "PValue", "FDR")]
  colnames(top) <- c("logFC", "mean", "pvalue", "qvalue")
  if (add.avg) {
    cpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
    avg <- tapply(
      colnames(counts), y,
      function(k) rowMeans(cpm[, k, drop = FALSE])
    )
    avg <- do.call(cbind, avg)
    avg <- avg[, c(ref, setdiff(y, ref))] ## reorder ref first
    avg <- avg[rownames(top), ]
    top <- cbind(top, avg)
  }
  top
}

res <- "notact"
stats.DESeq2 <- function(counts, y, ref = NULL, test = "Wald", add.avg = TRUE, ...) {
  if (length(unique(y[!is.na(y)])) > 2) stop("only 2 levels allowed")
  fy <- factor(as.character(y))
  if (!is.null(ref)) fy <- relevel(fy, ref = as.character(ref))
  ref <- levels(fy)[1]
  message("[stats.DESeq2] using refence level: ", ref)
  message(paste("[stats.DESeq2] comparison: ", levels(fy)[2], "vs", levels(fy)[1]))
  coldata <- data.frame(y = fy)
  counts <- round(counts)
  mode(counts) <- "integer"
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData = coldata,
    design = ~y
  )
  if (test == "Wald") {
    dds <- DESeq2::DESeq(dds, quiet = TRUE, parallel = TRUE)
    ## dds <- DESeq2::DESeq(dds, quiet=TRUE, parallel=TRUE, ...)
  } else if (test == "LRT") {
    dds <- DESeq2::DESeq(dds, quiet = TRUE, test = "LRT", reduced = ~1, fitType = "glmGamPoi")
  } else {
    stop("[stats.DESeq2] unknown test type", test)
  }

  DESeq2::resultsNames(dds)
  ## top <- DESeq2::lfcShrink(dds, coef=2, type="normal")  ## recommended with shrinkage
  top <- DESeq2::results(dds)
  top <- top[rownames(counts), ]
  top$logMean <- log2(1 + top$baseMean)
  top <- top[, c("log2FoldChange", "logMean", "pvalue", "padj")]
  colnames(top) <- c("logFC", "mean", "pvalue", "qvalue")
  top <- top[rownames(counts), ]
  if (add.avg) {
    vsd <- DESeq2::vst(dds, blind = FALSE)
    vx <- SummarizedExperiment::assay(vsd)
    avg <- tapply(
      colnames(vx), y,
      function(k) rowMeans(vx[, k, drop = FALSE])
    )
    avg <- do.call(cbind, avg)
    avg <- avg[, c(ref, setdiff(y, ref))] ## REF first
    avg <- avg[rownames(top), ]
    top <- cbind(top, avg)
  }
  top
}

stats.enrichment <- function(F, pattern = "^go|pathway") {
  ## Gene sets
  if (NCOL(F) == 1 && !"Matrix" %in% class(F)) {
    F <- cbind(F)
  }
  is.symbol <- (mean(rownames(F) %in% colnames(playdata::GSETxGENE)) > 0.2)
  if (!is.symbol) {
    stop("features are not SYMBOL")
    return(NULL)
  }
  gmt <- playdata::GSETxGENE
  if (!is.null(pattern)) {
    sel <- grep(pattern, rownames(playdata::GSETxGENE), value = TRUE, ignore.case = TRUE)
    gmt <- gmt[sel, ]
  }
  gsa <- playbase::gset.rankcor(F, Matrix::t(gmt), compute.p = TRUE)
  names(gsa)
  if (NCOL(F) == 1) {
    gsa <- data.frame(rho = gsa$rho[, 1], p.value = gsa$p.value[, 1], q.value = gsa$q.value[, 1])
    gsa <- gsa[order(gsa$p.value), ]
  } else {
    ii <- order(-rowMeans(-log(gsa$p.value)))
    gsa$rho <- gsa$rho[ii, , drop = FALSE]
    gsa$p.value <- gsa$p.value[ii, , drop = FALSE]
    gsa$q.value <- gsa$q.value[ii, , drop = FALSE]
  }
  gsa
}

## lfc=1;q=0.05;set.na = NULL;trend = TRUE;verbose = TRUE
stats.numsig <- function(X, y, lfc = 1, q = 0.05, set.na = NULL,
                         trend = TRUE, verbose = TRUE, gs.method = "rankcor") {
  if (!is.null(set.na)) y[y == set.na] <- NA

  ## select non-missing
  sel <- !is.na(y)
  y <- y[sel]
  X <- X[, sel, drop = FALSE]

  ## Genes
  res <- gx.limmaF(X, y, fdr = 1, lfc = 0, trend = trend, verbose = 0)
  #  res1 <- gx.limma(X, y, fdr=1, lfc=0, trend=trend, verbose=0)
  #  res2 <- stats.limma(X, y, ref=NULL, trend=trend)

  res <- res[order(res$P.Value), ]
  avx <- res[, grep("AveExpr", colnames(res))]
  ldiff <- apply(avx, 1, function(x) diff(range(x)))
  sig <- (abs(ldiff) > lfc & res$adj.P.Val < q)
  sig.genes <- rownames(res)[which(sig)]
  fc0 <- array(ldiff, dimnames = list(rownames(res)))
  names(fc0) <- toupper(names(fc0))

  ## Gene sets
  is.symbol <- (mean(names(fc0) %in% colnames(playdata::GSETxGENE)) > 0.2)
  sig.gsets <- c()
  gsa <- list(rho = NULL, p.value = NULL, q.value = NULL)
  if (is.symbol) {
    sel <- grep("^go|pathway|geo_human", rownames(playdata::GSETxGENE),
      value = TRUE, ignore.case = TRUE
    )
    gmt <- playdata::GSETxGENE[sel, ]
    if (gs.method == "rankcor") {
      gsa <- playbase::gset.rankcor(cbind(fc0), Matrix::t(gmt), compute.p = TRUE)
      gsa <- data.frame(rho = gsa$rho[, 1], p.value = gsa$p.value[, 1], q.value = gsa$q.value[, 1])
    } else if (gs.method == "fgsea") {
      gsets <- playbase::mat2gmt(t(gmt))
      gsa <- fgsea::fgsea(gsets, fc0)
      gsa <- data.frame(rho = gsa$NES, p.value = gsa$pval, q.value = gsa$padj)
    } else {
      stop("unknown geneset method: ", gs.method)
    }
    gsa <- gsa[order(gsa$p.value), ]
    sig.gs <- (gsa$q.value < q)
    sig.gsets <- rownames(gsa)[which(sig.gs)]
  }

  if (verbose) {
    cat("nsig.genes = ", length(sig.genes), "\n")
    cat("nsig.gsets = ", length(sig.gsets), "\n")
  }

  list(
    genes = sig.genes,
    gsets = sig.gsets,
    fc = fc0,
    p.genes = res$P.Value,
    p.gsets = gsa$p.value
  )
}
