##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

## -----------------------------------------------------------------------------
## -------------------- FIT ALL CONTRASTS --------------------------------------
## -----------------------------------------------------------------------------


#' @title ngs.fitContrastsWithAllMethods
#'
#' @description This function fits contrasts using multiple differential expression analysis methods on count data.
#' It applies the methods specified in the \code{methods} argument to the count data in \code{counts}, design matrix \code{design}, and contrast matrix \code{contr.matrix}.
#'
#' @param counts A matrix of counts, with genes in rows and samples in columns.
#' @param X Covariates to include in the design matrix. Default is NULL.
#' @param samples A vector of sample names that match the columns in \code{counts}.
#' @param design The design matrix, with samples in columns.
#' @param covariates Variables to be regressed out. Only valid for linear model-based DGE tests. Default is NULL.
#' @param contr.matrix The contrasts matrix, with contrasts in rows.
#' @param genes A vector of gene names that match the rows in \code{counts}. Default is NULL.
#' @param prior.cpm Prior counts per million. Default is 1.
#' @param prune.samples Whether to filter low count samples. Default is FALSE.
#' @param conform.output Whether to conform the output. Default is TRUE.
#' @param do.filter Whether to filter genes. Default is TRUE.
#' @param cpm.scale Counts per million scaling factor. Default is 1e6.
#' @param remove.batch Whether to remove batch effects. Default is TRUE.
#' @param methods Methods to apply. Can be one or more of:
#'   \itemize{
#'     \item{\code{ttest}}{t-test}
#'     \item{\code{ttest.welch}}{Welch's t-test}
#'     \item{\code{voom.limma}}{voom + limma}
#'     \item{\code{trend.limma}}{limma with trend}
#'     \item{\code{notrend.limma}}{limma without trend}
#'     \item{\code{deseq2.wald}}{DESeq2 Wald test}
#'     \item{\code{deseq2.lrt}}{DESeq2 LRT}
#'     \item{\code{edger.qlf}}{edgeR QLF}
#'     \item{\code{edger.lrt}}{edgeR LRT}
#'   }
#' @param correct.AveExpr Whether to correct for average expression. Default is TRUE.
#' @param custom Custom differential expression method. Default is NULL.
#' @param custom.name Name for custom method. Default is NULL.
#' @param time time character vector to be passed for time series analysis. Default NULL.
#'
#' @details This function provides a convenient wrapper to run multiple differential expression methods on the same data.
#' It runs the methods specified in \code{methods} on the provided count data and design.
#' The results are returned in a list, with each element containing the results from one method.
#' Filtering, batch effect removal, and output conforming can be controlled with arguments.
#' Custom methods can also be added.
#'
#' @return A list with results from each differential expression method.
#'
#' @examples
#' \dontrun{
#' # TODO
#' }
#' @export
ngs.fitContrastsWithAllMethods <- function(counts,
                                           X = NULL,
                                           samples,
                                           design,
                                           covariates = NULL,
                                           contr.matrix,
                                           genes = NULL,
                                           prior.cpm = 1,
                                           cpm.scale = 1e6,
                                           prune.samples = TRUE,
                                           conform.output = TRUE,
                                           do.filter = TRUE,
                                           remove.batch = TRUE,
                                           methods = c(
                                             "ttest", "ttest.welch", "wilcoxon.ranksum",
                                             "voom.limma", "trend.limma", "notrend.limma",
                                             "deseq2.wald", "deseq2.lrt", "edger.qlf", "edger.lrt"
                                           ),
                                           correct.AveExpr = TRUE,
                                           custom = NULL,
                                           custom.name = NULL,
                                           timeseries = FALSE) {

  ## Don't test fully missing features. Put them back in toptable.
  counts <- counts[which(rowMeans(is.na(counts)) < 1), ]
  if (!is.null(X)) X <- X[which(rowMeans(is.na(X)) < 1), ]

  if (methods[1] == "*") {
    methods <- c(
      "ttest", "ttest.welch", "wilcoxon.ranksum",
      "voom.limma", "trend.limma", "notrend.limma",
      "deseq2.wald", "deseq2.lrt", "edger.qlf", "edger.lrt"
    )
  }
  methods <- intersect(methods, c(
    "ttest", "ttest.welch", "wilcoxon.ranksum",
    "voom.limma", "trend.limma", "notrend.limma",
    "deseq2.wald", "deseq2.lrt", "edger.qlf", "edger.lrt"
  ))

  ## If degenerate set design to NULL
  if (!is.null(design) && ncol(design) >= ncol(X)) {
    ## "no-replicate" design!!!
    cat("WARNING: degenerate design. setting design to NULL\n")
    contr.matrix <- design %*% contr.matrix
    design <- NULL
  }

  counts <- pmax(counts, 0)

  ## -----------------------------------------------------------------------------
  ## Time series: determine variable 'time'
  ## -----------------------------------------------------------------------------
  if (timeseries) {
    time.var <- playbase::get_timevars()
    sel <- grep(time.var, colnames(samples), ignore.case = TRUE)
    timeseries <- as.character(samples[, sel[1]])
    names(timeseries) <- rownames(samples)
  } else {
    timeseries <- NULL
  }

  ## -----------------------------------------------------------------
  ## Covariates
  ## -----------------------------------------------------------------
  if (!is.null(covariates)) {
    kk <- intersect(covariates, colnames(samples))
    if (length(kk) > 0) {
      covariates <- samples[, kk, drop = FALSE]
    } else {
      covariates <- NULL
    }
  }
  
  ## ------------------------------------------------------------------
  ## define transformation methods: log2CPM for counts
  ## ------------------------------------------------------------------
  if (is.null(X)) {
    message("[ngs.fitContrastsWithAllMethods] prior CPM counts =", prior.cpm)
    message("[ngs.fitContrastsWithAllMethods] CPM scale =", cpm.scale)
    X <- log2(t(t(counts) / Matrix::colSums(counts)) * cpm.scale + prior.cpm) ## log2CPM
    X <- limma::normalizeQuantiles(X) ## in linear space
  } else {
    message("[ngs.fitContrastsWithAllMethods] Using input log-expression matrix X as is")
  }

  ## ------------------------------------------------------------------
  ## get main grouping variable for modeling
  ## ------------------------------------------------------------------
  group <- NULL
  if (!is.null(design)) {
    group <- colnames(design)[max.col(design)]
    if (nrow(design) == ncol(design) &&
      all(rownames(design) == colnames(design))) {
      group <- NULL
    }
  }

  timings <- list()
  outputs <- list()

  ## Skip tests that do not tolerate NAs. Inform the user.
  nmissing.counts <- sum(is.na(counts))
  nmissing <- sum(is.na(X))
  if (nmissing.counts > 0 & nmissing == 0) counts <- pmax(2**X - 1, 0)

  ## ---------------- t-test methods -------------------
  ttest.mtds <- c("ttest", "ttest.rank", "ttest.welch")
  ttest.mdls <- c("equalvar", "equalvar", "welch")
  cm.mtds <- intersect(methods, ttest.mtds)
  if (length(cm.mtds) > 0) {
    for (i in 1:length(cm.mtds)) {
      message("[ngs.fitContrastsWithAllMethods] Fitting using ", cm.mtds[i])
      X1 <- X
      if (cm.mtds[i] == "ttest.rank") {
        X1 <- scale(apply(X1, 2, rank, na.last = "keep"))
      }
      mdl <- ttest.mdls[match(cm.mtds[i], ttest.mtds)]
      timings[[cm.mtds[i]]] <- system.time(
        outputs[[cm.mtds[i]]] <- ngs.fitContrastsWithTTEST(
          X1, contr.matrix, design,
          method = mdl,
          conform.output = conform.output
        )
      )
    }
  }

  ## ---------- Wilcoxon test -----------
  if ("wilcoxon.ranksum" %in% methods) {
    message("[ngs.fitContrastsWithAllMethods] Fitting using Wilcoxon rank sum test")
    timings[["wilcoxon.ranksum"]] <- system.time(
      outputs[["wilcoxon.ranksum"]] <- ngs.fitContrastsWithWILCOXON(
        X, contr.matrix, design,
        conform.output = conform.output
      )
    )
  }
  
  ## ---------------- LIMMA methods -------------------
  limma.mtds <- c("trend.limma", "notrend.limma", "voom.limma")
  limma.mdls <- c("limma", "limma", "voom")
  cm.mtds <- intersect(methods, limma.mtds)
  if (length(cm.mtds) > 0) {
    for (i in 1:length(cm.mtds)) {
      X1 <- X
      mdl <- limma.mdls[match(cm.mtds[i], limma.mtds)]
      trend <- ifelse(grepl("notrend|voom", cm.mtds[i]), FALSE, TRUE)
      message("[ngs.fitContrastsWithAllMethods] Fitting using ", cm.mtds[i])
      if (cm.mtds[i] == "voom.limma" && nmissing > 0) {
        message("[ngs.fitContrastsWithAllMethods] Missing values detected. Cannot perform voom.limma")
        next
      } else {
        if (!is.null(timeseries) && cm.mtds[i] == "trend.limma") {
          time_var <- timeseries
        } else {
          time_var <- NULL
        }
        tt <- system.time(
          outputs[[cm.mtds[i]]] <- ngs.fitContrastsWithLIMMA(
            X1, contr.matrix,
            design, covariates = covariates,
            method = mdl,
            trend = trend, prune.samples = prune.samples,
            conform.output = conform.output, plot = FALSE, timeseries = time_var
          )
        )
        timings[[cm.mtds[i]]] <- round(as.numeric(tt), digits = 4)
      }
    }
  }

  ## ---------------- DESEQ2 methods -------------------
  deseq2.mtds <- c("deseq2.wald", "deseq2.lrt")
  deseq2.mdls <- c("Wald", "LRT")
  cm.mtds <- intersect(methods, deseq2.mtds)
  if (length(cm.mtds) > 0) {
    for (i in 1:length(cm.mtds)) {
      X1 <- X
      mdl <- deseq2.mdls[match(cm.mtds[i], deseq2.mtds)]
      message("[ngs.fitContrastsWithAllMethods] Fitting using ", cm.mtds[i])
      if (nmissing.counts > 0) {
        message("[ngs.fitContrastsWithAllMethods] Missing values detected. Cannot perform ", cm.mtds[i])
        next
      } else {
        if (!is.null(timeseries)) {
          time_var <- timeseries
        } else {
          time_var <- NULL
        }
        timings[[cm.mtds[i]]] <- system.time(
          outputs[[cm.mtds[i]]] <- ngs.fitContrastsWithDESEQ2(
            counts, group, contr.matrix,
            design, covariates = covariates,
            X = X1, genes = genes, test = mdl,
            prune.samples = prune.samples,
            conform.output = conform.output, timeseries = time_var
          )
        )
      }
    }
  }

  ## ---------------- EdgeR methods -------------------
  edger.mtds <- c("edger.qlf", "edger.lrt")
  edger.mdls <- c("qlf", "lrt")
  cm.mtds <- intersect(methods, edger.mtds)
  if (length(cm.mtds) > 0) {
    for (i in 1:length(cm.mtds)) {
      X1 <- X
      mdl <- edger.mdls[match(cm.mtds[i], edger.mtds)]
      message("[ngs.fitContrastsWithAllMethods] Fitting using ", cm.mtds[i])
      if (nmissing.counts > 0) {
        message("[ngs.fitContrastsWithAllMethods] Missing values detected. Cannot perform edgeR QL-F test or LRT.")
        next
      } else {
        time_var <- NULL
        if (!is.null(timeseries)) {
          time_var <- timeseries
        }
        timings[[cm.mtds[i]]] <- system.time(
          outputs[[cm.mtds[i]]] <- ngs.fitContrastsWithEDGER(
            counts, group, contr.matrix,
            design, covariates = covariates,
            method = mdl,
            X = X1, prune.samples = prune.samples,
            conform.output = conform.output, plot = FALSE, timeseries = time_var
          )
        )
      }
    }
  }

  if (!is.null(custom)) {
    message("[ngs.fitContrastsWithAllMethods] adding custom results table")
    if (is.null(custom.name)) custom.name <- "custom"
    if (length(outputs) == 0) {
      compare_output <- ngs.fitContrastsWithTTEST(
        X, contr.matrix, design,
        method = "equalvar",
        conform.output = conform.output
      )
    } else {
      compare_output <- outputs[[1]]
    }
    test_names <- unique(sub("\\..*", "", colnames(custom)))
    test_names <- test_names[match(names(compare_output$tables), test_names)]
    custom_tables <- list()
    for (test in test_names) {
      logFC_col <- paste0(test, ".logFC")
      pval_col <- paste0(test, ".P.Value")
      adjp_col <- paste0(test, ".adj.P.Val")
      test_df <- data.frame(
        logFC = custom[, logFC_col],
        P.Value = custom[, pval_col],
        adj.P.Val = custom[, adjp_col]
      )
      rownames(test_df) <- rownames(custom)
      custom_tables[[test]] <- test_df
    }
    custom <- custom_tables
    missing_rows <- rownames(compare_output$tables[[1]])[which(!rownames(compare_output$tables[[1]]) %in% rownames(custom[[1]]))]
    if (length(missing_rows) > 0) {
      missing_data <- matrix(NA,
        nrow = length(missing_rows),
        ncol = ncol(custom[[1]]),
        dimnames = list(missing_rows, colnames(custom[[1]]))
      )
      for (test in names(custom)) {
        custom[[test]] <- rbind(custom[[test]], missing_data)
        custom[[test]] <- custom[[test]][rownames(compare_output$tables[[1]]), ]
      }
    }
    outputs[[custom.name]]$tables <- custom
    timings[["custom"]] <- system.time(0)
  }

  ##-------temp
  saveRDS(list(outputs=outputs, contr.matrix=contr.matrix), "~/Desktop/oo.RDS")
  ##--------
  
  
  ## ----------------------------------------------------------------------
  ## "corrections" ...
  ## ----------------------------------------------------------------------
  correct.AveExpr <- TRUE
  if (correct.AveExpr) {
    message("[ngs.fitContrastsWithAllMethods] correcting AveExpr values...")

    ## EdgeR and Deseq2 compute weird normalized expr. matrix.
    ## We need to "correct" for those.

    exp.matrix <- contr.matrix

    if (!is.null(design)) exp.matrix <- (design %*% contr.matrix)
    samplesX <- lapply(apply(exp.matrix != 0, 2, which, simplify = FALSE), function(i) rownames(exp.matrix)[i])
    samples1 <- lapply(apply(exp.matrix > 0, 2, which, simplify = FALSE), function(i) rownames(exp.matrix)[i])
    samples0 <- lapply(apply(exp.matrix < 0, 2, which, simplify = FALSE), function(i) rownames(exp.matrix)[i])

    avgX <- sapply(samplesX, function(s) rowMeans(X[, s, drop = FALSE], na.rm = TRUE))
    avg.1 <- sapply(samples1, function(s) rowMeans(X[, s, drop = FALSE], na.rm = TRUE))
    avg.0 <- sapply(samples0, function(s) rowMeans(X[, s, drop = FALSE], na.rm = TRUE))

    i <- j <- 1
    for (i in 1:length(outputs)) {
      for (j in 1:length(outputs[[i]]$tables)) {
        outputs[[i]]$tables[[j]]$AveExpr <- avgX[, j]
        outputs[[i]]$tables[[j]]$AveExpr1 <- avg.1[, j]
        outputs[[i]]$tables[[j]]$AveExpr0 <- avg.0[, j]
      }
    }
  }

  ## -----------------------------------------------------------------------
  ## Put "IA:*" contrasts as last columns in outputs limma/DeSeq2/EdgeR. NEEDED??
  ## -----------------------------------------------------------------------
  ## if (!is.null(timeseries)) {
  ##   i <- 1
  ##   for (i in 1:length(outputs)) {
  ##     contr.names <- names(outputs[[i]]$tables)
  ##     chk1 <- grep("^IA:*", contr.names)
  ##     if (any(chk1)) {
  ##       ss <- unique(c(contr.names[-chk1], contr.names[chk1]))
  ##       jj <- match(ss, names(outputs[[i]]$tables))
  ##       outputs[[i]]$tables <- outputs[[i]]$tables[jj]
  ##       names(outputs[[i]]$tables) <- ss
  ##     }
  ##   }
  ## }

  ## ----------------------------------------------------------------------
  ## add some statistics
  ## ----------------------------------------------------------------------
  message("[ngs.fitContrastsWithAllMethods] calculating statistics...")

  for (i in 1:length(outputs)) {
    res <- outputs[[i]]
    M <- sapply(res$tables, function(x) x[, "AveExpr"]) ## !!!! edgeR and Deseq2 are weird!!!
    M0 <- sapply(res$tables, function(x) x[, "AveExpr0"])
    M1 <- sapply(res$tables, function(x) x[, "AveExpr1"])
    Q <- sapply(res$tables, function(x) x[, "adj.P.Val"])
    P <- sapply(res$tables, function(x) x[, "P.Value"])
    logFC <- sapply(res$tables, function(x) x[, "logFC"])
    # colnames(M) <- colnames(logFC) <- colnames(P) <- colnames(Q) <- colnames(contr.matrix)
    rownames(M) <- rownames(logFC) <- rownames(P) <- rownames(Q) <- rownames(res$tables[[1]])
    rownames(M0) <- rownames(M1) <- rownames(res$tables[[1]])

    ## count significant terms
    qvalues <- c(1e-16, 10**seq(-8, -2, 2), 0.05, 0.1, 0.2, 0.5, 1)
    lfc <- 1
    sig.both <- sapply(qvalues, function(q) Matrix::colSums((Q <= q) * (abs(logFC) > lfc), na.rm = TRUE))
    sig.up <- sapply(qvalues, function(q) Matrix::colSums((Q <= q) * (logFC > lfc), na.rm = TRUE))
    sig.down <- sapply(qvalues, function(q) Matrix::colSums((Q <= q) * (logFC < -lfc), na.rm = TRUE))
    sig.notsig <- sapply(qvalues, function(q) Matrix::colSums(Q > q | (abs(logFC) < lfc), na.rm = TRUE))
    if (NCOL(Q) == 1) {
      sig.both <- matrix(sig.both, nrow = 1)
      sig.up <- matrix(sig.up, nrow = 1)
      sig.down <- matrix(sig.down, nrow = 1)
      sig.notsig <- matrix(sig.notsig, nrow = 1)
      rownames(sig.both) <- rownames(sig.up) <- colnames(Q)
      rownames(sig.down) <- rownames(sig.notsig) <- colnames(Q)
    }
    colnames(sig.both) <- colnames(sig.notsig) <- qvalues
    colnames(sig.up) <- colnames(sig.down) <- qvalues

    res$sig.counts <- list(both = sig.both, up = sig.up, down = sig.down, notsig = sig.notsig)

    ## need this? takes space!!!
    res$p.value <- P
    res$q.value <- Q
    res$logFC <- logFC
    res$aveExpr <- M
    res$aveExpr0 <- M0
    res$aveExpr1 <- M1

    outputs[[i]] <- res
  }

  ## --------------------------------------------------------------
  ## Reshape matrices by comparison
  ## --------------------------------------------------------------
  message("[ngs.fitContrastsWithAllMethods] Reshape matrices...")

  tests <- colnames(outputs[[1]]$p.value)
  ntest <- length(tests)
  P <- lapply(1:ntest, function(i) sapply(outputs, function(x) x$p.value[, i]))
  Q <- lapply(1:ntest, function(i) sapply(outputs, function(x) x$q.value[, i]))
  logFC <- lapply(1:ntest, function(i) sapply(outputs, function(x) x$logFC[, i]))
  aveExpr <- lapply(1:ntest, function(i) sapply(outputs, function(x) x$aveExpr[, i]))
  aveExpr0 <- lapply(1:ntest, function(i) sapply(outputs, function(x) x$aveExpr0[, i]))
  aveExpr1 <- lapply(1:ntest, function(i) sapply(outputs, function(x) x$aveExpr1[, i]))

  sig.up <- lapply(1:ntest, function(i) {
    do.call(rbind, lapply(outputs, function(x) x$sig.counts[["up"]][i, ]))
  })

  sig.down <- lapply(1:ntest, function(i) {
    do.call(rbind, lapply(outputs, function(x) x$sig.counts[["down"]][i, ]))
  })

  sig.notsig <- lapply(1:ntest, function(i) {
    do.call(rbind, lapply(outputs, function(x) x$sig.counts[["notsig"]][i, ]))
  })

  sig.both <- lapply(1:ntest, function(i) {
    do.call(rbind, lapply(outputs, function(x) x$sig.counts[["both"]][i, ]))
  })

  names(P) <- names(Q) <- names(logFC) <- names(aveExpr) <- tests
  names(sig.up) <- names(sig.down) <- names(sig.both) <- names(sig.notsig) <- tests
  sig.counts <- list(up = sig.up, down = sig.down, both = sig.both, notsig = sig.notsig)

  ## --------------------------------------------------
  ## meta analysis, aggregate p-values
  ## --------------------------------------------------
  message("[ngs.fitContrastsWithAllMethods] aggregating p-values...")

  all.meta <- list()
  for (i in 1:ntest) {
    pv <- P[[i]]
    qv <- Q[[i]]
    fc <- logFC[[i]]
    mx <- aveExpr[[i]]
    mx0 <- aveExpr0[[i]]
    mx1 <- aveExpr1[[i]]

    ## avoid strange values
    fc[is.infinite(fc) | is.nan(fc)] <- NA
    pv <- pmax(pv, 1e-99)
    pv[is.na(pv)] <- 1 # ???
    qv[is.na(qv)] <- 1 # ????

    ## !!!!!!!!!!!!!!!!!!!!!!!! NEED RETHINK !!!!!!!!!!!!!!!!!!!!!!!!
    meta.p <- apply(pv, 1, function(p) exp(mean(log(p), na.rm = TRUE))) ## geometric mean
    meta.q <- stats::p.adjust(meta.p, method = "BH")
    meta.fx <- rowMeans(fc, na.rm = TRUE)
    meta.avg <- rowMeans(mx, na.rm = TRUE)
    meta.avg0 <- rowMeans(mx0, na.rm = TRUE)
    meta.avg1 <- rowMeans(mx1, na.rm = TRUE)

    meta <- data.frame(fx = meta.fx, p = meta.p, q = meta.q)
    avg <- data.frame(avg.0 = meta.avg0, avg.1 = meta.avg1)
    rownames(meta) <- rownames(logFC[[i]])
    rownames(avg) <- rownames(logFC[[i]])
    rownames(fc) <- NULL ## saves memory
    rownames(pv) <- NULL
    rownames(qv) <- NULL
    all.meta[[i]] <- data.frame(meta = meta, avg, fc = I(fc), p = I(pv), q = I(qv))
    if (!is.null(genes)) all.meta[[i]] <- cbind(genes, all.meta[[i]])
  }
  names(all.meta) <- tests

  timings0 <- do.call(rbind, timings)
  colnames(timings0) <- names(timings[[1]])
  colnames(timings0) <- c("user.self", "sys.self", "elapsed", "user.child", "sys.child")

  res <- list(
    outputs = outputs,
    meta = all.meta,
    sig.counts = sig.counts,
    timings = timings0,
    X = X
  )

  return(res)

}

## --------------------------------------------------------------------------------------------
## ----------------------------------- FIT ALL CONTRASTS --------------------------------------
## --------------------------------------------------------------------------------------------

#' @describeIn ngs.fitContrastsWithAllMethods Fits contrasts using t-test
#' @export
ngs.fitContrastsWithTTEST <- function(X,
                                      contr.matrix,
                                      design, method = "welch",
                                      conform.output = 0) {

  tables <- list()
  exp.matrix <- contr.matrix
  if (!is.null(design)) exp.matrix <- (design %*% contr.matrix)
  for (i in 1:ncol(exp.matrix)) {
    j1 <- which(exp.matrix[, i] > 0)
    j0 <- which(exp.matrix[, i] < 0)
    if (method == "welch") {
      suppressWarnings(rt <- matrixTests::row_t_welch(X[, j1, drop = FALSE], X[, j0, drop = FALSE]))
    } else {
      suppressWarnings(rt <- matrixTests::row_t_equalvar(X[, j1, drop = FALSE], X[, j0, drop = FALSE]))
    }
    kk <- c("mean.x", "mean.y", "mean.diff", "stderr", "df", "statistic", "pvalue")
    rt <- rt[, kk]
    rt$qvalue <- stats::p.adjust(rt[, "pvalue"], method = "BH")
    rt$mean.value <- (rt[, "mean.x"] + rt[, "mean.y"]) / 2
    tables[[i]] <- rt
  }
  names(tables) <- colnames(contr.matrix)
  if (conform.output == TRUE) {
    for (i in 1:length(tables)) {
      k1 <- c("mean.diff", "mean.value", "statistic", "pvalue", "qvalue", "mean.y", "mean.x")
      k2 <- c("logFC", "AveExpr", "statistic", "P.Value", "adj.P.Val", "AveExpr0", "AveExpr1")
      tables[[i]] <- tables[[i]][, k1]
      colnames(tables[[i]]) <- k2
    }
  }
  res <- list(tables = tables)
  return(res)
}

#' @describeIn ngs.fitContrastsWithAllMethods Fits contrasts using Wilcoxon Rank Sum test
#' @export
ngs.fitContrastsWithWILCOXON <- function(X, contr.matrix, design, conform.output = 0) {
  exp.matrix <- contr.matrix
  if (!is.null(design)) exp.matrix <- (design %*% contr.matrix)

  tables <- list()
  for (i in 1:ncol(exp.matrix)) {
    j1 <- which(exp.matrix[, i] > 0)
    j0 <- which(exp.matrix[, i] < 0)
    suppressWarnings(rt <- matrixTests::row_wilcoxon_twosample(
      X[, j1, drop = FALSE],
      X[, j0, drop = FALSE]
    ))

    rt$mean.x <- rowMeans(X[, j1, drop = FALSE], na.rm = TRUE)
    rt$mean.y <- rowMeans(X[, j0, drop = FALSE], na.rm = TRUE)
    rt$mean.diff <- rt$mean.x - rt$mean.y

    varx1 <- apply(X[, j1, drop = FALSE], 1, function(x) var(x) / length(x))
    varx0 <- apply(X[, j0, drop = FALSE], 1, function(x) var(x) / length(x))
    rt$stderr <- sqrt(varx1 + varx0)

    kk <- c("mean.x", "mean.y", "mean.diff", "stderr", "df", "statistic", "pvalue")
    kk <- intersect(kk, colnames(rt))
    rt <- rt[, kk, drop = FALSE]
    rt$qvalue <- stats::p.adjust(rt[, "pvalue"], method = "BH")
    rt$mean.value <- (rt[, "mean.x"] + rt[, "mean.y"]) / 2
    tables[[i]] <- rt
  }

  names(tables) <- colnames(contr.matrix)
  if (conform.output == TRUE) {
    for (i in 1:length(tables)) {
      k1 <- c("mean.diff", "mean.value", "statistic", "pvalue", "qvalue", "mean.y", "mean.x")
      k2 <- c("logFC", "AveExpr", "statistic", "P.Value", "adj.P.Val", "AveExpr0", "AveExpr1")
      tables[[i]] <- tables[[i]][, k1]
      colnames(tables[[i]]) <- k2
    }
  }

  res <- list(tables = tables)
  return(res)
}

#' @title ngs.fitContrastsWithLIMMA
#' @param X log2- and normalized expression matrix
#' @description ngs.fitContrastsWithAllMethods Fits contrasts using LIMMA
#' @return limma top tables
#' @export
#' @name ngs.fitContrastsWithLIMMA
ngs.fitContrastsWithLIMMA <- function(X,
                                      contr.matrix,                                    
                                      design,
                                      covariates = NULL,
                                      method = c("voom", "limma"),
                                      trend = TRUE,
                                      robust = TRUE,
                                      prune.samples = FALSE,
                                      conform.output = FALSE,
                                      plot = FALSE,
                                      timeseries = NULL) {

  ## Design (grouping): perform LIMMA on the entire contrast matrix.
  ## No design (no grouping): perform LIMMA per contrast one-by-one.

  ## Do not test full NA features. Put them back in topTable.
  if (!is.null(X)) X <- X[which(rowMeans(is.na(X)) < 1), ]
  method <- method[1]

  if (!is.null(design)) {
    message("[ngs.fitContrastsWithLIMMA] fitting LIMMA contrasts using design matrix")
    exp0 <- design %*% contr.matrix
    kk <- rownames(exp0)
    if (prune.samples) kk <- rownames(exp0)[which(rowSums(abs(exp0), na.rm = TRUE) > 0)]
    design1 <- design[kk, , drop = FALSE]
    X1 <- X[, kk, drop = FALSE]
    contr1 <- contr.matrix
    if (method == "voom") {
      v <- limma::voom(2**X1, design1, plot = plot, normalize.method = "none")
      vfit <- limma::lmFit(v, design1)
      trend <- FALSE ## no need
    } else {
      vfit <- limma::lmFit(X1, design1)
    }

    vfit <- limma::contrasts.fit(vfit, contrasts = contr1)
    efit <- limma::eBayes(vfit, trend = trend, robust = robust)
    if (plot == TRUE) limma::plotSA(efit)
    tables <- list()
    exp1 <- (design1 %*% contr1)
    for (i in 1:ncol(contr1)) {
      top <- limma::topTable(efit, coef = i, sort.by = "none", number = Inf, adjust.method = "BH")
      j1 <- which(exp1[, i] > 0)
      j2 <- which(exp1[, i] < 0)
      mean1 <- rowMeans(X1[, j1, drop = FALSE], na.rm = TRUE)
      mean0 <- rowMeans(X1[, j2, drop = FALSE], na.rm = TRUE)
      top <- top[rownames(X1), ]
      top <- cbind(top, "AveExpr0" = mean0, "AveExpr1" = mean1)
      tables[[i]] <- top
      names(tables)[i] <- colnames(contr1)[i]
    }
  } else {
    message("[ngs.fitContrastsWithLIMMA] Fitting LIMMA contrasts *without* design")
    exp0 <- contr.matrix ## sample-wise contrasts...
    tables <- list()
    for (i in 1:ncol(exp0)) {
      kk <- 1:nrow(exp0)
      if (prune.samples) kk <- which(!is.na(exp0[, i]) & exp0[, i] != 0)
      ct <- exp0[kk, i]
      y <- factor(c("neg", "o", "pos")[2 + sign(ct)])
      X1 <- X[, kk, drop = FALSE]
      if (grepl("^IA:*", colnames(exp0)[i]) && !is.null(timeseries)) {
        top <- try(ngs.fitContrastsWithLIMMA.timeseries(X1, y, timeseries, trend = trend), silent = TRUE)
      } else {
        design1 <- stats::model.matrix(~ 0 + y)
        if (method == "voom") {
          v <- limma::voom(2**X1, design1, plot = FALSE)
          suppressMessages(vfit <- limma::lmFit(v, design1))
          trend <- FALSE
        } else {
          suppressMessages(vfit <- limma::lmFit(X1, design1))
        }
        contr1 <- matrix(c(-1, 0, 1), nrow = 3)
        rownames(contr1) <- c("yneg", "yo", "ypos")
        colnames(contr1) <- "pos_vs_neg"
        contr1 <- contr1[colnames(vfit), , drop = FALSE]
        vfit <- limma::contrasts.fit(vfit, contrasts = contr1)
        efit <- try(limma::eBayes(vfit, trend = trend, robust = robust), silent = TRUE)
        if ("try-error" %in% class(efit)) {
          efit <- try(limma::eBayes(vfit, trend = trend, robust = FALSE), silent = TRUE)
          if ("try-error" %in% class(efit)) {
            efit <- try(limma::eBayes(vfit, trend = FALSE, robust = FALSE), silent = TRUE)
          }
        }
        top <- try(limma::topTable(efit, coef = 1, sort.by = "none", number = Inf, adjust.method = "BH"), silent = TRUE)
        if ("try-error" %in% class(top)) {
          top <- data.frame(matrix(NA, nrow = nrow(X1), ncol = 6))
          rownames(top) <- rownames(X1)
          colnames(top) <- c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
        } else {
          top <- top[rownames(X1), , drop = FALSE]        
        }

        ##---If covariates specified, regress them out and get p.value
        cov.pval <- list()
        if (!is.null(covariates)) {
          covariates1 <- covariates[kk, , drop = FALSE]
          k=1
          for (k in 1:ncol(covariates1)) {
            cov.name <- colnames(covariates1)[k]
            cov.val <- covariates1[, k]
            if (is.character(cov.val)) {
              cov.val <- factor(cov.val)
            } else if (is.numeric(cov.val)) {
              nn <- unique(cov.val[!is.na(cov.val)])
              if (length(nn) <= 4) { ## edge case: likely categorical
                cov.val <- factor(cov.val)
              } else {
                cov.val <- as.numeric(cov.val) ## no arbitrary binning
              }
            }
            model_data <- data.frame(y = y, cov = cov.val)
            design1_cov <- stats::model.matrix(~ 0 + y + cov, data = model_data)
            suppressMessages(vfit_cov <- limma::lmFit(X1, design1_cov))
            contr1_cov <- matrix(0, nrow = ncol(design1_cov), ncol = 1)
            rownames(contr1_cov) <- colnames(design1_cov)
            colnames(contr1_cov) <- "pos_vs_neg"
            y_cols <- intersect(c("yneg", "yo", "ypos"), colnames(design1_cov))
            contr1_cov[y_cols, 1] <- c(-1, 0, 1)[match(y_cols, c("yneg", "yo", "ypos"))]
            vfit_cov <- limma::contrasts.fit(vfit_cov, contrasts = contr1_cov)
            efit_cov <- try(limma::eBayes(vfit_cov, trend = trend_cov, robust = robust), silent = TRUE)
            if ("try-error" %in% class(efit_cov)) {
              efit_cov <- try(limma::eBayes(vfit_cov, trend = trend_cov, robust = FALSE), silent = TRUE)
              if ("try-error" %in% class(efit_cov)) {
                efit_cov <- try(limma::eBayes(vfit_cov, trend = FALSE, robust = FALSE), silent = TRUE)
              }
            }
            top_cov <- try(limma::topTable(efit_cov, coef = 1, sort.by = "none", number = Inf), silent = TRUE)
            if ("try-error" %in% class(top_cov)) {
              top_cov <- data.frame(matrix(NA, nrow = nrow(X1), ncol = 1))
              rownames(top_cov) <- rownames(X1)
              colnames(top_cov) <- paste0("P.Value.", cov.name)          
              cov.pval[[cov.name]] <- top_cov
            } else {
              top_cov <- top_cov[rownames(X1), "P.Value", drop = FALSE]
              colnames(top_cov) <- paste0("P.Value.", cov.name) 
              cov.pval[[cov.name]] <- top_cov
            } 
          }

        }
        ##--end covariates regression
      }

      j1 <- which(ct > 0)
      j0 <- which(ct < 0)
      mean1 <- rowMeans(X1[, j1, drop = FALSE], na.rm = TRUE)
      mean0 <- rowMeans(X1[, j0, drop = FALSE], na.rm = TRUE)

      top <- cbind(top, "AveExpr0" = mean0, "AveExpr1" = mean1)
      if (length(cov.pval) > 0) { ## add covariates' regression p.values    
        top_cov <- do.call(cbind, cov.pval)
        cm <- intersect(rownames(top_cov), rownames(top))
        top <- cbind(top[cm, , drop = FALSE], top_cov[cm, , drop = FALSE])
      }

      tables[[i]] <- top
      names(tables)[i] <- colnames(exp0)[i]

    }

  }

  if (conform.output) {
    for (i in 1:length(tables)) {
      k1 <- c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "AveExpr0", "AveExpr1")
      k2 <- c("logFC", "AveExpr", "statistic", "P.Value", "adj.P.Val", "AveExpr0", "AveExpr1")
      ec <- grep("^P\\.Value\\.", colnames(tables[[i]]), value = TRUE)
      tables[[i]] <- tables[[i]][, c(k1, ec)]
      colnames(tables[[i]])[1:length(k2)] <- k2
    }
  }

  return(list(tables = tables))

}


#' @describeIn ngs.fitContrastsWithLIMMA Fits contrasts using LIMMA with no design. For time-series analysis.
#' @export
ngs.fitContrastsWithLIMMA.timeseries <- function(X,
                                                 y,
                                                 timeseries,
                                                 trend = TRUE,
                                                 use.spline = NULL) {

  message("[ngs.fitContrastsWithLIMMA.timeseries] Fitting Limma with no design; time series analysis...")
  if (!all(colnames(X) %in% names(timeseries))) {
    stop("[ngs.fitContrastsWithLIMMA.timeseries] X and timeseries vector contain different set of samples.")
  }
  
  jj <- match(colnames(X), names(timeseries))
  time0 <- as.character(unname(timeseries[jj]))
  time0 <- gsub("\\D", "", time0)
  if (is.null(use.spline)) {
    use.spline <- !(length(unique(time0)) == 1 && unique(time0)[1] == "")
  }

  if (use.spline) {
    message("[ngs.fitContrastsWithLIMMA.timeseries]: Limma timeseries with interaction term & spline.")
    time0 <- as.numeric(time0)
  } else {
    message("[ngs.fitContrastsWithLIMMA.timeseries]: Limma timeseries with interaction term.")
    time0 <- as.character(unname(timeseries[jj]))
  }

  y <- factor(y)

  if (use.spline) {
    ## Iterate across df. Pick first valid run.
    idx <- 1:length(unique(time0))
    for (i in 1:length(idx)) {
      ndf <- length(unique(time0)) - idx[i]
      time.spline <- try(splines::ns(time0, df = ndf), silent = TRUE)
      if ("try-error" %in% class(time.spline)) next
      design <- model.matrix(~ y * time.spline)
      fit <- try(limma::lmFit(X, design), silent = TRUE)
      fit <- try(limma::eBayes(fit, trend = trend, robust = TRUE), silent = TRUE)
      if ("try-error" %in% class(fit)) {
        fit <- try(limma::eBayes(fit, trend = trend, robust = FALSE), silent = TRUE)
        if ("try-error" %in% class(fit)) {
          fit <- try(limma::eBayes(fit, trend = FALSE, robust = FALSE), silent = TRUE)
        }
      }
      if ("try-error" %in% class(fit)) next()

      coefs <- apply(fit$t, 2, function(x) sum(is.na(x)))
      est.coefs <- names(coefs[coefs != nrow(fit$t)])
      est.coefs <- est.coefs[grep(":time.spline*", est.coefs)]
      sel <- match(est.coefs, names(coefs))
      top <- try(limma::topTable(fit, coef = sel, sort.by = "none", number = Inf, adjust.method = "BH"),
        silent = TRUE
      )
      if ("try-error" %in% class(top) || nrow(top) == 0) next else break
    }

    if ("try-error" %in% class(top) || nrow(top) == 0) {
      top <- data.frame(matrix(NA, nrow = nrow(X), ncol = 5))
      rownames(top) <- rownames(X)
      colnames(top) <- c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val")
    } else {
      top0 <- data.frame(
        logFC = NA, AveExpr = top$AveExpr, t = NA,
        P.Value = top$P.Value, adj.P.Val = top$adj.P.Val, row.names = rownames(top)
      )
      i <- 1
      SEL <- list()
      for (i in 1:length(est.coefs)) {
        sel <- match(est.coefs[i], names(coefs))
        SEL[[est.coefs[i]]] <- limma::topTable(fit, coef = sel, sort.by = "none", number = Inf, adjust.method = "BH")
      }
      top <- do.call(cbind, SEL)
      index <- c("logFC", "t")
      for (i in 1:length(index)) {
        idx <- paste0(names(SEL), ".", index[i])
        sel <- match(idx, colnames(top))
        top0[, index[i]] <- apply(top[, sel, drop = FALSE], 1, function(x) x[which.max(abs(x))])
      }
      top <- top0
      rm(top0)
      ## kk <- c("logFC", "AveExpr", "F", "P.Value", "adj.P.Val")
      kk <- c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val")
      top <- top[, kk, drop = FALSE]
    }
  } else {
    time0 <- as.factor(time0)
    design <- stats::model.matrix(~ 0 + y + time0 + y:time0)
    fit <- limma::lmFit(X, design)
    fit <- limma::eBayes(fit, trend = trend)
    sel <- grep("*:time0*", colnames(coef(fit)))
    top <- limma::topTable(fit, coef = sel, number = nrow(X))
    hh <- grep(".time0", colnames(top))
    if (any(hh)) {
      logFC <- apply(top[, hh, drop = FALSE], 1, function(x) x[which.max(abs(x))])
      top$logFC <- logFC
      top <- top[, -hh, drop = FALSE]
      kk <- c("logFC", "AveExpr", "F", "P.Value", "adj.P.Val")
      top <- top[, kk, drop = FALSE]
    } else {
      kk <- c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val")
      top <- top[, kk, drop = FALSE]
    }
    colnames(top) <- c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val")
  }

  return(top)
}


#' @describeIn ngs.fitContrastsWithAllMethods Fit contrasts with EdgeR
#' @export
ngs.fitContrastsWithEDGER <- function(counts,
                                      group,
                                      contr.matrix,
                                      design,
                                      covariates = NULL,
                                      method = c("qlf", "lrt"),
                                      prune.samples = FALSE,
                                      X = NULL,
                                      conform.output = FALSE,
                                      robust = TRUE,
                                      calc.tmm = TRUE,
                                      recalc.fc = TRUE,
                                      plot = TRUE,
                                      timeseries = NULL) {
  method <- method[1]

  ## ps: EdgeR/Deseq2 tests will not be available for proteomics data.
  ## Therefore, this autoscaling will very rarely be run.
  counts <- playbase::counts.autoScaling(counts)$counts

  exp0 <- contr.matrix
  if (!is.null(design)) exp0 <- design %*% contr.matrix

  if (prune.samples) {
    kk <- which(rowSums(exp0 != 0, na.rm = TRUE) > 0)
    exp0 <- exp0[kk, , drop = FALSE]
    counts <- counts[, kk, drop = FALSE]
    if (!is.null(design)) design <- design[kk, ]
    if (is.null(design)) contr.matrix <- contr.matrix[kk, , drop = FALSE]
    if (!is.null(group)) group <- group[kk]
  }

  dge <- edgeR::DGEList(round(counts), group = NULL) ## we like integer counts...
  dge$samples$group <- group
  if (calc.tmm) dge <- edgeR::normLibSizes(dge, method = "TMM")

  if (is.null(design) && !prune.samples) { ## This will never run
    message("[ngs.fitContrastsWithEDGER] fitting EDGER contrasts *without* design, no pruning ")
    res <- .ngs.fitContrastsWithEDGER.nodesign(
      dge = dge, contr.matrix = contr.matrix, method = method,
      conform.output = conform.output, robust = robust, calc.tmm = calc.tmm,
      recalc.fc = recalc.fc, plot = plot, timeseries = NULL
    )
    return(res)
  }

  if (is.null(design) && prune.samples) {
    message("[ngs.fitContrastsWithEDGER] fitting EDGER contrasts *without* design, with pruning")
    res <- .ngs.fitContrastsWithEDGER.nodesign.pruned(
      counts = counts, contr.matrix = contr.matrix, method = method, group = group,
      conform.output = conform.output, robust = robust, calc.tmm = calc.tmm,
      recalc.fc = recalc.fc, plot = plot, timeseries = timeseries
    )
    return(res)
  }

  message("[ngs.fitContrastsWithEDGER] fitting EDGER contrasts using design matrix")
  dge <- try(dge <- edgeR::estimateDisp(dge, design = design, robust = robust), silent = TRUE)
  if ("try-error" %in% class(dge)) {
    message("[ngs.fitContrastsWithEDGER] retrying with trend.method 'loess'...")
    dge <- edgeR::estimateDisp(dge, design = design, robust = robust, trend.method = "loess")
  }

  if (is.null(X)) X <- edgeR::cpm(counts, log = TRUE)

  if (method == "qlf") {
    fit <- edgeR::glmQLFit(dge, design, robust = robust)
  } else if (method == "lrt") {
    fit <- edgeR::glmFit(dge, design, robust = robust)
  } else {
    stop("unknown method")
  }

  ## get top table and calculate means
  exp.matrix <- (design %*% contr.matrix)
  tables <- list()
  for (i in 1:ncol(contr.matrix)) {
    cntr <- contr.matrix[, i]
    if (method == "qlf") {
      ct <- edgeR::glmQLFTest(fit, contrast = cntr)
    } else if (method == "lrt") {
      ct <- edgeR::glmLRT(fit, contrast = cntr)
    } else {
      stop("unknown method")
    }
    top <- edgeR::topTags(ct, n = Inf, sort.by = "none")$table
    top <- data.frame(top[rownames(X), ])

    ## calculate means
    j1 <- which(exp.matrix[, i] > 0)
    j2 <- which(exp.matrix[, i] < 0)
    mean1 <- rowMeans(X[, j1, drop = FALSE], na.rm = TRUE)
    mean0 <- rowMeans(X[, j2, drop = FALSE], na.rm = TRUE)
    ## logFC of edgeR is not really reliable..
    if (recalc.fc) top$logFC <- (mean1 - mean0)
    top <- cbind(top, "AveExpr0" = mean0, "AveExpr1" = mean1)
    tables[[i]] <- top
  }
  names(tables) <- colnames(contr.matrix)

  if (conform.output == TRUE) {
    for (i in 1:length(tables)) {
      if (method == "qlf") {
        k1 <- c("logFC", "logCPM", "F", "PValue", "FDR", "AveExpr0", "AveExpr1")
      } else if (method == "lrt") {
        k1 <- c("logFC", "logCPM", "LR", "PValue", "FDR", "AveExpr0", "AveExpr1")
      } else {
        stop("switch method error")
      }
      k2 <- c("logFC", "AveExpr", "statistic", "P.Value", "adj.P.Val", "AveExpr0", "AveExpr1")
      tables[[i]] <- tables[[i]][, k1]
      colnames(tables[[i]]) <- k2
    }
  }
  res <- list(tables = tables)
  return(res)
}


#' @describeIn ngs.fitContrastsWithAllMethods contrasts with EdgeR without a design
#' @export
.ngs.fitContrastsWithEDGER.nodesign <- function(dge,
                                                contr.matrix,
                                                method = c("qlf", "lrt"),
                                                X = NULL,
                                                conform.output = FALSE,
                                                robust = TRUE,
                                                calc.tmm = TRUE,
                                                recalc.fc = TRUE,
                                                plot = TRUE,
                                                timeseries = NULL) {

  ## With no design matrix, we must do EdgeR per contrast
  ## one-by-one. Warning this can become very slow.

  if (class(dge) != "DGEList") stop("dge must be a DGEList object")
  method <- method[1]

  if (is.null(X)) X <- edgeR::cpm(dge$counts, log = TRUE)
  dge1 <- dge
  if (calc.tmm) dge1 <- edgeR::normLibSizes(dge1, method = "TMM")
  dge1 <- edgeR::estimateDisp(dge1, design = NULL, robust = robust)

  contr0 <- matrix(c(-1, 0, 1), nrow = 3)
  rownames(contr0) <- c("yneg", "yo", "ypos")
  colnames(contr0) <- "pos_vs_neg"

  tables <- list()
  for (i in 1:ncol(contr.matrix)) {
    ct <- contr.matrix[, i]
    y <- factor(c("neg", "o", "pos")[2 + sign(ct)])

    if (grepl("^IA:*", colnames(contr.matrix)[i]) && !is.null(timeseries)) {
      top <- .ngs.fitContrastsWithEDGER.nodesign.timeseries(
        dge = dge,
        counts = as.matrix(dge$counts),
        X = X,
        y = y,
        method = method,
        timeseries = timeseries,
        robust = robust
      )
    } else {
      design1 <- stats::model.matrix(~ 0 + y)
      if (method == "qlf") {
        fit <- edgeR::glmQLFit(dge1, design1, robust = robust)
        ctx <- contr0[colnames(stats::coef(fit)), ]
        res <- edgeR::glmQLFTest(fit, contrast = ctx)
      } else if (method == "lrt") {
        fit <- edgeR::glmFit(dge1, design1, robust = robust)
        ctx <- contr0[colnames(stats::coef(fit)), ]
        res <- edgeR::glmLRT(fit, contrast = ctx)
      } else {
        stop("unknown method: ", method)
      }
      top <- edgeR::topTags(res, n = 1e9)$table
      top <- data.frame(top[rownames(X), ])
    }

    j1 <- which(contr.matrix[, i] > 0)
    j0 <- which(contr.matrix[, i] < 0)
    mean1 <- rowMeans(X[, j1, drop = FALSE], na.rm = TRUE)
    mean0 <- rowMeans(X[, j0, drop = FALSE], na.rm = TRUE)
    ## logFC of edgeR is not really reliable..
    if (recalc.fc) top$logFC <- (mean1 - mean0)
    top <- cbind(top, "AveExpr0" = mean0, "AveExpr1" = mean1)
    tables[[i]] <- top
    names(tables)[i] <- colnames(contr.matrix)[i]
  }

  if (conform.output == TRUE) {
    for (i in 1:length(tables)) {
      if (method == "qlf") {
        k1 <- c("logFC", "logCPM", "F", "PValue", "FDR", "AveExpr0", "AveExpr1")
      } else if (method == "lrt") {
        k1 <- c("logFC", "logCPM", "LR", "PValue", "FDR", "AveExpr0", "AveExpr1")
      } else {
        stop("switch method error")
      }
      k2 <- c("logFC", "AveExpr", "statistic", "P.Value", "adj.P.Val", "AveExpr0", "AveExpr1")
      tables[[i]] <- tables[[i]][, k1]
      colnames(tables[[i]]) <- k2
    }
  }

  res <- list(tables = tables)
  return(res)

}


#' @describeIn ngs.fitContrastsWithAllMethods Fit contrasts with EdgeR without a design matrix
#' @export
.ngs.fitContrastsWithEDGER.nodesign.pruned <- function(counts,
                                                       contr.matrix,
                                                       group = NULL,
                                                       method = c("qlf", "lrt"),
                                                       X = NULL,
                                                       conform.output = FALSE,
                                                       robust = TRUE,
                                                       calc.tmm = TRUE,
                                                       recalc.fc = TRUE,
                                                       plot = TRUE,
                                                       timeseries = NULL) {
  method <- method[1]

  ## ps: EdgeR/Deseq2 tests will not be available for proteomics data.
  ## Therefore, this autoscaling will very rarely be run.
  counts <- playbase::counts.autoScaling(counts)$counts

  tables <- list()
  for (i in 1:NCOL(contr.matrix)) {
    kk <- which(!is.na(contr.matrix[, i]) & contr.matrix[, i] != 0)
    counts1 <- counts[, kk, drop = FALSE]
    # X1 <- NULL
    if (!is.null(X)) {
      X1 <- X[, kk, drop = FALSE]
    } else {
      X1 <- edgeR::cpm(counts1, log = TRUE)
    }

    group1 <- group
    if (!is.null(group)) group1 <- group[kk]
    ct <- contr.matrix[kk, i]
    y <- factor(c("neg", "o", "pos")[2 + sign(ct)])
    dge1 <- edgeR::DGEList(round(counts1), group = group1)
    if (calc.tmm) dge1 <- edgeR::normLibSizes(dge1, method = "TMM")

    if (grepl("^IA:*", colnames(contr.matrix)[i]) && !is.null(timeseries)) {
      top <- .ngs.fitContrastsWithEDGER.nodesign.timeseries(
        dge = dge1,
        counts = counts1,
        X = X1,
        y = y,
        method = method,
        timeseries = timeseries,
        robust = robust
      )
    } else {
      design1 <- stats::model.matrix(~y)
      dge1 <- edgeR::estimateDisp(dge1, design = design1, robust = robust)
      if (method == "qlf") {
        fit <- edgeR::glmQLFit(dge1, design1, robust = robust)
        res <- edgeR::glmQLFTest(fit, coef = 2)
      } else if (method == "lrt") {
        fit <- edgeR::glmFit(dge1, design1, robust = robust)
        res <- edgeR::glmLRT(fit, coef = 2)
      } else {
        stop("unknown method: ", method)
      }
      top <- edgeR::topTags(res, n = 1e9)$table
      top <- data.frame(top[rownames(X1), ])
    }

    contr1 <- contr.matrix[kk, i]
    j1 <- which(contr1 > 0)
    j0 <- which(contr1 < 0)
    mean1 <- rowMeans(X1[, j1, drop = FALSE], na.rm = TRUE)
    mean0 <- rowMeans(X1[, j0, drop = FALSE], na.rm = TRUE)
    ## logFC of edgeR is not really reliable..
    if (recalc.fc) top$logFC <- (mean1 - mean0)
    top <- cbind(top, "AveExpr0" = mean0, "AveExpr1" = mean1)
    Matrix::head(top)
    tables[[i]] <- top
    names(tables)[i] <- colnames(contr.matrix)[i]
  }

  if (conform.output == TRUE) {
    for (i in 1:length(tables)) {
      if (method %in% c("qlf", "qlf2")) {
        k1 <- c("logFC", "logCPM", "F", "PValue", "FDR", "AveExpr0", "AveExpr1")
      } else if (method %in% c("lrt", "lrt2")) {
        k1 <- c("logFC", "logCPM", "LR", "PValue", "FDR", "AveExpr0", "AveExpr1")
      } else {
        stop("switch method error")
      }
      k2 <- c("logFC", "AveExpr", "statistic", "P.Value", "adj.P.Val", "AveExpr0", "AveExpr1")
      tables[[i]] <- tables[[i]][, k1]
      colnames(tables[[i]]) <- k2
    }
  }

  res <- list(tables = tables)
  return(res)
}

## #' @describeIn ngs.fitContrastsWithAllMethods Fits time-series contrasts using edgeR QLF or LRT
## #' @export
.ngs.fitContrastsWithEDGER.nodesign.timeseries <- function(dge,
                                                           counts,
                                                           X,
                                                           y,
                                                           method,
                                                           timeseries,
                                                           use.spline = NULL,
                                                           robust = TRUE) {

  if (!all(colnames(counts) %in% names(timeseries))) {
    message("[ngs.fitConstrastsWithEDGER.nodesign.timeseries] Counts and timeseries vector contain different set of samples.")
  }

  if (!method %in% c("lrt", "qlf")) {
    stop("[ngs.fitConstrastsWithEDGER.nodesign.timeseries] EdgeR test unrecognized. Must be LRT or QLF.")
  }

  jj <- match(colnames(counts), names(timeseries))
  time0 <- as.character(timeseries[jj])
  time0 <- gsub("\\D", "", unname(time0))

  if (is.null(use.spline)) {
    use.spline <- !(length(unique(time0)) == 1 && unique(time0)[1] == "")
  }

  ## possible ranges
  if (any(grepl("_", unname(as.character(timeseries[jj]))))) {
    use.spline <- FALSE 
  }

  if (use.spline) {
    message("[ngs.fitConstrastsWithEDGER.nodesign.timeseries]: EdgeR timeseries with interaction & spline.")
    time0 <- as.numeric(time0)
  } else {
    message("[ngs.fitConstrastsWithEDGER.nodesign.timeseries]: EdgeR timeseries with interaction.")
    time0 <- as.character(timeseries[jj])
  }

  y <- as.factor(as.character(y))
  dge0 <- dge
  res <- NULL
  
  if (use.spline) {
    ## Iterate across df. Pick first valid run.
    idx <- 1:length(unique(time0))
    for (i in 1:length(idx)) {
      ndf <- length(unique(time0)) - idx[i]
      design <- try(stats::model.matrix(~ 0 + y * splines::ns(time0, df = ndf)), silent = TRUE)
      if ("try-error" %in% class(design)) next
      message("[ngs.fitConstrastsWithEDGER.nodesign.timeseries] Using splines with ", ndf, " degrees of freedom.")
      dge <- try(edgeR::estimateDisp(dge, design = design, robust = robust), silent = TRUE)
      if ("try-error" %in% class(dge)) next
      sel <- grep("*:splines::ns*", colnames(design))
      if (method == "qlf") {
        fit <- try(edgeR::glmQLFit(dge, design, robust = robust), silent = TRUE)
        # colnames(stats::coef(fit))[sel]
        if ("try-error" %in% class(fit)) next
        res <- try(edgeR::glmQLFTest(fit, coef = sel), silent = TRUE)
        if ("try-error" %in% class(res)) next else break
      } else if (method == "lrt") {
        fit <- try(edgeR::glmFit(dge, design, robust = robust), silent = TRUE)
        if ("try-error" %in% class(fit)) next
        res <- try(edgeR::glmLRT(fit, coef = sel), silent = TRUE)
        if ("try-error" %in% class(res)) next else break
      }
    }
  }

  if (!use.spline | is.null(res)) {
    design <- model.matrix(~ y * time0)
    dge0 <- try(edgeR::estimateDisp(dge0, design = design, robust = robust), silent = TRUE)
    # Test interaction terms directly
    sel <- grep("*:time0*", colnames(design)) ## ???? unclear.
    if (method == "qlf") {
      fit <- try(edgeR::glmQLFit(dge0, design, robust = robust), silent = TRUE)
      res <- try(edgeR::glmQLFTest(fit, coef = sel), silent = TRUE)
    } else if (method == "lrt") {
      fit <- try(edgeR::glmFit(dge0, design, robust = robust), silent = TRUE)
      res <- try(edgeR::glmLRT(fit, coef = sel), silent = TRUE)
    }
  }

  if ("try-error" %in% class(res)) {
    top <- data.frame(matrix(NA, nrow = nrow(X), ncol = 5))
    rownames(top) <- rownames(X)
    kk <- c("logFC", "logCPM", "LR", "PValue", "FDR")
    if (method == "qlf") kk[which(kk == "LR")] <- "F"
    colnames(top) <- kk
  } else {
    top <- edgeR::topTags(res, n = 1e9)$table
    top <- data.frame(top[rownames(X), ])
    sel <- grep("logFC.*", colnames(top))
    logFC <- as.numeric(apply(top[, sel, drop = FALSE], 1, function(x) x[which.max(abs(x))]))
    top <- cbind(logFC = logFC, top[, -sel, drop = FALSE])
  }

  return(top)

}

#' @describeIn ngs.fitContrastsWithAllMethods Fits contrasts using DESeq2 differential expression
#' analysis on count data
#' @export
ngs.fitContrastsWithDESEQ2 <- function(counts,
                                       group,
                                       contr.matrix,
                                       design,
                                       covariates = NULL,
                                       X = NULL,
                                       genes = NULL,
                                       test = "Wald",
                                       fitType = "mean",
                                       prune.samples = FALSE,
                                       conform.output = FALSE,
                                       timeseries = NULL) {
  ## ps: EdgeR/Deseq2 tests will not be available for proteomics data.
  ## Therefore, this autoscaling will very rarely be run.
  counts <- playbase::counts.autoScaling(counts)$counts

  exp0 <- contr.matrix
  if (!is.null(design)) exp0 <- design %*% contr.matrix

  if (prune.samples) {
    kk <- which(rowSums(exp0 != 0, na.rm = TRUE) > 0)
    exp0 <- exp0[kk, , drop = FALSE]
    counts <- counts[, kk, drop = FALSE]
    if (!is.null(design)) design <- design[kk, ]
    if (is.null(design)) contr.matrix <- contr.matrix[kk, , drop = FALSE]
    if (!is.null(group)) group <- group[kk]
  }

  if (is.null(design)) {
    message("[ngs.fitContrastsWithDESEQ2] fitting DESEQ2  *without* design")
    out <- .ngs.fitContrastsWithDESEQ2.nodesign(
      counts = counts,
      contr.matrix = contr.matrix,
      test = test,
      timeseries = timeseries,
      prune.samples = prune.samples,
      conform.output = conform.output
    )
    return(out)
  }
  message("[ngs.fitContrastsWithDESEQ2] fitting DESEQ2 using design matrix")

  design.formula <- stats::formula(" ~ 0 + group")
  message("[ngs.fitContrastsWithDESEQ2] using model design: ", as.character(design.formula))

  rownames.counts <- rownames(counts)
  counts <- round(counts) ## WARNING!!!

  if (all(rowSums(counts == 0, na.rm = TRUE) > 0)) {
    ## ERROR: 'every gene contains at least one zero, cannot compute log
    ## geometric means' so we fix it villager-style
    jmax <- which.max(rowSums(counts, na.rm = TRUE))
    counts[jmax, ] <- pmax(counts[jmax, ], 1)
  }

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    design = design.formula,
    colData = data.frame(group)
  )

  rownames(counts) <- rownames.counts

  ## Run DESeq : Modeling counts with generic 'group'
  # fitTypes <- c("parametric", "local", "mean", "glmGamPoi")
  if (test == "glmGamPoi" || fitType == "glmGamPoi") {
    dds <- try(DESeq2::DESeq(dds, fitType = "glmGamPoi", test = "LRT", reduced = ~1))
  } else if (test == "LRT") {
    dds <- try(DESeq2::DESeq(dds, fitType = fitType, test = "LRT", reduced = ~1))
  } else {
    dds <- try(DESeq2::DESeq(dds, fitType = fitType, test = "Wald"))
  }

  ## Sometime DESEQ2 fails and we resort to gene-wise estimates. (IK:
  ## what is this???)
  if ("try-error" %in% class(dds)) {
    message("[.ngs.fitContrastsWithDESEQ2] retrying DESEQ2 with gene-wise estimates...")
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = counts,
      design = design.formula,
      colData = data.frame(group)
    )
    dds <- DESeq2::estimateSizeFactors(dds)
    disp.type <- ifelse(test == "glmGamPoi", "glmGamPoi", "DESeq2")
    dds <- DESeq2::estimateDispersionsGeneEst(dds, type = disp.type)
    DESeq2::dispersions(dds) <- GenomicRanges::mcols(dds)$dispGeneEst
    suppressWarnings({
      if (test == "LRT") {
        dds <- try(DESeq2::nbinomLRT(dds), silent = TRUE)
      } else if (test == "glmGamPoi") {
        dds <- try(DESeq2::nbinomLRT(dds, type = "glmGamPoi"), silent = TRUE)
      } else {
        dds <- try(DESeq2::nbinomWaldTest(dds), silent = TRUE)
      }
    })
  }

  ## we add the gene annotation here (not standard...)
  if (!is.null(genes)) SummarizedExperiment::rowData(dds)$genes <- genes ## does this work??

  ## logCPM for calculating means
  if (is.null(X)) X <- edgeR::cpm(counts, log = TRUE)
  exp.matrix <- contr.matrix
  if (!is.null(design)) exp.matrix <- (design %*% contr.matrix)

  tables <- list()
  for (i in 1:ncol(contr.matrix)) {
    ## manually specify contrast vector. See also https://support.bioconductor.org/p/69104/
    contr <- contr.matrix[, i]
    contr[is.na(contr)] <- 0
    contr <- (contr > 0) / sum(contr > 0) - (contr < 0) / sum(contr < 0) ## mean zero, signed sum to one.
    DESeq2::resultsNames(dds)
    if (any(grepl("group", DESeq2::resultsNames(dds)))) {
      grp.contr <- contr
      names(grp.contr) <- paste0("group", names(contr))
      contr <- rep(0, length(DESeq2::resultsNames(dds)))
      names(contr) <- DESeq2::resultsNames(dds)
      contr[names(grp.contr)] <- grp.contr
    }
    ## do no set p values to NA
    resx <- DESeq2::results(dds, contrast = contr, cooksCutoff = FALSE, independentFiltering = FALSE)
    pos.samples <- which(exp.matrix[, i] > 0)
    neg.samples <- which(exp.matrix[, i] < 0)
    resx$AveExpr1 <- rowMeans(X[, pos.samples, drop = FALSE], na.rm = TRUE)
    resx$AveExpr0 <- rowMeans(X[, neg.samples, drop = FALSE], na.rm = TRUE)
    resx$log2BaseMean <- log2(0.0001 + resx$baseMean)
    if (conform.output) resx$log2FoldChange <- (resx$AveExpr1 - resx$AveExpr0) ## recompute
    tables[[i]] <- data.frame(resx)
    names(tables)[i] <- colnames(contr.matrix)[i]
  }
  names(tables) <- colnames(contr.matrix)

  if (conform.output == TRUE) {
    for (i in 1:length(tables)) {
      k1 <- c("log2FoldChange", "log2BaseMean", "stat", "pvalue", "padj", "AveExpr0", "AveExpr1")
      k2 <- c("logFC", "AveExpr", "statistic", "P.Value", "adj.P.Val", "AveExpr0", "AveExpr1")
      tables[[i]] <- tables[[i]][, k1]
      colnames(tables[[i]]) <- k2
    }
  }

  res <- list(tables = tables)
  return(res)
}


#' @describeIn ngs.fitContrastsWithAllMethods Fits contrasts using DESeq2 differential expression
#' @export
.ngs.fitContrastsWithDESEQ2.nodesign <- function(counts,
                                                 contr.matrix,
                                                 test = "Wald",
                                                 fitType = "mean",
                                                 prune.samples = FALSE,
                                                 conform.output = FALSE,
                                                 X = NULL,
                                                 timeseries = NULL) {
  ## ps: EdgeR/Deseq2 tests will not be available for proteomics data.
  ## Therefore, this autoscaling will very rarely be run.
  counts <- playbase::counts.autoScaling(counts)$counts

  counts <- round(counts)
  if (is.null(X)) X <- edgeR::cpm(counts, log = TRUE)

  if (nrow(contr.matrix) != ncol(X)) {
    stop("ngs.fitContrastsWithDESEQ2.nodesign:: contrast matrix must be by sample")
  }

  exp.matrix <- contr.matrix
  tables <- list()
  for (i in 1:ncol(exp.matrix)) {
    ## manual design matrix (CHECK THIS!!!)
    kk <- 1:nrow(exp.matrix)
    if (prune.samples) {
      kk <- which(!is.na(exp.matrix[, i]) & exp.matrix[, i] != 0)
    }
    ct <- exp.matrix[kk, i]
    y <- factor(c("neg", "zero", "pos")[2 + sign(ct)], levels = c("neg", "zero", "pos"))
    counts1 <- counts[, kk, drop = FALSE]
    if (grepl("^IA:*", colnames(exp.matrix)[i]) && !is.null(timeseries)) {
      resx <- .ngs.fitContrastsWithDESEQ2.nodesign.timeseries(counts1, y, timeseries, test = test)
    } else {
      colData <- data.frame(y, row.names = colnames(counts1))
      ## sample-wise model matrix (does this work???)
      colnames(counts1) <- NULL
      design.formula <- stats::formula("~ 0 + y")
      dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = counts1,
        design = design.formula,
        colData = colData
      )
      ## fitType <- "mean"
      suppressWarnings({
        if (test == "glmGamPoi" || fitType == "glmGamPoi") {
          dds <- try(DESeq2::DESeq(dds, fitType = "glmGamPoi", test = "LRT", reduced = ~1))
        } else if (test == "LRT") {
          dds <- try(DESeq2::DESeq(dds, fitType = fitType, test = "LRT", reduced = ~1))
        } else {
          dds <- try(DESeq2::DESeq(dds, fitType = fitType, test = "Wald"))
        }
      })
      ## sometime DESEQ2 fails
      if ("try-error" %in% class(dds)) {
        message("[.ngs.fitContrastsWithDESEQ2.nodesign] retrying DESEQ2 with gene-wise estimates...")
        dds <- DESeq2::DESeqDataSetFromMatrix(
          countData = counts1,
          design = design.formula,
          colData = data.frame(y)
        )
        dds <- DESeq2::estimateSizeFactors(dds)
        disp.type <- ifelse(test == "glmGamPoi", "glmGamPoi", "DESeq2")
        dds <- DESeq2::estimateDispersionsGeneEst(dds, type = disp.type)
        DESeq2::dispersions(dds) <- GenomicRanges::mcols(dds)$dispGeneEst
        suppressWarnings({
          if (test == "LRT") {
            dds <- try(DESeq2::nbinomLRT(dds))
          } else if (test == "glmGamPoi") {
            dds <- try(DESeq2::nbinomLRT(dds, type = "glmGamPoi"))
          } else {
            dds <- try(DESeq2::nbinomWaldTest(dds))
          }
        })
      }
      ## DESeq2::resultsNames(dds)
      ctx <- c("yneg" = -1, "yzero" = 0, "ypos" = 1)[DESeq2::resultsNames(dds)]
      resx <- DESeq2::results(dds, contrast = ctx, cooksCutoff = FALSE, independentFiltering = FALSE)
      ## we add the gene annotation here (not standard...)
      rownames(resx) <- rownames(SummarizedExperiment::rowData(dds))
    }
    X1 <- X[, kk, drop = FALSE]
    pos.samples <- which(exp.matrix[kk, i] > 0)
    neg.samples <- which(exp.matrix[kk, i] < 0)
    resx$AveExpr1 <- rowMeans(X1[, pos.samples, drop = FALSE], na.rm = TRUE)
    resx$AveExpr0 <- rowMeans(X1[, neg.samples, drop = FALSE], na.rm = TRUE)
    resx$log2BaseMean <- log2(0.0001 + resx$baseMean)
    if (conform.output) {
      ## For time-series we keep DESeq2 log2FC from LRT/Wald.
      if (!grepl("^IA:*", colnames(exp.matrix)[i])) {
        resx$log2FoldChange <- (resx$AveExpr1 - resx$AveExpr0) ## recompute
      }
    }
    tables[[i]] <- data.frame(resx)
    names(tables)[i] <- colnames(exp.matrix)[i]
  }

  if (conform.output) {
    for (i in 1:length(tables)) {
      k1 <- c("log2FoldChange", "log2BaseMean", "stat", "pvalue", "padj", "AveExpr0", "AveExpr1")
      k2 <- c("logFC", "AveExpr", "statistic", "P.Value", "adj.P.Val", "AveExpr0", "AveExpr1")
      tables[[i]] <- tables[[i]][, k1]
      colnames(tables[[i]]) <- k2
    }
  }

  return(list(tables = tables))
}

## Q: does the condition induces a change in gene expression at
## any time point after the reference level time point (time 0)?
## https://support.bioconductor.org/p/62684/
## Features showing a consistent difference from time 0 onward will not have a small p value.
## This is imporant because differences between groups at time 0 should be controlled for,
## e.g. random differences between the individuals chosen for each treatment group which
## are observable before the treatment takes effect.
# design.formula <- stats::formula("~ y + time + y:time")
# dds <- DESeq2::DESeqDataSetFromMatrix(counts, colData, design)
# dds <- DESeq2::DESeq(dds, test = "LRT", reduced = ~ y + time)
# resx <- DESeq2::results(dds, cooksCutoff = FALSE, independentFiltering = FALSE)

#' @describeIn ngs.fitContrastsWithAllMethods Fits time-series contrasts using DESeq2 LRT
#' @export
.ngs.fitContrastsWithDESEQ2.nodesign.timeseries <- function(counts,
                                                            y,
                                                            timeseries,
                                                            test = "LRT",
                                                            fitType = "mean",
                                                            use.spline = NULL) {

  ## ps: EdgeR/Deseq2 tests will not be available for proteomics data.
  ## Therefore, this autoscaling will very rarely be run.
  counts <- playbase::counts.autoScaling(counts)$counts

  if (!all(colnames(counts) %in% names(timeseries))) {
    stop("[ngs.fitConstrastsWithDESEQ2.nodesign.timeseries] counts and time contain different set of samples.")
  }

  if (!test %in% c("LRT", "Wald", "glmGamPoi")) {
    stop("[.ngs.fitConstrastsWithDESEQ2.nodesign.timeseries] DESeq test unrecognized. Must be LRT or Wald.")
  }

  jj <- match(colnames(counts), names(timeseries))
  time0 <- as.character(timeseries[jj])
  time0 <- gsub("\\D", "", unname(time0))
 
  if (is.null(use.spline)) {
    use.spline <- !(length(unique(time0)) == 1 && unique(time0)[1] == "")
  }

  ## possible ranges
  if (any(grepl("_", unname(as.character(timeseries[jj]))))) {
    use.spline <- FALSE 
  }

  if (use.spline) {
    time0 <- as.numeric(time0)
    # !!: splines::ns need data in range c(0,60). Else breaks
    # !![AZ]. IK: this is weird.
    if (max(range(time0, na.rm = TRUE)) > 60) time0 <- time0 / 60
    message("[ngs.fitConstrastsWithDESEQ2.nodesign.timeseries]: DESeq2 timeseries with interaction term & spline.")
  } else {
    time0 <- as.character(timeseries[jj])
    message("[ngs.fitConstrastsWithDESEQ2.nodesign.timeseries]: DESeq2 timeseries with interaction term. No spline.")
  }

  y <- as.factor(as.character(y))
  colData <- data.frame(y = y, time = factor(time0))

  if (use.spline) {

    ## Iterate across df. Pick first valid run.
    idx <- 1:length(unique(time0))
    for (i in 1:length(idx)) {
      ndf <- length(unique(time0)) - idx[i]
      sp <- splines::ns(time0, df = ndf)
      if ("try-error" %in% class(sp)) next
      design <- try(stats::model.matrix(~ y + sp + y:sp), silent = TRUE)
      if ("try-error" %in% class(design)) next
      red.design <- stats::model.matrix(~ y + sp)

      dds <- try(DESeq2::DESeqDataSetFromMatrix(counts, colData, design), silent = TRUE)
      if ("try-error" %in% class(dds)) next

      if (test == "glmGamPoi" || fitType == "glmGamPoi") {
        dds <- try(DESeq2::DESeq(dds, fitType = "glmGamPoi", test = "LRT", reduced = red.design), silent = TRUE)
      } else if (test == "LRT") {
        dds <- try(DESeq2::DESeq(dds, fitType = fitType, test = "LRT", reduced = red.design), silent = TRUE)
      } else if (test == "Wald") {
        dds <- try(DESeq2::DESeq(dds, fitType = fitType, test = "Wald"), silent = TRUE)
      } else {
        stop("[.ngs.fitContrastsWithDESEQ2.nodesign.timeseries] DESeq2::DESeq test unrecognized")
      }

      if ("try-error" %in% class(dds)) {
        dds <- DESeq2::DESeqDataSetFromMatrix(counts, colData, design)
        dds <- DESeq2::estimateSizeFactors(dds)
        disp.type <- ifelse(test == "glmGamPoi", "glmGamPoi", "DESeq2")
        dds <- try(DESeq2::estimateDispersionsGeneEst(dds, type = disp.type), silent = TRUE)
        if ("try-error" %in% class(dds)) next
        DESeq2::dispersions(dds) <- GenomicRanges::mcols(dds)$dispGeneEst
        if (test == "LRT") {
          dds <- try(DESeq2::nbinomLRT(dds, reduced = red.design), silent = TRUE)
        } else if (test == "glmGamPoi") {
          dds <- try(DESeq2::nbinomLRT(dds, type = "glmGamPoi"))
        } else {
          dds <- try(DESeq2::nbinomWaldTest(dds), silent = TRUE)
        }
        if ("try-error" %in% class(dds)) next
      }
      message("[ngs.fitConstrastsWithDESEQ2.nodesign.timeseries] Using splines with ", ndf, " degrees of freedom.")
      break ## stop for loop
    }

  } else {

    design <- stats::model.matrix(~ y + time0 + y:time0)
    red.design <- stats::model.matrix(~ y + time0)
    dds <- try(DESeq2::DESeqDataSetFromMatrix(counts, colData, design), silent = TRUE)
    if (test == "glmGamPoi" || fitType == "glmGamPoi") {
      dds <- try(DESeq2::DESeq(dds, fitType = "glmGamPoi", test = "LRT", reduced = red.design), silent = TRUE)
    } else if (test == "LRT") {
      dds <- try(DESeq2::DESeq(dds, test = "LRT", fitType = "mean", reduced = red.design), silent = TRUE)
    } else if (test == "Wald") {
      dds <- try(DESeq2::DESeq(dds, test = "Wald", fitType = "mean"), silent = TRUE)
    } else {
      stop("[.ngs.fitContrastsWithDESEQ2.nodesign.timeseries] DESeq2::DESeq test unrecognized")
    }
    if ("try-error" %in% class(dds)) {
      dds <- try(DESeq2::DESeqDataSetFromMatrix(counts, colData, design), silent = TRUE)
      dds <- try(DESeq2::estimateSizeFactors(dds), silent = TRUE)
      disp.type <- ifelse(test == "glmGamPoi", "glmGamPoi", "DESeq2")
      dds <- try(DESeq2::estimateDispersionsGeneEst(dds, type = disp.type), silent = TRUE)
      disp.values <- try(GenomicRanges::mcols(dds), silent = TRUE)
      c1 <- (! "try-error" %in% class(dds))
      c2 <- (! "try-error" %in% class(disp.values))
      if (c1 & c2) DESeq2::dispersions(dds) <- disp.values$dispGeneEst
      if (test == "LRT") {
        dds <- try(DESeq2::nbinomLRT(dds, reduced = red.design), silent = TRUE)
      } else if (test == "glmGamPoi") {
        dds <- try(DESeq2::nbinomLRT(dds, type = "glmGamPoi"), silent = TRUE)
      } else {
        dds <- try(DESeq2::nbinomWaldTest(dds), silent = TRUE)
      }
    }

  }

  resx <- try(DESeq2::results(dds, cooksCutoff = FALSE, independentFiltering = FALSE), silent = TRUE)
  if ("try-error" %in% class(resx)) {
    resx <- data.frame(matrix(NA, nrow = nrow(counts), ncol = 6))
    rownames(resx) <- rownames(counts)
    colnames(resx) <- c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
  } else {
    rownames(resx) <- rownames(SummarizedExperiment::rowData(dds))
  }

  return(resx)

}

## =====================================================================================
## =========================== END OF FILE =============================================
## =====================================================================================
