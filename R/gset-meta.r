##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

#' Fit Contrasts with All Gene Set Enrichment Methods
#'
#' This function fits contrasts using multiple gene set enrichment methods, such as fisher, ssgsea, gsva, spearman, camera, fry, gsea.permPH, gsea.permGS, gseaPR, and fgsea.
#'
#' @param gmt The gene set matrix.
#' @param X The gene expression matrix.
#' @param Y The phenotype data matrix.
#' @param G The gene annotation matrix.
#' @param design The experimental design matrix.
#' @param contr.matrix The contrast matrix.
#' @param methods A character vector specifying the gene set enrichment methods to use.
#' @param mc.threads The number of threads to use for parallel computing.
#' @param mc.cores The number of CPU cores to use for parallel computing.
#' @param batch.correct Logical indicating whether to correct for batch effects.
#'
#' @return A list containing the results of the gene set enrichment analyses for each method.
#'
#' @export
gset.fitContrastsWithAllMethods <- function(gmt,
                                            X,
                                            Y,
                                            G,
                                            design,
                                            contr.matrix,
                                            methods,
                                            mc.threads = 1,
                                            mc.cores = NULL,
                                            batch.correct = TRUE) {
  ALL.GENESET.METHODS <- c(
    "fisher", "ssgsea", "gsva", "spearman", "camera", "fry",
    "gsea.permPH", "gsea.permGS", "gseaPR", "fgsea"
  )

  timings <- c()

  if (is.null(mc.cores)) {
    mc.cores <- round(0.25 * parallel::detectCores(all.tests = TRUE, logical = FALSE)) ## max 25% of cores
    mc.cores <- pmax(mc.cores, 1) ## min 1 core
    mc.cores <- pmin(mc.cores, 16) ## max 16 cores
  }
  message("using ", mc.cores, " number of cores")
  message("using ", mc.threads, " number of threads")

  if (methods[1] == "*") {
    methods <- ALL.GENESET.METHODS
  }
  methods <- intersect(methods, ALL.GENESET.METHODS)
  message("calculating methods:", paste(methods,collapse=" "))

  ## If degenerate set design to NULL
  if (!is.null(design) && ncol(design) >= ncol(X)) {
    ## "no-replicate" design!!!
    message("WARNING: degenerate design. setting design to NULL")
    contr.matrix <- design %*% contr.matrix
    design <- NULL
  }

  ## experiment matrix
  if (!is.null(design)) {
    exp.matrix <- (design %*% contr.matrix)[colnames(X), , drop = FALSE]
  } else {
    exp.matrix <- contr.matrix[colnames(X), , drop = FALSE]
  }
  Y <- Y[colnames(X), , drop = FALSE]

  ## some "normalization" for single-sample methods
  my.normalize <- function(zx) {
    if (nrow(zx) <= 10) {
      return(zx)
    }
    zx <- scale(limma::normalizeQuantiles(zx))
    return(zx)
  }

  all.results <- list()
  ## pre-compute matrices
  zx.gsva  <- zx.ssgsea <- zx.rnkcorr <- NULL
  res.gsva <- res.ssgsea <- res.rnkcorr <- NULL

  G <- G[rownames(X), names(gmt), drop = FALSE]

  if ("spearman" %in% methods) {
    message("fitting contrasts using spearman/limma... ")

    ## single-sample gene set enrichment using (fast) rank correlation
    xx1 <- X - rowMeans(X, na.rm = TRUE) ## center it...
    xx1 <- apply(xx1, 2, rank, na.last = "keep") ## rank correlation (like spearman)
    jj <- intersect(rownames(G), rownames(xx1))
    tt <- system.time({
      zx.rnkcorr <- qlcMatrix::corSparse(G[jj, , drop = FALSE], xx1[jj, ]) ## superfast
      rownames(zx.rnkcorr) <- colnames(G)
      colnames(zx.rnkcorr) <- colnames(X)

      ## row-wise (per feature) scaling is 'good practice', see
      ## tests comparing rankcor and ssGSEA/gsva

      ## additional batch correction and NNM???
      zx.rnkcorr <- my.normalize(zx.rnkcorr)
      zx.rnkcorr <- zx.rnkcorr[names(gmt), colnames(X), drop = FALSE] ## make sure..

      ## compute LIMMA
      all.results[["spearman"]] <- playbase::gset.fitContrastsWithLIMMA(
        zx.rnkcorr,
        contr.matrix,
        design = design,
        trend = TRUE,
        conform.output = TRUE
      )
    })
    timings <- rbind(timings, c("spearman", tt))
  }

  if ("gsva" %in% methods) {
    message("fitting contrasts using GSVA/limma... ")

    ## check if we have the new version of GSVA
    new.gsva <- exists('gsvaParam', where=asNamespace('GSVA'), mode='function')
    new.gsva    
    tt <- system.time({
      zx.gsva <- NULL
      if(new.gsva) {
        zx.gsva <- try({
          bpparam <- BiocParallel::MulticoreParam(mc.cores)
          GSVA::gsva(GSVA::gsvaParam(as.matrix(X), gmt), BPPARAM = bpparam)
        })
      } else {
        zx.gsva <- try({
          GSVA::gsva(as.matrix(X), gmt, method = "gsva", parallel.sz = mc.cores, verbose = FALSE)
        })
      }

      if (is.null(zx.gsva) || "try-error" %in% class(zx.gsva)) {
        ## switch to single core...
        warning("WARNING: GSVA ERROR: retrying single core... ")
        if(new.gsva) {
          zx.gsva <- try(GSVA::gsva(GSVA::gsvaParam(as.matrix(X), gmt)))
        } else {
          zx.gsva <- try(GSVA::gsva(as.matrix(X), gmt, method = "gsva", parallel.sz = 1, verbose = FALSE ))
        }
      }

      if (!"try-error" %in% class(zx.gsva)) {
        zx.gsva <- my.normalize(zx.gsva)
        jj <- match(names(gmt), rownames(zx.gsva))
        zx.gsva <- zx.gsva[jj, colnames(X), drop = FALSE] ## make sure..
        zx.gsva[is.na(zx.gsva)] <- 0
        all.results[["gsva"]] <- playbase::gset.fitContrastsWithLIMMA(
          zx.gsva,
          contr.matrix,
          design = design,
          trend = TRUE,
          conform.output = TRUE
        )
      }
    })
    timings <- rbind(timings, c("gsva", tt))
  }

  if ("ssgsea" %in% methods) {
    message("fitting contrasts using ssGSEA/limma... ")
    tt <- system.time({
      zx.ssgsea <- try(GSVA::gsva(as.matrix(X), gmt[],
        method = "ssgsea",
        parallel.sz = mc.cores, verbose = FALSE
      ))
      if (!"try-error" %in% class(zx.ssgsea)) {
        zx.ssgsea <- my.normalize(zx.ssgsea)
        jj <- match(names(gmt), rownames(zx.ssgsea))
        zx.ssgsea <- zx.ssgsea[jj, colnames(X), drop = FALSE] ## make sure..
        zx.ssgsea[is.na(zx.ssgsea)] <- 0
        all.results[["ssgsea"]] <- playbase::gset.fitContrastsWithLIMMA(
          zx.ssgsea,
          contr.matrix,
          design,
          trend = TRUE,
          conform.output = TRUE
        )
      }
    })
    timings <- rbind(timings, c("ssgsea", tt))
  }


  ## --------------------------------------------------------------
  ## Fit remaining methods
  ## --------------------------------------------------------------

  k <- 1
  fitThisContrastWithMethod <- function(method, k) {
    jj <- which(exp.matrix[, k] != 0)
    yy <- 1 * (exp.matrix[jj, k] > 0)
    xx <- X[, jj]
    ref <- 0
    timings <- c()
    res <- list()

    ## Standard Fisher exact test
    if ("fisher" %in% method) {
      ## calculate significant genes with LIMMA (we need all genes for GSEA-PR)
      lfc <- 0
      lfc05 <- 0.2
      fdr <- 0.25 ## OLD thresholds
      lfc05 <- 0.0
      fdr <- 0.05 ## NEW thresholds (since oct2021)
      suppressWarnings(suppressMessages(
        limma0 <- playbase::gx.limma(xx, yy,
          fdr = 1.0, lfc = 0,
          ref = ref, trend = TRUE, verbose = 0
        ) ## trend true for NGS
      ))
      which.up <- which(limma0[, "adj.P.Val"] <= fdr & limma0[, "logFC"] > lfc05)
      which.dn <- which(limma0[, "adj.P.Val"] <= fdr & limma0[, "logFC"] < -lfc05)
      genes.up <- rownames(limma0)[which.up]
      genes.dn <- rownames(limma0)[which.dn]

      ## Always take at least first 100.. (HACK??!!!)
      if (length(genes.dn) < 100) {
        genes.dn0 <- rownames(limma0)[order(limma0[, "logFC"])]
        genes.dn <- utils::head(unique(c(genes.dn, genes.dn0)), 100)
      }
      if (length(genes.up) < 100) {
        genes.up0 <- rownames(limma0)[order(-limma0[, "logFC"])]
        genes.up <- utils::head(unique(c(genes.up, genes.up0)), 100)
      }

      tt <- system.time({
        output <- playbase::gset.fisher2(genes.up, genes.dn,
          genesets = gmt, fdr = 1.0,
          background = rownames(X), check.background = FALSE,
          min.genes = 0, max.genes = 99999,
          common.genes = FALSE, verbose = 0
        )
      })

      timings <- rbind(timings, c("fisher", tt))
      output <- output[match(names(gmt), rownames(output)), ]
      rownames(output) <- names(gmt)
      output <- output[, c("sign", "p.value", "q.value", "odd.ratio", "overlap")]
      colnames(output) <- c("score", "p.value", "q.value", "odd.ratio", "overlap")
      res[["fisher"]] <- output
    }

    ## ----------------------------------------------------
    ## GSVA-limma
    ## ----------------------------------------------------

    LIMMA.TREND <- TRUE
    if (FALSE && "ssgsea" %in% method) {
      zx <- zx.ssgsea[, colnames(xx)]
      gs <- intersect(names(gmt), rownames(zx))
      tt <- system.time(
        output <- playbase::gx.limma(zx[gs, ], yy, fdr = 1, lfc = 0, ref = ref, trend = LIMMA.TREND, verbose = 0) ## ssgsea
      )
      timings <- rbind(timings, c("ssgsea", tt))

      output <- output[match(names(gmt), rownames(output)), ]
      rownames(output) <- names(gmt)
      output <- output[, c("logFC", "P.Value", "adj.P.Val", "0", "1")]
      colnames(output) <- c("score", "p.value", "q.value", "AveExpr0", "AveExpr1")
      res[["ssgsea"]] <- output
    }

    if (FALSE && "gsva" %in% method) {
      zx <- zx.gsva[, colnames(xx)]
      gs <- intersect(names(gmt), rownames(zx))
      tt <- system.time({
        output <- playbase::gx.limma(zx[gs, ], yy,
          fdr = 1, lfc = 0, ref = ref,
          trend = LIMMA.TREND, verbose = 0
        ) ## ssgsea
      })
      timings <- rbind(timings, c("gsva", tt))
      output <- output[match(names(gmt), rownames(output)), ]
      rownames(output) <- names(gmt)
      output <- output[, c("logFC", "P.Value", "adj.P.Val", "0", "1")]
      colnames(output) <- c("score", "p.value", "q.value", "AveExpr0", "AveExpr1")
      res[["gsva"]] <- output
    }

    if ("spearman" %in% method && !is.null(zx.rnkcorr)) {
      tt <- system.time(
        output <- playbase::gx.limma(zx.rnkcorr[, ], yy, fdr = 1, lfc = 0, ref = ref, trend = LIMMA.TREND, verbose = 0) ## ssgsea
      )
      timings <- rbind(timings, c("spearman", tt))


      output <- output[match(names(gmt), rownames(output)), ]
      rownames(output) <- names(gmt)
      output <- output[, c("logFC", "P.Value", "adj.P.Val", "0", "1")]
      colnames(output) <- c("score", "p.value", "q.value", "AveExpr0", "AveExpr1")
      res[["spearman"]] <- output
    }

    ## ----------------------------------------------------
    ## LIMMA methods
    ## ----------------------------------------------------
    if ("camera" %in% method) {
      cdesign <- cbind(Intercept = 1, Group = yy)
      tt <- system.time({
        suppressWarnings(suppressMessages(
          output <- limma::camera(xx, gmt, cdesign, contrast = 2)
        ))
      })
      timings <- rbind(timings, c("camera", tt))
      ## note: camera does not provide any score!!!!
      output$score <- c(-1, 1)[1 + 1 * (output$Direction == "Up")] * -log10(output$PValue)
      output <- output[match(names(gmt), rownames(output)), ]
      rownames(output) <- names(gmt)
      output <- output[, c("score", "PValue", "FDR", "NGenes", "Direction")]
      colnames(output) <- c("score", "p.value", "q.value", "NGenes", "Direction")
      res[["camera"]] <- output
    }
    if ("fry" %in% method) {
      cdesign <- cbind(Intercept = 1, Group = yy)
      tt <- system.time(
        output <- limma::fry(xx, gmt, cdesign, contrast = 2)
      )
      timings <- rbind(timings, c("fry", tt))
      ## note: camera does not provide any logFC!!!!
      output$score <- c(-1, 1)[1 + 1 * (output$Direction == "Up")] * -log10(output$PValue)
      output <- output[match(names(gmt), rownames(output)), ]
      rownames(output) <- names(gmt)
      output <- output[, c("score", "PValue", "FDR", "NGenes", "Direction")]
      colnames(output) <- c("score", "p.value", "q.value", "NGenes", "Direction")
      res[["fry"]] <- output
    }

    ## fast GSEA
    if ("fgsea" %in% method) {
      rnk <- rowMeans(xx[, which(yy == 1), drop = FALSE]) - rowMeans(xx[, which(yy == 0), drop = FALSE])
      rnk <- rnk + 1e-8 * stats::rnorm(length(rnk))

      tt <- system.time(
        output <- fgsea::fgseaSimple(gmt, rnk,
          nperm = 10000,
          minSize = 1, maxSize = 9999, nproc = 1
        ) ## nproc0 fails!!!
      )
      timings <- rbind(timings, c("fgsea", tt))
      output <- as.data.frame(output)
      output <- output[match(names(gmt), output$pathway), ]
      rownames(output) <- names(gmt)
      output <- output[, c("NES", "pval", "padj", "size", "leadingEdge")]
      colnames(output) <- c("score", "p.value", "q.value", "SIZE", "LEADING GENES")
      res[["fgsea"]] <- output
    }

    res2 <- list(results = res, timings = timings)
    return(res2)
  }


  fitContrastsWithMethod <- function(method) {
    message("fitting contrasts using ", method, "... ")
    results <- list()
    timings <- c()
    k <- 1
    ncontrasts <- ncol(contr.matrix)
    for (k in 1:ncontrasts) {
      res <- fitThisContrastWithMethod(method = method, k)
      results[[k]] <- res$results[[1]]
      timings <- rbind(timings, res$timings)
      names(results)[k] <- colnames(contr.matrix)[k]
    }
    return(list(results = results, timings = timings))
  }

  names(all.results)
  methods2 <- setdiff(methods, names(all.results))
  methods2 <- setdiff(methods2, c("gsva","ssgsea")) ## not do

  for (m in methods2) {
    res <- fitContrastsWithMethod(method = m)
    all.results[[m]] <- res$results
    timings <- rbind(timings, res$timings)
  }

  ## --------------------------------------------------------------
  ## Reshape matrices by comparison
  ## --------------------------------------------------------------

  message("[gset.fitContrastsWithAllMethods] length(all.results) = ", length(all.results))
  tests <- names(all.results[[1]])
  ntest <- length(tests)

  P <- lapply(tests, function(k) do.call(cbind, lapply(all.results, function(x) x[[k]][, "p.value"])))
  Q <- lapply(tests, function(k) do.call(cbind, lapply(all.results, function(x) x[[k]][, "q.value"])))
  S <- lapply(tests, function(k) do.call(cbind, lapply(all.results, function(x) x[[k]][, "score"])))
  for (i in 1:ntest) {
    rownames(P[[i]]) <- names(gmt)
    rownames(Q[[i]]) <- names(gmt)
    rownames(S[[i]]) <- names(gmt)
  }
  names(P) <- names(Q) <- names(S) <- tests

  ## --------------------------------------------------------------
  ## Compute sig counts (by method)
  ## --------------------------------------------------------------
  message("computing sigcounts... ")

  methods <- colnames(Q[[1]])
  nmethod <- length(methods)
  pv.list <- sort(c(1e-16, 10**seq(-8, -2, 2), 0.05, 0.1, 0.2, 0.5, 1))
  sig.counts <- list()
  i <- 1
  fdr <- 0.05
  lfc <- 1e-3
  for (i in 1:nmethod) {
    q0 <- sapply(Q, function(x) x[, i])
    s0 <- sapply(S, function(x) x[, i])
    if (nrow(Q[[1]]) == 1) {
      q0 <- matrix(q0, nrow = nrow(Q[[1]]), dimnames = list(rownames(Q[[1]]), names(Q)))
      s0 <- matrix(s0, nrow = nrow(S[[1]]), dimnames = list(rownames(S[[1]]), names(S)))
    }
    q0[is.na(q0)] <- 1
    s0[is.na(s0)] <- 0
    up0 <- sapply(pv.list, function(p) Matrix::colSums(q0 <= p & s0 > 0, na.rm = TRUE))
    dn0 <- sapply(pv.list, function(p) Matrix::colSums(q0 <= p & s0 < 0, na.rm = TRUE))
    ns0 <- sapply(pv.list, function(p) Matrix::colSums(q0 > p, na.rm = TRUE))
    both0 <- sapply(pv.list, function(p) Matrix::colSums(q0 <= p & abs(s0) > 0, na.rm = TRUE))
    if (ncol(q0) == 1) {
      up0 <- matrix(up0, nrow = 1)
      dn0 <- matrix(dn0, nrow = 1)
      ns0 <- matrix(ns0, nrow = 1)
      both0 <- matrix(both0, nrow = 1)
      rownames(up0) <- rownames(dn0) <- rownames(ns0) <- rownames(both0) <- colnames(q0)[1]
    }
    colnames(up0) <- colnames(dn0) <- colnames(ns0) <- colnames(both0) <- pv.list
    m <- methods[i]
    sig.counts[[m]] <- list(both = both0, up = up0, down = dn0, notsig = ns0)
  }

  ## --------------------------------------------------
  ## meta analysis, aggregate p-values
  ## --------------------------------------------------
  message("computing meta-p values... \n")

  all.meta <- list()
  i <- 1
  for (i in 1:ntest) {
    pv <- P[[i]]
    qv <- Q[[i]]
    fc <- S[[i]]
    meta.p <- apply(pv, 1, max, na.rm = TRUE) ## maximum p-statistic (simple & fast)
    meta.q <- apply(qv, 1, max, na.rm = TRUE) ## maximum q-statistic (simple & fast)
    ss.rank <- function(x) scale(sign(x) * rank(abs(x), na.last = "keep"), center = FALSE)
    if (nrow(S[[i]]) == 1) {
      meta.fx <- S[[i]]
    } else {
      meta.fx <- rowMeans(apply(S[[i]], 2, ss.rank), na.rm = TRUE)
    }
    meta <- data.frame(fx = meta.fx, p = meta.p, q = meta.q)
    rownames(fc) <- NULL ## saves memory...
    rownames(pv) <- NULL
    rownames(qv) <- NULL
    all.meta[[i]] <- data.frame(meta = meta, fc = I(fc), p = I(pv), q = I(qv))
    rownames(all.meta[[i]]) <- rownames(S[[i]])
  }
  names(all.meta) <- tests

  ## --------------------------------------------------
  ## Add meta matrices (this becomes quite large...)
  ## --------------------------------------------------

  message("computing meta-matrix... ")
  m <- list(gsva = zx.gsva, ssgsea = zx.ssgsea, rnkcorr = zx.rnkcorr)
  m <- m[which(!sapply(m, is.null))]

  ## average expression of geneset members
  ng <- Matrix::colSums(G != 0)
  meta.matrix <- as.matrix(Matrix::t(G != 0) %*% X) / ng

  m[["meta"]] <- meta.matrix

  timings <- as.matrix(timings)
  rownames(timings) <- timings[, 1]
  timings0 <- matrix(timings[, -1], nrow = nrow(timings))
  timings0 <- matrix(as.numeric(timings0), nrow = nrow(timings0))
  rownames(timings0) <- rownames(timings)
  if (nrow(timings0) > 1 && sum(duplicated(rownames(timings0)) > 0)) {
    timings0 <- do.call(rbind, tapply(1:nrow(timings0), rownames(timings0), function(i) colSums(timings0[i, , drop = FALSE])))
  }

  res <- list(
    meta = all.meta, sig.counts = sig.counts, outputs = all.results,
    matrices = m, timings = timings0
  )

  return(res)
}


#' Fit Contrasts with LIMMA
#'
#' This function fits contrasts using the LIMMA package and performs statistical analysis on gene expression data.
#'
#' @param gsetX A gene expression dataset.
#' @param contr.matrix Contrast matrix specifying the contrasts of interest.
#' @param design Design matrix specifying the experimental design. Default is \code{NULL}.
#' @param trend Logical indicating whether to use trend estimation in the empirical Bayes smoothing. Default is \code{TRUE}.
#' @param conform.output Logical indicating whether to conform the output table columns. Default is \code{FALSE}.
#'
#' @return A list of tables containing the results of statistical analysis for each contrast.
#'
#' @export
gset.fitContrastsWithLIMMA <- function(gsetX, contr.matrix, design,
                                       trend = TRUE, conform.output = FALSE) {
  if (!is.null(design)) {
    message("fitting gset.LIMMA contrasts with design matrix... ")
    vfit <- limma::lmFit(gsetX, design)
    vfit <- limma::contrasts.fit(vfit, contrasts = contr.matrix)
    efit <- limma::eBayes(vfit, trend = trend, robust = TRUE)

    tables <- list()
    i <- 1
    exp.matrix <- (design %*% contr.matrix)
    for (i in 1:ncol(contr.matrix)) {
      top <- limma::topTable(efit, coef = i, sort.by = "none", number = Inf, adjust.method = "BH")
      rownames(top) <- rownames(gsetX)
      j1 <- which(exp.matrix[, i] > 0)
      j0 <- which(exp.matrix[, i] < 0)
      mean1 <- rowMeans(gsetX[, j1, drop = FALSE], na.rm = TRUE)
      mean0 <- rowMeans(gsetX[, j0, drop = FALSE], na.rm = TRUE)
      ## top <- top[rownames(gsetX), , drop = FALSE]
      top <- cbind(top, "AveExpr0" = mean0, "AveExpr1" = mean1)
      tables[[i]] <- top
    }
    names(tables) <- colnames(contr.matrix)
  } else {
    message("fitting gset.LIMMA contrasts without design.... ")
    tables <- list()
    i <- 1
    for (i in 1:ncol(contr.matrix)) {
      design0 <- cbind(1, contr.matrix[, i])
      colnames(design0) <- c("ref", colnames(contr.matrix)[i])
      vfit <- limma::lmFit(gsetX, design0)
      efit <- limma::eBayes(vfit, trend = trend, robust = TRUE)
      top <- limma::topTable(efit, coef = 2, sort.by = "none", number = Inf, adjust.method = "BH")
      j1 <- which(contr.matrix[, i] > 0)
      j0 <- which(contr.matrix[, i] < 0)
      mean1 <- rowMeans(gsetX[, j1, drop = FALSE], na.rm = TRUE)
      mean0 <- rowMeans(gsetX[, j0, drop = FALSE], na.rm = TRUE)
      top <- top[rownames(gsetX), ]
      top <- cbind(top, "AveExpr0" = mean0, "AveExpr1" = mean1)
      Matrix::head(top, 10)
      tables[[i]] <- top
    }
    names(tables) <- colnames(contr.matrix)
  }

  if (conform.output == TRUE) {
    for (i in 1:length(tables)) {
      jj <- match(rownames(gsetX), rownames(tables[[i]]))
      k1 <- c("logFC", "P.Value", "adj.P.Val", "AveExpr0", "AveExpr1")
      k2 <- c("score", "p.value", "q.value", "AveExpr0", "AveExpr1")
      tables[[i]] <- tables[[i]][jj, k1]
      colnames(tables[[i]]) <- k2
    }
  }

  return(tables)
}


#' @export
shortstring <- function(s, n) {
  s <- as.character(s)
  ifelse(nchar(s) <= n, s, paste0(substring(s, 1, n), "..."))
}


## ======================================================================
## ======================================================================
## ======================================================================
