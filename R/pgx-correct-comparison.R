# This file is part of the Omics Playground project.
# Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
# Batch-correction method comparison and evaluation.
# This file owns multi-method execution, scoring, and method selection.
# Comparison routines delegate correction arithmetic to established implementations.

## =============================================================================
## Run/compare multiple batch-correction methods
## =============================================================================

#' @export
runBatchCorrectionMethods <- function(X, batch, y, controls = NULL, ntop = 2000,
                                      sc = FALSE, prefix = "",
                                      methods = NULL, remove.failed = TRUE) {
  if (0) {
    controls <- NULL
    ntop <- 2000
    sc <- FALSE
    prefix <- ""
    methods <- NULL
    remove.failed <- TRUE
  }

  if (is.null(y)) {
    if ("uncorrected" %in% methods) {
      if (ntop < Inf) {
        X <- head(X[order(-matrixStats::rowSds(X, na.rm = TRUE)), ], ntop) ## faster
      }
      xlist[["uncorrected"]] <- X
    } else {
      return(NULL)
    }
  }

  mod <- model.matrix(~ factor(y))
  nlevel <- length(unique(y[!is.na(y)]))
  if (ntop < Inf) {
    X <- head(X[order(-matrixStats::rowSds(X, na.rm = TRUE)), ], ntop) ## faster
  }

  if (is.null(methods)) {
    methods <- c(
      "uncorrected", "normalized_to_control",
      "ComBat", "limma", "ComBat.no_mod", "limma.no_mod",
      "superBC", "PCA", "RUV", "SVA", "NPM", "MNN", "Harmony"
    )
  }

  xlist <- list()

  if ("uncorrected" %in% methods) {
    xlist[["uncorrected"]] <- X
  }

  ## --------------------------------------------------------------
  ## SUPERVISED METHODS (need batch parameter)
  ## --------------------------------------------------------------

  ## normalize to control
  if (!is.null(controls) && "normalized_to_control" %in% methods) {
    nX <- normalizeToControls(X, batch, y, controls)
    xlist[["normalized_to_control"]] <- nX
  }

  ## limma -------------------------------------------------------
  if ("limma" %in% methods && is.null(batch)) {
    xlist[["limma"]] <- X
  }
  if ("limma" %in% methods && !is.null(batch)) {
    cX <- try(playbase.preprocess::pp.batchCorrect(
      X,
      target = y,
      batch = batch,
      method = "limma"
    ))
    xlist[["limma"]] <- cX
  }
  if ("limma.no_mod" %in% methods && is.null(batch)) {
    xlist[["limma.no_mod"]] <- X
  }
  if ("limma.no_mod" %in% methods && !is.null(batch)) {
    cX <- try(playbase.preprocess::pp.batchCorrect(
      X,
      batch = batch,
      method = "limma"
    ))
    xlist[["limma.no_mod"]] <- cX
  }

  ## ComBat ------------------------------------------------------
  if ("ComBat" %in% methods && is.null(batch)) {
    xlist[["ComBat"]] <- X
  }
  if ("ComBat" %in% methods && !is.null(batch)) {
    if (max(table(batch), na.rm = TRUE) > 1) {
      bX <- try(playbase.preprocess::pp.batchCorrect(
        X,
        target = y,
        batch = batch,
        method = "ComBat"
      ))
      xlist[["ComBat"]] <- bX
    }
  }
  if ("ComBat.no_mod" %in% methods && !is.null(batch)) {
    bX <- try(playbase.preprocess::pp.batchCorrect(
      X,
      batch = batch,
      method = "ComBat"
    ))
    xlist[["ComBat.no_mod"]] <- bX
  }
  if ("ComBat.no_mod" %in% methods && is.null(batch)) {
    xlist[["ComBat.no_mod"]] <- X
  }

  ## superbatchcorrect
  if ("superBC" %in% methods) {
    df <- data.frame(y = y)
    if (!is.null(batch)) df <- cbind(df, batch = batch)
    xlist[["superBC"]] <- pgx.superBatchCorrect(X, df, model.par = "y", batch.par = "*")$X
  }

  ## --------------------------------------------------------------
  ## UNSUPERVISED METHODS (need pheno vector)
  ## --------------------------------------------------------------

  ## PCA
  if ("PCA" %in% methods) {
    xlist[["PCA"]] <- try(pcaCorrect(X, y = y, p.notsig = 0.20))
  }

  ## RUV and SVA
  if ("RUV" %in% methods) {
    xlist[["RUV"]] <- try(playbase.preprocess::pp.batchCorrect(
      X,
      target = y,
      method = "RUV"
    ))
  }

  if ("SVA" %in% methods) {
    xlist[["SVA"]] <- try(playbase.preprocess::pp.batchCorrect(
      X,
      target = y,
      method = "SVA"
    ))
  }

  if ("NPM" %in% methods) {
    xlist[["NPM"]] <- try(playbase.preprocess::pp.batchCorrect(
      X,
      target = y,
      method = "NPM"
    ))
    ## xlist[["NNM2"]] <- nnmCorrect2(X, y, use.design = TRUE)
    ##    xlist[["NNM.no_mod"]] <- nnmCorrect2(X, y, use.design = FALSE)
  }

  ## --------------------------------------------------------------
  ## scRNA-seq methods
  ## --------------------------------------------------------------

  if (sc) {
    ## MNN ---------------------------------------------------------
    if ("MNN" %in% methods) {
      xlist[["MNN"]] <- try(MNNcorrect(X, batch))

      ## restricted MNN ----------------------------------------------
      if (!is.null(controls)) {
        kk <- (y %in% controls)
        xlist[["rMNN"]] <- try(MNNcorrect(X, batch, controls = kk))
      }
    }

    ## Harmony ----------------------------------------------------
    if ("Harmony" %in% methods) {
      res <- try(runHarmony(X, batch = batch))
      if (!"try-error" %in% class(res)) {
        xlist[["Harmony"]] <- as.matrix(res$corrected)
      }
    }
  }

  if (remove.failed) {
    is.error <- sapply(xlist, function(x) ("try-error" %in% class(x)))
    is.nullrow <- sapply(sapply(xlist, nrow), is.null)
    is.xnull <- sapply(xlist, is.null)
    xlist <- xlist[which(!is.xnull & !is.nullrow & !is.error)]
  }

  names(xlist) <- paste0(prefix, names(xlist))
  xlist
}

#' @export
runTechCorrectionMethods <- function(X, samples, y, p.pca = 0.5, p.pheno = 0.05, nv = 1,
                                     xrank = NULL, force = FALSE,
                                     remove.failed = TRUE, ntop = Inf) {
  ##  p.pca = 0.5;p.pheno = 0.05;nv = 2;remove.failed = TRUE;ntop = Inf

  mod <- model.matrix(~y)
  nlevel <- length(unique(y[!is.na(y)]))
  X <- head(X[order(-matrixStats::rowSds(X, na.rm = TRUE)), ], ntop) ## faster

  xlist <- list()
  xlist[["uncorrected"]] <- X

  params <- c("lib", "gender", "mito", "ribo", "cellcycle")

  xlist[["lib"]] <- removeTechnicalEffects(
    X, samples,
    y = y, params = c("lib"),
    p.pheno = p.pheno, p.pca = p.pca, nv = nv,
    xrank = xrank, force = force
  )

  xlist[["gender"]] <- removeTechnicalEffects(
    X, samples,
    y = y, params = c("gender"),
    p.pheno = p.pheno, p.pca = p.pca, nv = nv,
    xrank = xrank, force = force
  )

  xlist[["mito"]] <- removeTechnicalEffects(
    X, samples,
    y = y, params = c("mito"),
    p.pheno = p.pheno, p.pca = p.pca, nv = nv,
    xrank = xrank, force = force
  )

  xlist[["ribo"]] <- removeTechnicalEffects(
    X, samples,
    y = y, params = c("ribo"),
    p.pheno = p.pheno, p.pca = p.pca, nv = nv,
    xrank = xrank, force = force
  )

  xlist[["cellcycle"]] <- removeTechnicalEffects(
    X, samples,
    y = y, params = c("cellcycle"),
    p.pheno = p.pheno, p.pca = p.pca, nv = nv,
    xrank = xrank, force = force
  )

  xlist[["<all>"]] <- removeTechnicalEffects(
    X, samples,
    y = y,
    params = c("lib", "gender", "mito", "ribo", "cellcycle"),
    p.pheno = p.pheno, p.pca = p.pca, nv = nv,
    xrank = xrank, force = force
  )

  if (remove.failed) {
    xlist <- xlist[!sapply(sapply(xlist, nrow), is.null)]
  }

  xlist
}

#' @export
bc.evaluateResults <- function(xlist, pheno, lfc = 0.2, q = 0.2, pos = NULL,
                               add.sil = TRUE, plot = TRUE, trend = TRUE,
                               ref = "uncorrected", clust = "tsne") {
  if (0) {
    lfc <- 0.2
    q <- 0.2
    pos <- NULL
    add.sil <- TRUE
    plot <- TRUE
    trend <- TRUE
    clust <- "tsne"
  }

  if (!ref %in% names(xlist)) ref <- names(xlist)[1]

  ## compute and make table
  message("[bc.evaluateResults] computing statistics...")
  numsig <- lapply(xlist, stats.numsig,
    y = pheno, lfc = lfc, q = q,
    trend = trend, verbose = FALSE
  )

  res <- t(sapply(numsig, function(r) {
    c(sapply(r[1:2], length), avg.fc = mean(abs(r[[3]]), na.rm = TRUE))
  }))

  sdx <- sapply(xlist, function(x) mean(matrixStats::rowSds(x, na.rm = TRUE)))
  snr <- res[, "avg.fc"] / sdx
  res <- cbind(res, avg.sd = sdx, SNR = snr)

  ## compute relative genes/geneset overlap
  message("[bc.evaluateResults] computing overlap...")
  g1 <- numsig[[ref]]$genes
  n1 <- sapply(numsig, function(s) length(intersect(s$genes, g1)))
  ##  n2 <- sapply(numsig, function(s) length(union(s$genes, g1)))
  n2 <- sapply(numsig, function(s) length(g1))
  ##  res <- cbind(res, r.genes=r1, r.gsets=r2, s.genes=s1, s.gsets=s2)
  r.genes <- n1 / (1e-3 + n2)
  res <- cbind(res, r.genes)

  any.gsets <- any(sapply(numsig, function(s) length(s$gsets) > 0))
  if (any.gsets) {
    s1 <- numsig[[ref]]$gsets
    m1 <- sapply(numsig, function(s) length(intersect(s$gsets, s1)))
    ##    m2 <- sapply(numsig, function(s) length(union(s$gsets, s1)))
    m2 <- sapply(numsig, function(s) length(s1))
    r.gsets <- m1 / (1e-3 + m2)
    res <- cbind(res, r.gsets)
  }

  ## centered top
  xlist1 <- lapply(xlist, function(x) {
    x <- head(x[order(-matrixStats::rowSds(x, na.rm = TRUE)), ], 1000)
    x <- as.matrix(x)
    (x - rowMeans(x, na.rm = TRUE))
  })

  message("[bc.evaluateResults] computing silhouette scores...")
  silhouette <- rep(1, nrow(res))
  if (add.sil) {
    if (is.null(pos)) {
      if (clust == "tsne") {
        nb <- max(0.33, min(30, round(ncol(xlist[[1]]) / 5)))
        ## CLUSTFUN <- function(x) uwot::tumap(scale(t(x), scale = FALSE), n_neighbors = nb)
        if (ncol(xlist[[1]]) <= 6) {
          nb <- 0.5
        }
        CLUSTFUN <- function(x) {
          Rtsne::Rtsne(scale(t(x)),
            check_duplicates = FALSE,
            perplexity = nb
          )$Y
        }
      } else {
        CLUSTFUN <- function(x) svd(scale(t(x), scale = FALSE))$u[, 1:2]
      }
      pos <- lapply(xlist1, function(x) CLUSTFUN(x))
    }
    pheno0 <- as.character(pheno)
    pheno0[is.na(pheno0)] <- "NA"
    silhouette <- sapply(pos, function(p) {
      score <- cluster::silhouette(as.integer(factor(pheno0)), stats::dist(p))
      mean(score[, "sil_width"], na.rm = TRUE)
    })
    silhouette <- pmax(silhouette, 1e-4)

    ## PCA score
    nu <- max(2, min(10, dim(xlist[[1]]) / 4))
    pca10 <- lapply(xlist1, function(x) {
      ## svd(scale(t(x), scale = FALSE), nu=nu, nv=0)$u
      svd(t(x), nu = nu, nv = 0)$u
    })
    Y <- model.matrix(~pheno)[, -1]
    rho <- lapply(pca10, function(x) cor(x, Y))
    rho <- lapply(rho, function(x) rowMeans(abs(x), na.rm = TRUE))
    pc1.ratio <- sapply(rho, function(r) abs(r[1]) / sum(abs(r)))

    res <- cbind(res, silhouette, pc1.ratio)
  }

  ## use only these for score
  sel <- c("genes", "gsets", "SNR", "pc1.ratio", "silhouette")
  sel <- intersect(sel, colnames(res))

  ##  score <- res.score * (silhouette / silhouette[1])**1
  overall.score <- t(t(1e-4 + res[, sel]) / (1e-4 + res[ref, sel]))
  overall.score[, "silhouette"] <- overall.score[, "silhouette"]**2 ## give more weight
  overall.score <- exp(rowMeans(log(overall.score), na.rm = TRUE)) ## geometric mean

  res1 <- cbind(score = overall.score, res)
  res1 <- res1[order(-res1[, "score"]), ]
  pos <- pos[rownames(res1)]

  if (plot) {
    nc <- ceiling(1.2 * sqrt(length(pos)))
    nr <- ceiling(length(pos) / nc)
    i <- 1
    xdim <- nrow(pos[[1]])
    cex1 <- cut(xdim, breaks = c(0, 20, 100, 400, 1000, 999999), c(1.8, 1.5, 1.2, 0.9, 0.6))
    cex1 <- as.numeric(as.character(cex1))

    par(mfrow = c(nr, nc))
    for (i in 1:length(pos)) {
      plot(pos[[i]][, 1:2],
        col = factor(pheno), pch = 20, cex = cex1,
        main = names(pos)[i], cex.main = 1.6
      )
      tt <- paste("score = ", round(res1[i, "score"], 3))
      legend("topright", legend = tt, cex = 1.1)
    }
  }

  p.genes <- lapply(numsig, function(s) s$p.genes)
  p.gsets <- lapply(numsig, function(s) s$p.gsets)

  list(scores = res1, pos = pos, p.genes = p.genes, p.gsets = p.gsets)
}

#' @export
compare_batchcorrection_methods <- function(X,
                                            samples,
                                            pheno,
                                            contrasts,
                                            methods = c(
                                              "uncorrected",
                                              "ComBat", "limma", "RUV", "SVA", "NPM"
                                            ),
                                            batch.pars = "<autodetect>",
                                            clust.method = "tsne",
                                            ntop = 4000,
                                            xlist.init = list(),
                                            ref = NULL,
                                            evaluate = TRUE,
                                            npc = 2) {
  if (is.null(pheno) && is.null(contrasts)) {
    stop("must give either pheno vector or contrasts matrix")
  }

  pars <- get_model_parameters(X, samples, pheno = pheno, contrasts = contrasts)
  if (length(batch.pars) && batch.pars[1] %in% c("autodetect", "<autodetect>")) {
    batch.pars <- pars$batch.pars
  }

  batch.pars <- intersect(batch.pars, colnames(samples))
  if (!is.null(batch.pars) && length(batch.pars)) {
    B <- samples[, batch.pars, drop = FALSE]
  } else {
    B <- NULL
  }

  nmissing <- sum(is.na(X))
  if (nmissing) message("WARNING: missing values in X. some methods may fail")

  message("Running batch-correction methods...")
  xlist <- runBatchCorrectionMethods(
    X = X,
    batch = B,
    y = pars$pheno,
    controls = NULL,
    methods = methods,
    ntop = ntop,
    sc = FALSE,
    remove.failed = TRUE
  )

  if (length(xlist.init) > 0) xlist <- c(xlist.init, xlist)
  common_rows <- Reduce(intersect, lapply(xlist, rownames))
  xlist <- c(lapply(xlist, function(x) x[common_rows, , drop = FALSE]))
  xlist <- xlist[order(names(xlist))]

  ## PCA is faster than UMAP
  pos <- NULL
  pca.varexp <- NULL
  loadings <- NULL
  pheno.cor <- NULL
  t2 <- double_center_scale_fast
  if (clust.method == "tsne" && nmissing == 0) {
    message("Computing t-SNE clustering...")
    nb <- max(0.33, round(min(30, dim(X) / 5)))
    if (ncol(X) <= 6) nb <- 0.5
    pos <- lapply(xlist, function(x) {
      Rtsne::Rtsne(t2(x), perplexity = nb, check_duplicates = FALSE)$Y
    })
  } else {
    message("Computing PCA clustering...")
    npc_eff <- max(2, min(npc, min(sapply(xlist, function(x) min(dim(x)))) - 1))
    ## numeric-code sample annotations once (factors -> level codes, as
    ## PCAtools::eigencorplot does); drop constant columns that cannot correlate
    meta_num <- NULL
    if (!is.null(samples) && ncol(samples)) {
      ## Build the numeric annotation matrix to correlate against the PCs.
      ## Numeric columns are kept as-is; categorical columns are one-hot encoded
      ## (one binary indicator per level) so each level gets its own point-
      ## biserial correlation / arrow, rather than arbitrary integer level codes
      ## that would impose a fake ordering on unordered categories. Indexing
      ## columns by name (not sapply over the whole object) also keeps this
      ## correct when samples is a character matrix rather than a data.frame.
      onehot_or_num <- function(cn) {
        v <- samples[, cn]
        ## numeric if it already is, or if every non-missing value coerces
        ## cleanly -- covers numeric columns stored as strings when samples is a
        ## character matrix rather than a data.frame
        vn <- suppressWarnings(as.numeric(as.character(v)))
        if (is.numeric(v) || !any(is.na(vn) & !is.na(v))) {
          return(matrix(vn, ncol = 1, dimnames = list(NULL, cn)))
        }
        v <- as.factor(v)
        ## skip pure-identifier columns (every sample is its own level)
        if (nlevels(v) >= sum(!is.na(v))) {
          return(NULL)
        }
        m <- vapply(levels(v), function(l) as.numeric(v == l), numeric(length(v)))
        colnames(m) <- paste0(cn, "=", levels(v))
        m
      }
      meta_num <- do.call(cbind, lapply(colnames(samples), onehot_or_num))
      if (!is.null(meta_num)) {
        rownames(meta_num) <- rownames(samples)
        keep <- apply(meta_num, 2, function(x) length(unique(x[!is.na(x)])) > 1)
        meta_num <- meta_num[, keep, drop = FALSE]
      }
    }
    for (i in 1:length(xlist)) {
      set.seed(1234)
      M <- t2(xlist[[i]])
      pca <- irlba::irlba(M, nu = npc_eff, nv = npc_eff)
      U <- pca$u[, seq_len(npc_eff), drop = FALSE]
      rownames(U) <- colnames(xlist[[i]])
      pos[[names(xlist)[i]]] <- U
      ## % variance explained relative to TOTAL variance (sum(M^2) == sum of all
      ## eigenvalues), not just the top npc_eff PCs irlba returns, so the values
      ## and their cumulative sum are honest (else they inflate and cum hits 100)
      pca.varexp[[names(xlist)[i]]] <- (pca$d^2 / sum(M^2)) * 100
      v <- pca$v[, seq_len(npc_eff), drop = FALSE]
      rownames(v) <- rownames(xlist[[i]])
      loadings[[names(xlist)[i]]] <- v
      ## correlation of each numeric-coded annotation with each PC, for the
      ## phenotype-projection biplot in the UI (annotations x npc)
      if (!is.null(meta_num) && ncol(meta_num)) {
        mm <- meta_num[rownames(U), , drop = FALSE]
        pheno.cor[[names(xlist)[i]]] <-
          t(suppressWarnings(stats::cor(U, mm, use = "pairwise.complete.obs")))
      }
    }
  }

  for (i in 1:length(pos)) {
    rownames(pos[[i]]) <- colnames(X)
  }

  res <- NULL
  best.method <- ref

  if (evaluate) {
    ## compare results using scoring
    res <- bc.evaluateResults(
      xlist,
      pheno = pars$pheno,
      lfc = 0.2,
      q = 0.05,
      pos = pos,
      add.sil = TRUE,
      plot = FALSE,
      trend = TRUE
    )
    ## shiny::removeModal()
    score <- res$scores[, "score"]
    if (is.null(ref)) ref <- names(xlist)[1]
    best.method <- names(which.max(score))

    ## if the improvement is small, we rather choose the uncorrected solution
    score.ratio <- score[best.method] / score[ref]
    best.method <- ifelse(score.ratio < 1.20, ref, best.method)
    message("[select_batchcorrect_method] best.method = ", best.method)
  }

  list(
    xlist = xlist,
    pos = pos,
    pca.varexp = pca.varexp,
    loadings = loadings,
    pheno.cor = pheno.cor,
    scores = res$scores,
    pheno = pars$pheno,
    pars = pars,
    best.method = best.method
  )
}
