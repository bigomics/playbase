# This file is part of the Omics Playground project.
# Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
# Batch-correction diagnostics and technical-effect analysis.
# This file owns confounding checks, effect detection, and vector estimation.
# Diagnostic routines must not implement production correction kernels.

#' @title Check for confounders
#'
#' @description
#' Performs a correlation analysis on a phenotype matrix to detect possible confounding factors.
#'
#' @param pheno Dataframe containing sample metadata/covariates.
#' @param model.par Vector of column names of covariates of interest \code{pheno}.
#' @param max.rho Maximum allowed correlation.
#'
#' @details This function performs a correlation analysis on the factor with respect to the factor
#' of interest. Correlation is compute on each level of all factors and a factor is deemed confounding
#' if any level has a correlation larger than max.rho. For example, such factors should not be used
#' for batch correction as they may remove 'genuine' signal correlated with the factor of interest.
#'
#' @return List of confounding and not-confounding factors. Correlation matrix rho.
#'
#' @export
checkConfounders <- function(pheno, model.par, max.rho = 0.3) {
  getModelMatrix <- function(v) {
    y <- as.character(pheno[, v])
    y[is.na(y)] <- "NA" ## or impute???
    m1 <- stats::model.matrix(~y)[, -1, drop = FALSE]
    colnames(m1) <- sub("^y", paste0(v, "="), colnames(m1))
    m1
  }

  mod1 <- do.call(cbind, lapply(model.par, getModelMatrix))
  rownames(mod1) <- rownames(pheno)
  mod1

  ## --------------------------------------------------------------------
  ## Check confounding
  ## --------------------------------------------------------------------
  batch.pars <- setdiff(colnames(pheno), model.par)
  if (is.null(batch.pars) || length(batch.pars) == 0) {
    return(c())
  }

  mod0 <- do.call(cbind, lapply(batch.pars, getModelMatrix))
  rho <- stats::cor(mod0, mod1)
  rho[is.na(rho)] <- 0
  confounding.pars <- c()

  if (max(abs(rho), na.rm = TRUE) > max.rho) {
    idx <- which(abs(rho) > max.rho, arr.ind = TRUE)
    idx
    for (i in 1:nrow(idx)) {
      v0 <- colnames(mod0)[idx[i, 1]]
      v1 <- colnames(mod1)[idx[i, 2]]
      dbg(paste0(
        "WARNING:: '", v0, "' is confounded with '", v1, "' ",
        ": rho= ", round(rho[idx[i, 1], idx[i, 2]], 3), "\n"
      ))
    }
    confounding <- colnames(mod0)[idx[, 1]]
    confounding <- unique(gsub("=.*", "", confounding))
    dbg("WARNING:: removing confounding batch factors:", confounding, "\n")
    confounding.pars <- c(confounding.pars, confounding)
    batch.pars <- setdiff(batch.pars, confounding.pars)
  }

  list(
    confounding = confounding.pars,
    not.confounding = batch.pars,
    rho = rho
  )
}

#' Correlate principal components with phenotypes
#'
#' @param X Expression matrix
#' @param Y Data frame of sample phenotypes
#' @param nv Number of principal components to use
#' @param stat Statistic for categorical phenotypes ("F" or "t").
#' @param plot Logical, whether to generate a PCA plot.
#' @param main Title for PCA plot.
#'
#' @return Named vector of correlation coefficients
#'
#' @description
#' Calculates the correlation between principal components of the expression matrix X
#' and sample phenotypes in Y.
#'
#' @details
#' This function calculates the top nv principal components of the expression matrix X
#' using irlba::irlba. It then correlates each PC with each phenotype in the Y data frame.
#'
#' For categorical phenotypes, it converts to a factor and calculates the correlation with
#' the model matrix. For numeric phenotypes it calculates the standard correlation coefficient.
#'
#' It returns a named vector of correlation coefficients, with names corresponding to
#' the phenotypes.
#'
#' @export
pgx.PC_correlation <- function(X, Y, nv = 3, stat = "F",
                               expand = FALSE, collapse = TRUE,
                               plot = TRUE, horiz = FALSE,
                               main = NULL, text.cex = 1) {
  getF <- function(x, y) {
    x <- t(scale(t(x))) ## rowscale
    ii <- which(!is.na(y))
    y1 <- y[ii]
    if (inherits(y1, c("factor", "character", "logical"))) {
      y1 <- factor(as.character(y1))
    } else {
      y1 <- y1 + 1e-8 * stats::rnorm(length(y1))
      y1 <- (y1 > stats::median(y1, na.rm = TRUE))
    }
    design <- stats::model.matrix(~y1)
    fit <- limma::lmFit(x[, ii], design)
    suppressWarnings(fit <- try(limma::eBayes(fit, trend = FALSE)))
    if (inherits(fit, "try-error")) {
      return(NULL)
    }
    suppressMessages(top <- limma::topTableF(fit, number = nrow(x), sort.by = "none"))
    list(val = top$F, pv = top$P.Value)
  }
  getCor <- function(x, y) {
    ii <- which(!is.na(y))
    y1 <- y[ii]
    if (inherits(y1, "factor")) y1 <- factor(as.character(y1))
    design <- stats::model.matrix(~ 0 + y1)
    r1 <- stats::cor(t(x[, ii]), design, use = "pairwise")
    r1 <- rowMeans(abs(r1), na.rm = TRUE)
    pv <- cor.pvalue(r1, length(y1))
    list(val = r1, pv = pv)
  }

  if (expand) {
    Y <- expandPhenoMatrix(Y, drop.ref = FALSE)
  }
  X <- X - rowMeans(X, na.rm = TRUE) ## center features
  X[is.na(X)] <- mean(X, na.rm = TRUE) ## no missing allowed
  nv <- min(nv, dim(X) - 1)
  if (nv < ncol(X) / 2) {
    set.seed(1234)
    V <- irlba::irlba(X, nv = nv, nu = 0)$v
  } else {
    V <- svd(X, nv = nv, nu = 0)$v
  }
  rownames(V) <- colnames(X)
  colnames(V) <- paste0("PC", 1:ncol(V))

  rho <- list()
  pv <- list()
  p <- "Chemotherapy"
  for (p in c("<random>", colnames(Y))) {
    if (p == "<random>") {
      y <- sample(c("a", "b"), ncol(X), replace = TRUE)
    } else {
      y <- Y[, p]
    }
    nlevels <- length(unique(y[!is.na(y)]))
    if (nlevels > 1) {
      if (stat == "cor") {
        res <- getCor(x = t(V), y)
        rho[[p]] <- res$val
        pv[[p]] <- res$pv
      }
      if (stat == "F") {
        res <- getF(x = t(V), y)
        rho[[p]] <- res$val
        pv[[p]] <- res$pv
      }
    }
  }
  R <- do.call(rbind, rho)
  P <- do.call(rbind, pv)
  colnames(R) <- paste0("PC", 1:ncol(R))
  colnames(P) <- paste0("PC", 1:ncol(P))

  if (!horiz && stat == "F") {
    R <- t(t(R) / colMeans(R, na.rm = TRUE))
  }

  if (collapse) {
    params <- gsub("[=:].*", "", rownames(R))
    #    params <- sub("^[.]","_",params)
    #    params <- sub("[.].*","",params)

    rr <- tapply(1:nrow(R), params, function(i) colMeans(R[i, , drop = FALSE]))
    R <- do.call(rbind, rr)
    ##    pp <- tapply(1:nrow(P), params, function(i) apply(P[i,,drop=FALSE],2,min))
    pp <- tapply(1:nrow(P), params, function(i) exp(colMeans(log(1e-20 + P[i, , drop = FALSE]))))
    P <- do.call(rbind, pp)
  }

  if (plot) {
    stat0 <- c("correlation", "F-statistic")[1 + 1 * (stat == "F")]
    tt0 <- c("PC correlation", "PC variation")[1 + 1 * (stat == "F")]
    if (is.null(main)) main <- tt0
    if (horiz) {
      plt <- plot_ggbarplot((R),
        ylab = "", srt = 0, horiz = TRUE,
        legend.cex = 1.0 * text.cex,
        label.cex = 1.15 * text.cex,
        axis.cex = 1.05 * text.cex,
        group.name = ""
      )
      plt <- plt + ggplot2::theme(
        plot.margin = ggplot2::margin(t = 0, r = 4, b = 0, l = 8, "pt"),
        plot.title = ggplot2::element_text(size = 13 * text.cex)
      ) +
        ggplot2::xlab(stat0) + ggplot2::ggtitle(main)
    } else {
      plt <- plot_ggbarplot(t(R),
        ylab = stat0, srt = 45, horiz = FALSE,
        legend.cex = 1.0 * text.cex,
        label.cex = 1.15 * text.cex,
        axis.cex = 1.05 * text.cex,
        group.name = ""
      ) +
        ggplot2::theme(
          plot.margin = ggplot2::margin(t = 0, r = 4, b = 0, l = 8, "pt"),
          plot.title = ggplot2::element_text(size = 13 * text.cex)
        ) +
        ggplot2::xlab("") + ggplot2::ggtitle(main)
    }
    plt
    return(plt)
  }
  list(R = R, P = P, V = V)
}

#' @title Estimate technical effects variation
#'
#' @param X Gene expression matrix, with genes in rows and samples in columns
#' @param is.count Logical indicating if X contains counts (TRUE) or log-expression values (FALSE)
#'
#' @return List containing:
#' \itemize{
#'  \item{pct.mito}{Percent mitochondrial genes}
#'  \item{pct.ribo}{Percent ribosomal genes}
#'  \item{biological}{Biological coefficient of variation}
#' }
#'
#' @description
#' Estimates the biological variation and fraction of mitochondrial and ribosomal genes from a gene expression matrix.
#'
#' @details
#' This function calculates the biological variation (BCV) for each gene as the coefficient of variation of expression across samples.
#'
#' It also calculates the percentage of mitochondrial and ribosomal genes based on gene symbols.
#'
#' If the input matrix X contains counts, it will be transformed to log2-CPM.
#' If X contains log-expression values, it will be shifted to start at the 1% quantile.
#' @examples
#' \dontrun{
#' data(sample.ExpressionSet)
#' results <- pgx.computeTechnicalEffects(sample.ExpressionSet)
#' head(results$biological)
#' }
#' @export
pgx.computeTechnicalEffects <- function(X, is.count = FALSE, nmin = 3, nv = 1) {
  ## estimate biological variation
  ##
  ## X:     log-expression matrix
  ##

  ##  is.count = FALSE; nmin = 3;nv=2
  ##  nmin=3:nv=2

  message("[pgx.computeTechnicalEffects] estimating technical effects...")
  nv <- min(nmin, nv)
  nv <- min(nv, ncol(X) / 2)

  ## shift zero to 1% percentile
  if (!is.count) {
    counts <- pmax(2**X - 1, 0, na.rm = TRUE) ## counts
  } else {
    counts <- X
    X <- log2(counts + 1e-8)
  }

  ## technical covariates
  nfeature <- Matrix::colSums(counts > 0, na.rm = TRUE) + 1
  libsize <- Matrix::colSums(counts, na.rm = TRUE)
  libmedian <- apply(X, 2, median, na.rm = TRUE)
  sel.big3 <- head(order(-rowMeans(counts, na.rm = TRUE)), 3)
  counts.big3 <- colSums(counts[sel.big3, , drop = FALSE], na.rm = TRUE)
  big3 <- log2((1 + counts.big3) / (1 + libsize))

  ## mito/ribo genes
  mt.genes <- grep("^MT-|^ATP[68]$|^COX[1-2]$|^ND[1-6]$|^CYTB$",
    rownames(X),
    ignore.case = TRUE, value = TRUE
  )
  rb.genes <- grep("^RP[SL]", rownames(X), ignore.case = TRUE, value = TRUE)
  mito <- ribo <- NA
  pct.mito <- pct.ribo <- NA

  pheno <- data.frame(
    lib.size = log2(libsize + 1e-8),
    lib.median = libmedian,
    lib.big3 = big3,
    #    mito = mito.pc,
    #    ribo = ribo.pc,
    check.names = FALSE
  )
  ## colnames(pheno) <- c("libsize",colnames(mito.pc),colnames(ribo.pc))
  mt.genes
  rb.genes

  if (length(mt.genes) >= nmin) {
    ## mt.genes <- mt.genes[order(-apply(X[mt.genes, , drop = FALSE], 1, stats::sd, na.rm = TRUE))]
    mt.genes <- mt.genes[order(-matrixStats::rowSds(X[mt.genes, , drop = FALSE], na.rm = TRUE))]
    mt.genes <- head(mt.genes, 20)
    mito <- Matrix::colMeans(X[mt.genes, , drop = FALSE], na.rm = TRUE)
    pct.mito <- Matrix::colSums(counts[mt.genes, , drop = FALSE], na.rm = TRUE) / libsize

    mito.pc <- svd(X[mt.genes, ])$v[, 1:nv]
    if (NCOL(mito.pc) > 1) {
      colnames(mito.pc) <- paste0("PC", 1:ncol(mito.pc))
      rownames(mito.pc) <- colnames(X)
    }
    pheno <- cbind(pheno, mito = mito.pc)
  }

  if (length(rb.genes) >= nmin) {
    ## rb.genes <- rb.genes[order(-apply(X[rb.genes, , drop = FALSE], 1, stats::sd, na.rm = TRUE))]
    rb.genes <- rb.genes[order(-matrixStats::rowSds(X[rb.genes, , drop = FALSE], na.rm = TRUE))]
    rb.genes <- head(rb.genes, 20)
    ribo <- Matrix::colMeans(X[rb.genes, , drop = FALSE], na.rm = TRUE)
    pct.ribo <- Matrix::colSums(counts[rb.genes, , drop = FALSE], na.rm = TRUE) / libsize

    ribo.pc <- svd(X[rb.genes, ])$v[, 1:nv]
    if (NCOL(ribo.pc) > 1) {
      colnames(ribo.pc) <- paste0("PC", 1:ncol(ribo.pc))
      rownames(ribo.pc) <- colnames(X)
    }
    pheno <- cbind(pheno, ribo = ribo.pc)
  }

  cc.score <- try(pgx.scoreCellCycle(counts))
  if (!any(class(cc.score) == "try-error")) {
    cc.score <- cc.score[, c("s_score", "g2m_score")]
    colnames(cc.score) <- c("S", "G2M")
    if (nv == 1) cc.score <- rowMeans(cc.score, na.rm = TRUE)
    pheno <- cbind(pheno, cellcycle = cc.score)
  }

  ## create gender model matrix
  gender. <- pgx.inferGender(counts)
  if (length(table(gender.)) > 1) {
    gender.[is.na(gender.)] <- "NA"
    mod.gender <- model.matrix(~ 0 + gender.)
    mod.gender <- mod.gender[, setdiff(colnames(mod.gender), "gender.NA"), drop = FALSE]
    pheno <- cbind(pheno, mod.gender)
  }

  ## take out constant columns
  sel <- which(apply(pheno, 2, sd, na.rm = TRUE) > 0)
  pheno <- pheno[, sel, drop = FALSE]

  return(pheno)
}

#' @export
detectBatchEffects <- function(X, samples, pheno, contrasts = NULL,
                               params = c("statistical", "technical", "pca"),
                               p.pca = 0.5, p.pheno = 0.05, force = FALSE,
                               k.pca = 10, nv = 1, xrank = NULL) {
  if (0) {
    p.pca <- 0.5
    p.pheno <- 0.05
    k.pca <- 10
    nv <- 2
    contrasts <- NULL
    params <- c("statistical", "technical", "pca")
    params <- c("statistical")
    params <- c("technical")
  }

  if (force == TRUE) {
    p.pca <- 1
    p.pheno <- 0
    ## params = c("statistical", "technical", "pca")
  }

  if (!all(params %in% c("statistical", "technical", "pca"))) {
    params1 <- setdiff(params, c("statistical", "technical", "pca"))
    stop("[detectBatchEffects] unknown parameter type: ", params1)
  }

  Y <- samples[, 0]
  if ("statistical" %in% params) {
    Y <- samples
  }
  if ("technical" %in% params) {
    B <- pgx.computeTechnicalEffects(X, nv = nv)
    Y <- cbind(Y, B)
  }

  V <- NULL
  if ("pca" %in% params) {
    ## add PC components
    cX <- X - rowMeans(X, na.rm = TRUE)
    k.pca <- ceiling(min(k.pca, dim(cX) - 1))
    if (k.pca < min(dim(cX)) / 3) {
      set.seed(1234)
      V <- irlba::irlba(cX, nv = k.pca, nu = 0)$v
    } else {
      V <- svd(cX, nv = k.pca, nu = 0)$v
    }
    rownames(V) <- colnames(cX)
    colnames(V) <- paste0("pca.PC", 1:ncol(V))
    Y <- cbind(Y, V)
  }
  dim(Y)

  if (!is.null(contrasts) && is.null(pheno)) {
    pheno <- contrasts2pheno(contrasts, samples)
  } else if (length(pheno) == 1 && pheno[1] %in% colnames(samples)) {
    pheno <- samples[, pheno]
  } else if (length(pheno) == nrow(samples)) {} else {
    stop("invalid pheno argument type: need pheno vector or contrast matrix")
  }

  ## add pheno vector
  dpheno <- model.matrix(~ 0 + pheno)
  colnames(dpheno) <- sub("^pheno", "*pheno*=", colnames(dpheno))
  dY <- cbind(dpheno, expandPhenoMatrix(Y))
  Y <- cbind("*pheno*" = pheno, Y)

  ## detect possible batch covariates (not correlated with phenotype)
  ## determine batch covariates
  dY <- scale(dY)

  res <- gx.limmaF(t(dY), pheno, fdr = 1, lfc = 0, compute.means = FALSE, verbose = 0)
  param <- sub("=.*", "", rownames(res))
  pv.pheno <- tapply(res$P.Value, param, min)

  ## compute correlation with PC components
  pc <- pgx.PC_correlation(X, dY, nv = k.pca, collapse = FALSE, plot = FALSE)
  pca.pars <- sub("=.*", "", rownames(pc$P))
  P0 <- apply(pc$P, 2, function(x) tapply(x, pca.pars, min))
  pp <- intersect(rownames(P0), names(pv.pheno))
  P0 <- P0[pp, , drop = FALSE]

  ## determine phenotype rank
  if (is.null(xrank)) {
    xrank <- which.min(P0["*pheno*", ])
  }
  xrank <- min(xrank, ncol(P0))
  xrank

  ##  pv.pca.min <- apply(cbind(1, P0)[, 1:xrank, drop = FALSE], 1, min, na.rm = TRUE)
  pv.pca.min <- apply(P0[, 1:xrank, drop = FALSE], 1, min, na.rm = TRUE)
  pv.pca.min

  ## select parameters with significant correlation with PC and
  ## not-significant with phenotype.
  pp <- setdiff(pp, "*pheno*")
  P <- cbind(p.pca = pv.pca.min[pp], p.pheno = pv.pheno[pp])
  P

  params <- names(which(P[, 1] < p.pca & P[, 2] > p.pheno))
  params

  batch.matrix <- NULL
  batch.vec <- NULL
  if (length(params) > 0) {
    batch.matrix <- Y[, params, drop = FALSE]
    batch.matrix <- expandPhenoMatrix(batch.matrix, keep.numeric = TRUE)
    ## determine batch covariates
    B <- 1 * expandPhenoMatrix(batch.matrix)
    B[is.na(B)] <- "x"
    batch.vec <- apply(B, 1, paste, collapse = "")
    table(batch.vec)
  }

  ## divide parameters by technical and statistical
  if (length(params)) {
    stats.params <- intersect(params, colnames(samples))
    pca.params <- grep("^pca.PC", params, value = TRUE)
    tech.params <- setdiff(params, c(pca.params, stats.params))
    params <- list(
      "statistical" = stats.params,
      "technical" = tech.params,
      "pca" = pca.params
    )
  }

  ## covariate plus phenotype
  #  Y2 <- expandPhenoMatrix(Y, keep.numeric = TRUE)
  M <- expandPhenoMatrix(cbind(pheno))
  batch.matrix_plus <- cbind(batch.matrix, pheno = M)

  list(
    params = params,
    batch = batch.vec,
    covariates = batch.matrix,
    covariates_plus = batch.matrix_plus,
    pheno = pheno,
    p.values = P,
    p.pca = P0,
    v.pca = pc$V,
    xrank = xrank,
    Y = Y,
    pc = pc
  )
}

#' @export
get_model_parameters <- function(X, samples, pheno = NULL, contrasts = NULL) {
  if (is.null(pheno) && is.null(contrasts)) {
    stop("must give either pheno vector or contrasts matrix")
  }

  ## if a contrasts is given, create a virtual phenotype that spans
  ## all comparison groups.
  if (!is.null(contrasts)) {
    pheno <- contrasts2pheno(contrasts, samples)
  }

  bc <- detectBatchEffects(X, samples, pheno,
    params = "statistical",
    k.pca = 10, p.pca = 0.5, p.pheno = 0.05, xrank = 10
  )
  bc$p.values

  ## Check for any parameter that is highly correlated with the
  ## 'phenotype' as defined by the contrasts.
  p.pheno <- bc$p.values[, "p.pheno"]
  if (nrow(bc$p.values) == 1) names(p.pheno)[1] <- rownames(bc$p.values)
  p.pheno
  pheno.pars <- names(which(p.pheno == min(p.pheno, na.rm = TRUE) | p.pheno < 1e-80))
  pheno.pars <- pheno.pars[order(p.pheno[pheno.pars])]
  pheno.pars

  batch.pars <- NULL
  if ("statistical" %in% names(bc$params)) {
    batch.pars <- bc$params$statistical
    batch.pars
    batch.pars <- setdiff(batch.pars, pheno.pars)
    ## if there are columns names batch, we add them to batch parameters
    batch.pars2 <- grep("batch", colnames(samples), ignore.case = TRUE, value = TRUE)
    if (length(batch.pars2)) batch.pars <- c(batch.pars, batch.pars2)
    batch.pars <- sort(unique(batch.pars))
  }
  batch.pars

  ## we need to use the discretized phenotype matrix because we do not
  ## want create continuous levels.
  dsamples <- expandPhenoMatrix(samples)
  dsamples.pars <- sub("=.*", "", colnames(dsamples))

  batch.vec <- NULL
  if (length(batch.pars)) {
    ## batch.vec <- apply(samples[, batch.pars, drop = FALSE], 1, paste, collapse = "_")
    sel <- which(dsamples.pars %in% batch.pars)
    pheno.vec <- apply(dsamples[, sel, drop = FALSE], 1, paste, collapse = "_")
  }

  pheno.vec <- NULL
  if (length(pheno.pars)) {
    ## pheno.vec <- apply(samples[, pheno.pars, drop = FALSE], 1, paste, collapse = "_")
    sel <- which(dsamples.pars %in% pheno.pars)
    pheno.vec <- apply(dsamples[, sel, drop = FALSE], 1, paste, collapse = "_")
  }

  list(
    batch.pars = batch.pars,
    pheno.pars = pheno.pars,
    batch = batch.vec,
    pheno = pheno.vec
  )
}

#' Estimate batch correction vectors from corrected cX and uncorrected
#' matrix X.
#'
#' @export
estimateBatchCorrectionVectors <- function(cX, X, k = NULL, threshold = 0.8) {
  res <- svd(X - cX)
  cumcv <- (cumsum(res$d**2) / sum(res$d**2))
  cumcv
  if (is.null(k)) {
    k <- min(which(cumcv >= threshold))
  }
  ## return batch vectors
  res$V[, 1:k, drop = FALSE]
}
