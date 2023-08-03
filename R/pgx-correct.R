##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

## Batch correction methods
##
##


#' @title Supervised batch correction
#'
#' @description
#' Performs supervised batch correction on a gene expression matrix, using known technical factors and biological covariates.
#'
#' @param X Gene expression matrix, genes in rows, samples in columns.
#' @param batch Factor or matrix indicating batch/technical groups for each sample.
#' @param model.par Vector of column names of biological covariates in \code{pheno}.
#' @param pheno Dataframe containing sample metadata/covariates. Must match colnames of \code{X}.
#' @param method Batch correction method to apply ("combat", "limma", "sva", etc).
#'
#' @details This function performs supervised batch correction, using known batch groups and biological covariates.
#' It constructs a design matrix containing batch groups and covariates.
#' The \code{method} batch correction is then applied to the expression matrix \code{X}, using this design matrix.
#'
#' Technical effects are estimated from the data as follows:
#' \enumerate{
#' \item Compute average expression within each batch
#' \item Perform PCA on the batch means
#' \item Extract PCs explaining at least 80\% variance
#' }
#'
#' Biological effects are included using the specified \code{model.par} covariates.
#'
#' @return The batch corrected gene expression matrix.
#'
#' @export
pgx.superBatchCorrect <- function(X, pheno, model.par, partype = NULL,
                                  batch.par = "*",
                                  lib.correct = TRUE,
                                  bio.correct = c("mito", "ribo", "cell_cycle", "gender"),
                                  sva.correct = TRUE, pca.correct = TRUE, hc.correct = TRUE,
                                  mnn.correct = NULL, nnm.correct = TRUE,
                                  max.rho = 0.3, max.iter = 10, hc.top = 50) {
  getModelMatrix <- function(v) {
    y <- as.character(pheno[, v])
    y[is.na(y)] <- "NA" ## or impute???
    m1 <- stats::model.matrix(~y)[, -1, drop = FALSE]
    colnames(m1) <- sub("^y", paste0(v, "="), colnames(m1))
    m1
  }
  if (is.null(model.par) && is.null(batch.par)) {
    stop("ERROR:: model and batch cannot be both NULL")
  }

  ## tidy up pheno matrix?? get correct parameter types

  pheno <- utils::type.convert(pheno, as.is = TRUE)

  message("[pgx.superBatchCorrect] 1 : dim.pheno = ", paste(dim(pheno), collapse = "x"))
  message("[pgx.superBatchCorrect] 1 : model.par = ", paste(model.par, collapse = "x"))
  message("[pgx.superBatchCorrect] 1 : colnames.pheno = ", paste(colnames(pheno), collapse = "x"))

  ## setup model matrix
  mod1 <- NULL
  if (!is.null(model.par) && length(model.par) > 0) {
    model.par <- intersect(model.par, colnames(pheno))
    mod1 <- do.call(cbind, lapply(model.par, getModelMatrix))
    rownames(mod1) <- rownames(pheno)
  }
  model.par

  ## get technical/biological effects
  Y <- pgx.computeBiologicalEffects(X)
  colnames(Y) <- paste0(".", colnames(Y))

  message("[pgx.superBatchCorrect] 2 : dim.Y = ", paste(dim(Y), collapse = "x"))

  ## add to phenotype matrix
  pheno <- cbind(pheno, Y)
  not.na <- colMeans(is.na(pheno)) < 1
  nlev <- apply(pheno, 2, function(x) length(unique(x[!is.na(x)])))
  pheno <- pheno[, which(nlev > 1 & not.na), drop = FALSE]
  partype <- sapply(pheno, class)

  message("[pgx.superBatchCorrect] 3 : dim.pheno = ", paste(dim(pheno), collapse = "x"))

  ## --------------------------------------------------------------------
  ## select parameters
  ## --------------------------------------------------------------------


  ## select all non-model variables
  if (!is.null(batch.par) && batch.par[1] == "*") {
    batch.par <- setdiff(colnames(pheno), model.par)
  }

  if ("mito" %in% bio.correct) {
    b1 <- grep("mito", colnames(pheno), value = TRUE)
    batch.par <- c(batch.par, b1)
  }
  if ("ribo" %in% bio.correct) {
    b1 <- grep("ribo", colnames(pheno), value = TRUE)
    batch.par <- c(batch.par, b1)
  }
  if ("cell_cycle" %in% bio.correct) {
    b1 <- grep("cc[.].*score$", colnames(pheno), value = TRUE) ## only s.score and g2m.score
    batch.par <- c(batch.par, b1)
  }
  if ("gender" %in% bio.correct) {
    b1 <- grep("gender", colnames(pheno), value = TRUE)
    batch.par <- c(batch.par, b1)
  }

  batch.par <- intersect(batch.par, colnames(pheno))
  batch.par <- setdiff(batch.par, model.par)
  batch.par <- setdiff(batch.par, c("group", "cluster", "condition")) ## never???

  ## --------------------------------------------------------------------
  ## guess parameter type
  ## --------------------------------------------------------------------

  ## select which are (continuous) covariates or (discrete) factors
  sel1 <- which(partype %in% c("factor", "character", "discrete", "logical"))
  batch.prm <- intersect(batch.par, names(partype[sel1]))
  sel2 <- which(partype %in% c("integer", "numeric"))
  batch.cov <- intersect(batch.par, names(partype[sel2]))


  model.par <- intersect(model.par, colnames(pheno))
  batch.prm <- intersect(batch.prm, colnames(pheno))
  batch.cov <- intersect(batch.cov, colnames(pheno))
  if (length(model.par) == 0) model.par <- NULL
  if (length(batch.prm) == 0) batch.prm <- NULL
  if (length(batch.cov) == 0) batch.cov <- NULL

  ## --------------------------------------------------------------------
  ## Check confounding
  ## --------------------------------------------------------------------
  if (!is.null(batch.prm) && !is.null(mod1)) {
    mod0 <- do.call(cbind, lapply(batch.prm, getModelMatrix))
    rho <- stats::cor(mod0, mod1)
    rho[is.na(rho)] <- 0
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
      confounding.pars <- colnames(mod0)[idx[, 1]]
      confounding.pars <- unique(gsub("=.*", "", confounding.pars))
      dbg("WARNING:: removing confounding batch factors:", confounding.pars, "\n")
      batch.prm <- setdiff(batch.prm, confounding.pars)
    }
  }

  if (!is.null(batch.cov) && !is.null(mod1)) {
    cvar <- data.matrix(pheno[, batch.cov])
    rho1 <- stats::cor(cvar, mod1, use = "pairwise")
    rho1
    rho1[is.na(rho1)] <- 0
    if (max(abs(rho1), na.rm = TRUE) > max.rho) {
      idx <- which(abs(rho1) > max.rho, arr.ind = TRUE)
      for (i in 1:nrow(idx)) {
        v0 <- colnames(cvar)[idx[i, 1]]
        v1 <- colnames(mod1)[idx[i, 2]]
        dbg(paste0(
          "WARNING:: '", v0, "' is confounded with '", v1, "' ",
          ": rho= ", round(rho1[idx[i, 1], idx[i, 2]], 3), "\n"
        ))
      }
      confounding.cov <- colnames(cvar)[idx[, 1]]
      confounding.cov <- unique(gsub("=.*", "", confounding.cov))
      dbg("WARNING:: removing confounding batch covariates:", confounding.cov, "\n")
      batch.cov <- setdiff(batch.cov, confounding.cov)
    }
  }

  dbg("[pgx.superBatchCorrect] model.par=", model.par, "\n")
  dbg("[pgx.superBatchCorrect] batch.prm=", batch.prm, "\n")
  dbg("[pgx.superBatchCorrect] batch.cov=", batch.cov, "\n")

  cX <- X
  mod1x <- matrix(1, ncol(cX), 1)
  if (!is.null(mod1)) mod1x <- cbind(1, mod1)

  B <- mod1x[, 0] ## accumulate batch-correction matrix

  ## --------------------------------------------------------------------
  ## Remove (unwanted) technical experiment effects (libsize, nfeature, etc.)
  ## --------------------------------------------------------------------
  if (lib.correct) {
    sel <- grep("libsize|nfeature", colnames(pheno), value = TRUE)
    if (length(sel)) {
      dbg("[pgx.superBatchCorrect] Correcting for unwanted library effects:", sel, "\n")
      exp.pheno <- as.matrix(pheno[, sel, drop = FALSE])
      exp.pheno <- apply(exp.pheno, 2, function(x) {
        x[is.na(x)] <- stats::median(x, na.rm = TRUE)
        x
      })
      cX <- limma::removeBatchEffect(cX, covariates = exp.pheno, design = mod1x)
      B <- cbind(B, exp.pheno)
    }
  }

  ## --------------------------------------------------------------------
  ## Remove (unwanted) biological effects
  ## --------------------------------------------------------------------
  if (!is.null(bio.correct) && length(bio.correct) > 0 && bio.correct[1] != FALSE) {
    p1 <- intersect(batch.prm, colnames(Y))
    dbg("[pgx.superBatchCorrect] Correcting for unwanted biological factors:", p1, "\n")
    if (length(p1)) {
      i <- 1
      for (i in 1:length(p1)) {
        b1 <- as.character(pheno[, p1[i]])
        b1[is.na(b1)] <- "NA" ## NA is third group?? better to impute??
        cX <- limma::removeBatchEffect(cX, batch = b1, design = mod1x)
        b1x <- stats::model.matrix(~b1)[, -1, drop = FALSE]
        colnames(b1x) <- sub("^b1", paste0(p1[i], "."), colnames(b1x))
        B <- cbind(B, b1x)
      }
    }

    p2 <- intersect(batch.cov, colnames(Y))
    if (length(p2)) {
      dbg("[pgx.superBatchCorrect] Correcting for unwanted biological covariates:", p2, "\n")
      b2 <- as.matrix(pheno[, p2, drop = FALSE])
      b2 <- apply(b2, 2, function(x) {
        x[is.na(x)] <- stats::median(x, na.rm = TRUE)
        x
      })
      cX <- limma::removeBatchEffect(cX, covariates = b2, design = mod1x)
      B <- cbind(B, b2)
    }

    ## out <- pgx.removeBiologicalEffect(cX, pheno, model.par=model.par,
    ##                                  correct=bio.correct, force=force)
  }

  ## --------------------------------------------------------------------
  ## batch correct other parameters with limma
  ## --------------------------------------------------------------------
  if (!is.null(batch.prm) && length(batch.prm) > 0) {
    batch.prm1 <- setdiff(batch.prm, colnames(Y))
    dbg("[pgx.superBatchCorrect] Batch correction for factors:", batch.prm1, "\n")
    b <- batch.prm1[1]
    for (b in batch.prm1) {
      batch <- as.character(pheno[, b])
      nna <- sum(is.na(batch))
      if (nna > 0) {
        ## impute missing values
        batch[is.na(batch)] <- sample(batch[!is.na(batch)], nna, replace = TRUE)
      }
      mod1x <- matrix(1, ncol(cX), 1)
      if (!is.null(mod1)) mod1x <- cbind(1, mod1)
      cX <- limma::removeBatchEffect(cX, batch = batch, design = mod1x)

      b1x <- stats::model.matrix(~batch)[, -1, drop = FALSE]
      colnames(b1x) <- sub("^batch", paste0(b, "."), colnames(b1x))
      B <- cbind(B, b1x)
    }
  }

  if (!is.null(batch.cov) && length(batch.cov) > 0) {
    batch.cov
    batch.cov1 <- setdiff(batch.cov, colnames(Y))
    dbg("[pgx.superBatchCorrect] Batch correction for covariates:", batch.cov1, "\n")
    for (b in batch.cov1) {
      batch <- as.numeric(pheno[, b])

      nna <- sum(is.na(batch))
      if (nna > 0) {
        batch[is.na(batch)] <- sample(batch[!is.na(batch)], nna, replace = TRUE)
      }
      mod1x <- matrix(1, ncol(cX), 1)
      if (!is.null(mod1)) mod1x <- cbind(1, mod1)
      cX <- limma::removeBatchEffect(cX, covariates = batch, design = mod1x)
      B <- cbind(B, batch)
    }
  }

  ## --------------------------------------------------------------------
  ## MNN correction (e.g. for single-cell)
  ## --------------------------------------------------------------------
  if (!is.null(mnn.correct)) {
    mnn.correct <- intersect(mnn.correct, colnames(pheno))
    if (length(mnn.correct) == 0) mnn.correct <- NULL
  }
  if (!is.null(mnn.correct)) {
    dbg("[pgx.superBatchCorrect] Mutual Nearest Neighbour (MNN) correction on", mnn.correct, "\n")
    b <- pheno[, mnn.correct]
    out <- mnnCorrect(cX, batch = b, cos.norm.out = FALSE)
    cX <- out@assays@data[["corrected"]]
  }

  ## --------------------------------------------------------------------
  ## Nearest-neighbour matching (NNM)
  ## --------------------------------------------------------------------

  if (1 && nnm.correct) {
    dbg("[pgx.superBatchCorrect] Correcting with nearest-neighbour matching (NNM)")
    dbg("[pgx.superBatchCorrect] NNM :: model.par = ", model.par)
    y1 <- pheno[, model.par, drop = FALSE]
    y1 <- apply(y1, 1, paste, collapse = ":")
    cX <- gx.nnmcorrect(cX, y1, center.x = TRUE, center.m = TRUE)$X
  }

  ## --------------------------------------------------------------------
  ## SVA correction (removing unwanted variation)
  ## --------------------------------------------------------------------
  if (sva.correct && !is.null(mod1)) {
    message("[pgx.superBatchCorrect] Calculating SVA...")

    ##
    ## This is a combination of methods from SVA and SmartSVA
    ## because of speed.
    ##

    mod1x <- cbind(1, mod1)
    mod0x <- mod1x[, 1, drop = FALSE] ## just ones...

    ## fast method using SmartSVA
    pp <- paste0(model.par, collapse = "+")
    lm.expr <- paste0("lm(t(cX) ~ ", pp, ", data=pheno)")
    X.r <- t(resid(eval(parse(text = lm.expr))))
    n.sv <- isva::EstDimRMT(X.r, FALSE)$dim + 1

    cX1 <- Matrix::head(cX[order(-apply(cX, 1, sd)), ], 1000) ## top 1000 genes only (faster)
    sv <- try(sva::sva(cX1, mod1x, mod0 = mod0x, n.sv = n.sv)$sv)

    if (any(class(sv) == "try-error")) {
      ## try again with little bit of noise...
      a <- 0.01 * mean(apply(cX, 1, sd))
      cX1 <- cX + a * matrix(rnorm(length(cX)), nrow(cX), ncol(cX))
      cX1 <- Matrix::head(cX1[order(-apply(cX1, 1, sd)), ], 1000) ## top 1000 genes only (faster)
      sv <- try(sva::sva(cX1, mod1x, mod0 = mod0x, n.sv = pmax(n.sv - 1, 1))$sv)
    }
    if (!any(class(sv) == "try-error")) {
      message("[pgx.superBatchCorrect] Performing SVA correction...")

      rownames(sv) <- colnames(cX)
      colnames(sv) <- paste0("SV.", 1:ncol(sv))
      cX <- limma::removeBatchEffect(cX, covariates = sv, design = mod1x)

      B <- cbind(B, sv)
    }
  }


  ## --------------------------------------------------------------------
  ## PCA correction: remove remaining batch effect using PCA
  ## (iteratively, only SV larger than max correlated SV)
  ## --------------------------------------------------------------------
  if (pca.correct && !is.null(mod1)) {
    ii <- 1:99
    niter <- 0
    nremove <- 0
    pX <- NULL
    while (length(ii) > 0 && niter < max.iter) {
      nv <- min(10, ncol(cX) - 1)
      suppressWarnings(suppressMessages(
        pc <- irlba::irlba(cX, nv = nv)$v
      ))
      pc.rho <- stats::cor(pc, mod1)
      pc.rho <- apply(abs(pc.rho), 1, max)
      ii <- which(pc.rho < max.rho)
      ii <- ii[ii < which.max(pc.rho)]
      if (length(ii) > 0) {
        mod1x <- cbind(1, mod1)
        cX <- limma::removeBatchEffect(cX, covariates = pc[, ii], design = mod1x)
        pX <- cbind(pX, pc[, ii, drop = FALSE])
        nremove <- nremove + 1
      }
      niter <- niter + 1
    }
    niter
    if (niter == max.iter) {
      dbg("WARNING:: PCA correction did not converge after", nremove, "iterations\n")
    } else {
      dbg("Performed", nremove, "iterations of PCA batch correction\n")
    }
    if (!is.null(pX)) {
      colnames(pX) <- paste0("PC.", 1:ncol(pX))
      B <- cbind(B, pX) ## update batch correction matrix
    }
  }

  ## --------------------------------------------------------------------
  ## HC correction: remove remaining batch effect iteratively using
  ## hclust
  ## --------------------------------------------------------------------
  if (hc.correct && !is.null(mod1)) {
    ii <- 1:99
    niter <- 0
    nremove <- 0
    pX <- NULL
    while (length(ii) > 0 && niter < max.iter) {
      xx <- Matrix::head(cX[order(-apply(cX, 1, sd)), ], hc.top)
      hc <- stats::cutree(fastcluster::hclust(dist(t(xx)), method = "ward.D2"), 2)
      hc.rho <- stats::cor(hc, mod1)
      hc.rho <- apply(abs(hc.rho), 1, max)
      ii <- which(hc.rho < max.rho)
      if (length(ii) > 0) {
        mod1x <- cbind(1, mod1)
        hc <- scale(hc)
        cX <- limma::removeBatchEffect(cX, covariates = hc, design = mod1x)
        pX <- cbind(pX, hc)
        nremove <- nremove + 1
      }
      niter <- niter + 1
    }
    if (niter == max.iter) {
      dbg("WARNING:: HC correction did not converge after", nremove, "iterations\n")
    } else {
      dbg("Performed", nremove, "iterations of HC batch correction\n")
    }
    if (!is.null(pX)) B <- cbind(B, pX) ## update batch correction matrix
  }

  ## --------------------------------------------------------------------
  ## important: means seems to be affected!!! regressed out??
  ## --------------------------------------------------------------------
  cX <- cX - rowMeans(cX) + rowMeans(X)

  dbg("[pgx.superBatchCorrect] almost done!\n")

  ## matrix B contains the active batch correction vectors

  res <- list(X = cX, Y = pheno, B = B)

  dbg("[pgx.superBatchCorrect] done!\n")

  return(res)
}


#' Correlate principal components with phenotypes
#'
#' @param X Expression matrix
#' @param pheno Data frame of sample phenotypes
#' @param nv Number of principal components to use
#' @param stat Statistic for categorical phenotypes ("F" or "t").
#' @param plot Logical, whether to generate a PCA plot.
#' @param main Title for PCA plot.
#'
#' @return Named vector of correlation coefficients
#'
#' @description
#' Calculates the correlation between principal components of the expression matrix X
#' and sample phenotypes in pheno.
#'
#' @details
#' This function calculates the top nv principal components of the expression matrix X
#' using irlba::irlba. It then correlates each PC with each phenotype in the pheno data frame.
#'
#' For categorical phenotypes, it converts to a factor and calculates the correlation with
#' the model matrix. For numeric phenotypes it calculates the standard correlation coefficient.
#'
#' It returns a named vector of correlation coefficients, with names corresponding to
#' the phenotypes.
#'
#' @export
pgx.PC_correlation <- function(X, pheno, nv = 3, stat = "F", plot = TRUE, main = NULL) {
  getF <- function(x, y) {
    x <- t(scale(t(x))) ## rowscale
    ii <- which(!is.na(y))
    y1 <- y[ii]
    if (inherits(y1, c("factor", "character", "logical"))) {
      y1 <- factor(as.character(y1))
    } else {
      y1 <- y1 + 1e-8 * stats::rnorm(length(y1))
      y1 <- (y1 > stats::median(y1))
    }
    design <- stats::model.matrix(~ 1 + y1)
    fit <- limma::lmFit(x[, ii], design)
    suppressWarnings(fit <- try(limma::eBayes(fit, trend = FALSE)))

    if (inherits(fit, "try-error")) {
      return(NULL)
    }
    suppressMessages(top <- limma::topTableF(fit, number = nrow(x)))

    return(top$F)
  }
  getCor <- function(x, y) {
    ii <- which(!is.na(y))
    y1 <- y[ii]
    if (inherits(y1, "factor")) y1 <- factor(as.character(y1))
    design <- stats::model.matrix(~ 0 + y1)

    r1 <- stats::cor(t(x[, ii]), design)
    rowMeans(abs(r1))
  }


  X <- X - rowMeans(X) ## center features
  V <- irlba::irlba(X, nv = nv)$v
  rho <- list()
  px <- pheno
  p <- "Chemotherapy"
  for (p in c("<random>", colnames(px))) {
    if (p == "<random>") {
      y <- sample(c("a", "b"), ncol(X), replace = TRUE)
    } else {
      y <- px[, p]
    }
    nlevels <- length(unique(y[!is.na(y)]))
    nrep <- max(table(y))
    if (nlevels > 1 && nrep >= 2) {
      if (stat == "cor") {
        rho[[p]] <- getCor(x = t(V), y)
      }
      if (stat == "F") {
        rho[[p]] <- getF(x = t(V), y)
      }
    }
  }

  R <- do.call(rbind, rho)
  colnames(R) <- paste0("PC", 1:ncol(R))
  if (stat == "F") R <- t(t(R) / colMeans(R))

  if (plot) {
    stat0 <- c("correlation", "F-statistic")[1 + 1 * (stat == "F")]
    tt0 <- c("PC correlation", "PC variation")[1 + 1 * (stat == "F")]
    if (is.null(main)) main <- tt0

    plt <- plot_ggbarplot(t(R), ylab = stat0, srt = 45, group.name = "") +
      ggplot2::theme(
        plot.margin = ggplot2::margin(2, 2, 0, 2, "mm"),
        plot.title = ggplot2::element_text(size = 12)
      ) +
      ggplot2::xlab("") + ggplot2::ggtitle(main)
    ## plt
    return(plt)
  }
  R
}


NORMALIZATION.METHODS <- c("none", "mean", "scale", "NC", "CPM", "TMM", "RLE", "RLE2", "quantile")

#' @title Normalize counts with TMM method
#'
#' @param counts Count matrix with genes in rows and samples in columns
#' @param log Logical, whether to log transform normalized counts
#' @param method Normalization method, default "TMM"
#'
#' @return Matrix of normalized counts
#'
#' @description
#' Normalizes a count matrix using the TMM (Trimmed Mean of M-values) method.
#'
#' @details
#' This function takes a count matrix and normalizes the counts using the TMM method
#' implemented in the edgeR package.
#'
#' It converts the count matrix to a DGEList object. Then it calculates normalization
#' factors using the calcNormFactors function with default TMM method.
#'
#' The normalized counts are returned as CPM values, optionally log-transformed.
#'
#' @examples
#' \dontrun{
#' counts <- matrix(rnbinom(100, mu = 10, size = 1), ncol = 10)
#' norm_counts <- normalizeTMM(counts)
#' }
#' @export
normalizeTMM <- function(counts, log = FALSE, method = "TMM") {
  dge <- edgeR::DGEList(as.matrix(counts), group = NULL)
  dge <- edgeR::calcNormFactors(dge, method = method)
  edgeR::cpm(dge, log = log)
}


#' Normalize count data using RLE method
#'
#' @description This function normalizes count data using the Relative Log Expression (RLE) method.
#' It supports normalization using either the edgeR or DESeq2 package.
#'
#' @param counts A numeric matrix of count data, where rows represent features and columns represent samples.
#' @param log Logical, whether to return log-transformed normalized counts.
#' @param use Character, the package to use for normalization, either "edger" or "deseq2".
#'
#' @return A numeric matrix of normalized count data.
#'
#' @examples
#' \dontrun{
#' counts <- matrix(rnbinom(100, mu = 10, size = 1), ncol = 10)
#' norm_counts <- normalizeRLE(counts)
#' }
#' @export
normalizeRLE <- function(counts, log = FALSE, use = "deseq2") {
  outx <- NULL
  if (use == "edger") {
    dge <- edgeR::DGEList(as.matrix(counts), group = NULL)
    dge <- edgeR::calcNormFactors(dge, method = "RLE")
    outx <- edgeR::cpm(dge, log = log)
  } else if (use == "deseq2") {
    cts <- counts
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = cts,
      colData = cbind(colnames(cts), 1),
      design = ~1
    )
    disp <- DESeq2::estimateSizeFactors(dds)
    outx <- DESeq2::counts(disp, normalized = TRUE)
  } else {
    stop("unknown method")
  }
  return(outx)
}


#' @title Normalize count data
#'
#' @description
#' Normalizes a matrix of RNA-seq count data using different methods.
#'
#' @param x Matrix of count data, with genes in rows and samples in columns.
#' @param methods Normalization method(s) to use. Options are "none", "scale", "CPM", "TMM", "RLE".
#' @param keep.zero Logical indicating whether to retain zero counts. Default is TRUE.
#'
#' @details
#' This function normalizes a matrix of RNA-seq count data using different normalization methods:
#'
#' - "none": No normalization
#' - "scale": Scale normalization by column means
#' - "CPM": Counts per million
#' - "TMM": Trimmed mean of M-values (edgeR)
#' - "RLE": Relative log expression (DESeq2)
#'
#' Zero counts can be retained or set to NA after normalization based on the keep.zero parameter.
#'
#' @return
#' The normalized count matrix.
#'
#' @examples
#' \dontrun{
#' counts <- matrix(rnbinom(10000, mu = 10, size = 1), 100, 100)
#' normalized <- pgx.countNormalization(counts, c("TMM", "RLE"))
#' }
#' @export
pgx.countNormalization <- function(x, methods, keep.zero = TRUE) {
  ## Column-wise normalization (along samples).
  ##
  ## x:        counts (linear)
  ## method:   single method

  methods <- methods[1]
  which.zero <- which(x == 0, arr.ind = TRUE)

  for (m in methods) {
    if (m == "none") {
      ## normalization on individual mean
      x <- x
    } else if (m == "scale") {
      ## normalization on individual mean
      mx <- mean(x, na.rm = TRUE)
      x <- t(t(x) / colMeans(x, na.rm = TRUE)) * mx
    } else if (m == "CPM") {
      x <- t(t(x) / Matrix::colSums(x, na.rm = TRUE)) * 1e6
    } else if (m == "TMM") {
      ## normalization on total counts (linear scale)
      x <- normalizeTMM(x, log = FALSE) ## does TMM on counts (edgeR)
    } else if (m == "RLE") {
      ## normalization on total counts (linear scale)
      x <- normalizeRLE(x, log = FALSE, use = "deseq2") ## does RLE on counts (Deseq2)
    } else if (m == "RLE2") {
      ## normalization on total counts (linear scale)
      x <- normalizeRLE(x, log = FALSE, use = "edger") ## does RLE on counts (Deseq2)
      ##        } else if(m %in% c("upperquartile")) {
      ##            ## normalization on total counts (linear scale)
    } else if (m == "quantile") {
      new.x <- 0.01 * limma::normalizeQuantiles(as.matrix(100 * x)) ## shift to avoid clipping
      rownames(new.x) <- rownames(x)
      colnames(new.x) <- colnames(x)
      x <- new.x
    }
  } ## end of for method

  x <- pmax(x, 0) ## prevent negative values
  ## put back zeros as zeros
  if (keep.zero && nrow(which.zero) > 0) {
    x[which.zero] <- 0
  }

  return(x)
}








#' @title Estimate biological variation
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
#' results <- pgx.computeBiologicalEffects(sample.ExpressionSet)
#' head(results$biological)
#' }
#' @export
pgx.computeBiologicalEffects <- function(X, is.count = FALSE) {
  ## estimate biological variation
  ##
  ## X:     log-expression matrix
  ##

  message("[pgx.computeBiologicalEffects] estimating biological effects...")

  ## shift zero to 1% percentile
  if (!is.count) {
    q0 <- stats::quantile(X[X > 0], probs = 0.01, na.rm = TRUE)

    tx <- pmax(X - q0, 0) ## log expression
    cx <- pmax(2**tx - 1, 0) ## counts
  } else {
    cx <- X
    tx <- log2(cx + 1)
  }
  nfeature <- Matrix::colSums(cx > 0) + 1
  libsize <- Matrix::colSums(cx)

  mt.genes <- grep("^MT-", rownames(X), ignore.case = TRUE, value = TRUE)
  rb.genes <- grep("^RP[SL]", rownames(X), ignore.case = TRUE, value = TRUE)
  mito <- ribo <- NA
  pct.mito <- pct.ribo <- NA

  if (length(mt.genes) >= 10) {
    mito <- Matrix::colMeans(tx[mt.genes, , drop = FALSE])
    pct.mito <- Matrix::colSums(cx[mt.genes, , drop = FALSE], na.rm = TRUE) / libsize
  }
  if (length(rb.genes) >= 10) {
    ii <- rb.genes[order(-apply(tx[rb.genes, , drop = FALSE], 1, sd, na.rm = TRUE))]
    sel20 <- Matrix::head(ii, 20)
    ribo <- Matrix::colMeans(tx[rb.genes, , drop = FALSE])
    ribo20 <- Matrix::colMeans(tx[sel20, , drop = FALSE])
    pct.ribo <- Matrix::colSums(cx[rb.genes, , drop = FALSE], na.rm = TRUE) / libsize
  }
  pheno <- data.frame(
    mito = mito,
    ribo = ribo,
    libsize = log2(libsize + 1),
    check.names = FALSE
  )

  cc.score <- try(pgx.scoreCellCycle(cx))
  if (!any(class(cc.score) == "try-error")) {
    cc.score <- cc.score[, c("s_score", "g2m_score")]
    colnames(cc.score) <- paste0("cc.", colnames(cc.score))
    pheno <- cbind(pheno, cc.score)
  }
  pheno$gender <- pgx.inferGender(cx)

  return(pheno)
}




## ================================================================================
## Automagic batch correction by trying all combinations of batch parameters and
## optimizing the number of significant genes.
## ================================================================================




#' @title Compute number of significant genes
#'
#' @param ngs Object containing NGS analysis data
#' @param X Expression matrix
#' @param contrast Contrast to test
#' @param fc Fold change cutoff
#' @param qv FDR cutoff
#'
#' @return Named vector containing number of significant genes for each parameter combination
#'
#' @description Counts number of significant genes for given parameters.
#'
#' @details This function fits the specified contrast and counts the number of genes
#' significant at the given FDR and fold change cutoffs.
#'
#' It requires an NGS analysis object, expression matrix, and contrast. Batch
#' correction is performed using parameters defined in the NGS object. The contrast is
#' then fit using limma-voom and number of significant genes counted.
#'
#' Fold change and FDR cutoffs can also be specified. The number of significant genes
#' passing these thresholds is returned for each parameter combination.
#'
#' @export
pgx._computeNumSig <- function(ngs, X, contrast = NULL, fc = 0, qv = 0.05) {
  samples <- colnames(X)
  design <- ngs$model.parameters$design[samples, ]
  contr.matrix <- ngs$model.parameters$contr.matrix
  if (is.null(contrast)) contrast <- colnames(contr.matrix)
  contr.matrix <- contr.matrix[, contrast, drop = FALSE]
  res <- ngs.fitContrastsWithLIMMA(
    X, contr.matrix, design,
    method = "limma", trend = TRUE,
    conform.output = FALSE, plot = FALSE
  )

  fc0 <- sapply(res$tables, function(x) x$logFC)
  qv0 <- sapply(res$tables, function(x) x$adj.P.Val)
  numsig <- mean(Matrix::colSums(abs(fc0) >= fc & qv0 <= qv, na.rm = TRUE))

  return(numsig)
}

## =====================================================================================
## =========================== END OF FILE =============================================
## =====================================================================================
