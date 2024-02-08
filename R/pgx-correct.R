##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

## Batch correction methods
##
##


#' @title Supervised batch correction
#'
#' @description Performs supervised batch correction on a gene
#'   expression matrix, using known technical factors and biological
#'   covariates.
#'
#' @param X Gene expression matrix, genes in rows, samples in columns.
#' @param pheno Dataframe containing sample metadata/covariates. Must match colnames of \code{X}.
#' @param model.par Vector of column names of biological covariates in \code{pheno}.
#' @param batch.par Vector of column names of batch covariates in \code{pheno}.
#'
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
pgx.superBatchCorrect <- function(X, pheno,
                                  model.par,
                                  partype = NULL,
                                  batch.par = "*",
                                  lib.correct = TRUE,
                                  bio.correct = c("mito", "ribo", "cell_cycle", "gender"),
                                  sva.correct = TRUE,
                                  pca.correct = TRUE,
                                  hc.correct = TRUE,
                                  mnn.correct = NULL,
                                  nnm.correct = TRUE,
                                  bc.methods = NULL,
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

  ## override old-style
  if (!is.null(bc.methods)) {
    bc.methods <- tolower(bc.methods)
  } else {
    bc <- c()
    if (!is.null(mnn.correct)) bc <- c(bc, "mnn")
    if (nnm.correct) bc <- c(bc, "nnm")
    if (sva.correct) bc <- c(bc, "sva")
    if (pca.correct) bc <- c(bc, "pca")
    if (hc.correct) bc <- c(bc, "hc")
    bc.methods <- bc
  }

  ## tidy up pheno matrix?? get correct parameter types
  pheno <- utils::type.convert(pheno, as.is = TRUE)

  ## setup model matrix
  mod1 <- NULL
  if (!is.null(model.par) && length(model.par) > 0) {
    model.par <- intersect(model.par, colnames(pheno))
    mod1 <- do.call(cbind, lapply(model.par, getModelMatrix))
    rownames(mod1) <- rownames(pheno)
  }
  model.par

  ## get technical/biological effects
  Y <- pgx.computeBiologicalEffects.DEPRECATED(X)
  colnames(Y) <- paste0(".", colnames(Y))

  ## add to phenotype matrix
  pheno <- cbind(pheno, Y)
  not.na <- colMeans(is.na(pheno)) < 1
  nlev <- apply(pheno, 2, function(x) length(unique(x[!is.na(x)])))
  pheno <- pheno[, which(nlev > 1 & not.na), drop = FALSE]
  partype <- sapply(pheno, class)

  ## --------------------------------------------------------------------
  ## select parameters
  ## --------------------------------------------------------------------

  ## select all non-model variables
  if (!is.null(batch.par) && batch.par[1] == "*") {
    batch.par <- setdiff(colnames(pheno), model.par)
  }

  if ("mito" %in% bio.correct) {
    b1 <- grep("^mito", colnames(pheno), value = TRUE)
    batch.par <- c(batch.par, b1)
  }
  if ("ribo" %in% bio.correct) {
    b1 <- grep("^ribo", colnames(pheno), value = TRUE)
    batch.par <- c(batch.par, b1)
  }
  if ("cell_cycle" %in% bio.correct) {
    b1 <- grep("^cc[.]|^cellcycle", colnames(pheno), value = TRUE) ## only s.score and g2m.score
    batch.par <- c(batch.par, b1)
  }
  if ("gender" %in% bio.correct) {
    b1 <- grep("^gender", colnames(pheno), value = TRUE)
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

  ## iterate over bc.methods in this order
  for (bc in bc.methods) {
    if (bc == "mnn") {
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
    }

    if (bc == "nnm") {
      ## --------------------------------------------------------------------
      ## Nearest-neighbour matching (NNM)
      ## --------------------------------------------------------------------
      dbg("[pgx.superBatchCorrect] Correcting with nearest-neighbour matching (NNM)")
      dbg("[pgx.superBatchCorrect] model.par = ", model.par)
      y1 <- pheno[, model.par, drop = FALSE]
      y1 <- apply(y1, 1, paste, collapse = ":")
      cX <- gx.nnmcorrect(cX, y1, center.x = TRUE, center.m = TRUE)$X
    }

    if (bc == "sva") {
      ## --------------------------------------------------------------------
      ## SVA correction (removing unwanted variation)
      ## --------------------------------------------------------------------
      if (!is.null(mod1)) {
        message("[pgx.superBatchCorrect] Calculating SVA...")
        ##
        ## This is a combination of methods from SVA and SmartSVA
        ## because of speed.
        mod1x <- cbind(1, mod1)
        mod0x <- mod1x[, 1, drop = FALSE] ## just ones...

        ## fast method using SmartSVA
        pp <- paste0(model.par, collapse = "+")
        lm.expr <- paste0("lm(t(cX) ~ ", pp, ", data=pheno)")
        X.r <- t(stats::resid(eval(parse(text = lm.expr))))
        n.sv <- isva::EstDimRMT(X.r, FALSE)$dim + 1

        cX1 <- Matrix::head(cX[order(-apply(cX, 1, stats::sd)), ], 1000) ## top 1000 genes only (faster)
        sv <- try(sva::sva(cX1, mod1x, mod0 = mod0x, n.sv = n.sv)$sv)

        if (any(class(sv) == "try-error")) {
          ## try again with little bit of noise...
          a <- 0.01 * mean(apply(cX, 1, stats::sd, na.rm = TRUE), na.rm = TRUE)
          cX1 <- cX + a * matrix(stats::rnorm(length(cX)), nrow(cX), ncol(cX))
          cX1 <- Matrix::head(cX1[order(-apply(cX1, 1, stats::sd)), ], 1000) ## top 1000 genes only (faster)
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
    }


    if (bc == "pca") {
      ## --------------------------------------------------------------------
      ## PCA correction: remove remaining batch effect using PCA
      ## (iteratively, only SV larger than max correlated SV)
      ## --------------------------------------------------------------------
      if (!is.null(mod1)) {
        ii <- 1:99
        niter <- 0
        nremoved <- 0
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
            nremoved <- nremoved + 1
          }
          niter <- niter + 1
        }
        niter
        if (niter == max.iter) {
          dbg("WARNING:: PCA correction did not converge after", nremoved, "iterations\n")
        } else {
          dbg("PCA batch correction: removed", nremoved, "principal components\n")
        }
        if (!is.null(pX)) {
          colnames(pX) <- paste0("PC.", 1:ncol(pX))
          B <- cbind(B, pX) ## update batch correction matrix
        }
      }
    }

    if (bc == "hc") {
      ## --------------------------------------------------------------------
      ## HC correction: remove remaining batch effect iteratively using
      ## hclust
      ## --------------------------------------------------------------------
      if (!is.null(mod1)) {
        ii <- 1:99
        niter <- 0
        nremoved <- 0
        pX <- NULL
        while (length(ii) > 0 && niter < max.iter) {
          xx <- Matrix::head(cX[order(-apply(cX, 1, stats::sd)), ], hc.top)
          hc <- stats::cutree(fastcluster::hclust(stats::dist(t(xx)), method = "ward.D2"), 2)
          hc.rho <- stats::cor(hc, mod1)
          hc.rho <- apply(abs(hc.rho), 1, max)
          ii <- which(hc.rho < max.rho)
          if (length(ii) > 0) {
            mod1x <- cbind(1, mod1)
            hc <- scale(hc)
            cX <- limma::removeBatchEffect(cX, covariates = hc, design = mod1x)
            pX <- cbind(pX, hc)
            nremoved <- nremoved + 1
          }
          niter <- niter + 1
        }
        if (niter == max.iter) {
          dbg("WARNING:: HC correction did not converge after", nremoved, "iterations\n")
        } else {
          dbg("Performed", nremoved, "iterations of HC batch correction\n")
        }
        if (!is.null(pX)) B <- cbind(B, pX) ## update batch correction matrix
      }
    }
  } ## end of for bc.methods

  ## --------------------------------------------------------------------
  ## important: means seems to be affected!!! regressed out??
  ## --------------------------------------------------------------------
  cX <- cX - rowMeans(cX, na.rm = TRUE) + rowMeans(X, na.rm = TRUE)

  ## matrix B contains the active batch correction vectors
  res <- list(X = cX, Y = pheno, B = B)

  dbg("[pgx.superBatchCorrect] done!\n")

  return(res)
}


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
    r1 <- rowMeans(abs(r1))
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
pgx.computeBiologicalEffects.DEPRECATED <- function(X, is.count = FALSE) {
  message("[pgx.computeBiologicalEffects] estimating biological effects...")

  ## shift zero to 1% percentile
  if (!is.count) {
    q0 <- quantile(X[X > 0], probs = 0.01, na.rm = TRUE)
    q0
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
  mt.genes
  rb.genes
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
    ## ribo20 = ribo20,
    ## pct.mito = pct.mito,
    ## pct.ribo = pct.ribo,
    libsize = log2(libsize + 1),
    ## nfeature = log2(nfeature+1),
    check.names = FALSE
  )

  cc.score <- try(pgx.scoreCellCycle(cx))
  Matrix::head(cc.score)
  if (!any(class(cc.score) == "try-error")) {
    ## cc.score <- cc.score[,c("s_score","g2m_score","diff_score")]
    cc.score <- cc.score[, c("s_score", "g2m_score")]
    colnames(cc.score) <- paste0("cc.", colnames(cc.score))
    pheno <- cbind(pheno, cc.score)
  }
  pheno$gender <- pgx.inferGender(cx)

  Matrix::head(pheno)
  return(pheno)
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
pgx.computeTechnicalEffects <- function(X, is.count = FALSE, nmin = 3, nv = 2) {
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
    mt.genes <- mt.genes[order(-apply(X[mt.genes, , drop = FALSE], 1, stats::sd, na.rm = TRUE))]
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
    rb.genes <- rb.genes[order(-apply(X[rb.genes, , drop = FALSE], 1, stats::sd, na.rm = TRUE))]
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
    if (nv == 1) cc.score <- rowMeans(cc.score)
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
                               k.pca = 10, nv = 2, xrank = NULL) {
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
  } else if (length(pheno) == nrow(samples)) {
  } else {
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
bc.plotCovariateHeatmap <- function(bc.res) {
  ## bc <- detectBatchEffects(X, samples, pheno, contrasts = NULL,
  ##                          params = c("statistical", "technical", "pca"),
  ##                          p.pca = 0.5, p.pheno = 0.05,
  ##                          k.pca = 10, nv = 2, xrank = NULL)
  B <- bc.res$covariates_plus
  rho <- cor(apply(B, 2, rank))
  colnames(rho) <- rep("", ncol(rho))
  gx.heatmap(rho,
    sym = TRUE, mar = c(1, 15), keysize = 0.4, cexCol = 0.0001,
    scale = "none", key = FALSE
  )
}

## p.pheno=0.05;p.pca=0.5;nmin=3;nv=2
## params = c("lib","gender","mito","ribo","cellcycle")

#' @export
removeTechnicalEffects <- function(X, samples, y, p.pheno = 0.05, p.pca = 0.5,
                                   params = c("lib", "mito", "ribo", "cellcycle", "gender"),
                                   force = FALSE, nv = 2, k.pca = 10, xrank = NULL) {
  ##  p.pheno = 0.05;p.pca = 0.5;force = FALSE; nv = 2;k.pca = 10;xrank = NULL
  ##  params = c("lib","mito","ribo","cellcycle","gender")

  X1 <- X
  X1 <- svdImpute2(X1) ## temporary hack. need to refactor: allowing for real NA!!
  if (force) {
    bc <- detectBatchEffects(X1, samples, y,
      params = "technical",
      p.pca = 1, p.pheno = 0, k.pca = k.pca, nv = nv, xrank = 999
    )
  } else {
    bc <- detectBatchEffects(X1, samples, y,
      params = "technical",
      p.pca = p.pca, p.pheno = p.pheno, k.pca = k.pca,
      nv = nv, xrank = xrank
    )
  }
  bc$params

  message("[removeTechnicalEffect] params = ", paste(params, collapse = " "))
  message("[removeTechnicalEffect] length(bc.params) = ", length(bc$params))
  message("[removeTechnicalEffect] bc.params = ", paste(unlist(bc$params), collapse = " "))

  if (!is.null(params)) {
    B <- bc$covariates
    sel <- lapply(params, function(p) grep(paste0("^", p, "[.]"), colnames(B)))
    sel <- sort(unique(unlist(sel)))
    bc$covariates <- B[, sel, drop = FALSE]
    dim(bc$covariates)
  }

  if (!is.null(bc$covariates)) {
    ## perform regression
    B <- scale(bc$covariates)
    B[is.nan(B) | is.na(B)] <- 0
    B[is.infinite(B)] <- 0
    bb <- paste(colnames(B), collapse = "+")
    message("[removeTechnicalEffect] correcting for parameters: ", bb)
    design <- model.matrix(~y)
    bX <- limma::removeBatchEffect(X1, batch = NULL, covariates = B, design = design)
  } else {
    message("[removeTechnicalEffect] no significant technical batch effects. correction skipped.")
    bX <- X
  }

  ## put back missing values
  bX[is.na(X)] <- NA

  ## put back on original row means
  bX <- bX - rowMeans(bX, na.rm = TRUE) + rowMeans(X, na.rm = TRUE)
  bX
}

## ================================================================================
## Run/compare multiple batch-correction methods
## ================================================================================


#' @export
runBatchCorrectionMethods <- function(X, batch, y, controls = NULL, ntop = 2000,
                                      combatx = FALSE, sc = FALSE, prefix = "",
                                      methods = NULL, remove.failed = TRUE) {
  if (0) {
    controls <- NULL
    ntop <- 2000
    combatx <- FALSE
    sc <- FALSE
    prefix <- ""
    methods <- NULL
    remove.failed <- TRUE
  }

  mod <- model.matrix(~y)
  nlevel <- length(unique(y[!is.na(y)]))
  if (ntop < Inf) {
    X <- head(X[order(-matrixStats::rowSds(X, na.rm = TRUE)), ], ntop) ## faster
  }

  if (is.null(methods)) {
    methods <- c(
      "uncorrected", "normalized_to_control", "ComBat",
      "limma", "superBC", "PCA", "RUV", "SVA", "NNM", "MNN", "Harmony"
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
  if ("limma" %in% methods) {
    if (is.null(batch)) {
      xlist[["limma"]] <- X
    } else {
      cX <- try(limma::removeBatchEffect(X,
        batch = batch, covariates = NULL, mod = mod1
      ))
      xlist[["limma"]] <- cX
      cX <- try(limma::removeBatchEffect(X,
        batch = batch, covariates = NULL
      ))
      xlist[["limma.no_mod"]] <- cX
    }
  }

  ## ComBat ------------------------------------------------------
  if ("ComBat" %in% methods) {
    if (is.null(batch)) {
      xlist[["ComBat"]] <- X
    } else {
      if (max(table(batch)) > 1) {
        mod1 <- model.matrix(~ factor(y))
        bX <- try(sva::ComBat(X, batch = batch, mod = mod1, par.prior = TRUE))
        xlist[["ComBat"]] <- bX
        bX <- try(sva::ComBat(X, batch = batch, mod = NULL, par.prior = TRUE))
        xlist[["ComBat.no_mod"]] <- bX
      }
    }

    ## ComBatX
    if (combatx && !is.null(batch)) {
      xlist[["ComBatX.a"]] <- ComBatX(X, batch = batch, y = NULL, controls = NULL)
      xlist[["ComBatX.b"]] <- ComBatX(X, batch = batch, y = y, controls = NULL, bc.dim = 0)
      xlist[["ComBatX.c"]] <- ComBatX(X, batch = batch, y = y, controls = NULL, bc.dim = -1)
      ##  xlist[['ComBatX.d']] <- ComBatX(X, batch=batch, y=y, controls=NULL, bc.dim=10)
      xlist[["ComBatX.e"]] <- ComBatX(X, batch = batch, y = NULL, controls = controls)
      xlist[["ComBatX.f"]] <- ComBatX(X, batch = batch, y = y, controls = controls, bc.dim = 0)
      xlist[["ComBatX.g"]] <- ComBatX(X, batch = batch, y = y, controls = controls, bc.dim = -1)
    }
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
    message("[runBatchCorrectionMethods] correcting with PCA")
    ##  xlist[["PCA"]] <- try(pcaCorrect.OLD(X, y = y, max.rho=0.3))
    xlist[["PCA"]] <- try(pcaCorrect2(X, y = y, p.notsig = 0.20))
  }

  ## RUV and SVA
  if ("RUV" %in% methods) {
    message("[runBatchCorrectionMethods] correcting with RUV3")
    xlist[["RUV3"]] <- try(ruvCorrect(X, y, k = NULL, type = "III"))
    ##  xlist[["RUVg"]] <- try(ruvCorrect(X, y, k = NULL, type = "g"))
  }

  if ("SVA" %in% methods) {
    message("[runBatchCorrectionMethods] correcting with SVA")
    xlist[["SVA"]] <- try(svaCorrect(X, y))
  }

  if ("NNM" %in% methods) {
    ## xlist[["NNM"]] <- gx.nnmcorrect(X, y)$X
    xlist[["NNM1"]] <- nnmCorrect(X, y, use.design = TRUE)
    xlist[["NNM2"]] <- nnmCorrect2(X, y, use.design = TRUE)      
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
    xlist <- xlist[!is.null(unlist(sapply(xlist, nrow)))]
  }

  names(xlist) <- paste0(prefix, names(xlist))
  xlist
}

#' @export
runTechCorrectionMethods <- function(X, samples, y, p.pca = 0.5, p.pheno = 0.05, nv = 2,
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
  message("computing statistics...")
  numsig <- lapply(xlist, stats.numsig,
    y = pheno, lfc = lfc, q = q,
    trend = trend, verbose = FALSE
  )

  res <- t(sapply(numsig, function(r) {
    c(sapply(r[1:2], length), avg.fc = mean(abs(r[[3]])))
  }))
  sdx <- sapply(xlist, function(x) mean(matrixStats::rowSds(x, na.rm = TRUE)))
  snr <- res[, "avg.fc"] / sdx
  res <- cbind(res, avg.sd = sdx, SNR = snr)

  ## compute relative genes/geneset overlap
  message("computing overlap...")
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
    (x - rowMeans(x))
  })

  message("computing silhouette scores...")
  silhouette <- rep(1, nrow(res))
  if (add.sil) {
    if (is.null(pos)) {
      if (clust == "tsne") {
        nb <- max(1, min(30, round(ncol(xlist[[1]]) / 5)))
        ## CLUSTFUN <- function(x) uwot::tumap(scale(t(x), scale = FALSE), n_neighbors = nb)
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
      mean(score[, "sil_width"])
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
    rho <- lapply(rho, function(x) rowMeans(abs(x)))
    pc1.ratio <- sapply(rho, function(r) abs(r[1]) / sum(abs(r)))

    res <- cbind(res, silhouette, pc1.ratio)
  }

  ## use only these for score
  sel <- c("genes", "gsets", "SNR", "pc1.ratio", "silhouette")
  sel <- intersect(sel, colnames(res))

  ##  score <- res.score * (silhouette / silhouette[1])**1
  overall.score <- t(t(1e-4 + res[, sel]) / (1e-4 + res[ref, sel]))
  overall.score[, "silhouette"] <- overall.score[, "silhouette"]**2 ## give more weight
  overall.score <- exp(rowMeans(log(overall.score))) ## geometric mean

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
bc.plotResults <- function(X, xlist, pos, pheno, samples = NULL, scores = NULL,
                           type = "umap", nmax = 1000, cex = 1, text.cex = 1,
                           ncol = NULL, par = TRUE) {
  ## samples=NULL;scores = NULL;type='umap';nmax=1000;cex=1;text.cex = 1
  if (par) {
    if (is.null(ncol)) {
      ncol <- ceiling(sqrt(length(xlist)))
    }
    nr <- ceiling(length(xlist) / ncol)
    par(mfrow = c(nr, ncol))
  }

  methods <- names(xlist)
  if (!is.null(pos)) {
    methods <- intersect(methods, names(pos))
  }

  if (!is.null(scores)) {
    methods <- intersect(methods, rownames(scores))
    m.score <- scores[methods, "score"]
    methods <- methods[order(-m.score)]
    scores <- scores[methods, ]
    xlist <- xlist[methods]
    pos <- pos[methods]
  }

  cex1 <- ifelse(length(pheno) > 20, 3, 4)
  cex1 <- ifelse(length(pheno) > 100, 2.5, cex1)
  cex1 <- ifelse(length(pheno) > 400, 2, cex1)
  cex1 <- ifelse(length(pheno) > 1000, 1, cex1)
  cex1 <- 0.7 * cex * cex1

  if (tolower(type) == "umap") {
    par(mar = c(2.4, 3, 2, 1))
    i <- 1
    for (m in methods) {
      plot(pos[[m]], col = factor(pheno), pch = 20, cex = cex1)
      title(main = m, cex.main = 1.4, line = 0.5)
      ## tt <- paste("score = ",round(res$results[m,"score"],2))
      ## legend("topleft", legend=tt, cex=1.4)
    }
  }

  if (type == "heatmap") {
    par(mar = c(2, 3, 1.8, 1))
    i <- 1
    for (m in methods) {
      xx <- xlist[[m]]
      xx <- head(xx[order(-apply(xx, 1, sd)), ], nmax)
      xx <- xx - rowMeans(xx)
      xx <- abs(xx)**0.5 * sign(xx)
      gx.imagemap(xx, main = m, cex.main = 1.4, cex = 0)
      mtext("samples", 1, line = 0.5, las = 1)
      mtext("genes", 2, line = 0.5, las = 3)
    }
  }

  if (tolower(type) %in% c("pc", "pc2")) {
    if (is.null(samples)) message("samples in NULL!")
    B <- pgx.computeTechnicalEffects(X, nv = 2)
    bcat <- sub("[.].*", "", colnames(B))
    colnames(B) <- paste0(bcat, ":", colnames(B)) ## for collapsing
    if (!is.null(samples)) B <- cbind(B, samples)
    horiz <- ifelse(tolower(type) == "pc2", TRUE, FALSE)

    plist <- list()
    i <- 1
    for (m in methods) {
      xx <- xlist[[m]]
      plist[[m]] <- pgx.PC_correlation(xx, B,
        nv = 3, stat = "F",
        plot = TRUE, main = m, expand = FALSE, collapse = TRUE,
        horiz = horiz, text.cex = text.cex
      )
    }

    gridExtra::grid.arrange(grobs = plist, ncol = ncol, padding = unit(0.0, "line"))
  }

  if (type == "hist") {
    par(mar = c(3, 3, 4, 2))
    i <- 1
    for (m in methods) {
      xx <- xlist[[m]]
      hist(xx, breaks = 100, main = m, cex.main = 1.8)
    }
  }

  if (type == "scores" && !is.null(scores)) {
    m <- 1
    plt <- list()
    ylabs <- c(
      "score" = "overall score",
      "genes" = "significant genes",
      "gsets" = "significant genesets",
      "avg.fc" = "average abs.logFC",
      "avg.sd" = "average SD",
      "r.genes" = "gene.coverage",
      "r.gsets" = "gset coverage",
      "SNR" = "signal-to-noise",
      "pc1.ratio" = "PC1 ratio",
      "silhouette" = "silhoutte score"
    )
    for (i in 1:ncol(scores)) {
      nn <- colnames(scores)[i]
      xx <- scores[, i]
      plt[[i]] <- plot_ggbarplot(
        t(xx),
        ylab = ylabs[nn],
        srt = 45,
        legend.cex = 1.2 * text.cex,
        label.cex = 1.2 * text.cex,
        axis.cex = 1.2 * text.cex,
        group.name = ""
      ) +
        ggplot2::theme(
          plot.margin = ggplot2::margin(t = 0, r = 4, b = 0, l = 8, "pt"),
          plot.title = ggplot2::element_text(size = 15 * text.cex)
        ) +
        ggplot2::xlab("") + ggplot2::ggtitle(nn)
    }

    gridExtra::grid.arrange(grobs = plt, ncol = ncol, padding = unit(0.0, "line"))
  }
}

#' @export
bc.CovariateAnalysisPlot <- function(bc.results, k = 1:3, par = TRUE, col = 1) {
  bc <- bc.results
  if (par == TRUE) par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))
  pp <- intersect(rownames(bc$p.pca), rownames(bc$p.values))
  pp <- grep("^pca", pp, value = TRUE, invert = TRUE)
  k <- k[which(k <= ncol(bc$p.pca))]
  pxx <- bc$p.pca[pp, k, drop = FALSE]
  py <- bc$p.value[pp, 2]
  x1 <- -log10(1e-04 + py)
  col <- rep(col, length(k))
  for (i in c(0, k)) {
    if (i == 0) {
      plot.new()
      abline(h = 0.5, lty = 2)
      abline(v = 0.5, lty = 2)
      mtext("correlation with PC   →", side = 2, line = 1.3, cex = 0.85)
      mtext("correlation with phenotype   →", side = 1, line = 1.3, cex = 0.85)
      ##      axis(side=1, tick='n', cex.axis=0.001)
      ##      axis(side=2, tick='n', cex.axis=0.001)
      text(
        x = 0.2, y = 0.80, adj = 0.5,
        labels = "strong batch-effects\nor\nstratification factors"
      )
      text(
        x = 0.75, y = 0.80, adj = 0.5,
        labels = "well designed model-parameters\nor\nstrong confouders"
      )
      text(
        x = 0.2, y = 0.20, adj = 0.5,
        labels = 'nuisance parameters\nor\n"noise"'
      )
      text(
        x = 0.75, y = 0.20, adj = 0.5,
        labels = "weak model-parameters\nor\nweak confouders"
      )
    } else {
      y1 <- -log10(1e-04 + pxx[, i])
      ylim <- c(-0.1 * max(y1), 1.1 * max(y1))
      col1 <- col[i]
      plot(x1, y1,
        pch = 20, cex = 1.5, col = col1,
        xlab = "significance with phenotype (-log10p)",
        ylab = "significance with PC (-log10p)",
        xlim = c(-0.4, 4.4), ylim = ylim
      )
      title(paste0("PC", i), cex.main = 1.4)
      text(x1, y1, pp, pos = c(1:4), cex = 1.3, col = col1)
    }
  }
}

#' @export
get_model_parameters <- function(X, samples, pheno = NULL, contrasts = NULL) {
  if (is.null(pheno) && is.null(contrasts)) {
    stop("must give either pheno vector or contrasts matrix")
  }
  if (!is.null(contrasts)) {
    pheno <- playbase::contrasts2pheno(contrasts, samples)
  }

  bc <- playbase::detectBatchEffects(X, samples, pheno,
    params = "statistical",
    k.pca = 10, p.pca = 0.5, p.pheno = 0.05, xrank = 10
  )

  bc$p.values

  ## pheno.pars <- names(which(bc$p.values[,'p.pheno'] < 1e-20))
  p.pheno <- bc$p.values[, "p.pheno"]
  if (nrow(bc$p.values) == 1) names(p.pheno)[1] <- rownames(bc$p.values)
  p.pheno
  pheno.pars <- names(which(p.pheno == min(p.pheno) | p.pheno < 1e-80))
  pheno.pars <- pheno.pars[order(p.pheno[pheno.pars])]
  pheno.pars

  batch.pars <- NULL
  if ("statistical" %in% names(bc$params)) {
    batch.pars <- bc$params$statistical
    batch.pars
    batch.pars <- setdiff(batch.pars, pheno.pars)
    batch.pars2 <- grep("batch", colnames(samples), ignore.case = TRUE, value = TRUE)
    if (length(batch.pars2)) batch.pars <- c(batch.pars, batch.pars2)
    batch.pars <- sort(batch.pars)
  }

  batch.vec <- NULL
  if (length(batch.pars)) {
    batch.vec <- apply(samples[, batch.pars, drop = FALSE], 1, paste, collapse = "_")
  }
  pheno.vec <- NULL
  if (length(pheno.pars)) {
    pheno.vec <- apply(samples[, pheno.pars, drop = FALSE], 1, paste, collapse = "_")
  }

  list(
    batch.pars = batch.pars,
    pheno.pars = pheno.pars,
    batch = batch.vec,
    pheno = pheno.vec
  )
}


#' @export
compare_batchcorrection_methods <- function(X, samples, pheno, contrasts,
                                            methods = c(
                                              "uncorrected", "ComBat",
                                              "limma", "RUV", "SVA", "NNM"
                                            ),
                                            ntop = 4000, xlist.init = list(),
                                            ref = NULL) {
  ## methods <- c("uncorrected","ComBat", "limma","RUV","SVA","NNM")
  ## ntop = 4000; xlist.init = list()
  batch <- NULL
  pars <- get_model_parameters(X, samples, pheno = pheno, contrasts = contrasts)

  message("Running batch-correction methods...")
  xlist <- runBatchCorrectionMethods(
    X = X,
    batch = pars$batch,
    y = pars$pheno,
    controls = NULL,
    methods = methods,
    combatx = FALSE,
    ntop = ntop,
    sc = FALSE,
    remove.failed = TRUE
  )
  names(xlist)
  if (length(xlist.init) > 0) xlist <- c(xlist.init, xlist)


  dbg("[compare_batchcorrection_methods] names.xlist = ", names(xlist))
  dbg("[compare_batchcorrection_methods] xlist.nrow = ", sapply(xlist, nrow))
  dbg("[compare_batchcorrection_methods] xlist.ncol = ", sapply(xlist, ncol))

  ## PCA is faster than UMAP
  pos <- list()
  t2 <- function(x) t(scale(t(scale(t(x), scale = FALSE))))
  nb <- max(2, round(min(30, dim(X) / 5)))

  ##  incProgress( amount = 0.1, "Computing PCA clustering...")
  message("Computing PCA clustering...")
  pos[["pca"]] <- lapply(xlist, function(x) {
    irlba::irlba(t2(x), nu = 2, nv = 2)$u[, 1:2]
  })
  for (i in 1:length(pos[["pca"]])) {
    rownames(pos[["pca"]][[i]]) <- colnames(X)
  }

  ##  incProgress( amount = 0.1, "Computing t-SNE clustering...")
  message("Computing t-SNE clustering...")
  message("[Comuputing] nb = ", nb)
  pos[["tsne"]] <- lapply(xlist, function(x) {
    Rtsne::Rtsne(t2(x), perplexity = nb, check_duplicates = FALSE)$Y
  })
  for (i in 1:length(pos[["tsne"]])) {
    rownames(pos[["tsne"]][[i]]) <- colnames(X)
  }

  ##  incProgress( amount = 0.1, "Comparing results...")
  message("Comparing results...")
  res <- playbase::bc.evaluateResults(
    xlist,
    pheno = pars$pheno,
    lfc = 0.2,
    q = 0.05,
    pos = pos[["tsne"]],
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

  list(
    xlist = xlist,
    pos = pos,
    scores = res$scores,
    pheno = pars$pheno,
    pars = pars,
    best.method = best.method
  )
}

## ================================================================================
## Single batch-correction methods wrappers
## ================================================================================


#' @export
superBC2 <- function(X, samples, y, batch = NULL,
                     ## methods = c("technical","batch","statistical","pca","sva","nnm"),
                     methods = c("batch", "technical", "statistical", "sva", "nnm2"),
                     p.pca = 0.5, p.pheno = 0.05, k.pca = 10, nv = 2,
                     xrank = NULL, use.design = TRUE) {
  if (y[1] %in% colnames(samples)) {
    y <- samples[, y[1]]
  }
  if (!is.null(batch) && grepl("batch", colnames(samples), ignore.case = TRUE)) {
    batch.col <- grep("batch", colnames(samples), ignore.case = TRUE)[1]
    message("[superBC2] found batch column in sample info: ", colnames(samples)[batch.col])
    batch <- samples[, batch.col]
  }
  if (!is.null(batch) && batch[1] %in% colnames(samples)) {
    message("[superBC2] using batch column in sample info: ", batch[1])
    batch <- samples[, batch[1]]
  }

  cX <- X
  methods <- intersect(methods, c(
    "technical", "batch", "statistical", "pca", "sva",
    "ruv", "nnm", "nnm2"
  ))

  for (m in methods) {
    ## correct explicit batch effect
    if (!is.null(batch) && m == "batch") {
      message("[superBC2] correcting for: batch")
      if (use.design) {
        mod1 <- model.matrix(~y)
        cX <- limma::removeBatchEffect(cX, batch = batch, design = mod1)
      } else {
        cX <- limma::removeBatchEffect(cX, batch = batch)
      }
    }

    ## this removes typical batch effects
    if (m %in% c("statistical", "technical", "pca")) {
      bc <- detectBatchEffects(
        ## cX,
        X, ## on original matrix??
        samples, y,
        params = m,
        p.pca = p.pca, p.pheno = p.pheno, k.pca = k.pca,
        nv = nv, xrank = xrank
      )
      if (!is.null(bc$covariates)) {
        message("[superBC2] correcting for: ", m)
        B <- scale(bc$covariates)
        B[is.nan(B) | is.na(B)] <- 0
        B[is.infinite(B)] <- 0
        if (use.design) {
          mod1 <- model.matrix(~ bc$pheno)
          cX <- limma::removeBatchEffect(cX, covariates = B, design = mod1)
        } else {
          cX <- limma::removeBatchEffect(cX, covariates = B)
        }
      }
    }

    ## 5. additional unsupervised correction
    if (m == "sva") {
      message("[superBC2] correcting for: SVA")
      cX <- svaCorrect(cX, y)
    }
    if (m == "ruv") {
      message("[superBC2] correcting for: RUV")
      cX <- ruvCorrect(cX, y)
    }
    if (m == "nnm") {
      message("[superBC2] correcting for: NNM")
      cX <- gx.nnmcorrect(cX, y, use.design = use.design)$X
    }
    if (m == "nnm2") {
      message("[superBC2] correcting for: NNM2")
      cX <- gx.nnmcorrect2(cX, y, r = 0.35, use.design = use.design)$X
    }
  }

  cX <- cX - rowMeans(cX, na.rm = TRUE) + rowMeans(X, na.rm = TRUE)
  cX
}


#' @export
svaCorrect <- function(X, y) {
  ##
  ## This is a combination of methods from SVA and SmartSVA
  ## because of speed.

  if (any(is.na(X))) {
    stop("[svaCorrect] cannot handle missing values in X")
  }
  dimx <- ncol(X)

  ## sva doesn't like if the dimension is too small
  if (ncol(X) < 10) {
    X <- cbind(X, X, X)
    y <- rep(y, 3)
  }

  mod1x <- model.matrix(~ 1 + y)
  mod0x <- mod1x[, 1, drop = FALSE] ## just ones...

  ## fast method using SmartSVA
  #  pp <- paste0(model.par, collapse = "+")
  #  lm.expr <- paste0("lm(t(X) ~ ", pp, ", data=pheno)")
  X.r <- t(stats::resid(lm(t(X) ~ y)))
  n.sv <- isva::EstDimRMT(X.r, FALSE)$dim + 1
  ## top 1000 genes only (faster)
  X1 <- Matrix::head(X[order(-matrixStats::rowSds(X, na.rm = TRUE)), ], 1000)
  ## add a little bit of noise to avoid singular error
  a <- 0.01 * mean(apply(X1, 1, stats::sd, na.rm = TRUE), na.rm = TRUE)
  X1 <- X1 + a * matrix(stats::rnorm(length(X1)), nrow(X1), ncol(X1))
  sv <- try(sva::sva(X1, mod1x, mod0 = mod0x, n.sv = pmax(n.sv - 1, 1))$sv)
  class(sv)
  if (!any(class(sv) == "try-error")) {
    message("[svaCorrect] Performing SVA correction...")
    rownames(sv) <- colnames(X)
    colnames(sv) <- paste0("SV.", 1:ncol(sv))
    X <- limma::removeBatchEffect(X, covariates = sv, design = mod1x)
  } else {
    message("[svaCorrect] WARNING could not get covariates. no correction.")
  }
  X <- X[, 1:dimx]
  X
}

#' @export
ruvCorrect <- function(X, y, k = NULL, type = c("III", "g"), controls = 0.10) {
  if (any(is.na(X))) {
    stop("[ruvCorrect] cannot handle missing values")
  }

  ## F-test using limma just variables
  if (length(controls) == 1 && is.numeric(controls[1])) {
    ii <- which(!duplicated(rownames(X)))
    F <- gx.limmaF(X[ii, ], y, lfc = 0, fdr = 1, method = 1, sort.by = "none", compute.means = FALSE, verbose = 0)
    nc <- pmax(nrow(X) * as.numeric(controls), 1)
    sel <- head(order(-F$P.Value), nc)
    controls <- rownames(F)[sel]
    controls <- intersect(controls, rownames(X))
  }
  message(paste("[ruvCorrect] Number of control features:", length(controls)))

  if (is.null(k)) {
    ## this is from SVA
    ##    X.r <- t(stats::resid(lm(t(X) ~ y)))
    ##    k <- isva::EstDimRMT(X.r, FALSE)$dim + 1
    fit <- lm(t(X[controls, , drop = FALSE]) ~ y)
    X.r <- t(stats::resid(fit))
    k <- try(isva::EstDimRMT(X.r, FALSE)$dim + 1)
    if ("try-error" %in% class(k)) {
      k0 <- isva::EstDimRMT(X, FALSE)$dim
      k <- ncol(X) - k0
    }
  }
  message(paste("[ruvCorrect] Number of significant surrogate variables is:", k))

  M <- model.matrix(~ 0 + y)
  rownames(M) <- colnames(X)

  # Actually run correction
  # Set number of threads for CRAN
  ctl <- which(controls %in% rownames(X))
  type <- type[1]
  if (type == "III") {
    ruvX <- t(ruv::RUVIII(Y = t(X), M = M, ctl = ctl, k = k, eta = NULL))
    ruvX <- ruvX - mean(ruvX) + mean(X) ## ??
  } else if (type == "g") {
    ruvX <- RUVSeq::RUVg(x = X, cIdx = ctl, k = k, isLog = TRUE)$normalizedCounts
  } else {
    stop("[ruvCorrect] unknown RUV type", type)
  }
  ruvX
}


#' @export
pcaCorrect.OLD <- function(X, y, k = NULL, max.rho = 0.3) {
  ## --------------------------------------------------------------------
  ## PCA correction: remove remaining batch effect using PCA
  ## (iteratively, only SV larger than max correlated SV)
  ## --------------------------------------------------------------------
  mod1 <- model.matrix(~ 0 + y)
  ii <- 1:99
  niter <- 0
  if (is.null(k)) k <- ncol(X)
  nremoved <- 0
  cX <- X
  while (length(ii) > 0 && niter < k) {
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
      nremoved <- nremoved + 1
    }
    niter <- niter + 1
  }
  niter
  if (niter == k) {
    dbg("WARNING:: PCA correction did not converge after", nremoved, "iterations\n")
  } else {
    dbg("[pcaCorrect] removed", nremoved, "principal components\n")
  }

  ## bring back mean
  cX <- cX - rowMeans(cX, na.rm = TRUE) + rowMeans(X, na.rm = TRUE)
  cX
}

#' @export
pcaCorrect2 <- function(X, y, k = 10, p.notsig = 0.20) {
  ## --------------------------------------------------------------------
  ## PCA correction: remove remaining batch effect using PCA
  ## (iteratively, only SV larger than max correlated SV)
  ## --------------------------------------------------------------------
  mod1 <- model.matrix(~ 0 + y)
  k <- min(k, ncol(X) - 1)
  suppressWarnings(suppressMessages({
    if (k < min(dim(X)) / 3) {
      V <- irlba::irlba(X, nv = k, nu = 0)$v
    } else {
      V <- svd(X, nv = k, nu = 0)$v
    }
  }))
  rownames(V) <- colnames(X)
  colnames(V) <- paste0("PC", 1:ncol(V))
  res <- gx.limmaF(t(V), y, lfc = 0, fdr = 1, sort.by = "none", compute.means = FALSE, verbose = 0)
  res <- res[colnames(V), ]
  round(res$P.Value, 4)

  ## select pc comps that are not correlated with pheno and are
  ## smaller than strongest pheno PC component
  xrank <- which.min(res$P.Value)
  sel <- which(res$P.Value > p.notsig & 1:nrow(res) < xrank)
  sel
  cX <- X
  if (length(sel)) {
    V <- V[, sel, drop = FALSE]
    ## remove batch-suspected PC vectors
    mod1 <- model.matrix(~ 0 + y)
    cX <- limma::removeBatchEffect(X, covariates = V, design = mod1)
    nremoved <- ncol(V)
    dbg("[pcaCorrect2] removed", nremoved, "principal components\n")
  } else {
    dbg("[pcaCorrect2] no correction\n")
  }

  ## bring back mean
  cX <- cX - rowMeans(cX, na.rm = TRUE) + rowMeans(X, na.rm = TRUE)
  cX
}

#' @export
pcaCorrect3 <- function(X, y, k = 10, xrank = NULL, p.notsig = 0.20) {
  ## this removes typical batch effects
  cX <- X
  bc <- detectBatchEffects(cX, samples, y,
    params = "pca",
    p.pca = 0.5, p.pheno = p.notsig, k.pca = k, xrank = xrank
  )
  if (!is.null(bc$covariates)) {
    mod1 <- model.matrix(~y)
    B <- scale(bc$covariates)
    B[is.nan(B) | is.na(B)] <- 0
    B[is.infinite(B)] <- 0
    cX <- limma::removeBatchEffect(cX, covariates = B, design = mod1)
  }
  cX
}

#' @export
runHarmony <- function(X, batch) {
  library(Seurat)
  library(harmony)
  nx <- ncol(X)
  cn <- colnames(X)
  X1 <- X
  if (ncol(X) < 10) {
    X1 <- cbind(X, X, X)
    batch <- rep(batch, 3)
  }
  colnames(X1) <- paste0("col", 1:ncol(X1), "_", colnames(X1))
  M <- data.frame(batch = batch)
  rownames(M) <- colnames(X1)
  if (is.null(rownames(X1))) rownames(X1) <- paste0("row", 1:nrow(X1))
  obj <- CreateSeuratObject(2**X1, meta.data = M)
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj, verbose = FALSE)
  sX <- t(scale(t(X1)))
  hvg <- VariableFeatures(obj)
  sX <- sX[hvg, ]
  obj[["RNA"]]@layers[["scale.data"]] <- sX
  npcs <- min(30L, ncol(X1) / 2)
  npcs
  nn <- npcs
  obj <- RunPCA(obj, npcs = npcs, verbose = FALSE)
  obj <- RunUMAP(obj, reduction = "pca", dims = 1:npcs, n.neighbors = nn)
  pos0 <- obj@reductions[["umap"]]@cell.embeddings
  ## DimPlot(obj, reduction="umap", group.by='batch',pt.size=5)
  hm.obj <- harmony::RunHarmony(obj, "batch", verbose = TRUE, sigma = 0.1)
  hm.obj <- RunUMAP(hm.obj, reduction = "harmony", dims = 1:npcs, n.neighbors = nn)
  ## DimPlot(hm.obj, reduction="umap", group.by="batch", pt.size=4)
  pos1 <- hm.obj@reductions[["umap"]]@cell.embeddings
  ## reconstruct corrected matrix
  u2 <- Loadings(hm.obj, "harmony")
  u2[is.nan(u2) | is.na(u2)] <- 0 ## sometimes...
  u2 <- t(t(u2) / sqrt(1e-4 + colSums(u2**2)))
  v2 <- Embeddings(hm.obj, "harmony")
  X2 <- as.matrix(u2 %*% t(v2))

  X2 <- X2[, 1:nx]
  pos0 <- pos0[1:nx, ]
  pos1 <- pos1[1:nx, ]
  colnames(X2) <- rownames(pos0) <- rownames(pos0) <- cn
  list(corrected = X2, umap = pos1, umap.orig = pos0)
}

#' @export
ComBatX <- function(X, batch, y = NULL, controls = NULL, b = 50,
                    recenter = TRUE, add.star = TRUE, bc.dim = 3,
                    bc.method = "combat") {
  if (0) {
    y <- NULL
    controls <- NULL
    b <- 50
    recenter <- TRUE
    add.star <- TRUE
    bc.dim <- 3
    bc.method <- "combat"
  }

  batch <- paste0("b_", as.character(batch))
  ## Get anchors: phenotypes that are in more than one batches
  if (is.null(y)) {
    if (!is.null(controls)) controls <- NULL
    y <- rep("all.samples", ncol(X))
  }
  anchors <- names(which(table(y) > 1))
  if (!is.null(batch)) anchors <- names(which(rowSums(table(y, batch) > 0) > 1))
  anchors <- setdiff(anchors, NA)
  if (!is.null(controls)) {
    if (!any(controls %in% y)) {
      stop("[ComBatX] controls do not match phenotype")
    }
    anchors <- intersect(anchors, controls)
  }
  message("Found ", length(anchors), " anchor conditions")
  message("anchors: ", paste(anchors, sep = "", collapse = " "))

  if (!is.null(batch) && add.star) anchors <- c("*", anchors)

  ## For each anchor determine correction vector
  cat("Computing correction vectors\n")
  ref <- anchors[1]
  dx.list <- list()
  message("[ComBatX] using batch vector")
  message("[ComBatX] bc.dim = ", bc.dim)
  for (ref in anchors) {
    if (ref == "*") {
      ii <- 1:ncol(X)
    } else {
      ii <- which(y == ref)
    }
    X0 <- X[, ii, drop = FALSE]
    dX <- NULL
    has.batch <- length(unique(batch[ii])) > 1
    if (!is.null(batch) && has.batch) {
      bb <- batch[ii]
      if (bc.method == "combat") {
        if (ncol(X0) == 2) {
          mx <- rowMeans(X0)
          X1 <- cbind(mx, mx)
        } else {
          X1 <- suppressWarnings(suppressMessages(
            sva::ComBat(X0, batch = bb)
          ))
        }
      } else {
        X1 <- limma::removeBatchEffect(X0, batch = bb)
      }
      dX <- X1 - X0
    } else if (NCOL(X0) > 1) {
      dX <- X0 - rowMeans(X0, na.rm = TRUE)
      bb <- rep(1, ncol(dX))
    }

    if (!is.null(dX)) {
      ## get correction vectors
      if (bc.dim < 0) {
        dxx <- dX
      }
      if (bc.dim == 0) {
        ## determine delta vectors in each batch
        dxx <- tapply(1:ncol(dX), bb, function(ii) {
          if (length(ii) == 1) ii <- c(ii, ii)
          rowMeans(dX[, ii])
        })
        dxx <- do.call(cbind, dxx)
      }
      if (bc.dim > 0) {
        ## determine delta vectors in each batch
        if (0) {
          dxx <- tapply(1:ncol(dX), bb, function(ii) {
            if (length(ii) == 1) ii <- c(ii, ii)
            nu <- min(bc.dim, length(ii) - 1)
            res <- irlba::irlba(dX[, ii], nu = nu)
            res$u[, 1:nu]
          })
        } else {
          dxx <- list()
          for (b1 in unique(bb)) {
            ii <- which(bb == b1)
            if (length(ii) == 1) {
              dxx[[b1]] <- dX[, ii]
            } else {
              nu <- min(bc.dim, length(ii) - 1)
              res <- irlba::irlba(dX[, ii], nu = nu)
              dxx[[b1]] <- res$u[, 1:nu]
            }
          }
        }
        dxx <- do.call(cbind, dxx)
      }
    }
    dx.list[[ref]] <- dxx
  }

  dX <- do.call(cbind, dx.list)
  dim(dX)
  cat("Found", ncol(dX), "correction vectors\n")

  ## combat sometimes give
  sum(is.na(dX))
  dX[is.na(dX)] <- 0

  ## Correct original data with all correction vectors
  # dX <- irlba::irlba(dX, nv=round(0.9*ncol(dX)))$u
  nv <- min(b, ncol(dX) - 1)
  cat("Factorizing differences using", nv, "dimensions\n")
  sv <- try(irlba::irlba(dX, nv = nv))
  if ("try-error" %in% class(sv)) {
    sv <- svd(dX, nv = nv)
  }
  dX <- sv$u[, 1:nv, drop = FALSE]

  cat("Correcting input matrix\n")
  cX <- t(limma::removeBatchEffect(t(X), covariates = dX))

  ## recenter on original mean?
  if (recenter) {
    cX <- cX - rowMeans(cX) + rowMeans(X)
  }

  cX
}


#' @export
MNNcorrect <- function(X, batch, controls = NULL) {
  if (is.null(controls)) {
    hvg <- head(rownames(X)[order(-apply(X, 1, sd, na.rm = TRUE))], 1000)
    res <- batchelor::fastMNN(X, batch = batch, subset.row = hvg, correct.all = TRUE)
  } else {
    xx <- tapply(1:ncol(X), batch, function(i) X[, i, drop = FALSE])
    rr <- tapply(1:ncol(X), batch, function(i) which(controls[i]))
    res <- batchelor::fastMNN(xx, restrict = rr)
  }
  cX <- as.matrix(SummarizedExperiment::assay(res))
  cX
}

#' @export
normalizeToControls <- function(X, batch, y, controls) {
  ii <- which(y %in% controls)
  batch.ctl <- tapply(ii, batch[ii], function(k) rowMeans(X[, k, drop = FALSE]))
  batch.ctl <- do.call(cbind, batch.ctl)
  nX <- (X - batch.ctl[, batch]) + rowMeans(X)
  ## nX <- limma::normalizeQuantiles(X)
  nX
}


##' @name bbknn State Matrix
##' @rdname bbknn
##'
##' @title Run bbknn clustering algorithm
##'
##' @description Implements the bbknn clustering algorithm in R using
##'   reticulate to run the Python version. Requires the python
##'   "bbknn" and "igraph" modules to be installed. Returns a vector
##'   of partition indices. From: https://rdrr.io/github/TomKellyGenetics/bbknn/src/R/bbknn.R
##'
##' @param data_matrix A matrix (genes x samples or cells) for expression data
##' @param batch An integer vector of batches to correct for (converts factors or numeric vectors)
##' @param pca whether to compute pca (defaults to TRUE) or apply correction to the raw matrix (FALSE)
##' @param compute_pca whether to compute PCA in Python (defaults to TRUE, requires scanpy library) or with R functions (FALSE)
##' @param nPcs number of principal components to compute (defaults to 50 if more than 50 genes)
##'
##' @return returns a list with the following components
##'   \item{corrected matrix}{matrix of data corrected by the BBKNN
##'   (batch based K nearest neighbours)}\item{pca}{principal
##'   components(matrix with row for every sample and column for each
##'   component)}\item{tsne}{t-distributed stochastic neighbour
##'   embedding (matrix with row for every sample)}\item{umap}{uniform
##'   manifold approximation and projection (matrix with row for every
##'   sample)}
##'
##' @keywords graph network igraph mvtnorm simulation
##' @import reticulate
##' @importFrom stats prcomp
##' @export
bbknn <- function(data_matrix, batch, pca = TRUE, compute_pca = "python", nPcs = NULL) {
  ## reticulate::py_install("anndata")
  ## reticulate::py_install("bbknn")

  # import python modules with reticulate
  if (!is.matrix(data_matrix)) {
    warning("matrix expected for data_matrix")
    data_matrix <- as.matrix(data_matrix)
  }
  if (is.null(nPcs)) {
    nPcs <- min(50, nrow(data_matrix), ncol(data_matrix))
  }
  if (nPcs > nrow(data_matrix)) {
    warning("number of genes less than nPcs")
    print(paste("using", nrow(data_matrix), "components"))
    ## nPcs <- nrow(data_matrix)
    nPcs <- min(nPcs, dim(data_matrix) - 1)
  }
  # reticulate::use_python("/usr/local/bin/python3")
  ##  reticulate::py_install("anndata")
  ##  reticulate::py_install("bbknn")
  ##  reticulate::py_install("scanpy")

  anndata <- reticulate::import("anndata", convert = FALSE)
  bbknn <- reticulate::import("bbknn", convert = FALSE)
  ##  sc <- reticulate::import("scanpy.api",convert=FALSE)
  sc <- reticulate::import("scanpy", convert = FALSE)

  # set up annotation data for batches
  if (is.character(batch)) batch <- as.factor(batch)
  if (is.factor(batch)) batch <- as.numeric(batch)
  if (is.numeric(batch)) batch <- as.integer(batch)

  # perform PCA
  if (pca) {
    # sc$tl$pca(adata)
    if (compute_pca == "python") {
      # use PCA computed in Python
      pca <- sc$pp$pca(t(data_matrix))
    } else if (compute_pca != "python") {
      # use PCA computed in R
      print("test")
      pca <- reticulate::r_to_py(t(prcomp(data_matrix)$x[1:nPcs, ]))
    }
    adata <- anndata$AnnData(X = pca, obs = batch)
    sc$tl$pca(adata)
    adata$obsm$X_pca <- pca
  } else {
    # use full matrix
    adata <- anndata$AnnData(X = t(data_matrix), obs = batch)
    sc$tl$pca(adata)
  }
  # perform BBKNN to derive corrected components
  bbknn$bbknn(adata, batch_key = 0)
  corrected_matrix <- t(reticulate::py_to_r(adata$X))
  sc$tl$pca(adata)
  pca <- reticulate::py_to_r(adata$obsm["X_pca"])
  sc$tl$tsne(adata)
  tsne <- reticulate::py_to_r(adata$obsm["X_tsne"])
  sc$tl$umap(adata)
  umap <- reticulate::py_to_r(adata$obsm["X_umap"])
  output <- list(corrected_matrix = corrected_matrix, pca = pca, tsne = tsne, umap = umap)
  return(output)
}


#' Nearest neighbor matching batch correction
#'
#' Correct for batch effects in a gene expression matrix using nearest neighbor matching.
#'
#' @param X Numeric matrix of gene expression values (genes in rows, samples in columns).
#' @param y Factor vector indicating batch for each sample.
#' @param dist.method Distance metric to use for matching ('cor' or 'euclidean').
#' @param center.x Logical for whether to center gene expression by row means.
#' @param center.m Logical for whether to center expression by batch means.
#' @param sdtop Number of top variable genes to use for correlation.
#'
#' @return List containing:
#' \itemize{
#'   \item X - Batch corrected gene expression matrix
#'   \item pairings - Matrix of sample pairings used for correction
#' }
#'
#' @details This function performs batch correction using the following steps:
#' \enumerate{
#'   \item Compute distance matrix between all samples
#'   \item Find nearest neighbor matches between batches
#'   \item Construct full paired dataset with matches
#'   \item Apply limma batch correction to paired data
#'   \item Average paired samples back to original samples
#' }
#'
#'
#' @seealso
#' \code{\link[limma]{removeBatchEffect}} for the batch correction method used.
#'
#' @examples
#' # TODO
#'
#' @export
nnmCorrect <- function(X, y, dist.method = "cor", center.x = TRUE, center.m = TRUE,
                       knn = 1, sdtop = 2000, return.B = FALSE,
                       use.design = TRUE, use.cov = FALSE) {
  ## Nearest-neighbour matching for batch correction. This
  ## implementation creates a fully paired dataset with nearest
  ## matching neighbours when pairs are missing.

  ## use.design=TRUE;dist.method="cor";center.x=TRUE;center.m=TRUE;sdtop=1000;knn=2

  ## compute distance matrix for NNM-pairing
  y1 <- paste0("y=", y)
  dX <- X

  ## reduce for speed
  sdx <- apply(dX, 1, stats::sd)
  ii <- Matrix::head(order(-sdx), sdtop)
  dX <- dX[ii, ]

  if (center.x) {
    dX <- dX - rowMeans(dX, na.rm = TRUE)
  }
  if (center.m) {
    ## center per condition group (takes out batch differences)
    mX <- tapply(1:ncol(dX), y1, function(i) rowMeans(dX[, i, drop = FALSE]))
    mX <- do.call(cbind, mX)
    dX <- dX - mX[, y1]
  }

  if (dist.method == "cor") {
    message("[nnmCorrect] computing correlation matrix D...")
    ## D <- 1 - crossprod(scale(dX)) / (nrow(dX) - 1) ## faster
    D <- 1 - cor(dX)
  } else {
    message("[nnmCorrect] computing distance matrix D...\n")
    D <- as.matrix(stats::dist(t(dX)))
  }
  ## remove(dX)
  D[is.na(D)] <- 0 ## might have NA

  ## find neighbours
  if (knn > 1) {
    message(paste0("[nnmCorrect] finding ",knn,"-nearest neighbours..."))
    bb <- apply(D, 1, function(r) tapply(r, y1, function(s) head(names(sort(s)), knn)))
    B <- do.call(rbind, lapply(bb, function(x) unlist(x)))
    colnames(B) <- unlist(mapply(rep, names(bb[[1]]), sapply(bb[[1]],length)),use.names=FALSE)
  } else {
    message("[nnmCorrect] finding nearest neighbours...")       
    B <- t(apply(D, 1, function(r) tapply(r, y1, function(s) names(which.min(s)))))
  }
  rownames(B) <- colnames(X)

  ## ensure sample is always present in own group
  idx <- cbind(1:nrow(B), match(y1, colnames(B)))
  B[idx] <- rownames(B)
##  B <- cbind(rownames(B), B)

  ## imputing full paired data set
  kk <- match(as.vector(B), rownames(B))
  full.y <- y1[kk]
  full.pairs <- rep(rownames(B), ncol(B))
  full.X <- X[, kk]
  dim(full.X)

  ## remove pairing effect
  message("[nnmCorrect] correcting for pairing effects...")
  if (use.cov == FALSE) {
    design <- stats::model.matrix(~full.y)
    if (!use.design) design <- matrix(1, ncol(full.X), 1)
    full.X <- limma::removeBatchEffect(
      full.X,
      batch = full.pairs,
      design = design
    )
  } else {
    V <- model.matrix(~ 0 + full.pairs)
    design <- stats::model.matrix(~full.y)
    if (!use.design) design <- matrix(1, ncol(full.X), 1)
    full.X <- limma::removeBatchEffect(
      full.X,
      covariates = scale(V),
      design = design
    )
  }

  ## now contract to original samples
  message("[nnmCorrect] matching result...")
  full.idx <- rownames(B)[kk]
  cX <- do.call(cbind, tapply(
    1:ncol(full.X), full.idx,
    function(i) rowMeans(full.X[, i, drop = FALSE])
  ))
  cX <- cX[, colnames(X)]

  ## retain original row means
  cX <- cX - rowMeans(cX, na.rm = TRUE) + rowMeans(X, na.rm = TRUE)
  res <- cX
  if (return.B) {
    res <- list(X = cX, pairings = B)
  }
  return(res)
}


#' @export
nnmCorrect2 <- function(X, y, r = 0.35, center.x = TRUE, center.m = TRUE,
                        scale.x = FALSE, center.y = TRUE, mode = "sym",
                        knn = 1, sdtop = 2000, return.B = FALSE, use.design = TRUE) {

##  center.x=TRUE;center.m=TRUE;scale.x=FALSE;sdtop=1000;r=0.35;knn=5
    
  ## compute distance matrix for NNM-pairing
  y1 <- paste0("y=", y)
  dX <- X

  ## reduce for speed
  sdx <- apply(dX, 1, stats::sd)
  ii <- Matrix::head(order(-sdx), sdtop)
  dX <- dX[ii, ]
  if (center.x) {
    dX <- dX - rowMeans(dX, na.rm = TRUE)
  }
  if (scale.x) {
    row.sdx <- matrixStats::rowSds(dX, na.rm = TRUE)
    dX <- dX / (row.sdx + 1e-4 * mean(row.sdx, na.rm = TRUE))
  }
  if (center.m) {
    ## center per condition group (takes out batch differences)
    mX <- tapply(1:ncol(dX), y1, function(i) rowMeans(dX[, i, drop = FALSE]))
    mX <- do.call(cbind, mX)
    dX <- dX - mX[, y1]
  }
  if (center.y) {
    ## finally sample scaling
    dX <- scale(dX)
    dX[is.na(dX)] <- 0 ## zero SD can cause NA
  }

  ## find neighbours using fast KNN search
  message("[nnmCorrect2] finding nearest neighbours...")
  a <- y1[1]
  bb <- list()
  nn <- list()
  for (a in sort(unique(y1))) {
    x1 <- dX[, which(y1 == a), drop = FALSE]
    knn1 <- min(knn, ncol(x1))
    res <- FNN::get.knnx(t(x1), query = t(dX), k = knn1)
    bb[[a]] <- apply(res$nn.index, 2, function(i) colnames(x1)[i])
    nn[[a]] <- knn1
  }
  B <- do.call(cbind, bb)
  colnames(B) <- as.vector(unlist(mapply(rep, names(bb), nn)))
  rownames(B) <- colnames(dX)

  if(1) {
    ## delete neighbours in own group (???)
    jj <- lapply(y1, function(a) which(colnames(B)==a)) 
    ii <- mapply(rep, 1:length(y1), sapply(jj,length))
    idx <- cbind( as.vector(unlist(ii)), as.vector(unlist(jj)))
    B[idx] <- NA
  }
    
  ## ensure sample is always present in own group
  idx <- cbind(1:nrow(B), match(y1, colnames(B)))
  B[idx] <- rownames(B)

  ## create pairing design matrix manually
  idx <- apply(B, 1, function(x) match(x, rownames(B)))
  jj <- as.vector(idx)
  ii <- as.vector(mapply(rep, 1:ncol(idx), nrow(idx)))
  ii <- ii[!is.na(jj)]
  jj <- jj[!is.na(jj)]    
  P <- Matrix::sparseMatrix(
    i = jj, j = ii, x = rep(1, length(ii)),
    dims = c(nrow(B), nrow(B))
  )
  P <- as.matrix(P)    
  P <- 1 * (P > 0) ## dupcliated got summed
    
  ## correct for pairing effect
  message("[nnmCorrect2] correcting for pairing effects...")
  P1 <- P
  if (mode == "sym") P1 <- P + Matrix::t(P) ## make symmetric
  if (mode == "tr") P1 <- Matrix::t(P) ## transposed

  ## take out duplicate columns  
  P1 <- P1[,!duplicated.matrix( t(as.matrix(P1))), drop=FALSE]  
  dim(P1)

  if (r < 1) {
    k <- round(min(r * dim(P), dim(P) - 1)) ## critical
    k <- max(k, 1)
    if (r > 0.2) {
      sv <- svd(P1, nu = k, nv = k)
    } else {
      sv <- irlba::irlba(P1, nu = k, nv = k)
    }
    P1 <- sv$u
  }

  design <- stats::model.matrix(~y1)
  if (!use.design) design <- matrix(1, ncol(X), 1)
  cX <- limma::removeBatchEffect(X, covariates = scale(P1), design = design)

  ## retain original row means
  cX <- cX - rowMeans(cX, na.rm = TRUE) + rowMeans(X, na.rm = TRUE)
  res <- cX
  if (return.B) {
    res <- list(X = cX, pairings = B)
  }
  return(res)
}

#' Nearest neighbor matching batch correction
#'
#' Correct for batch effects in a gene expression matrix using nearest neighbor matching.
#'
#' @param x Numeric matrix of gene expression values (genes in rows, samples in columns).
#' @param y Factor vector indicating batch for each sample.
#' @param k Number of nearest neighbors to use (default 3).
#'
#' @return Batch corrected gene expression matrix.
#'
#' @details This function performs batch correction by matching each sample to its
#' k-nearest neighbors from other batches based on expression profile correlation.
#' The batch effect for a sample is estimated as the mean difference between it and its
#' matched neighbors. This difference is subtracted from the sample's expression profile.
#'
#' @seealso
#' \code{\link[limma]{removeBatchEffect}} for an alternative batch correction method
#'
#' @examples
#' \dontrun{
#' x <- matrix(rnorm(100 * 30), 100, 30) # random expression matrix
#' y <- gl(3, 10) # 3 batches of 10 samples each
#' xcorr <- nnmCorrect.SAVE(x, y)
#' }
#' @export
nnmCorrect.SIMPLE <- function(x, y, k = 3) {
  ## -----------------------------------------------------
  ## nearest-neighbour matching for batch correction
  ## -----------------------------------------------------
  xcor <- stats::cor(x)
  diag(xcor) <- 0
  nx <- x
  j <- 1
  for (j in 1:ncol(x)) {
    nj <- which(y != y[j])
    nn <- intersect(order(-xcor[j, ]), nj)
    nn <- Matrix::head(nn, k)
    nx[, j] <- x[, j] - rowMeans(x[, nn, drop = FALSE])
  }
  nx <- nx + rowMeans(x)
  return(nx)
}


## compatibility

#' @export
gx.nnmcorrect <- function(...) nnmCorrect(..., return.B = TRUE)

#' @export
gx.nnmcorrect2 <- function(...) nnmCorrect2(..., return.B = TRUE)


#' @export
bcKNN <- function(X, batch, y = NULL) {
  if (is.null(batch)) stop("batch must be provided")
  bb.list <- list()
  batch <- as.character(batch)
  b <- batch[1]
  for (b in sort(unique(batch))) {
    X1 <- X[, which(batch == b), drop = FALSE]
    res <- FNN::get.knnx(t(X1), query = t(X), k = 1)
    xx <- X1[, res$nn.index[, 1]]
    bb.list[[b]] <- xx
  }
  if (!is.null(y)) {
    bb.list <- lapply(bb.list, function(x) t(rowsum(t(x), y)))
  } else {
    bb.list <- lapply(bb.list, function(x) rowMeans(x))
  }
  bb.mean <- Reduce("+", bb.list) / length(bb.list)
  bb.list <- lapply(bb.list, function(x) (x - bb.mean))
  B <- do.call(cbind, bb.list)
  dim(B)
  bX <- t(limma::removeBatchEffect(t(X), covariates = scale(B))) ## model??
  bX
}

## =====================================================================================
## =========================== END OF FILE =============================================
## =====================================================================================
