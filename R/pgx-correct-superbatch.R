# This file is part of the Omics Playground project.
# Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
# Superbatch correction implementations retained by playbase.
# This file owns supervised, technical-effect, PCA, and nearest-neighbour methods.
# Keep these implementations behavior-compatible with their existing callers.

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
    if (nnm.correct) bc <- c(bc, "npm")
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
  Y <- pgx.computeTechnicalEffects(X)
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
        ## out <- mnnCorrect(cX, batch = b, cos.norm.out = FALSE)
        ## cX <- out@assays@data[["corrected"]]
        cX <- MNNcorrect(cX, batch = b)
      }
    }

    if (bc == "npm") {
      ## --------------------------------------------------------------------
      ## Nearest-neighbour matching (NNM)
      ## --------------------------------------------------------------------
      dbg("[pgx.superBatchCorrect] Correcting with nearest-neighbour matching (NPM)")
      dbg("[pgx.superBatchCorrect] model.par = ", model.par)
      y1 <- pheno[, model.par, drop = FALSE]
      y1 <- apply(y1, 1, paste, collapse = ":")
      cX <- playbase.preprocess::pp.batchCorrect(
        cX,
        target = y1,
        method = "NPM",
        center_x = TRUE,
        center_m = TRUE
      )
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

        ## cX1 <- Matrix::head(cX[order(-apply(cX, 1, stats::sd, na.rm = TRUE)), ], 1000) ## top 1000 genes only (faster)
        cX1 <- Matrix::head(cX[order(-matrixStats::rowSds(cX, na.rm = TRUE)), ], 1000)
        sv <- try(sva::sva(cX1, mod1x, mod0 = mod0x, n.sv = n.sv)$sv)

        if (any(class(sv) == "try-error")) {
          ## try again with little bit of noise...
          ## a <- 0.01 * mean(apply(cX, 1, stats::sd, na.rm = TRUE), na.rm = TRUE)
          a <- 0.01 * mean(matrixStats::rowSds(cX, na.rm = TRUE), na.rm = TRUE)
          cX1 <- cX + a * matrix(stats::rnorm(length(cX)), nrow(cX), ncol(cX))
          ## cX1 <- Matrix::head(cX1[order(-apply(cX1, 1, stats::sd, na.rm = TRUE)), ], 1000) ## top 1000 genes only (faster)
          cX1 <- Matrix::head(cX1[order(-matrixStats::rowSds(cX1, na.rm = TRUE)), ], 1000)
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
          set.seed(1234)
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
          ## xx <- Matrix::head(cX[order(-apply(cX, 1, stats::sd, na.rm = TRUE)), ], hc.top)
          xx <- Matrix::head(cX[order(-matrixStats::rowSds(cX, na.rm = TRUE)), ], hc.top)
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

  return(res)
}

#' @export
removeTechnicalEffects <- function(X, samples, y, p.pheno = 0.05, p.pca = 0.5,
                                   params = c("lib", "mito", "ribo", "cellcycle", "gender"),
                                   force = FALSE, nv = 1, k.pca = 10, xrank = NULL) {
  ##  p.pheno = 0.05;p.pca = 0.5;force = FALSE; nv = 2;k.pca = 10;xrank = NULL
  ##  params = c("lib","mito","ribo","cellcycle","gender")

  X1 <- X
  X1 <- .pgx_impute_svd2(X1)
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

## =============================================================================
## Single batch-correction methods wrappers
## =============================================================================

#' @export
pcaCorrect <- function(X, y, k = 10, p.notsig = 0.20) {
  ## --------------------------------------------------------------------
  ## PCA correction: remove remaining batch effect using PCA
  ## (iteratively, only SV larger than max correlated SV)
  ## --------------------------------------------------------------------
  mod1 <- model.matrix(~ 0 + y)
  k <- min(k, ncol(X) - 1)
  suppressWarnings(suppressMessages({
    if (k < min(dim(X)) / 3) {
      set.seed(1234)
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
    dbg("[pcaCorrect] removed", nremoved, "principal components\n")
  } else {
    dbg("[pcaCorrect] no correction\n")
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
normalizeToControls <- function(X, batch, y, controls) {
  if (!is.null(ncol(batch))) batch <- apply(batch, 1, paste, collapse = "_")
  ii <- which(y %in% controls)
  batch.ctl <- tapply(ii, batch[ii], function(k) rowMeans(X[, k, drop = FALSE], na.rm = TRUE))
  batch.ctl <- do.call(cbind, batch.ctl)
  nX <- (X - batch.ctl[, batch]) + rowMeans(X, na.rm = TRUE)
  nX
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
nnmCorrect2 <- function(X, y, r = 0.35, center.x = TRUE, center.m = TRUE,
                        scale.x = FALSE, center.y = TRUE, mode = "sym",
                        knn = 1, sdtop = 2000, return.B = FALSE,
                        use.design = TRUE, delete.kin = TRUE) {
  ##  center.x=TRUE;center.m=TRUE;scale.x=FALSE;sdtop=1000;r=0.35;knn=5

  ## compute distance matrix for NNM-pairing
  y1 <- paste0("y=", y)
  dX <- X

  ## reduce for speed
  ## sdx <- apply(dX, 1, stats::sd, na.rm = TRUE)
  sdx <- matrixStats::rowSds(dX, na.rm = TRUE)
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
    mX <- tapply(1:ncol(dX), y1, function(i) rowMeans(dX[, i, drop = FALSE], na.rm = TRUE))
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

  if (delete.kin) {
    ## delete neighbours in own group (???)
    jj <- lapply(y1, function(a) which(colnames(B) == a))
    ii <- mapply(rep, 1:length(y1), sapply(jj, length))
    idx <- cbind(as.vector(unlist(ii)), as.vector(unlist(jj)))
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
  P1 <- P1[, !duplicated.matrix(t(as.matrix(P1))), drop = FALSE]
  dim(P1)

  if (r < 1) {
    k <- round(min(r * dim(P), dim(P) - 1)) ## critical
    k <- max(k, 1)
    if (r > 0.2) {
      sv <- svd(P1, nu = k, nv = k)
    } else {
      set.seed(1234)
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
    nx[, j] <- x[, j] - rowMeans(x[, nn, drop = FALSE], na.rm = TRUE)
  }
  nx <- nx + rowMeans(x, na.rm = TRUE)
  return(nx)
}
