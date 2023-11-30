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

  message("[pgx.superBatchCorrect] 1 : dim.pheno = ", paste(dim(pheno), collapse = "x"), "\n")
  message("[pgx.superBatchCorrect] 1 : model.par = ", paste(model.par, collapse = "x"), "\n")
  message("[pgx.superBatchCorrect] 1 : colnames.pheno = ", paste(colnames(pheno), collapse = "x"), "\n")

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
      y1 <- (y1 > stats::median(y1, na.rm = TRUE))
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

    r1 <- stats::cor(t(x[, ii]), design, use = "pairwise")
    rowMeans(abs(r1))
  }


  X <- X - rowMeans(X, na.rm = TRUE) ## center features
  X[is.na(X)] <- mean(X, na.rm = TRUE) ## no missing allowed
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
  if (stat == "F") R <- t(t(R) / colMeans(R, na.rm = TRUE))

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
    tx <- pmax(X - q0, 0, na.rm = TRUE) ## log expression
    cx <- pmax(2**tx - 1, 0, na.rm = TRUE) ## counts
  } else {
    cx <- X
    tx <- log2(cx + 1)
  }
  nfeature <- Matrix::colSums(cx > 0, na.rm = TRUE) + 1
  libsize <- Matrix::colSums(cx, na.rm = TRUE)

  mt.genes <- grep("^MT-", rownames(X), ignore.case = TRUE, value = TRUE)
  rb.genes <- grep("^RP[SL]", rownames(X), ignore.case = TRUE, value = TRUE)
  mito <- ribo <- NA
  pct.mito <- pct.ribo <- NA

  if (length(mt.genes) >= 10) {
    mito <- Matrix::colMeans(tx[mt.genes, , drop = FALSE], na.rm = TRUE)
    pct.mito <- Matrix::colSums(cx[mt.genes, , drop = FALSE], na.rm = TRUE) / libsize
  }
  if (length(rb.genes) >= 10) {
    ii <- rb.genes[order(-apply(tx[rb.genes, , drop = FALSE], 1, stats::sd, na.rm = TRUE))]
    sel20 <- Matrix::head(ii, 20)
    ribo <- Matrix::colMeans(tx[rb.genes, , drop = FALSE], na.rm = TRUE)
    ribo20 <- Matrix::colMeans(tx[sel20, , drop = FALSE], na.rm = TRUE)
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


#' @export
svaCorrect <- function(X, y) {
  ##
  ## This is a combination of methods from SVA and SmartSVA
  ## because of speed.

  if (any(is.na(X))) {
    stop("[svaCorrect] cannot handle missing values in X")
  }

  mod1x <- model.matrix(~ 1 + y)
  mod0x <- mod1x[, 1, drop = FALSE] ## just ones...

  ## fast method using SmartSVA
  #  pp <- paste0(model.par, collapse = "+")
  #  lm.expr <- paste0("lm(t(X) ~ ", pp, ", data=pheno)")
  X.r <- t(stats::resid(lm(t(X) ~ y)))
  n.sv <- isva::EstDimRMT(X.r, FALSE)$dim + 1
  X1 <- Matrix::head(X[order(-matrixStats::rowSds(X, na.rm = TRUE)), ], 1000) ## top 1000 genes only (faster)
  sv <- try(sva::sva(X1, mod1x, mod0 = mod0x, n.sv = n.sv)$sv)

  if (any(class(sv) == "try-error")) {
    ## try again with little bit of noise...
    a <- 0.01 * mean(apply(X, 1, stats::sd, na.rm = TRUE), na.rm = TRUE)
    X1 <- X + a * matrix(stats::rnorm(length(X)), nrow(X), ncol(X))
    X1 <- Matrix::head(X1[order(-matrixStats::rowSds(X, na.rm = TRUE)), ], 1000) ## top 1000 genes only (faster)
    sv <- try(sva::sva(X1, mod1x, mod0 = mod0x, n.sv = pmax(n.sv - 1, 1))$sv)
  }
  if (!any(class(sv) == "try-error")) {
    message("[svaCorrect] Performing SVA correction...")
    rownames(sv) <- colnames(X)
    colnames(sv) <- paste0("SV.", 1:ncol(sv))
    X <- limma::removeBatchEffect(X, covariates = sv, design = mod1x)
  } else {
    message("[svaCorrect] WARNING could not get covariates. no correction.")
  }
  X
}

#' @export
ruvCorrect <- function(X, y, k = NULL, type = c("III", "g"), controls = NULL) {
  if (any(is.na(X))) {
    stop("[ruvCorrect] cannot handle missing values")
  }

  ## F-test using limma just variables
  if (is.null(controls)) {
    F <- gx.limmaF(X, y, lfc = 0, fdr = 1, method = 1, compute.means = FALSE)
    F <- F[rownames(X), ]
    sel <- head(order(-F$P.Value), nrow(X) / 4)
    controls <- rownames(F)[sel]
  }

  if (is.null(k)) {
    ## this is from SVA
    ##    X.r <- t(stats::resid(lm(t(X) ~ y)))
    ##    k <- isva::EstDimRMT(X.r, FALSE)$dim + 1
    X.r <- t(stats::resid(lm(t(X[controls, ]) ~ y)))
    k <- isva::EstDimRMT(X.r, FALSE)$dim + 1
  }
  message(paste("[ruvCorrect] Number of significant surrogate variables is:", k))

  M <- model.matrix(~ y - 1)
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


runHarmony <- function(X, batch) {
  library(Seurat)
  library(harmony)
  M <- data.frame(batch = batch)
  rownames(M) <- colnames(X)
  table(batch)
  if (is.null(rownames(X))) rownames(X) <- paste0("row", 1:nrow(X))
  obj <- CreateSeuratObject(2**X, meta.data = M)
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj, verbose = FALSE)
  sX <- t(scale(t(X)))
  hvg <- VariableFeatures(obj)
  sX <- sX[hvg, ]
  obj[["RNA"]]@layers[["scale.data"]] <- sX
  obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
  obj <- RunUMAP(obj, reduction = "pca", dims = 1:30)
  pos0 <- obj@reductions[["umap"]]@cell.embeddings
  ## DimPlot(obj, reduction="umap", group.by='batch',pt.size=5)
  hm.obj <- harmony::RunHarmony(obj, "batch", verbose = TRUE)
  hm.obj <- RunUMAP(hm.obj, reduction = "harmony", dims = 1:30)
  ## DimPlot(hm.obj, reduction="umap", group.by="batch", pt.size=4)
  pos1 <- hm.obj@reductions[["umap"]]@cell.embeddings
  ## reconstruct corrected matrix
  u2 <- Loadings(hm.obj, "harmony")
  u2[is.nan(u2) | is.na(u2)] <- 0 ## sometimes...
  u2 <- t(t(u2) / sqrt(colSums(u2**2)))
  v2 <- Embeddings(hm.obj, "harmony")
  X2 <- as.matrix(u2 %*% t(v2))
  list(corrected = X2, umap = pos1, umap.orig = pos0)
}


ComBatX <- function(X, batch, pheno = NULL, controls = NULL, b = 50,
                    recenter = TRUE, add.star = TRUE, bc.dim = 0,
                    bc.method = "combat") {
  ## Get anchors: phenotypes that are in more than one batches
  if (is.null(pheno)) {
    if (!is.null(controls)) controls <- NULL
    pheno <- rep("all.samples", ncol(X))
  }
  anchors <- names(which(table(pheno) > 1))
  if (!is.null(batch)) anchors <- names(which(rowSums(table(pheno, batch) > 0) > 1))
  anchors <- setdiff(anchors, NA)
  if (!is.null(controls)) {
    if (!any(controls %in% pheno)) {
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
      ii <- which(pheno == ref)
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


MNNcorrect <- function(X, batch, controls = NULL) {
  if (is.null(controls)) {
    hvg <- head(rownames(X)[order(-apply(X, 1, sd, na.rm = TRUE))], 1000)
    res <- batchelor::fastMNN(X, batch = batch, subset.row = hvg, correct.all = TRUE)
  } else {
    xx <- tapply(1:ncol(X), batch, function(i) X[, i, drop = FALSE])
    rr <- tapply(1:ncol(X), batch, function(i) which(pheno[i] %in% controls))
    res <- batchelor::fastMNN(xx, restrict = rr)
  }
  cX <- as.matrix(SummarizedExperiment::assay(res))
  cX
}

normalizeToControls <- function(X, pheno, controls) {
  ii <- which(pheno %in% controls)
  ctl.means <- tapply(colnames(X)[ii], batch[ii], function(k) rowMeans(X[, k, drop = FALSE]))
  ctl.means <- do.call(cbind, ctl.means)
  nX <- (X - ctl.means[, batch]) + rowMeans(X)
  ## nX <- limma::normalizeQuantiles(X)
  nX
}


runBatchCorrectionMethods <- function(X, batch, pheno, controls, combinations = TRUE) {
  ## these are methods commonly used in scRNA-seq.
  ##

  table(pheno, batch)
  mod <- model.matrix(~pheno)
  nlevel <- length(unique(pheno[!is.na(pheno)]))

  xlist <- list()
  xlist[["uncorrected"]] <- X

  ## normalize to control
  if (!is.null(controls)) {
    nX <- normalizeToControls(X, pheno, controls)
    xlist[["normalized_to_control"]] <- nX
  }

  ## ComBat ------------------------------------------------------
  mod1 <- model.matrix(~pheno)
  bX <- sva::ComBat(X, batch = batch, mod = mod1)
  xlist[["ComBat"]] <- bX

  ## limma -------------------------------------------------------
  cX <- limma::removeBatchEffect(X, batch = batch, mod = mod1)
  xlist[["limma"]] <- cX

  ## MNN ---------------------------------------------------------
  xlist[["fastMNN"]] <- MNNcorrect(X, batch)

  ## restricted MNN ----------------------------------------------
  if (!is.null(controls)) {
    xlist[["restrictedMNN"]] <- MNNcorrect(X, batch, controls = controls)
  }

  ## Harmony ----------------------------------------------------
  res <- try(runHarmony(X, batch = batch))
  if (!"try-error" %in% class(res)) {
    xlist[["Harmony"]] <- as.matrix(res$corrected)
  }

  ## RUV and SVA
  xlist[["RUV3"]] <- ruvCorrect(X, pheno, k = NULL, type = "III")
  xlist[["RUVg"]] <- ruvCorrect(X, pheno, k = NULL, type = "g")
  xlist[["SVA"]] <- svaCorrect(X, pheno)

  ## ComBatX
  xlist[["ComBatX.a"]] <- ComBatX(X, batch = batch, pheno = NULL, controls = NULL)
  xlist[["ComBatX.b"]] <- ComBatX(X, batch = batch, pheno = pheno, controls = NULL, bc.dim = 0)
  xlist[["ComBatX.c"]] <- ComBatX(X, batch = batch, pheno = pheno, controls = NULL, bc.dim = -1)
  ##  xlist[['ComBatX.d']] <- ComBatX(X, batch=batch, pheno=pheno, controls=NULL, bc.dim=10)
  xlist[["ComBatX.e"]] <- ComBatX(X, batch = batch, pheno = NULL, controls = controls)
  xlist[["ComBatX.f"]] <- ComBatX(X, batch = batch, pheno = pheno, controls = controls, bc.dim = 0)
  xlist[["ComBatX.g"]] <- ComBatX(X, batch = batch, pheno = pheno, controls = controls, bc.dim = -1)

  if (nlevel <= 4) {
    xlist[["NNM"]] <- gx.nnmcorrect(X, pheno)$X
    if (combinations) xlist[["CB+NNM"]] <- gx.nnmcorrect(bX, pheno)$X
  } else {
    message("warning: too many levels. skipping NNM.")
  }

  if (combinations) {
    xlist[["RUV3+SVA"]] <- svaCorrect(xlist[["RUV3"]], pheno)
    xlist[["SVA+RUV3"]] <- ruvCorrect(xlist[["SVA"]], pheno)

    xlist[["CB+SVA"]] <- svaCorrect(bX, pheno)
    xlist[["CB+RUV3"]] <- ruvCorrect(bX, pheno)

    xlist[["CB+SVA+RUV3"]] <- ruvCorrect(xlist[["CB+SVA"]], pheno)
    xlist[["CB+RUV3+SVA"]] <- svaCorrect(xlist[["CB+RUV3"]], pheno)

    zX2 <- ruvCorrect(xlist[["ComBatX.c"]], pheno)
    xlist[["ComBatX.c+RUV3+SVA"]] <- svaCorrect(zX2, pheno)
  }

  xlist
}


bc.evaluateResults <- function(xlist, pheno, ref, pos = NULL) {
  ## determine strongest signal
  ct <- makeDirectContrasts(data.frame(pheno), ref = ref)
  Y <- ct$exp.matrix
  Y[Y == 0] <- NA
  colnames(Y)

  if (ncol(Y) > 1) {
    i <- 14
    nsig <- c()
    for (i in 1:ncol(Y)) {
      res <- stats.numsig(xlist[[1]], y = Y[, i], ref = -1, lfc = 0.1, q = 0.2)
      nsig <- rbind(nsig, c(length(res$genes), length(res$gsets)))
    }
    rownames(nsig) <- colnames(Y)
    nsig
    nsig <- nsig[order(-exp(rowMeans(log(nsig)))), ]
    head(nsig)
    colnames(Y)
    strongest.signal <- rownames(nsig)[1]
    strongest.signal
    y <- Y[, strongest.signal]
  } else {
    y <- Y[, 1]
  }

  ## compute and make table
  numsig <- lapply(xlist, stats.numsig, y = y, ref = -1, lfc = 0.2, q = 0.05)
  res <- t(sapply(numsig, function(r) {
    c(sapply(r[1:2], length), avg.fc = mean(abs(r[[3]])))
  }))

  str(numsig)
  ## r1 <- sapply( numsig, function(s) mean(numsig[[1]]$genes %in% s$genes))
  ## r2 <- sapply( numsig, function(s) mean(numsig[[1]]$gsets %in% s$gsets))
  ## s1 <- sapply( numsig, function(s) mean(s$genes %in% numsig[[1]]$genes))
  ## s2 <- sapply( numsig, function(s) mean(s$gsets %in% numsig[[1]]$gsets))
  n1 <- sapply(numsig, function(s) length(intersect(s$genes, numsig[[1]]$genes)))
  n2 <- sapply(numsig, function(s) length(union(s$genes, numsig[[1]]$genes)))
  m1 <- sapply(numsig, function(s) length(intersect(s$gsets, numsig[[1]]$gsets)))
  m2 <- sapply(numsig, function(s) length(union(s$gsets, numsig[[1]]$gsets)))
  ##  res <- cbind(res, r.genes=r1, r.gsets=r2, s.genes=s1, s.gsets=s2)
  res <- cbind(res, r.genes = n1 / n2, r.gsets = m1 / m2)
  res.score <- exp(rowMeans(log(t(t(res) / res[1, ]))))

  if (is.null(pos)) {
    CLUSTFUN <- function(x) uwot::tumap(scale(t(x), scale = FALSE))
    pos <- lapply(xlist, function(x) CLUSTFUN(x))
  }
  sil_score <- sapply(pos, function(p) {
    score <- cluster::silhouette(as.integer(factor(pheno)), stats::dist(p))
    mean(score[, "sil_width"])
  })
  sil_score1 <- pmax(sil_score, 1e-4)
  better <- res.score * (sil_score1 / sil_score1[1])**0.5

  res1 <- cbind(res, sil_score, better)
  res1 <- res1[order(-res1[, "better"]), ]
  res1
}





## =====================================================================================
## =========================== END OF FILE =============================================
## =====================================================================================
