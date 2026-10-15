##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##


## ----------------------------------------------------------------------
## Biomarker selection: data preparation and method registry
## ----------------------------------------------------------------------

# Phenotype, sample filter, imputation, feature filter and top-SD cut shared by
# pgx.compute_importance() and the biomarker pipeline. Returns NULL (with a
# message) when too few classes, samples or features remain.
.biomarker_prepare <- function(pgx, pheno, level = "genes",
                               filter_features = NULL,
                               select_features = NULL,
                               select_samples = NULL,
                               nfeatures = 50) {
  ft <- ifelse(is.null(filter_features), "<all>", filter_features)
  sel <- select_features

  if (!(pheno %in% colnames(pgx$samples))) {
    message("ERROR. pheno not in pgx$samples")
    return(NULL)
  }

  ## WARNING. this converts any phenotype to discrete
  ## A phenotype is only ever read next to an expression profile -- `X` is
  ## subset by names(y) below -- so `y` is defined on pgx$X's samples. Outlier
  ## removal may have cut those below the upload pgx$samples spans (D-24).
  ss <- intersect(rownames(pgx$samples), colnames(pgx$X))
  y <- as.character(pgx$samples[ss, pheno])
  names(y) <- ss
  if (!is.null(select_samples)) {
    y <- y[names(y) %in% select_samples]
  }
  y <- y[!is.na(y)]

  if (length(unique(y)) < 2) {
    message("[pgx.compute_importance] Less than 2 phenotype levels after sample filtering. Returning NULL.")
    return(NULL)
  }
  if (length(y) < 4) {
    message("[pgx.compute_importance] Not enough samples after filtering (n=", length(y), "). Returning NULL.")
    return(NULL)
  }

  ## -------------------------------------------
  ## select features
  ## -------------------------------------------
  if (FALSE && level == "geneset") {
    X <- pgx$gsetX
    if (any(is.na(X))) {
      X <- X[complete.cases(X), , drop = FALSE]
    }
  } else {
    X <- pgx$X
    if (any(is.na(X))) {
      X <- .pgx_impute_svd2(X)
    }
  }

  ## ----------- filter with selected features
  info("[pgx.compute_importance] Filtering features...")

  is.family <- (ft %in% c(names(pgx$families), names(playdata::iGSETS)))
  if (ft == "<custom>" && !is.null(sel) && length(sel) > 0) {
    ## ------------- filter with user selection
    if (sel[1] != "") {
      features <- playbase::filterProbes(pgx$genes, sel)
      pp <- which(rownames(X) %in% features)
      X <- X[pp, , drop = FALSE]
    }
  } else if (is.family) {
    pp <- rownames(X)
    if (ft == "<all>") {
      pp <- rownames(X)
    } else if (ft %in% names(pgx$families)) {
      gg <- pgx$families[[ft]]
      pp <- filterProbes(pgx$genes, gg)
    } else if (ft %in% names(playdata::iGSETS)) {
      gg <- unlist(playdata::getGSETS(ft))
      pp <- filterProbes(pgx$genes, gg)
    }
    pp <- intersect(pp, rownames(X))
    X <- X[pp, , drop = FALSE]
  }

  if (nrow(X) == 0) {
    message("[pgx.compute_importance] No features remaining after filtering. Returning NULL.")
    return(NULL)
  }

  ## ----------- restrict to top SD -----------
  if (pgx$datatype == "multi-omics") {
    X <- mofa.topSD(X, 10 * nfeatures)
  } else {
    sdx <- matrixStats::rowSds(X, na.rm = TRUE)
    X <- head(X[order(-sdx), , drop = FALSE], 10 * nfeatures)
  }

  list(X = X, y = y, class_counts = c(table(y)))
}


## ----------------------------------------------------------------------
## Method registry
## ----------------------------------------------------------------------
##
## The one place that says what a selection method does. Every entry has
##
##   importance(X, y)            non-negative score per feature, named by
##                               rownames(X); higher is more important
##   fit_predict(Xtr, ytr, Xte)  class probabilities, ncol(Xte) x levels(ytr),
##                               rows summing to 1
##
## X matrices are features x samples and y is a class vector (factor or
## character). "model" methods refit their own classifier; "filter" methods
## only score features, so they predict with a ridge logistic regression.
## Nothing here seeds the RNG: the caller does, and every draw comes from R's.

BIOMARKER_METHODS <- list(
  xgboost = list(
    label = "XGBoost",
    kind = "model",
    importance = function(X, y) {
      bst <- .biomarker_xgboost(X, as.factor(y))
      gain <- xgboost::xgb.importance(model = bst)
      score <- gain$Gain[match(.biomarker_xgboost_names(X), gain$Feature)]
      score[is.na(score)] <- 0 ## features the trees never split on
      stats::setNames(score, rownames(X))
    },
    fit_predict = function(Xtr, ytr, Xte) {
      ytr <- as.factor(ytr)
      bst <- .biomarker_xgboost(Xtr, ytr)
      p <- stats::predict(bst, .biomarker_xgboost_data(Xte))
      ## binary gives P(second level); multiclass a samples x classes matrix
      ## (a row-major vector before xgboost 2) of single precision sums
      if (nlevels(ytr) == 2) {
        p <- cbind(1 - p, p)
      } else {
        p <- matrix(p, ncol = nlevels(ytr), byrow = is.null(dim(p)))
        p <- p / rowSums(p)
      }
      dimnames(p) <- list(colnames(Xte), levels(ytr))
      p
    }
  ),
  elastic_net = list(
    label = "Elastic net",
    kind = "model",
    importance = function(X, y) {
      net <- .biomarker_glmnet(X, as.factor(y), alpha = 0.5)
      beta <- stats::coef(net$fit, s = net$s)
      if (is.list(beta)) beta <- do.call(cbind, beta) ## multinomial: one per class
      beta <- as.matrix(beta)[-1, , drop = FALSE]
      ## scale by SD so the coefficients compare across features
      rowMeans(abs(beta)) * matrixStats::rowSds(X)
    },
    fit_predict = function(Xtr, ytr, Xte) {
      .biomarker_glmnet_predict(Xtr, ytr, Xte, alpha = 0.5)
    }
  ),
  random_forest = list(
    label = "Random forest",
    kind = "model",
    importance = function(X, y) {
      fit <- randomForest::randomForest(t(X), as.factor(y), ntree = 500)
      fit$importance[, "MeanDecreaseGini"]
    },
    fit_predict = function(Xtr, ytr, Xte) {
      fit <- randomForest::randomForest(t(Xtr), as.factor(ytr), ntree = 500)
      p <- stats::predict(fit, t(Xte), type = "prob")
      matrix(p, nrow = ncol(Xte), dimnames = list(colnames(Xte), colnames(p)))
    }
  ),
  splsda = list(
    label = "sPLS-DA",
    kind = "model",
    importance = function(X, y) {
      fit <- .biomarker_splsda(X, as.factor(y))
      rowMeans(fit$loadings$X[rownames(X), , drop = FALSE]^2)
    },
    fit_predict = function(Xtr, ytr, Xte) {
      ytr <- as.factor(ytr)
      fit <- .biomarker_splsda(Xtr, ytr)
      ## predicted class indicators on the last component; they are not
      ## bounded, so clip at zero and renormalise
      p <- stats::predict(fit, t(Xte))$predict[, , fit$ncomp]
      p <- pmax(matrix(p, nrow = ncol(Xte)), 0)
      p[rowSums(p) == 0, ] <- 1
      p <- p / rowSums(p)
      dimnames(p) <- list(colnames(Xte), levels(ytr))
      p
    }
  ),
  correlation = list(
    label = "Correlation",
    kind = "filter",
    importance = function(X, y) {
      dy <- stats::model.matrix(~ 0 + y, data.frame(y = as.factor(y)))
      rho <- suppressWarnings(stats::cor(t(X), dy)) ## constant features give NA
      rho[is.na(rho)] <- 0
      sqrt(rowMeans(rho^2))
    },
    fit_predict = function(Xtr, ytr, Xte) {
      .biomarker_glmnet_predict(Xtr, ytr, Xte, alpha = 0)
    }
  ),
  ftest = list(
    label = "F-test",
    kind = "filter",
    importance = function(X, y) {
      res <- gx.limmaF(X, as.character(y), lfc = 0, fdr = 1, compute.means = FALSE, verbose = 0)
      stats::setNames(-log10(res[rownames(X), "P.Value"]), rownames(X))
    },
    fit_predict = function(Xtr, ytr, Xte) {
      .biomarker_glmnet_predict(Xtr, ytr, Xte, alpha = 0)
    }
  )
)

# Stratified fold ids: each class is dealt over the k folds in random order.
.biomarker_folds <- function(y, k) {
  folds <- integer(length(y))
  for (idx in split(seq_along(y), y)) {
    folds[idx] <- rep_len(seq_len(k), length(idx))[sample.int(length(idx))]
  }
  folds
}

# glmnet with lambda.min from a stratified inner CV. With 2 samples in a class
# no inner split leaves glmnet the 2 it needs, so the path's least penalised
# lambda is used instead.
.biomarker_glmnet <- function(X, y, alpha) {
  family <- if (nlevels(y) > 2) "multinomial" else "binomial"
  nfolds <- min(5L, min(table(y)))
  ## glmnet warns below 8 samples per class, which is the norm for biomarkers
  withCallingHandlers(
    if (nfolds < 3) {
      fit <- glmnet::glmnet(t(X), y, family = family, alpha = alpha)
      list(fit = fit, s = min(fit$lambda))
    } else {
      ## grouped = FALSE: folds hold too few samples for per-fold deviance
      fit <- glmnet::cv.glmnet(t(X), y,
        family = family, alpha = alpha,
        foldid = .biomarker_folds(y, nfolds), grouped = FALSE
      )
      list(fit = fit, s = fit$lambda.min)
    },
    warning = function(w) {
      if (grepl("fewer than 8", conditionMessage(w))) invokeRestart("muffleWarning")
    }
  )
}

# Class probabilities from a .biomarker_glmnet() fit on the training samples.
.biomarker_glmnet_predict <- function(Xtr, ytr, Xte, alpha) {
  ytr <- as.factor(ytr)
  net <- .biomarker_glmnet(Xtr, ytr, alpha)
  p <- stats::predict(net$fit, t(Xte), s = net$s, type = "response")
  ## binomial gives P(second level); multinomial a samples x classes x 1 array
  if (nlevels(ytr) == 2) p <- cbind(1 - p, p)
  matrix(p, nrow = ncol(Xte), dimnames = list(colnames(Xte), levels(ytr)))
}

# xgboost input with positional feature names, as xgboost rejects some
# characters found in feature ids.
.biomarker_xgboost_names <- function(X) paste0("f", seq_len(nrow(X)))

.biomarker_xgboost_data <- function(X, y = NULL) {
  x <- t(X)
  colnames(x) <- .biomarker_xgboost_names(X)
  xgboost::xgb.DMatrix(x, label = if (!is.null(y)) as.integer(y) - 1L)
}

# Boosted trees; one thread and no subsampling keep the fit deterministic.
.biomarker_xgboost <- function(X, y) {
  params <- list(
    booster = "gbtree", eta = 0.1, max_depth = 3, subsample = 1, nthread = 1,
    objective = if (nlevels(y) > 2) "multi:softprob" else "binary:logistic"
  )
  if (nlevels(y) > 2) params$num_class <- nlevels(y)
  xgboost::xgb.train(params, .biomarker_xgboost_data(X, y), nrounds = 100, verbose = 0)
}

# Sparse PLS-DA on one component per class contrast, keeping at most 25
# features per component.
.biomarker_splsda <- function(X, y) {
  ncomp <- nlevels(y) - 1L
  mixOmics::splsda(t(X), y, ncomp = ncomp, keepX = rep(min(25L, nrow(X)), ncomp))
}
