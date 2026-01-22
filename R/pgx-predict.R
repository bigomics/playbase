##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


## ----------------------------------------------------------------------
## Variable importance functions
## ----------------------------------------------------------------------

#' Main compute driver function. Compute variable importance for one
#' specific contrast.
#'
#' @export
pgx.compute_importance <- function(pgx, pheno, level = "genes",
                                   filter_features = NULL,
                                   select_features = NULL,
                                   select_samples = NULL,
                                   nfeatures = 60,
                                   multiomics = 2,
                                   niter = 1) {
  if (0) {
    level <- "genes"
    filter_features <- NULL
    select_features <- NULL
    select_samples <- NULL
    nfeatures <- 60
    multiomics <- 2
  }

  ft <- ifelse(is.null(filter_features), "<all>", filter_features)
  sel <- select_features

  if (!(pheno %in% colnames(pgx$samples))) {
    message("ERROR. pheno not in pgx$samples")
    return(NULL)
  }

  ## WARNING. this converts any phenotype to discrete
  y <- as.character(pgx$samples[, pheno])
  names(y) <- rownames(pgx$samples)
  if (!is.null(select_samples)) {
    y <- y[names(y) %in% select_samples]
  }
  y <- y[!is.na(y)]

  ## -------------------------------------------
  ## select features
  ## -------------------------------------------
  if (FALSE && level == "geneset") {
    X <- pgx$gsetX ## NB: this will augment
    if (any(is.na(X))) {
      X <- X[complete.cases(X), , drop = FALSE]
    }
  } else {
    X <- pgx$X ## NB: this will augment
    is.mox <- is.multiomics(rownames(X))
    if (any(is.na(X))) {
      if (is.mox) {
        X <- imputeMissing.mox(X, method = "SVD2")
      } else {
        X <- imputeMissing(X, method = "SVD2")
      }
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

  ## ----------- restrict to top SD -----------
  is.multiomics <- (pgx$datatype == "multi-omics")
  if (is.multiomics) {
    X <- mofa.topSD(X, 10 * nfeatures) ## top 100
  } else {
    sdx <- matrixStats::rowSds(X, na.rm = TRUE)
    X <- head(X[order(-sdx), , drop = FALSE], 10 * nfeatures) ## top 100
  }

  ## ----------------------------------------
  ## augment
  ## ----------------------------------------
  ## augment to at least 100 samples per level :)
  ii <- tapply(1:length(y), y, function(ii) {
    sample(c(ii, ii), size = 100, replace = TRUE)
  })
  y <- y[unlist(ii)]
  X <- X[, names(y), drop = FALSE]

  ## -------------------------------------------
  ## compute importance values
  ## -------------------------------------------

  ## determine is dataset is multi-omics
  has.mofa <- "mofa" %in% names(pgx)
  is.multiomics <- (pgx$datatype == "multi-omics" && has.mofa)

  P <- NULL
  if (is.multiomics && multiomics > 0) {
    ## compute variable importance for MULTI-OMICS. We need not to
    ## use the augmented data.
    message("Computing multi-omics MOFA variable importance")
    kernels <- c(
      "mofa", "pca", "nmf", "nmf2", "mcia", "wgcna", "diablo", "rgcca",
      "rgcca,rgcca", "rgcca.rgccda", "rgcca.mcoa"
    )
    P <- pgx.compute_mofa_importance(
      pgx, pheno,
      numfactors = 8, use.sdwt = TRUE, kernels = kernels
    )

    if (multiomics == 2) {
      ## compute variable importance for 2nd pass using top scoring
      ## from mofa_importance.
      sel <- head(intersect(rownames(P), rownames(X)), 4 * nfeatures) ## TUNE TOP
      X <- X[sel, ]
    } else {
      ## Only single pass with MOFA is not using augmented data.
      X <- pgx$X
      y <- pgx$samples[, pheno]
      names(y) <- rownames(pgx$samples)
    }
  }

  if (!is.multiomics || multiomics != 1) {
    ## compute variable importance using ML methods. We use the
    ## augmented data.

    ## add some noise
    sdx0 <- matrixStats::rowSds(X, na.rm = TRUE)
    sdx1 <- 0.5 * sdx0 + 0.5 * mean(sdx0, na.rm = TRUE)

    methods <- c(
      "glmnet", "randomforest", "xgboost", "splsda",
      "correlation", "ftest"
    )
    message("Computing ML variable importance...")
    P <- 0
    i <- 1
    for (i in 1:niter) {
      X1 <- X + 0.25 * sdx1 * matrix(rnorm(length(X)), nrow(X), ncol(X))
      y1 <- y
      names(y1) <- colnames(X1) <- paste0("x", 1:ncol(X))
      res <- pgx.variableImportance(
        X1,
        y1,
        methods = methods,
        reduce = 1000,
        resample = 0,
        scale = FALSE,
        add.noise = 0
      )
      P <- P + res$importance
    }
    P <- P / niter
    remove(X1, y1)
  }

  ## normalize importance measures
  P <- abs(P) ## sometimes negative according to sign
  P[is.na(P)] <- 0
  P[is.nan(P)] <- 0

  ## Convert to elevated rank. take top features
  R <- P
  if (nrow(R) > 1) {
    R <- (apply(P, 2, rank) / nrow(P))**4 ## exponent???
    R <- R[order(-rowSums(R, na.rm = TRUE)), , drop = FALSE]
  }
  sel <- head(rownames(R), nfeatures) ## top features
  R <- R[sel, , drop = FALSE]
  X <- X[sel, , drop = FALSE]

  ## reduce R and y (from augmented set)
  kk <- names(y)[which(!duplicated(names(y)))]
  X <- X[, kk, drop = FALSE]
  y <- y[kk]

  ## make partition tree
  rf <- makePartitionTree(X, y, add.splits = 0)

  ## rf = NULL
  res <- list(R = R, y = y, X = X, rf = rf)

  return(res)
}


#' @title Variable importance for survival models
#'
#' @param X Numeric matrix of predictor variables, rows are variables
#' @param time Vector of event/censoring times
#' @param status Vector indicating event (1) or censoring (0)
#' @param methods Methods for computing variable importance. Options are "glmnet", "randomforest", "boruta", "xgboost", "pls".
#'
#' @return Named list with variable importance scores by method.
#'
#' @description Compute variable importance scores for predictors in a survival analysis model using different methods.
#'
#' @details This function calculates variable importance scores using a variety of methods suitable for survival analysis.
#' The input data consists of a predictor matrix \code{X}, a vector of event/censoring times \code{time} and a status indicator vector \code{status}.
#'
#' The following methods can be selected via the \code{methods} parameter:
#' \itemize{
#' \item glmnet: Absolute value of coefficients from elastic net Cox model
#' \item randomForest: Variable importance from random survival forest
#' \item boruta: Boruta variable selection algorithm
#' \item xgboost: Importance scores from XGBoost survival model
#' \item pls: Absolute coefficients from partial least squares Cox model
#' }
#'
#' Variable importance scores are returned for each method in a named list.
#' These scores can be used to select important predictors for survival modeling.
#'
#' @export
pgx.survivalVariableImportance <- function(X, time, status,
                                           methods = c("glmnet", "randomforest", "boruta", "xgboost", "pls")) {
  ## ----------------------------------------------------------------------
  ## multi-class version
  ## ----------------------------------------------------------------------

  imp <- list()
  xnames <- rownames(X)
  ## sdx <- apply(X, 1, stats::sd, na.rm = TRUE)
  sdx <- matrixStats::rowSds(X, na.rm = TRUE)
  if (nrow(X) == 1) X <- rbind(X, X)

  if (!inherits(status, "logical") && all(status %in% c(0, 1, NA))) {
    stop("status must be logical or 0/1")
  }

  y <- survival::Surv(time + 0.0001, status)

  if ("glmnet" %in% methods) {
    fam <- "cox"

    NFOLD <- 5
    out0 <- glmnet::cv.glmnet(t(X), y, alpha = 0, family = fam, standardize = TRUE, nfolds = NFOLD)
    cf0 <- glmnet::coef.glmnet(out0, s = "lambda.min")[, 1]
    out1 <- glmnet::cv.glmnet(t(X), y, alpha = 1, family = fam, standardize = TRUE, nfolds = NFOLD)
    cf1 <- glmnet::coef.glmnet(out1, s = "lambda.min")[, 1]
    out0a <- glmnet::cv.glmnet(t(X), y, alpha = 0, family = fam, standardize = FALSE, nfolds = NFOLD)
    cf0a <- glmnet::coef.glmnet(out0a, s = "lambda.min")[, 1]
    out1a <- glmnet::cv.glmnet(t(X), y, alpha = 1, family = fam, standardize = FALSE, nfolds = NFOLD)
    cf1a <- glmnet::coef.glmnet(out1a, s = "lambda.min")[, 1]
    imp[["coxnet.a0"]] <- (cf0 / max(abs(cf0)) + cf0a / max(abs(cf0a))) * sdx[names(cf0)]
    imp[["coxnet.a1"]] <- (cf1 / max(abs(cf1)) + cf1a / max(abs(cf1a))) * sdx[names(cf1)]
  }

  if ("randomforest" %in% methods) {
    df <- data.frame(time = time, status = status, t(X))
    fit_rf <- randomForestSRC::rfsrc(survival::Surv(time, status) ~ ., data = df)
    vimp <- randomForestSRC::vimp(fit_rf)$importance
    imp[["randomForest"]] <- vimp
  }

  if ("boruta" %in% methods) {
    imp4 <- rep(0, nrow(X))
    niter <- 4
    for (k in 1:niter) {
      jj <- sample(ncol(X), ncol(X) * 0.9)
      out3 <- Boruta(t(X[, jj, drop = FALSE]), y[jj])
      fd <- factor(out3$finalDecision, levels = c("Rejected", "Tentative", "Confirmed"))
      fd <- (as.integer(fd) - 1) / 2
      imp4 <- imp4 + fd / niter
    }

    names(imp4) <- rownames(X)
    imp[["Boruta"]] <- imp4
  }

  if ("xgboost" %in% methods) {
    yy <- ifelse(!status, -time, time)
    bst <- xgboost::xgboost(
      data = t(X), label = yy, booster = "gbtree",
      max_depth = 2, eta = 1, nthread = 2, nrounds = 2,
      verbose = 0, objective = "survival:cox"
    )

    xgmat <- xgboost::xgb.importance(model = bst)
    imp5 <- xgmat$Gain**0.33
    imp5 <- imp5[match(rownames(X), xgmat$Feature)]
    names(imp5) <- rownames(X)
    imp5[which(is.na(imp5))] <- 0
    Matrix::tail(sort(imp5))
    imp[["xgboost"]] <- imp5

    ## linear model
    bst2 <- xgboost::xgboost(
      data = t(X), label = yy, booster = "gblinear",
      max_depth = 2, eta = 1, nthread = 2, nrounds = 2,
      verbose = 0, objective = "survival:cox"
    )
    xgmat <- xgboost::xgb.importance(model = bst2)
    imp6 <- xgmat$Weight
    names(imp6) <- xgmat$Feature
    imp6 <- imp6[match(rownames(X), names(imp6))]
    Matrix::tail(sort(imp6))
    imp[["xgboost.lin"]] <- imp6
  }

  if ("pls" %in% methods) {
    res <- plsRcox::plsRcox(t(X), time = time, event = status, nt = 5)
    cf <- res$Coeffs[, 1]
    cf[is.na(cf)] <- 0
    cf <- cf * sdx[names(cf)] ## really?
    imp[["pls.cox"]] <- abs(cf) / max(abs(cf), na.rm = TRUE)
  }

  P <- do.call(cbind, imp)
  P <- abs(P) ## always positive??
  P[is.na(P)] <- 0
  P <- P[xnames, , drop = FALSE]
  return(P)
}


CARET.METHODS <- c("svmLinear", "rpart", "glm", "glmStepAIC", "pda", "sparseLDA", "pls", "plsda")
CARET.METHODS <- c("pls", "plsda")

## for compatibility before v3.3.0
#' @export
pgx.multiclassVariableImportance <- function(X, y, methods) {
  res <- pgx.variableImportance(X = X, y = y, methods = methods)
  res$importance
}

#' @describeIn pgx.survivalVariableImportance Calculates variable
#'   importance scores for predictors of a multiclass response using
#'   various methods.
#'
#' @param y Multiclass factor response variable. Contains the class
#'   labels for each sample
#' @export
pgx.variableImportance <- function(X, y,
                                   methods = c(
                                     "glmnet", "randomforest",
                                     "xgboost", "splsda", "correlation", "ftest",
                                     "boruta", CARET.METHODS
                                   )[1:6],
                                   scale = TRUE, reduce = 1000, resample = 0,
                                   add.noise = 0) {
  ## variables
  imp <- list()
  runtime <- list()

  if (nrow(X) == 1) X <- rbind(X, X)
  xnames <- rownames(X)

  ## drop missing??
  sel <- which(!is.na(y) & y != "")
  X <- X[, sel, drop = FALSE]
  y <- y[sel]

  ## convert to factor
  y <- factor(as.character(y))

  ## resample to minimum size to balance groups
  if (resample > 0) {
    NSIZE <- resample
    jj <- tapply(1:length(y), y, function(ii) {
      if (length(ii) < NSIZE) ii <- sample(c(ii, ii), NSIZE, replace = TRUE)
      return(ii)
    })
    jj <- unlist(jj)
    X <- X[, jj]
    y <- y[jj]
  }


  ## add noise
  if (add.noise > 0) {
    sdx0 <- matrixStats::rowSds(X, na.rm = TRUE)
    sdx1 <- add.noise * (0.5 * sdx0 + 0.5 * mean(sdx0, na.rm = TRUE))
    X <- X + sdx1 * matrix(rnorm(length(X)), nrow(X), ncol(X)) ## add some noise
  }

  if (scale) {
    sdx <- matrixStats::rowSds(X, na.rm = TRUE)
    X <- (X - rowMeans(X, na.rm = TRUE)) / (sdx + 1e-8)
  }

  ## reduce features using fast correlation/cov
  if (reduce > 0 && reduce < nrow(X)) {
    dy <- model.matrix(~ 0 + y)
    rho <- cov(t(X), dy) ## note: using covariance
    imp1 <- sqrt(rowMeans(rho**2))
    sel <- head(names(sort(-imp1)), reduce)
    X <- X[sel, ]
  }

  if ("correlation" %in% methods) {
    runtime[["correlation"]] <- system.time({
      dy <- model.matrix(~ 0 + y)
      rho <- cov(t(X), dy)
    })
    imp[["correlation"]] <- sqrt(rowMeans(rho**2))
  }

  if ("ftest" %in% methods) {
    runtime[["ftest"]] <- system.time({
      res <- gx.limma(X, y, lfc = 0, fdr = 1, f.test = TRUE)
    })
    imp1 <- -log10(res[rownames(X), "P.Value"])
    names(imp1) <- rownames(X)
    imp[["ftest"]] <- imp1
  }

  if ("glmnet" %in% methods) {
    fam <- "multinomial"
    NFOLD <- 5
    runtime[["glmnet.a0"]] <- system.time({
      out0 <- glmnet::cv.glmnet(t(X), y, alpha = 0, family = fam, standardize = TRUE, nfold = NFOLD)
      cf0 <- Matrix::rowMeans(do.call(cbind, coef(out0, s = "lambda.min"))[-1, ]**2)
      out0a <- glmnet::cv.glmnet(t(X), y, alpha = 0, family = fam, standardize = FALSE, nfold = NFOLD)
      cf0a <- Matrix::rowMeans(do.call(cbind, coef(out0a, s = "lambda.min"))[-1, ]**2)
    })

    runtime[["glmnet.a1"]] <- system.time({
      out1 <- glmnet::cv.glmnet(t(X), y, alpha = 0.99, family = fam, standardize = TRUE, nfold = NFOLD)
      cf1 <- Matrix::rowMeans(do.call(cbind, coef(out1, s = "lambda.min"))[-1, ]**2)
      out1a <- glmnet::cv.glmnet(t(X), y, alpha = 0.99, family = fam, standardize = FALSE, nfold = NFOLD)
      cf1a <- Matrix::rowMeans(do.call(cbind, coef(out1a, s = "lambda.min"))[-1, ]**2)
    })

    sdx <- matrixStats::rowSds(X, na.rm = TRUE)
    imp[["glmnet.a0"]] <- (cf0 / max(abs(cf0)) + cf0a / max(abs(cf0a))) * sdx[names(cf0)]
    imp[["glmnet.a1"]] <- (cf1 / max(abs(cf1)) + cf1a / max(abs(cf1a))) * sdx[names(cf1)]
  }

  if ("randomforest" %in% methods) {
    runtime[["randomForest"]] <- system.time({
      fit_rf <- randomForest::randomForest(t(X), factor(y))
    })
    imp[["randomForest"]] <- fit_rf$importance[, 1]
  }

  if ("boruta" %in% methods) {
    imp4 <- rep(0, nrow(X))
    niter <- 4
    runtime[["Boruta"]] <- system.time({
      for (k in 1:niter) {
        jj <- sample(ncol(X), ncol(X) * 0.9)
        out3 <- Boruta::Boruta(t(X[, jj, drop = FALSE]), y[jj])
        fd <- factor(out3$finalDecision,
          levels = c("Rejected", "Tentative", "Confirmed")
        )
        fd <- (as.integer(fd) - 1) / 2
        imp4 <- imp4 + fd / niter
      }
    })
    names(imp4) <- rownames(X)
    imp[["Boruta"]] <- imp4
  }
  
  if ("xgboost" %in% methods) {

    ny <- length(table(y))
    yy <- as.character(y)
    obj <- "binary:logistic"
    if (length(unique(yy)) > 2) obj <- "multi:softprob" 

    runtime[["xgboost"]] <- system.time({
      bst <- xgboost::xgboost(
        x = t(X), y = yy, booster = "gbtree",
        max_depth = 2, nthread = 2,
        nrounds = 2, objective = obj
      )
    })

    xgmat <- xgboost::xgb.importance(model = bst)
    imp5 <- xgmat$Gain**0.2
    imp5 <- imp5[match(rownames(X), xgmat$Feature)]
    names(imp5) <- rownames(X)
    imp[["xgboost"]] <- imp5

    runtime[["xgboost.lin"]] <- system.time({
      bst2 <- xgboost::xgboost(
        x = t(X), y = yy, booster = "gblinear",
        max_depth = 2, nthread = 2,
        nrounds = 2, objective = obj
      )
    })

    xgmat <- xgboost::xgb.importance(model = bst2)
    imp6 <- xgmat$Weight
    names(imp6) <- xgmat$Feature
    imp6 <- imp6[match(rownames(X), names(imp6))]
    imp[["xgboost.lin"]] <- imp6

  }

  if ("splsda" %in% methods) {
    n <- min(25, nrow(X))
    colnames(X) <- names(y) <- paste0("sample", 1:length(y))
    rownames(X) <- make.unique(rownames(X))
    runtime[["splsda"]] <- system.time({
      res <- mixOmics::splsda(t(X), y, keepX = c(n, n))
    })
    impx <- rowMeans(res$loadings$X[rownames(X), ]**2)
    imp[["splsda"]] <- impx
  }

  #  CARET.METHODS = c("glmnet",'svmLinear','rpart','xgbLinear','glm','glmStepAIC',
  #  'pda','sparseLDA','pls','plsda')
  caret.methods <- intersect(methods, CARET.METHODS)
  if (length(caret.methods)) {
    repeatedcv <- caret::trainControl(method = "repeatedcv", number = 25, repeats = 1)
    m <- "pls"
    for (m in caret.methods) {
      message("computing variable importance for ", m)
      t0 <- NULL
      if (m == "plsda") {
        runtime[["plsda"]] <- system.time({
          fit <- try(caret::plsda(x = t(X), y = y, ncomp = min(20, nrow(X) / 3)))
        })
        if (!"try-error" %in% class(fit)) {
          imp[[m]] <- sqrt(rowMeans(caret::varImp(fit)**2))
        }
      } else {
        runtime[[m]] <- system.time({
          fit <- try(caret::train(
            x = t(X), y = y, method = m,
            trControl = repeatedcv, trace = FALSE
          ))
        })
        if (!"try-error" %in% class(fit)) {
          imp[[m]] <- sqrt(rowMeans(caret::varImp(fit)$importance**2))
        }
      }
    } ## for caret.methods
  }
  names(imp)

  imp <- imp[!sapply(imp, is.null)]
  imp <- lapply(imp, function(v) v[match(xnames, names(v))])
  P <- do.call(cbind, imp)
  colnames(P) <- names(imp)
  rownames(P) <- xnames
  P <- abs(P) ## always positive??
  P[is.na(P)] <- 0
  ## P <- P[xnames, , drop = FALSE]

  ## normalize values by ranking
  ii <- which(P == 0)
  P <- 1 / apply(-P, 2, rank)
  P[ii] <- 0

  ## parse runtimes
  tt <- unlist(sapply(runtime[colnames(P)], function(t) t[1]))
  names(tt) <- sub(".user.self", "", names(tt))
  tt <- tt[match(colnames(P), names(tt))]
  names(tt) <- colnames(P)
  tt
  sort(tt)

  res <- list(importance = P, runtime = tt)
  return(res)
}


#' @export
plotImportance <- function(P, p.sign = NULL, top = 50, runtime = NULL) {
  P[is.na(P)] <- 0
  if (is.null(p.sign)) p.sign <- rep(1, nrow(P))
  p.sign <- sign(p.sign)
  P1 <- 1 / apply(-abs(P), 2, rank) * p.sign
  P1 <- P1[order(-rowMeans(P1, na.rm = TRUE)), ]
  if (any(p.sign < 0)) {
    sel <- unique(head(rownames(P1), top / 2), tail(rownames(P1), top / 2))
  } else {
    sel <- head(rownames(P1), top)
  }
  P1 <- P1[sel, ]
  P1 <- P1[order(rowMeans(P1, na.rm = TRUE)), ]

  par(mfrow = c(2, 2), mar = c(5, 15, 4, 4))
  if (is.null(runtime)) par(mfrow = c(1, 2))
  barplot(t(P1),
    las = 1, horiz = TRUE,
    xlab = "cumulative importance", main = "variable importance"
  )
  legend("bottomright", legend = colnames(P1), fill = grey.colors(ncol(P1)))

  par(mar = c(8, 4, 4, 15))
  gx.imagemap(scale(P1, center = FALSE))

  if (!is.null(runtime)) {
    par(mar = c(5, 15, 4, 4))
    runtime <- sort(runtime, decreasing = TRUE)
    runtime <- pmax(runtime, 0.1 * min(runtime[runtime > 0], na.rm = TRUE))
    log.runtime <- log10(runtime)
    log.runtime <- log.runtime - min(log10(runtime))
    barplot(log.runtime,
      las = 1, horiz = TRUE,
      xlab = "time (log10)", main = "algorithm runtime"
    )
  }
}


makePartitionTree <- function(X, y, add.splits = 0) {
  ## formula wants clean names, so save original names
  tx <- t(X)
  colnames(tx) <- gsub("[: +-.,]", "_", colnames(tx))
  colnames(tx) <- gsub("[')(]", "", colnames(tx))
  colnames(tx) <- gsub("\\[|\\]", "", colnames(tx))

  ## create partition tree
  do.survival <- FALSE
  if (do.survival) {
    time <- abs(y)
    status <- (y > 0) ## dead if positive time
    df <- data.frame(time = time + 0.001, status = status, tx)
    rf <- rpart::rpart(survival::Surv(time, status) ~ ., data = df)
  } else {
    df <- data.frame(y = y, tx)
    ## rf <- rpart::rpart(y ~ ., data = df)
    rf <- rpart::rpart(y ~ .,
      data = df, method = "class",
      control = rpart::rpart.control(minsplit = 2, maxdepth = 10)
    )
  }
  rf$orig.names <- rownames(X)
  names(rf$orig.names) <- colnames(tx)

  ## prune the tree
  rf.nsplit <- rf$cptable[, "nsplit"]
  MAXSPLIT <- length(unique(y)) + 1 + add.splits ## maximum N+1 groups
  if (max(rf.nsplit) > MAXSPLIT) {
    cp.idx <- max(which(rf.nsplit <= MAXSPLIT))
    cp0 <- rf$cptable[cp.idx, "CP"]
    rf <- rpart::prune(rf, cp = cp0)
  }
  rf
}


#' @export
plotDecisionTreeFromImportance <- function(imp, add.splits = 0, rf = NULL,
                                           type = c("fancy", "simple", "extended")) {
  if (is.null(imp) && is.null(rf)) {
    message("ERROR: must provide imp or rf")
    return(NULL)
  }

  if (is.null(rf) && !is.null(imp)) {
    kk <- which(!duplicated(names(imp$y)))
    y <- imp$y[kk]
    sel <- intersect(rownames(imp$X), rownames(imp$R))
    X <- imp$X[sel, kk, drop = FALSE]
    rf <- makePartitionTree(X, y, add.splits = add.splits)
  } ## end-if-null rf

  ## plot the tree
  is.surv <- grepl("Surv", rf$call)[2]
  is.surv
  if (is.surv) {
    pkrf <- partykit::as.party(rf)
    partykit::plot.party(pkrf)
  } else {
    if (type == "fancy") {
      ## rattle::fancyRpartPlot(rf, caption = NULL, type=4)
      rpart.plot::rpart.plot(rf)
    } else if (type %in% c("simple", "extended")) {
      pk <- partykit::as.party(rf)
      plot(pk, type = type)
    } else {
      message("[plotDecisionTreeFromImportance] ERROR. unknown type: ", type)
      return(NULL)
    }
  }
}



## ===============================================================================
## ===================== END OF FILE =============================================
## ===============================================================================
