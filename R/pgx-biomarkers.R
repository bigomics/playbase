##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##


## ----------------------------------------------------------------------
## Biomarker selection: data preparation and method registry
## ----------------------------------------------------------------------

# Phenotype, sample filter, imputation, feature filter and top-SD cut shared by
# pgx.compute_importance() and the biomarker pipeline. Returns NULL (with a
# message) when too few classes, samples or features remain. drop_empty also
# drops samples whose phenotype is "" (the legacy path keeps them as a class).
.biomarker_prepare <- function(pgx, pheno, level = "genes",
                               filter_features = NULL,
                               select_features = NULL,
                               select_samples = NULL,
                               nfeatures = 50,
                               drop_empty = FALSE) {
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
  if (drop_empty) y <- y[y != ""]

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


## ----------------------------------------------------------------------
## Selection, stability and validation
## ----------------------------------------------------------------------

# Every ticked method on the given samples. S holds each method's score as a
# rank percentile (features x methods, in (0, 1], higher is better); the
# consensus is the rank product (geometric mean of the methods' ranks) as a
# percentile; rank holds integer ranks, 1 = best, for methods and consensus.
# scored (same columns as rank) is TRUE where the raw score is above 0, for
# the consensus in any method: sparse methods (xgboost, elastic net, sPLS-DA)
# score most features exactly 0, and those ties rank by row, i.e. by SD.
biomarker_select <- function(X, y, methods) {
  y <- as.factor(y)
  d <- nrow(X)
  scores <- vapply(methods, function(m) {
    BIOMARKER_METHODS[[m]]$importance(X, y)[rownames(X)]
  }, numeric(d))
  scores <- matrix(scores, d, dimnames = list(rownames(X), methods))
  S <- scores
  S[] <- apply(scores, 2, rank) / d
  product <- exp(rowMeans(log(d + 1 - d * S)))
  consensus <- rank(-product) / d
  rank <- apply(-cbind(S, consensus = consensus), 2, rank, ties.method = "first")
  scored <- !is.na(scores) & scores > 0 ## an F-test can give NA on a constant row
  scored <- cbind(scored, consensus = rowSums(scored) > 0)
  list(S = S, consensus = consensus, rank = rank, scored = scored)
}

# Selection on B stratified subsamples (frac of each class, or a stratified
# bootstrap when that leaves a class under 2): how often each feature is
# selected, its rank median and IQR, and the Nogueira index of the selections.
# A feature is selected in a round when it ranks in the top_n and is scored
# (biomarker_select), so a 0 score never counts however it ranks.
biomarker_stability <- function(X, y, methods, B, frac = 0.8, top_n, seed,
                                progress = NULL) {
  y <- as.factor(y)
  by_class <- split(seq_along(y), y)
  size <- floor(frac * lengths(by_class))
  bootstrap <- any(size < 2)
  if (bootstrap) size <- lengths(by_class)
  rounds <- .biomarker_rounds(B, seed, function(b) {
    idx <- unlist(Map(function(i, n) i[sample.int(length(i), n, replace = bootstrap)], by_class, size))
    selection <- biomarker_select(X[, idx, drop = FALSE], y[idx], methods)
    list(rank = selection$rank, selected = selection$rank <= top_n & selection$scored)
  }, "Stability round", progress)
  ## features x (methods + consensus) x rounds
  ranks <- simplify2array(lapply(rounds, `[[`, "rank"))
  in_top <- simplify2array(lapply(rounds, `[[`, "selected"))
  list(
    B = B,
    bootstrap = bootstrap,
    frequency = apply(in_top, c(1, 2), mean),
    rank_median = apply(ranks, c(1, 2), stats::median),
    rank_iqr_lo = apply(ranks, c(1, 2), stats::quantile, 0.25, names = FALSE),
    rank_iqr_hi = apply(ranks, c(1, 2), stats::quantile, 0.75, names = FALSE),
    index = apply(in_top, 2, .biomarker_nogueira)
  )
}

# Nogueira et al. (2018) stability of selections (features x rounds, logical;
# sizes may vary): 1 when every round picks the same set, about 0 for random
# sets. NA when the rounds select nothing or everything.
.biomarker_nogueira <- function(selected) {
  d <- nrow(selected)
  k <- mean(colSums(selected))
  if (k == 0 || k >= d) {
    return(NA_real_)
  }
  p <- rowMeans(selected)
  variance <- ncol(selected) / (ncol(selected) - 1) * p * (1 - p)
  1 - mean(variance) / (k / d * (1 - k / d))
}

# Runs fun(i) for i in 1..n, each under its own seed, seed + i wrapped into
# R's integer range, in batches of getOption("mc.cores") forks, so results
# depend on neither scheduling nor core count. progress(message) is called
# from this process per finished round.
.biomarker_rounds <- function(n, seed, fun, label, progress = NULL) {
  cores <- getOption("mc.cores", 1L)
  out <- vector("list", n)
  for (batch in split(seq_len(n), ceiling(seq_len(n) / cores))) {
    out[batch] <- parallel::mclapply(batch, function(i) {
      set.seed(as.integer((as.numeric(seed) + i) %% .Machine$integer.max))
      fun(i)
    }, mc.cores = cores)
    failed <- Filter(function(x) inherits(x, "try-error"), out[batch])
    if (length(failed)) stop(attr(failed[[1]], "condition"))
    if (!is.null(progress)) for (i in batch) progress(sprintf("%s %d/%d", label, i, n))
  }
  out
}

# Repeated stratified K-fold CV with selection and fitting inside each training
# fold, at each panel size k; then n_perm label permutations of one repeat at
# null_k. Set select_in_fold = FALSE only to measure selection leakage: it
# selects once on all samples, which inflates the AUC. CV repeats draw their
# seeds from seed + 1e6 and permutations from seed + 2e6, apart from the
# stability rounds' seed + i.
biomarker_cv <- function(X, y, methods, k_grid = c(5, 10, 20, 50), folds = 5L,
                         repeats, n_perm, seed, null_k = 10L,
                         select_in_fold = TRUE, progress = NULL) {
  y <- as.factor(y)
  K <- min(folds, min(table(y)))
  k_grid <- unique(pmin(as.integer(k_grid), nrow(X)))
  null_k <- min(as.integer(null_k), nrow(X))
  ranking <- if (!select_in_fold) biomarker_select(X, y, methods)$rank
  passes <- .biomarker_rounds(repeats, seed + 1e6, function(r) {
    .biomarker_cv_pass(X, y, methods, k_grid, K, ranking)
  }, "Cross-validation repeat", progress)
  null_passes <- .biomarker_rounds(n_perm, seed + 2e6, function(p) {
    yp <- stats::setNames(sample(y), names(y))
    list(y = yp, pass = .biomarker_cv_pass(X, yp, methods, null_k, K, ranking))
  }, "Permutation", progress)

  binary <- nlevels(y) == 2
  cells <- expand.grid(method = c(methods, "consensus"), k = k_grid, stringsAsFactors = FALSE)
  per_cell <- function(fun) {
    out <- do.call(rbind, unname(Map(fun, cells$method, as.character(cells$k))))
    rownames(out) <- NULL
    out
  }
  interval <- function(auc) stats::quantile(auc, c(0.025, 0.975), names = FALSE)

  auc <- per_cell(function(m, k) {
    a <- vapply(passes, function(p) .biomarker_prob_auc(p$prob[, , m, k], y), 0)
    data.frame(method = m, k = as.integer(k), auc = mean(a), lo = interval(a)[1], hi = interval(a)[2])
  })
  roc <- if (binary) {
    per_cell(function(m, k) {
      score <- unlist(lapply(passes, function(p) p$prob[, 1, m, k]))
      data.frame(method = m, k = as.integer(k), .biomarker_roc(score, rep(y == levels(y)[1], repeats)))
    })
  }
  ## held-out calls summed over repeats, then the mean per repeat
  confusion_of <- function(m, k) {
    Reduce(`+`, lapply(passes, function(p) {
      table(truth = y, predicted = factor(levels(y)[p$call[, m, k]], levels(y)))
    })) / repeats
  }
  confusion <- per_cell(function(m, k) {
    counts <- .biomarker_round_rows(confusion_of(m, k))
    data.frame(method = m, k = as.integer(k), as.data.frame(counts, responseName = "count", stringsAsFactors = FALSE))
  })
  metrics <- per_cell(function(m, k) {
    data.frame(method = m, k = as.integer(k), t(.biomarker_metrics(confusion_of(m, k), binary)))
  })
  null <- if (n_perm > 0) {
    do.call(rbind, lapply(c(methods, "consensus"), function(m) {
      a <- vapply(null_passes, function(p) .biomarker_prob_auc(p$pass$prob[, , m, 1], p$y), 0)
      data.frame(method = m, k = null_k, auc = mean(a), lo = interval(a)[1], hi = interval(a)[2])
    }))
  }

  list(
    folds = K, repeats = repeats, n_perm = n_perm, k_grid = k_grid, null_k = null_k,
    positive_class = if (binary) levels(y)[1] else NA_character_,
    auc = auc, roc = roc, confusion = confusion, metrics = metrics, null = null
  )
}

# Rounds each row of a table to integers that keep the row's (integer) sum:
# rows are rounded down, and the units lost go to the largest remainders.
.biomarker_round_rows <- function(x) {
  counts <- floor(x)
  for (i in seq_len(nrow(x))) {
    short <- round(sum(x[i, ]) - sum(counts[i, ]))
    up <- order(counts[i, ] - x[i, ])[seq_len(short)]
    counts[i, up] <- counts[i, up] + 1
  }
  storage.mode(counts) <- "integer"
  counts
}

# One stratified K-fold pass: held-out probabilities (samples x classes x
# (methods + consensus) x k) and class calls (level index) for every sample.
# A method's panel is its top k. The consensus panel is the consensus top k:
# every method refits on it and the consensus is their soft vote (mean
# probability). One method is its own consensus and is not refitted.
.biomarker_cv_pass <- function(X, y, methods, k_grid, K, ranking = NULL) {
  folds <- .biomarker_folds(y, K)
  cols <- c(methods, "consensus")
  dims <- list(names(y), levels(y), cols, as.character(k_grid))
  size <- c(length(y), lengths(dims[-1]))
  prob <- array(NA_real_, size, dims)
  call <- array(NA_integer_, size[-2], dims[-2])
  for (f in seq_len(K)) {
    train <- folds != f
    rank <- ranking
    if (is.null(rank)) rank <- biomarker_select(X[, train, drop = FALSE], y[train], methods)$rank
    ## predict every sample: training rows set the threshold, held-out rows are kept
    fit_predict <- function(m, panel) {
      BIOMARKER_METHODS[[m]]$fit_predict(X[panel, train, drop = FALSE], y[train], X[panel, , drop = FALSE])
    }
    for (k in dims[[4]]) {
      p <- lapply(methods, function(m) fit_predict(m, rank[, m] <= as.integer(k)))
      consensus <- if (length(methods) == 1L) {
        p[[1]]
      } else {
        panel <- rank[, "consensus"] <= as.integer(k)
        Reduce(`+`, lapply(methods, fit_predict, panel)) / length(methods)
      }
      p <- c(p, list(consensus))
      for (j in seq_along(cols)) {
        prob[!train, , j, k] <- p[[j]][!train, ]
        call[!train, j, k] <- .biomarker_call(p[[j]][train, , drop = FALSE], y[train], p[[j]][!train, , drop = FALSE])
      }
    }
  }
  list(prob = prob, call = call)
}

# Class calls for held-out samples. Two classes: the first class when its
# probability reaches the Youden-optimal cut on the training samples' own
# (resubstituted) probabilities. More classes: the most probable class.
.biomarker_call <- function(p_train, y_train, p_test) {
  if (nlevels(y_train) > 2) {
    return(max.col(p_test, ties.method = "first"))
  }
  cut <- .biomarker_youden(p_train[, 1], y_train == levels(y_train)[1])
  ifelse(p_test[, 1] >= cut, 1L, 2L)
}

# Cut maximising sensitivity + specificity - 1, placed midway between the
# score it admits last and the next lower one so it does not sit on a sample.
.biomarker_youden <- function(score, positive) {
  cuts <- sort(unique(score), decreasing = TRUE)
  j <- vapply(cuts, function(t) mean(score[positive] >= t) - mean(score[!positive] >= t), 0)
  i <- which.max(j)
  mean(c(cuts[i], cuts[i + 1L]), na.rm = TRUE)
}

# Mann-Whitney AUC of score for the positive (TRUE) samples.
.biomarker_auc <- function(score, positive) {
  n1 <- sum(positive)
  n0 <- length(positive) - n1
  (sum(rank(score)[positive]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# AUC of class probabilities: of the first class for two classes, else the
# one-vs-rest macro average.
.biomarker_prob_auc <- function(prob, y) {
  classes <- if (nlevels(y) == 2) levels(y)[1] else levels(y)
  mean(vapply(classes, function(cl) .biomarker_auc(prob[, cl], y == cl), 0))
}

# ROC as the true-positive rate on a fixed grid of false-positive rates, so a
# pooled curve stays small and curves align across methods.
.biomarker_roc <- function(score, positive, grid = seq(0, 1, by = 0.01)) {
  o <- order(score, decreasing = TRUE)
  last <- !duplicated(score[o], fromLast = TRUE) ## one point per distinct cut
  tpr <- c(0, cumsum(positive[o])[last] / sum(positive))
  fpr <- c(0, cumsum(!positive[o])[last] / sum(!positive))
  data.frame(fpr = grid, tpr = tpr[findInterval(grid, fpr)])
}

# Sensitivity, specificity, balanced accuracy, PPV and NPV from a truth x
# predicted table: of the first class for two classes, else macro averages of
# the one-vs-rest values. Undefined ratios are NA.
.biomarker_metrics <- function(confusion, binary) {
  classes <- if (binary) rownames(confusion)[1] else rownames(confusion)
  per_class <- vapply(classes, function(cl) {
    tp <- confusion[cl, cl]
    fn <- sum(confusion[cl, ]) - tp
    fp <- sum(confusion[, cl]) - tp
    tn <- sum(confusion) - tp - fn - fp
    c(sensitivity = tp / (tp + fn), specificity = tn / (tn + fp), ppv = tp / (tp + fp), npv = tn / (tn + fn))
  }, numeric(4))
  m <- rowMeans(per_class)
  m <- c(m[1:2], balanced_accuracy = mean(m[1:2]), m[3:4])
  m[is.nan(m)] <- NA
  m
}

# Descriptive single-feature AUC on all samples: each feature's AUC for the
# class it is highest in (vs the rest), named by that class.
.biomarker_single_auc <- function(X, y) {
  y <- as.factor(y)
  means <- vapply(levels(y), function(cl) rowMeans(X[, y == cl, drop = FALSE]), numeric(nrow(X)))
  direction <- levels(y)[max.col(matrix(means, nrow(X)), ties.method = "first")]
  auc <- vapply(seq_len(nrow(X)), function(i) .biomarker_auc(X[i, ], y == direction[i]), 0)
  data.frame(feature = rownames(X), auc = auc, direction = direction)
}


## ----------------------------------------------------------------------
## Driver
## ----------------------------------------------------------------------

#' Seeded biomarker selection with stability and held-out validation
#'
#' @description Selects biomarkers for a discrete phenotype with one or more
#' methods on the original samples, measures how stable the selection is under
#' resampling, and estimates held-out classification performance against a
#' permutation null. The same inputs and seed give an identical result.
#'
#' @param pgx A pgx object.
#' @param pheno Column of \code{pgx$samples} to predict. It is discretised as
#'   in \code{pgx.compute_importance()}; samples whose value is \code{NA} or
#'   \code{""} are dropped.
#' @param methods Selection methods: any of \code{"xgboost"},
#'   \code{"elastic_net"}, \code{"random_forest"}, \code{"splsda"},
#'   \code{"correlation"} and \code{"ftest"}.
#' @param resampling \code{"standard"}, \code{"fast"} or \code{"thorough"}:
#'   (stability rounds, CV repeats, permutations) of (50, 10, 20), (10, 2, 5)
#'   and (200, 20, 100).
#' @param seed Integer seed.
#' @param top_n Number of biomarkers to report. Ten times as many features,
#'   by SD, are screened.
#' @param ... Passed to the data preparation: \code{level},
#'   \code{filter_features}, \code{select_features}, \code{select_samples}.
#' @param validate Run the nested cross-validation and the permutation null.
#' @param progress \code{NULL}, or \code{function(percent, message)} called
#'   after every stability round, CV repeat and permutation.
#'
#' @details Steps, in order:
#' \enumerate{
#'   \item Selection on all samples. Each method's scores become rank
#'     percentiles; the consensus is their rank product (geometric mean of
#'     ranks). A decision tree is grown on the consensus top features and on
#'     each method's.
#'   \item Stability: selection on B stratified 80\% subsamples (a stratified
#'     bootstrap when a class is too small), giving per-feature selection
#'     frequency, rank median and IQR, and the Nogueira stability index. A
#'     feature counts as selected in a round when it ranks in the top
#'     \code{top_n} and its raw score is above 0 (for the consensus, in any
#'     method): sparse methods score most features exactly 0, and their rank
#'     among those ties only reflects SD order. Ranks are reported unchanged.
#'   \item Validation: repeated stratified K-fold CV, K = min(5, smallest class),
#'     with selection and fitting inside each training fold at panel sizes
#'     5, 10, 20 and 50. Each method is fitted on its own top k. The consensus
#'     validates the consensus top k: every method is refitted on it and the
#'     consensus predicts by the mean of their probabilities (a single method
#'     is its own consensus). Two-class calls use a Youden threshold chosen on
#'     the training fold; multiclass calls take the most probable class. Labels
#'     are permuted to give a null AUC at panel size 10 (\code{cv$null_k}),
#'     with the same definitions.
#' }
#' Rounds run under \code{parallel::mclapply} with
#' \code{getOption("mc.cores", 1)} cores, each with its own seed, so results do
#' not depend on the core count. Stability rounds, CV repeats and permutations
#' draw seeds from separate ranges, wrapped into R's integer range, so any
#' seed in [0, 2^31 - 1] works. The caller's RNG kind is restored on exit.
#' Multi-omics MOFA importance is not supported here; see
#' \code{pgx.compute_importance()}.
#'
#' @return \code{NULL} when too few classes, samples or features remain, as
#' \code{pgx.compute_importance()}. Otherwise a list whose feature rows are
#' the union of the consensus and each method's top \code{top_n}, in consensus
#' order (the first \code{top_n} rows are the consensus panel):
#' \item{R}{features x methods rank percentiles, the legacy importance matrix.}
#' \item{y}{named phenotype vector of the samples used.}
#' \item{X}{features x samples expression.}
#' \item{rf}{rpart decision tree on the consensus panel, as \code{trees$consensus}.}
#' \item{trees}{named list of rpart decision trees, \code{consensus} and one per
#'   method, each grown on that ranking's top \code{top_n} features.}
#' \item{pipeline}{\code{"v2"}.}
#' \item{methods}{the methods run.}
#' \item{params}{resampling, B, repeats, n_perm, seed and top_n.}
#' \item{class_counts}{samples per class.}
#' \item{n_features_screened}{features entering selection.}
#' \item{consensus}{consensus rank percentile per feature.}
#' \item{rank}{features x (methods, consensus) integer ranks, 1 = best.}
#' \item{stability}{B, bootstrap, and features x (methods, consensus)
#'   frequency, rank_median, rank_iqr_lo, rank_iqr_hi; index per column.}
#' \item{cv}{\code{NULL} unless \code{validate}: folds, repeats, n_perm,
#'   k_grid, null_k (the null's panel size), positive_class, and data frames
#'   auc, roc (two classes only), confusion (integer mean counts per repeat,
#'   each row adding up to its class size), metrics and null.}
#' \item{single_auc}{data frame of feature, all-sample AUC and the class the
#'   feature is highest in (descriptive, not validated).}
#' \item{playbase_version}{package version.}
#'
#' An error of class \code{biomarker_small_class} (carrying
#' \code{class_counts}) is raised when a class has fewer than 3 samples, and
#' one of class \code{biomarker_too_few_features} (carrying
#' \code{n_features}) when fewer than 5 features remain after filtering.
#'
#' @export
pgx.compute_biomarkers <- function(pgx, pheno, methods = "xgboost",
                                   resampling = c("standard", "fast", "thorough"),
                                   seed = 42L, top_n = 50L, ...,
                                   validate = TRUE, progress = NULL) {
  methods <- match.arg(methods, names(BIOMARKER_METHODS), several.ok = TRUE)
  resampling <- match.arg(resampling)
  prep <- .biomarker_prepare(pgx, pheno, nfeatures = top_n, drop_empty = TRUE, ...)
  if (is.null(prep)) {
    return(NULL)
  }
  if (nrow(prep$X) < 5) {
    stop(errorCondition(
      paste0(
        "At least 5 features are needed to select biomarkers; ",
        nrow(prep$X), " remain after filtering"
      ),
      class = "biomarker_too_few_features", n_features = nrow(prep$X)
    ))
  }
  small <- prep$class_counts < 3
  if (any(small)) {
    stop(errorCondition(
      paste0(
        "A class needs at least 3 samples to validate: ",
        paste0(names(small)[small], " (n=", prep$class_counts[small], ")", collapse = ", ")
      ),
      class = "biomarker_small_class", class_counts = prep$class_counts
    ))
  }

  rng_kind <- RNGkind()
  on.exit(do.call(RNGkind, as.list(rng_kind)), add = TRUE)
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)

  preset <- list(fast = c(10L, 2L, 5L), standard = c(50L, 10L, 20L), thorough = c(200L, 20L, 100L))[[resampling]]
  params <- list(
    resampling = resampling, B = preset[1], repeats = preset[2], n_perm = preset[3],
    seed = seed, top_n = top_n
  )
  ## prepare cuts y to select_samples but keeps X's columns, as
  ## pgx.compute_importance aligns them itself afterwards
  X <- prep$X[, names(prep$y), drop = FALSE]
  y <- factor(prep$y)

  step <- 0
  steps <- params$B + if (validate) params$repeats + params$n_perm else 0
  tick <- if (!is.null(progress)) {
    function(message) {
      step <<- step + 1
      progress(100 * step / steps, message)
    }
  }

  ## the reported ranking and its tree come first, off the global seed, so the
  ## rounds below (each under its own seed) cannot shift them
  selection <- biomarker_select(X, y, methods)
  keep <- rownames(X)[order(selection$rank[, "consensus"])]
  keep <- keep[rowSums(selection$rank[keep, , drop = FALSE] <= top_n) > 0]
  grow_tree <- function(panel) {
    tree <- makePartitionTree(X[panel, , drop = FALSE], prep$y, add.splits = 0)
    ## the formula would otherwise carry the builder's frame, data included,
    ## into the payload and make two equal runs non-identical()
    environment(tree$terms) <- baseenv()
    tree
  }
  ## one tree per ranking, grown serially in a fixed order; a single method's
  ## consensus is that method's ranking, so its tree is shared, not regrown
  trees <- list(consensus = grow_tree(utils::head(keep, top_n)))
  for (method in methods) {
    trees[[method]] <- if (length(methods) == 1L) {
      trees$consensus
    } else {
      grow_tree(utils::head(keep[order(selection$rank[keep, method])], top_n))
    }
  }

  stability <- biomarker_stability(X, y, methods, params$B, top_n = top_n, seed = seed, progress = tick)
  by_feature <- c("frequency", "rank_median", "rank_iqr_lo", "rank_iqr_hi")
  stability[by_feature] <- lapply(stability[by_feature], function(m) m[keep, , drop = FALSE])
  cv <- if (validate) {
    biomarker_cv(X, y, methods,
      repeats = params$repeats, n_perm = params$n_perm, seed = seed, progress = tick
    )
  }

  list(
    R = selection$S[keep, , drop = FALSE],
    y = prep$y,
    X = X[keep, , drop = FALSE],
    rf = trees$consensus,
    trees = trees,
    pipeline = "v2",
    methods = methods,
    params = params,
    class_counts = prep$class_counts,
    n_features_screened = nrow(X),
    consensus = selection$consensus[keep],
    rank = selection$rank[keep, , drop = FALSE],
    stability = stability,
    cv = cv,
    single_auc = .biomarker_single_auc(X[keep, , drop = FALSE], y),
    playbase_version = unname(getNamespaceVersion("playbase"))
  )
}
