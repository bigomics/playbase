## The biomarker method registry, BIOMARKER_METHODS (pgx-biomarkers.R): every
## method scores every feature and predicts proper class probabilities, on two
## and on three classes, at the 6-samples-per-class size the board meets.

## features x samples, 6 samples per class, the first 3 features shifted by class
biomarker_fixture <- function(n_class, n_per_class = 6, n_features = 40, seed = 1) {
  withr::local_seed(seed)
  y <- factor(rep(paste0("class", seq_len(n_class)), each = n_per_class))
  X <- matrix(
    stats::rnorm(n_features * length(y)),
    nrow = n_features,
    dimnames = list(paste0("gene", seq_len(n_features)), paste0("s", seq_along(y)))
  )
  X[1:3, ] <- X[1:3, ] + 2 * as.integer(y)[col(X[1:3, ])]
  ## two samples per class held out for prediction
  test <- unlist(lapply(split(seq_along(y), y), utils::head, 2))
  list(X = X, y = y, test = test)
}

fixtures <- list("2-class" = biomarker_fixture(2), "3-class" = biomarker_fixture(3))
methods <- playbase:::BIOMARKER_METHODS

test_that("the registry holds the six methods with their kinds", {
  expect_named(
    methods,
    c("xgboost", "elastic_net", "random_forest", "splsda", "correlation", "ftest")
  )
  kinds <- vapply(methods, `[[`, "", "kind")
  expect_equal(unname(kinds), rep(c("model", "filter"), c(4, 2)))
  for (m in methods) {
    expect_type(m$label, "character")
    expect_true(is.function(m$importance) && is.function(m$fit_predict))
  }
})

for (name in names(methods)) {
  for (fx in names(fixtures)) {
    test_that(paste(name, "scores and predicts on the", fx, "fixture"), {
      method <- methods[[name]]
      X <- fixtures[[fx]]$X
      y <- fixtures[[fx]]$y
      test <- fixtures[[fx]]$test

      withr::local_seed(42)
      score <- method$importance(X, y)
      expect_length(score, nrow(X))
      expect_named(score, rownames(X))
      expect_true(all(is.finite(score)))
      expect_true(all(score >= 0))
      expect_true(names(which.max(score)) %in% c("gene1", "gene2", "gene3"))

      prob <- method$fit_predict(X[, -test], y[-test], X[, test])
      expect_true(is.matrix(prob))
      expect_identical(dimnames(prob), list(colnames(X)[test], levels(y)))
      expect_true(all(prob >= 0 & prob <= 1))
      expect_equal(unname(rowSums(prob)), rep(1, length(test)), tolerance = 1e-8)
    })
  }
}

test_that("every method is deterministic for a given seed", {
  X <- fixtures[["3-class"]]$X
  y <- fixtures[["3-class"]]$y
  test <- fixtures[["3-class"]]$test
  run <- function(method) {
    withr::local_seed(7)
    list(
      method$importance(X, y),
      method$fit_predict(X[, -test], y[-test], X[, test])
    )
  }
  for (name in names(methods)) {
    expect_identical(run(methods[[name]]), run(methods[[name]]), label = name)
  }
})

test_that("glmnet methods fit when a training class has only 2 samples", {
  fx <- fixtures[["2-class"]]
  train <- c(1:2, 7:8)
  for (name in c("elastic_net", "correlation")) {
    withr::local_seed(3)
    prob <- methods[[name]]$fit_predict(fx$X[, train], fx$y[train], fx$X[, -train])
    expect_equal(unname(rowSums(prob)), rep(1, ncol(fx$X) - 4), tolerance = 1e-8)
  }
})


## ----------------------------------------------------------------------
## Selection, stability and validation
## ----------------------------------------------------------------------

## features x samples of pure noise, with optionally the first 3 features
## shifted by `signal` SDs per class step
signal_fixture <- function(n_per_class, n_features, signal = 0, n_class = 2, seed = 1) {
  withr::local_seed(seed)
  y <- factor(rep(paste0("class", seq_len(n_class)), each = n_per_class))
  X <- matrix(
    stats::rnorm(n_features * length(y)),
    nrow = n_features,
    dimnames = list(paste0("gene", seq_len(n_features)), paste0("s", seq_along(y)))
  )
  X[1:3, ] <- X[1:3, ] + signal * (as.integer(y) - 1)[col(X[1:3, ])]
  list(X = X, y = y)
}

test_that("biomarker_select ranks by percentile and rank product", {
  fx <- signal_fixture(10, 40, signal = 3)
  withr::local_seed(1)
  one <- playbase:::biomarker_select(fx$X, fx$y, "correlation")
  expect_identical(dimnames(one$S), list(rownames(fx$X), "correlation"))
  expect_true(all(one$S > 0 & one$S <= 1))
  ## one method: the consensus is that method
  expect_equal(one$consensus, one$S[, "correlation"])
  expect_identical(one$rank[, "consensus"], one$rank[, "correlation"])

  two <- playbase:::biomarker_select(fx$X, fx$y, c("correlation", "ftest"))
  expect_identical(colnames(two$rank), c("correlation", "ftest", "consensus"))
  expect_setequal(two$rank[, "consensus"], seq_len(nrow(fx$X)))
  expect_setequal(names(sort(two$rank[, "consensus"]))[1:3], c("gene1", "gene2", "gene3"))
})

test_that("biomarker_stability is reproducible and finds a planted signal", {
  fx <- signal_fixture(10, 60, signal = 2.5)
  run <- function() {
    playbase:::biomarker_stability(fx$X, fx$y, "correlation", B = 10, top_n = 5, seed = 42)
  }
  res <- run()
  expect_identical(res, run())
  expect_false(res$bootstrap)
  expect_identical(dimnames(res$frequency), list(rownames(fx$X), c("correlation", "consensus")))
  expect_true(all(res$frequency[c("gene1", "gene2", "gene3"), "consensus"] >= 0.9))
  expect_true(all(res$rank_iqr_lo <= res$rank_median & res$rank_median <= res$rank_iqr_hi))
  expect_true(res$index[["consensus"]] > 0 && res$index[["consensus"]] <= 1)

  ## rounds carry their own seeds, so the core count does not matter
  skip_on_os("windows")
  withr::local_options(mc.cores = 2L)
  expect_true(identical(run(), res)) ## base identical(), not waldo
})

test_that("biomarker_stability falls back to a bootstrap for a tiny class", {
  fx <- signal_fixture(2, 30, signal = 3)
  res <- playbase:::biomarker_stability(fx$X, fx$y, "correlation", B = 3, top_n = 5, seed = 1)
  expect_true(res$bootstrap)
  expect_true(all(is.finite(res$frequency)))
})

## Held-out AUC on pure noise, averaged over datasets: a single noise dataset
## sits anywhere in about 0.5 +/- 0.25 at this size, the mean of six does not.
noise_auc <- function(select_in_fold = TRUE) {
  rowMeans(vapply(1:6, function(seed) {
    fx <- signal_fixture(30, 200, seed = seed)
    cv <- playbase:::biomarker_cv(fx$X, fx$y, "correlation",
      k_grid = 10, repeats = 3, n_perm = 2, seed = seed,
      select_in_fold = select_in_fold
    )
    c(cv = cv$auc[cv$auc$method == "consensus", "auc"], null = cv$null[cv$null$method == "consensus", "auc"])
  }, numeric(2)))
}

test_that("biomarker_cv on pure noise stays at chance", {
  auc <- noise_auc()
  expect_lt(abs(auc[["cv"]] - 0.5), 0.1)
  expect_lt(abs(auc[["null"]] - 0.5), 0.1)
})

test_that("biomarker_cv selection inside the folds is what keeps noise at chance", {
  ## the scientific acceptance gate: selecting on all samples before the CV
  ## leaks the held-out labels into the panel and inflates the AUC
  honest <- noise_auc(select_in_fold = TRUE)[["cv"]]
  leaky <- noise_auc(select_in_fold = FALSE)[["cv"]]
  expect_lt(abs(honest - 0.5), 0.1)
  expect_gt(leaky, 0.75)
})

test_that("biomarker_cv validates a planted signal on held-out samples", {
  fx <- signal_fixture(12, 60, signal = 2)
  cv <- playbase:::biomarker_cv(fx$X, fx$y, c("correlation", "elastic_net"),
    repeats = 2, n_perm = 2, seed = 42
  )
  expect_identical(cv$folds, 5L)
  expect_identical(cv$positive_class, "class1")
  expect_equal(cv$k_grid, c(5, 10, 20, 50))
  cells <- 3 * 4 ## (2 methods + consensus) x 4 panel sizes
  expect_equal(nrow(cv$auc), cells)
  expect_true(all(cv$auc$lo <= cv$auc$auc & cv$auc$auc <= cv$auc$hi))
  expect_gt(cv$auc[cv$auc$method == "consensus" & cv$auc$k == 5, "auc"], 0.85)
  expect_lt(cv$null[cv$null$method == "consensus", "auc"], 0.7)

  expect_equal(nrow(cv$roc), cells * 101)
  expect_true(all(diff(cv$roc$tpr[1:101]) >= 0))
  expect_equal(nrow(cv$confusion), cells * 4)
  ## mean counts per repeat: every row adds up to its class size
  sums <- tapply(cv$confusion$count, cv$confusion[c("method", "k", "truth")], sum)
  expect_true(all(sums == 12))
  expect_named(cv$metrics, c(
    "method", "k", "sensitivity", "specificity", "balanced_accuracy", "ppv", "npv"
  ))
})

test_that("biomarker_cv on three classes reports macro AUC and argmax calls", {
  fx <- signal_fixture(6, 40, signal = 3, n_class = 3)
  cv <- playbase:::biomarker_cv(fx$X, fx$y, "correlation",
    k_grid = c(5, 10), repeats = 2, n_perm = 2, seed = 1
  )
  expect_identical(cv$positive_class, NA_character_)
  expect_null(cv$roc)
  expect_gt(cv$auc[cv$auc$method == "consensus" & cv$auc$k == 5, "auc"], 0.8)
  expect_equal(nrow(cv$confusion), 2 * 2 * 9)
  expect_setequal(cv$confusion$predicted, levels(fx$y))
})

test_that("single-feature AUC names the class a feature is highest in", {
  fx <- signal_fixture(10, 20, signal = 3)
  res <- playbase:::.biomarker_single_auc(fx$X, fx$y)
  expect_identical(res$feature, rownames(fx$X))
  expect_true(all(res$direction[1:3] == "class2"))
  expect_true(all(res$auc[1:3] > 0.9))
})


## ----------------------------------------------------------------------
## pgx.compute_biomarkers
## ----------------------------------------------------------------------

test_that("pgx.compute_biomarkers is reproducible and keeps the legacy shape", {
  pgx <- outlier_pgx()
  run <- function(...) {
    suppressMessages(playbase::pgx.compute_biomarkers(
      pgx, "activated",
      resampling = "fast", top_n = 10, ...
    ))
  }
  calls <- character()
  kind <- RNGkind()
  res <- run(progress = function(percent, message) calls <<- c(calls, sprintf("%.0f %s", percent, message)))
  expect_identical(RNGkind(), kind)
  expect_true(identical(run(), res)) ## base identical(), not waldo

  expect_identical(res$pipeline, "v2")
  expect_identical(res$methods, "xgboost")
  expect_identical(res$params[c("B", "repeats", "n_perm")], list(B = 10L, repeats = 2L, n_perm = 5L))
  expect_length(calls, 10 + 2 + 5)
  expect_identical(utils::tail(calls, 1), "100 Permutation 5/5")

  ## legacy readers: R, X and y over the panel, ordered by consensus
  top <- rownames(res$R)
  expect_identical(colnames(res$R), "xgboost")
  expect_length(top, 10)
  expect_identical(rownames(res$X), top)
  expect_identical(names(res$y), colnames(res$X))
  expect_s3_class(res$rf, "rpart")
  ## one method: its tree is the consensus tree, shared rather than regrown
  expect_identical(names(res$trees), c("consensus", "xgboost"))
  expect_identical(res$trees$consensus, res$rf)
  expect_identical(res$trees$xgboost, res$rf)
  expect_identical(unname(res$rank[top, "consensus"]), 1:10)
  expect_identical(rownames(res$stability$frequency), top)
  expect_identical(res$single_auc$feature, top)
  expect_identical(sum(res$class_counts), length(res$y))
  expect_false(is.null(res$cv))
})

test_that("pgx.compute_biomarkers grows one tree per ranking, on its own panel", {
  pgx <- outlier_pgx()
  run <- function() {
    suppressMessages(playbase::pgx.compute_biomarkers(
      pgx, "activated",
      methods = c("correlation", "ftest"), resampling = "fast", top_n = 10,
      validate = FALSE
    ))
  }
  res <- run()
  expect_true(identical(run(), res)) ## base identical(), not waldo
  expect_true(identical(withr::with_options(list(mc.cores = 2L), run()), res))

  expect_identical(names(res$trees), c("consensus", "correlation", "ftest"))
  expect_identical(res$trees$consensus, res$rf)
  for (ranking in names(res$trees)) {
    tree <- res$trees[[ranking]]
    expect_s3_class(tree, "rpart")
    panel <- rownames(res$rank)[order(res$rank[, ranking])][1:10]
    expect_identical(unname(tree$orig.names), panel, label = ranking)
    expect_identical(environment(tree$terms), baseenv(), label = ranking)
  }
})

test_that("pgx.compute_biomarkers refuses a class under 3 samples", {
  pgx <- outlier_pgx()
  y <- pgx$samples[colnames(pgx$X), "activated"]
  keep <- c(
    utils::head(colnames(pgx$X)[y == "act"], 4),
    utils::head(colnames(pgx$X)[y == "notact"], 2)
  )
  err <- expect_error(
    suppressMessages(playbase::pgx.compute_biomarkers(pgx, "activated", select_samples = keep)),
    class = "biomarker_small_class"
  )
  expect_identical(err$class_counts[["notact"]], 2L)
})

test_that("pgx.compute_biomarkers scores only the selected samples", {
  pgx <- outlier_pgx()
  y <- pgx$samples[colnames(pgx$X), "activated"]
  keep <- c(
    utils::head(colnames(pgx$X)[y == "act"], 4),
    utils::head(colnames(pgx$X)[y == "notact"], 4)
  )
  res <- suppressMessages(playbase::pgx.compute_biomarkers(
    pgx, "activated", select_samples = keep, resampling = "fast", top_n = 10L
  ))
  expect_setequal(colnames(res$X), keep)
  expect_identical(colnames(res$X), names(res$y))
})
