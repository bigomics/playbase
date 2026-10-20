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

test_that("xgboost ranks many features on a small cohort", {
  ## 9 per class and one feature splitting the classes: every tree used to
  ## split on it alone, leaving the rest of the ranking to ties at 0
  for (seed in 1:3) {
    fx <- signal_fixture(9, 200, signal = 4, seed = seed)
    withr::local_seed(42)
    score <- methods$xgboost$importance(fx$X, fx$y)
    label <- paste("seed", seed)
    expect_gte(sum(score > 0), 20, label = label)
    ## top 5, not top 3: on Linux seed 3 puts a noise feature above gene3
    expect_true(all(c("gene1", "gene2", "gene3") %in% names(sort(-score))[1:5]), label = label)
  }
})

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

test_that("biomarker_stability does not count 0-score ties as selected", {
  ## pure noise, rows in SD order as the driver screens them: sparse methods
  ## score most features exactly 0, and those ties rank by row. Were they
  ## counted, the highest-SD rows would fill every round's top 50.
  for (seed in 1:2) {
    fx <- signal_fixture(15, 500, seed = seed)
    X <- fx$X[order(-matrixStats::rowSds(fx$X)), ]
    stability <- function(method) {
      playbase:::biomarker_stability(X, fx$y, method, B = 10, top_n = 50, seed = seed)
    }
    baseline <- stability("correlation")$index[["correlation"]]
    for (method in c("xgboost", "elastic_net", "splsda")) {
      res <- stability(method)
      label <- paste(method, "seed", seed)
      expect_lte(res$index[[method]], baseline + 0.15, label = label)
      ## not the run of top-SD rows that tie-breaking would select every round
      expect_lte(sum(res$frequency[1:50, method] == 1), 3, label = label)
      expect_identical(res$frequency[, method], res$frequency[, "consensus"], label = label)
    }
  }
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

test_that("biomarker_cv confusion counts are integers adding up to each class", {
  ## largest remainders: the row (0.2, 1.6, 1.6, 1.6) of a class of 5 rounds
  ## to (0, 2, 2, 1), where rounding then fixing the diagonal gave -1
  row <- matrix(c(0.2, 1.6, 1.6, 1.6), 1)
  expect_identical(playbase:::.biomarker_round_rows(row), matrix(c(0L, 2L, 2L, 1L), 1))

  withr::local_seed(3)
  y <- factor(rep(c("a", "b", "c"), c(7, 9, 5)))
  X <- matrix(stats::rnorm(30 * length(y)), 30, dimnames = list(paste0("gene", 1:30), paste0("s", seq_along(y))))
  cv <- playbase:::biomarker_cv(X, y, "correlation", k_grid = 5, repeats = 20, n_perm = 0, seed = 1)
  expect_type(cv$confusion$count, "integer")
  expect_true(all(cv$confusion$count >= 0))
  sums <- with(cv$confusion, tapply(count, list(method, truth), sum))
  expect_identical(sums["consensus", ], c(a = 7L, b = 9L, c = 5L))
  expect_identical(sums["correlation", ], c(a = 7L, b = 9L, c = 5L))
})

test_that("the consensus is validated on the consensus top k, refitted by every method", {
  fx <- signal_fixture(10, 30, signal = 2)
  ## a fixed ranking with three different top-5 panels: correlation genes
  ## 1-5, ftest 6-10, consensus 11-15
  ranking <- cbind(correlation = 1:30, ftest = c(6:10, 1:5, 11:30), consensus = c(6:15, 1:5, 16:30))
  rownames(ranking) <- rownames(fx$X)
  selection <- list(rank = ranking, scored = ranking > 0)
  panels <- character()
  spy <- lapply(playbase:::BIOMARKER_METHODS, function(m) {
    fit_predict <- m$fit_predict
    m$fit_predict <- function(Xtr, ytr, Xte) {
      panels <<- c(panels, paste(rownames(Xtr), collapse = ","))
      fit_predict(Xtr, ytr, Xte)
    }
    m
  })
  local_mocked_bindings(BIOMARKER_METHODS = spy, .package = "playbase")
  genes <- function(i) paste0("gene", i, collapse = ",")

  withr::local_seed(1)
  playbase:::.biomarker_cv_pass(fx$X, fx$y, c("correlation", "ftest"), 5L, 5L, selection)
  ## per fold: each method on its own panel, then both on the consensus panel
  expect_identical(panels, rep(c(genes(1:5), genes(6:10), genes(11:15), genes(11:15)), 5))

  ## one method is its own consensus: not refitted
  panels <- character()
  one <- lapply(selection, function(m) m[, c("correlation", "consensus")])
  pass <- playbase:::.biomarker_cv_pass(fx$X, fx$y, "correlation", 5L, 5L, one)
  expect_identical(panels, rep(genes(1:5), 5))
  expect_identical(pass$prob[, , "consensus", ], pass$prob[, , "correlation", ])
})

test_that("a validation panel holds only the features its ranking scored", {
  fx <- signal_fixture(10, 30, signal = 2)
  rank <- cbind(elastic_net = 1:30, consensus = 1:30)
  rownames(rank) <- rownames(fx$X)
  panels <- character()
  spy <- lapply(playbase:::BIOMARKER_METHODS, function(m) {
    fit_predict <- m$fit_predict
    m$fit_predict <- function(Xtr, ytr, Xte) {
      panels <<- c(panels, paste(rownames(Xtr), collapse = ","))
      fit_predict(Xtr, ytr, Xte)
    }
    m
  })
  local_mocked_bindings(BIOMARKER_METHODS = spy, .package = "playbase")

  ## one scored feature: k = 5 and k = 20 both fit on it alone, unpadded by
  ## the unscored features ranked after it
  withr::local_seed(1)
  pass <- playbase:::.biomarker_cv_pass(fx$X, fx$y, "elastic_net", c(5L, 20L), 5L, list(rank = rank, scored = rank <= 1))
  expect_identical(panels, rep("gene1", 10))
  expect_true(all(is.finite(pass$prob)))

  ## nothing scored: no fit, the training fold's class frequencies
  panels <- character()
  pass <- playbase:::.biomarker_cv_pass(fx$X, fx$y, "elastic_net", 5L, 5L, list(rank = rank, scored = rank < 0))
  expect_length(panels, 0)
  expect_equal(unique(as.vector(pass$prob)), 0.5)
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
  expect_identical(res$cv$null_k, 10L)
  expect_identical(unique(res$cv$null$k), 10L)

  ## rounds carry their own seeds, so the core count does not matter
  skip_on_os("windows")
  expect_true(identical(withr::with_options(list(mc.cores = 2L), run()), res))
})

test_that("pgx.compute_biomarkers runs at the largest seed", {
  ## per-round seeds wrap into R's integer range instead of overflowing
  res <- suppressMessages(playbase::pgx.compute_biomarkers(
    outlier_pgx(), "activated",
    resampling = "fast", top_n = 10, seed = .Machine$integer.max
  ))
  expect_identical(res$params$seed, .Machine$integer.max)
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

test_that("pgx.compute_biomarkers reports only the features a ranking scored", {
  ## xgboost and sPLS-DA leave most of the 500 screened features at 0, and
  ## those ties rank by SD; sPLS-DA keeps 25 features on its one component.
  ## mixOmics warns of a feature constant within a stability subsample.
  res <- suppressWarnings(suppressMessages(playbase::pgx.compute_biomarkers(
    outlier_pgx(), "activated",
    methods = c("xgboost", "splsda"), resampling = "fast", top_n = 50L, validate = FALSE
  )))
  in_panel <- res$rank <= 50 & res$scored
  expect_identical(dimnames(res$scored), dimnames(res$rank))
  expect_true(all(rowSums(in_panel) > 0)) ## no row is filler
  expect_identical(sum(in_panel[, "splsda"]), 25L)
  expect_lt(sum(in_panel[, "xgboost"]), 50L)
  expect_identical(rownames(res$stability$frequency), rownames(res$X))
  ## the consensus panel leads, and each tree is grown on its own panel
  expect_true(all(in_panel[seq_len(sum(in_panel[, "consensus"])), "consensus"]))
  for (ranking in names(res$trees)) {
    panel <- rownames(res$rank)[in_panel[, ranking]]
    panel <- panel[order(res$rank[panel, ranking])]
    expect_identical(unname(res$trees[[ranking]]$orig.names), panel, label = ranking)
  }
})

test_that("pgx.compute_biomarkers refuses a run where no method scores", {
  zero <- playbase:::BIOMARKER_METHODS
  zero$correlation$importance <- function(X, y) stats::setNames(numeric(nrow(X)), rownames(X))
  local_mocked_bindings(BIOMARKER_METHODS = zero, .package = "playbase")
  expect_error(
    suppressMessages(playbase::pgx.compute_biomarkers(outlier_pgx(), "activated", methods = "correlation")),
    class = "biomarker_none_scored"
  )
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

test_that("pgx.compute_biomarkers refuses fewer than 5 features", {
  pgx <- outlier_pgx()
  pgx$X <- pgx$X[1:4, ]
  err <- expect_error(
    suppressMessages(playbase::pgx.compute_biomarkers(pgx, "activated")),
    class = "biomarker_too_few_features"
  )
  expect_identical(err$n_features, 4L)
})

test_that("pgx.compute_biomarkers drops samples with an empty phenotype", {
  pgx <- outlier_pgx()
  blank <- colnames(pgx$X)[1:2]
  pgx$samples[blank, "activated"] <- ""
  res <- suppressMessages(playbase::pgx.compute_biomarkers(
    pgx, "activated", resampling = "fast", top_n = 10L, validate = FALSE
  ))
  expect_false(any(blank %in% names(res$y)))
  expect_setequal(names(res$class_counts), c("act", "notact"))
  expect_identical(sum(res$class_counts), length(res$y))
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
