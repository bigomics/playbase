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
