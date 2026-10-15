## Regression cover for the Biomarker board's compute against the counts/X
## sample split. `pgx.compute_importance()` named its phenotype vector after
## `pgx$samples` and then subset `pgx$X` by those names, which is out of bounds
## once outlier removal has taken a sample out of X (D-24).
## The fixture is helper-outlier-pgx.R's `outlier_pgx()`; it really drops two.

test_that("pgx.compute_importance runs on X's samples when an outlier was dropped", {
  pgx <- outlier_pgx()
  dropped <- setdiff(rownames(pgx$samples), colnames(pgx$X))
  expect_gt(length(dropped), 0)

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  res <- suppressMessages(utils::capture.output(
    imp <- playbase::pgx.compute_importance(pgx, pheno = "activated", nfeatures = 20)
  ))
  expect_false(is.null(imp))

  ## the phenotype is defined on X's samples, so the dropped ones are not scored
  expect_true(all(names(imp$y) %in% colnames(pgx$X)))
  expect_false(any(dropped %in% names(imp$y)))
  expect_equal(length(imp$y), ncol(pgx$X))
  expect_equal(
    as.character(imp$y),
    as.character(pgx$samples[names(imp$y), "activated"])
  )
  expect_gt(nrow(imp$R), 0)
})

## The same object shape on the `multiomics = 1` branch, which nothing exercised.
## That branch rebuilt `y` from `pgx$samples` instead of keeping the vector the
## top of the function had already discretised, cut to pgx$X's samples and cut
## to `select_samples` -- so it crashed on a dropped outlier, returned a factor
## where the function promises a discrete character, and ignored select_samples.
## The fixture is `outlier_mofa_pgx()`; the drop is real and at the shipped
## default outlier_threshold = 3.

test_that("pgx.compute_importance(multiomics = 1) keeps the phenotype it built", {
  pgx <- outlier_mofa_pgx()
  dropped <- setdiff(rownames(pgx$samples), colnames(pgx$X))
  expect_gt(length(dropped), 0)

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  suppressMessages(suppressWarnings(utils::capture.output(
    imp <- playbase::pgx.compute_importance(
      pgx,
      pheno = "activated", nfeatures = 20, multiomics = 1
    )
  )))
  expect_false(is.null(imp))

  ## (1) defined on X's samples, so the dropped one is not scored
  expect_true(all(names(imp$y) %in% colnames(pgx$X)))
  expect_false(any(dropped %in% names(imp$y)))
  expect_equal(length(imp$y), ncol(pgx$X))

  ## (2) discrete, as the function's own contract says, and the right values
  expect_type(imp$y, "character")
  expect_equal(
    as.character(imp$y),
    as.character(pgx$samples[names(imp$y), "activated"])
  )

  ## the matrix it hands back is pgx$X itself, on those samples
  expect_identical(colnames(imp$X), names(imp$y))
  expect_identical(imp$X, pgx$X[rownames(imp$X), names(imp$y), drop = FALSE])
  expect_gt(nrow(imp$R), 0)
})

test_that("pgx.compute_importance(multiomics = 1) honours select_samples", {
  pgx <- outlier_mofa_pgx()
  keep <- head(colnames(pgx$X), 10)

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  suppressMessages(suppressWarnings(utils::capture.output(
    imp <- playbase::pgx.compute_importance(
      pgx,
      pheno = "activated", nfeatures = 20,
      multiomics = 1, select_samples = keep
    )
  )))

  ## (3) the branch used to drop this filter on the floor and score all of them
  expect_setequal(names(imp$y), keep)
  expect_identical(colnames(imp$X), names(imp$y))
})

## Seeded snapshot of the whole result, recorded before the data preparation
## moved into `.biomarker_prepare()` (pgx-biomarkers.R). Every random draw in
## the function -- the augmentation `sample()`, the `rnorm()` noise, the
## cv.glmnet folds, randomForest and rpart's cross-validation -- comes from R's
## RNG and runs serially, so one explicit seed and RNG kind fix the result. The
## tree is reduced to its structure; the rpart object also holds a call and
## terms that say nothing about behaviour. Recorded on macOS: xgboost, glmnet
## and randomForest builds can rank ties differently elsewhere, so CI skips it.

test_that("pgx.compute_importance output is unchanged (seeded snapshot)", {
  skip_on_ci()
  pgx <- outlier_pgx()

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  withr::local_seed(
    2024,
    .rng_kind = "Mersenne-Twister",
    .rng_normal_kind = "Inversion",
    .rng_sample_kind = "Rejection"
  )
  suppressMessages(utils::capture.output(
    imp <- playbase::pgx.compute_importance(pgx, pheno = "activated", nfeatures = 20)
  ))

  expect_named(imp, c("R", "y", "X", "rf"))
  expect_snapshot_value(
    list(
      R = imp$R,
      y = imp$y,
      X = imp$X,
      rf = list(
        frame = imp$rf$frame[, c("var", "n", "wt", "dev", "yval")],
        splits = imp$rf$splits,
        orig.names = imp$rf$orig.names
      )
    ),
    style = "deparse",
    tolerance = 1e-6
  )
})
