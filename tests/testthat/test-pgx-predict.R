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
