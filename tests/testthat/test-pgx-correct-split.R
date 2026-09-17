# Batch-correction split smoke coverage.
# These tests keep live diagnostic and comparison entry points loadable.
# They avoid optional correction methods and exercise stable public behavior.

test_that("split diagnostics remain callable", {
  pheno <- data.frame(
    group = c("A", "A", "B", "B"),
    batch = c("x", "y", "x", "y"),
    row.names = paste0("s", 1:4)
  )

  result <- checkConfounders(pheno, model.par = "group", max.rho = 0.99)

  expect_named(result, c("confounding", "not.confounding", "rho"))
  expect_null(result$confounding)
  expect_identical(result$not.confounding, "batch")
  expect_equal(unname(result$rho), matrix(0, nrow = 1, ncol = 1))
})

test_that("split comparison runner preserves the uncorrected path", {
  X <- matrix(
    seq_len(12),
    nrow = 3,
    dimnames = list(paste0("g", 1:3), paste0("s", 1:4))
  )

  result <- runBatchCorrectionMethods(
    X,
    batch = NULL,
    y = c("A", "A", "B", "B"),
    ntop = Inf,
    methods = "uncorrected"
  )

  expect_named(result, "uncorrected")
  expect_identical(result$uncorrected, X)
})
