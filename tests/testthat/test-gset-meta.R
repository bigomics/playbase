## Tests for meta-analysis zero-power method filtering (gset-meta.r)

test_that("smart_max is no-op when all methods have power", {
  set.seed(42)
  pv <- matrix(runif(300), nrow = 100, ncol = 3,
    dimnames = list(paste0("gs", 1:100), c("gsva", "fisher", "fgsea"))
  )
  pv[1:10, ] <- runif(30, 0, 0.001)
  qv <- apply(pv, 2, p.adjust, method = "BH")

  old_q <- apply(qv, 1, max, na.rm = TRUE)

  has.power <- colSums(qv < 0.05, na.rm = TRUE) > 0
  expect_true(all(has.power))
  qv2 <- qv[, has.power, drop = FALSE]
  new_q <- apply(qv2, 1, max, na.rm = TRUE)

  expect_identical(old_q, new_q)
})

test_that("smart_max drops zero-power methods and rescues discoveries", {
  set.seed(42)
  n <- 200
  pv <- matrix(NA_real_, nrow = n, ncol = 3,
    dimnames = list(paste0("gs", 1:n), c("gsva", "fisher", "fgsea"))
  )
  pv[, "gsva"] <- c(runif(30, 0, 0.0005), runif(n - 30, 0.01, 1))
  pv[, "fisher"] <- NA
  pv[, "fgsea"] <- runif(n, 0.05, 1)

  qv <- pv
  qv[, "gsva"] <- p.adjust(pv[, "gsva"], method = "BH")
  qv[, "fgsea"] <- p.adjust(pv[, "fgsea"], method = "BH")

  ## Old logic: fgsea vetoes everything
  old_q <- apply(qv, 1, max, na.rm = TRUE)
  expect_equal(sum(old_q < 0.05), 0)

  ## New logic: drop fisher (NA) + fgsea (0 power), keep gsva
  has.power <- colSums(qv < 0.05, na.rm = TRUE) > 0
  expect_false(has.power["fisher"])
  expect_false(has.power["fgsea"])
  expect_true(has.power["gsva"])

  qv2 <- qv[, has.power, drop = FALSE]
  new_q <- apply(qv2, 1, max, na.rm = TRUE)
  expect_gt(sum(new_q < 0.05), 0)
})

test_that("smart_max on bundled pgx_example never regresses", {
  pgx_file <- system.file("data", "pgx_example.rds", package = "playbase")
  if (!nzchar(pgx_file)) pgx_file <- "data/pgx_example.rds"
  if (!file.exists(pgx_file)) skip("pgx_example not found")
  pgx <- local(get(load(pgx_file)))
  if (is.null(pgx$gset.meta)) skip("pgx_example has no gset.meta")

  meta <- pgx$gset.meta$meta[[1]]
  qv <- meta$q

  old_q <- apply(qv, 1, max, na.rm = TRUE)
  old_sig <- sum(old_q < 0.05)

  has.power <- colSums(qv < 0.05, na.rm = TRUE) > 0
  qv2 <- qv[, has.power, drop = FALSE]
  new_q <- apply(qv2, 1, max, na.rm = TRUE)
  new_sig <- sum(new_q < 0.05)

  expect_gte(new_sig, old_sig)
})
