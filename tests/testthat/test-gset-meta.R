library(Matrix)

#' Fixture for gene-set testing: a small expression matrix with three groups,
#' four gene sets and two contrasts. Deterministic, so the reference
#' statistics recorded below stay comparable.
gset_test_fixture <- function() {
  set.seed(42)

  genes <- paste0("G", 1:200)
  samples <- paste0("S", 1:6)
  grp <- factor(c("A", "A", "B", "B", "C", "C"))

  X <- matrix(
    rnorm(200 * 6), nrow = 200, ncol = 6,
    dimnames = list(genes, samples)
  )
  X[1:30, 3:6] <- X[1:30, 3:6] + 2
  X[101:140, 5:6] <- X[101:140, 5:6] - 1.5

  G <- matrix(
    0L, nrow = 200, ncol = 4,
    dimnames = list(genes, c("SET_A", "SET_B", "SET_C", "SET_D"))
  )
  G[1:40, "SET_A"] <- 1L
  G[21:90, "SET_B"] <- sample(c(0L, 1L), 70, replace = TRUE)
  G[101:160, "SET_C"] <- 1L
  G[c(5:15, 150:199), "SET_D"] <- 1L
  G <- Matrix(G, sparse = TRUE)

  # cell-means design, as produced by the playbase model parameters
  design <- model.matrix(~ 0 + grp)
  rownames(design) <- samples
  colnames(design) <- levels(grp)
  contr.matrix <- matrix(
    c(-1, 1, 0, -1, 0, 1), nrow = 3, ncol = 2,
    dimnames = list(levels(grp), c("B_vs_A", "C_vs_A"))
  )

  list(
    X = X, G = G, gmt = playbase::mat2gmt(G),
    design = design, contr.matrix = contr.matrix
  )
}

run_gset_test <- function(fix, methods) {
  playbase::gset.fitContrastsWithAllMethods(
    gmt = fix$gmt,
    X = fix$X,
    G = fix$G,
    design = fix$design,
    contr.matrix = fix$contr.matrix,
    methods = methods,
    use.replaid = FALSE,
    mc.cores = 1
  )
}

#' Issue #563: a dataset without gene-set test methods (gset.methods: [])
#' crashed in gset.fitContrastsWithAllMethods instead of skipping the test.
test_that("no gene-set methods returns empty results instead of failing", {
  fix <- gset_test_fixture()

  expect_no_error(res <- run_gset_test(fix, methods = character(0)))

  expect_named(res, c("meta", "sig.counts", "outputs", "matrices", "timings"))
  expect_length(res$meta, 0)
  expect_length(res$sig.counts, 0)
  expect_length(res$outputs, 0)
  expect_length(res$matrices, 0)
  expect_equal(dim(res$timings), c(0L, 5L))
})

#' Methods that are all unsupported leave nothing to test either.
test_that("unsupported method names are handled like no methods", {
  fix <- gset_test_fixture()

  expect_no_error(res <- run_gset_test(fix, methods = "not_a_method"))

  expect_length(res$meta, 0)
  expect_length(res$sig.counts, 0)
  expect_length(res$outputs, 0)
})

#' The compute step must not test or store anything when no methods are given.
test_that("compute_testGenesets skips the genesets test without methods", {
  fix <- gset_test_fixture()
  pgx <- list(
    X = fix$X,
    GMT = fix$G,
    timings = c(),
    model.parameters = list(
      design = fix$design,
      contr.matrix = fix$contr.matrix
    )
  )

  expect_no_error(res <- playbase::compute_testGenesets(
    pgx,
    test.methods = character(0)
  ))

  expect_null(res$gset.meta)
  expect_null(res$gsetX)
  expect_identical(dim(res$GMT), dim(fix$G))
  expect_identical(res$X, fix$X)
})

#' The empty-methods guard must not touch the case where methods are given:
#' these numbers were recorded from playbase 4.1.3 before the fix.
test_that("gene-set statistics are unchanged when methods are given", {
  fix <- gset_test_fixture()
  res <- run_gset_test(fix, methods = c("camera", "fry"))

  expect_named(res$meta, c("B_vs_A", "C_vs_A"))
  expect_named(res$sig.counts, c("camera", "fry"))
  expect_named(res$matrices, "meta")
  expect_equal(dim(res$timings), c(2L, 5L))

  psums <- unname(vapply(res$meta, function(m) sum(unlist(m$p)), numeric(1)))
  qsums <- unname(vapply(res$meta, function(m) sum(unlist(m$q)), numeric(1)))
  expect_equal(psums, c(2.96604634073192, 1.86057815753857), tolerance = 0.1)
  expect_equal(qsums, c(3.62245660837170, 1.92073672563712), tolerance = 0.1)

  both <- unname(vapply(res$sig.counts, function(s) sum(s$both), numeric(1)))
  up <- unname(vapply(res$sig.counts, function(s) sum(s$up), numeric(1)))
  down <- unname(vapply(res$sig.counts, function(s) sum(s$down), numeric(1)))
  expect_equal(both, c(35, 28), tolerance = 0.1)
  expect_equal(up, c(23, 22), tolerance = 0.1)
  expect_equal(down, c(12, 6), tolerance = 0.1)

  # the gene-set meta matrix is centred, so it sums to zero
  expect_equal(dim(res$matrices$meta), c(4L, 6L))
  expect_equal(sum(res$matrices$meta), 0, tolerance = 1e-10)

  res2 <- run_gset_test(fix, methods = "spearman")
  expect_named(res2$sig.counts, "spearman")
  expect_named(res2$matrices, c("rnkcorr", "meta"))
  expect_equal(dim(res2$timings), c(1L, 5L))
  expect_equal(
    unname(vapply(res2$meta, function(m) sum(unlist(m$p)), numeric(1))),
    c(1.86161384412956, 0.730440106312101), tolerance = 0.1
  )
  expect_equal(sum(res2$matrices$rnkcorr), 0.0872114636966989, tolerance = 0.1)
})
