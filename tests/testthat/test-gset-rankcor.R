#' Test for gset.rankcor
set.seed(1746)
ranks <- sample(1:10000, 1000, replace = TRUE)
names(ranks) <- replicate(1000, paste(sample(LETTERS, 4, replace = TRUE), collapse = ""))
genesets <- Matrix(rnorm(1000 * 20), ncol = 20)
rownames(genesets) <- names(ranks)

result <- playbase::gset.rankcor(ranks, genesets)

test_that("gset.rankcor returns a list", {
  expect_type(result, "list")
})

test_that("gset.rankcor matrix has correct dimensions", {
  expect_equal(dim(result$rho), c(20, 1))
  expect_equal(dim(result$rho), c(20, 1))
})

test_that("gset.rankcor constant rho", {
  # Mean correlation between random val should be 0
  expect_equal(mean(result$rho), 0, tolerance = 1)
})


test_that("gset.averageFC averages only the members that have a value", {
  F <- cbind(c1 = c(a = 1, b = NA, c = 3, d = NA))
  G <- Matrix::Matrix(cbind(
    s1 = c(1, 1, 1, 0), s2 = c(0, 1, 0, 1), s3 = c(0, 0, 0, 0)
  ), sparse = TRUE)
  rownames(G) <- rownames(F)
  fc <- playbase::gset.averageFC(F, G)
  expect_equal(fc["s1", "c1"], 2, tolerance = 1e-6)
  ## a set whose members are all NA stays NA; an empty set stays 0
  expect_true(is.na(fc["s2", "c1"]))
  expect_equal(fc["s3", "c1"], 0, tolerance = 1e-6)
})
