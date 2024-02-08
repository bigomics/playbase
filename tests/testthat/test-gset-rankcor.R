#' Test for gset.rankcor
set.seed(1746)
ranks <- sample(1:10000, 1000, replace = TRUE)
names(ranks) <- replicate(1000, paste(sample(LETTERS, 4, replace = TRUE), collapse = ""))
genesets <- matrix(rnorm(1000 * 20), ncol = 20)
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

# ' Test for gx.nnmcorrect.SAVE
# ' deprecated?
# x <- matrix(rnorm(100 * 30), 100, 30)
# y <- gl(3, 10)
# test_that("gx.nnmcorrect.SAVE returns a matrix of proper dim", {
#   xcorr <- playbase::gx.nnmcorrect.SAVE(x, y)

#   expect_type(xcorr, "double")
#   expect_equal(dim(xcorr), dim(x))
# })
