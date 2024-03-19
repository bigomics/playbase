testthat::test_that("read_counts works", {
  counts <- playbase::read_counts(
    playbase::example_file("counts.csv"),
    convert_names = FALSE
  )
  testthat::expect_equal(nrow(counts), 7439)
})

testthat::test_that("read_samples works", {
  samples <- playbase::read_samples(playbase::example_file("samples.csv"))
  testthat::expect_equal(nrow(samples), 18)
})

testthat::test_that("read_contrasts works", {
  contrasts <- playbase::read_contrasts(playbase::example_file("contrasts.csv"))
  testthat::expect_equal(nrow(contrasts), 6)
})


test_that("read.as_matrix works correctly", {
  # Test that read.as_matrix correctly reads a simple matrix
  #write.table(matrix(1:4, nrow = 2), file = ".//tests/data/test1.csv", sep = ";", row.names = c("a", "b"), col.names = c("c","d") )
  expect_equal(as.numeric(read.as_matrix("./tests/data/test1.csv")), c(1,2,3,4))

  # Test that read.as_matrix correctly reads a matrix with extreme values
  #write.table(matrix(c(1e100, 1e-100, -1e100, -1e-100), nrow = 2), file = ".//tests/data/test2.csv", sep = ";", row.names = c("a", "b"), col.names = c("c","d") )
  expect_equal(as.numeric(read.as_matrix("./tests/data/test2.csv")), c(1e100, 1e-100, -1e100, -1e-100))

})