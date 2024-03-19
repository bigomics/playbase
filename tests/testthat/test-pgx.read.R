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


# Load the testthat package
library(testthat)

# Define the tests
test_that("read.as_matrix works correctly", {
  # Test that read.as_matrix correctly reads a simple matrix
  write.table(matrix(1:4, nrow = 2), file = "test1.txt", sep = "\t", row.names = FALSE, col.names = FALSE)
  expect_equal(read.as_matrix("test1.txt"), matrix(1:4, nrow = 2))

  # Test that read.as_matrix correctly reads a matrix with extreme values
  write.table(matrix(c(1e100, 1e-100, -1e100, -1e-100), nrow = 2), file = "test2.txt", sep = "\t", row.names = FALSE, col.names = FALSE)
  expect_equal(read.as_matrix("test2.txt"), matrix(c(1e100, 1e-100, -1e100, -1e-100), nrow = 2))

  # Add more tests as needed
})

# Run the tests
test_dir("path/to/your/tests")