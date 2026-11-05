testthat::test_that("read_counts works", {
  counts <- playbase::read_counts(
    playbase::example_file("counts.csv")
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
  # write.table(matrix(1:4, nrow = 2), file = ".//tests/data/test1.csv", sep = ";", row.names = c("a", "b"), col.names = c("c","d") )
  temp_file1 <- here::here("tests", "data", "test1.csv")
  expect_equal(as.numeric(read.as_matrix(temp_file1)), c(1, 2, 3, 4))

  # Test that read.as_matrix correctly reads a matrix with extreme values
  # write.table(matrix(c(1e100, 1e-100, -1e100, -1e-100), nrow = 2), file = ".//tests/data/test2.csv", sep = ";", row.names = c("a", "b"), col.names = c("c","d") )
  temp_file2 <- here::here("tests", "data", "test2.csv")
  expect_equal(as.numeric(read.as_matrix(temp_file2)), c(1e100, 1e-100, -1e100, -1e-100))

  temp_file3 <- here::here("tests", "data", "large_integers.csv")
  expect_equal(
    as.numeric(read.as_matrix(temp_file3)),
    c(395000000, 895050000, 84760000000, 4760700000, 2390000000, 1290000000, 4680000000, 4680000000)
  )
})

test_that("read.as_matrix reads a CP1252 (Windows Excel) csv", {
  ## 0xB5 is "µ" in CP1252 and an invalid byte in UTF-8
  f <- tempfile(fileext = ".csv")
  writeBin(charToRaw("Gene;Gd 10 \xb5gL 1;Ctrl 1\r\ntpma;1.5;2\r\n\xb5x;3;4\r\n"), f)
  x <- read.as_matrix(f, as.char = FALSE)
  expect_equal(colnames(x), c("Gd 10 µgL 1", "Ctrl 1"))
  expect_equal(rownames(x), c("tpma", "µx"))
  expect_equal(as.numeric(x), c(1.5, 3, 2, 4))
})
