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
