test_that("read_counts works", {
  counts <- read_counts(playbase::example_file("counts.csv"),
    convert_names = FALSE
  )
  expect_equal(nrow(counts), 7439)
})

test_that("read_samples works", {
  samples <- read_samples(playbase::example_file("samples.csv"))
  expect_equal(nrow(samples), 18)
})

test_that("read_contrasts works", {
  contrasts <- read_contrasts(playbase::example_file("contrasts.csv"))
  expect_equal(nrow(contrasts), 6)
})
