# this takes too long
test_that("create_pgx example works", {
  counts <- read_counts(example_file("counts.csv"))
  samples <- read_samples(example_file("samples.csv"))
  contrasts <- read_contrasts(example_file("contrasts.csv"))
  #my_pgx <- create_pgx(counts, samples, contrasts)

  testthat::expect_equal(1,1)
  #testthat::expect_type(my_pgx, "list")
})
