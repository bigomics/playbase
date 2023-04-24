if (!identical(Sys.getenv("NOT_CRAN"), "true")) return()

# this takes too long
test_that("createPGX and computePGX works on default example", {
  counts <- read_counts(example_file("counts.csv"))
  samples <- read_samples(example_file("samples.csv"))
  contrasts <- read_contrasts(example_file("contrasts.csv"))

  # create pgx from default counts, samples, contrasts
  testthat::expect_no_error(
    my_pgx <- pgx.createPGX(counts, samples, contrasts)
  )

  # compute pgx from the created pgx with default methods
  testthat::expect_no_error(
    my_pgx <- pgx.computePGX(my_pgx)
  )

})

