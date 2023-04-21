if (!identical(Sys.getenv("NOT_CRAN"), "true")) return()

# this takes too long
test_that("create_pgx and compute example works", {
  counts <- read_counts(example_file("counts.csv"))
  samples <- read_samples(example_file("samples.csv"))
  contrasts <- read_contrasts(example_file("contrasts.csv"))

  testthat::expect_no_error(
    my_pgx2 <- pgx.createPGX(counts, samples, contrasts)
  )

})

