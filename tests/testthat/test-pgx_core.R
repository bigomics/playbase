if (!identical(Sys.getenv("NOT_CRAN"), "true")) return()

# this takes too long
testthat::test_that("createPGX and computePGX works on default example", {
  # create pgx from default counts, samples, contrasts
  testthat::expect_no_error(
    my_pgx <- pgx <- playbase::pgx.createPGX(
      counts = playbase::COUNTS,
      samples = playbase::SAMPLES,
      contrasts = playbase::CONTRASTS[,1, drop = FALSE]
    )
  )
  # compute pgx from the created pgx with default methods
  testthat::expect_no_error(
    my_pgx <- playbase::pgx.computePGX(
      my_pgx,
      gx.methods = c("ttest.welch"),
      gset.methods = c("fisher"))
  )
})

