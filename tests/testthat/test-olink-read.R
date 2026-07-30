testthat::test_that("OlinkAnalyze::read_NPX works with parquet file", {
  parquet_file <- system.file("npx_data_ext.parquet", package = "playbase")
  testthat::expect_true(file.exists(parquet_file))
  npx <- OlinkAnalyze::read_NPX(parquet_file)
  testthat::expect_s3_class(npx, "data.frame")
  testthat::expect_gt(nrow(npx), 0)
})

testthat::test_that("OlinkAnalyze::read_NPX works with CSV file", {
  csv_file <- system.file("abundance_NPX_Data_3K.csv", package = "playbase")
  testthat::expect_true(file.exists(csv_file))
  npx <- OlinkAnalyze::read_NPX(csv_file)
  testthat::expect_s3_class(npx, "data.frame")
  testthat::expect_gt(nrow(npx), 0)
})
