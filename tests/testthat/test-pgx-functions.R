testthat::test_that("selectSamplesFromSelectedLevels returns correct samples", {

  # Test data
  Y <- data.frame(
    Pheno1 = c("A", "B", "C"),
    Pheno2 = c(1, 2, 3)
  )
  rownames(Y) <- c("sample_1", "sample_2", "sample_3")
  # Case 1: No levels
  testthat::expect_equal(selectSamplesFromSelectedLevels(Y, levels = NULL), rownames(Y))

  # Case 2: Empty levels
  testthat::expect_equal(selectSamplesFromSelectedLevels(Y, levels = ""), rownames(Y))

  # Case 3: Single level
  levels <- "Pheno1=B"
  testthat::expect_equal(selectSamplesFromSelectedLevels(Y, levels), "sample_2")

  # Case 4: Multiple levels
  levels <- c("Pheno1=A", "Pheno2=3")
  testthat::expect_equal(selectSamplesFromSelectedLevels(Y, levels), c("sample_1", "sample_3"))

})