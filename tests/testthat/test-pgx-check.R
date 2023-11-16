#' Test for pgx.checkINPUT
# Load the necessary libraries

# Define the common test data
df <- data.frame(matrix(rnorm(20), nrow=10))
colnames(df) <- c("sample1", "sample2")
rownames(df) <- paste0("gene", 1:10)

# Test 1: 
test_that("checkINPUT function returns a list", {
  expect_is(playbase::pgx.checkINPUT(df, "SAMPLES"), "list")
})

# Test 2: 
test_that("checkINPUT function handles non-existent datatype", {
  expect_error(playbase::pgx.checkINPUT(df, "NON_EXISTENT"))
})

# Test 3: 
test_that("checkINPUT function handles duplicated column names", {
  df_dup <- cbind(df, df)
  check_result <- playbase::pgx.checkINPUT(df_dup, "COUNTS")
  expect_false(check_result$PASS)
  expect_equal(check_result$checks$e6, c("sample1", "sample2"))
})

# Test 4: 
test_that("checkINPUT function handles duplicated row names", {
  df_dup <- as.matrix(df)
  rownames(df_dup)[1] <- rownames(df_dup)[2]

  # No pass for samples
  samples_check_result <- playbase::pgx.checkINPUT(df_dup, "SAMPLES")
  expect_false(samples_check_result$PASS)
  expect_equal(samples_check_result$checks$e1, "gene2")

  # Pass for counts
  counts_check_result <- playbase::pgx.checkINPUT(df_dup, "COUNTS")
  expect_true(counts_check_result$PASS)
  expect_equal(counts_check_result$checks$e7, "gene2")

})

# Test 5: 
test_that("checkINPUT function handles zero count rows", {
  df_zero <- df
  df_zero[c(1,5),] <- 0
  check_result <- playbase::pgx.checkINPUT(df_zero, "COUNTS")
  expect_true(check_result$PASS)  
  expect_equal(check_result$checks$e9, c("gene1", "gene5"))
})

# Test 6: 
test_that("checkINPUT function handles zero count columns", {
  df_zero <- df
  df_zero[,1] <- 0
  counts_zero_col_check <- playbase::pgx.checkINPUT(df_zero, "COUNTS")
  expect_true(counts_zero_col_check$PASS)
  expect_equal(counts_zero_col_check$checks$e10, "sample1")
})

# Test 7: 
test_that("checkINPUT function handles valid contrast names", {
  df_contrasts <- df
  colnames(df_contrasts) <- c("sample1_sample2", "sample2_sample3")
  expect_error(pgx.checkINPUT(df_contrasts, "CONTRASTS"))
})

# Test 8: 
test_that("checkINPUT function handles invalid contrast names", {
  df_contrasts <- df
  colnames(df_contrasts) <- c("sample1_sample2", "sample2_sample3")
  expect_error(playbase::pgx.checkINPUT(df_contrasts, "CONTRASTS"))
})

#' Test for pgx.crosscheckINPUT
# Define the test data
COUNTS <- data.frame(matrix(rnorm(40), nrow=10))
colnames(COUNTS) <- c("sample1", "sample2", "sample3", "sample4")
rownames(COUNTS) <- paste0("gene", 1:10)
SAMPLES <- data.frame(type  = c("A", "A", "B", "B"),
                      group = c("high", "low", "high", "low"))
rownames(SAMPLES) <- colnames(COUNTS)
CONTRASTS <- data.frame(`type:A_vs_B` = c("A", "A", "B", "B"),
                        `group:high_vs_low` = c("high", "low", "high", "low"))

# Test 2: Check if the function returns a list
test_that("crosscheckINPUT function returns a list", {
  result <- playbase::pgx.crosscheckINPUT(SAMPLES, COUNTS, CONTRASTS)
  expect_type(result, "list")
  expect_equal(result$checks, list())
  expect_true(result$PASS)

})

# Test 3: Check if the function handles non-matching samples
test_that("crosscheckINPUT function handles non-matching sample and count names", {
  SAMPLES_mismatch <- SAMPLES
  rownames(SAMPLES_mismatch)[1] <- "mismatch"
  suppressWarnings(result <- playbase::pgx.crosscheckINPUT(SAMPLES_mismatch, COUNTS, CONTRASTS))
  expect_true(result$PASS)
  expect_equal(result$checks$e19, c("mismatch", "sample1"))
  expect_equal(result$checks$e17, c("sample2", "sample3", "sample4", "1", "2", "3", "4"))

})

# Test 6: Check if the function handles non-matching order of sample and count names
test_that("crosscheckINPUT function handles non-matching order of sample and count names", {
  COUNTS_mismatch_order <- COUNTS[, 2:1]
  result <- playbase::pgx.crosscheckINPUT(SAMPLES, COUNTS_mismatch_order, CONTRASTS)
  expect_true(result$PASS)
  expect_equal(result$checks$e19, c("sample3", "sample4"))
  expect_equal(result$checks$e17, c("sample1", "sample2", "1", "2", "3", "4"))

})

#' Test for contrasts_conversion_check