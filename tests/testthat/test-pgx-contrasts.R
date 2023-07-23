
#' Test for contrastAsLabels
test_that("contrastAsLabels works as intended", {

  contrast <- playbase::CONTRASTS
  
  # Run function
  result <- playbase::contrastAsLabels(contrast)
  
  # Check class
  expect_s3_class(result, "data.frame") 
  
  # Check dimensions
  expect_equal(dim(result), dim(z))

  # Check preserve rownames 
  expect_true(all(rownames(contrast) == rownames(z)))
  expect_true(all(colnames(contrast) == colnames(z)))

  # Check 0 become NAs
  expect_equal(sum(contrast == 0), sum(is.na(z)))

})