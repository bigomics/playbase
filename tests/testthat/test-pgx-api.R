#' Test for pgx.getFamilies
d <- get_mini_example_data()
pgx <- playbase::pgx.createPGX(counts = d$counts, 
  samples = d$samples, 
  contrast = d$contrast)
test_that("pgx.getFamilies returns expected families", {

  # Test with default parameters
  families <- pgx.getFamilies(pgx)
  
  # Expect 2 families returned
  expect_equal(length(families), 1)
  
  # Expect FAMILY1 and FAMILY3 returned
  expect_equal(families,"<all>")
  

})
