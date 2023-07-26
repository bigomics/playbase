#' Test for seq_limma
test_that("seq_limma runs without errors", {

    # Generate mock data
    set.seed(151)
    gx <- matrix(sample(1:100, 100*20, replace = TRUE), 100, 20)
    rownames(gx) <- replicate(100, sample(paste0(LETTERS, 1:50), 1)) 
    colnames(gx) <- sample(letters, 20)
    y <- sample(letters[1:4], 20, replace = TRUE)
    result <- playbase::seq_limma(gx, y)
  
    expect_type(result, "double")
    expect_equal(class(result), c("matrix", "array"))
    expect_equal(dim(result), c(100, 5))
    expect_equal(sum(result), 6511.413, tolerance = 1e-3)
})

