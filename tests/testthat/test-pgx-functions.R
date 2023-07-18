#' Test for pgx.phenoMatrix
#' 
#' 

#' Test for text_repel.NOTWORKING
#' 
#' 

#' Test for pos.compact
test_that("pos.compact works", {
     # Sample input
     set.seed(123)
     pos <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), ncol = 2)

     
     # Expected output
     expected <- matrix(c(1.00, 1.04, 1.08, 1.12, 1.16, 
                        6.00, 6.04, 6.08, 6.12, 6.16), ncol = 2)
    # Run function
    result <- pos.compact(pos)
    expect_equal(class(result), c("matrix", "array"))
    expect_equal(dim(result), dim(pos))
    expect_equal(result, expected, tolerance = 1e-3)

})

#' Test for util.findboxes
#' 
#' 

#' Test for selectSamplesFromSelectedLevels
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


#' Test for star.symbols
test_that("star.symbols generates correct symbols", {

  # Test cases
  expect_equal(star.symbols(0), "")
  expect_equal(star.symbols(1), "\u2605") 
  expect_equal(star.symbols(5), "\u2605\u2605\u2605\u2605\u2605")
  
  # Default pch
  expect_equal(star.symbols(3), "\u2605\u2605\u2605")
  
  # Custom pch
  expect_equal(star.symbols(2, pch = "X"), "XX")
  
})

#' Test for search_path
test_that("search_path finds file in paths", {
  # Create temp dirs and files
  dir1 <- tempdir()
  dir2 <- tempdir()
  file <- "test.txt"
  
  file.create(file.path(dir1, file))
  
  # Test
  path <- search_path(c(dir1, dir2), file)
  
  # Assertions
  expect_equal(path, file.path(dir1, file))
  expect_true(file.exists(path))
  
  # Clean up
  unlink(file.path(dir1, file))
})

test_that("search_path returns NULL if file not found", {

  # Create temp dir
  dir1 <- tempfile()
  
  # Test
  path <- search_path(dir1, "not_found.txt")
  
  # Assertion
  expect_null(path)
  
})


#' Test for rowscale
test_that("rowscale scales rows correctly", {

  # Generate test matrix
  m <- matrix(c(1,2,3,4,5,
                6,7,8,9,10,
                11,12,13,14,
                15,16,17,18,19,20), nrow = 5)
  
  # Expected result
  expected_vector <- rep(c(-1.341617, -0.4472056, 0.4472056, 1.341617), 5)
  expected <- matrix(expected_vector, nrow = 5, byrow = TRUE)
  
  # Test rowscale 
  result <- rowscale(m)

  # Check class
  expect_equal(class(result), c("matrix", "array"))

  # Check dimensions
  expect_equal(dim(result), dim(m))

  # Check values
  expect_equal(result, expected, tolerance = 1e-6)

})

# Test for add_opacity
test_that("add_opacity adds opacity correctly", {

  # Input vectors
  hexcol <- c("#FFFFFF", "#000000", NA)
  opacity <- 0.5
  
  # Expected output
  expected <- c(rgb(1, 1, 1, 0.5), rgb(0, 0, 0, 0.5), NA)
  
  # Run function
  result <- add_opacity(hexcol, opacity)
  
  # Check class
  expect_equal(class(result), "character")
  
  # Check length 
  expect_equal(length(result), length(hexcol))

  # Check values
  expect_equal(result, expected)
  
})

#' Test for logCPM
#' 

#' Test for pgx.checkObject 
#' 

#' Test for matGroupMeans
#' 
#' 


#' Test for knnMedianFilter
#' 
#' 

#' Test for nmfImpute
#' 
#' 

#' Test for knnImputeMissing
#' 
#' 

#' Test for randomImputeMissing
#' 
#' 

#' Test for human2mouse.SLLOWWW
#' 
#' 

#' Test for human2mouse
#' 
#' 


#' Test for mouse2human
#' 
#' 

#' Test for probe2symbol
#' 
#' 

#' Test for trimsame
#'
#'

#' Test for trimsame.ends 
#' 
#' 

#' Test for trimsame0
#' 
#' 

#' Test for dbg.BAK
#' 
#' 

#' Test for read.csv3.BAK 
#' 
#' 

#' Test for read.csv3
#' 
#' 

#' Test for read.as_matrix.SAVE 
#' 
#' 

#' Test for read.as_matrix 
test_that("read.as_matrix reads file as matrix", {
  # Create temp file
  tmp <- tempfile()
  writeLines(c(",1,2,3", "gene1,1,2,3", "gene2,4,5,6"), tmp)
  
  # Expected matrix
  expected <- matrix(c(1,2,3,4,5,6), ncol=3, byrow=TRUE,
                     dimnames=list(c("gene1", "gene2"), NULL))
  colnames(expected) <- c(1, 2, 3)
  # Test function
  result <- read.as_matrix(tmp)

  # Check class
  expect_equal(class(result), c("matrix", "array"))

  # Check dimensions
  expect_equal(dim(result), dim(expected)) 
  
  # Check values
  expect_equal(result, expected)

  # Clean up
  unlink(tmp)
})


#' Test for fread.csv
#' 
#' 

#' Test for tagDuplicates
test_that("tagDuplicates tags duplicates correctly", {

  # Test input
  s <- c("A", "B", "C", "A", "B")
  
  # Expected output
  expected <- c("A", "B", "C", "A ", "B ")
  
  # Run function
  result <- tagDuplicates(s)

  # Check class
  expect_equal(class(result), "character")

  # Check length 
  expect_equal(length(result), length(s))

  # Check tagged values
  expect_equal(result, expected)
  
})

#' Test for wrapHyperLink
#' 
#' 

#' Test for reverse.AvsB
test_that("reverse.AvsB reverses A vs B contrasts", {

  comps <- c("A_vs_B", "grp1_vs_grp2", "cond1:A_vs_B@time1")
  
  expected <- c("B_vs_A", "grp2_vs_grp1", "cond1:B_vs_condA@time1")
  
  result <- reverse.AvsB(comps)
  
  expect_equal(class(result), "character")
  expect_equal(length(result), length(comps))
  expect_equal(result, expected)
  
})

#' Test for is.POSvsNEG
#' 
#' 

#' Test for getLevels
#' 
#' 


#' Test for getMyGeneInfo
#' 
#' 

#' Test for getHSGeneInfo
#' 
#' 

#' Test for pgx.getGeneFamilies
#' 
#' 

#' Test for pgx.getGeneSetCollections
#' 
#' 

#' Test for filterFamily
#' 

#' Test for filterProbes
#'
#' 

#' Test for computeFeatureScore
#' 
#' 

#' Test for getKeggID
#' 
#' 

#' Test for is.Date
test_that("is.Date correctly identifies date columns", {
  # Test data
  x1 <- c("2022-01-01", "2022-02-01", "2022-03-01")
  x2 <- c("01/01/2022", "02/01/2022", "03/01/2022")
  x3 <- c("2022/01/01", "2022/02/01", "2022/03/01")
  x4 <- c("not a date", "also not a date", "still not a date")
  x5 <- c("abc", "123")
  x6 <- c(NA, NA, NA)

  # Case 1: Date column in YYYY-MM-DD format
  expect_true(is.Date(x1))
  
  # Case 2: Date column in DD/MM/YYYY format
  expect_true(is.Date(x2))
  
  # Case 3: Date column in YYYY/MM/DD format
  expect_true(is.Date(x3))
  
  # Case 4: Non-date column
  expect_false(is.Date(x4))
  expect_false(is.Date(x5))
  expect_false(is.Date(x6))


})

#' Test for averageByGroup
#' 
#' 

#' Test for makeAcronym
test_that("makeAcronym generates expected acronyms", {

  # Test single word
  expect_equal(makeAcronym("hello"), "he")
  
  # Test multiple words
  expect_equal(makeAcronym("hello_world"), "HW")
  
  # Test hyphenated words
  expect_equal(makeAcronym("hello-world"), "HW")

  # Test already uppercase
  expect_equal(makeAcronym("HELLO"), "HE")

# Test data
  x1 <- "hello"
  x2 <- "hello_world"
  x3 <- "hello-world"
  x4 <- c("hello", "world")
  x5 <- "HELLO"
  x6 <- c("United Nations", "European Union")
  
  # Expected results
    expected1 <- "he"
    expected2 <- "HW"
    expected3 <- "HW"
    expected4 <- c("he", "wo")
    expected5 <- "HE" 
    expected6 <- c("UN", "EU")
  
  # Test function
  result1 <- makeAcronym(x1)
  result2 <- makeAcronym(x2)
  result3 <- makeAcronym(x3)
  result4 <- makeAcronym(x4)
  result5 <- makeAcronym(x5)
  result6 <- makeAcronym(x6)
  
  # Check class
  expect_equal(class(result1), "character")
  expect_equal(class(result2), "character")
  expect_equal(class(result3), "character")
  expect_equal(class(result4), "character")
  expect_equal(class(result5), "character")
  expect_equal(class(result6), "character")
  
  # Check length
  expect_equal(length(result1), length(expected1))
  expect_equal(length(result2), length(expected2))
  expect_equal(length(result3), length(expected3))
  expect_equal(length(result4), length(expected4))
  expect_equal(length(result5), length(expected5))
  expect_equal(length(result6), length(expected6))
  
  # Check values
  ## Single word
  expect_equal(result1, expected1)
  ## Single word (underscore)
  expect_equal(result2, expected2)
  ## Single word (dash)
  expect_equal(result3, expected3)
  ## Multiple words
  expect_equal(result4, expected4)
  ## Single words capital letters
  expect_equal(result5, expected5)
  ## Multiple items vector, multiple words 
  expect_equal(result6, expected6)
  
})

#' Test for relevelFactorFirst
test_that("relevelFactorFirst relevels factor with first level", {

  # Input factor
  f <- factor(c("B", "C", "A"))
  
  # Expected releveled factor
  expected <- factor(c("B", "C", "A"), levels = c("B", "C", "A"))

  # Test function
  result <- relevelFactorFirst(f)

  # Check class
  expect_equal(class(result), "factor")

  # Check levels
  expect_equal(levels(result), levels(expected))

  # Check values
  expect_equal(as.character(result), as.character(expected))
  
})

#' Test for extremeCorrelation
#' 
#' 

#' Test for alias2hugo
#' 
#' 

#' Test for breakstringBROKEN
#' 
#' 

#' Test for breakstring
#' 
#' 

#' Test for breakstring2
#' 
#' 

#' Test for shortstring
#' 
#' 

#' Test for shortstring0
#' 
#' 

#' Test for psort 
test_that("psort sorts data frame by p-value column", {

  # Create test data frame
  df <- data.frame(
    Gene = c("A", "B", "C"),
    pvalue = c(0.04, 0.02, 0.01),
    FC = c(1.5, 2.0, 1.2)
  )

  # Expected order
  expected_order <- c("C", "B", "A")

  # Sort
  result <- psort(df)

  # Test class
  expect_equal(class(result), "data.frame")

  # Test sorting
  expect_equal(result$Gene, expected_order)

  # Test sorting by specified column
  result2 <- psort(df, p.col = "pvalue")
  expect_equal(result2$Gene, expected_order)

})

#' Test for tidy.dataframe
#' 
#' 

#' Test for param.class 
#' 
#' 

#' Test for is.num

test_that("is.num correctly identifies numeric vectors", {
  # Numeric vector
  expect_true(is.num(c(1, 2, 3)))

  # Character vector
  expect_false(is.num(c("a", "b", "c")))

  # Mixed vector
  expect_true(is.num(c(1, 2, "a")))

  # All NA
  expect_false(is.num(rep(NA, 10)))

  # Mostly unique 
  expect_false(is.num(sample(letters, 20, replace=FALSE)))

  # Mostly duplicates
  expect_false(is.num(sample(1:3, 20, replace=TRUE)))
})

#' Test for isanumber
test_that("isanumber correctly identifies numeric vectors", {

  numeric_vector <- c(1, 2, 3)
  expect_true(isanumber(numeric_vector))

  character_vector <- c("a", "b", "c") 
  expect_false(isanumber(character_vector))

  mixed_vector <- c(1, 2, "a")
  expect_true(isanumber(mixed_vector))

  mostly_na <- c(NA, NA, 1, 2)
  expect_true(isanumber(mostly_na))

  empty <- c()
  expect_false(isanumber(empty))
  
})

#' Test for expandAnnotationMatrix
#' 
#' 

#' Test for expandAnnotationMatrixSAVE
#' 
#' 

#' Test for expandPhenoMatrix
#' 
#' 

#' Test for correctMarchSeptemberGenes
#' 
#' 

#' Test for cor.pvalue
test_that("cor.pvalue calculates p-values correctly", {

  # Sample input
  x <- 0.5 
  n <- 100
  
  # Expected p-value
  expected_p <- cor.pvalue(x, n)
  
  # Calculate p-value
  result_p <- cor.pvalue(x, n)
  
  # Test class
  expect_equal(class(result_p), "numeric")
  
  # Test length 
  expect_equal(length(result_p), 1)

  # Test value
  expect_equal(result_p, expected_p, tolerance = 1e-5)

  # Test different inputs
  x <- c(-0.8, 0.3)
  n <- c(50, 150)
  expected_p <- cor.pvalue(x, n)
  result_p <- cor.pvalue(x, n)
  
  expect_equal(length(result_p), 2)
  expect_equal(result_p, expected_p, tolerance = 1e-5)
  
})

#' Test for getGSETS_playbase
#' 
#' 

#' Test for getGSETS_playbase.SAVE
#' 
#' 
