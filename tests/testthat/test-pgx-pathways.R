test_that("pathbankview works correctly", {

  # Test 1: Valid PathBank ID and coloring node
  test_that("Valid PathBank ID with coloring", {
    pb <- "SMP0080852"
    val <- c("PW_C000064" = 2.0)
    result <- pathbankview(pb, val)
    expect_type(result, "character")  # Expect a character string (path to SVG)
    expect_true(file.exists(result))  # Ensure the returned path exists
  })

  # Test 2: Valid PathBank ID with no matching nodes
  test_that("Valid PathBank ID with no matching nodes for coloring", {
    pb <- "SMP0080852"
    val <- c("nonexistent_node" = 1.5)
    result <- pathbankview(pb, val)
    expect_type(result, "character")  # Expect a character string
    expect_true(file.exists(result))  # Ensure the returned path exists
  })

  # Test 3: Invalid PathBank ID
  test_that("Invalid PathBank ID returns NULL", {
    pb <- "INVALID_ID"
    val <- c("PW_C000064" = 2.0)
    result <- pathbankview(pb, val)
    expect_null(result)  # Expect NULL when the ID is invalid
  })

  # Test 4: Valid PathBank ID with NULL val
  test_that("Valid PathBank ID with NULL val", {
    pb <- "SMP0080852"
    val <- NULL
    result <- pathbankview(pb, val)
    expect_type(result, "character")  # Expect a character string
    expect_true(file.exists(result))  # Ensure the returned path exists
  })

  # Test 5: Valid PathBank ID with mixed values in val
  test_that("Valid PathBank ID with mixed values", {
    pb <- "SMP0080852"
    val <- c("PW_C000064" = 1.5, "PW_C000065" = -2.0)
    result <- pathbankview(pb, val)
    expect_type(result, "character")  # Expect a character string
    expect_true(file.exists(result))  # Ensure the returned path exists
  })

  # Test 6: Empty val vector
  test_that("Valid PathBank ID with empty val", {
    pb <- "SMP0080852"
    val <- c()
    result <- pathbankview(pb, val)
    expect_type(result, "character")  # Expect a character string
    expect_true(file.exists(result))  # Ensure the returned path exists
  })
})

test_that("wikipath (modern) works correctly", {

  # Test 1: Valid WikiPathway ID and coloring node
  test_that("Valid WikiPathway ID with coloring", {
    wp <- "WP5405"
    val <- c("CRYL1" = 2.0)
    result <- wikipathview(wp, val)
    expect_type(result, "character")  # Expect a character string (path to SVG)
    expect_true(file.exists(result))  # Ensure the returned path exists
  })

  # Test 2: Valid WikiPathway ID with no matching nodes
  test_that("Valid WikiPathway ID with no matching nodes for coloring", {
    wp <- "WP5405"
    val <- c("node" = 2.0)
    result <- wikipathview(wp, val)
    expect_type(result, "character")  # Expect a character string
    expect_true(file.exists(result))  # Ensure the returned path exists
  })

  # Test 3: Invalid WikiPathway ID
  test_that("Invalid WikiPathway ID returns NULL", {
    wp <- "invalid"
    val <- c("CRYL1" = 2.0)
    result <- wikipathview(wp, val)
    expect_null(result)  # Expect NULL when the ID is invalid
  })

  # Test 4: Valid WikiPathway ID with NULL val
  test_that("Valid WikiPathway ID with NULL val", {
    wp <- "WP5405"
    val <- NULL
    result <- wikipathview(wp, val)
    expect_type(result, "character")  # Expect a character string
    expect_true(file.exists(result))  # Ensure the returned path exists
  })

  # Test 5: Valid WikiPathway ID with mixed values in val
  test_that("Valid WikiPathway ID with mixed values", {
    wp <- "WP5405"
    val <- c("CRYL1" = 2.0, "CRYL2" = 1.0)
    result <- wikipathview(wp, val)
    expect_type(result, "character")  # Expect a character string
    expect_true(file.exists(result))  # Ensure the returned path exists
  })

  # Test 6: Empty val vector
  test_that("Valid WikiPathway ID with empty val", {
    wp <- "WP5405"
    val <- c()
    result <- wikipathview(wp, val)
    expect_type(result, "character")  # Expect a character string
    expect_true(file.exists(result))  # Ensure the returned path exists
  })
})

test_that("wikipath (classic) works correctly", {

  # Test 1: Valid WikiPathway ID and coloring node
  test_that("Valid WikiPathway ID with coloring", {
    wp <- "WP4876"
    val <- c("TRAF3" = 2.0)
    result <- wikipathview(wp, val)
    expect_type(result, "character")  # Expect a character string (path to SVG)
    expect_true(file.exists(result))  # Ensure the returned path exists
  })

  # Test 2: Valid WikiPathway ID with no matching nodes
  test_that("Valid WikiPathway ID with no matching nodes for coloring", {
    wp <- "WP4876"
    val <- c("node" = 2.0)
    result <- wikipathview(wp, val)
    expect_type(result, "character")  # Expect a character string
    expect_true(file.exists(result))  # Ensure the returned path exists
  })

  # Test 3: Invalid WikiPathway ID
  test_that("Invalid WikiPathway ID returns NULL", {
    wp <- "invalid"
    val <- c("TRAF3" = 2.0)
    result <- wikipathview(wp, val)
    expect_null(result)  # Expect NULL when the ID is invalid
  })

  # Test 4: Valid WikiPathway ID with NULL val
  test_that("Valid WikiPathway ID with NULL val", {
    wp <- "WP4876"
    val <- NULL
    result <- wikipathview(wp, val)
    expect_type(result, "character")  # Expect a character string
    expect_true(file.exists(result))  # Ensure the returned path exists
  })

  # Test 5: Valid WikiPathway ID with mixed values in val
  test_that("Valid WikiPathway ID with mixed values", {
    wp <- "WP4876"
    val <- c("TRAF3" = 2.0, "TRAF4" = 1.0)
    result <- wikipathview(wp, val)
    expect_type(result, "character")  # Expect a character string
    expect_true(file.exists(result))  # Ensure the returned path exists
  })

  # Test 6: Empty val vector
  test_that("Valid WikiPathway ID with empty val", {
    wp <- "WP4876"
    val <- c()
    result <- wikipathview(wp, val)
    expect_type(result, "character")  # Expect a character string
    expect_true(file.exists(result))  # Ensure the returned path exists
  })
})
