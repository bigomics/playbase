is_ci <- function() {
  Sys.getenv("CI") == "true"
}

## TRUE when PathBank serves a diagram for `pb`; it sometimes answers with an
## empty SVG document.
pathbank_serves <- function(pb) {
  url <- paste0("https://www.pathbank.org/view/", pb, "/download?type=simple_large_font_vector_image")
  destfile <- tempfile(fileext = ".svg")
  tryCatch(
    {
      utils::download.file(url, destfile, quiet = TRUE)
      length(xml2::xml_children(xml2::read_xml(destfile))) > 0
    },
    error = function(e) FALSE
  )
}

test_that("pathbankview works correctly", {
  skip_if(is_ci(), "Skipping PathBank tests in CI environment (GitHub Actions randomyl fails to download SVGs)")
  skip_if_not(pathbank_serves("SMP0080852"), "PathBank serves no diagram for SMP0080852")

  # Test 1: Valid PathBank ID and coloring node
  test_that("Valid PathBank ID with coloring", {
    pb <- "SMP0080852"
    val <- c("PW_C000064" = 2.0)
    result <- pathbankview(pb, val)
    expect_type(result, "character") # Expect a character string (path to SVG)
    expect_true(file.exists(result)) # Ensure the returned path exists
  })

  # Test 2: Valid PathBank ID with no matching nodes
  test_that("Valid PathBank ID with no matching nodes for coloring", {
    pb <- "SMP0080852"
    val <- c("nonexistent_node" = 1.5)
    result <- pathbankview(pb, val)
    expect_type(result, "character") # Expect a character string
    expect_true(file.exists(result)) # Ensure the returned path exists
  })

  # Test 3: Invalid PathBank ID
  test_that("Invalid PathBank ID returns NULL", {
    pb <- "INVALID_ID"
    val <- c("PW_C000064" = 2.0)
    result <- pathbankview(pb, val)
    expect_null(result) # Expect NULL when the ID is invalid
  })

  # Test 4: Valid PathBank ID with NULL val
  test_that("Valid PathBank ID with NULL val", {
    pb <- "SMP0080852"
    val <- NULL
    result <- pathbankview(pb, val)
    expect_type(result, "character") # Expect a character string
    expect_true(file.exists(result)) # Ensure the returned path exists
  })

  # Test 5: Valid PathBank ID with mixed values in val
  test_that("Valid PathBank ID with mixed values", {
    pb <- "SMP0080852"
    val <- c("PW_C000064" = 1.5, "PW_C000065" = -2.0)
    result <- pathbankview(pb, val)
    expect_type(result, "character") # Expect a character string
    expect_true(file.exists(result)) # Ensure the returned path exists
  })

  # Test 6: Empty val vector
  test_that("Valid PathBank ID with empty val", {
    pb <- "SMP0080852"
    val <- c()
    result <- pathbankview(pb, val)
    expect_type(result, "character") # Expect a character string
    expect_true(file.exists(result)) # Ensure the returned path exists
  })
})

test_that("wikipath (modern) works correctly", {
  skip_if(is_ci(), "Skipping WikiPathways modern tests in CI environment")

  # Test 1: Valid WikiPathway ID and coloring node
  test_that("Valid WikiPathway ID with coloring", {
    wp <- "WP5405"
    val <- c("CRYL1" = 2.0)
    result <- wikipathview(wp, val)
    expect_type(result, "character") # Expect a character string (path to SVG)
    expect_true(file.exists(result)) # Ensure the returned path exists
  })

  # Test 2: Valid WikiPathway ID with no matching nodes
  test_that("Valid WikiPathway ID with no matching nodes for coloring", {
    wp <- "WP5405"
    val <- c("node" = 2.0)
    result <- wikipathview(wp, val)
    expect_type(result, "character") # Expect a character string
    expect_true(file.exists(result)) # Ensure the returned path exists
  })

  # Test 3: Invalid WikiPathway ID
  test_that("Invalid WikiPathway ID returns NULL", {
    wp <- "invalid"
    val <- c("CRYL1" = 2.0)
    result <- wikipathview(wp, val)
    expect_null(result) # Expect NULL when the ID is invalid
  })

  # Test 4: Valid WikiPathway ID with NULL val
  test_that("Valid WikiPathway ID with NULL val", {
    wp <- "WP5405"
    val <- NULL
    result <- wikipathview(wp, val)
    expect_type(result, "character") # Expect a character string
    expect_true(file.exists(result)) # Ensure the returned path exists
  })

  # Test 5: Valid WikiPathway ID with mixed values in val
  test_that("Valid WikiPathway ID with mixed values", {
    wp <- "WP5405"
    val <- c("CRYL1" = 2.0, "CRYL2" = 1.0)
    result <- wikipathview(wp, val)
    expect_type(result, "character") # Expect a character string
    expect_true(file.exists(result)) # Ensure the returned path exists
  })

  # Test 6: Empty val vector
  test_that("Valid WikiPathway ID with empty val", {
    wp <- "WP5405"
    val <- c()
    result <- wikipathview(wp, val)
    expect_type(result, "character") # Expect a character string
    expect_true(file.exists(result)) # Ensure the returned path exists
  })
})

test_that("wikipath (classic) works correctly", {
  skip_if(is_ci(), "Skipping WikiPathways classic tests in CI environment")

  # Test 1: Valid WikiPathway ID and coloring node
  test_that("Valid WikiPathway ID with coloring", {
    wp <- "WP4876"
    val <- c("TRAF3" = 2.0)
    result <- wikipathview(wp, val)
    expect_type(result, "character") # Expect a character string (path to SVG)
    expect_true(file.exists(result)) # Ensure the returned path exists
  })

  # Test 2: Valid WikiPathway ID with no matching nodes
  test_that("Valid WikiPathway ID with no matching nodes for coloring", {
    wp <- "WP4876"
    val <- c("node" = 2.0)
    result <- wikipathview(wp, val)
    expect_type(result, "character") # Expect a character string
    expect_true(file.exists(result)) # Ensure the returned path exists
  })

  # Test 3: Invalid WikiPathway ID
  test_that("Invalid WikiPathway ID returns NULL", {
    wp <- "invalid"
    val <- c("TRAF3" = 2.0)
    result <- wikipathview(wp, val)
    expect_null(result) # Expect NULL when the ID is invalid
  })

  # Test 4: Valid WikiPathway ID with NULL val
  test_that("Valid WikiPathway ID with NULL val", {
    wp <- "WP4876"
    val <- NULL
    result <- wikipathview(wp, val)
    expect_type(result, "character") # Expect a character string
    expect_true(file.exists(result)) # Ensure the returned path exists
  })

  # Test 5: Valid WikiPathway ID with mixed values in val
  test_that("Valid WikiPathway ID with mixed values", {
    wp <- "WP4876"
    val <- c("TRAF3" = 2.0, "TRAF4" = 1.0)
    result <- wikipathview(wp, val)
    expect_type(result, "character") # Expect a character string
    expect_true(file.exists(result)) # Ensure the returned path exists
  })

  # Test 6: Empty val vector
  test_that("Valid WikiPathway ID with empty val", {
    wp <- "WP4876"
    val <- c()
    result <- wikipathview(wp, val)
    expect_type(result, "character") # Expect a character string
    expect_true(file.exists(result)) # Ensure the returned path exists
  })
})

test_that("getReactomeSVG serves the lib-x diagram, resolving sub-pathways", {
  skip_if(is.null(omicspgx::libx_dir()))
  svg <- getReactomeSVG("R-HSA-1059683", cache.dir = tempfile("reactome-svg"))
  expect_equal(basename(svg), paste0(reactomeDiagramId("R-HSA-1059683"), ".svg"))
  expect_gt(file.info(svg)$size, 0)
})

test_that("getReactomeSVG returns NULL without a lib-x release", {
  old.opt <- options(omicspgx.libx_dir = NULL)
  old.env <- Sys.getenv("LIBX_DIR")
  on.exit({
    options(old.opt)
    Sys.setenv(LIBX_DIR = old.env)
  })
  Sys.setenv(LIBX_DIR = "")
  expect_length(reactomeDiagrams(), 0)
  expect_null(getReactomeSVG("R-HSA-73843"))
})
