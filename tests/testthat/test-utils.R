#' Test for mem.proc
test_that("mem.proc returns default value on non-Linux systems", {
  # Skip test if running on Linux
  if (Sys.info()["sysname"] == "Linux") {
    skip("Test not applicable on Linux systems")
  }

  # Test function
  result <- mem.proc()

  # Check result
  expect_equal(result, "[? MB]")
})

#' Test for info
#'
#'

#' Test for dbg
