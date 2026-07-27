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


#' Test for reorder_uniprots
test_that("reorder_uniprots tolerates spaces around ';' (client bug)", {
  # Swiss-Prot (Q02978) must outrank TrEMBL (I3L1P8) regardless of spacing.
  expect_equal(reorder_uniprots("I3L1P8;Q02978")$feature, "Q02978;I3L1P8")
  expect_equal(reorder_uniprots("I3L1P8; Q02978")$feature, "Q02978;I3L1P8")
  expect_equal(reorder_uniprots("I3L1P8;   Q02978")$feature, "Q02978;I3L1P8")
  expect_equal(reorder_uniprots(" I3L1P8 ; Q02978 ")$feature, "Q02978;I3L1P8")
})
