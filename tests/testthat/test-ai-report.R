##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

test_that(".ai_report_normheadings delegates heading spacing to omicsai", {
  input <- paste(
    "Opening text",
    "## **Discussion**",
    "Discussion text",
    sep = "\n"
  )
  expected <- paste(
    "Opening text",
    "",
    "## **Discussion**",
    "",
    "Discussion text",
    sep = "\n"
  )

  expect_equal(.ai_report_normheadings(input), expected)
})
