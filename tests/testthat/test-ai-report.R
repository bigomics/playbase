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

test_that(".pathways_meta_fx recovers direction when meta.fx is NaN", {
  # Reproduces the GSEA "no enrichment" bug: meta.fx is NaN (member-gene average
  # undefined for sparse/multi-omics sets) while per-method effects exist.
  mx <- list(
    meta.fx = c(GS1 = NaN, GS2 = NaN, GS3 = NaN),
    meta.q  = c(GS1 = 0.001, GS2 = 0.01, GS3 = 0.2),
    fc = matrix(c(1.2, -0.8, NA, 0.9, -1.1, NA), nrow = 3,
                dimnames = list(c("GS1", "GS2", "GS3"), c("gsva", "fgsea")))
  )
  eff <- .pathways_meta_fx(mx, rep(TRUE, 3))

  expect_true(is.finite(eff[["GS1"]]) && eff[["GS1"]] > 0)  # elevated
  expect_true(is.finite(eff[["GS2"]]) && eff[["GS2"]] < 0)  # suppressed
  expect_true(is.nan(eff[["GS3"]]))                         # no per-method fc

  # The report's significance filter now counts sets previously dropped.
  significant <- is.finite(eff) & is.finite(mx$meta.q) & mx$meta.q < 0.05
  expect_equal(sum(significant), 2L)
})

test_that(".pathways_meta_fx leaves finite meta.fx untouched", {
  mx <- list(meta.fx = c(A = 1.5, B = -2.0), fc = NULL)
  expect_equal(.pathways_meta_fx(mx, c(TRUE, TRUE)), c(A = 1.5, B = -2.0))
})
