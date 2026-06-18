##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

if (file.exists("R/ai-copilot.R")) {
  source("R/ai-copilot.R")
}

test_that("ai.build_report_context reads reports from pgx$ai", {
  pgx <- list(
    name = "example-data",
    description = "T-cell activation",
    datatype = "proteomics",
    organism = "Homo sapiens",
    samples = data.frame(group = c("resting", "activated")),
    contrasts = matrix(
      c(1, -1),
      nrow = 2,
      dimnames = list(c("resting", "activated"), "activated_vs_resting")
    ),
    settings = list(method = "test"),
    drugs = list(
      "L1000_ACTIVITYS_N20D1011" = list(),
      "GDSC/sensitivity" = list()
    ),
    ai = list(
      de = list(report = "DE report"),
      pathways = list(report = "Pathway report"),
      drugs_L1000_ACTIVITYS_N20D1011 = list(report = "L1000 report"),
      drugs_GDSC_sensitivity = list(report = "GDSC report")
    ),
    wgcna = list(report = list(report = "old WGCNA report"))
  )

  txt <- ai.build_report_context(pgx, collate = TRUE)

  expect_true(grepl("## Description", txt, fixed = TRUE))
  expect_true(grepl("T-cell activation", txt, fixed = TRUE))
  expect_true(grepl("## Dataset Info", txt, fixed = TRUE))
  expect_true(grepl("## Differential Expression", txt, fixed = TRUE))
  expect_true(grepl("DE report", txt, fixed = TRUE))
  expect_true(grepl("## Geneset Enrichment", txt, fixed = TRUE))
  expect_true(grepl("Pathway report", txt, fixed = TRUE))
  expect_true(grepl("## L1000 Activity", txt, fixed = TRUE))
  expect_true(grepl("L1000 report", txt, fixed = TRUE))
  expect_true(grepl("## GDSC/sensitivity", txt, fixed = TRUE))
  expect_true(grepl("GDSC report", txt, fixed = TRUE))
  expect_false(grepl("old WGCNA report", txt, fixed = TRUE))
})

test_that("ai.build_report_context honors sections and output shape", {
  pgx <- list(
    name = "example-data",
    ai = list(
      de = list(report = "DE report"),
      combined = list(report = "Summary report")
    )
  )

  txt <- ai.build_report_context(
    pgx,
    sections = c("description", "differential_expression"),
    collate = TRUE
  )

  expect_true(grepl("## Description", txt, fixed = TRUE))
  expect_true(grepl("DE report", txt, fixed = TRUE))
  expect_false(grepl("Summary report", txt, fixed = TRUE))

  parts <- ai.build_report_context(
    pgx,
    sections = c("description", "differential_expression"),
    collate = FALSE
  )
  expect_equal(names(parts), c("description", "differential_expression"))

  short <- ai.build_report_context(pgx, max_chars = 20)
  expect_true(nchar(short) > 20)
  expect_true(grepl("[Report context truncated]", short, fixed = TRUE))
})
