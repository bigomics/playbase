##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

test_that("ai_infographic_set/get round-trip on a plain list", {
  pgx <- list(ai = list(combined = list(report = "## Summary\n\nSome text.")))

  tmp <- tempfile(fileext = ".png")
  writeBin(as.raw(c(1, 2, 3)), tmp)
  result <- list(path = tmp, prompt = "a prompt", metadata = list(model = "img-model"))

  pgx <- ai_infographic_set(pgx, "combined", result, status = "done", style = "bigomics")

  img <- ai_infographic_get(pgx, "combined")
  expect_false(is.null(img))
  expect_equal(img$status, "done")
  expect_equal(img$bytes, as.raw(c(1, 2, 3)))
  expect_equal(img$content_type, "image/png")
  expect_equal(img$model, "img-model")
  expect_equal(img$style, "bigomics")

  unlink(tmp)
})

test_that("ai_infographic_get returns NULL for a slot with no infographic", {
  pgx <- list(ai = list(combined = list(report = "text")))
  expect_null(ai_infographic_get(pgx, "combined"))
  expect_null(ai_infographic_get(pgx, "missing_slot"))
})

test_that("ai_infographic_get surfaces stored error entries", {
  pgx <- list(ai = list(combined = list(report = "text")))
  pgx <- ai_infographic_set(pgx, "combined", result = NULL,
    status = "error", error = "boom")

  img <- ai_infographic_get(pgx, "combined")
  expect_false(is.null(img))
  expect_equal(img$status, "error")
  expect_true(nzchar(img$error))
})

# -----------------------------------------------------------------------------
# pgx.update_infographics() -- omicsai::omicsai_gen_image() mocked out
# -----------------------------------------------------------------------------

.stub_gen_image <- function(template, params = NULL, config = NULL, filename = NULL) {
  tmp <- tempfile(fileext = ".png")
  writeBin(as.raw(c(10, 20, 30)), tmp)
  list(path = tmp, prompt = "stub prompt", metadata = list(model = "stub-model"))
}

.stub_pgx <- function() {
  list(
    organism = "human",
    ai = list(
      combined = list(report = "## Summary\n\nSome report text.")
    )
  )
}

test_that("pgx.update_infographics populates bytes for a slot with a report", {
  testthat::local_mocked_bindings(
    omicsai_gen_image = .stub_gen_image,
    omicsai_image_config = function(...) list(...),
    omicsai_strip_report_noise = function(report) report,
    omicsai_image_species_visual = function(organism) organism,
    image_prompt = function(...) list(...),
    build_prompt = function(prompt) list(system = "sys", board = "board"),
    frag = function(...) NULL,
    .package = "omicsai"
  )

  pgx <- .stub_pgx()
  pgx <- pgx.update_infographics(pgx, ai = list(select = "combined"))

  expect_equal(pgx$ai$combined$infographic$status, "done")
  expect_equal(pgx$ai$combined$infographic$bytes, as.raw(c(10, 20, 30)))
})

test_that("pgx.update_infographics skips a slot with an existing infographic unless forced", {
  testthat::local_mocked_bindings(
    omicsai_gen_image = .stub_gen_image,
    omicsai_image_config = function(...) list(...),
    omicsai_strip_report_noise = function(report) report,
    omicsai_image_species_visual = function(organism) organism,
    image_prompt = function(...) list(...),
    build_prompt = function(prompt) list(system = "sys", board = "board"),
    frag = function(...) NULL,
    .package = "omicsai"
  )

  pgx <- .stub_pgx()
  pgx$ai$combined$infographic <- list(status = "done", bytes = as.raw(99), date = "prior")

  pgx_unforced <- pgx.update_infographics(pgx, ai = list(select = "combined"))
  expect_equal(pgx_unforced$ai$combined$infographic$bytes, as.raw(99))

  pgx_forced <- pgx.update_infographics(pgx, ai = list(select = "combined", force = TRUE))
  expect_equal(pgx_forced$ai$combined$infographic$bytes, as.raw(c(10, 20, 30)))
})

test_that("pgx.update_infographics skips a slot with no report", {
  testthat::local_mocked_bindings(
    omicsai_gen_image = .stub_gen_image,
    omicsai_image_config = function(...) list(...),
    omicsai_strip_report_noise = function(report) report,
    omicsai_image_species_visual = function(organism) organism,
    image_prompt = function(...) list(...),
    build_prompt = function(prompt) list(system = "sys", board = "board"),
    frag = function(...) NULL,
    .package = "omicsai"
  )

  pgx <- list(organism = "human", ai = list(wgcna = list()))
  pgx <- pgx.update_infographics(pgx, ai = list(select = "wgcna"))

  expect_null(pgx$ai$wgcna$infographic)
})

test_that("pgx.update_infographics returns pgx unchanged when ai is NULL", {
  pgx <- .stub_pgx()
  expect_identical(pgx.update_infographics(pgx, ai = NULL), pgx)
})
