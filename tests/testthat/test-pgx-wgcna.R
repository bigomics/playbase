## Filter exemption in wgcna.computeModuleEnrichment: which geneset families
## survive the reference-geneset name filter. Mocks WGCNAplus to capture the
## filter that is actually handed down.

.filter_passed_to_wgcnaplus <- function(gsets, filter = "PATHWAY|HALLMARK|^GO|^C[1-9]") {
  GMT <- matrix(0, nrow = 1, ncol = length(gsets), dimnames = list("A1BG", gsets))
  seen <- NULL
  testthat::local_mocked_bindings(
    computeModuleEnrichment = function(..., filter) {
      seen <<- filter
      NULL
    },
    .package = "WGCNAplus"
  )
  wgcna.computeModuleEnrichment(wgcna = NULL, GMT = GMT, filter = filter)
  seen
}

test_that("CUSTOM and TEST genesets are always exempt from the filter", {
  f <- .filter_passed_to_wgcnaplus(c("HALLMARK:APOPTOSIS", "CUSTOM:my_set", "TEST:set1"))
  expect_true(all(grepl(f, c("HALLMARK:APOPTOSIS", "CUSTOM:my_set", "TEST:set1"))))
})

test_that("METABOLITE ontology/chemclass sets are dropped when pathways match", {
  f <- .filter_passed_to_wgcnaplus(
    c("METABOLITE_PATHWAY:glycolysis", "METABOLITE_ONTOLOGY:x", "METABOLITE_CHEMCLASS:y")
  )
  expect_true(grepl(f, "METABOLITE_PATHWAY:glycolysis"))
  expect_false(grepl(f, "METABOLITE_ONTOLOGY:x"))
  expect_false(grepl(f, "METABOLITE_CHEMCLASS:y"))
})

test_that("METABOLITE sets are kept as fallback when nothing matches the filter", {
  f <- .filter_passed_to_wgcnaplus(c("METABOLITE_ONTOLOGY:x", "METABOLITE_CHEMCLASS:y"))
  expect_true(all(grepl(f, c("METABOLITE_ONTOLOGY:x", "METABOLITE_CHEMCLASS:y"))))
})

test_that("filter = NULL is passed through untouched", {
  expect_null(.filter_passed_to_wgcnaplus("HALLMARK:APOPTOSIS", filter = NULL))
})

test_that("consensus WGCNA takes its reference sets from WGCNAplus::getLibxGMT", {
  libx_gmt <- Matrix::sparseMatrix(
    i = 1:2, j = c(1, 1), x = 1, dims = c(2, 1),
    dimnames = list(c("A1BG", "A2M"), "MOCK:set1")
  )
  seen <- NULL
  testthat::local_mocked_bindings(
    getLibxGMT = function() libx_gmt,
    runConsensusWGCNA = function(..., GMT) {
      seen <<- GMT
      NULL
    },
    .package = "WGCNAplus"
  )
  wgcna.runConsensusWGCNA(exprList = list(), phenoData = list(), compute.enrichment = TRUE)
  expect_true("MOCK:set1" %in% colnames(seen))
})

test_that("WGCNAplus::getLibxGMT returns a non-empty genes x sets matrix", {
  skip_if(is.null(tryCatch(omicspgx::libx_dir(), error = function(e) NULL)))
  GMT <- WGCNAplus::getLibxGMT()
  expect_true(nrow(GMT) > 0 && ncol(GMT) > 0)
  expect_true(any(grepl("^[A-Za-z]", rownames(GMT))))
})
