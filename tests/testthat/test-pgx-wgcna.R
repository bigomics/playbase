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
