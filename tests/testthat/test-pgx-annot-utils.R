## Optional Plasmodium OrgDb

test_that("missing optional Pf OrgDb reports the existing unsupported result", {
  queried <- character()
  testthat::local_mocked_bindings(
    AnnotationHub = function(...) list(),
    query = function(ah, pattern, ...) {
      queried <<- pattern
      list()
    },
    .package = "AnnotationHub"
  )
  get_db <- playbase:::.getOrgDb
  environment(get_db) <- list2env(
    list(requireNamespace = function(package, quietly) {
      expect_identical(package, "org.Pf.plasmo.db")
      expect_true(quietly)
      FALSE
    }),
    parent = environment(get_db)
  )

  for (flag in list(NULL, FALSE, TRUE)) {
    expect_message(
      expect_null(get_db("Plasmodium falciparum", use.ah = flag)),
      "WARNING: organism 'Plasmodium falciparum' not in AnnotationHub.",
      fixed = TRUE
    )
    expect_identical(queried, c("Plasmodium falciparum", "OrgDb"))
  }
})
