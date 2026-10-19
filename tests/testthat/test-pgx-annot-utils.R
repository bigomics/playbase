## Mouse and rat OrgDb from lib-x

# the configured lib-x release (gates set LIBX_DIR; tests skip without one)
libx_next <- Sys.getenv("LIBX_DIR", unset = "")

orgdb_organism <- function(db) {
  meta <- AnnotationDbi::metadata(db)
  meta$value[meta$name == "ORGANISM"]
}

test_that(".getOrgDb takes the mouse and rat OrgDb from lib-x", {
  withr::local_envvar(LIBX_DIR = libx_next)
  skip_if(is.null(tryCatch(omicspgx::libx_orgdb("Mus musculus"), error = function(e) NULL)))

  for (organism in c("Mus musculus", "Rattus norvegicus")) {
    db <- playbase:::.getOrgDb(organism)
    expect_s4_class(db, "OrgDb")
    expect_equal(orgdb_organism(db), organism)
  }
})

test_that(".getOrgDb falls back to AnnotationHub without lib-x", {
  skip_if_not_installed("AnnotationHub")
  withr::local_envvar(LIBX_DIR = NA)
  withr::local_options(omicspgx.libx_dir = NULL)
  testthat::local_mocked_bindings(
    AnnotationHub = function(...) stop("offline"),
    .package = "AnnotationHub"
  )
  for (organism in c("Mus musculus", "Rattus norvegicus")) {
    suppressWarnings(
      expect_error(playbase:::.getOrgDb(organism), "offline")
    )
  }
})

## Human OrgDb

test_that(".getOrgDb returns the human OrgDb from lib-x first", {
  sentinel <- structure(list(), class = "sentinel_orgdb")
  testthat::local_mocked_bindings(
    libx_orgdb = function(organism, ...) sentinel,
    .package = "omicspgx"
  )
  expect_identical(playbase:::.getOrgDb("Homo sapiens"), sentinel)
})

test_that(".getOrgDb gives an OrgDb for human with lib-x", {
  skip_if(libx_next == "")
  withr::local_envvar(LIBX_DIR = libx_next)
  skip_if(is.null(tryCatch(omicspgx::libx_orgdb("Homo sapiens"), error = function(e) NULL)))
  db <- playbase:::.getOrgDb("Homo sapiens")
  expect_s4_class(db, "OrgDb")
  expect_equal(orgdb_organism(db), "Homo sapiens")
})

test_that(".getOrgDb uses the org.Hs.eg.db package object without lib-x", {
  skip_if_not_installed("org.Hs.eg.db")
  testthat::local_mocked_bindings(
    libx_orgdb = function(organism, ...) stop("no lib-x"),
    .package = "omicspgx"
  )
  pkg_db <- org.Hs.eg.db::org.Hs.eg.db
  skip_if(inherits(tryCatch(AnnotationDbi::keytypes(pkg_db), error = function(e) e), "error"))
  expect_identical(playbase:::.getOrgDb("Homo sapiens"), pkg_db)
})
