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
