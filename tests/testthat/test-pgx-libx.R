## Reference data read from the lib-x release, the code that replaced
## reference objects, and the provenance stamp of createPGX.

test_that("BLUERED returns the blue-grey-red palette", {
  expect_identical(BLUERED(7), c(
    "#3A5FCD", "#738BD5", "#ACB8DD", "#E5E5E5", "#DDB5B5", "#D58585",
    "#CD5555"
  ))
  expect_identical(BLUERED(32), c(
    "#3A5FCD", "#4567CE", "#5070D0", "#5B78D1", "#6681D3", "#718AD4",
    "#7C92D6", "#879BD7", "#92A4D9", "#9DACDA", "#A8B5DC", "#B3BEDE",
    "#BEC6DF", "#C9CFE1", "#D4D8E2", "#DFE0E4", "#E4E0E0", "#E2D7D7",
    "#E1CDCD", "#DFC4C4", "#DEBBBB", "#DCB1B1", "#DAA8A8", "#D99F9F",
    "#D79696", "#D68C8C", "#D48383", "#D37A7A", "#D17070", "#D06767",
    "#CE5E5E", "#CD5555"
  ))
  expect_identical(BLUERED(33), c(
    "#3A5FCD", "#4467CE", "#4F6FD0", "#5A78D1", "#6480D3", "#6F88D4",
    "#7A91D6", "#8499D7", "#8FA2D9", "#9AAADA", "#A4B2DC", "#AFBBDD",
    "#BAC3DF", "#C4CBE0", "#CFD4E2", "#DADCE3", "#E5E5E5", "#E3DCDC",
    "#E2D3D3", "#E0CACA", "#DFC1C1", "#DDB8B8", "#DCAFAF", "#DAA6A6",
    "#D99D9D", "#D79494", "#D68B8B", "#D48282", "#D37979", "#D17070",
    "#D06767", "#CE5E5E", "#CD5555"
  ))
})

test_that(".pb_libx errors on an unknown name", {
  skip_if(is.null(omicspgx::libx_dir()))
  expect_error(.pb_libx("no_such_object"), "unknown lib-x name")
  expect_error(.pb_libx("fam"), "unknown lib-x name")
})

test_that("getGSETS returns the genes of each named gene set", {
  skip_if(is.null(omicspgx::libx_dir()))
  igsets <- .pb_libx("igsets")
  gs <- names(igsets)[c(1, 100, 1000)]
  res <- getGSETS(gs)
  expect_type(res, "list")
  expect_named(res, gs)
  expect_identical(lengths(res, use.names = FALSE), lengths(igsets[gs], use.names = FALSE))
  expect_type(res[[1]], "character")
})

test_that("pgx.createPGX stamps the lib-x provenance", {
  skip_if(is.null(omicspgx::libx_dir()))
  ## The cell-cycle and gender stage stands in for any lib-x read during createPGX.
  testthat::local_mocked_bindings(
    compute_cellcycle_gender = function(pgx, ...) {
      .pb_libx("families")
      pgx
    },
    .package = "playbase"
  )
  pgx <- suppressMessages(playbase::pgx.createPGX(
    counts = as.matrix(playbase::COUNTS),
    samples = playbase::SAMPLES,
    contrasts = playbase::CONTRASTS,
    organism = "Human",
    datatype = "RNA-seq",
    preprocess = list(remove_outliers = FALSE, impute = FALSE),
    add.gmt = FALSE
  ))
  man <- omicspgx::libx_manifest()
  expect_identical(pgx$versions$libx_release, attr(man, "release"))
  expect_identical(pgx$versions$libx_manifest_sha256, attr(man, "manifest_sha256"))
  expect_true(.pb_LIBX_PATHS[["families"]] %in% pgx$versions$libx_assets$path)
  expect_identical(pgx$versions$reference_source, "libx")
  expect_null(pgx$versions$playdata_version)
})

test_that("a second stamp keeps the release and appends the assets", {
  assets <- function(path) {
    data.frame(path = path, sha256 = "s", source = "libx", stringsAsFactors = FALSE)
  }
  prov <- function(release, path) {
    list(release = release, manifest_sha256 = paste0("m", release),
         reference_source = "libx", assets = assets(path))
  }
  pgx <- .pb_libx_stamp(list(versions = list()), prov("r1", "a"))
  same <- .pb_libx_stamp(pgx, prov("r1", "b"))
  expect_identical(same$versions$libx_release, "r1")
  expect_identical(same$versions$libx_assets$path, c("a", "b"))
  expect_identical(same$versions$reference_source, "libx")
  other <- .pb_libx_stamp(pgx, prov("r2", "b"))
  expect_identical(other$versions$libx_release, "r1")
  expect_identical(other$versions$libx_manifest_sha256, "mr1")
  expect_identical(other$versions$reference_source, "mixed")
})

test_that("phosphosite detection needs no lib-x release", {
  withr::local_options(omicspgx.libx_dir = "")
  expect_false(annotate_phospho_residue(c("P04637_S15", "P38398_T20"), detect.only = TRUE))
  expect_false(annotate_phospho_residue(c("TP53", "BRCA1"), detect.only = TRUE))
  expect_identical(annotate_phospho_residue(c("TP53", "BRCA1")), c("TP53", "BRCA1"))
})

test_that("drug enrichment without a lib-x release errors instead of skipping", {
  withr::local_options(omicspgx.libx_dir = "")
  expect_error(compute_drugActivityEnrichment(list()), "no lib-x release")
  expect_error(compute_drugSensitivityEnrichment(list()), "no lib-x release")
})

test_that("pgx.computeDrugEnrichment keeps its X = NULL default", {
  expect_null(formals(pgx.computeDrugEnrichment)$X)
})

test_that("the default drug matrix is one named asset, missing with a typed error", {
  withr::local_options(omicspgx.libx_dir = "")
  expect_error(
    pgx.computeDrugEnrichment(list(gx.meta = list(), genes = data.frame())),
    class = "pgx_error_libx_missing_asset"
  )
  skip_if_not(nzchar(Sys.getenv("LIBX_DIR")))
  man <- omicspgx::libx_manifest(Sys.getenv("LIBX_DIR"))
  expect_true(.pb_LIBX_PATHS[["l1000_activity"]] %in% man$path)
})

test_that("compute_extra takes sigdb from the lib-x release, never a plain folder", {
  release <- withr::local_tempdir()
  dir.create(file.path(release, "sigdb"))
  file.create(file.path(release, "sigdb", "sigdb-a.h5"))
  writeLines(
    paste("publish", strrep("0", 64), "0", "sigdb/sigdb-a.h5", sep = "\t"),
    file.path(release, "manifest.tsv")
  )
  legacy <- withr::local_tempdir()
  dir.create(file.path(legacy, "sigdb"))
  file.create(file.path(legacy, "sigdb", "sigdb-old.h5"))
  scored <- character()
  testthat::local_mocked_bindings(
    pgx.computeConnectivityScores = function(pgx, db, ...) {
      scored <<- c(scored, basename(db))
      NULL
    },
    .package = "playbase"
  )
  run <- function(libx.dir) {
    suppressMessages(compute_extra(list(),
      extra = "connectivity", pgx.dir = NULL,
      libx.dir = libx.dir, user_input_dir = tempdir()
    ))
  }

  withr::local_options(omicspgx.libx_dir = release)
  run(release)
  run(legacy) # not a release: the configured one is used
  expect_identical(scored, c("sigdb-a.h5", "sigdb-a.h5"))
})

test_that("a NULL lib-x directory resolves to the configured release", {
  skip_if(is.null(omicspgx::libx_dir()))
  expect_identical(.pb_libx_release(NULL), omicspgx::libx_dir())
  expect_error(.pb_libx_release(tempdir()), class = "pgx_error_libx_bad_release")
})

test_that("assets read by pattern are recorded in the provenance", {
  skip_if(is.null(omicspgx::libx_dir()))
  dir <- omicspgx::libx_dir()
  omicspgx::libx_provenance_begin("test")
  assets <- .pb_libx_assets("^refmat/LM22\\.rds$", dir)
  prov <- omicspgx::libx_provenance_end()
  expect_named(assets, "LM22.rds")
  expect_identical(assets[[1]], .pb_libx("lm22"))
  expect_identical(prov$assets$path, "refmat/LM22.rds")
  expect_identical(prov$reference_source, "libx")
})
