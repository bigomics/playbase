## playbase keeps these entry points as forwarders to playbase.annotation:
## same signature as the package function, same result.
forwarders <- c(
  getProbeAnnotation = "getProbeAnnotation",
  getHumanOrtholog = "getHumanOrtholog",
  detect_probetype = "detect_probetype",
  check_species_probetype = "check_species_probetype",
  pgx.getFeatureInfo = "pgx.getFeatureInfo",
  pgx.add_GMT = "pgx.addGMT"
)

for (fn in names(forwarders)) {
  test_that(paste(fn, "has the formals of its playbase.annotation target"), {
    target <- getExportedValue("playbase.annotation", forwarders[[fn]])
    expect_identical(formals(get(fn, envir = asNamespace("playbase"))), formals(target))
  })
}

test_that("detect_probetype gives the playbase.annotation result", {
  skip_if_not_installed("org.Hs.eg.db")
  probes <- c("TP53", "BRCA1", "EGFR", "MYC", "KRAS")
  expect_identical(
    detect_probetype("Human", probes),
    playbase.annotation::detect_probetype("Human", probes)
  )
})

test_that("mx.convert_probe detects the probe type through playbase.annotation", {
  probes <- c("8874", "15903")
  expect_identical(mx.convert_probe(probes), mx.convert_probe(probes, "ChEBI"))
})
