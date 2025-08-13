# example of metabolite ids for testing
ChEBI <- c("8874", "UnrecognizedId", "1445644", "49603", "15903", "", NA)
HMDB <- c("HMDB0015464", "HMDB0250234", "HMDB0000190", NA, "HMDB0000466", "UnrecognizedId", "HMDB999999")
KEGG <- c("C06958", "C06849", "C18044", "UnrecognizedId", "C16545")

# Test cases
test_that("probe_type is correctly validated and converted for metabolomics", {
  # Test valid probe types
  expect_equal(playbase::mx.convert_probe(probes = ChEBI, probe_type = "ChEBI"), c("8874", NA, "1445644", "49603", "15903", NA, NA))
  expect_equal(playbase::mx.convert_probe(probes = HMDB, probe_type = "HMDB"), c("10093", NA, "16651", NA, "9171", NA, NA))
  expect_equal(playbase::mx.convert_probe(probes = KEGG, "KEGG"), c("101278", "3082", "35420", NA, "63825"))

  # chebi handles duplicated values
  expect_equal(playbase::mx.convert_probe(c("8874", "8874"), "ChEBI"), c("8874", "8874"))

  # HMDB handles duplicated values
  expect_equal(playbase::mx.convert_probe(c("HMDB0015464", "HMDB0015464"), "HMDB"), c("10093", "10093"))

  # KEGG handles duplicated values
  expect_equal(playbase::mx.convert_probe(c("C06958", "C06958"), "KEGG"), c("101278", "101278"))

  # Test invalid probe type returns error
  expect_null(playbase::mx.convert_probe(ChEBI, "InvalidType"))

  # Test empty string is converted to NA
  expect_equal(playbase::mx.convert_probe(c("", "HMDB0015464"), "HMDB"), c(NA, "10093"))
})
