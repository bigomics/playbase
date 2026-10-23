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

test_that("mx.lipid_characteristics parses class, sum composition and chains", {
  lc <- playbase::mx.lipid_characteristics(
    c("PC 16:0_18:1", "PE O-38:4", "Cer 32:0;2", "CE 14:0", "not a lipid", "-")
  )
  expect_equal(nrow(lc), 6)
  expect_equal(lc$lipid_class[1:4], c("PC", "PE O-", "Cer", "CE"))
  expect_equal(lc$lipid_category[1], "Glycerophospholipids")
  expect_equal(lc$total_c[1:2], c(34L, 38L))
  expect_equal(lc$total_db[1:2], c(1L, 4L))
  expect_equal(lc$fa_chains[1], "16:0, 18:1")
  expect_true(is.na(lc$fa_chains[2])) ## species level: chains unknown
  expect_true(all(is.na(lc[5:6, "lipid_class"])))
})

test_that("mx.create_lipid_sets builds the five families from parsed lipids", {
  annot <- data.frame(symbol = c("PC 16:0_18:1", "PC 18:0_20:4", "PE 38:4", "x"))
  annot <- cbind(annot, playbase::mx.lipid_characteristics(annot$symbol))
  sets <- playbase::mx.create_lipid_sets(annot)
  expect_setequal(sets[["LIPID_CLASS:PC"]], c("PC 16:0_18:1", "PC 18:0_20:4"))
  expect_setequal(sets[["LIPID_UNSAT:4+ double bonds"]], c("PC 18:0_20:4", "PE 38:4"))
  expect_setequal(sets[["LIPID_CHAIN:38 carbons"]], c("PC 18:0_20:4", "PE 38:4"))
  expect_equal(sets[["LIPID_FA:FA 20:4"]], "PC 18:0_20:4")
  expect_false(any(grepl("^LIPID_FA:", names(sets)) & sapply(sets, function(s) "PE 38:4" %in% s)))
  ## nothing parsed: no lipid sets, the class-set fallback takes over
  expect_length(playbase::mx.create_lipid_sets(annot[4, , drop = FALSE]), 0)
})
