# Sample data frame containing all IDs

# devtools::load_all()

ChEBI <- c("15377", "UnrecognizedId", "145677", "15377", "15377", NA)
HMDB <- c("HMDB0015464", "HMDB0250234", "HMDB0000190", NA, "HMDB00004", "UnrecognizedId")
KEGG <- c("C00001", "C00002", "C00003", "UnrecognizedId", "C00004")


# Test cases
test_that("probe_type is correctly validated and converted for metabolomics", {
    # Test valid probe types
    expect_equal(convert_probe_to_chebi(ChEBI, "ChEBI"), "ADD TEST HERE")
    expect_equal(convert_probe_to_chebi(ids = HMDB, probe_type = "HMDB"), "ADD TEST HERE")
    expect_equal(convert_probe_to_chebi(ChEBI, "ChEBI"), "ADD TEST HERE")


    # Test invalid probe type
    expect_error(convert_to_chebi(ids, "InvalidType"), "arg should be one of")
})

# Run the tests
test_file("test_convert_to_chebi.R")
