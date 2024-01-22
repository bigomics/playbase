# Generate shared objects
counter <- 0
while (!exists("ensembl_human")) {
  Sys.sleep(60 * counter)
  ensembl <- biomaRt::useEnsembl(biomart = "genes", version = 110)
  ensembl_human <- biomaRt::useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)
  counter <- counter + 1
}

#' Test for guess_probetype
test_that("guess_probetype can detect ensembl IDs", {
  # Create input data with <- for reuse
  probes <- c(
    "ENSG00000230915.1", "ENSG00000275728.1", "ENSG00000277599.1",
    "ENSG00000186163.9", "ENSG00000164823.11", "ENSG00000274234.1",
    "ENSG00000282461.1", "ENSG00000283056.1", "ENSG00000239021.1",
    "ENSG00000214268.2", "ENSG00000206687.1", "ENSG00000171148.14",
    "ENSG00000250027.1", "ENSG00000244217.1", "ENSG00000103502.14",
    "ENSG00000213178.3", "ENSG00000235059.5", "ENSG00000204555.3",
    "ENSG00000221044.2", "ENSG00000267162.1"
  )

  # Run function
  # Use while to prevent crash on ensembl calls
  counter <- 0
  while (!exists("type")) {
    type <- playbase::guess_probetype()(probes, mart = ensembl_human)
    counter <- counter + 1
  }

  # Check output
  expect_equal(type, "ensembl_gene_id_version")
})

test_that("guess_probetype stops without mart", {
  expect_error(playbase::guess_probetype(probes))
})


#' Test for ngs.getGeneAnnotation
test_that("ngs.getGeneAnnotation returns annotation for genes", {
  # Input data
  probes <- c(
    "ENSG00000230915.1", "ENSG00000275728.1", "ENSG00000277599.1",
    "ENSG00000186163.9", "ENSG00000164823.11", "ENSG00000274234.1",
    "ENSG00000282461.1", "ENSG00000283056.1", "ENSG00000239021.1",
    "ENSG00000214268.2", "ENSG00000206687.1", "ENSG00000171148.14",
    "ENSG00000250027.1", "ENSG00000244217.1", "ENSG00000103502.14",
    "ENSG00000213178.3", "ENSG00000235059.5", "ENSG00000204555.3",
    "ENSG00000221044.2", "ENSG00000267162.1"
  )

  # Run function <- for reuse
  counter <- 0
  while (!exists("result")) {
    Sys.sleep(60 * counter)
    result <- playbase::ngs.getGeneAnnotation(
      probes = probes,
      organism = "Human",
      probe_type = "ensembl_gene_id_version",
      mart = ensembl_human
    )
    counter <- counter + 1
  }
  # Check class
  expect_s3_class(result, "data.frame")

  # Check number of rows
  expect_equal(nrow(result), length(probes))

  # Check presence of expected columns
  expected_cols <- c("gene_name", "gene_title", "gene_biotype", "chr", "pos", "tx_len", "map")
  expect_true(all(expected_cols %in% colnames(result)))

  # Check gene names match
  expect_equal(result$gene_name, probes)
})


#' Test for probe2symbol
test_that("probe2symbol returns expected output", {
  # Input data
  probes <- c(
    "ENSG00000230915.1", "ENSG00000275728.1", "ENSG00000277599.1",
    "ENSG00000186163.9", "ENSG00000164823.11", "ENSG00000274234.1",
    "ENSG00000282461.1", "ENSG00000283056.1", "ENSG00000239021.1",
    "ENSG00000214268.2", "ENSG00000206687.1", "ENSG00000171148.14",
    "ENSG00000250027.1", "ENSG00000244217.1", "ENSG00000103502.14",
    "ENSG00000213178.3", "ENSG00000235059.5", "ENSG00000204555.3",
    "ENSG00000221044.2", "ENSG00000267162.1"
  )

  # Run function
  counter <- 0
  while (!exists("result")) {
    Sys.sleep(60 * counter)
    result <- playbase::ngs.getGeneAnnotation(
      probes = probes,
      organism = "Human",
      probe_type = "ensembl_gene_id_version",
      mart = ensembl_human
    )
    counter <- counter + 1
  }
  # Run function
  symbol <- playbase::probe2symbol(probes, result, query = "symbol", fill_na = FALSE)

  # Default parameters
  expect_equal(length(symbol), length(probes))
  expect_type(symbol, "character")

  # Test handling NAs with fill_na = FALSE
  symbol_na <- playbase::probe2symbol(probes, result, query = "symbol", fill_na = FALSE)
  expect_type(symbol_na, "character")
  expect_true(sum(symbol_na == "") == 4)

  # Test handling NAs with fill_na = TRUE
  symbol_na <- playbase::probe2symbol(probes, result, query = "symbol", fill_na = TRUE)
  expect_type(symbol_na, "character")
  expect_true(sum(symbol_na == "") == 0)
})


#' Test for guess_probetype
test_that("detects ENSEMBL for human probes", {
  probes <- c("ENSG00000136997", "ENSG00000241860")
  expect_equal(guess_probetype(probes, "Human"), "ENSEMBL")

  # UNIPROT genes
  uniprot_genes <- c("P31749", "P04637", "Q9Y6K9", "O15111", "Q9UM73")
  expect_equal(guess_probetype(uniprot_genes, "Human"), "UNIPROT")

  # Fake genes
  probes <- c("ENSG00088136997", "ENSG00099241860")
  expect_error(guess_probetype(probes, "Human"))
})

# Test with valid mouse probes
test_that("detects Ensembl for mouse probes", {
  probes <- c("ENSMUSG00000051951", "ENSMUSG00000033845")
  expect_equal(guess_probetype(probes, "Mouse"), "ENSEMBL")

  probes <- c(
    "NM_001081979", "NM_001081980", "NM_001081981", "NM_001081982",
    "NM_001081983"
  )
  expect_equal(guess_probetype(probes, "Mouse"), "REFSEQ")
})


#' Test for guess_probetype
test_that("ngs.getGeneAnnotation_ORGDB function works correctly", {
  # Test 1: Check that the function returns the correct annotation for a known human gene
  expect_equal(rownames(ngs.getGeneAnnotation_ORGDB("ENSG00000141510", "ENSEMBL", "Human"))[1], "ENSG00000141510")

  # Test 2: Check that the function returns the correct annotation for a known mouse gene
  expect_equal(rownames(ngs.getGeneAnnotation_ORGDB("ENSMUSG00000051951", "ENSEMBL", "Mouse"))[1], "ENSMUSG00000051951")

  # Test 3: Check that the function handles multiple probes correctly
  probes <- c("ENSG00000141510", "ENSG00000139618")
  expect_equal(nrow(ngs.getGeneAnnotation_ORGDB(probes, "ENSEMBL", "Human")), length(probes))

  # Test 4: Check that the function handles an unknown organism correctly
  expect_error(ngs.getGeneAnnotation_ORGDB("ENSG00000141510", "ENSEMBL", "Unknown"))

  # Test 5: Check that the function handles an unknown probe correctly
  expect_error(ngs.getGeneAnnotation_ORGDB("Unknown", "ENSEMBL", "Human"))

  # Test 6: Check that the function handles a NULL probe correctly
  expect_error(ngs.getGeneAnnotation_ORGDB(NULL, "ENSEMBL", "Human"))

  # Test 7: Check that the function handles a NULL organism correctly
  expect_error(ngs.getGeneAnnotation_ORGDB("ENSG00000141510", "ENSEMBL", NULL))

  # Test 8: Check that the function handles an empty string probe correctly
  expect_error(ngs.getGeneAnnotation_ORGDB("", "ENSEMBL", "Human"))

  # Test 9: Check that the function handles an empty string organism correctly
  expect_error(ngs.getGeneAnnotation_ORGDB("ENSG00000141510", "ENSEMBL", ""))
})
