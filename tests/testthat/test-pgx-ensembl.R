# Generate shared objects
counter <- 0
while (!exists("ensembl_human")) {
  Sys.sleep(60 * counter)
  ensembl <- biomaRt::useEnsembl(biomart="genes", version = 110)
  ensembl_human <- biomaRt::useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)
  counter <- counter + 1
}

#' Test for detect_probe
test_that("detect_probe can detect ensembl IDs", {
  # Create input data with <- for reuse
  probes <- c("ENSG00000230915.1", "ENSG00000275728.1", "ENSG00000277599.1",
    "ENSG00000186163.9", "ENSG00000164823.11", "ENSG00000274234.1",
    "ENSG00000282461.1", "ENSG00000283056.1", "ENSG00000239021.1",
    "ENSG00000214268.2", "ENSG00000206687.1", "ENSG00000171148.14",
    "ENSG00000250027.1", "ENSG00000244217.1", "ENSG00000103502.14",
    "ENSG00000213178.3", "ENSG00000235059.5", "ENSG00000204555.3",
    "ENSG00000221044.2", "ENSG00000267162.1")

  # Run function
  # Use while to prevent crash on ensembl calls
  counter <- 0
  while (!exists("type")) {  
      type <- playbase::detect_probe(probes, mart = ensembl_human)
      counter <- counter + 1
  }

  # Check output
  expect_equal(type, "ensembl_gene_id_version")
  
  })

test_that("detect_probe stops without mart", {
  expect_error(playbase::detect_probe(probes))
})


#' Test for ngs.getGeneAnnotation
test_that("ngs.getGeneAnnotation returns annotation for genes", {

  # Input data
  probes <- c("ENSG00000230915.1", "ENSG00000275728.1", "ENSG00000277599.1",
    "ENSG00000186163.9", "ENSG00000164823.11", "ENSG00000274234.1",
    "ENSG00000282461.1", "ENSG00000283056.1", "ENSG00000239021.1",
    "ENSG00000214268.2", "ENSG00000206687.1", "ENSG00000171148.14",
    "ENSG00000250027.1", "ENSG00000244217.1", "ENSG00000103502.14",
    "ENSG00000213178.3", "ENSG00000235059.5", "ENSG00000204555.3",
    "ENSG00000221044.2", "ENSG00000267162.1")
  
  # Run function <- for reuse
  counter <- 0
  while (!exists("result")) {
    Sys.sleep(60 * counter)
    result <- playbase::ngs.getGeneAnnotation(probes = probes, 
                                              organism = "Human", 
                                              probe_type = "ensembl_gene_id_version", 
                                              mart = ensembl_human)
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
  probes <- c("ENSG00000230915.1", "ENSG00000275728.1", "ENSG00000277599.1",
    "ENSG00000186163.9", "ENSG00000164823.11", "ENSG00000274234.1",
    "ENSG00000282461.1", "ENSG00000283056.1", "ENSG00000239021.1",
    "ENSG00000214268.2", "ENSG00000206687.1", "ENSG00000171148.14",
    "ENSG00000250027.1", "ENSG00000244217.1", "ENSG00000103502.14",
    "ENSG00000213178.3", "ENSG00000235059.5", "ENSG00000204555.3",
    "ENSG00000221044.2", "ENSG00000267162.1")

  # Run function
  counter <- 0
  while (!exists("result")) {
    Sys.sleep(60 * counter)
    result <- playbase::ngs.getGeneAnnotation(probes = probes, 
                                              organism = "Human", 
                                              probe_type = "ensembl_gene_id_version", 
                                              mart = ensembl_human)
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