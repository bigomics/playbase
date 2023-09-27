#' Test for pgx.createFromFiles
#'
#'


#' Test for pgx.createPGX
test_that("pgx.createPGX runs without errors", {
  # Create mock input data
  pgx_data <- playbase:::get_mini_example_data()

  # Run function
  pgx <- playbase::pgx.createPGX(
    samples = pgx_data$samples,
    counts = pgx_data$counts,
    contrasts = pgx_data$contrast
  )

  # Create expected outputs
  expected_tests <- c("name", "date", "datatype", "description", "samples", "counts", "contrasts",
                      "X", "total_counts", "counts_multiplier", "genes", "all_genes", "symbol",
                      "tsne2d", "tsne3d", "cluster", "cluster.genes")
  total_counts <- apply(pgx_data$counts, 2, sum)
  new_genes <- data.table::data.table(
    hgnc_symbol = c("A1BG", "ACOT9", "AGAP3", "AGAP3", "ALG12", "ARF5"),
    external_gene_name = c("A1BG", "ACOT9", "AGAP3", "AGAP3", "ALG12", "ARF5"),
    description = c("alpha-1-B glycoprotein [Source:HGNC Symbol;Acc:HGNC:5]",
                    "acyl-CoA thioesterase 9 [Source:HGNC Symbol;Acc:HGNC:17152]",
                    "ArfGAP with GTPase domain, ankyrin repeat and PH domain 3 [Source:HGNC Symbol;Acc:HGNC:16923]",
                    "ArfGAP with GTPase domain, ankyrin repeat and PH domain 3 [Source:HGNC Symbol;Acc:HGNC:16923]",
                    "ALG12 alpha-1,6-mannosyltransferase [Source:HGNC Symbol;Acc:HGNC:19358]",
                    "ADP ribosylation factor 5 [Source:HGNC Symbol;Acc:HGNC:658]"),
    gene_biotype = rep("protein_coding", 6),
    chromosome_name = c(19, "X", 7, 7, 22, 7),
    transcript_start = c(58345178, 23703705, 151089712, 151120016, 49903650, 127588411),
    transcript_length = c(2134, 1682, 1019, 413, 674, 1032),
    band = c("q13.43", "p22.11", "q36.1", "q36.1", "q13.33", "q32.1")
  )
  # genes <- data.frame(
  #   gene_name = c("A1BG", "ACOT9", "ALG12", "ANXA7", "ARF5", "ARPC4"),
  #   gene_title = c(
  #     "alpha-1-B glycoprotein", "acyl-CoA thioesterase 9",
  #     "ALG12 alpha-1,6-mannosyltransferase", "annexin A7",
  #     "ADP ribosylation factor 5", "actin related protein 2/3 complex subunit 4"
  #   ),
  #   gene_biotype = rep("protein_coding", 6),
  #   chr = c(19, "X", 22, 10, 7, 3),
  #   pos = c("58345182", "23701055", "49900228", "73375100", "127588410", "9793081"),
  #   tx_len = c("2592", "1427", "2373", "1734", "1022", "1624"),
  #   map = c("19q13.43", "Xp22.11", "22q13.33", "10q22.2", "7q32.1", "3p25.3")
  # )

  # Check output
  ## Check all te test present
  expect_true(all(names(pgx) == expected_tests))

  ## Check the total counts in each sample
  expect_equal(pgx$total_counts, total_counts)

  ## Check multiplier
  expect_equal(pgx$counts_multiplier, 1)

  ## Check contrast
  expect_equal(dim(pgx$contrasts), c(7, 3))


  ## Check that the gene info is generated correctly
  expect_equal(pgx$genes[c(1, 10, 20, 30, 40, 50)], new_genes)
  expect_equal(
    apply(pgx$genes, 2, class),
    apply(new_genes, 2, class)
  )

  ## Check cluster.genes
  expect_equal(dim(pgx$cluster$pos$pca2d), c(ncol(pgx_data$counts), 2))
  expect_equal(dim(pgx$cluster$pos$tsne3d), c(ncol(pgx_data$counts), 3))
  expect_equal(dim(pgx$cluster$pos$umap2d), c(ncol(pgx_data$counts), 2))
  expect_equal(dim(pgx$cluster$pos$umap3d), c(ncol(pgx_data$counts), 3))
})


#' Test for pgx.computePGX
test_that("pgx.computePGX runs without errors", {
  # Create mock input data
  pgx <- playbase::pgx.createPGX(
    samples = playbase::SAMPLES,
    counts = playbase::COUNTS,
    contrasts = playbase::CONTRASTS
  )

  # Run function
  pgx_comp <- playbase::pgx.computePGX(pgx)

  # Create expected outputs
  expected_tests <- c("name", "date", "datatype", "description", "samples", "counts", "contrasts",
                      "X", "total_counts", "counts_multiplier", "genes", "all_genes", "symbol",
                      "tsne2d", "tsne3d", "cluster", "cluster.genes", "model.parameters",
                      "filtered", "timings", "gx.meta", "gset.meta", "gsetX", "GMT",
                      "cluster.gsets", "meta.go")
  # Check output
  expect_true(all(names(pgx_comp) == expected_tests))
})
