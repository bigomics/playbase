#' Test for pgx.createFromFiles
#'
#'

#' Test for pgx.createPGX
test_that("pgx.createPGX produce all pgx slots", {
  # Call example data
  # pgx_data <- playbase::get_mini_example_data()

  # Run function
  # Use while to prevent crash on ensembl calls
  suppressWarnings(pgx <- playbase::pgx.createPGX(
    samples = playbase::SAMPLES,
    counts = playbase::COUNTS,
    contrasts = playbase::CONTRASTS,
    organism = "Human"
  ))
  # For every function that uses biomaRt, we need to wait 60 seconds
  # Sys.sleep(60)
  # Create expected outputs
  expected_tests <- c(
    "name", "organism", "version", "date", "creator", "datatype",
    "description", "samples", "counts", "contrasts", "X",
    "total_counts", "counts_multiplier", "genes", "all_genes",
    "probe_type", "tsne2d", "tsne3d", "cluster", "cluster.genes"
  )
  total_counts <- apply(playbase::COUNTS, 2, sum)

  gene_table <- data.frame(
    symbol = c("A1BG", "AGAP2", "ANXA4", "ARPC1A", "BATF", "C19orf53"),
    gene_title = c(
      "alpha-1-B glycoprotein", "ArfGAP with GTPase domain, ankyrin repeat and PH domain 2",
      "annexin A4", "actin related protein 2/3 complex subunit 1A", "basic leucine zipper ATF-like transcription factor",
      "chromosome 19 open reading frame 53"
    ),
    gene_biotype = rep("protein_coding", 6),
    chr = c("19", "12", "2", "7", "14", "19"),
    pos = c(58345178, 57723761, 69644425, 99325898, 75522455, 13774456),
    tx_len = c(2134, 5388, 905, 1582, 617, 897),
    map = c("q13.43", "q14.1", "p13.3", "q22.1", "q24.3", "p13.13"),
    source = c(
      "Source:HGNC Symbol;Acc:HGNC:5", "Source:HGNC Symbol;Acc:HGNC:16921", "Source:HGNC Symbol;Acc:HGNC:542",
      "Source:HGNC Symbol;Acc:HGNC:703", "Source:HGNC Symbol;Acc:HGNC:958", "Source:HGNC Symbol;Acc:HGNC:24991"
    ),
    gene_name = c("A1BG", "AGAP2", "ANXA4", "ARPC1A", "BATF", "C19orf53"),
    feature = c("A1BG", "AGAP2", "ANXA4", "ARPC1A", "BATF", "C19orf53")
  )
  rownames(gene_table) <- gene_table$symbol

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
  expect_equal(pgx$genes[c(1, 10, 20, 30, 40, 50), , drop = FALSE], gene_table)
  expect_equal(
    apply(pgx$genes, 2, class),
    apply(gene_table, 2, class)
  )
shiny::observe({
        shiny::req(pgx$X)
        ct <- colnames(pgx$model.parameters$contr.matrix)
        ct <- sort(ct)
        selected_ct <- head(ct, 7)
        shiny::updateSelectInput(
          session,
          "selected_contrasts",
          choices = ct,
          selected = selected_ct)
      })
  ## Check cluster.genes
  expect_equal(dim(pgx$cluster$pos$pca2d), c(ncol(playbase::COUNTS), 2))
  expect_equal(dim(pgx$cluster$pos$tsne3d), c(ncol(playbase::COUNTS), 3))
  expect_equal(dim(pgx$cluster$pos$umap2d), c(ncol(playbase::COUNTS), 2))
  expect_equal(dim(pgx$cluster$pos$umap3d), c(ncol(playbase::COUNTS), 3))
})


#' Test for pgx.computePGX
test_that("pgx.computePGX runs without errors", {
  # Run function
  suppressWarnings(pgx_comp <- playbase::pgx.computePGX(playbase::PGX_CREATE))

  # Expected outputs
  expected_slots <- c(
    "name", "organism", "version", "date", "datatype", "description", "samples",
    "counts", "contrasts", "X", "total_counts", "counts_multiplier", "genes",
    "all_genes", "probe_type", "tsne2d", "tsne3d", "cluster", "cluster.genes",
    "model.parameters", "filtered", "timings", "gx.meta", "gset.meta", "gsetX",
    "GMT", "cluster.gsets", "meta.go"
  )
  # Check output
  expect_equal(names(pgx_comp), expected_slots)
})
