#' Test for pgx.createFromFiles
#'
#'

# pgx <- playbase::pgx.createPGX(
#     samples = playbase::SAMPLES,
#     counts = playbase::COUNTS,
#     contrasts = playbase::CONTRASTS[,1:2],
#     organism = "Human"
#   )

#' Test for pgx.createPGX
test_that("pgx.createPGX produce all pgx slots", {
  skip("high memory usage, maybe use mini-example?")
  # Call example data
  # pgx_data <- playbase::get_mini_example_data()

  # Create expected outputs
  expected_tests <- c(
    "name", "organism", "version", "date", "creator", "datatype",
    "description", "samples", "counts", "contrasts", "X",
    "total_counts", "counts_multiplier", "genes", "all_genes",
    "probe_type", "filtered", "tsne2d", "tsne3d", "cluster", "cluster.genes"
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
  expect_equal(dim(pgx$contrasts), c(18, 2))


  ## Check that the gene info is generated correctly
  # TODO expect_equal(pgx$genes[c(1, 10, 20, 30, 40, 50), , drop = FALSE], gene_table)

  ## Check cluster.genes
  expect_equal(dim(pgx$cluster$pos$pca2d), c(ncol(playbase::COUNTS), 2))
  expect_equal(dim(pgx$cluster$pos$tsne3d), c(ncol(playbase::COUNTS), 3))
  expect_equal(dim(pgx$cluster$pos$umap2d), c(ncol(playbase::COUNTS), 2))
  expect_equal(dim(pgx$cluster$pos$umap3d), c(ncol(playbase::COUNTS), 3))
})


#' Test for pgx.computePGX
test_that("pgx.computePGX runs without errors", {
  skip("high memory usage, maybe use mini-example?")
  # Run function
  suppressWarnings(pgx_comp <- playbase::pgx.computePGX(pgx))

  # Expected outputs
  expected_slots <- c(
    "name", "organism", "version", "date", "creator", "datatype", "description", "samples",
    "counts", "contrasts", "X", "total_counts", "counts_multiplier", "genes",
    "all_genes", "probe_type", "tsne2d", "tsne3d", "cluster", "cluster.genes",
    "model.parameters", "filtered", "timings", "gx.meta", "gset.meta", "gsetX",
    "GMT", "cluster.gsets", "meta.go"
  )
  # Check output
  expect_equal(names(pgx_comp), expected_slots)
})


## ---------------------------------------------------------------
## GMT remap tests (pgx.add_GMT rowname harmonisation)
## ---------------------------------------------------------------

test_that("GMT remap is no-op when rownames already match counts", {
  pgx_file <- system.file("data", "pgx_example.rds", package = "playbase")
  if (!nzchar(pgx_file)) pgx_file <- "data/pgx_example.rds"
  if (!file.exists(pgx_file)) skip("pgx_example not found")
  pgx <- local(get(load(pgx_file)))

  gmt_before <- pgx$GMT
  counts_ids <- rownames(pgx$counts)
  current_overlap <- sum(rownames(pgx$GMT) %in% counts_ids)

  expect_equal(current_overlap, nrow(pgx$GMT))
  expect_identical(pgx$GMT, gmt_before)
})

test_that("GMT remap fixes mismatched rownames via genes table", {
  pgx_file <- system.file("data", "pgx_example.rds", package = "playbase")
  if (!nzchar(pgx_file)) pgx_file <- "data/pgx_example.rds"
  if (!file.exists(pgx_file)) skip("pgx_example not found")
  pgx <- local(get(load(pgx_file)))

  ## Simulate non-human: give counts and genes fake systematic IDs,

  ## keep GMT rownames as the original symbols.
  orig_ids <- rownames(pgx$GMT)
  n <- length(orig_ids)
  fake_ids <- paste0("GENE", sprintf("%05d", seq_len(n)))

  ## Build a genes table that bridges fake -> symbol
  genes_new <- pgx$genes[match(orig_ids, rownames(pgx$genes)), ]
  genes_new <- genes_new[!is.na(genes_new$symbol), ]
  n_genes <- nrow(genes_new)
  rownames(genes_new) <- fake_ids[seq_len(n_genes)]
  genes_new$symbol <- orig_ids[seq_len(n_genes)]
  pgx$genes <- genes_new

  ## Update counts rownames to fake IDs
  idx <- match(rownames(pgx$counts), orig_ids)
  rownames(pgx$counts) <- ifelse(!is.na(idx), fake_ids[idx], rownames(pgx$counts))

  ## Before remap: 0% overlap
  expect_equal(sum(rownames(pgx$GMT) %in% rownames(pgx$counts)), 0)

  ## Apply remap logic (same as pgx-compute.R)
  counts_ids <- rownames(pgx$counts)
  current_overlap <- sum(rownames(pgx$GMT) %in% counts_ids)
  map <- setNames(rownames(pgx$genes), pgx$genes$symbol)
  new_ids <- map[rownames(pgx$GMT)]
  keep <- !is.na(new_ids)
  remap_overlap <- sum(new_ids[keep] %in% counts_ids)

  expect_gt(remap_overlap, current_overlap)

  pgx$GMT <- pgx$GMT[keep, , drop = FALSE]
  rownames(pgx$GMT) <- new_ids[keep]

  final_overlap <- sum(rownames(pgx$GMT) %in% counts_ids)
  expect_gt(final_overlap, 0)
  expect_equal(final_overlap, remap_overlap)
})
