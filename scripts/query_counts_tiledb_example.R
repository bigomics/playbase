#!/usr/bin/env Rscript
#' Example: Query gene counts from TileDB database
#'
#' This script demonstrates how to query the TileDB database
#' created by build_counts_tiledb.R

library(playbase)

# Path to your TileDB database (adjust as needed)
tiledb_path <- "./counts_tiledb"

# Check if database exists
if (!dir.exists(tiledb_path)) {
  stop("TileDB database not found at: ", tiledb_path,
       "\nPlease run build_counts_tiledb.R first to create the database.")
}

# Get database info
cat("=== Database Info ===\n")
playbase::pgx.infoTileDB(tiledb_path)

# List available genes (first 20)
cat("\n=== Available genes (first 20) ===\n")
genes <- playbase::pgx.listGenesTileDB(tiledb_path)
cat(paste(head(genes, 20), collapse = ", "), "...\n")
cat("Total genes:", length(genes), "\n")

# List available samples (first 10)
cat("\n=== Available samples (first 10) ===\n")
samples <- playbase::pgx.listSamplesTileDB(tiledb_path)
cat(paste(head(samples, 10), collapse = "\n"), "\n")
cat("Total samples:", length(samples), "\n")

# Query a single gene
cat("\n=== Query single gene (TP53) ===\n")
if ("TP53" %in% genes) {
  counts_tp53 <- playbase::pgx.queryTileDB(tiledb_path, genes = "TP53")
  cat("Dimensions:", dim(counts_tp53), "\n")
  cat("First 5 samples:\n")
  print(head(counts_tp53[1, counts_tp53[1, ] > 0], 5))
} else {
  cat("TP53 not found in database. Trying first available gene...\n")
  counts_first <- playbase::pgx.queryTileDB(tiledb_path, genes = genes[1])
  cat("Gene:", genes[1], "\n")
  cat("Dimensions:", dim(counts_first), "\n")
}

# Query multiple genes
cat("\n=== Query multiple genes ===\n")
query_genes <- head(genes, 5)
counts_multi <- playbase::pgx.queryTileDB(tiledb_path, genes = query_genes)
cat("Genes queried:", paste(query_genes, collapse = ", "), "\n")
cat("Matrix dimensions:", dim(counts_multi), "\n")

# Get as data.frame (long format)
cat("\n=== Query as data.frame (first 10 rows) ===\n")
df <- playbase::pgx.queryTileDB(tiledb_path, genes = query_genes[1], as_matrix = FALSE)
print(head(df, 10))

