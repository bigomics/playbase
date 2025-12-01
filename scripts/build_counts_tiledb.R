#!/usr/bin/env Rscript
#' Build TileDB database from PGX files
#'
#' This script scans a folder for .pgx files, extracts their counts matrices,
#' and stores them in a TileDB database for efficient querying.
#'
#' Usage:
#'   Rscript build_counts_tiledb.R <pgx_folder> <tiledb_output_path> [--overwrite]
#'
#' Example:
#'   Rscript build_counts_tiledb.R /path/to/pgx/files /path/to/output/counts_db
#'   Rscript build_counts_tiledb.R /path/to/pgx/files /path/to/output/counts_db --overwrite

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  cat("Usage: Rscript build_counts_tiledb.R <pgx_folder> <tiledb_output_path> [--overwrite]\n")
  cat("\n")
  cat("Arguments:\n")
  cat("  pgx_folder        Path to folder containing .pgx files\n")
  cat("  tiledb_output_path Path where TileDB database will be created\n")
  cat("  --overwrite       Optional: overwrite existing database\n")
  cat("\n")
  cat("Example:\n")
  cat("  Rscript build_counts_tiledb.R /home/bigomics/Playground/omicsplayground/data ./counts_tiledb\n")
  quit(status = 1)
}

pgx_folder <- args[1]
tiledb_path <- args[2]
overwrite <- "--overwrite" %in% args

# Check dependencies
cat("Checking dependencies...\n")

if (!requireNamespace("tiledb", quietly = TRUE)) {
  stop("Package 'tiledb' is required. Install with: install.packages('tiledb')")
}

if (!requireNamespace("playbase", quietly = TRUE)) {
  stop("Package 'playbase' is required. Install from GitHub.")
}

library(playbase)
library(tiledb)

cat("Dependencies OK\n\n")

# Validate input folder
if (!dir.exists(pgx_folder)) {
  stop("PGX folder does not exist: ", pgx_folder)
}

# Count .pgx files
pgx_files <- list.files(pgx_folder, pattern = "\\.pgx$", full.names = TRUE)
cat("Found", length(pgx_files), ".pgx files in:", pgx_folder, "\n\n")

if (length(pgx_files) == 0) {
  stop("No .pgx files found in the specified folder")
}

# Build the database
cat("Building TileDB database...\n")
cat("Output path:", tiledb_path, "\n")
cat("Overwrite:", overwrite, "\n\n")

start_time <- Sys.time()

metadata <- playbase::pgx.buildTileDB(
  pgx_folder = pgx_folder,
  tiledb_path = tiledb_path,
  overwrite = overwrite,
  verbose = TRUE
)

end_time <- Sys.time()
elapsed <- difftime(end_time, start_time, units = "mins")

cat("\n")
cat("========================================\n")
cat("TileDB database built successfully!\n")
cat("========================================\n")
cat("Time elapsed:", round(as.numeric(elapsed), 2), "minutes\n")
cat("Total genes:", metadata$n_genes, "\n")
cat("Total samples:", metadata$n_samples, "\n")
cat("Total files processed:", metadata$n_files, "\n")
cat("\nDatabase location:", tiledb_path, "\n")
cat("Metadata file:", paste0(tiledb_path, "_metadata.rds"), "\n")

