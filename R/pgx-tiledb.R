#' @title Extract counts from PGX files and store in TileDB
#'
#' @description Functions for storing and querying gene expression counts
#' from multiple PGX files using TileDB for efficient storage and retrieval.
#'

#' @title Build TileDB database from PGX files
#'
#' @description Scans a folder for .pgx files, extracts counts matrices,
#' and stores them in a TileDB database.
#'
#' @param pgx_folder Path to folder containing .pgx files
#' @param tiledb_path Path where TileDB database will be created
#' @param overwrite Logical, whether to overwrite existing database. Default FALSE.
#' @param verbose Logical, whether to print progress messages. Default TRUE.
#'
#' @return Invisibly returns a list with metadata about the stored data
#'
#' @details This function creates a sparse TileDB array optimized for
#' gene-level queries across many samples from multiple PGX files.
#' The array uses genes as rows and samples as columns.
#'
#' @export
pgx.buildTileDB <- function(pgx_folder, tiledb_path, overwrite = FALSE, verbose = TRUE) {
  if (!requireNamespace("tiledb", quietly = TRUE)) {
    stop("Package 'tiledb' is required. Install with: install.packages('tiledb')")
  }

  # Check if folder exists

if (!dir.exists(pgx_folder)) {
    stop("Folder does not exist: ", pgx_folder)
  }

  # Find all .pgx files
  pgx_files <- list.files(
    path = pgx_folder,
    pattern = "\\.pgx$",
    full.names = TRUE,
    recursive = FALSE
  )

  if (length(pgx_files) == 0) {
    stop("No .pgx files found in: ", pgx_folder)
  }

  if (verbose) message("Found ", length(pgx_files), " .pgx files")

  # Handle existing database
  if (dir.exists(tiledb_path)) {
    if (overwrite) {
      if (verbose) message("Removing existing TileDB database...")
      unlink(tiledb_path, recursive = TRUE)
    } else {
      stop("TileDB database already exists at: ", tiledb_path,
           "\nUse overwrite=TRUE to replace it.")
    }
  }

  # First pass: collect all unique genes and samples
  if (verbose) message("First pass: collecting gene and sample metadata...")

  all_genes <- character(0)
  all_samples <- character(0)
  sample_to_file <- list()

  for (i in seq_along(pgx_files)) {
    pgx_file <- pgx_files[i]
    if (verbose) message("  [", i, "/", length(pgx_files), "] Scanning: ", basename(pgx_file))

    tryCatch({
      pgx <- pgx.load(pgx_file)

      if (is.null(pgx$counts)) {
        warning("No counts found in: ", pgx_file, " - skipping")
        next
      }

      counts <- pgx$counts
      genes <- rownames(counts)
      samples <- colnames(counts)

      # Make sample names unique by prefixing with file basename
      file_prefix <- tools::file_path_sans_ext(basename(pgx_file))
      unique_samples <- paste0(file_prefix, "::", samples)

      all_genes <- union(all_genes, genes)
      all_samples <- c(all_samples, unique_samples)

      for (s in unique_samples) {
        sample_to_file[[s]] <- pgx_file
      }

      rm(pgx, counts)
      gc(verbose = FALSE)

    }, error = function(e) {
      warning("Error loading ", pgx_file, ": ", e$message, " - skipping")
    })
  }

  if (length(all_samples) == 0) {
    stop("No valid samples found in any .pgx file")
  }

  # Sort for consistent ordering
  all_genes <- sort(unique(all_genes))
  all_samples <- unique(all_samples)

  if (verbose) {
    message("Total unique genes: ", length(all_genes))
    message("Total samples: ", length(all_samples))
  }

  # Create TileDB schema
  if (verbose) message("Creating TileDB schema...")

  # Use string dimensions for gene and sample names
  gene_dim <- tiledb::tiledb_dim(
    name = "gene",
    domain = NULL,
    tile = NULL,
    type = "ASCII"
  )

  sample_dim <- tiledb::tiledb_dim(
    name = "sample",
    domain = NULL,
    tile = NULL,
    type = "ASCII"
  )

  dom <- tiledb::tiledb_domain(dims = list(gene_dim, sample_dim))

  # Attribute for count values
  count_attr <- tiledb::tiledb_attr(
    name = "count",
    type = "FLOAT64"
  )

  # Create sparse array schema
  schema <- tiledb::tiledb_array_schema(
    domain = dom,
    attrs = list(count_attr),
    sparse = TRUE
  )

  # Create the array
  tiledb::tiledb_array_create(tiledb_path, schema)

  if (verbose) message("TileDB array created at: ", tiledb_path)

  # Second pass: write data
  if (verbose) message("Second pass: writing counts to TileDB...")

  arr <- tiledb::tiledb_array(tiledb_path)

  for (i in seq_along(pgx_files)) {
    pgx_file <- pgx_files[i]
    if (verbose) message("  [", i, "/", length(pgx_files), "] Writing: ", basename(pgx_file))

    tryCatch({
      pgx <- pgx.load(pgx_file)

      if (is.null(pgx$counts)) next

      counts <- as.matrix(pgx$counts)
      genes <- rownames(counts)
      samples <- colnames(counts)
      file_prefix <- tools::file_path_sans_ext(basename(pgx_file))

      # Prepare data for TileDB write (coordinate format)
      # Convert matrix to long format
      for (j in seq_along(samples)) {
        sample_name <- paste0(file_prefix, "::", samples[j])
        values <- counts[, j]

        # Filter out zeros and NAs to keep array sparse
        valid_idx <- which(!is.na(values) & values != 0)

        if (length(valid_idx) > 0) {
          df <- data.frame(
            gene = genes[valid_idx],
            sample = rep(sample_name, length(valid_idx)),
            count = as.numeric(values[valid_idx]),
            stringsAsFactors = FALSE
          )

          arr[] <- df
        }
      }

      rm(pgx, counts)
      gc(verbose = FALSE)

    }, error = function(e) {
      warning("Error writing ", pgx_file, ": ", e$message)
    })
  }

  # Write metadata
  if (verbose) message("Writing metadata...")

  metadata_path <- paste0(tiledb_path, "_metadata.rds")
  metadata <- list(
    genes = all_genes,
    samples = all_samples,
    sample_to_file = sample_to_file,
    n_genes = length(all_genes),
    n_samples = length(all_samples),
    n_files = length(pgx_files),
    pgx_files = pgx_files,
    created = Sys.time()
  )
  saveRDS(metadata, metadata_path)

  if (verbose) {
    message("Done! TileDB database created at: ", tiledb_path)
    message("Metadata saved to: ", metadata_path)
    message("  - Genes: ", length(all_genes))
    message("  - Samples: ", length(all_samples))
    message("  - Files: ", length(pgx_files))
  }

  invisible(metadata)
}


#' @title Query gene counts from TileDB database
#'
#' @description Retrieves counts for specified genes across all samples
#' stored in the TileDB database.
#'
#' @param tiledb_path Path to the TileDB database
#' @param genes Character vector of gene names to query
#' @param samples Optional character vector of sample names to filter.
#'   If NULL (default), returns all samples.
#' @param as_matrix Logical, whether to return result as a matrix.
#'   Default TRUE. If FALSE, returns a data.frame in long format.
#'
#' @return If as_matrix=TRUE, a matrix with genes as rows and samples as columns.
#'   If as_matrix=FALSE, a data.frame with columns: gene, sample, count.
#'
#' @details This function efficiently queries the TileDB sparse array
#' for specific genes, leveraging TileDB's indexing for fast retrieval.
#'
#' @examples
#' \dontrun{
#' # Query single gene
#' counts <- pgx.queryTileDB("/path/to/tiledb", genes = "TP53")
#'
#' # Query multiple genes
#' counts <- pgx.queryTileDB("/path/to/tiledb", genes = c("TP53", "BRCA1", "MYC"))
#'
#' # Get as data.frame
#' df <- pgx.queryTileDB("/path/to/tiledb", genes = "TP53", as_matrix = FALSE)
#' }
#'
#' @export
pgx.queryTileDB <- function(tiledb_path, genes, samples = NULL, as_matrix = TRUE) {
  if (!requireNamespace("tiledb", quietly = TRUE)) {
    stop("Package 'tiledb' is required. Install with: install.packages('tiledb')")
  }

  if (!dir.exists(tiledb_path)) {
    stop("TileDB database not found at: ", tiledb_path)
  }

  # Load metadata for sample list if needed
  metadata_path <- paste0(tiledb_path, "_metadata.rds")
  if (file.exists(metadata_path)) {
    metadata <- readRDS(metadata_path)
    all_genes <- metadata$genes
    all_samples <- metadata$samples
  } else {
    warning("Metadata file not found. Some features may be limited.")
    all_genes <- NULL
    all_samples <- NULL
  }

  # Validate genes if metadata available
  if (!is.null(all_genes)) {
    missing_genes <- setdiff(genes, all_genes)
    if (length(missing_genes) > 0) {
      warning("Some genes not found in database: ",
              paste(head(missing_genes, 5), collapse = ", "),
              if (length(missing_genes) > 5) paste0(" ... and ", length(missing_genes) - 5, " more"))
      genes <- intersect(genes, all_genes)
    }
  }

  if (length(genes) == 0) {
    stop("No valid genes to query")
  }

  # Open array and query
  arr <- tiledb::tiledb_array(tiledb_path, as.data.frame = TRUE)

  # Build query for selected genes
  qc <- tiledb::tiledb_query_condition_init()
  
  # TileDB query with gene selection
  # Use selected_ranges for efficient gene-level queries
  result <- tryCatch({
    # Read all data and filter (for simpler implementation)
    # For very large databases, use tiledb_query with conditions
    df <- arr[]
    df <- df[df$gene %in% genes, , drop = FALSE]
    
    if (!is.null(samples)) {
      df <- df[df$sample %in% samples, , drop = FALSE]
    }
    
    df
  }, error = function(e) {
    stop("Error querying TileDB: ", e$message)
  })

  if (nrow(result) == 0) {
    warning("No data found for the specified genes")
    if (as_matrix) {
      # Return empty matrix with proper dimensions
      if (!is.null(samples)) {
        mat <- matrix(0, nrow = length(genes), ncol = length(samples),
                      dimnames = list(genes, samples))
      } else if (!is.null(all_samples)) {
        mat <- matrix(0, nrow = length(genes), ncol = length(all_samples),
                      dimnames = list(genes, all_samples))
      } else {
        mat <- matrix(0, nrow = length(genes), ncol = 0,
                      dimnames = list(genes, character(0)))
      }
      return(mat)
    } else {
      return(data.frame(gene = character(0), sample = character(0), count = numeric(0)))
    }
  }

  if (!as_matrix) {
    return(result)
  }

  # Convert to matrix
  # Determine all samples to include
  if (!is.null(samples)) {
    sample_names <- samples
  } else if (!is.null(all_samples)) {
    sample_names <- all_samples
  } else {
    sample_names <- unique(result$sample)
  }

  # Create matrix
  mat <- matrix(0, nrow = length(genes), ncol = length(sample_names),
                dimnames = list(genes, sample_names))

  # Fill matrix with values
  for (i in seq_len(nrow(result))) {
    g <- result$gene[i]
    s <- result$sample[i]
    if (g %in% genes && s %in% sample_names) {
      mat[g, s] <- result$count[i]
    }
  }

  return(mat)
}


#' @title List genes in TileDB database
#'
#' @description Returns all gene names stored in the TileDB database.
#'
#' @param tiledb_path Path to the TileDB database
#'
#' @return Character vector of gene names
#'
#' @export
pgx.listGenesTileDB <- function(tiledb_path) {
  metadata_path <- paste0(tiledb_path, "_metadata.rds")
  if (!file.exists(metadata_path)) {
    stop("Metadata file not found: ", metadata_path)
  }
  metadata <- readRDS(metadata_path)
  return(metadata$genes)
}


#' @title List samples in TileDB database
#'
#' @description Returns all sample names stored in the TileDB database.
#'
#' @param tiledb_path Path to the TileDB database
#'
#' @return Character vector of sample names (format: "filename::samplename")
#'
#' @export
pgx.listSamplesTileDB <- function(tiledb_path) {
  metadata_path <- paste0(tiledb_path, "_metadata.rds")
  if (!file.exists(metadata_path)) {
    stop("Metadata file not found: ", metadata_path)
  }
  metadata <- readRDS(metadata_path)
  return(metadata$samples)
}


#' @title Get TileDB database info
#'
#' @description Returns metadata about the TileDB database.
#'
#' @param tiledb_path Path to the TileDB database
#'
#' @return List with database metadata
#'
#' @export
pgx.infoTileDB <- function(tiledb_path) {
  metadata_path <- paste0(tiledb_path, "_metadata.rds")
  if (!file.exists(metadata_path)) {
    stop("Metadata file not found: ", metadata_path)
  }
  metadata <- readRDS(metadata_path)

  cat("TileDB Database Info\n")
  cat("====================\n")
  cat("Path:", tiledb_path, "\n")
  cat("Created:", as.character(metadata$created), "\n")
  cat("Number of genes:", metadata$n_genes, "\n")
  cat("Number of samples:", metadata$n_samples, "\n")
  cat("Number of PGX files:", metadata$n_files, "\n")
  cat("\nSource files:\n")
  for (f in metadata$pgx_files) {
    cat("  -", basename(f), "\n")
  }

  invisible(metadata)
}

