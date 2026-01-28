##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2025 BigOmics Analytics SA. All rights reserved.
##


## ============================================================================
## TileDB COUNTS DATABASE - DATA MODEL
## ============================================================================
##
## This module stores gene expression counts from multiple PGX files in a
## TileDB sparse array for efficient cross-dataset queries.
##
## FILE STRUCTURE
## --------------
##
##   <tiledb_path>/              <- TileDB sparse array (counts)
##   <tiledb_path>_metadata.rds  <- R metadata file
##
##
## COUNTS ARRAY (TileDB Sparse)
## ----------------------------
##
##   Dimensions:  gene (ASCII), sample (ASCII)
##   Attributes:  count (FLOAT64), zscore (FLOAT64)
##
##   Stored as coordinate (COO) format:
##
##       gene     | sample                   | count   | zscore
##       ---------+--------------------------+---------+---------
##       TP53     | TCGA_BRCA::TCGA-A1-A0SB  | 1523.5  | 1.23
##       TP53     | GSE12345::sample_001     | 892.3   | -0.45
##       BRCA1    | TCGA_BRCA::TCGA-A1-A0SB  | 234.1   | 0.78
##       ...      | ...                      | ...     | ...
##
##   - Genes are human orthologs (converted via pgx$genes$human_ortholog)
##   - Samples are prefixed with dataset name: "dataset::original_sample_id"
##   - Zeros ARE stored (distinguish from missing data which is NA)
##   - Z-scores are computed per-dataset from log-scaled X: (X - gene_mean) / gene_sd
##   - Z-scores are NA for genes with zero variance or single-sample datasets
##
##
## METADATA RDS FILE
## -----------------
##
##   list(
##     genes        = c("TP53", "BRCA1", ...),     # all gene symbols
##     samples      = c("ds1::s1", "ds1::s2", ...), # all sample IDs
##     phenotypes   = c("age", "sex", ...),        # all phenotype column names
##     phenotypes_by_dataset = list(               # phenotypes per dataset
##       "TCGA_BRCA" = c("age", "stage", ...),
##       "GSE12345" = c("treatment", ...)
##     ),
##     phenotype_data = data.frame(                # long format
##       sample    = c("ds1::s1", "ds1::s1", ...),
##       phenotype = c("age", "sex", ...),
##       value     = c("45", "male", ...)          # all values as strings
##     ),
##     n_genes, n_samples, n_phenotypes, n_files,
##     pgx_files,                                  # source file paths
##     created                                     # timestamp
##   )
##
##
## SAMPLE ID FORMAT
## ----------------
##
##   Full ID:    "dataset_name::original_sample_id"
##   Dataset:    getDatasetFromSample()  -> "dataset_name"
##   Short name: getSampleShortName()    -> "original_sample_id"
##
##   Dataset name = PGX filename without .pgx extension
##
## ============================================================================


#' @title Build TileDB database from PGX files
#'
#' @description Scans a folder for .pgx files, extracts counts matrices
#' and sample phenotypes, and stores them in a TileDB database.
#'
#' @param pgx_folder Path to folder containing .pgx files
#' @param tiledb_path Path where TileDB database will be created
#' @param overwrite Logical, whether to overwrite existing database. Default FALSE.
#' @param verbose Logical, whether to print progress messages. Default TRUE.
#'
#' @return Invisibly returns a list with metadata about the stored data
#'
#' @details This function creates:
#' \itemize{
#'   \item \code{tiledb_path}: TileDB sparse array with counts (genes x samples)
#'   \item \code{tiledb_path_metadata.rds}: RDS file with metadata and phenotype data
#' }
#'
#' Sample phenotypes are stored in the metadata file (not TileDB) for simplicity.
#' Since different datasets may have different phenotype columns, all values are
#' stored as strings for flexibility.
#'
#' @export
buildTileDB <- function(pgx_folder, tiledb_path, overwrite = FALSE, verbose = TRUE) {
  if (!requireNamespace("tiledb", quietly = TRUE)) {
    stop("Package 'tiledb' is required. Install with: install.packages('tiledb')")
  }

  if (!dir.exists(pgx_folder)) {
    stop("Folder does not exist: ", pgx_folder)
  }

  pgx_files <- list.files(pgx_folder, pattern = "\\.pgx$", full.names = TRUE)
  if (length(pgx_files) == 0) {
    stop("No .pgx files found in: ", pgx_folder)
  }

  if (verbose) message("Found ", length(pgx_files), " .pgx files")

  phenotypes_path <- paste0(tiledb_path, "_phenotypes")

  if (dir.exists(tiledb_path) || dir.exists(phenotypes_path)) {
    if (overwrite) {
      if (verbose) message("Removing existing TileDB databases...")
      unlink(tiledb_path, recursive = TRUE)
      unlink(phenotypes_path, recursive = TRUE)
    } else {
      stop("TileDB database already exists at: ", tiledb_path,
           "\nUse overwrite=TRUE to replace it.")
    }
  }

  if (verbose) message("First pass: collecting metadata...")

  all_genes <- character(0)
  all_samples <- character(0)
  all_phenotypes <- character(0)
  phenotypes_by_dataset <- list()

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
      counts <- rename_by2(counts, pgx$genes, new_id = "human_ortholog", keep.prefix = FALSE)

      file_prefix <- tools::file_path_sans_ext(basename(pgx_file))
      unique_samples <- paste0(file_prefix, "::", colnames(counts))

      all_genes <- union(all_genes, rownames(counts))
      all_samples <- c(all_samples, unique_samples)

      if (!is.null(pgx$samples) && is.data.frame(pgx$samples)) {
        pheno_cols <- colnames(pgx$samples)
        pheno_cols <- pheno_cols[!grepl("^\\.", pheno_cols)]
        all_phenotypes <- union(all_phenotypes, pheno_cols)
        phenotypes_by_dataset[[file_prefix]] <- pheno_cols
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

  all_genes <- sort(unique(all_genes))
  all_samples <- unique(all_samples)
  all_phenotypes <- sort(unique(all_phenotypes))

  if (verbose) {
    message("Total unique genes: ", length(all_genes))
    message("Total samples: ", length(all_samples))
    message("Total unique phenotype columns: ", length(all_phenotypes))
  }

  if (verbose) message("Creating TileDB schema...")

  gene_dim <- tiledb::tiledb_dim(name = "gene", domain = NULL, tile = NULL, type = "ASCII")
  sample_dim <- tiledb::tiledb_dim(name = "sample", domain = NULL, tile = NULL, type = "ASCII")
  dom <- tiledb::tiledb_domain(dims = list(gene_dim, sample_dim))
  count_attr <- tiledb::tiledb_attr(name = "count", type = "FLOAT64")
  zscore_attr <- tiledb::tiledb_attr(name = "zscore", type = "FLOAT64")
  schema <- tiledb::tiledb_array_schema(domain = dom, attrs = list(count_attr, zscore_attr), sparse = TRUE)
  tiledb::tiledb_array_create(tiledb_path, schema)

  if (verbose) message("TileDB array created at: ", tiledb_path)

  if (verbose) message("Second pass: writing data...")

  arr <- tiledb::tiledb_array(tiledb_path)
  all_phenotype_data <- list()

  for (i in seq_along(pgx_files)) {
    pgx_file <- pgx_files[i]
    if (verbose) message("  [", i, "/", length(pgx_files), "] Writing: ", basename(pgx_file))

    tryCatch({
      pgx <- pgx.load(pgx_file)
      if (is.null(pgx$counts)) next

      file_prefix <- tools::file_path_sans_ext(basename(pgx_file))

      ## Prepare and write counts/z-scores data
      data <- tiledb.prepareData(pgx)
      tiledb.writeData(arr, data, file_prefix)

      ## Extract phenotype data
      pheno <- tiledb.extractPhenotypes(pgx, file_prefix)
      all_phenotype_data <- c(all_phenotype_data, pheno$phenotype_data)

      rm(pgx, data)
      gc(verbose = FALSE)
    }, error = function(e) {
      warning("Error writing ", pgx_file, ": ", e$message)
    })
  }

  phenotype_data <- if (length(all_phenotype_data) > 0) {
    do.call(rbind, all_phenotype_data)
  } else {
    data.frame(sample = character(0), phenotype = character(0), value = character(0))
  }

  if (verbose) message("Writing metadata...")

  metadata_path <- paste0(tiledb_path, "_metadata.rds")
  metadata <- list(
    genes = all_genes,
    samples = all_samples,
    phenotypes = all_phenotypes,
    phenotypes_by_dataset = phenotypes_by_dataset,
    phenotype_data = phenotype_data,
    n_genes = length(all_genes),
    n_samples = length(all_samples),
    n_phenotypes = length(all_phenotypes),
    n_files = length(pgx_files),
    pgx_files = pgx_files,
    created = Sys.time()
  )
  saveRDS(metadata, metadata_path)

  if (verbose) {
    message("Done! TileDB database created at: ", tiledb_path)
    message("  - Genes: ", length(all_genes))
    message("  - Samples: ", length(all_samples))
    message("  - Phenotype columns: ", length(all_phenotypes))
    message("  - Files: ", length(pgx_files))
  }

  invisible(metadata)
}


#' @title Query gene counts or z-scores from TileDB database
#'
#' @param tiledb_path Path to the TileDB database
#' @param genes Character vector of gene names to query
#' @param samples Optional character vector of sample names to filter
#' @param as_matrix Return as matrix (TRUE) or data.frame (FALSE)
#' @param value Which value to return: "count" for raw counts, "zscore" for
#'   per-dataset standardized values. Default is "count".
#'
#' @return Matrix (genes x samples) or data.frame with columns: gene, sample, and the selected value
#'
#' @export
queryTileDB <- function(tiledb_path, genes, samples = NULL, as_matrix = TRUE,
                        value = c("count", "zscore")) {
  if (!requireNamespace("tiledb", quietly = TRUE)) {
    stop("Package 'tiledb' is required. Install with: install.packages('tiledb')")
  }

  value <- match.arg(value)

  if (!dir.exists(tiledb_path)) {
    stop("TileDB database not found at: ", tiledb_path)
  }

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

  arr <- tiledb::tiledb_array(tiledb_path, return_as = "data.frame")

  df <- arr[]
  df <- df[df$gene %in% genes, , drop = FALSE]
  if (!is.null(samples)) {
    df <- df[df$sample %in% samples, , drop = FALSE]
  }

  if (nrow(df) == 0) {
    warning("No data found for the specified genes")
    if (as_matrix) {
      ncols <- if (!is.null(samples)) length(samples) else if (!is.null(all_samples)) length(all_samples) else 0
      colnms <- if (!is.null(samples)) samples else if (!is.null(all_samples)) all_samples else character(0)
      return(matrix(NA_real_, nrow = length(genes), ncol = ncols, dimnames = list(genes, colnms)))
    } else {
      return(data.frame(gene = character(0), sample = character(0), value = numeric(0)))
    }
  }

  ## Select the requested value column
  df$value <- df[[value]]

  if (!as_matrix) {
    return(df[, c("gene", "sample", value), drop = FALSE])
  }

  sample_names <- if (!is.null(samples)) samples else if (!is.null(all_samples)) all_samples else unique(df$sample)

  mat <- matrix(NA_real_, nrow = length(genes), ncol = length(sample_names),
                dimnames = list(genes, sample_names))

  for (i in seq_len(nrow(df))) {
    g <- df$gene[i]
    s <- df$sample[i]
    if (g %in% genes && s %in% sample_names) {
      mat[g, s] <- df$value[i]
    }
  }

  mat
}


#' @title List genes in TileDB database
#' @param tiledb_path Path to the TileDB database
#' @return Character vector of gene names
#' @export
listGenesTileDB <- function(tiledb_path) {
  metadata <- readRDS(paste0(tiledb_path, "_metadata.rds"))
  metadata$genes
}


#' @title List samples in TileDB database
#' @param tiledb_path Path to the TileDB database
#' @return Character vector of sample names (format: "filename::samplename")
#' @export
listSamplesTileDB <- function(tiledb_path) {
  metadata <- readRDS(paste0(tiledb_path, "_metadata.rds"))
  metadata$samples
}


#' @title List datasets in TileDB database
#' @param tiledb_path Path to the TileDB database
#' @return Character vector of dataset names
#' @export
listDatasetsTileDB <- function(tiledb_path) {
  metadata <- readRDS(paste0(tiledb_path, "_metadata.rds"))
  sort(tools::file_path_sans_ext(basename(metadata$pgx_files)))
}


#' @title Get dataset information table from TileDB database
#'
#' @description Returns a data.frame with metadata for datasets in the TileDB database.
#' Cross-references with datasets-info.csv if provided for additional metadata.
#'
#' @param tiledb_path Path to the TileDB database
#' @param datasets_info_file Optional path to datasets-info.csv file for metadata lookup
#'
#' @return Data.frame with columns: dataset, nsamples, and optionally organism, datatype, description, date
#' @export
getDatasetInfoTileDB <- function(tiledb_path, datasets_info_file = NULL) {
  metadata <- readRDS(paste0(tiledb_path, "_metadata.rds"))

  datasets <- tools::file_path_sans_ext(basename(metadata$pgx_files))

  ## Count samples per dataset
  sample_datasets <- getDatasetFromSample(metadata$samples)
  samples_per_dataset <- table(sample_datasets)

  ## Build base info data.frame

  df <- data.frame(
    dataset = datasets,
    nsamples = as.integer(samples_per_dataset[datasets]),
    stringsAsFactors = FALSE
  )

  ## Cross-reference with datasets-info.csv if available
  if (!is.null(datasets_info_file) && file.exists(datasets_info_file)) {
    info_csv <- utils::read.csv(datasets_info_file, stringsAsFactors = FALSE)
    ## Normalize dataset names (remove .pgx extension if present)
    if ("dataset" %in% colnames(info_csv)) {
      info_csv$dataset <- sub("[.]pgx$", "", info_csv$dataset)
    }
    ## Select useful columns (including metadata_* columns)
    base_cols <- c("dataset", "organism", "datatype", "description", "date", "nfeatures", "conditions")
    metadata_cols <- grep("^metadata_", colnames(info_csv), value = TRUE)
    info_cols <- c(base_cols, metadata_cols)
    info_cols <- intersect(info_cols, colnames(info_csv))
    if (length(info_cols) > 1) {
      info_csv <- info_csv[, info_cols, drop = FALSE]
      df <- merge(df, info_csv, by = "dataset", all.x = TRUE)
    }
  }

  df <- df[order(df$dataset), ]
  rownames(df) <- NULL

  df
}


#' @title List datasets in TileDB database
#'
#' @description Returns all unique dataset names stored in the TileDB database.
#'
#' @param tiledb_path Path to the TileDB database
#'
#' @return Character vector of dataset names (sorted alphabetically)
#'
#' @export
listDatasetsTileDB <- function(tiledb_path) {
  metadata_path <- paste0(tiledb_path, "_metadata.rds")
  if (!file.exists(metadata_path)) {
    stop("Metadata file not found: ", metadata_path)
  }
  metadata <- readRDS(metadata_path)
  datasets <- tools::file_path_sans_ext(basename(metadata$pgx_files))
  return(sort(datasets))
}


#' @title Get TileDB database info
#' @param tiledb_path Path to the TileDB database
#' @return List with database metadata (invisibly)
#' @export
infoTileDB <- function(tiledb_path) {
  metadata <- readRDS(paste0(tiledb_path, "_metadata.rds"))

  cat("TileDB Database Info\n")
  cat("====================\n")
  cat("Path:", tiledb_path, "\n")
  cat("Created:", as.character(metadata$created), "\n")
  cat("Genes:", metadata$n_genes, "\n")
  cat("Samples:", metadata$n_samples, "\n")
  cat("Phenotypes:", metadata$n_phenotypes, "\n")
  cat("Files:", metadata$n_files, "\n")
  cat("\nSource files:\n")
  for (f in metadata$pgx_files) cat("  -", basename(f), "\n")
  if (length(metadata$phenotypes) > 0) {
    cat("\nPhenotype columns:\n  ", paste(head(metadata$phenotypes, 10), collapse = ", "))
    if (length(metadata$phenotypes) > 10) cat(" ... and", length(metadata$phenotypes) - 10, "more")
    cat("\n")
  }

  invisible(metadata)
}


#' @title List phenotype columns in TileDB database
#' @param tiledb_path Path to the TileDB database
#' @return Character vector of phenotype column names
#' @export
listPhenotypesTileDB <- function(tiledb_path) {
  metadata <- readRDS(paste0(tiledb_path, "_metadata.rds"))
  if (is.null(metadata$phenotypes)) return(character(0))
  metadata$phenotypes
}


#' @title List phenotype columns by dataset
#' @param tiledb_path Path to the TileDB database
#' @param dataset Optional dataset name. If NULL, returns list for all datasets.
#' @return Named list of phenotype columns per dataset, or character vector if dataset specified
#' @export
listPhenotypesByDatasetTileDB <- function(tiledb_path, dataset = NULL) {
  metadata <- readRDS(paste0(tiledb_path, "_metadata.rds"))

  if (is.null(metadata$phenotypes_by_dataset)) return(list())
  if (is.null(dataset)) return(metadata$phenotypes_by_dataset)
  if (!dataset %in% names(metadata$phenotypes_by_dataset)) {
    warning("Dataset not found: ", dataset)
    return(character(0))
  }
  metadata$phenotypes_by_dataset[[dataset]]
}


#' @title Query phenotype data from TileDB database
#'
#' @param tiledb_path Path to the TileDB database
#' @param samples Optional sample names to query
#' @param phenotypes Optional phenotype column names to retrieve
#' @param datasets Optional dataset names to filter by
#' @param as_wide Return as wide data.frame (TRUE) or long format (FALSE)
#'
#' @return Data.frame with samples as rows and phenotypes as columns (wide),
#'   or columns: sample, phenotype, value (long)
#'
#' @export
queryPhenotypesTileDB <- function(tiledb_path, samples = NULL, phenotypes = NULL,
                                       datasets = NULL, as_wide = TRUE) {
  metadata <- readRDS(paste0(tiledb_path, "_metadata.rds"))

  if (is.null(metadata$phenotype_data) || nrow(metadata$phenotype_data) == 0) {
    warning("No phenotype data found in metadata")
    if (as_wide) return(data.frame(row.names = character(0)))
    return(data.frame(sample = character(0), phenotype = character(0), value = character(0)))
  }

  all_samples <- metadata$samples
  all_phenotypes <- metadata$phenotypes
  result <- metadata$phenotype_data

  if (!is.null(datasets)) {
    sample_datasets <- getDatasetFromSample(all_samples)
    samples_in_datasets <- all_samples[sample_datasets %in% datasets]
    samples <- if (is.null(samples)) samples_in_datasets else intersect(samples, samples_in_datasets)
  }

  if (!is.null(phenotypes)) {
    missing_pheno <- setdiff(phenotypes, all_phenotypes)
    if (length(missing_pheno) > 0) {
      warning("Some phenotypes not found: ", paste(head(missing_pheno, 5), collapse = ", "))
      phenotypes <- intersect(phenotypes, all_phenotypes)
    }
  }

  if (!is.null(samples)) result <- result[result$sample %in% samples, , drop = FALSE]
  if (!is.null(phenotypes)) result <- result[result$phenotype %in% phenotypes, , drop = FALSE]

  if (nrow(result) == 0) {
    warning("No phenotype data found for the specified query")
    if (as_wide) return(data.frame(row.names = character(0)))
    return(data.frame(sample = character(0), phenotype = character(0), value = character(0)))
  }

  if (!as_wide) return(result)

  sample_names <- if (!is.null(samples)) samples else unique(result$sample)
  pheno_names <- if (!is.null(phenotypes)) phenotypes else unique(result$phenotype)

  wide_df <- data.frame(row.names = sample_names, stringsAsFactors = FALSE)
  for (pheno in pheno_names) wide_df[[pheno]] <- NA_character_

  for (i in seq_len(nrow(result))) {
    s <- result$sample[i]
    p <- result$phenotype[i]
    if (s %in% sample_names && p %in% pheno_names) {
      wide_df[s, p] <- result$value[i]
    }
  }

  for (col in colnames(wide_df)) {
    vals <- wide_df[[col]]
    numeric_vals <- suppressWarnings(as.numeric(vals))
    if (!all(is.na(numeric_vals[!is.na(vals)])) && all(is.na(vals) | !is.na(numeric_vals))) {
      wide_df[[col]] <- numeric_vals
    }
  }

  wide_df$dataset <- getDatasetFromSample(rownames(wide_df))
  wide_df$sample_short <- getSampleShortName(rownames(wide_df))
  wide_df
}


#' @title Get unique values for a phenotype column
#' @param tiledb_path Path to the TileDB database
#' @param phenotype Name of the phenotype column
#' @param datasets Optional dataset names to filter by
#' @return Character vector of unique values for the phenotype
#' @export
getPhenotypeValuesTileDB <- function(tiledb_path, phenotype, datasets = NULL) {
  metadata <- readRDS(paste0(tiledb_path, "_metadata.rds"))

  if (is.null(metadata$phenotype_data) || nrow(metadata$phenotype_data) == 0) {
    return(character(0))
  }

  result <- metadata$phenotype_data[metadata$phenotype_data$phenotype == phenotype, , drop = FALSE]
  if (nrow(result) == 0) return(character(0))

  if (!is.null(datasets)) {
    sample_datasets <- getDatasetFromSample(result$sample)
    result <- result[sample_datasets %in% datasets, , drop = FALSE]
  }

  values <- unique(result$value)
  sort(values[!is.na(values)])
}


#' @title Extract dataset name from sample identifiers
#' @param samples Character vector of sample identifiers ("dataset::sample")
#' @return Character vector of dataset names
#' @export
getDatasetFromSample <- function(samples) {
  sub("::.*", "", samples)
}


#' @title Extract short sample name from sample identifiers
#' @param samples Character vector of sample identifiers ("dataset::sample")
#' @return Character vector of short sample names
#' @export
getSampleShortName <- function(samples) {
  sub(".*::", "", samples)
}


#' @title Convert TileDB query result to plotting data.frame
#'
#' @param result Matrix or data.frame from queryTileDB
#' @param gene_name Optional gene name for single-gene matrix input
#' @param tiledb_path Optional path to TileDB database for merging phenotype data
#' @param phenotypes Optional phenotype columns to include
#'
#' @return Data.frame with columns: sample, count, dataset, sample_short, gene, phenotypes
#'
#' @export
tiledbToPlotDF <- function(result, gene_name = NULL, tiledb_path = NULL, phenotypes = NULL) {
  if (is.matrix(result) || inherits(result, "Matrix")) {
    if (nrow(result) == 1) {
      df <- data.frame(
        sample = colnames(result),
        count = as.numeric(result[1, ]),
        gene = if (!is.null(gene_name)) gene_name else rownames(result)[1],
        stringsAsFactors = FALSE
      )
    } else {
      df <- data.frame(
        gene = rep(rownames(result), ncol(result)),
        sample = rep(colnames(result), each = nrow(result)),
        count = as.numeric(result),
        stringsAsFactors = FALSE
      )
    }
  } else if (is.data.frame(result)) {
    df <- result
  } else {
    stop("Input must be a matrix or data.frame")
  }

  df$dataset <- getDatasetFromSample(df$sample)
  df$sample_short <- getSampleShortName(df$sample)

  if (!is.null(tiledb_path) && !is.null(phenotypes)) {
    pheno_df <- queryPhenotypesTileDB(
      tiledb_path, samples = unique(df$sample), phenotypes = phenotypes, as_wide = TRUE
    )
    if (nrow(pheno_df) > 0) {
      pheno_df$dataset <- NULL
      pheno_df$sample_short <- NULL
      pheno_df$sample <- rownames(pheno_df)
      df <- merge(df, pheno_df, by = "sample", all.x = TRUE)
    }
  }

  df
}


## ============================================================================
## INTERNAL HELPER FUNCTIONS
## ============================================================================


#' @title Extract counts and z-scores from PGX object
#' @description Prepares counts matrix and computes z-scores from log-scaled X matrix
#' @param pgx PGX object with counts and X matrices
#' @return List with counts, zscores, genes, and samples
#' @keywords internal
tiledb.prepareData <- function(pgx) {
  if (is.null(pgx$counts)) {
    stop("No counts found in PGX object")
  }

  ## Extract and transform counts
  counts <- as.matrix(pgx$counts)
  counts <- rename_by2(counts, pgx$genes, new_id = "human_ortholog", keep.prefix = FALSE)

  ## Compute z-scores from log-scaled X matrix
  X <- as.matrix(pgx$X)
  X <- rename_by2(X, pgx$genes, new_id = "human_ortholog", keep.prefix = FALSE)
  gene_means <- rowMeans(X, na.rm = TRUE)
  gene_sds <- apply(X, 1, sd, na.rm = TRUE)
  gene_sds[gene_sds == 0] <- NA  ## Avoid division by zero
  zscores <- (X - gene_means) / gene_sds

  list(
    counts = counts,
    zscores = zscores,
    genes = rownames(counts),
    samples = colnames(counts)
  )
}


#' @title Write counts and z-scores to TileDB array
#' @description Writes data for all samples to TileDB sparse array
#' @param arr TileDB array object
#' @param data List from tiledb.prepareData()
#' @param file_prefix Dataset name prefix for sample IDs
#' @keywords internal
tiledb.writeData <- function(arr, data, file_prefix) {
  genes <- data$genes
  samples <- data$samples
  counts <- data$counts
  zscores <- data$zscores

  for (j in seq_along(samples)) {
    sample_name <- paste0(file_prefix, "::", samples[j])
    count_values <- counts[, j]
    zscore_values <- zscores[, j]
    valid_idx <- which(!is.na(count_values))

    if (length(valid_idx) > 0) {
      df <- data.frame(
        gene = genes[valid_idx],
        sample = rep(sample_name, length(valid_idx)),
        count = as.numeric(count_values[valid_idx]),
        zscore = as.numeric(zscore_values[valid_idx]),
        stringsAsFactors = FALSE
      )
      arr[] <- df
    }
  }
}


#' @title Extract phenotype data from PGX object
#' @description Extracts phenotype data in long format
#' @param pgx PGX object with samples data.frame
#' @param file_prefix Dataset name prefix for sample IDs
#' @return List with phenotype_data (list of data.frames) and phenotypes (column names)
#' @keywords internal
tiledb.extractPhenotypes <- function(pgx, file_prefix) {
  phenotype_data <- list()
  phenotypes <- character(0)

  if (!is.null(pgx$samples) && is.data.frame(pgx$samples)) {
    pheno_data <- pgx$samples
    pheno_cols <- colnames(pheno_data)
    pheno_cols <- pheno_cols[!grepl("^\\.", pheno_cols)]
    phenotypes <- pheno_cols

    if (length(pheno_cols) > 0) {
      for (j in seq_len(nrow(pheno_data))) {
        sample_id <- rownames(pheno_data)[j]
        sample_name <- paste0(file_prefix, "::", sample_id)

        for (col in pheno_cols) {
          val <- pheno_data[j, col]
          val_str <- if (is.na(val)) NA_character_ else as.character(val)
          phenotype_data[[length(phenotype_data) + 1]] <- data.frame(
            sample = sample_name,
            phenotype = col,
            value = val_str,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  list(
    phenotype_data = phenotype_data,
    phenotypes = phenotypes
  )
}


## ============================================================================
## INCREMENTAL UPDATE FUNCTIONS
## ============================================================================


#' @title Add a dataset to existing TileDB database
#'
#' @description Adds a single PGX file to an existing TileDB database.
#' If the dataset already exists, it can be overwritten. If the TileDB
#' database doesn't exist, it will be created.
#'
#' @param tiledb_path Path to the TileDB database
#' @param pgx_file Path to .pgx file to add
#' @param overwrite If TRUE and dataset already exists, remove old data first. Default TRUE.
#' @param verbose Logical, whether to print progress messages. Default TRUE.
#'
#' @return Invisibly returns updated metadata
#'
#' @details This function enables incremental updates to the TileDB database
#' without rebuilding from scratch. It:
#' \itemize{
#'   \item Loads the pgx file and extracts counts/phenotypes
#'   \item Appends new data to the TileDB sparse array
#'   \item Updates the metadata RDS file
#' }
#'
#' @export
tiledb.addDataset <- function(tiledb_path, pgx_file, overwrite = TRUE, verbose = TRUE) {
  if (!requireNamespace("tiledb", quietly = TRUE)) {
    stop("Package 'tiledb' is required. Install with: install.packages('tiledb')")
  }

  ## Normalize pgx_file path
  if (!file.exists(pgx_file)) {
    stop("PGX file does not exist: ", pgx_file)
  }

  metadata_path <- paste0(tiledb_path, "_metadata.rds")

  ## If TileDB doesn't exist, create it with just this file

  if (!dir.exists(tiledb_path)) {
    if (verbose) message("TileDB database doesn't exist. Creating new database...")
    pgx_folder <- dirname(pgx_file)
    ## Create with just this single file by building from folder but filtering
    return(tiledb.addDataset.create(tiledb_path, pgx_file, verbose))
  }

  ## Load existing metadata
  if (!file.exists(metadata_path)) {
    stop("Metadata file not found: ", metadata_path)
  }
  metadata <- readRDS(metadata_path)

  ## Get dataset name from filename
  dataset_name <- tools::file_path_sans_ext(basename(pgx_file))

  ## Check if dataset already exists
  existing_datasets <- tools::file_path_sans_ext(basename(metadata$pgx_files))
  dataset_exists <- dataset_name %in% existing_datasets

  if (dataset_exists) {
    if (!overwrite) {
      if (verbose) message("Dataset '", dataset_name, "' already exists. Use overwrite=TRUE to replace.")
      return(invisible(metadata))
    }
    if (verbose) message("Removing existing data for dataset: ", dataset_name)
    metadata <- tiledb.removeDataset.internal(tiledb_path, dataset_name, metadata, verbose = FALSE)
  }

  ## Load the PGX file
  if (verbose) message("Loading PGX file: ", basename(pgx_file))

  pgx <- tryCatch({
    pgx.load(pgx_file)
  }, error = function(e) {
    stop("Error loading PGX file: ", e$message)
  })

  ## Prepare and write counts/z-scores data
  file_prefix <- dataset_name
  data <- tiledb.prepareData(pgx)

  if (verbose) message("Adding ", length(data$samples), " samples, ", length(data$genes), " genes")

  arr <- tiledb::tiledb_array(tiledb_path)
  tiledb.writeData(arr, data, file_prefix)

  ## Extract phenotype data
  pheno <- tiledb.extractPhenotypes(pgx, file_prefix)

  ## Update metadata
  new_samples <- paste0(file_prefix, "::", data$samples)

  metadata$genes <- sort(unique(c(metadata$genes, data$genes)))
  metadata$samples <- c(metadata$samples, new_samples)
  metadata$phenotypes <- sort(unique(c(metadata$phenotypes, pheno$phenotypes)))
  metadata$phenotypes_by_dataset[[file_prefix]] <- pheno$phenotypes

  if (length(pheno$phenotype_data) > 0) {
    new_pheno_df <- do.call(rbind, pheno$phenotype_data)
    metadata$phenotype_data <- rbind(metadata$phenotype_data, new_pheno_df)
  }

  metadata$pgx_files <- c(metadata$pgx_files, normalizePath(pgx_file))
  metadata$n_genes <- length(metadata$genes)
  metadata$n_samples <- length(metadata$samples)
  metadata$n_phenotypes <- length(metadata$phenotypes)
  metadata$n_files <- length(metadata$pgx_files)
  metadata$updated <- Sys.time()

  ## Save updated metadata
  saveRDS(metadata, metadata_path)

  if (verbose) {
    message("Dataset '", dataset_name, "' added successfully!")
    message("  Total genes: ", metadata$n_genes)
    message("  Total samples: ", metadata$n_samples)
    message("  Total datasets: ", metadata$n_files)
  }

  rm(pgx, data)
  gc(verbose = FALSE)

  invisible(metadata)
}


#' @title Create new TileDB database with single file (internal)
#' @keywords internal
tiledb.addDataset.create <- function(tiledb_path, pgx_file, verbose = TRUE) {
  if (!requireNamespace("tiledb", quietly = TRUE)) {
    stop("Package 'tiledb' is required. Install with: install.packages('tiledb')")
  }

  ## Load the PGX file
  if (verbose) message("Loading PGX file: ", basename(pgx_file))

  pgx <- tryCatch({
    pgx.load(pgx_file)
  }, error = function(e) {
    stop("Error loading PGX file: ", e$message)
  })

  ## Prepare data
  file_prefix <- tools::file_path_sans_ext(basename(pgx_file))
  data <- tiledb.prepareData(pgx)

  if (verbose) message("Creating TileDB with ", length(data$samples), " samples, ", length(data$genes), " genes")

  ## Create TileDB schema
  gene_dim <- tiledb::tiledb_dim(name = "gene", domain = NULL, tile = NULL, type = "ASCII")
  sample_dim <- tiledb::tiledb_dim(name = "sample", domain = NULL, tile = NULL, type = "ASCII")
  dom <- tiledb::tiledb_domain(dims = list(gene_dim, sample_dim))
  count_attr <- tiledb::tiledb_attr(name = "count", type = "FLOAT64")
  zscore_attr <- tiledb::tiledb_attr(name = "zscore", type = "FLOAT64")
  schema <- tiledb::tiledb_array_schema(domain = dom, attrs = list(count_attr, zscore_attr), sparse = TRUE)
  tiledb::tiledb_array_create(tiledb_path, schema)

  ## Write data
  arr <- tiledb::tiledb_array(tiledb_path)
  tiledb.writeData(arr, data, file_prefix)

  ## Extract phenotype data
  pheno <- tiledb.extractPhenotypes(pgx, file_prefix)

  phenotype_data <- if (length(pheno$phenotype_data) > 0) {
    do.call(rbind, pheno$phenotype_data)
  } else {
    data.frame(sample = character(0), phenotype = character(0), value = character(0))
  }

  ## Create metadata
  all_samples <- paste0(file_prefix, "::", data$samples)
  phenotypes_by_dataset <- list()
  phenotypes_by_dataset[[file_prefix]] <- pheno$phenotypes

  metadata_path <- paste0(tiledb_path, "_metadata.rds")
  metadata <- list(
    genes = sort(unique(data$genes)),
    samples = all_samples,
    phenotypes = sort(unique(pheno$phenotypes)),
    phenotypes_by_dataset = phenotypes_by_dataset,
    phenotype_data = phenotype_data,
    n_genes = length(unique(data$genes)),
    n_samples = length(all_samples),
    n_phenotypes = length(unique(pheno$phenotypes)),
    n_files = 1,
    pgx_files = normalizePath(pgx_file),
    created = Sys.time()
  )
  saveRDS(metadata, metadata_path)

  if (verbose) {
    message("TileDB database created at: ", tiledb_path)
    message("  Genes: ", metadata$n_genes)
    message("  Samples: ", metadata$n_samples)
    message("  Datasets: 1")
  }

  rm(pgx, data)
  gc(verbose = FALSE)

  invisible(metadata)
}


#' @title Remove dataset from TileDB (internal helper)
#' @keywords internal
tiledb.removeDataset.internal <- function(tiledb_path, dataset_name, metadata, verbose = TRUE) {
  ## Remove samples belonging to this dataset from metadata
  sample_datasets <- getDatasetFromSample(metadata$samples)
  keep_samples <- metadata$samples[sample_datasets != dataset_name]

  ## Remove phenotype data for this dataset
  if (!is.null(metadata$phenotype_data) && nrow(metadata$phenotype_data) > 0) {
    pheno_datasets <- getDatasetFromSample(metadata$phenotype_data$sample)
    metadata$phenotype_data <- metadata$phenotype_data[pheno_datasets != dataset_name, , drop = FALSE]
  }

  ## Remove from phenotypes_by_dataset
  metadata$phenotypes_by_dataset[[dataset_name]] <- NULL

  ## Remove from pgx_files
  pgx_basenames <- tools::file_path_sans_ext(basename(metadata$pgx_files))
  keep_files <- metadata$pgx_files[pgx_basenames != dataset_name]

  ## Update metadata
  metadata$samples <- keep_samples
  metadata$pgx_files <- keep_files
  metadata$n_samples <- length(keep_samples)
  metadata$n_files <- length(keep_files)

  ## Note: We don't actually delete from TileDB sparse array - old data remains

  ## but is effectively orphaned since metadata no longer references it.
  ## A full rebuild would clean this up if storage becomes an issue.

  if (verbose) message("Removed dataset '", dataset_name, "' from metadata")

  metadata
}


#' @title Check if TileDB database needs update
#'
#' @description Checks if the TileDB database is missing any .pgx files
#' from the specified directory.
#'
#' @param pgx_dir Directory containing .pgx files
#' @param tiledb_path Path to TileDB database. If NULL, defaults to pgx_dir/counts_tiledb
#'
#' @return Logical indicating if update is needed (TRUE if TileDB doesn't exist
#'   or is missing datasets)
#'
#' @export
tiledb.needUpdate <- function(pgx_dir, tiledb_path = NULL) {
  if (!dir.exists(pgx_dir)) {
    return(FALSE)
  }

  ## Get pgx files in directory - if none, nothing to update
  pgx_files <- list.files(pgx_dir, pattern = "\\.pgx$", full.names = FALSE)
  if (length(pgx_files) == 0) {
    return(FALSE)
  }

  if (is.null(tiledb_path)) {
    tiledb_path <- file.path(pgx_dir, "counts_tiledb")
  }

  metadata_path <- paste0(tiledb_path, "_metadata.rds")

  ## If TileDB doesn't exist, needs update
  if (!dir.exists(tiledb_path) || !file.exists(metadata_path)) {
    return(TRUE)
  }
  pgx_datasets <- tools::file_path_sans_ext(pgx_files)

  ## Get datasets in TileDB
  metadata <- tryCatch({
    readRDS(metadata_path)
  }, error = function(e) {
    return(TRUE)  ## Can't read metadata, needs rebuild
  })

  if (is.null(metadata) || is.null(metadata$pgx_files)) {
    return(TRUE)
  }

  tiledb_datasets <- tools::file_path_sans_ext(basename(metadata$pgx_files))

  ## Check if any pgx files are missing from TileDB
  missing <- setdiff(pgx_datasets, tiledb_datasets)

  return(length(missing) > 0)
}


#' @title Update TileDB database for a folder
#'
#' @description Updates the TileDB database to include any new or modified
#' .pgx files. If a specific file is provided via new_pgx, only that file
#' is added/updated.
#'
#' @param pgx_dir Directory containing .pgx files
#' @param tiledb_path Path to TileDB database. If NULL, defaults to pgx_dir/counts_tiledb
#' @param new_pgx Optional. Specific .pgx filename that was just added/modified.
#'   If provided, only this file is processed.
#' @param verbose Logical, whether to print progress messages. Default TRUE.
#'
#' @return Invisibly returns the updated metadata, or NULL if no update needed
#'
#' @details This function mirrors the behavior of pgxinfo.updateDatasetFolder()
#' for TileDB databases. It:
#' \itemize{
#'   \item Creates the TileDB database if it doesn't exist
#'   \item Adds specific new datasets when new_pgx is provided
#'   \item Or scans for all missing datasets and adds them
#' }
#'
#' @export
tiledb.updateDatasetFolder <- function(pgx_dir, tiledb_path = NULL, new_pgx = NULL, verbose = TRUE) {
  if (!dir.exists(pgx_dir)) {
    if (verbose) message("[tiledb.updateDatasetFolder] Directory does not exist: ", pgx_dir)
    return(invisible(NULL))
  }

  if (is.null(tiledb_path)) {
    tiledb_path <- file.path(pgx_dir, "counts_tiledb")
  }

  metadata_path <- paste0(tiledb_path, "_metadata.rds")

  ## If specific file provided, just add that one
  if (!is.null(new_pgx)) {
    new_pgx <- sub("\\.pgx$", "", new_pgx)  ## strip extension if present
    pgx_file <- file.path(pgx_dir, paste0(new_pgx, ".pgx"))

    if (!file.exists(pgx_file)) {
      if (verbose) message("[tiledb.updateDatasetFolder] PGX file not found: ", pgx_file)
      return(invisible(NULL))
    }

    if (verbose) message("[tiledb.updateDatasetFolder] Adding dataset: ", new_pgx)
    result <- tiledb.addDataset(tiledb_path, pgx_file, overwrite = TRUE, verbose = verbose)
    return(invisible(result))
  }

  ## Otherwise, check what's missing and add all
  pgx_files <- list.files(pgx_dir, pattern = "\\.pgx$", full.names = TRUE)
  if (length(pgx_files) == 0) {
    if (verbose) message("[tiledb.updateDatasetFolder] No .pgx files found in: ", pgx_dir)
    return(invisible(NULL))
  }

  pgx_datasets <- tools::file_path_sans_ext(basename(pgx_files))

  ## Get existing datasets in TileDB (if exists)
  tiledb_datasets <- character(0)
  if (dir.exists(tiledb_path) && file.exists(metadata_path)) {
    metadata <- tryCatch({
      readRDS(metadata_path)
    }, error = function(e) NULL)

    if (!is.null(metadata) && !is.null(metadata$pgx_files)) {
      tiledb_datasets <- tools::file_path_sans_ext(basename(metadata$pgx_files))
    }
  }

  ## Find missing datasets
  missing_datasets <- setdiff(pgx_datasets, tiledb_datasets)

  if (length(missing_datasets) == 0) {
    if (verbose) message("[tiledb.updateDatasetFolder] TileDB is up to date")
    return(invisible(NULL))
  }

  if (verbose) message("[tiledb.updateDatasetFolder] Adding ", length(missing_datasets), " missing dataset(s)")

  ## Add each missing dataset
  result <- NULL
  for (ds in missing_datasets) {
    pgx_file <- file.path(pgx_dir, paste0(ds, ".pgx"))
    if (file.exists(pgx_file)) {
      if (verbose) message("  Adding: ", ds)
      result <- tryCatch({
        tiledb.addDataset(tiledb_path, pgx_file, overwrite = TRUE, verbose = FALSE)
      }, error = function(e) {
        warning("Error adding dataset '", ds, "': ", e$message)
        NULL
      })
    }
  }

  if (verbose && !is.null(result)) {
    message("[tiledb.updateDatasetFolder] Done! Total datasets: ", result$n_files)
  }

  invisible(result)
}

