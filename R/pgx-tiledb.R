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
##   Attributes:  count (FLOAT64)
##
##   Stored as coordinate (COO) format:
##
##       gene     | sample                   | count
##       ---------+--------------------------+--------
##       TP53     | TCGA_BRCA::TCGA-A1-A0SB  | 1523.5
##       TP53     | GSE12345::sample_001     | 892.3
##       BRCA1    | TCGA_BRCA::TCGA-A1-A0SB  | 234.1
##       ...      | ...                      | ...
##
##   - Genes are human orthologs (converted via pgx$genes$human_ortholog)
##   - Samples are prefixed with dataset name: "dataset::original_sample_id"
##   - Zeros ARE stored (distinguish from missing data which is NA)
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
##   Dataset:    pgx.getDatasetFromSample()  -> "dataset_name"
##   Short name: pgx.getSampleShortName()    -> "original_sample_id"
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
pgx.buildTileDB <- function(pgx_folder, tiledb_path, overwrite = FALSE, verbose = TRUE) {
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
  schema <- tiledb::tiledb_array_schema(domain = dom, attrs = list(count_attr), sparse = TRUE)
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

      counts <- as.matrix(pgx$counts)
      counts <- rename_by2(counts, pgx$genes, new_id = "human_ortholog", keep.prefix = FALSE)

      genes <- rownames(counts)
      samples <- colnames(counts)
      file_prefix <- tools::file_path_sans_ext(basename(pgx_file))

      for (j in seq_along(samples)) {
        sample_name <- paste0(file_prefix, "::", samples[j])
        values <- counts[, j]
        valid_idx <- which(!is.na(values))

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

      if (!is.null(pgx$samples) && is.data.frame(pgx$samples)) {
        pheno_data <- pgx$samples
        pheno_cols <- colnames(pheno_data)
        pheno_cols <- pheno_cols[!grepl("^\\.", pheno_cols)]

        if (length(pheno_cols) > 0) {
          for (j in seq_len(nrow(pheno_data))) {
            sample_id <- rownames(pheno_data)[j]
            sample_name <- paste0(file_prefix, "::", sample_id)

            for (col in pheno_cols) {
              val <- pheno_data[j, col]
              val_str <- if (is.na(val)) NA_character_ else as.character(val)
              all_phenotype_data[[length(all_phenotype_data) + 1]] <- data.frame(
                sample = sample_name,
                phenotype = col,
                value = val_str,
                stringsAsFactors = FALSE
              )
            }
          }
        }
      }

      rm(pgx, counts)
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


#' @title Query gene counts from TileDB database
#'
#' @param tiledb_path Path to the TileDB database
#' @param genes Character vector of gene names to query
#' @param samples Optional character vector of sample names to filter
#' @param as_matrix Return as matrix (TRUE) or data.frame (FALSE)
#'
#' @return Matrix (genes x samples) or data.frame with columns: gene, sample, count
#'
#' @export
pgx.queryTileDB <- function(tiledb_path, genes, samples = NULL, as_matrix = TRUE) {
  if (!requireNamespace("tiledb", quietly = TRUE)) {
    stop("Package 'tiledb' is required. Install with: install.packages('tiledb')")
  }

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

  arr <- tiledb::tiledb_array(tiledb_path, as.data.frame = TRUE)

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
      return(data.frame(gene = character(0), sample = character(0), count = numeric(0)))
    }
  }

  if (!as_matrix) return(df)

  sample_names <- if (!is.null(samples)) samples else if (!is.null(all_samples)) all_samples else unique(df$sample)

  mat <- matrix(NA_real_, nrow = length(genes), ncol = length(sample_names),
                dimnames = list(genes, sample_names))

  for (i in seq_len(nrow(df))) {
    g <- df$gene[i]
    s <- df$sample[i]
    if (g %in% genes && s %in% sample_names) {
      mat[g, s] <- df$count[i]
    }
  }

  mat
}


#' @title List genes in TileDB database
#' @param tiledb_path Path to the TileDB database
#' @return Character vector of gene names
#' @export
pgx.listGenesTileDB <- function(tiledb_path) {
  metadata <- readRDS(paste0(tiledb_path, "_metadata.rds"))
  metadata$genes
}


#' @title List samples in TileDB database
#' @param tiledb_path Path to the TileDB database
#' @return Character vector of sample names (format: "filename::samplename")
#' @export
pgx.listSamplesTileDB <- function(tiledb_path) {
  metadata <- readRDS(paste0(tiledb_path, "_metadata.rds"))
  metadata$samples
}


#' @title List datasets in TileDB database
#' @param tiledb_path Path to the TileDB database
#' @return Character vector of dataset names
#' @export
pgx.listDatasetsTileDB <- function(tiledb_path) {
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
pgx.getDatasetInfoTileDB <- function(tiledb_path, datasets_info_file = NULL) {
  metadata <- readRDS(paste0(tiledb_path, "_metadata.rds"))

  datasets <- tools::file_path_sans_ext(basename(metadata$pgx_files))

  ## Count samples per dataset
  sample_datasets <- pgx.getDatasetFromSample(metadata$samples)
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
pgx.listDatasetsTileDB <- function(tiledb_path) {
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
pgx.infoTileDB <- function(tiledb_path) {
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
pgx.listPhenotypesTileDB <- function(tiledb_path) {
  metadata <- readRDS(paste0(tiledb_path, "_metadata.rds"))
  if (is.null(metadata$phenotypes)) return(character(0))
  metadata$phenotypes
}


#' @title List phenotype columns by dataset
#' @param tiledb_path Path to the TileDB database
#' @param dataset Optional dataset name. If NULL, returns list for all datasets.
#' @return Named list of phenotype columns per dataset, or character vector if dataset specified
#' @export
pgx.listPhenotypesByDatasetTileDB <- function(tiledb_path, dataset = NULL) {
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
pgx.queryPhenotypesTileDB <- function(tiledb_path, samples = NULL, phenotypes = NULL,
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
    sample_datasets <- pgx.getDatasetFromSample(all_samples)
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

  wide_df$dataset <- pgx.getDatasetFromSample(rownames(wide_df))
  wide_df$sample_short <- pgx.getSampleShortName(rownames(wide_df))
  wide_df
}


#' @title Get unique values for a phenotype column
#' @param tiledb_path Path to the TileDB database
#' @param phenotype Name of the phenotype column
#' @param datasets Optional dataset names to filter by
#' @return Character vector of unique values for the phenotype
#' @export
pgx.getPhenotypeValuesTileDB <- function(tiledb_path, phenotype, datasets = NULL) {
  metadata <- readRDS(paste0(tiledb_path, "_metadata.rds"))

  if (is.null(metadata$phenotype_data) || nrow(metadata$phenotype_data) == 0) {
    return(character(0))
  }

  result <- metadata$phenotype_data[metadata$phenotype_data$phenotype == phenotype, , drop = FALSE]
  if (nrow(result) == 0) return(character(0))

  if (!is.null(datasets)) {
    sample_datasets <- pgx.getDatasetFromSample(result$sample)
    result <- result[sample_datasets %in% datasets, , drop = FALSE]
  }

  values <- unique(result$value)
  sort(values[!is.na(values)])
}


#' @title Extract dataset name from sample identifiers
#' @param samples Character vector of sample identifiers ("dataset::sample")
#' @return Character vector of dataset names
#' @export
pgx.getDatasetFromSample <- function(samples) {
  sub("::.*", "", samples)
}


#' @title Extract short sample name from sample identifiers
#' @param samples Character vector of sample identifiers ("dataset::sample")
#' @return Character vector of short sample names
#' @export
pgx.getSampleShortName <- function(samples) {
  sub(".*::", "", samples)
}


#' @title Convert TileDB query result to plotting data.frame
#'
#' @param result Matrix or data.frame from pgx.queryTileDB
#' @param gene_name Optional gene name for single-gene matrix input
#' @param tiledb_path Optional path to TileDB database for merging phenotype data
#' @param phenotypes Optional phenotype columns to include
#'
#' @return Data.frame with columns: sample, count, dataset, sample_short, gene, phenotypes
#'
#' @export
pgx.tiledbToPlotDF <- function(result, gene_name = NULL, tiledb_path = NULL, phenotypes = NULL) {
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

  df$dataset <- pgx.getDatasetFromSample(df$sample)
  df$sample_short <- pgx.getSampleShortName(df$sample)

  if (!is.null(tiledb_path) && !is.null(phenotypes)) {
    pheno_df <- pgx.queryPhenotypesTileDB(
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

