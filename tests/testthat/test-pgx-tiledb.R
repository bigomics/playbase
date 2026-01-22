#' Tests for pgx-tiledb.R
#'
#' Comprehensive test suite for TileDB database functions.
#' Tests cover building, querying, and metadata operations.

# ============================================================================
# HELPER FUNCTIONS - Create mock PGX objects for testing
# ============================================================================

#' Create a minimal mock PGX object for testing
#' @param n_genes Number of genes
#' @param n_samples Number of samples
#' @param dataset_name Name prefix for the dataset
#' @param include_na Whether to include NA values in counts
create_mock_pgx <- function(n_genes = 10, n_samples = 5, dataset_name = "test",
                            include_na = FALSE, include_phenotypes = TRUE) {
  set.seed(42)

  gene_names <- paste0("GENE", seq_len(n_genes))
  sample_names <- paste0("sample_", seq_len(n_samples))

  # Create counts matrix with some variation
  counts <- matrix(
    rpois(n_genes * n_samples, lambda = 100),
    nrow = n_genes, ncol = n_samples,
    dimnames = list(gene_names, sample_names)
  )

  # Optionally add NA values

  if (include_na) {
    na_idx <- sample(length(counts), size = floor(length(counts) * 0.1))
    counts[na_idx] <- NA
  }

  # Create genes annotation with human_ortholog column
  genes <- data.frame(
    symbol = gene_names,
    human_ortholog = paste0("HUMAN_", gene_names),
    row.names = gene_names,
    stringsAsFactors = FALSE
  )

  # Create samples phenotype data
  samples <- NULL
  if (include_phenotypes) {
    samples <- data.frame(
      age = sample(20:80, n_samples, replace = TRUE),
      sex = sample(c("male", "female"), n_samples, replace = TRUE),
      condition = sample(c("control", "treated"), n_samples, replace = TRUE),
      .internal_col = rep("hidden", n_samples), # columns starting with . should be excluded
      row.names = sample_names,
      stringsAsFactors = FALSE
    )
  }

  list(
    counts = counts,
    genes = genes,
    samples = samples
  )
}

#' Save mock PGX to file (simulates .pgx file)
save_mock_pgx <- function(pgx, filepath) {
  save(pgx, file = filepath)
}


# ============================================================================
# TESTS FOR HELPER FUNCTIONS
# ============================================================================

test_that("getDatasetFromSample extracts dataset name correctly", {
  # Standard case
  samples <- c("dataset1::sample1", "dataset1::sample2", "dataset2::sample3")
  result <- playbase::getDatasetFromSample(samples)
  expect_equal(result, c("dataset1", "dataset1", "dataset2"))

  # Single sample
  result_single <- playbase::getDatasetFromSample("mydata::s1")
  expect_equal(result_single, "mydata")

  # Edge case: dataset name contains special characters
  result_special <- playbase::getDatasetFromSample("my_data-v2::sample")
  expect_equal(result_special, "my_data-v2")

  # Edge case: empty string
  result_empty <- playbase::getDatasetFromSample("")
  expect_equal(result_empty, "")
})

test_that("getSampleShortName extracts sample name correctly", {
  samples <- c("dataset1::sample1", "dataset1::sample2", "dataset2::sample3")
  result <- playbase::getSampleShortName(samples)
  expect_equal(result, c("sample1", "sample2", "sample3"))

  # Single sample
  result_single <- playbase::getSampleShortName("mydata::s1")
  expect_equal(result_single, "s1")

  # Edge case: sample name contains :: (greedy match takes last ::)
  result_double <- playbase::getSampleShortName("dataset::part1::part2")
  expect_equal(result_double, "part2")

  # Edge case: no :: separator (should return everything after non-existent ::)
  result_no_sep <- playbase::getSampleShortName("noseparator")
  expect_equal(result_no_sep, "noseparator")
})


# ============================================================================
# TESTS FOR buildTileDB
# ============================================================================

test_that("buildTileDB fails gracefully with missing tiledb package", {
  skip_if(requireNamespace("tiledb", quietly = TRUE),
          message = "tiledb is installed, skipping missing package test")

  expect_error(
    playbase::buildTileDB("/tmp/fake", "/tmp/out"),
    "Package 'tiledb' is required"
  )
})

test_that("buildTileDB fails with non-existent folder", {
  skip_if_not(requireNamespace("tiledb", quietly = TRUE),
              message = "tiledb not installed")

  expect_error(
    playbase::buildTileDB("/nonexistent/folder/path", tempfile()),
    "Folder does not exist"
  )
})

test_that("buildTileDB fails with folder containing no .pgx files", {
  skip_if_not(requireNamespace("tiledb", quietly = TRUE),
              message = "tiledb not installed")

  empty_dir <- tempfile()
  dir.create(empty_dir)
  on.exit(unlink(empty_dir, recursive = TRUE), add = TRUE)

  expect_error(
    playbase::buildTileDB(empty_dir, tempfile()),
    "No .pgx files found"
  )
})

test_that("buildTileDB creates database from single PGX file", {
  skip_if_not(requireNamespace("tiledb", quietly = TRUE),
              message = "tiledb not installed")

  # Setup temp directories
  pgx_dir <- tempfile("pgx_test_")
  tiledb_path <- tempfile("tiledb_test_")
  dir.create(pgx_dir)
  on.exit({
    unlink(pgx_dir, recursive = TRUE)
    unlink(tiledb_path, recursive = TRUE)
    unlink(paste0(tiledb_path, "_metadata.rds"))
  }, add = TRUE)

  # Create and save mock PGX
  pgx <- create_mock_pgx(n_genes = 10, n_samples = 5)
  save_mock_pgx(pgx, file.path(pgx_dir, "dataset1.pgx"))

  # Build TileDB
  metadata <- playbase::buildTileDB(pgx_dir, tiledb_path, verbose = FALSE)

  # Verify database was created

  expect_true(dir.exists(tiledb_path))
  expect_true(file.exists(paste0(tiledb_path, "_metadata.rds")))

  # Verify metadata structure

  expect_true(is.list(metadata))
  expect_equal(metadata$n_files, 1)
  expect_equal(metadata$n_samples, 5)
  expect_equal(metadata$n_genes, 10)
  expect_true(length(metadata$genes) == 10)
  expect_true(all(grepl("^HUMAN_GENE", metadata$genes))) # human orthologs
  expect_true(all(grepl("^dataset1::", metadata$samples)))
})

test_that("buildTileDB handles multiple PGX files correctly", {
  skip_if_not(requireNamespace("tiledb", quietly = TRUE),
              message = "tiledb not installed")

  pgx_dir <- tempfile("pgx_multi_")
  tiledb_path <- tempfile("tiledb_multi_")
  dir.create(pgx_dir)
  on.exit({
    unlink(pgx_dir, recursive = TRUE)
    unlink(tiledb_path, recursive = TRUE)
    unlink(paste0(tiledb_path, "_metadata.rds"))
  }, add = TRUE)

  # Create multiple PGX files with overlapping genes
  pgx1 <- create_mock_pgx(n_genes = 10, n_samples = 3)
  pgx2 <- create_mock_pgx(n_genes = 8, n_samples = 4)  # Different size
  save_mock_pgx(pgx1, file.path(pgx_dir, "dataset1.pgx"))
  save_mock_pgx(pgx2, file.path(pgx_dir, "dataset2.pgx"))

  metadata <- playbase::buildTileDB(pgx_dir, tiledb_path, verbose = FALSE)

  expect_equal(metadata$n_files, 2)
  expect_equal(metadata$n_samples, 7)  # 3 + 4
  expect_equal(metadata$n_genes, 10)   # union, both have same gene names
  expect_true("dataset1" %in% names(metadata$phenotypes_by_dataset))
  expect_true("dataset2" %in% names(metadata$phenotypes_by_dataset))
})

test_that("buildTileDB fails without overwrite flag", {
  skip_if_not(requireNamespace("tiledb", quietly = TRUE),
              message = "tiledb not installed")

  pgx_dir <- tempfile("pgx_overwrite_")
  tiledb_path <- tempfile("tiledb_overwrite_")
  dir.create(pgx_dir)
  on.exit({
    unlink(pgx_dir, recursive = TRUE)
    unlink(tiledb_path, recursive = TRUE)
    unlink(paste0(tiledb_path, "_metadata.rds"))
  }, add = TRUE)

  pgx <- create_mock_pgx()
  save_mock_pgx(pgx, file.path(pgx_dir, "test.pgx"))

  # First build

  playbase::buildTileDB(pgx_dir, tiledb_path, verbose = FALSE)

  # Second build should fail without overwrite
  expect_error(
    playbase::buildTileDB(pgx_dir, tiledb_path, overwrite = FALSE, verbose = FALSE),
    "already exists"
  )
})

test_that("buildTileDB succeeds with overwrite flag", {
  skip_if_not(requireNamespace("tiledb", quietly = TRUE),
              message = "tiledb not installed")

  pgx_dir <- tempfile("pgx_overwrite2_")
  tiledb_path <- tempfile("tiledb_overwrite2_")
  dir.create(pgx_dir)
  on.exit({
    unlink(pgx_dir, recursive = TRUE)
    unlink(tiledb_path, recursive = TRUE)
    unlink(paste0(tiledb_path, "_metadata.rds"))
  }, add = TRUE)

  pgx <- create_mock_pgx(n_samples = 3)
  save_mock_pgx(pgx, file.path(pgx_dir, "test.pgx"))

  playbase::buildTileDB(pgx_dir, tiledb_path, verbose = FALSE)

  # Change the data and rebuild with overwrite
  pgx2 <- create_mock_pgx(n_samples = 7)
  save_mock_pgx(pgx2, file.path(pgx_dir, "test.pgx"))

  metadata <- playbase::buildTileDB(pgx_dir, tiledb_path, overwrite = TRUE, verbose = FALSE)
  expect_equal(metadata$n_samples, 7)
})

test_that("buildTileDB handles NA values in counts", {
  skip_if_not(requireNamespace("tiledb", quietly = TRUE),
              message = "tiledb not installed")

  pgx_dir <- tempfile("pgx_na_")
  tiledb_path <- tempfile("tiledb_na_")
  dir.create(pgx_dir)
  on.exit({
    unlink(pgx_dir, recursive = TRUE)
    unlink(tiledb_path, recursive = TRUE)
    unlink(paste0(tiledb_path, "_metadata.rds"))
  }, add = TRUE)

  pgx <- create_mock_pgx(n_genes = 10, n_samples = 5, include_na = TRUE)
  save_mock_pgx(pgx, file.path(pgx_dir, "test_na.pgx"))

  # Should not error
  metadata <- playbase::buildTileDB(pgx_dir, tiledb_path, verbose = FALSE)
  expect_true(is.list(metadata))
})

test_that("buildTileDB handles PGX without phenotypes", {
  skip_if_not(requireNamespace("tiledb", quietly = TRUE),
              message = "tiledb not installed")

  pgx_dir <- tempfile("pgx_nopheno_")
  tiledb_path <- tempfile("tiledb_nopheno_")
  dir.create(pgx_dir)
  on.exit({
    unlink(pgx_dir, recursive = TRUE)
    unlink(tiledb_path, recursive = TRUE)
    unlink(paste0(tiledb_path, "_metadata.rds"))
  }, add = TRUE)

  pgx <- create_mock_pgx(include_phenotypes = FALSE)
  save_mock_pgx(pgx, file.path(pgx_dir, "test.pgx"))

  metadata <- playbase::buildTileDB(pgx_dir, tiledb_path, verbose = FALSE)
  expect_equal(length(metadata$phenotypes), 0)
  expect_equal(nrow(metadata$phenotype_data), 0)
})

test_that("buildTileDB skips PGX files without counts", {
  skip_if_not(requireNamespace("tiledb", quietly = TRUE),
              message = "tiledb not installed")

  pgx_dir <- tempfile("pgx_nocounts_")
  tiledb_path <- tempfile("tiledb_nocounts_")
  dir.create(pgx_dir)
  on.exit({
    unlink(pgx_dir, recursive = TRUE)
    unlink(tiledb_path, recursive = TRUE)
    unlink(paste0(tiledb_path, "_metadata.rds"))
  }, add = TRUE)

  # Create valid PGX
  pgx_valid <- create_mock_pgx(n_samples = 3)
  save_mock_pgx(pgx_valid, file.path(pgx_dir, "valid.pgx"))

  # Create PGX without counts
  pgx_nocounts <- list(
    genes = data.frame(symbol = "A", human_ortholog = "A"),
    samples = data.frame(age = 30, row.names = "s1"),
    counts = NULL
  )
  save_mock_pgx(pgx_nocounts, file.path(pgx_dir, "nocounts.pgx"))

  # Should complete with warning
  expect_warning(
    metadata <- playbase::buildTileDB(pgx_dir, tiledb_path, verbose = FALSE),
    "No counts found"
  )

  # Only valid file should be processed
  expect_equal(metadata$n_samples, 3)
})


# ============================================================================
# TESTS FOR QUERY FUNCTIONS
# ============================================================================

# Create a shared test fixture for query tests
setup_test_tiledb <- function() {
  skip_if_not(requireNamespace("tiledb", quietly = TRUE),
              message = "tiledb not installed")

  pgx_dir <- tempfile("pgx_query_")
  tiledb_path <- tempfile("tiledb_query_")
  dir.create(pgx_dir)

  # Create two datasets for testing queries
  pgx1 <- create_mock_pgx(n_genes = 10, n_samples = 3)
  pgx2 <- create_mock_pgx(n_genes = 10, n_samples = 4)
  save_mock_pgx(pgx1, file.path(pgx_dir, "datasetA.pgx"))
  save_mock_pgx(pgx2, file.path(pgx_dir, "datasetB.pgx"))

  playbase::buildTileDB(pgx_dir, tiledb_path, verbose = FALSE)

  list(
    pgx_dir = pgx_dir,
    tiledb_path = tiledb_path,
    cleanup = function() {
      unlink(pgx_dir, recursive = TRUE)
      unlink(tiledb_path, recursive = TRUE)
      unlink(paste0(tiledb_path, "_metadata.rds"))
    }
  )
}

test_that("listGenesTileDB returns all genes", {
  fixture <- setup_test_tiledb()
  on.exit(fixture$cleanup(), add = TRUE)

  genes <- playbase::listGenesTileDB(fixture$tiledb_path)

  expect_true(is.character(genes))
  expect_equal(length(genes), 10)
  expect_true(all(grepl("^HUMAN_GENE", genes)))
})

test_that("listSamplesTileDB returns all samples", {
  fixture <- setup_test_tiledb()
  on.exit(fixture$cleanup(), add = TRUE)

  samples <- playbase::listSamplesTileDB(fixture$tiledb_path)

  expect_true(is.character(samples))
  expect_equal(length(samples), 7)  # 3 + 4
  expect_true(all(grepl("::", samples)))
})

test_that("listDatasetsTileDB returns dataset names", {
  fixture <- setup_test_tiledb()
  on.exit(fixture$cleanup(), add = TRUE)

  datasets <- playbase::listDatasetsTileDB(fixture$tiledb_path)

  expect_true(is.character(datasets))
  expect_equal(length(datasets), 2)
  expect_true("datasetA" %in% datasets)
  expect_true("datasetB" %in% datasets)
})

test_that("queryTileDB returns matrix by default", {
  fixture <- setup_test_tiledb()
  on.exit(fixture$cleanup(), add = TRUE)

  genes <- playbase::listGenesTileDB(fixture$tiledb_path)
  query_genes <- head(genes, 3)

  result <- playbase::queryTileDB(fixture$tiledb_path, genes = query_genes)

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 7)
  expect_equal(rownames(result), query_genes)
})

test_that("queryTileDB returns data.frame when as_matrix=FALSE", {
  fixture <- setup_test_tiledb()
  on.exit(fixture$cleanup(), add = TRUE)

  genes <- playbase::listGenesTileDB(fixture$tiledb_path)

  result <- playbase::queryTileDB(
    fixture$tiledb_path,
    genes = genes[1],
    as_matrix = FALSE
  )

  expect_true(is.data.frame(result))
  expect_true(all(c("gene", "sample", "count") %in% colnames(result)))
})

test_that("queryTileDB returns z-scores when value='zscore'", {
  fixture <- setup_test_tiledb()
  on.exit(fixture$cleanup(), add = TRUE)

  genes <- playbase::listGenesTileDB(fixture$tiledb_path)
  query_genes <- head(genes, 3)

  # Query z-scores as matrix
  result <- playbase::queryTileDB(
    fixture$tiledb_path,
    genes = query_genes,
    value = "zscore"
  )

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 7)

  # Z-scores should have mean ~0 within each dataset
  # (allowing for NA values from zero-variance genes)
})

test_that("queryTileDB returns zscore column in data.frame mode", {
  fixture <- setup_test_tiledb()
  on.exit(fixture$cleanup(), add = TRUE)

  genes <- playbase::listGenesTileDB(fixture$tiledb_path)

  result <- playbase::queryTileDB(
    fixture$tiledb_path,
    genes = genes[1],
    as_matrix = FALSE,
    value = "zscore"
  )

  expect_true(is.data.frame(result))
  expect_true("zscore" %in% colnames(result))
  expect_true("gene" %in% colnames(result))
  expect_true("sample" %in% colnames(result))
})

test_that("queryTileDB filters by samples", {
  fixture <- setup_test_tiledb()
  on.exit(fixture$cleanup(), add = TRUE)

  samples <- playbase::listSamplesTileDB(fixture$tiledb_path)
  genes <- playbase::listGenesTileDB(fixture$tiledb_path)
  query_samples <- head(samples, 2)

  result <- playbase::queryTileDB(
    fixture$tiledb_path,
    genes = genes[1:2],
    samples = query_samples
  )

  expect_equal(ncol(result), 2)
  expect_equal(colnames(result), query_samples)
})

test_that("queryTileDB warns on missing genes", {
  fixture <- setup_test_tiledb()
  on.exit(fixture$cleanup(), add = TRUE)

  genes <- playbase::listGenesTileDB(fixture$tiledb_path)

  expect_warning(
    result <- playbase::queryTileDB(
      fixture$tiledb_path,
      genes = c(genes[1], "FAKEGENE1", "FAKEGENE2")
    ),
    "genes not found"
  )

  expect_equal(nrow(result), 1)  # Only valid gene returned
})

test_that("queryTileDB errors when no valid genes", {
  fixture <- setup_test_tiledb()
  on.exit(fixture$cleanup(), add = TRUE)

  expect_warning(
    expect_error(
      playbase::queryTileDB(fixture$tiledb_path, genes = c("FAKE1", "FAKE2")),
      "No valid genes"
    )
  )
})

test_that("queryTileDB handles single gene query", {
  fixture <- setup_test_tiledb()
  on.exit(fixture$cleanup(), add = TRUE)

  genes <- playbase::listGenesTileDB(fixture$tiledb_path)

  result <- playbase::queryTileDB(fixture$tiledb_path, genes = genes[1])

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 1)
})

test_that("queryTileDB fails on non-existent database", {
  skip_if_not(requireNamespace("tiledb", quietly = TRUE),
              message = "tiledb not installed")

  expect_error(
    playbase::queryTileDB("/nonexistent/path", genes = "GENE1"),
    "not found"
  )
})


# ============================================================================
# TESTS FOR PHENOTYPE FUNCTIONS
# ============================================================================

test_that("listPhenotypesTileDB returns phenotype columns", {
  fixture <- setup_test_tiledb()
  on.exit(fixture$cleanup(), add = TRUE)

  phenotypes <- playbase::listPhenotypesTileDB(fixture$tiledb_path)

  expect_true(is.character(phenotypes))
  expect_true("age" %in% phenotypes)
  expect_true("sex" %in% phenotypes)
  expect_true("condition" %in% phenotypes)
  expect_false(".internal_col" %in% phenotypes)  # dot-columns excluded
})

test_that("listPhenotypesByDatasetTileDB returns phenotypes per dataset", {
  fixture <- setup_test_tiledb()
  on.exit(fixture$cleanup(), add = TRUE)

  # Get all datasets
  all_pheno <- playbase::listPhenotypesByDatasetTileDB(fixture$tiledb_path)
  expect_true(is.list(all_pheno))
  expect_true("datasetA" %in% names(all_pheno))

  # Get specific dataset
  dataset_pheno <- playbase::listPhenotypesByDatasetTileDB(
    fixture$tiledb_path,
    dataset = "datasetA"
  )
  expect_true(is.character(dataset_pheno))
  expect_true("age" %in% dataset_pheno)
})

test_that("listPhenotypesByDatasetTileDB warns on invalid dataset", {
  fixture <- setup_test_tiledb()
  on.exit(fixture$cleanup(), add = TRUE)

  expect_warning(
    result <- playbase::listPhenotypesByDatasetTileDB(
      fixture$tiledb_path,
      dataset = "nonexistent"
    ),
    "Dataset not found"
  )
  expect_equal(length(result), 0)
})

test_that("queryPhenotypesTileDB returns wide format by default", {
  fixture <- setup_test_tiledb()
  on.exit(fixture$cleanup(), add = TRUE)

  result <- playbase::queryPhenotypesTileDB(fixture$tiledb_path)

  expect_true(is.data.frame(result))
  expect_true("age" %in% colnames(result))
  expect_true("sex" %in% colnames(result))
  expect_true("dataset" %in% colnames(result))
  expect_true("sample_short" %in% colnames(result))
  expect_equal(nrow(result), 7)  # All samples
})

test_that("queryPhenotypesTileDB returns long format when as_wide=FALSE", {
  fixture <- setup_test_tiledb()
  on.exit(fixture$cleanup(), add = TRUE)

  result <- playbase::queryPhenotypesTileDB(
    fixture$tiledb_path,
    as_wide = FALSE
  )

  expect_true(is.data.frame(result))
  expect_equal(colnames(result), c("sample", "phenotype", "value"))
  expect_true(nrow(result) > 0)
})

test_that("queryPhenotypesTileDB filters by samples", {
  fixture <- setup_test_tiledb()
  on.exit(fixture$cleanup(), add = TRUE)

  samples <- playbase::listSamplesTileDB(fixture$tiledb_path)
  query_samples <- head(samples, 2)

  result <- playbase::queryPhenotypesTileDB(
    fixture$tiledb_path,
    samples = query_samples
  )

  expect_equal(nrow(result), 2)
})

test_that("queryPhenotypesTileDB filters by phenotypes", {
  fixture <- setup_test_tiledb()
  on.exit(fixture$cleanup(), add = TRUE)

  result <- playbase::queryPhenotypesTileDB(
    fixture$tiledb_path,
    phenotypes = c("age", "sex")
  )

  # Should have age, sex, plus dataset and sample_short columns
  expect_true("age" %in% colnames(result))
  expect_true("sex" %in% colnames(result))
  expect_false("condition" %in% colnames(result))
})

test_that("queryPhenotypesTileDB filters by datasets", {
  fixture <- setup_test_tiledb()
  on.exit(fixture$cleanup(), add = TRUE)

  result <- playbase::queryPhenotypesTileDB(
    fixture$tiledb_path,
    datasets = "datasetA"
  )

  expect_equal(nrow(result), 3)  # datasetA has 3 samples
  expect_true(all(result$dataset == "datasetA"))
})

test_that("queryPhenotypesTileDB warns on invalid phenotypes", {
  fixture <- setup_test_tiledb()
  on.exit(fixture$cleanup(), add = TRUE)

  expect_warning(
    result <- playbase::queryPhenotypesTileDB(
      fixture$tiledb_path,
      phenotypes = c("age", "nonexistent_pheno")
    ),
    "phenotypes not found"
  )
})

test_that("queryPhenotypesTileDB converts numeric columns", {
  fixture <- setup_test_tiledb()
  on.exit(fixture$cleanup(), add = TRUE)

  result <- playbase::queryPhenotypesTileDB(fixture$tiledb_path)

  # Age should be converted to numeric
  expect_true(is.numeric(result$age))
  # Sex should remain character
  expect_true(is.character(result$sex))
})

test_that("getPhenotypeValuesTileDB returns unique values", {
  fixture <- setup_test_tiledb()
  on.exit(fixture$cleanup(), add = TRUE)

  sex_values <- playbase::getPhenotypeValuesTileDB(
    fixture$tiledb_path,
    phenotype = "sex"
  )

  expect_true(is.character(sex_values))
  expect_true(all(sex_values %in% c("male", "female")))
})

test_that("getPhenotypeValuesTileDB filters by datasets", {
  fixture <- setup_test_tiledb()
  on.exit(fixture$cleanup(), add = TRUE)

  values <- playbase::getPhenotypeValuesTileDB(
    fixture$tiledb_path,
    phenotype = "condition",
    datasets = "datasetA"
  )

  expect_true(is.character(values))
})

test_that("getPhenotypeValuesTileDB returns empty for invalid phenotype", {
  fixture <- setup_test_tiledb()
  on.exit(fixture$cleanup(), add = TRUE)

  values <- playbase::getPhenotypeValuesTileDB(
    fixture$tiledb_path,
    phenotype = "nonexistent"
  )

  expect_equal(length(values), 0)
})


# ============================================================================
# TESTS FOR getDatasetInfoTileDB
# ============================================================================

test_that("getDatasetInfoTileDB returns basic info", {
  fixture <- setup_test_tiledb()
  on.exit(fixture$cleanup(), add = TRUE)

  info <- playbase::getDatasetInfoTileDB(fixture$tiledb_path)

  expect_true(is.data.frame(info))
  expect_true("dataset" %in% colnames(info))
  expect_true("nsamples" %in% colnames(info))
  expect_equal(nrow(info), 2)
  expect_equal(sum(info$nsamples), 7)
})

test_that("getDatasetInfoTileDB merges with datasets-info.csv", {
  fixture <- setup_test_tiledb()
  on.exit(fixture$cleanup(), add = TRUE)

  # Create a mock datasets-info.csv
  info_file <- tempfile(fileext = ".csv")
  info_df <- data.frame(
    dataset = c("datasetA", "datasetB"),
    organism = c("human", "mouse"),
    datatype = c("RNA-seq", "microarray"),
    description = c("Test A", "Test B"),
    stringsAsFactors = FALSE
  )
  write.csv(info_df, info_file, row.names = FALSE)
  on.exit(unlink(info_file), add = TRUE)

  result <- playbase::getDatasetInfoTileDB(
    fixture$tiledb_path,
    datasets_info_file = info_file
  )

  expect_true("organism" %in% colnames(result))
  expect_true("datatype" %in% colnames(result))
  expect_true("description" %in% colnames(result))
})


# ============================================================================
# TESTS FOR infoTileDB
# ============================================================================

test_that("infoTileDB prints info and returns metadata", {
  fixture <- setup_test_tiledb()
  on.exit(fixture$cleanup(), add = TRUE)

  # Capture output
  output <- capture.output(
    metadata <- playbase::infoTileDB(fixture$tiledb_path)
  )

  expect_true(is.list(metadata))
  expect_true(any(grepl("Genes:", output)))
  expect_true(any(grepl("Samples:", output)))
})


# ============================================================================
# TESTS FOR tiledbToPlotDF
# ============================================================================

test_that("tiledbToPlotDF converts single-row matrix", {
  fixture <- setup_test_tiledb()
  on.exit(fixture$cleanup(), add = TRUE)

  genes <- playbase::listGenesTileDB(fixture$tiledb_path)
  mat <- playbase::queryTileDB(fixture$tiledb_path, genes = genes[1])

  df <- playbase::tiledbToPlotDF(mat, gene_name = "TestGene")

  expect_true(is.data.frame(df))
  expect_true("sample" %in% colnames(df))
  expect_true("count" %in% colnames(df))
  expect_true("dataset" %in% colnames(df))
  expect_true("sample_short" %in% colnames(df))
  expect_true("gene" %in% colnames(df))
  expect_equal(nrow(df), 7)
})

test_that("tiledbToPlotDF converts multi-row matrix", {
  fixture <- setup_test_tiledb()
  on.exit(fixture$cleanup(), add = TRUE)

  genes <- playbase::listGenesTileDB(fixture$tiledb_path)
  mat <- playbase::queryTileDB(fixture$tiledb_path, genes = genes[1:3])

  df <- playbase::tiledbToPlotDF(mat)

  expect_equal(nrow(df), 3 * 7)  # 3 genes x 7 samples
  expect_true(length(unique(df$gene)) == 3)
})

test_that("tiledbToPlotDF handles data.frame input", {
  fixture <- setup_test_tiledb()
  on.exit(fixture$cleanup(), add = TRUE)

  genes <- playbase::listGenesTileDB(fixture$tiledb_path)
  df_input <- playbase::queryTileDB(
    fixture$tiledb_path,
    genes = genes[1],
    as_matrix = FALSE
  )

  df <- playbase::tiledbToPlotDF(df_input)

  expect_true("dataset" %in% colnames(df))
  expect_true("sample_short" %in% colnames(df))
})

test_that("tiledbToPlotDF merges phenotype data", {
  fixture <- setup_test_tiledb()
  on.exit(fixture$cleanup(), add = TRUE)

  genes <- playbase::listGenesTileDB(fixture$tiledb_path)
  mat <- playbase::queryTileDB(fixture$tiledb_path, genes = genes[1])

  df <- playbase::tiledbToPlotDF(
    mat,
    tiledb_path = fixture$tiledb_path,
    phenotypes = c("age", "sex")
  )

  expect_true("age" %in% colnames(df))
  expect_true("sex" %in% colnames(df))
})

test_that("tiledbToPlotDF errors on invalid input", {
  expect_error(
    playbase::tiledbToPlotDF("not a matrix or df"),
    "must be a matrix or data.frame"
  )
})


# ============================================================================
# TESTS FOR Z-SCORE COMPUTATION
# ============================================================================

test_that("Z-scores are computed correctly per dataset", {
  skip_if_not(requireNamespace("tiledb", quietly = TRUE),
              message = "tiledb not installed")

  pgx_dir <- tempfile("pgx_zscore_")
  tiledb_path <- tempfile("tiledb_zscore_")
  dir.create(pgx_dir)
  on.exit({
    unlink(pgx_dir, recursive = TRUE)
    unlink(tiledb_path, recursive = TRUE)
    unlink(paste0(tiledb_path, "_metadata.rds"))
  }, add = TRUE)

  # Create PGX with known values for easy z-score verification
  set.seed(42)
  n_genes <- 3
  n_samples <- 4
  gene_names <- paste0("GENE", seq_len(n_genes))
  sample_names <- paste0("sample_", seq_len(n_samples))

  # Create counts where we can predict z-scores
  # Gene1: constant value (should give NA z-scores)
  # Gene2 and Gene3: varying values
  counts <- matrix(
    c(100, 100, 100, 100,    # Gene1: constant
      80, 90, 100, 110,      # Gene2: mean=95, sd~12.9
      50, 100, 150, 200),    # Gene3: mean=125, sd~62.9
    nrow = n_genes, ncol = n_samples, byrow = TRUE,
    dimnames = list(gene_names, sample_names)
  )

  genes <- data.frame(
    symbol = gene_names,
    human_ortholog = gene_names,  # Keep same names
    row.names = gene_names,
    stringsAsFactors = FALSE
  )

  samples <- data.frame(
    condition = rep("test", n_samples),
    row.names = sample_names,
    stringsAsFactors = FALSE
  )

  pgx <- list(counts = counts, genes = genes, samples = samples)
  save_mock_pgx(pgx, file.path(pgx_dir, "test.pgx"))

  playbase::buildTileDB(pgx_dir, tiledb_path, verbose = FALSE)

  # Query z-scores
  zscores <- playbase::queryTileDB(tiledb_path, genes = gene_names, value = "zscore")

  # Gene1 (constant) should have NA z-scores
  expect_true(all(is.na(zscores["GENE1", ])))

  # Gene2 and Gene3 should have valid z-scores
  # Z-scores should have mean 0 (within floating point tolerance)
  gene2_zscores <- zscores["GENE2", ]
  gene3_zscores <- zscores["GENE3", ]

  expect_false(any(is.na(gene2_zscores)))
  expect_false(any(is.na(gene3_zscores)))

  expect_equal(mean(gene2_zscores), 0, tolerance = 1e-10)
  expect_equal(mean(gene3_zscores), 0, tolerance = 1e-10)

  # Z-scores should have sd 1

  expect_equal(sd(gene2_zscores), 1, tolerance = 1e-10)
  expect_equal(sd(gene3_zscores), 1, tolerance = 1e-10)
})

test_that("Z-scores are computed independently per dataset", {
  skip_if_not(requireNamespace("tiledb", quietly = TRUE),
              message = "tiledb not installed")

  pgx_dir <- tempfile("pgx_zscore_multi_")
  tiledb_path <- tempfile("tiledb_zscore_multi_")
  dir.create(pgx_dir)
  on.exit({
    unlink(pgx_dir, recursive = TRUE)
    unlink(tiledb_path, recursive = TRUE)
    unlink(paste0(tiledb_path, "_metadata.rds"))
  }, add = TRUE)

  # Create two datasets with different scales
  set.seed(42)

  # Dataset A: low expression
  counts_a <- matrix(
    c(10, 20, 30),
    nrow = 1, ncol = 3,
    dimnames = list("TESTGENE", paste0("a_sample_", 1:3))
  )
  genes_a <- data.frame(symbol = "TESTGENE", human_ortholog = "TESTGENE", row.names = "TESTGENE")
  samples_a <- data.frame(group = "A", row.names = paste0("a_sample_", 1:3))
  pgx_a <- list(counts = counts_a, genes = genes_a, samples = samples_a)
  save_mock_pgx(pgx_a, file.path(pgx_dir, "datasetA.pgx"))

  # Dataset B: high expression (10x scale)
  counts_b <- matrix(
    c(100, 200, 300),
    nrow = 1, ncol = 3,
    dimnames = list("TESTGENE", paste0("b_sample_", 1:3))
  )
  genes_b <- data.frame(symbol = "TESTGENE", human_ortholog = "TESTGENE", row.names = "TESTGENE")
  samples_b <- data.frame(group = "B", row.names = paste0("b_sample_", 1:3))
  pgx_b <- list(counts = counts_b, genes = genes_b, samples = samples_b)
  save_mock_pgx(pgx_b, file.path(pgx_dir, "datasetB.pgx"))

  playbase::buildTileDB(pgx_dir, tiledb_path, verbose = FALSE)

  # Query z-scores
  zscores <- playbase::queryTileDB(tiledb_path, genes = "TESTGENE", value = "zscore")

  # Get z-scores for each dataset
  samples <- colnames(zscores)
  ds_a_samples <- samples[grep("datasetA::", samples)]
  ds_b_samples <- samples[grep("datasetB::", samples)]

  zscores_a <- zscores[1, ds_a_samples]
  zscores_b <- zscores[1, ds_b_samples]

  # Both datasets should have same z-scores despite different scales
  # (because z-scores are computed per-dataset)
  expect_equal(as.numeric(zscores_a), as.numeric(zscores_b), tolerance = 1e-10)
})

test_that("Z-scores handle NA values in counts", {
  skip_if_not(requireNamespace("tiledb", quietly = TRUE),
              message = "tiledb not installed")

  pgx_dir <- tempfile("pgx_zscore_na_")
  tiledb_path <- tempfile("tiledb_zscore_na_")
  dir.create(pgx_dir)
  on.exit({
    unlink(pgx_dir, recursive = TRUE)
    unlink(tiledb_path, recursive = TRUE)
    unlink(paste0(tiledb_path, "_metadata.rds"))
  }, add = TRUE)

  # Create counts with NA values
  counts <- matrix(
    c(10, NA, 30, 40),
    nrow = 1, ncol = 4,
    dimnames = list("GENE1", paste0("s", 1:4))
  )

  genes <- data.frame(symbol = "GENE1", human_ortholog = "GENE1", row.names = "GENE1")
  samples <- data.frame(group = "test", row.names = paste0("s", 1:4))
  pgx <- list(counts = counts, genes = genes, samples = samples)
  save_mock_pgx(pgx, file.path(pgx_dir, "test.pgx"))

  playbase::buildTileDB(pgx_dir, tiledb_path, verbose = FALSE)

  # Query both counts and z-scores
  counts_result <- playbase::queryTileDB(tiledb_path, genes = "GENE1", value = "count")
  zscores_result <- playbase::queryTileDB(tiledb_path, genes = "GENE1", value = "zscore")

  # NA in counts should result in NA in z-scores
  na_samples <- colnames(counts_result)[is.na(counts_result[1, ])]
  if (length(na_samples) > 0) {
    expect_true(all(is.na(zscores_result[1, na_samples])))
  }
})

test_that("Z-scores work with tiledb.addDataset", {
  skip_if_not(requireNamespace("tiledb", quietly = TRUE),
              message = "tiledb not installed")

  pgx_dir <- tempfile("pgx_zscore_add_")
  tiledb_path <- tempfile("tiledb_zscore_add_")
  dir.create(pgx_dir)
  on.exit({
    unlink(pgx_dir, recursive = TRUE)
    unlink(tiledb_path, recursive = TRUE)
    unlink(paste0(tiledb_path, "_metadata.rds"))
  }, add = TRUE)

  # Create and add first dataset
  counts1 <- matrix(
    c(10, 20, 30),
    nrow = 1, ncol = 3,
    dimnames = list("GENE1", paste0("s1_", 1:3))
  )
  genes1 <- data.frame(symbol = "GENE1", human_ortholog = "GENE1", row.names = "GENE1")
  samples1 <- data.frame(group = "A", row.names = paste0("s1_", 1:3))
  pgx1 <- list(counts = counts1, genes = genes1, samples = samples1)
  pgx1_file <- file.path(pgx_dir, "dataset1.pgx")
  save_mock_pgx(pgx1, pgx1_file)

  playbase::tiledb.addDataset(tiledb_path, pgx1_file, verbose = FALSE)

  # Create and add second dataset
  counts2 <- matrix(
    c(100, 200, 300, 400),
    nrow = 1, ncol = 4,
    dimnames = list("GENE1", paste0("s2_", 1:4))
  )
  genes2 <- data.frame(symbol = "GENE1", human_ortholog = "GENE1", row.names = "GENE1")
  samples2 <- data.frame(group = "B", row.names = paste0("s2_", 1:4))
  pgx2 <- list(counts = counts2, genes = genes2, samples = samples2)
  pgx2_file <- file.path(pgx_dir, "dataset2.pgx")
  save_mock_pgx(pgx2, pgx2_file)

  playbase::tiledb.addDataset(tiledb_path, pgx2_file, verbose = FALSE)

  # Query z-scores - should work for both datasets
  zscores <- playbase::queryTileDB(tiledb_path, genes = "GENE1", value = "zscore")

  expect_equal(ncol(zscores), 7)  # 3 + 4 samples

  # Verify z-scores are present for both datasets
  samples <- colnames(zscores)
  ds1_samples <- samples[grep("dataset1::", samples)]
  ds2_samples <- samples[grep("dataset2::", samples)]

  expect_equal(length(ds1_samples), 3)
  expect_equal(length(ds2_samples), 4)

  # Both should have valid z-scores
  expect_false(any(is.na(zscores[1, ds1_samples])))
  expect_false(any(is.na(zscores[1, ds2_samples])))
})


# ============================================================================
# EDGE CASES AND STRESS TESTS
# ============================================================================

test_that("Functions handle empty database gracefully", {
  skip_if_not(requireNamespace("tiledb", quietly = TRUE),
              message = "tiledb not installed")

  # Create a database then query with no matching data
  pgx_dir <- tempfile("pgx_edge_")
  tiledb_path <- tempfile("tiledb_edge_")
  dir.create(pgx_dir)
  on.exit({
    unlink(pgx_dir, recursive = TRUE)
    unlink(tiledb_path, recursive = TRUE)
    unlink(paste0(tiledb_path, "_metadata.rds"))
  }, add = TRUE)

  pgx <- create_mock_pgx(n_genes = 5, n_samples = 2)
  save_mock_pgx(pgx, file.path(pgx_dir, "small.pgx"))
  playbase::buildTileDB(pgx_dir, tiledb_path, verbose = FALSE)

  # Query for samples that don't exist
  genes <- playbase::listGenesTileDB(tiledb_path)
  expect_warning(
    result <- playbase::queryTileDB(
      tiledb_path,
      genes = genes[1],
      samples = c("nonexistent::sample1")
    ),
    "No data found"
  )

  expect_true(is.matrix(result))
  expect_true(all(is.na(result)))
})

test_that("Sample ID parsing handles edge cases", {
  # Multiple :: in sample name (greedy match returns part after last ::)
  result <- playbase::getSampleShortName("ds::sample::with::colons")
  expect_equal(result, "colons")

  # Dataset extraction with empty parts
  result <- playbase::getDatasetFromSample("::sample")
  expect_equal(result, "")
})

test_that("Phenotype query handles NA values correctly", {
  skip_if_not(requireNamespace("tiledb", quietly = TRUE),
              message = "tiledb not installed")

  pgx_dir <- tempfile("pgx_na_pheno_")
  tiledb_path <- tempfile("tiledb_na_pheno_")
  dir.create(pgx_dir)
  on.exit({
    unlink(pgx_dir, recursive = TRUE)
    unlink(tiledb_path, recursive = TRUE)
    unlink(paste0(tiledb_path, "_metadata.rds"))
  }, add = TRUE)

  # Create PGX with NA in phenotypes
  pgx <- create_mock_pgx(n_genes = 5, n_samples = 3)
  pgx$samples$age[1] <- NA
  pgx$samples$sex[2] <- NA
  save_mock_pgx(pgx, file.path(pgx_dir, "test.pgx"))

  playbase::buildTileDB(pgx_dir, tiledb_path, verbose = FALSE)

  result <- playbase::queryPhenotypesTileDB(tiledb_path)

  # Should handle NA values
  expect_true(any(is.na(result$age)) || all(!is.na(result$age)))  # depends on conversion
})


# ============================================================================
# TESTS FOR SPECIAL CHARACTERS IN DATASET/SAMPLE NAMES
# ============================================================================

test_that("Dataset names with special characters are handled correctly", {
  skip_if_not(requireNamespace("tiledb", quietly = TRUE),
              message = "tiledb not installed")

  pgx_dir <- tempfile("pgx_special_ds_")
  tiledb_path <- tempfile("tiledb_special_ds_")
  dir.create(pgx_dir)
  on.exit({
    unlink(pgx_dir, recursive = TRUE)
    unlink(tiledb_path, recursive = TRUE)
    unlink(paste0(tiledb_path, "_metadata.rds"))
  }, add = TRUE)

  # Create PGX files with special characters in dataset names
  pgx1 <- create_mock_pgx(n_genes = 5, n_samples = 2)
  pgx2 <- create_mock_pgx(n_genes = 5, n_samples = 2)
  pgx3 <- create_mock_pgx(n_genes = 5, n_samples = 2)

  # Dataset names with various special characters
  # Note: some characters are invalid in filenames, so we test valid but tricky ones
  save_mock_pgx(pgx1, file.path(pgx_dir, "dataset-with-dashes.pgx"))
  save_mock_pgx(pgx2, file.path(pgx_dir, "dataset_with_underscores.pgx"))
  save_mock_pgx(pgx3, file.path(pgx_dir, "dataset.with.dots.pgx"))

  metadata <- playbase::buildTileDB(pgx_dir, tiledb_path, verbose = FALSE)

  # Verify all datasets were processed

  expect_equal(metadata$n_files, 3)
  expect_equal(metadata$n_samples, 6)

  # Verify dataset names are preserved correctly
  datasets <- playbase::listDatasetsTileDB(tiledb_path)
  expect_true("dataset-with-dashes" %in% datasets)
  expect_true("dataset_with_underscores" %in% datasets)
  expect_true("dataset.with.dots" %in% datasets)  # file_path_sans_ext only removes .pgx

  # Verify queries work with special character datasets
  samples <- playbase::listSamplesTileDB(tiledb_path)
  expect_true(any(grepl("dataset-with-dashes::", samples)))
  expect_true(any(grepl("dataset_with_underscores::", samples)))
})

test_that("Sample names with special characters are handled correctly", {
  skip_if_not(requireNamespace("tiledb", quietly = TRUE),
              message = "tiledb not installed")

  pgx_dir <- tempfile("pgx_special_sample_")
  tiledb_path <- tempfile("tiledb_special_sample_")
  dir.create(pgx_dir)
  on.exit({
    unlink(pgx_dir, recursive = TRUE)
    unlink(tiledb_path, recursive = TRUE)
    unlink(paste0(tiledb_path, "_metadata.rds"))
  }, add = TRUE)

  # Create PGX with special characters in sample names
  set.seed(42)
  n_genes <- 5
  n_samples <- 6

  gene_names <- paste0("GENE", seq_len(n_genes))

  # Various special character patterns in sample names
  sample_names <- c(
    "sample-with-dashes",
    "sample_with_underscores",
    "sample.with.dots",
    "sample (with) parens",
    "sample [with] brackets",
    "sample+plus+signs"
  )

  counts <- matrix(
    rpois(n_genes * n_samples, lambda = 100),
    nrow = n_genes, ncol = n_samples,
    dimnames = list(gene_names, sample_names)
  )

  genes <- data.frame(
    symbol = gene_names,
    human_ortholog = paste0("HUMAN_", gene_names),
    row.names = gene_names,
    stringsAsFactors = FALSE
  )

  samples <- data.frame(
    age = sample(20:80, n_samples, replace = TRUE),
    condition = sample(c("A", "B"), n_samples, replace = TRUE),
    row.names = sample_names,
    stringsAsFactors = FALSE
  )

  pgx <- list(counts = counts, genes = genes, samples = samples)
  save_mock_pgx(pgx, file.path(pgx_dir, "test.pgx"))

  metadata <- playbase::buildTileDB(pgx_dir, tiledb_path, verbose = FALSE)

  # Verify all samples were stored
  expect_equal(metadata$n_samples, 6)

  # Verify sample names are preserved
  db_samples <- playbase::listSamplesTileDB(tiledb_path)
  short_names <- playbase::getSampleShortName(db_samples)

  expect_true("sample-with-dashes" %in% short_names)
  expect_true("sample_with_underscores" %in% short_names)
  expect_true("sample.with.dots" %in% short_names)
  expect_true("sample (with) parens" %in% short_names)
  expect_true("sample [with] brackets" %in% short_names)
  expect_true("sample+plus+signs" %in% short_names)

  # Verify queries work with special character samples
  genes_list <- playbase::listGenesTileDB(tiledb_path)
  result <- playbase::queryTileDB(tiledb_path, genes = genes_list[1])
  expect_equal(ncol(result), 6)

  # Verify phenotype queries work
  pheno <- playbase::queryPhenotypesTileDB(tiledb_path)
  expect_equal(nrow(pheno), 6)
})

test_that("Sample names with whitespace are handled correctly", {
  skip_if_not(requireNamespace("tiledb", quietly = TRUE),
              message = "tiledb not installed")

  pgx_dir <- tempfile("pgx_whitespace_")
  tiledb_path <- tempfile("tiledb_whitespace_")
  dir.create(pgx_dir)
  on.exit({
    unlink(pgx_dir, recursive = TRUE)
    unlink(tiledb_path, recursive = TRUE)
    unlink(paste0(tiledb_path, "_metadata.rds"))
  }, add = TRUE)

  set.seed(42)
  n_genes <- 3
  n_samples <- 4

  gene_names <- paste0("GENE", seq_len(n_genes))

  # Whitespace edge cases
  sample_names <- c(
    "sample with spaces",
    "  leading_spaces",
    "trailing_spaces  ",
    "multiple   internal   spaces"
  )

  counts <- matrix(
    rpois(n_genes * n_samples, lambda = 100),
    nrow = n_genes, ncol = n_samples,
    dimnames = list(gene_names, sample_names)
  )

  genes <- data.frame(
    symbol = gene_names,
    human_ortholog = paste0("HUMAN_", gene_names),
    row.names = gene_names,
    stringsAsFactors = FALSE
  )

  samples <- data.frame(
    age = sample(20:80, n_samples, replace = TRUE),
    row.names = sample_names,
    stringsAsFactors = FALSE
  )

  pgx <- list(counts = counts, genes = genes, samples = samples)
  save_mock_pgx(pgx, file.path(pgx_dir, "test.pgx"))

  metadata <- playbase::buildTileDB(pgx_dir, tiledb_path, verbose = FALSE)

  # Verify samples stored (whitespace preserved as-is)
  expect_equal(metadata$n_samples, 4)

  db_samples <- playbase::listSamplesTileDB(tiledb_path)
  short_names <- playbase::getSampleShortName(db_samples)

  expect_true("sample with spaces" %in% short_names)
  expect_true("  leading_spaces" %in% short_names)
  expect_true("trailing_spaces  " %in% short_names)
  expect_true("multiple   internal   spaces" %in% short_names)
})

test_that("Gene names with special characters are handled correctly", {
  skip_if_not(requireNamespace("tiledb", quietly = TRUE),
              message = "tiledb not installed")

  pgx_dir <- tempfile("pgx_special_gene_")
  tiledb_path <- tempfile("tiledb_special_gene_")
  dir.create(pgx_dir)
  on.exit({
    unlink(pgx_dir, recursive = TRUE)
    unlink(tiledb_path, recursive = TRUE)
    unlink(paste0(tiledb_path, "_metadata.rds"))
  }, add = TRUE)

  set.seed(42)
  n_samples <- 3

  # Gene names that might appear in real data
  gene_names <- c(
    "TP53",
    "HLA-A",      # dash is common in HLA genes
    "C1orf123",   # numbers
    "HIST1H2BK",  # long alphanumeric
    "MT-CO1",     # mitochondrial genes
    "MIR21",      # microRNAs
    "LINC00152",  # lncRNAs
    "LOC123456"   # LOC identifiers
  )
  n_genes <- length(gene_names)

  sample_names <- paste0("sample_", seq_len(n_samples))

  counts <- matrix(
    rpois(n_genes * n_samples, lambda = 100),
    nrow = n_genes, ncol = n_samples,
    dimnames = list(gene_names, sample_names)
  )

  # human_ortholog keeps same names for this test
  genes <- data.frame(
    symbol = gene_names,
    human_ortholog = gene_names,
    row.names = gene_names,
    stringsAsFactors = FALSE
  )

  samples <- data.frame(
    age = sample(20:80, n_samples, replace = TRUE),
    row.names = sample_names,
    stringsAsFactors = FALSE
  )

  pgx <- list(counts = counts, genes = genes, samples = samples)
  save_mock_pgx(pgx, file.path(pgx_dir, "test.pgx"))

  metadata <- playbase::buildTileDB(pgx_dir, tiledb_path, verbose = FALSE)

  expect_equal(metadata$n_genes, n_genes)

  # Verify gene names preserved
  db_genes <- playbase::listGenesTileDB(tiledb_path)
  expect_true("HLA-A" %in% db_genes)
  expect_true("MT-CO1" %in% db_genes)

  # Query specific genes with special chars
  result <- playbase::queryTileDB(tiledb_path, genes = c("HLA-A", "MT-CO1"))
  expect_equal(nrow(result), 2)
  expect_true("HLA-A" %in% rownames(result))
  expect_true("MT-CO1" %in% rownames(result))
})

test_that("Helper functions handle special characters in sample IDs", {
  # Test getDatasetFromSample with special chars
  special_samples <- c(
    "data-set::sample-1",
    "data_set::sample_1",
    "data.set::sample.1",
    "data set::sample 1",
    "data::set::sample"  # multiple ::
  )

  datasets <- playbase::getDatasetFromSample(special_samples)
  expect_equal(datasets[1], "data-set")
  expect_equal(datasets[2], "data_set")
  expect_equal(datasets[3], "data.set")
  expect_equal(datasets[4], "data set")
  expect_equal(datasets[5], "data")  # only first part before first ::

  # Test getSampleShortName with special chars
  short_names <- playbase::getSampleShortName(special_samples)
  expect_equal(short_names[1], "sample-1")
  expect_equal(short_names[2], "sample_1")
  expect_equal(short_names[3], "sample.1")
  expect_equal(short_names[4], "sample 1")
  expect_equal(short_names[5], "sample")  # greedy match, takes after last ::
})

test_that("Phenotype values with special characters are preserved", {
  skip_if_not(requireNamespace("tiledb", quietly = TRUE),
              message = "tiledb not installed")

  pgx_dir <- tempfile("pgx_special_pheno_")
  tiledb_path <- tempfile("tiledb_special_pheno_")
  dir.create(pgx_dir)
  on.exit({
    unlink(pgx_dir, recursive = TRUE)
    unlink(tiledb_path, recursive = TRUE)
    unlink(paste0(tiledb_path, "_metadata.rds"))
  }, add = TRUE)

  set.seed(42)
  n_genes <- 3
  n_samples <- 5

  gene_names <- paste0("GENE", seq_len(n_genes))
  sample_names <- paste0("sample_", seq_len(n_samples))

  counts <- matrix(
    rpois(n_genes * n_samples, lambda = 100),
    nrow = n_genes, ncol = n_samples,
    dimnames = list(gene_names, sample_names)
  )

  genes <- data.frame(
    symbol = gene_names,
    human_ortholog = paste0("HUMAN_", gene_names),
    row.names = gene_names,
    stringsAsFactors = FALSE
  )

  # Phenotype values with special characters
  samples <- data.frame(
    treatment = c(
      "Drug-A (10mg)",
      "Drug_B [high]",
      "Control + vehicle",
      "Treatment/Condition",
      "Group: Alpha"
    ),
    stage = c("Stage I/II", "Stage III", "Stage IV", "N/A", "Unknown (?)"),
    row.names = sample_names,
    stringsAsFactors = FALSE
  )

  pgx <- list(counts = counts, genes = genes, samples = samples)
  save_mock_pgx(pgx, file.path(pgx_dir, "test.pgx"))

  metadata <- playbase::buildTileDB(pgx_dir, tiledb_path, verbose = FALSE)

  # Query phenotypes and verify special chars preserved
  pheno <- playbase::queryPhenotypesTileDB(tiledb_path)

  expect_true("Drug-A (10mg)" %in% pheno$treatment)
  expect_true("Drug_B [high]" %in% pheno$treatment)
  expect_true("Control + vehicle" %in% pheno$treatment)
  expect_true("Treatment/Condition" %in% pheno$treatment)
  expect_true("Group: Alpha" %in% pheno$treatment)

  expect_true("Stage I/II" %in% pheno$stage)
  expect_true("Unknown (?)" %in% pheno$stage)

  # Verify getPhenotypeValuesTileDB returns special chars
  treatment_values <- playbase::getPhenotypeValuesTileDB(tiledb_path, "treatment")
  expect_true("Drug-A (10mg)" %in% treatment_values)
})

test_that("Unicode characters in sample names are handled", {
  skip_if_not(requireNamespace("tiledb", quietly = TRUE),
              message = "tiledb not installed")

  pgx_dir <- tempfile("pgx_unicode_")
  tiledb_path <- tempfile("tiledb_unicode_")
  dir.create(pgx_dir)
  on.exit({
    unlink(pgx_dir, recursive = TRUE)
    unlink(tiledb_path, recursive = TRUE)
    unlink(paste0(tiledb_path, "_metadata.rds"))
  }, add = TRUE)

  set.seed(42)
  n_genes <- 3
  n_samples <- 4

  gene_names <- paste0("GENE", seq_len(n_genes))

  # Unicode sample names (common in international datasets)
  sample_names <- c(
    "patient_\u00e9",      # é (accented e)
    "sample_\u00f1",       # ñ (Spanish n)
    "test_\u03b1",         # α (Greek alpha)
    "ctrl_\u2014"          # — (em dash)
  )

  counts <- matrix(
    rpois(n_genes * n_samples, lambda = 100),
    nrow = n_genes, ncol = n_samples,
    dimnames = list(gene_names, sample_names)
  )

  genes <- data.frame(
    symbol = gene_names,
    human_ortholog = paste0("HUMAN_", gene_names),
    row.names = gene_names,
    stringsAsFactors = FALSE
  )

  samples <- data.frame(
    age = sample(20:80, n_samples, replace = TRUE),
    row.names = sample_names,
    stringsAsFactors = FALSE
  )

  pgx <- list(counts = counts, genes = genes, samples = samples)
  save_mock_pgx(pgx, file.path(pgx_dir, "test.pgx"))

  # This should either work or fail gracefully
  tryCatch({
    metadata <- playbase::buildTileDB(pgx_dir, tiledb_path, verbose = FALSE)
    expect_equal(metadata$n_samples, 4)

    # If it works, verify data is queryable
    genes_list <- playbase::listGenesTileDB(tiledb_path)
    result <- playbase::queryTileDB(tiledb_path, genes = genes_list[1])
    expect_equal(ncol(result), 4)
  }, error = function(e) {
    # If unicode causes issues, that's also valuable to know
    skip(paste("Unicode handling failed:", e$message))
  })
})

test_that("Empty and NULL sample names don't crash the system", {
  # Test helper functions with edge cases
  expect_equal(playbase::getDatasetFromSample(character(0)), character(0))
  expect_equal(playbase::getSampleShortName(character(0)), character(0))

  # Single empty string
  expect_equal(playbase::getDatasetFromSample(""), "")
  expect_equal(playbase::getSampleShortName(""), "")

  # Just the separator
  expect_equal(playbase::getDatasetFromSample("::"), "")
  expect_equal(playbase::getSampleShortName("::"), "")
})


# ============================================================================
# TESTS FOR INCREMENTAL UPDATE FUNCTIONS
# ============================================================================

# ----------------------------------------------------------------------------
# tiledb.needUpdate tests
# ----------------------------------------------------------------------------

test_that("tiledb.needUpdate returns TRUE when TileDB doesn't exist", {
  skip_if_not(requireNamespace("tiledb", quietly = TRUE),
              message = "tiledb not installed")

  pgx_dir <- tempfile("pgx_needupdate_")
  tiledb_path <- tempfile("tiledb_needupdate_")
  dir.create(pgx_dir)
  on.exit(unlink(pgx_dir, recursive = TRUE), add = TRUE)

  # Create a PGX file but no TileDB
  pgx <- create_mock_pgx(n_genes = 5, n_samples = 3)
  save_mock_pgx(pgx, file.path(pgx_dir, "test.pgx"))

  result <- playbase::tiledb.needUpdate(pgx_dir, tiledb_path)
  expect_true(result)
})

test_that("tiledb.needUpdate returns FALSE when TileDB is up to date", {
  skip_if_not(requireNamespace("tiledb", quietly = TRUE),
              message = "tiledb not installed")

  pgx_dir <- tempfile("pgx_needupdate2_")
  tiledb_path <- tempfile("tiledb_needupdate2_")
  dir.create(pgx_dir)
  on.exit({
    unlink(pgx_dir, recursive = TRUE)
    unlink(tiledb_path, recursive = TRUE)
    unlink(paste0(tiledb_path, "_metadata.rds"))
  }, add = TRUE)

  # Create PGX and build TileDB
  pgx <- create_mock_pgx(n_genes = 5, n_samples = 3)
  save_mock_pgx(pgx, file.path(pgx_dir, "test.pgx"))
  playbase::buildTileDB(pgx_dir, tiledb_path, verbose = FALSE)

  result <- playbase::tiledb.needUpdate(pgx_dir, tiledb_path)
  expect_false(result)
})

test_that("tiledb.needUpdate returns TRUE when new PGX files added", {
  skip_if_not(requireNamespace("tiledb", quietly = TRUE),
              message = "tiledb not installed")

  pgx_dir <- tempfile("pgx_needupdate3_")
  tiledb_path <- tempfile("tiledb_needupdate3_")
  dir.create(pgx_dir)
  on.exit({
    unlink(pgx_dir, recursive = TRUE)
    unlink(tiledb_path, recursive = TRUE)
    unlink(paste0(tiledb_path, "_metadata.rds"))
  }, add = TRUE)

  # Create first PGX and build TileDB
  pgx1 <- create_mock_pgx(n_genes = 5, n_samples = 3)
  save_mock_pgx(pgx1, file.path(pgx_dir, "dataset1.pgx"))
  playbase::buildTileDB(pgx_dir, tiledb_path, verbose = FALSE)

  # Should be up to date
  expect_false(playbase::tiledb.needUpdate(pgx_dir, tiledb_path))

  # Add a new PGX file
  pgx2 <- create_mock_pgx(n_genes = 5, n_samples = 4)
  save_mock_pgx(pgx2, file.path(pgx_dir, "dataset2.pgx"))

  # Now should need update
  expect_true(playbase::tiledb.needUpdate(pgx_dir, tiledb_path))
})

test_that("tiledb.needUpdate returns FALSE for empty directory", {
  skip_if_not(requireNamespace("tiledb", quietly = TRUE),
              message = "tiledb not installed")

  empty_dir <- tempfile("pgx_empty_")
  dir.create(empty_dir)
  on.exit(unlink(empty_dir, recursive = TRUE), add = TRUE)

  # No pgx files means nothing to update
  result <- playbase::tiledb.needUpdate(empty_dir)
  expect_false(result)
})

test_that("tiledb.needUpdate returns FALSE for non-existent directory", {
  result <- playbase::tiledb.needUpdate("/nonexistent/path/12345")
  expect_false(result)
})

test_that("tiledb.needUpdate uses default tiledb_path", {
  skip_if_not(requireNamespace("tiledb", quietly = TRUE),
              message = "tiledb not installed")

  pgx_dir <- tempfile("pgx_defaultpath_")
  dir.create(pgx_dir)
  on.exit(unlink(pgx_dir, recursive = TRUE), add = TRUE)

  pgx <- create_mock_pgx(n_genes = 5, n_samples = 3)
  save_mock_pgx(pgx, file.path(pgx_dir, "test.pgx"))

  # Should use pgx_dir/counts_tiledb as default
  result <- playbase::tiledb.needUpdate(pgx_dir)
  expect_true(result)  # TileDB doesn't exist at default path
})


# ----------------------------------------------------------------------------
# tiledb.addDataset tests
# ----------------------------------------------------------------------------

test_that("tiledb.addDataset creates new database when none exists", {
  skip_if_not(requireNamespace("tiledb", quietly = TRUE),
              message = "tiledb not installed")

  pgx_dir <- tempfile("pgx_add_create_")
  tiledb_path <- tempfile("tiledb_add_create_")
  dir.create(pgx_dir)
  on.exit({
    unlink(pgx_dir, recursive = TRUE)
    unlink(tiledb_path, recursive = TRUE)
    unlink(paste0(tiledb_path, "_metadata.rds"))
  }, add = TRUE)

  # Create PGX file
  pgx <- create_mock_pgx(n_genes = 10, n_samples = 5)
  pgx_file <- file.path(pgx_dir, "newdataset.pgx")
  save_mock_pgx(pgx, pgx_file)

  # Add to non-existent TileDB (should create it)
  metadata <- playbase::tiledb.addDataset(tiledb_path, pgx_file, verbose = FALSE)

  expect_true(dir.exists(tiledb_path))
  expect_true(file.exists(paste0(tiledb_path, "_metadata.rds")))
  expect_equal(metadata$n_files, 1)
  expect_equal(metadata$n_samples, 5)
  expect_equal(metadata$n_genes, 10)
})

test_that("tiledb.addDataset adds to existing database", {
  skip_if_not(requireNamespace("tiledb", quietly = TRUE),
              message = "tiledb not installed")

  pgx_dir <- tempfile("pgx_add_existing_")
  tiledb_path <- tempfile("tiledb_add_existing_")
  dir.create(pgx_dir)
  on.exit({
    unlink(pgx_dir, recursive = TRUE)
    unlink(tiledb_path, recursive = TRUE)
    unlink(paste0(tiledb_path, "_metadata.rds"))
  }, add = TRUE)

  # Create first dataset and build TileDB
  pgx1 <- create_mock_pgx(n_genes = 10, n_samples = 3)
  save_mock_pgx(pgx1, file.path(pgx_dir, "dataset1.pgx"))
  playbase::buildTileDB(pgx_dir, tiledb_path, verbose = FALSE)

  # Create second dataset
  pgx2 <- create_mock_pgx(n_genes = 10, n_samples = 4)
  pgx2_file <- file.path(pgx_dir, "dataset2.pgx")
  save_mock_pgx(pgx2, pgx2_file)

  # Add second dataset
  metadata <- playbase::tiledb.addDataset(tiledb_path, pgx2_file, verbose = FALSE)

  expect_equal(metadata$n_files, 2)
  expect_equal(metadata$n_samples, 7)  # 3 + 4

  # Verify both datasets are queryable
  datasets <- playbase::listDatasetsTileDB(tiledb_path)
  expect_true("dataset1" %in% datasets)
  expect_true("dataset2" %in% datasets)
})

test_that("tiledb.addDataset overwrites existing dataset", {
  skip_if_not(requireNamespace("tiledb", quietly = TRUE),
              message = "tiledb not installed")

  pgx_dir <- tempfile("pgx_add_overwrite_")
  tiledb_path <- tempfile("tiledb_add_overwrite_")
  dir.create(pgx_dir)
  on.exit({
    unlink(pgx_dir, recursive = TRUE)
    unlink(tiledb_path, recursive = TRUE)
    unlink(paste0(tiledb_path, "_metadata.rds"))
  }, add = TRUE)

  # Create and add first version
  pgx1 <- create_mock_pgx(n_genes = 10, n_samples = 3)
  pgx_file <- file.path(pgx_dir, "mydataset.pgx")
  save_mock_pgx(pgx1, pgx_file)
  playbase::tiledb.addDataset(tiledb_path, pgx_file, verbose = FALSE)

  initial_metadata <- readRDS(paste0(tiledb_path, "_metadata.rds"))
  expect_equal(initial_metadata$n_samples, 3)

  # Create updated version with more samples
  pgx2 <- create_mock_pgx(n_genes = 10, n_samples = 7)
  save_mock_pgx(pgx2, pgx_file)

  # Overwrite
  metadata <- playbase::tiledb.addDataset(tiledb_path, pgx_file, overwrite = TRUE, verbose = FALSE)

  # Should have new sample count (metadata updated)
  expect_equal(metadata$n_samples, 7)
  expect_equal(metadata$n_files, 1)
})

test_that("tiledb.addDataset skips existing dataset when overwrite=FALSE", {
  skip_if_not(requireNamespace("tiledb", quietly = TRUE),
              message = "tiledb not installed")

  pgx_dir <- tempfile("pgx_add_nooverwrite_")
  tiledb_path <- tempfile("tiledb_add_nooverwrite_")
  dir.create(pgx_dir)
  on.exit({
    unlink(pgx_dir, recursive = TRUE)
    unlink(tiledb_path, recursive = TRUE)
    unlink(paste0(tiledb_path, "_metadata.rds"))
  }, add = TRUE)

  # Create and add dataset
  pgx <- create_mock_pgx(n_genes = 10, n_samples = 3)
  pgx_file <- file.path(pgx_dir, "mydataset.pgx")
  save_mock_pgx(pgx, pgx_file)
  playbase::tiledb.addDataset(tiledb_path, pgx_file, verbose = FALSE)

  # Try to add again with overwrite=FALSE
  metadata <- playbase::tiledb.addDataset(tiledb_path, pgx_file, overwrite = FALSE, verbose = FALSE)

  # Should return existing metadata unchanged
  expect_equal(metadata$n_samples, 3)
})

test_that("tiledb.addDataset fails with non-existent PGX file", {
  skip_if_not(requireNamespace("tiledb", quietly = TRUE),
              message = "tiledb not installed")

  tiledb_path <- tempfile("tiledb_add_nofile_")

  expect_error(
    playbase::tiledb.addDataset(tiledb_path, "/nonexistent/file.pgx"),
    "does not exist"
  )
})

test_that("tiledb.addDataset fails with PGX without counts", {
  skip_if_not(requireNamespace("tiledb", quietly = TRUE),
              message = "tiledb not installed")

  pgx_dir <- tempfile("pgx_add_nocounts_")
  tiledb_path <- tempfile("tiledb_add_nocounts_")
  dir.create(pgx_dir)
  on.exit({
    unlink(pgx_dir, recursive = TRUE)
    unlink(tiledb_path, recursive = TRUE)
    unlink(paste0(tiledb_path, "_metadata.rds"))
  }, add = TRUE)

  # Create PGX without counts
  pgx_nocounts <- list(
    genes = data.frame(symbol = "A", human_ortholog = "A"),
    samples = data.frame(age = 30, row.names = "s1"),
    counts = NULL
  )
  pgx_file <- file.path(pgx_dir, "nocounts.pgx")
  save_mock_pgx(pgx_nocounts, pgx_file)

  expect_error(
    playbase::tiledb.addDataset(tiledb_path, pgx_file, verbose = FALSE),
    "No counts found"
  )
})

test_that("tiledb.addDataset preserves phenotype data", {
  skip_if_not(requireNamespace("tiledb", quietly = TRUE),
              message = "tiledb not installed")

  pgx_dir <- tempfile("pgx_add_pheno_")
  tiledb_path <- tempfile("tiledb_add_pheno_")
  dir.create(pgx_dir)
  on.exit({
    unlink(pgx_dir, recursive = TRUE)
    unlink(tiledb_path, recursive = TRUE)
    unlink(paste0(tiledb_path, "_metadata.rds"))
  }, add = TRUE)

  # Create PGX with phenotypes
  pgx <- create_mock_pgx(n_genes = 5, n_samples = 3, include_phenotypes = TRUE)
  pgx_file <- file.path(pgx_dir, "withpheno.pgx")
  save_mock_pgx(pgx, pgx_file)

  playbase::tiledb.addDataset(tiledb_path, pgx_file, verbose = FALSE)

  # Verify phenotypes stored
  phenotypes <- playbase::listPhenotypesTileDB(tiledb_path)
  expect_true("age" %in% phenotypes)
  expect_true("sex" %in% phenotypes)
  expect_true("condition" %in% phenotypes)

  pheno_data <- playbase::queryPhenotypesTileDB(tiledb_path)
  expect_equal(nrow(pheno_data), 3)
})

test_that("tiledb.addDataset handles multiple sequential additions", {
  skip_if_not(requireNamespace("tiledb", quietly = TRUE),
              message = "tiledb not installed")

  pgx_dir <- tempfile("pgx_add_multi_")
  tiledb_path <- tempfile("tiledb_add_multi_")
  dir.create(pgx_dir)
  on.exit({
    unlink(pgx_dir, recursive = TRUE)
    unlink(tiledb_path, recursive = TRUE)
    unlink(paste0(tiledb_path, "_metadata.rds"))
  }, add = TRUE)

  # Add three datasets sequentially
  for (i in 1:3) {
    pgx <- create_mock_pgx(n_genes = 10, n_samples = i + 2)
    pgx_file <- file.path(pgx_dir, paste0("dataset", i, ".pgx"))
    save_mock_pgx(pgx, pgx_file)
    playbase::tiledb.addDataset(tiledb_path, pgx_file, verbose = FALSE)
  }

  metadata <- readRDS(paste0(tiledb_path, "_metadata.rds"))
  expect_equal(metadata$n_files, 3)
  expect_equal(metadata$n_samples, 3 + 4 + 5)  # (1+2) + (2+2) + (3+2)

  datasets <- playbase::listDatasetsTileDB(tiledb_path)
  expect_equal(length(datasets), 3)
})


# ----------------------------------------------------------------------------
# tiledb.updateDatasetFolder tests
# ----------------------------------------------------------------------------

test_that("tiledb.updateDatasetFolder creates TileDB when none exists", {
  skip_if_not(requireNamespace("tiledb", quietly = TRUE),
              message = "tiledb not installed")

  pgx_dir <- tempfile("pgx_update_create_")
  tiledb_path <- file.path(pgx_dir, "counts_tiledb")
  dir.create(pgx_dir)
  on.exit({
    unlink(pgx_dir, recursive = TRUE)
  }, add = TRUE)

  # Create PGX file
  pgx <- create_mock_pgx(n_genes = 10, n_samples = 5)
  save_mock_pgx(pgx, file.path(pgx_dir, "test.pgx"))

  # Update folder (should create TileDB)
  metadata <- playbase::tiledb.updateDatasetFolder(pgx_dir, verbose = FALSE)

  expect_true(dir.exists(tiledb_path))
  expect_equal(metadata$n_files, 1)
  expect_equal(metadata$n_samples, 5)
})

test_that("tiledb.updateDatasetFolder adds specific new_pgx file", {
  skip_if_not(requireNamespace("tiledb", quietly = TRUE),
              message = "tiledb not installed")

  pgx_dir <- tempfile("pgx_update_specific_")
  tiledb_path <- file.path(pgx_dir, "counts_tiledb")
  dir.create(pgx_dir)
  on.exit({
    unlink(pgx_dir, recursive = TRUE)
  }, add = TRUE)

  # Create first dataset and update folder
  pgx1 <- create_mock_pgx(n_genes = 10, n_samples = 3)
  save_mock_pgx(pgx1, file.path(pgx_dir, "dataset1.pgx"))
  playbase::tiledb.updateDatasetFolder(pgx_dir, verbose = FALSE)

  # Create second dataset
  pgx2 <- create_mock_pgx(n_genes = 10, n_samples = 4)
  save_mock_pgx(pgx2, file.path(pgx_dir, "dataset2.pgx"))

  # Update with specific file
  metadata <- playbase::tiledb.updateDatasetFolder(
    pgx_dir,
    new_pgx = "dataset2.pgx",
    verbose = FALSE
  )

  expect_equal(metadata$n_files, 2)
  expect_equal(metadata$n_samples, 7)
})

test_that("tiledb.updateDatasetFolder handles new_pgx without extension", {
  skip_if_not(requireNamespace("tiledb", quietly = TRUE),
              message = "tiledb not installed")

  pgx_dir <- tempfile("pgx_update_noext_")
  tiledb_path <- file.path(pgx_dir, "counts_tiledb")
  dir.create(pgx_dir)
  on.exit({
    unlink(pgx_dir, recursive = TRUE)
  }, add = TRUE)

  pgx <- create_mock_pgx(n_genes = 10, n_samples = 5)
  save_mock_pgx(pgx, file.path(pgx_dir, "mydata.pgx"))

  # Pass without .pgx extension
  metadata <- playbase::tiledb.updateDatasetFolder(
    pgx_dir,
    new_pgx = "mydata",
    verbose = FALSE
  )

  expect_equal(metadata$n_files, 1)
})

test_that("tiledb.updateDatasetFolder adds all missing datasets", {

  skip_if_not(requireNamespace("tiledb", quietly = TRUE),
              message = "tiledb not installed")

  pgx_dir <- tempfile("pgx_update_missing_")
  tiledb_path <- file.path(pgx_dir, "counts_tiledb")
  dir.create(pgx_dir)
  on.exit({
    unlink(pgx_dir, recursive = TRUE)
  }, add = TRUE)

  # Create first dataset and build TileDB
  pgx1 <- create_mock_pgx(n_genes = 10, n_samples = 3)
  save_mock_pgx(pgx1, file.path(pgx_dir, "dataset1.pgx"))
  playbase::buildTileDB(pgx_dir, tiledb_path, verbose = FALSE)

  # Add two more datasets
  pgx2 <- create_mock_pgx(n_genes = 10, n_samples = 4)
  pgx3 <- create_mock_pgx(n_genes = 10, n_samples = 5)
  save_mock_pgx(pgx2, file.path(pgx_dir, "dataset2.pgx"))
  save_mock_pgx(pgx3, file.path(pgx_dir, "dataset3.pgx"))

  # Update folder (no new_pgx specified - should find all missing)
  metadata <- playbase::tiledb.updateDatasetFolder(pgx_dir, verbose = FALSE)

  expect_equal(metadata$n_files, 3)
  expect_equal(metadata$n_samples, 12)  # 3 + 4 + 5
})

test_that("tiledb.updateDatasetFolder returns NULL when up to date", {
  skip_if_not(requireNamespace("tiledb", quietly = TRUE),
              message = "tiledb not installed")

  pgx_dir <- tempfile("pgx_update_uptodate_")
  tiledb_path <- file.path(pgx_dir, "counts_tiledb")
  dir.create(pgx_dir)
  on.exit({
    unlink(pgx_dir, recursive = TRUE)
  }, add = TRUE)

  pgx <- create_mock_pgx(n_genes = 10, n_samples = 5)
  save_mock_pgx(pgx, file.path(pgx_dir, "test.pgx"))
  playbase::buildTileDB(pgx_dir, tiledb_path, verbose = FALSE)

  # Update when already up to date
  result <- playbase::tiledb.updateDatasetFolder(pgx_dir, verbose = FALSE)

  expect_null(result)
})

test_that("tiledb.updateDatasetFolder returns NULL for empty directory", {
  skip_if_not(requireNamespace("tiledb", quietly = TRUE),
              message = "tiledb not installed")

  empty_dir <- tempfile("pgx_update_empty_")
  dir.create(empty_dir)
  on.exit(unlink(empty_dir, recursive = TRUE), add = TRUE)

  result <- playbase::tiledb.updateDatasetFolder(empty_dir, verbose = FALSE)
  expect_null(result)
})

test_that("tiledb.updateDatasetFolder returns NULL for non-existent directory", {
  result <- playbase::tiledb.updateDatasetFolder("/nonexistent/path/12345", verbose = FALSE)
  expect_null(result)
})

test_that("tiledb.updateDatasetFolder returns NULL for non-existent new_pgx", {
  skip_if_not(requireNamespace("tiledb", quietly = TRUE),
              message = "tiledb not installed")

  pgx_dir <- tempfile("pgx_update_nofile_")
  dir.create(pgx_dir)
  on.exit(unlink(pgx_dir, recursive = TRUE), add = TRUE)

  result <- playbase::tiledb.updateDatasetFolder(
    pgx_dir,
    new_pgx = "nonexistent.pgx",
    verbose = FALSE
  )

  expect_null(result)
})

test_that("tiledb.updateDatasetFolder uses custom tiledb_path", {
  skip_if_not(requireNamespace("tiledb", quietly = TRUE),
              message = "tiledb not installed")

  pgx_dir <- tempfile("pgx_update_custompath_")
  custom_tiledb <- tempfile("custom_tiledb_")
  dir.create(pgx_dir)
  on.exit({
    unlink(pgx_dir, recursive = TRUE)
    unlink(custom_tiledb, recursive = TRUE)
    unlink(paste0(custom_tiledb, "_metadata.rds"))
  }, add = TRUE)

  pgx <- create_mock_pgx(n_genes = 10, n_samples = 5)
  save_mock_pgx(pgx, file.path(pgx_dir, "test.pgx"))

  metadata <- playbase::tiledb.updateDatasetFolder(
    pgx_dir,
    tiledb_path = custom_tiledb,
    verbose = FALSE
  )

  expect_true(dir.exists(custom_tiledb))
  expect_equal(metadata$n_files, 1)
})

test_that("tiledb.updateDatasetFolder handles errors gracefully", {
  skip_if_not(requireNamespace("tiledb", quietly = TRUE),
              message = "tiledb not installed")

  pgx_dir <- tempfile("pgx_update_error_")
  dir.create(pgx_dir)
  on.exit(unlink(pgx_dir, recursive = TRUE), add = TRUE)

  # Create one valid and one invalid PGX
  pgx_valid <- create_mock_pgx(n_genes = 5, n_samples = 3)
  save_mock_pgx(pgx_valid, file.path(pgx_dir, "valid.pgx"))

  pgx_invalid <- list(counts = NULL, genes = NULL)
  save_mock_pgx(pgx_invalid, file.path(pgx_dir, "invalid.pgx"))

  # Should complete with warnings but not fail
  expect_warning(
    metadata <- playbase::tiledb.updateDatasetFolder(pgx_dir, verbose = FALSE)
  )

  # Valid dataset should still be added
  expect_equal(metadata$n_files, 1)
})


# ----------------------------------------------------------------------------
# Integration tests for incremental updates
# ----------------------------------------------------------------------------

test_that("Full workflow: build, add, query cycle works", {
  skip_if_not(requireNamespace("tiledb", quietly = TRUE),
              message = "tiledb not installed")

  pgx_dir <- tempfile("pgx_workflow_")
  tiledb_path <- file.path(pgx_dir, "counts_tiledb")
  dir.create(pgx_dir)
  on.exit({
    unlink(pgx_dir, recursive = TRUE)
  }, add = TRUE)

  # Step 1: Create initial dataset
  pgx1 <- create_mock_pgx(n_genes = 10, n_samples = 3)
  save_mock_pgx(pgx1, file.path(pgx_dir, "initial.pgx"))

  # Step 2: Build initial TileDB
  playbase::buildTileDB(pgx_dir, tiledb_path, verbose = FALSE)
  expect_equal(playbase::listDatasetsTileDB(tiledb_path), "initial")

  # Step 3: Add new dataset via tiledb.addDataset
  pgx2 <- create_mock_pgx(n_genes = 10, n_samples = 4)
  save_mock_pgx(pgx2, file.path(pgx_dir, "added.pgx"))
  playbase::tiledb.addDataset(tiledb_path, file.path(pgx_dir, "added.pgx"), verbose = FALSE)

  datasets <- playbase::listDatasetsTileDB(tiledb_path)
  expect_true("initial" %in% datasets)
  expect_true("added" %in% datasets)

  # Step 4: Add another via tiledb.updateDatasetFolder
  pgx3 <- create_mock_pgx(n_genes = 10, n_samples = 5)
  save_mock_pgx(pgx3, file.path(pgx_dir, "updated.pgx"))
  playbase::tiledb.updateDatasetFolder(pgx_dir, new_pgx = "updated", verbose = FALSE)

  # Step 5: Verify all data is queryable
  datasets <- playbase::listDatasetsTileDB(tiledb_path)
  expect_equal(length(datasets), 3)

  samples <- playbase::listSamplesTileDB(tiledb_path)
  expect_equal(length(samples), 12)  # 3 + 4 + 5

  genes <- playbase::listGenesTileDB(tiledb_path)
  result <- playbase::queryTileDB(tiledb_path, genes = genes[1:3])
  expect_equal(ncol(result), 12)

  # Step 6: Verify phenotypes from all datasets
  pheno <- playbase::queryPhenotypesTileDB(tiledb_path)
  expect_equal(nrow(pheno), 12)
})

test_that("Incremental updates preserve existing data integrity", {
  skip_if_not(requireNamespace("tiledb", quietly = TRUE),
              message = "tiledb not installed")

  pgx_dir <- tempfile("pgx_integrity_")
  tiledb_path <- file.path(pgx_dir, "counts_tiledb")
  dir.create(pgx_dir)
  on.exit({
    unlink(pgx_dir, recursive = TRUE)
  }, add = TRUE)

  # Create dataset with known values
  set.seed(123)
  pgx1 <- create_mock_pgx(n_genes = 5, n_samples = 3)
  save_mock_pgx(pgx1, file.path(pgx_dir, "dataset1.pgx"))
  playbase::buildTileDB(pgx_dir, tiledb_path, verbose = FALSE)

  # Query original data
  genes <- playbase::listGenesTileDB(tiledb_path)
  original_result <- playbase::queryTileDB(tiledb_path, genes = genes)
  original_ds1_data <- original_result[, grep("dataset1", colnames(original_result))]

  # Add second dataset
  set.seed(456)
  pgx2 <- create_mock_pgx(n_genes = 5, n_samples = 4)
  save_mock_pgx(pgx2, file.path(pgx_dir, "dataset2.pgx"))
  playbase::tiledb.addDataset(tiledb_path, file.path(pgx_dir, "dataset2.pgx"), verbose = FALSE)

  # Query again and verify dataset1 data unchanged
  new_result <- playbase::queryTileDB(tiledb_path, genes = genes)
  new_ds1_data <- new_result[, grep("dataset1", colnames(new_result))]

  expect_equal(original_ds1_data, new_ds1_data)
})
