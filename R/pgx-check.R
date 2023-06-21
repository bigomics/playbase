#' Check input files for pgx.computePGX
#'
#' @param df data.frame. The data frame corresponding to the input file as
#'   described in Omics Playground documentation
#' @param type character. One of "SAMPLES", "COUNTS", "EXPRESSION", "CONTRASTS" - The type
#'  of data file to check.
#'
#' @return a list with two elements: `checks` which contains the status of the checks, and
#'  `PASS` which contains the overall status of the check.
#' @export
#'
#' @examples
pgx.checkPGX <- function(
    df,
    type = c("SAMPLES", "COUNTS", "EXPRESSION", "CONTRASTS")
) {
  datatype <- match.arg(type)
  df_clean <- df

  PASS = TRUE

  check_return <- list()

  if (datatype == "COUNTS" || datatype == "EXPRESSION") {
    feature_names <- rownames(df_clean)

    # check for duplicated rownames (but pass)
    ANY_DUPLICATED <- unique(feature_names[which(duplicated(feature_names))])

    if (length(x = ANY_DUPLICATED) > 0 && PASS) {
      check_return$e6 <- ANY_DUPLICATED
    }

    # check for zero count rows, remove them

    ANY_ROW_ZERO <- which(rowSums(df_clean)==0)

    if (length(ANY_ROW_ZERO) > 0 && PASS) {

      # get the row names with all zeros
      check_return$e9 <- names(ANY_ROW_ZERO)

      # remove the rownames with all zeros by using check_return$e9
      df_clean <- df_clean[!(rownames(df_clean) %in% check_return$e9), ]
    }

    # check for zero count columns, remove them
    ANY_COLUMN_ZERO <- which(colSums(df_clean)==0)

    if (length(ANY_COLUMN_ZERO) > 0 && PASS) {
      check_return$e10 <- names(ANY_COLUMN_ZERO)
      # remove the column names with all zeros by using check_return$e9
      df_clean <- df_clean[, !(colnames(df_clean) %in% check_return$e10)]
    }

  }

  if (datatype == "SAMPLES") {
    feature_names <- rownames(df_clean)

    # check for duplicated rownames

    ANY_DUPLICATED <- unique(feature_names[which(duplicated(feature_names))])

    if (length(x = ANY_DUPLICATED) > 0 && PASS) {
      check_return$e1 <- ANY_DUPLICATED
      PASS = FALSE
    }
  }

  if (datatype == "CONTRASTS") {
    feature_names <- rownames(df_clean)

    # check for duplicated rownames (but pass)
    ANY_DUPLICATED <- unique(feature_names[which(duplicated(feature_names))])

    if (length(x = ANY_DUPLICATED) > 0 && PASS) {
      check_return$e11 <- ANY_DUPLICATED
      PASS = FALSE
    }

  }

  # general checks for all data datatypes

  # check for empty df
  IS_DF_EMPTY <- dim(df_clean)[1] == 0 || dim(df_clean)[2] == 0

  if (!IS_DF_EMPTY) {
    df_clean <- as.matrix(df_clean)
  }

  if (IS_DF_EMPTY && PASS) {
    check_return$e15 <- "empty dataframe"
    pass = FALSE
  }

  return(
    list(
      df = df_clean,
      checks = check_return,
      PASS = PASS)
  )
}

#' Check all input files for pgx.computePGX
#'
#' @param SAMPLE data.frame. The data frame corresponding to the input file as in playbase::SAMPLES
#' @param COUNTS data.frame. The data frame corresponding to the input file as in playbase::COUNTS 
#' @param CONTRASTS data.frame. The data frame corresponding to the input file as in playbase::CONTRASTS
#'
#' @return a list with FIVE elements: SAMPLES, COUNTS and CONTRASTS that are the cleaned version of the
#'  input data frames, `checks` which contains the status of the checks, and
#'  `PASS` which contains the overall status of the check.
#' @export
#'
#' @examples
pgx.checkPGX_all <- function(
    SAMPLES = NULL,
    COUNTS = NULL,
    CONTRASTS = NULL
) {

  samples = SAMPLES
  counts = COUNTS
  contrasts = CONTRASTS
  PASS = TRUE

  check_return <- list()
  
  # Check that rownames(samples) match colnames(counts)
  SAMPLE_NAMES_NOT_MATCHING_COUNTS <- intersect(
    rownames(samples),
    colnames(uploaded$counts.csv)
  )

  
  if (length(SAMPLE_NAMES_NOT_MATCHING_COUNTS) == 0 && PASS) {
    check_return$e16 <- SAMPLE_NAMES_NOT_MATCHING_COUNTS
    pass = FALSE
  }

  # Check that rownames(samples) match colnames(counts)
  SAMPLE_NAMES_PARTIAL_MATCHING_COUNTS <- intersect(
    rownames(samples),
    colnames(counts)
  )

  nsamples <- max(ncol(counts, nrow(samples))
  ok.samples <- intersect(
            rownames(uploaded$samples.csv),
            colnames(uploaded$counts.csv)
          )


  if (n.ok > 0 && n.ok < nsamples && PASS) {
    check_return$e19 <- SAMPLE_NAMES_PARTIAL_MATCHING_COUNTS
    pass = TRUE
  }

  # Check that rownames(samples) match long contrast rownames.

  if(dim(samples)[1] == dim(contrasts)[1]){ # check that contrasts are in long format
    SAMPLE_NAMES_NOT_MATCHING_CONTRASTS <- intersect(
    rownames(samples),
    rownames(contrasts)
  )

    if (length(SAMPLE_NAMES_NOT_MATCHING_CONTRASTS) == 0 && PASS) {
      check_return$e17 <- SAMPLE_NAMES_NOT_MATCHING_CONTRASTS
      pass = FALSE
    }
  }

  # Check that counts have the same order as samples.

  MATCH_SAMPLES_COUNTS_ORDER <- all(diff(match(rownames(samples), colnames(counts)))>0)

  if(!MATCH_SAMPLES_COUNTS_ORDER && PASS){
    check_return$e18 <- "samples and counts do not have the same order"
    pass = TRUE
    counts <- counts[,match(rownames(samples), colnames(counts))]
  }

  return(
    list(
      SAMPLES = samples,
      COUNTS = counts,
      CONTRASTS = contrasts,
      checks = check_return,
      PASS = PASS)
  )
}