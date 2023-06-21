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

   if (!is.null(samples) && !is.null(counts)) {
    # Check that rownames(samples) match colnames(counts)
    SAMPLE_NAMES_NOT_MATCHING_COUNTS <- intersect(
      rownames(samples),
      colnames(counts)
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
    
    nsamples <- max(ncol(counts), nrow(samples))
    
    if (
      length(SAMPLE_NAMES_PARTIAL_MATCHING_COUNTS) > 0 && 
      length(SAMPLE_NAMES_PARTIAL_MATCHING_COUNTS) < nsamples &&
      PASS
      ) {
      TOTAL_SAMPLES_NAMES <- unique(c(rownames(samples),colnames(counts)))
      check_return$e19 <- TOTAL_SAMPLES_NAMES[!TOTAL_SAMPLES_NAMES %in% SAMPLE_NAMES_PARTIAL_MATCHING_COUNTS]
      samples <- samples[SAMPLE_NAMES_PARTIAL_MATCHING_COUNTS, , drop = FALSE]
      counts <- counts[, SAMPLE_NAMES_PARTIAL_MATCHING_COUNTS, drop = FALSE]
      pass = TRUE
    }
   }

   if (!is.null(samples) && !is.null(contrasts)) {
    contrasts <- contrasts_conversion(samples, contrasts)

    # Check that rownames(samples) match long contrast rownames.

    if(dim(samples)[1] == dim(contrasts)[1]){ # check that contrasts are in long format
      SAMPLE_NAMES_NOT_MATCHING_CONTRASTS <- intersect(
      rownames(samples),
      rownames(contrasts)
      )
    }

    if (length(SAMPLE_NAMES_NOT_MATCHING_CONTRASTS) == 0 && PASS) {
      check_return$e17 <- SAMPLE_NAMES_NOT_MATCHING_CONTRASTS
      pass = FALSE
    }
    # Check that counts have the same order as samples.

    MATCH_SAMPLES_COUNTS_ORDER <- all(diff(match(rownames(samples), colnames(counts)))>0)

    if(!MATCH_SAMPLES_COUNTS_ORDER && PASS){
      check_return$e18 <- "samples and counts do not have the same order"
      pass = TRUE
      counts <- counts[,match(rownames(samples), colnames(counts))]
      }
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

#' Convert contrasts for OPG
#'
#' @param SAMPLE data.frame. The data frame corresponding to the input file as in playbase::SAMPLES
#' @param CONTRASTS data.frame. The data frame corresponding to the input file as in playbase::CONTRASTS
#'
#' @return converted contrast df
#' @export
#'
#' @examples
contrasts_conversion <- function(SAMPLES, CONTRASTS){
  samples1 <- SAMPLES
  contrasts1 <- CONTRASTS
  group.col <- grep("group", tolower(colnames(samples1)))
  old1 <- (length(group.col) > 0 &&
    nrow(contrasts1) < nrow(samples1) &&
    all(rownames(contrasts1) %in% samples1[, group.col[1]])
  )
  old2 <- all(rownames(contrasts1) == rownames(samples1)) &&
    all(unique(as.vector(contrasts1)) %in% c(-1, 0, 1, NA))

  old.style <- (old1 || old2)
  if (old.style && old1) {
    message("[UploadModule] WARNING: converting old1 style contrast to new format")
    new.contrasts <- samples1[, 0]
    if (NCOL(contrasts1) > 0) {
      new.contrasts <- playbase::contrastAsLabels(contrasts1)
      grp <- as.character(samples1[, group.col])
      new.contrasts <- new.contrasts[grp, , drop = FALSE]
      rownames(new.contrasts) <- rownames(samples1)
    }
    contrasts1 <- new.contrasts
  }
  if (old.style && old2) {
    message("[UploadModule] WARNING: converting old2 style contrast to new format")
    new.contrasts <- samples1[, 0]
    if (NCOL(contrasts1) > 0) {
      new.contrasts <- playbase::contrastAsLabels(contrasts1)
      rownames(new.contrasts) <- rownames(samples1)
    }
    contrasts1 <- new.contrasts
  }

  dbg("[UploadModule] 1 : dim.contrasts1 = ", dim(contrasts1))
  dbg("[UploadModule] 1 : dim.samples1   = ", dim(samples1))

  ok.contrast <- length(intersect(rownames(samples1), rownames(contrasts1))) > 0
  if (ok.contrast && NCOL(contrasts1) > 0) {
    ## always clean up
    contrasts1 <- apply(contrasts1, 2, as.character)
    rownames(contrasts1) <- rownames(samples1)
    for (i in 1:ncol(contrasts1)) {
      isz <- (contrasts1[, i] %in% c(NA, "NA", "NA ", "", " ", "  ", "   ", " NA"))
      if (length(isz)) contrasts1[isz, i] <- NA
    }
  }

  return(contrasts1)
  
}

