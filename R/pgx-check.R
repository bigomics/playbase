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
pgx.checkINPUT <- function(
    df,
    type = c("SAMPLES", "COUNTS", "EXPRESSION", "CONTRASTS")) {
  datatype <- match.arg(type)
  df_clean <- df
  PASS <- TRUE
  check_return <- list()

  if (datatype == "COUNTS" || datatype == "EXPRESSION") {
    # remove special characters before converting to numeric (keep commas, dots) as they define the decimal separator
    df_clean <- apply(df_clean, 2, function(x) gsub("[^0-9,.]", "", x))
    # convert matrix from character to numeric, sometimes we receive character matrix from read.csv function
    df_clean <- apply(df_clean, 2, as.numeric, simplify = TRUE)
    rownames(df_clean) <- rownames(df)

    sample_names <- colnames(df_clean)

    # check if any value in df had a non-NA value that was converted to NA in df_clean
    ANY_NON_NUMERIC <- which(!is.na(df) & is.na(df_clean), arr.ind = TRUE)

    if (length(ANY_NON_NUMERIC) > 0 && PASS) {
      check_return$e27 <- paste("gene:", rownames(ANY_NON_NUMERIC), " and ", "sample:", colnames(df_clean)[ANY_NON_NUMERIC[, 2]])
    }

    # replace infinite values in counts by the maximum value of the feature + 10%
    # check if there are any infinite values
    ANY_INFINITE <- which(is.infinite(df_clean), arr.ind = TRUE)

    if (length(ANY_INFINITE) > 0 && PASS) {
      df_clean <- apply(df_clean, 1, function(x) {
        max_val <- max(x, na.rm = TRUE)
        ifelse(is.infinite(x), max_val * 1.1, x)
      })
      check_return$e28 <- paste("gene:", rownames(ANY_INFINITE), " and ", "sample:", colnames(df_clean)[ANY_INFINITE[, 2]])
    }

    # check for duplicated colnanes (gives error)
    ANY_DUPLICATED <- unique(sample_names[which(duplicated(sample_names))])

    if (length(ANY_DUPLICATED) > 0 && PASS) {
      PASS <- FALSE
      check_return$e6 <- ANY_DUPLICATED
    }

    feature_names <- rownames(df_clean)

    # check for duplicated rownames (but pass)
    ANY_DUPLICATED <- unique(feature_names[which(duplicated(feature_names))])

    if (length(ANY_DUPLICATED) > 0 && PASS) {
      check_return$e7 <- ANY_DUPLICATED
    }

    # check for zero count rows, remove them
    ANY_ROW_ZERO <- which(rowSums(df_clean, na.rm = TRUE) == 0, )

    if (length(ANY_ROW_ZERO) > 0 && PASS) {
      # get the row names with all zeros
      zero.rows <- names(ANY_ROW_ZERO)

      # remove the rownames with all zeros
      df_clean <- df_clean[!(rownames(df_clean) %in% zero.rows), , drop = FALSE]

      nzerorows <- length(ANY_ROW_ZERO)
      if (nzerorows < 10) {
        err.mesg <- zero.rows
      } else {
        err.mesg <- c(head(zero.rows, 10), "+more", paste("(total", nzerorows, "rows)"))
      }
      check_return$e9 <- err.mesg
    }

    # check for zero count columns, remove them
    ANY_COLUMN_ZERO <- which(colSums(df_clean) == 0)

    if (length(ANY_COLUMN_ZERO) > 0 && PASS) {
      check_return$e10 <- names(ANY_COLUMN_ZERO)
      # remove the column names with all zeros by using check_return$e9
      df_clean <- df_clean[, !(colnames(df_clean) %in% check_return$e10)]
    }

    # check if counts are log transformed
    check.log <- is_logged(df_clean)
    if (check.log) {
      check_return$e29 <- "Possible log transformed counts detected."
    }
  }

  if (datatype == "SAMPLES") {
    feature_names <- rownames(df_clean)

    # check for duplicated rownames
    ANY_DUPLICATED <- unique(feature_names[which(duplicated(feature_names))])

    if (length(x = ANY_DUPLICATED) > 0 && PASS) {
      check_return$e1 <- ANY_DUPLICATED
      PASS <- FALSE
    }
  }

  if (datatype == "CONTRASTS") {
    feature_names <- rownames(df_clean)

    # check that contrasts has at least one column

    COMPARISONS_WITHOUT_COLUMNS <- dim(df_clean)[2] == 0

    if (COMPARISONS_WITHOUT_COLUMNS && PASS) {
      check_return$e26 <- "No columns provided in comparisons."
      PASS <- FALSE
    }

    # check for duplicated rownames (but pass)
    ANY_DUPLICATED <- unique(feature_names[which(duplicated(feature_names))])

    if (length(x = ANY_DUPLICATED) > 0 && PASS) {
      check_return$e11 <- ANY_DUPLICATED
      PASS <- FALSE
    }

    ## check that numerator_vs_denominator is in the contrasts
    if (all(grepl(" vs ", colnames(df_clean)))) {
      colnames(df_clean) <- gsub(" vs ", "_vs_", colnames(df_clean))
    }
    has_no_vs <- which(!grepl("_vs_", colnames(df_clean)))
    if (length(has_no_vs) > 0 && PASS) {
      check_return$e24 <- colnames(df_clean)[has_no_vs]
      PASS <- FALSE
    }

    if (PASS) {
      # Split the column names at "_vs_"
      split_names <- strsplit(colnames(df_clean), "_vs_")

      # Get the numerators and denominators
      numerators <- sapply(split_names, "[[", 1)
      split_numerators <- strsplit(numerators, ":")

      # if colon is present in numerators, keep the elements after colon
      numerators <- sapply(split_numerators, function(x) {
        if (length(x) > 1) {
          x[2]
        } else {
          x[1]
        }
      })

      denominators <- sapply(split_names, "[[", 2)
      ## Check if all elements in the matrix are character
      all_numeric <- any(apply(df_clean, c(1, 2), is.numeric))

      if (!all_numeric && PASS) {
        ## only run if we have characters in matrix

        COLUMN_IN_GROUPS <- sapply(1:length(denominators), function(i) {
          vv <- setdiff(df_clean[, i], c(NA, "", " ", "NA"))
          all(grepl(paste0("^", numerators[i], "|^", denominators[i]), vv))
        })
        CONTRASTS_IN_GROUPS <- COLUMN_IN_GROUPS

        if (all(!CONTRASTS_IN_GROUPS) && PASS) {
          check_return$e23 <- "All comparisons were invalid."
          PASS <- FALSE
        }

        if (any(!CONTRASTS_IN_GROUPS) && PASS) {
          check_return$e22 <- colnames(df_clean)[!CONTRASTS_IN_GROUPS]
          df_clean <- df_clean[, CONTRASTS_IN_GROUPS, drop = FALSE]
        }
      }
    } ## if PASS
  }

  # general checks for all data datatypes

  # check for empty df
  IS_DF_EMPTY <- any(dim(df_clean) == 0)

  if (!IS_DF_EMPTY) {
    df_clean <- as.matrix(df_clean)
  }

  if (IS_DF_EMPTY && PASS) {
    check_return$e15 <- "empty dataframe"
    PASS <- FALSE
  }

  return(
    list(
      df = df_clean,
      checks = check_return,
      PASS = PASS
    )
  )
}

#' Cross check input files for pgx.computePGX
#'
#' @param SAMPLE data.frame. The data frame corresponding to the input file as in playbase::SAMPLES
#' @param COUNTS data.frame. The data frame corresponding to the input file as in playbase::COUNTS
#' @param CONTRASTS data.frame. The data frame corresponding to the input file as in playbase::CONTRASTS
#'
#' @return a list with FIVE elements: SAMPLES, COUNTS and CONTRASTS that are the cleaned version of the
#'  input data frames, `checks` which contains the status of the checks, and
#'  `PASS` which contains the overall status of the check.
#' @export
pgx.crosscheckINPUT <- function(
    SAMPLES = NULL,
    COUNTS = NULL,
    CONTRASTS = NULL,
    PASS = TRUE) {
  samples <- SAMPLES
  counts <- COUNTS
  contrasts <- CONTRASTS
  PASS <- PASS
  check_return <- list()

  if (!is.null(samples) && !is.null(counts)) {
    # Check that rownames(samples) match colnames(counts)
    COUNTS_NAMES_NOT_MATCHING_SAMPLES <- colnames(counts)[!colnames(counts) %in% rownames(samples)]

    # if there are not matches between samples and counts, return error e25
    if (length(COUNTS_NAMES_NOT_MATCHING_SAMPLES) == length(colnames(counts)) && PASS) {
      check_return$e25 <- COUNTS_NAMES_NOT_MATCHING_SAMPLES
      PASS <- FALSE
    }

    if (length(COUNTS_NAMES_NOT_MATCHING_SAMPLES) > 0 && PASS) {
      check_return$e21 <- COUNTS_NAMES_NOT_MATCHING_SAMPLES
      PASS <- TRUE
      # align counts columns with samples rownames
      samples_in_counts <- rownames(samples)[rownames(samples) %in% colnames(counts)]
      counts <- counts[, samples_in_counts, drop = FALSE]
    }
    SAMPLE_NAMES_NOT_MATCHING_COUNTS <- rownames(samples)[!rownames(samples) %in% colnames(counts)]

    if (length(SAMPLE_NAMES_NOT_MATCHING_COUNTS) > 0 && PASS) {
      check_return$e16 <- SAMPLE_NAMES_NOT_MATCHING_COUNTS
      PASS <- TRUE
      # align samples rows with counts colnames
      counts_in_samples <- colnames(counts)[colnames(counts) %in% rownames(samples)]
      samples <- samples[counts_in_samples, , drop = FALSE]
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
      TOTAL_SAMPLES_NAMES <- unique(c(rownames(samples), colnames(counts)))
      check_return$e19 <- TOTAL_SAMPLES_NAMES[!TOTAL_SAMPLES_NAMES %in% SAMPLE_NAMES_PARTIAL_MATCHING_COUNTS]
      samples <- samples[SAMPLE_NAMES_PARTIAL_MATCHING_COUNTS, , drop = FALSE]
      counts <- counts[, SAMPLE_NAMES_PARTIAL_MATCHING_COUNTS, drop = FALSE]
    }

    # Check that counts have the same order as samples.

    MATCH_SAMPLES_COUNTS_ORDER <- all(diff(match(rownames(samples), colnames(counts))) > 0)

    # in case no matches are found, we get an NA, which should be converted to FALSE
    if (is.na(MATCH_SAMPLES_COUNTS_ORDER)) {
      MATCH_SAMPLES_COUNTS_ORDER <- FALSE
    }

    if (!MATCH_SAMPLES_COUNTS_ORDER && PASS) {
      check_return$e18 <- "We will reorder your samples and counts."
      counts <- counts[, match(rownames(samples), colnames(counts))]
    }
  }

  if (!is.null(samples) && !is.null(contrasts)) {
    # Conver contrasts and Check that rows names of contrasts match rownames of samples.
    contrasts_check_results <- contrasts_conversion_check(samples, contrasts, PASS)

    if (contrasts_check_results$PASS == FALSE && PASS) {
      PASS <- FALSE
      check_return$e20 <- rownames(contrasts)[!rownames(contrasts) %in% rownames(samples)]
    }

    contrasts <- contrasts_check_results$CONTRASTS

    # Check that rownames(samples) match long contrast rownames.

    SAMPLE_NAMES_NOT_MATCHING_CONTRASTS <- NULL

    if (dim(contrasts)[1] > dim(samples)[1] && PASS) { # check that contrasts are in long format
      SAMPLE_NAMES_NOT_MATCHING_CONTRASTS <- c(
        setdiff(rownames(samples), rownames(contrasts)),
        setdiff(rownames(contrasts), rownames(samples))
      )
    }
    if (length(SAMPLE_NAMES_NOT_MATCHING_CONTRASTS) > 0 && PASS) {
      check_return$e17 <- SAMPLE_NAMES_NOT_MATCHING_CONTRASTS
    }
  }

  return(
    list(
      SAMPLES = samples,
      COUNTS = counts,
      CONTRASTS = contrasts,
      checks = check_return,
      PASS = PASS
    )
  )
}

#' Convert contrasts for OPG
#'
#' @param SAMPLE data.frame. The data frame corresponding to the input file as in playbase::SAMPLES
#' @param CONTRASTS data.frame. The data frame corresponding to the input file as in playbase::CONTRASTS
#' @param PASS boolean. The status of the checks.
#' @return converted contrast df
#' @export
contrasts_conversion_check <- function(SAMPLES, CONTRASTS, PASS) {
  samples1 <- SAMPLES
  contrasts1 <- contrasts.convertToLabelMatrix(CONTRASTS, SAMPLES)

  if (is.null(contrasts1)) {
    message("[contrasts_conversion_check] WARNING: could not convert contrasts!")
    return(list(CONTRASTS = CONTRASTS, PASS = FALSE))
  }

  ok.contrast <- length(intersect(rownames(samples1), rownames(contrasts1))) > 0
  if (ok.contrast && NCOL(contrasts1) > 0 && PASS) {
    # check that dimentions of contrasts match samples
    if (dim(contrasts1)[1] != dim(samples1)[1] && PASS) {
      message("[contrasts_conversion_check] WARNING: numrows of contrast1 do not match samples!")
      return(list(CONTRASTS = contrasts1, PASS = FALSE))
    }
    rownames(contrasts1) <- rownames(samples1)
    for (i in 1:ncol(contrasts1)) {
      isz <- (contrasts1[, i] %in% c(NA, "NA", "NA ", "", " ", "  ", "   ", " NA"))
      if (length(isz)) contrasts1[isz, i] <- NA
    }
  }
  return(list(CONTRASTS = contrasts1, PASS = PASS))
}

contrasts_conversion_check.SAVE <- function(SAMPLES, CONTRASTS, PASS) {
  samples1 <- SAMPLES
  contrasts1 <- CONTRASTS
  PASS <- PASS ## ???

  group.col <- grep("group", tolower(colnames(samples1)))
  is.numeric.contrast <- all(as.vector(unlist(contrasts1)) %in% c(-1, 0, 1, NA))
  is.numeric.contrast
  ## old1: group-wise -1/0/1 matrix
  old1 <- (length(group.col) > 0 &&
    nrow(contrasts1) < nrow(samples1) &&
    all(rownames(contrasts1) %in% samples1[, group.col[1]]) &&
    is.numeric.contrast)
  ## old2: sample-wise -1/0/1 matrix
  old2 <- (nrow(contrasts1) == nrow(samples1) &&
    all(rownames(contrasts1) == rownames(samples1)) &&
    is.numeric.contrast)
  ## old3: group-wise label matrix
  old3 <- (length(group.col) > 0 &&
    nrow(contrasts1) < nrow(samples1) &&
    all(rownames(contrasts1) %in% samples1[, group.col[1]]) &&
    !is.numeric.contrast)

  old1
  old2
  old3

  old.style <- (old1 || old2 || old3)
  if (old.style && old1) {
    message("[contrasts_conversion_check] WARNING: converting old1 style contrast to new format")
    new.contrasts <- samples1[, 0]
    if (NCOL(contrasts1) > 0) {
      new.contrasts <- contrastAsLabels(contrasts1)
      grp <- as.character(samples1[, group.col])
      new.contrasts <- new.contrasts[grp, , drop = FALSE]
      rownames(new.contrasts) <- rownames(samples1)
    }
    contrasts1 <- new.contrasts
  }
  if (old.style && old2) {
    message("[contrasts_conversion_check] WARNING: converting old2 style contrast to new format")
    new.contrasts <- samples1[, 0]
    if (NCOL(contrasts1) > 0) {
      new.contrasts <- contrastAsLabels(contrasts1)
      rownames(new.contrasts) <- rownames(samples1)
    }
    contrasts1 <- new.contrasts
  }
  if (old.style && old3) {
    message("[contrasts_conversion_check] WARNING: converting group-wise label contrast to new format")
    new.contrasts <- samples1[, 0]
    if (NCOL(contrasts1) > 0) {
      grp <- as.character(samples1[, group.col])
      new.contrasts <- contrasts1[grp, , drop = FALSE]
      rownames(new.contrasts) <- rownames(samples1)
    }
    contrasts1 <- new.contrasts
  }

  ok.contrast <- length(intersect(rownames(samples1), rownames(contrasts1))) > 0
  if (ok.contrast && NCOL(contrasts1) > 0 && PASS) {
    ## always clean up
    contrasts1 <- apply(contrasts1, 2, as.character)

    # check that dimentions of contrasts match samples
    if (dim(contrasts1)[1] != dim(samples1)[1] && PASS) {
      message("[contrasts_conversion_check] WARNING: numrows of contrast1 do not match samples!")
      PASS <- FALSE
      return(list(CONTRASTS = contrasts1, PASS = PASS))
    }
    dbg("[contrasts_conversion_check] dim.CONTRASTS =", dim(CONTRASTS))
    dbg("[contrasts_conversion_check] dim.contrasts1 =", dim(contrasts1))
    dbg("[contrasts_conversion_check] dim.samples1 =", dim(samples1))
    rownames(contrasts1) <- rownames(samples1)
    for (i in 1:ncol(contrasts1)) {
      isz <- (contrasts1[, i] %in% c(NA, "NA", "NA ", "", " ", "  ", "   ", " NA"))
      if (length(isz)) contrasts1[isz, i] <- NA
    }
  }
  return(list(CONTRASTS = contrasts1, PASS = PASS))
}
