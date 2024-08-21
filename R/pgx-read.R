#' Read data file as matrix
#'
#' @param file Path to input data file
#' @param skip_row_check (default `FALSE`) Flag to skip the removal
#' of empty rows
#'
#' @return Matrix object containing data from file
#'
#' @description Reads a tabular data file and returns a matrix object.
#' Automatically detects separator and allows duplicate row names.
#'
#' @details This function reads a tabular text file, automatically detecting the
#' separator (tab, comma, semicolon). It returns a matrix containing the data values,
#' using the first column as rownames (allowing duplicates). Blank rows and rows with
#' NA as rowname are skipped.
#'
#' @examples
#' \dontrun{
#' mymatrix <- read.as_matrix(mydata.csv)
#' }
#' @export
read.as_matrix <- function(file, skip_row_check = FALSE, row.names = 1) {
  ## determine if there are empty lines in header
  x0 <- data.table::fread(
    file = file,
    header = FALSE,
    nrow = 100
  )
  x0[is.na(x0)] <- ""
  skip.rows <- min(which(cumsum(rowMeans(x0 != "")) > 0)) - 1

  ## try to detect decimal separator
  dec <- detect_decimal(file)

  ## read delimited table automatically determine separator. allow
  ## duplicated rownames. This implements with faster fread.
  x0 <- data.table::fread(
    file = file,
    check.names = FALSE,
    header = TRUE,
    dec = dec,
    fill = TRUE,
    skip = skip.rows,
    blank.lines.skip = TRUE,
    stringsAsFactors = FALSE,
    integer64 = "numeric"
  )

  ## For csv with missing rownames field at (1,1) in the header,
  ## fill=TRUE will fail. Check header with slow read.csv() and
  ## correct if needed. fread is fast but is not so robust...
  ## detect delimiter/seperator
  sep <- detect_delim(file)
  hdr <- utils::read.csv(
    file = file, sep = sep, check.names = FALSE, na.strings = NULL,
    header = TRUE, nrows = 1, skip = skip.rows, row.names = NULL
  )
  if (NCOL(x0) > 0 && !all(colnames(x0) == colnames(hdr))) {
    message("read.as_matrix: warning correcting header")
    colnames(x0) <- colnames(hdr)
  }

  x <- NULL
  ## drop rows without rownames
  sel <- which(!as.character(x0[[1]]) %in% c("", " ", "NA", "na", NA))
  if (length(sel)) {
    ##
    x <- x0[sel, , drop = FALSE]
    # convert x from data.table to matrix. as.matrix means we do not
    # have mixed types (such as in dataframes). This is needed as
    # fread reads numeric as integer64, which is not supported by
    # matrix at this stage since we handle samples, counts and
    # contrasts, it is dangerous to force conversion to numeric hence,
    # the conversion from character to numeric for counts is handled
    # at pgx.checkINPUT
    colnames_x <- colnames(x)
    x <- sapply(1:ncol(x), function(i) {
      as.character(x[[i]]) ## always character!
    })
    colnames(x) <- colnames_x
  } else {
    return(NULL)
  }

  ## for character matrix, we strip whitespace
  if (is.character(x)) {
    x <- trimws(x)
  }

  ## set rownames
  if (!is.null(row.names) && !is.na(row.names) && row.names >= 1) {
    x <- x[, -row.names, drop = FALSE]
    rownames(x) <- x0[[row.names]]
  }

  ## some csv have trailing empty rows/cols at end of table
  if (NCOL(x) && !skip_row_check) { # bypass in case full NA rows
    empty.row <- (rowSums(is.na(x)) == ncol(x))
    if (tail(empty.row, 1)) {
      n <- which(!rev(empty.row))[1] - 1
      ii <- (nrow(x) - n + 1):nrow(x)
      x <- x[-ii, , drop = FALSE]
    }
    ## some csv have trailing empty columns at end of table
    empty.col <- (colSums(is.na(x)) == nrow(x))
    if (tail(empty.col, 1)) {
      n <- which(!rev(empty.col))[1] - 1
      ii <- (ncol(x) - n + 1):ncol(x)
      x <- x[, -ii, drop = FALSE]
    }
  }
  return(x)
}

#' Detect delimiter of text file from header (or first line)
#'
detect_delim <- function(file) {
  hdr <- readLines(file, n = 1)
  n_commas <- length(setdiff(gregexpr(",", hdr, fixed = TRUE)[[1]], -1))
  n_semicolons <- length(setdiff(gregexpr(";", hdr, fixed = TRUE)[[1]], -1))
  n_tabs <- length(setdiff(gregexpr("\t", hdr, fixed = TRUE)[[1]], -1))
  delim <- c(",", ";", "\t")[which.max(c(n_commas, n_semicolons, n_tabs))]
  delim
}

#' Detect delimiter of text file from first 10 lines. Assumes there is
#' a header and rownames column.
#'
detect_decimal <- function(file) {
  ff <- data.table::fread(file, header = TRUE, colClasses = "character")
  vals <- as.vector(as.matrix(ff[, -1]))
  n_commas <- length(grep(",", vals, fixed = TRUE))
  n_dots <- length(grep(".", vals, fixed = TRUE))
  dec <- c(".", ",")[which.max(c(n_dots, n_commas))]
  dec
}

#' Read CSV file into R efficiently
#'
#' @param file Path to CSV file
#' @param check.names Logical, should column names be checked for syntactic validity. Default is FALSE.
#' @param row.names Column to use for row names, default is 1 (first column).
#' @param sep Separator character, default is "auto" for automatic detection.
#' @param stringsAsFactors Logical, should character columns be converted to factors? Default is FALSE.
#' @param header Logical, does the file have a header row? Default is TRUE.
#' @param asMatrix Logical, should the result be returned as a matrix instead of a data frame? Default is TRUE.
#'
#' @return A data frame or matrix containing the parsed CSV data.
#'
#' @details This function efficiently reads a CSV file into R using \code{data.table::fread()}, then converts it into a regular data frame or matrix.
#' It is faster than \code{read.csv()} especially for large files.
#'
#' By default it converts the result to a matrix if all columns are numeric, character or integer. The row names are taken from the first column.
#' Factor conversion, column type checking, and header parsing can be controlled with parameters.
#'
#' @examples
#' \dontrun{
#' dat <- fread.csv("data.csv")
#' }
#' @export
fread.csv <- function(file, check.names = FALSE, row.names = 1, sep = ",",
                      stringsAsFactors = FALSE, header = TRUE, asMatrix = TRUE) {
  df <- data.table::fread(
    file = file, check.names = check.names, header = header, sep = sep, fill = TRUE
  )
  if (NCOL(df) == 1) {
    ## empty file, only rownames
    x <- matrix(NA, nrow(df), 0)
    rownames(x) <- df[[row.names]] ## allow dups if matrix
    return(x)
  }
  n0 <- ifelse(row.names == 0 || is.null(row.names), 1, 2)
  colnames(df) <- substring(colnames(df), 1, 1000) ## safety, avoid length overflow
  x <- data.frame(df[, n0:ncol(df)],
    stringsAsFactors = stringsAsFactors,
    check.names = check.names
  )
  ## check&correct for truncated header
  if (row.names == 0) {
    rn <- NULL
  } else {
    rn <- 1
  }
  hdr <- colnames(read.csv(file,
    nrow = 1, sep = sep, header = TRUE,
    row.names = rn, check.names = check.names
  ))
  if (!all(colnames(x) == hdr)) {
    colnames(x) <- hdr
  }
  is.num <- all(sapply(x, class) == "numeric")
  is.char <- all(sapply(x, class) == "character")
  is.int <- all(sapply(x, class) == "integer")
  if (asMatrix && (is.num || is.char || is.int)) x <- as.matrix(x)
  if (!is.null(row.names) && row.names != 0) {
    rownames(x) <- df[[row.names]] ## allow dups if matrix
  }
  return(x)
}

#' Read all CSV file in folder or current directory
#'
#' @export
read_files <- function(dir = ".", pattern = NULL) {
  ff <- dir(dir, pattern = pattern)
  file1 <- head(grep("count", ff, value = TRUE), 1)
  file2 <- head(grep("sample", ff, value = TRUE), 1)
  file3 <- head(grep("contrast|comparison", ff, value = TRUE), 1)
  counts <- samples <- contrasts <- NULL
  if (length(file1)) counts <- read_counts(file1)
  if (length(file2)) samples <- read_samples(file2)
  if (length(file3)) contrasts <- read_contrasts(file3)

  list(
    counts = counts,
    samples = samples,
    contrasts = contrasts
  )
}



#' Read counts data from file
#'
#' @param file string. path to file
#' @param drop_na_rows boolean. drop rows without rownames
#' @param first boolean. drop multiple feature names (separated by ;)
#' @param unique boolean. make duplicated rows unique by pasting a number
#'
#' @details This function reads a count matrix using \code{playbase::read.as_matrix()},
#' validates it with \code{validate_counts()}, and optionally converts row names from IDs to
#' gene symbols using \code{playbase::probe2symbol()}.
#'
#' It removes rows with NA, blank or invalid symbols, and collapses any duplicate symbols by
#' summing counts across rows.
#'
#' @examples
#' \dontrun{
#' counts <- read_counts(playbase::example_file("counts.csv"))
#' }
#' @export
read_counts <- function(file, first = FALSE, unique = TRUE, paste_char = "_") {
  if (is.character(file)) {
    df <- read.as_matrix(file)
  } else if (is.matrix(file) || is.data.frame(file)) {
    df <- file
  } else {
    stop("input error")
  }
  is_valid <- validate_counts(df)
  if (!is_valid) stop("Counts file is not valid.")

  ## determine column types
  df1 <- type.convert(data.frame(df, check.names = FALSE), as.is = TRUE)
  col.type <- sapply(df1, class)
  col.type
  char.cols <- which(col.type == "character")
  last.charcol <- tail(char.cols, 1)
  last.charcol

  ## if the rownames are not unique, and some more character columns
  ## exists, then search for best column and paste after rowname.
  if (length(char.cols) > 0 && sum(duplicated(rownames(df)))) {
    ndup <- sapply(char.cols, function(k) {
      sum(duplicated(paste0(rownames(df), "_", df[, k])))
    })
    ndup
    sel <- names(which.min(ndup))
    rownames(df) <- paste0(rownames(df), "_", df[, sel])
  }

  ## As expression values we take all columns after the last character
  ## column (if any).
  if (length(last.charcol)) {
    df <- df[, (last.charcol + 1):ncol(df)]
  }

  ## convert to numeric if needed (probably yes...)
  is.numeric.matrix <- all(apply(df, 2, is.numeric))
  if (!is.numeric.matrix) {
    message("warning: converting to numeric values")
    rn <- rownames(df)
    suppressWarnings(df <- apply(df, 2, as.numeric))
    rownames(df) <- rn
  }

  ## when multiple feature names, should we take only first feature?
  if (first) rownames(df) <- first_feature(rownames(df))

  ## if rownames are duplicated, we append a number behing
  if (unique) rownames(df) <- make_unique(rownames(df))
  return(df)
}

#' Read samples data from file
#'
#' @param file Path to input sample data file. Should be a matrix with samples as rows and metadata as columns.
#'
#' @return dataframe the file with the data
#' @details This function reads the sample matrix with the meta-data information
#' of the counts and converts it to a dataframe.
#'
#' @examples
#' samples <- read_samples(playbase::example_file("samples.csv"))
#' @export
read_samples <- function(file) {
  df <- read.as_matrix(file)
  is_valid <- validate_samples(df)
  if (!is_valid) stop("Samples file is not valid.")
  df <- as.data.frame(df)
  return(df)
}

#' Read contrasts data from file
#'
#' @param file string. path to file
#'
#' @return matrix. the file with the data
#'
#' @examples
#' contrasts <- read_contrasts(playbase::example_file("contrasts.csv"))
#' @export
read_contrasts <- function(file) {
  df <- read.as_matrix(file)

  is_valid <- validate_contrasts(df)
  if (!is_valid) stop("Contrasts file is not valid.")

  df
}

#' Read gene/probe annotation file
#'
#' @export
read_annot <- function(file, unique = TRUE) {
  if (is.character(file)) {
    ## we read without rownames because we want to retain full header
    df <- read.as_matrix(file, row.names = NULL)
    rownames(df) <- df[, 1]
  } else if (is.matrix(file) || is.data.frame(file)) {
    df <- file
  } else {
    stop("input error")
  }

  ## add column title
  if (colnames(df)[1] == "") colnames(df)[1] <- "row.names"

  ## determine last character column
  df1 <- type.convert(data.frame(df, check.names = FALSE), as.is = TRUE)
  col.type <- sapply(df1, class)
  col.type
  char.cols <- which(col.type == "character")
  char.cols
  last.charcol <- tail(char.cols, 1)
  last.charcol

  ## if the rownames are not unique, and some more character columns
  ## exists, then search for best column and paste after rownames
  ## (first column).
  if (sum(duplicated(df[, 1])) && length(char.cols) >= 2) {
    ndup <- sapply(char.cols[-1], function(k) {
      sum(duplicated(paste0(df[, 1], "_", df[, k])))
    })
    ndup
    sel <- names(which.min(ndup))[1]
    rownames(df) <- paste0(rownames(df), "_", df[, sel])
  }

  ## drop numerical columns (these can be intensities). We check if we
  ## have equal or more than two columns. First column are rownames.
  if (length(char.cols) && last.charcol >= 2) {
    df <- df[, 1:last.charcol, drop = FALSE]
    df <- data.frame(df, check.names = FALSE)
  } else {
    df <- NULL
  }
  return(df)
}

#' Validate counts data
#'
#' Counts data is valid if:
#'  - no duplicate rows
#'  - no empty rows
#'  - no duplicate cols
#'
#' @param data matrix.
#'
#' @return boolean. true if data is valid
#' @export
validate_counts <- function(data) {
  t1 <- check_duplicate_rows(data)
  t2 <- check_empty_rows(data)
  t3 <- check_duplicate_cols(data)
  return(t2 && t3)
}

#' Validate samples data
#'
#' Samples data is valid if:
#'  - no duplicate rows
#'  - no empty rows
#'  - no duplicate cols
#'  - contains less than the max samples allowed
#'
#' @param data matrix.
#'
#' @return boolean. true if data is valid
#' @export
validate_samples <- function(data) {
  t1 <- check_duplicate_rows(data)
  t2 <- check_empty_rows(data)
  t3 <- check_duplicate_cols(data)
  t4 <- check_max_samples(data)
  return(t1 & t2 & t3 & t4)
}

#' Validate contrasts data
#'
#' Contrasts data is valid if:
#'  - no duplicate rows
#'  - no empty rows
#'  - no duplicate cols
#'  - contains only cols with "_vs_" in names
#'
#' @param data matrix.
#'
#' @return boolean. true if data is valid
#' @export
validate_contrasts <- function(data) {
  t1 <- check_duplicate_rows(data)
  t2 <- check_empty_rows(data)
  t3 <- check_duplicate_cols(data)
  vs_names <- colnames(data)[-1]
  t4 <- all(grepl("_vs_", vs_names))
  return(t1 & t2 & t3 & t4)
}


#' @describeIn check_duplicate_cols check if there is any duplicate row in the input data
#' @export
check_duplicate_rows <- function(data) {
  rn <- setdiff(data[[1]], c("", "NA", NA))
  t1 <- sum(duplicated(rn)) == 0
  return(t1)
}


#' @describeIn check_duplicate_cols checks if there is any empty row in the data
#' @export
check_empty_rows <- function(data) {
  t1 <- nrow(data) > 0
  return(t1)
}


#' @title Input Checks
#'
#' @param data A data frame or matrix.
#'
#' @return Logical indicating output of the checks.
#'
#' @description Checks if a data frame or matrix contains duplicate column names.
#'
#' @details This function takes a data frame or matrix \code{data} as input and checks if it contains any duplicate column names.
#' It compares the column names against each other to look for duplicates.
#'
#' The output is a logical value indicating whether any duplicate names were found.
#' \code{TRUE} means duplicate names were detected, \code{FALSE} means no duplicates.
#'
#' This can be used to validate data before further analysis, to ensure no columns are duplicated.
#'
#' @examples
#' \dontrun{
#' data <- data.frame(A = 1:3, B = 4:6, A = 7:9)
#' check_duplicate_cols(data)
#' # Returns TRUE
#'
#' data <- data.frame(A = 1:3, B = 4:6, C = 7:9)
#' check_duplicate_cols(data)
#' # Returns FALSE
#' }
#' @export
check_duplicate_cols <- function(data) {
  t1 <- sum(duplicated(colnames(data))) == 0
  return(t1)
}


#' @describeIn check_duplicate_cols Checks if the number of sample is below the allowed maximum
#' @export
check_max_samples <- function(data, max_samples = 2000) {
  MAXSAMPLES <- as.integer(max_samples)
  t1 <- (ncol(data) - 1) <= MAXSAMPLES
  return(t1)
}
