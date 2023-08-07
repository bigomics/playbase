#' Read counts data from file
#'
#' @param file string. path to file
#' @param convert_names boolean.
#' @details This function reads a count matrix using \code{playbase::read.as_matrix()},
#' validates it with \code{validate_counts()}, and optionally converts row names from IDs to
#' gene symbols using \code{playbase::probe2symbol()}.
#'
#' It removes rows with NA, blank or invalid symbols, and collapses any duplicate symbols by
#' summing counts across rows.
#'
#' @examples
#' counts <- read_counts(playbase::example_file("counts.csv"))
#' @export
read_counts <- function(file, convert_names = FALSE) {
  df <- read.as_matrix(file)

  is_valid <- validate_counts(df)
  if (!is_valid) stop("Counts file is not valid.")

  # convert to gene names (needed for biological effects)
  if (convert_names) {
    pp <- rownames(df)
    rownames(df) <- probe2symbol(pp)
    sel <- !(rownames(df) %in% c(NA, "", "NA"))
    df <- df[sel, ]
    xx <- tapply(
      1:nrow(df), rownames(df),
      function(i) colSums(df[i, , drop = FALSE])
    )
    df <- do.call(rbind, xx)
  }

  return(df)
}

#' Read expression data from file
#'
#' Reads an expression data matrix from file and performs validation and preprocessing.
#'
#' @param file Path to input expression data file. Should be matrix with genes as rows and samples as columns.
#' @param convert_names Logical indicating whether to convert row names from IDs to gene symbols. Default is TRUE.
##'
#' @details This function reads an expression matrix using playbase::read.as_matrix(),
#' squares the values, validates with validate_counts(), and optionally converts row names
#' from IDs to gene symbols using playbase::probe2symbol().
#'
#' It removes rows with NA, blank or invalid symbols, and collapses duplicate symbols by
#' summing expression across rows.
#' @return matrix. the file with the data
#'
#' @examples
#' counts <- read_expression(playbase::example_file("counts.csv"))
#' @export
read_expression <- function(file, convert_names = TRUE) {
  df <- read.as_matrix(file)
  df <- df**2

  is_valid <- validate_counts(df)
  if (!is_valid) stop("Expression file is not valid.")

  # convert to gene names (needed for biological effects)
  if (convert_names) {
    pp <- rownames(df)
    rownames(df) <- probe2symbol(pp)
    sel <- !(rownames(df) %in% c(NA, "", "NA"))
    df <- df[sel, ]
    xx <- tapply(
      1:nrow(df), rownames(df),
      function(i) colSums(df[i, , drop = FALSE])
    )
    df <- do.call(rbind, xx)
  }

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
