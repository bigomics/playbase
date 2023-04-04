
#' Read counts data from file
#'
#' @param file string. path to file
#' @param convert_names boolean.
#'
#' @return matrix. the file with the data
#' @export
#'
#' @examples
#' counts <- read_counts(playbase::example_file("counts.csv"))
read_counts <- function(file, convert_names = FALSE) {
  df <- read_as_matrix(file)

  is_valid <- validate_counts(df)
  if (!is_valid) stop("Counts file is not valid.")

  # convert to gene names (needed for biological effects)
  if (convert_names) {
    pp <- rownames(df)
    rownames(df) <- probe_to_symbol(pp)
    sel <- !(rownames(df) %in% c(NA, "", "NA"))
    df <- df[sel, ]
    xx <- tapply(
      1:nrow(df), rownames(df),
      function(i) colSums(df[i, , drop = FALSE])
    )
    df <- do.call(rbind, xx)
  }

  df
}

#' Read expression data from file
#'
#' @param file string. path to file
#' @param convert_names boolean.
#'
#' @return matrix. the file with the data
#' @export
#'
#' @examples
#' counts <- read_expression(playbase::example_file("counts.csv"))
read_expression <- function(file, convert_names = TRUE) {
  df <- read_as_matrix(file)
  df <- df**2

  is_valid <- validate_counts(df)
  if (!is_valid) stop("Expression file is not valid.")

  # convert to gene names (needed for biological effects)
  if (convert_names) {
    pp <- rownames(df)
    rownames(df) <- probe_to_symbol(pp)
    sel <- !(rownames(df) %in% c(NA, "", "NA"))
    df <- df[sel, ]
    xx <- tapply(
      1:nrow(df), rownames(df),
      function(i) colSums(df[i, , drop = FALSE])
    )
    df <- do.call(rbind, xx)
  }

  df
}

#' Read samples data from file
#'
#' @param file string. path to file
#'
#' @return dataframe the file with the data
#' @export
#'
#' @examples
#' samples <- read_samples(playbase::example_file("samples.csv"))
read_samples <- function(file) {
  df <- read_as_matrix(file)

  is_valid <- validate_samples(df)
  if (!is_valid) stop("Samples file is not valid.")

  as.data.frame(df)
}

#' Read contrasts data from file
#'
#' @param file string. path to file
#'
#' @return matrix. the file with the data
#' @export
#'
#' @examples
#' contrasts <- read_contrasts(playbase::example_file("contrasts.csv"))
read_contrasts <- function(file) {
  df <- read_as_matrix(file)

  is_valid <- validate_contrasts(df)
  if (!is_valid) stop("Contrasts file is not valid.")

  df
}

#' Read a file as a matrix
#'
#' This function is used to read counts/samples/contrasts data from file.
#'
#' @param file string. filepath to data
#'
#' @return matrix
#'
#' @examples
#' x <- 1
read_as_matrix <- function(file) {
  x0 <- data.table::fread(
    file = file,
    check.names = FALSE,
    header = TRUE,
    blank.lines.skip = TRUE,
    stringsAsFactors = FALSE
  )
  x <- NULL
  sel <- which(!as.character(x0[[1]]) %in% c("", " ", "NA", "na", NA))
  if (length(sel)) {
    x <- as.matrix(x0[sel, -1, drop = FALSE])
    rownames(x) <- x0[[1]][sel]
  }
  return(x)
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
validate_contrasts <- function(data) {
  t1 <- check_duplicate_rows(data)
  t2 <- check_empty_rows(data)
  t3 <- check_duplicate_cols(data)
  vs_names <- colnames(data)[-1]
  t4 <- all(grepl("_vs_", vs_names))
  return(t1 & t2 & t3 & t4)
}

check_duplicate_rows <- function(data) {
  rn <- setdiff(data[[1]], c("", "NA", NA))
  t1 <- sum(duplicated(rn)) == 0
  return(t1)
}

check_empty_rows <- function(data) {
  t1 <- nrow(data) > 0
  return(t1)
}

check_duplicate_cols <- function(data) {
  t1 <- sum(duplicated(colnames(data))) == 0
  return(t1)
}

check_max_samples <- function(data, max_samples = 2000) {
  MAXSAMPLES <- as.integer(max_samples)
  t1 <- (ncol(data) - 1) <= MAXSAMPLES
  return(t1)
}
