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
read.as_matrix <- function(file, skip_row_check = FALSE, as.char = TRUE,
                           as.matrix = TRUE, row.names = 1) {
  ## determine if there are empty lines in header
  x0 <- data.table::fread(
    file = file,
    header = FALSE,
    nrow = 100
  )
  x0[is.na(x0)] <- ""
  skip.rows <- min(which(cumsum(rowMeans(x0 != "")) > 0)) - 1
  skip.colsN <- min(which(cumsum(rev(colMeans(x0 != ""))) > 0)) - 1
  keep.cols <- seq_len(ncol(x0) - skip.colsN)

  ## try to detect decimal separator
  sep <- detect_delim(file)
  dec <- detect_decimal(file)
  if (dec == "," && sep == ",") dec <- "." ## exception

  ## read delimited table automatically determine separator. allow
  ## duplicated rownames. This implements with faster fread.
  x0 <- data.table::fread(
    file = file,
    check.names = FALSE,
    header = TRUE,
    dec = dec,
    # fill = TRUE, ## fill=TRUE will fail for some datasets
    skip = skip.rows,
    select = keep.cols,
    blank.lines.skip = TRUE,
    stringsAsFactors = FALSE,
    integer64 = "numeric"
  )
  # Handle weird case where file contains quotes on each row
  # e.g. "field1 field2"
  #      "val1   val2"
  # where columns are read ok by fread
  # but first and last column maintain the "
  first_column <- x0[[1]] # Extract the first column
  last_column <- x0[[ncol(x0)]] # Extract the last column
  first_column <- iconv2utf8(first_column)
  last_column <- iconv2utf8(last_column)
  if (all(grepl('^"', first_column)) && all(grepl('"$', last_column))) {
    x0[[1]] <- gsub('^"', "", first_column)
    x0[[ncol(x0)]] <- as.numeric(gsub('"$', "", last_column))
  }

  ## see: https://github.com/Rdatatable/data.table/issues/2607
  coln.int64 <- names(which(sapply(x0, bit64::is.integer64)))
  if (length(coln.int64) > 0L) {
    x0[, c(coln.int64) := lapply(.SD, as.numeric), .SDcols = coln.int64]
  }

  ## For csv with missing rownames field at (1,1) in the header,
  ## fill=TRUE will fail. Check header with slow read.csv() and
  ## correct if needed. fread is fast but is not so robust...
  ## detect delimiter/seperator
  hdr <- utils::read.csv(
    file = file, sep = sep, check.names = FALSE, na.strings = NULL,
    header = TRUE, nrows = 1, skip = skip.rows, row.names = NULL
  )
  if (ncol(hdr) > 1) { # In some rare cases, hdr can be read as a single column (wrong), skip correction when that happens
    hdr <- hdr[, keep.cols]
    if (NCOL(x0) > 0 && !all(colnames(x0) == colnames(hdr))) {
      message("[read.as_matrix] correcting header")
      colnames(x0) <- colnames(hdr)
    }
  }

  ## give first column name
  if (colnames(x0)[1] == "") colnames(x0)[1] <- "row.names"

  ## drop rows without rownames
  sel <- which(!as.character(x0[[1]]) %in% c("", " ", "NA", "na", "-", NA))
  if (length(sel) == 0) {
    return(NULL)
  }
  x <- x0[sel, , drop = FALSE]

  # Convert x from data.table to matrix. With as.char = TRUE,
  # as.matrix() does not return mixed types (such as in dataframes).
  if (as.char) {
    ## drop duplicated columns
    ## otherwise as.char will crash
    x <- x[, which(duplicated(colnames(x))) := NULL]
    colnames_x <- colnames(x)[-1]
    x[, c(colnames_x) := lapply(.SD, as.character), .SDcols = colnames_x]
  }

  ## for character matrix, we strip whitespace
  which.char <- which(sapply(x, class) == "character")
  if (length(which.char)) {
    char.cols <- colnames(x)[which.char]
    x[, c(char.cols) := lapply(.SD, iconv2utf8), .SDcols = char.cols]
    x[, c(char.cols) := lapply(.SD, trimws), .SDcols = char.cols]
  }

  ## set rownames, convert to matrix.
  if (as.matrix) {
    ## this allow duplicated rownames
    x <- as.matrix(x, rownames = row.names) ## this can be slow!!
  } else {
    if (!is.null(row.names)) {
      rownamesx <- x[[1]]
      data.table::set(x, j = as.integer(row.names), value = NULL)
      ## duplicated rownames not allowed!
      if (sum(duplicated(rownamesx)) > 0) {
        message("warning: duplicated rownames will be made unique")
      }
      rownames(x) <- make_unique(rownamesx)
    }
  }

  ## some csv have trailing empty rows/cols at end of table
  last.row.empty <- mean(is.na(x[nrow(x), ])) == 1
  last.col.empty <- mean(is.na(x[, ncol(x)])) == 1

  if (last.row.empty && !skip_row_check) { # bypass in case full NA rows
    empty.row <- (rowSums(is.na(x) | x %in% c("", NA, "NA", " ")) == ncol(x))
    empty.row <- empty.row & rownames(x) %in% c(NA, "", " ")
    if (tail(empty.row, 1)) {
      n <- which(!rev(empty.row))[1] - 1
      ii <- (nrow(x) - n + 1):nrow(x)
      x <- x[-ii, , drop = FALSE]
    }
  }
  if (last.col.empty && !skip_row_check) { # bypass in case full NA rows
    ## some csv have trailing empty columns at end of table
    ## detect as empty column, columns with NA and weird characters
    empty.col <- (colSums(is.na(x) | x %in% c("", "-", ".", "NA", " ")) == nrow(x))
    empty.col <- empty.col & colnames(x) %in% c(NA, "", " ")
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
detect_delim <- function(file, delims = c(",", "\t", " ", "|", ":", ";")) {
  # Code extracted from vroom:::guess_delim (version 1.5.7)
  lines <- readLines(file, n = 10)

  # blank text within quotes
  lines <- gsub('"[^"]*"', "", lines)

  splits <- lapply(delims, strsplit, x = lines, useBytes = TRUE, fixed = TRUE)
  counts <- lapply(splits, function(x) table(lengths(x)))
  num_fields <- vapply(counts, function(x) as.integer(names(x)[[1]]), integer(1))
  num_lines <- vapply(counts, function(x) (x)[[1]], integer(1))
  top_lines <- 0
  top_idx <- 0
  for (i in seq_along(delims)) {
    if (num_fields[[i]] >= 2 && num_lines[[i]] > top_lines ||
      (top_lines == num_lines[[i]] && (top_idx <= 0 || num_fields[[top_idx]] < num_fields[[i]]))) {
      top_lines <- num_lines[[i]]
      top_idx <- i
    }
  }
  if (top_idx == 0) {
    return(",") # default to comma
  }

  delims[[top_idx]]
}

#' Detect delimiter of text file from first 10 lines. Assumes there is
#' a header and rownames column.
#'
detect_decimal <- function(file) {
  f1 <- data.table::fread(file, header = TRUE, nrows = 100)
  f2 <- data.table::fread(file, header = TRUE, colClasses = "character", nrows = 100)
  numcols <- which(sapply(f1, class) == "numeric")
  if (length(numcols) == 0) {
    return(".")
  } ## default
  numvals <- as.vector(as.matrix(f2[, ..numcols]))
  n_commas <- length(grep(",", numvals, fixed = TRUE))
  n_dots <- length(grep(".", numvals, fixed = TRUE))
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
  file1 <- head(grep("count|expression|abundance|concentration|intensity",
    ff,
    value = TRUE
  ), 1)
  file2 <- head(grep("sample", ff, value = TRUE), 1)
  file3 <- head(grep("contrast|comparison", ff, value = TRUE), 1)
  counts <- samples <- contrasts <- NULL
  if (length(file1)) counts <- read_counts(file.path(dir, file1))
  if (length(file2)) samples <- read_samples(file.path(dir, file2))
  if (length(file3)) contrasts <- read_contrasts(file.path(dir, file3))

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
read_counts <- function(file, first = FALSE, unique = FALSE, paste_char = "_") {
  if (is.character(file)) {
    df <- read.as_matrix(file, as.char = FALSE)
  } else if (is.matrix(file) || is.data.frame(file)) {
    df <- file
  } else {
    message("[read_counts] ERROR: Input must be filename or matrix")
    return(NULL)
  }

  ## is_valid <- validate_counts(df)
  ## if (!is_valid) {
  ##   message("[read_counts] WARNING: Counts file has errors")
  ##   ## return(NULL)
  ## }

  ## determine column types (NEED RETHINK!)
  df1 <- type.convert(data.frame(head(df, 20), check.names = FALSE), as.is = TRUE)
  col.type <- sapply(df1, class)
  xannot.names <- "gene|^id$$|metabolite|compound|position|phospo.*site" ## possible numeric annotations
  is.xannot <- grepl(xannot.names, tolower(colnames(df)))
  char.cols <- which(col.type == "character" | is.xannot)
  last.charcol <- tail(char.cols, 1)
  last.charcol

  ## NEED RETHINK. if the rownames are not unique, and some more
  ## character columns exists, then search for best column and paste
  ## after rownames.
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
    message("[read_counts] extra annotation columns = ", paste(1:last.charcol, collapse = " "))
    df <- df[, (last.charcol + 1):ncol(df)]
  }

  ## convert to numeric if needed (probably yes...)
  is.numeric.matrix <- all(apply(head(df, 20), 2, is.numeric))
  if (!is.numeric.matrix) {
    message("[read_counts] force to numeric values")
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
  if (!is_valid) message("[read_samples] WARNING: Samples file has errors")
  return(as.data.frame(df))
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
  if (!is_valid) message("[read_contrasts] WARNING: Contrasts file has errors")
  return(df)
}


#' Read Olink NPX data and create a counts matrix
#' @param NPX_data Path to input Olink NPX data file. Must be standard Olink NPX format as per OlinkAnalyze R package.
#' @return data matrix (features on rows; samples on columns)
#' @export
read_Olink_NPX <- function(NPX_data) {
  NPX <- try(OlinkAnalyze::read_NPX(NPX_data), silent = TRUE)
  if (!inherits(NPX, "try-error")) {
    NPX <- as.data.frame(NPX)
  } else {
    dbg("[read_Olink_NPX]. Uploaded proteomics data file is not Olink or does not adhere with standard Olink NPX format.")
    return(NULL)
  }

  hh1 <- grep("NPX", colnames(NPX), ignore.case = TRUE)
  if (any(hh1)) {
    npx.id <- colnames(NPX)[hh1[1]]
  } else {
    dbg("[read_Olink_NPX] The uploaded Olink NPX file does not contain the variable NPX")
    return(NULL)
  }

  feature.id <- NULL
  hh <- grepl("uniprot", tolower(colnames(NPX)))
  if (any(hh)) feature.id <- colnames(NPX)[hh][1]
  if (is.null(feature.id)) {
    hh <- grepl("assay$", tolower(colnames(NPX)))
    if (any(hh)) feature.id <- colnames(NPX)[hh][1]
  }
  if (is.null(feature.id)) {
    dbg("[read_Olink_NPX] The uploaded Olink NPX file does not contain the variables Uniprot or Assay")
    return(NULL)
  }

  sample.id <- NULL
  hh <- grepl("sampleid", tolower(colnames(NPX)))
  if (any(hh)) sample.id <- colnames(NPX)[hh][1]
  if (is.null(sample.id)) {
    dbg("[read_Olink_NPX] The uploaded Olink NPX file does not contain the variable SampleID")
    return(NULL)
  }

  # Convert to matrix
  fm <- as.formula(paste0(feature.id, "~", sample.id))
  counts.df <- reshape2::dcast(NPX, fm, value.var = npx.id, fun.aggregate = mean)
  counts <- as.matrix(counts.df[sapply(counts.df, is.numeric)])
  rownames(counts) <- counts.df[, feature.id]

  return(counts)
}

#' Read Olink NPX data and automatically create samples dataframe
#' @param NPX_data Path to input Olink NPX data file. Must be standard Olink NPX format as per OlinkAnalyze R package.
#' @return dataframe with the metadata
#' @export
read_Olink_samples <- function(NPX_data) {
  NPX <- try(OlinkAnalyze::read_NPX(NPX_data), silent = TRUE)
  if (!inherits(NPX, "try-error")) {
    NPX <- as.data.frame(NPX)
  } else {
    dbg("[read_Olink_NPX]. Uploaded proteomics data file is not Olink or does not adhere with standard Olink NPX format.")
    return(NULL)
  }

  samples.id <- NULL
  hh <- grep("SampleID", colnames(NPX), ignore.case = TRUE)
  if (any(hh)) samples.id <- unique(NPX[, hh[1]])
  if (is.null(samples.id)) {
    dbg("[read_Olink_NPX] The uploaded Olink NPX file does not contain the variable SampleID")
    return(NULL)
  }

  # https://cran.r-project.org/web/packages/OlinkAnalyze/vignettes/Vignett.html
  exclude.vars <- ""
  hh <- grepl("uniprot|olinkid|assay|npx|freq|lod", tolower(colnames(NPX)))
  if (any(hh)) exclude.vars <- colnames(NPX)[hh]
  samples <- matrix(NA, nrow = length(samples.id), ncol = sum(!hh))
  samples <- data.frame(samples, row.names = samples.id)
  colnames(samples) <- colnames(NPX)[!hh]

  hh <- grep("SampleID", colnames(NPX), ignore.case = TRUE)
  i <- 1
  for (i in 1:nrow(samples)) {
    t <- 1
    for (t in 1:ncol(samples)) {
      jj <- match(rownames(samples)[i], NPX[, hh[1]])
      jx <- match(colnames(samples)[t], colnames(NPX))
      if (!any(jj) | !any(jx)) next
      samples[i, t] <- NPX[jj, jx]
    }
  }
  samples <- samples[, -hh[1]]

  return(samples)
}


#' Read scRNA-seq counts matrix in h5 format.
#' Automatically 'infer' counts;features;cells from the h5 file structure.
#' Attempts multiple ways.
#' @export
read_h5_counts <- function(h5.file) {

  message("[playbase::read_h5_counts] Reading h5 file: ", h5.file)
  df <- NULL
  
  FF <- tryCatch( { rhdf5::h5ls(h5.file) }, error = function(w) { NULL } )

  if (!is.null(FF) && all(c("group","name") %in% colnames(FF))) {

    h5.ems <- paste0(FF[,"group"], "/", FF[,"name"])
    LL <- lapply(h5.ems, function(ems) rhdf5::h5read(file = h5.file, name = ems))

    dims <- unlist(lapply(LL, function(x) length(dim(x))))
    ll <- lapply(LL, length)
    
    if (any(dims == 2)) df <-  LL[[which(dims == 2)[1]]]

    if (!is.null(df)) {
      if (is.null(rownames(df)) && any(ll == nrow(df))) {
        rownames(df) <- as.character(LL[[which(ll == nrow(df))[1]]])
      }
      if (is.null(colnames(df)) && any(ll == ncol(df))) {
        colnames(df) <- as.character(LL[[which(ll == ncol(df))[1]]])
      }
    }
  }

  if (is.null(df))
    df <- tryCatch({ Seurat::Read10X_h5(h5.file) }, error = function(w) { NULL } )

  if (is.null(df))
    df <- tryCatch( { playbase::h5.readMatrix(h5.file) }, error = function(w) { NULL } )

  if (!is.null(df) & all(class(df) %in% c("matrix", "array"))) {
    message("[playbase::read_h5_counts] success!")
  } else {
    message("[playbase::read_h5_counts] df is null!")
  }

  return(df)

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
    message("[read_annot] ERROR: input not valid.")
    return(NULL)
  }

  ## add column title
  if (colnames(df)[1] == "") colnames(df)[1] <- "row.names"

  ## determine last character column
  df1 <- type.convert(data.frame(head(df, 20), check.names = FALSE), as.is = TRUE)
  col.type <- sapply(df1, class)
  xannot.names <- "gene|symbol|protein|compound|title|description|name|position"
  is.xannot <- grepl(xannot.names, tolower(colnames(df)))
  char.cols <- which(col.type == "character" | is.xannot)
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
  } else {
    df <- NULL
  }
  return(df)
}

getError <- function(e, what = "Description") {
  ERROR_MSG <- playbase::PGX_CHECKS
  if (!e %in% ERROR_MSG$error) {
    return(paste("unknown error", e))
  }
  ERROR_MSG[match(e, ERROR_MSG$error), what]
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
validate_counts <- function(df) {
  if (is.character(df) && is.null(dim(df))) {
    df <- read.as_matrix(df)
  }
  chk <- pgx.checkINPUT(df, "COUNTS")
  ERROR_MSG <- playbase::PGX_CHECKS
  err <- names(chk$checks)
  if (length(err)) {
    msg <- lapply(err, function(e) getError(e))
    msg <- paste(msg, collapse = "; ")
    message("WARNING: ", msg)
  }
  chk$PASS
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
validate_samples <- function(df) {
  if (is.character(df) && is.null(dim(df))) {
    df <- read.as_matrix(df)
  }
  chk <- pgx.checkINPUT(df, "SAMPLES")
  ERROR_MSG <- playbase::PGX_CHECKS
  err <- names(chk$checks)
  if (length(err)) {
    msg <- lapply(err, function(e) getError(e))
    msg <- paste(msg, collapse = "; ")
    message("WARNING: ", msg)
  }
  chk$PASS
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
validate_contrasts <- function(df) {
  if (is.character(df) && is.null(dim(df))) {
    df <- read.as_matrix(df)
  }
  chk <- pgx.checkINPUT(df, "CONTRASTS")
  ERROR_MSG <- playbase::PGX_CHECKS
  err <- names(chk$checks)
  if (length(err)) {
    msg <- lapply(err, function(e) getError(e))
    msg <- paste(msg, collapse = "; ")
    message("WARNING: ", msg)
  }
  chk$PASS
}

## --------------------------------------------------------------------
## --------------------------------------------------------------------
## --------------------------------------------------------------------
