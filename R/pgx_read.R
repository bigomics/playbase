
#' Read counts data from file
#'
#' @param file string. path to file
#'
#' @return matrix. the file with the data
#' @export
#'
#' @examples
#' counts <- read_counts(playbase::playbase_example('counts.csv'))
read_counts <- function(file, convert_names = TRUE) {
    df <- .read_as_matrix(file)

    # convert to gene names (needed for biological effects)
    if (convert_names) {
        pp <- rownames(df)
        rownames(df) <- .probe2symbol(pp)
        sel <- !(rownames(df) %in% c(NA, '', 'NA'))
        df <- df[sel, ]
        xx <- tapply(1:nrow(df), rownames(df),
                     function(i) colSums(df[i, , drop = FALSE]))
        df <- do.call(rbind, xx)
    }

    df
}

#' Read expression data from file
#'
#' @param file string. path to file
#'
#' @return matrix. the file with the data
#' @export
#'
#' @examples
#' counts <- read_expression(playbase::playbase_example('counts.csv'))
read_expression <- function(file, convert_names = TRUE) {
    df <- .read_as_matrix(file)
    df <- df ** 2

    # convert to gene names (needed for biological effects)
    if (convert_names) {
        pp <- rownames(df)
        rownames(df) <- .probe2symbol(pp)
        sel <- !(rownames(df) %in% c(NA, '', 'NA'))
        df <- df[sel, ]
        xx <- tapply(1:nrow(df), rownames(df),
                     function(i) colSums(df[i, , drop = FALSE]))
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
#' samples <- read_samples(playbase::playbase_example('samples.csv'))
read_samples <- function(file) {
    df <- .read_as_matrix(file)
    df <- as.data.frame(df)
    df
}

#' Read contrasts data from file
#'
#' @param file string. path to file
#'
#' @return matrix. the file with the data
#' @export
#'
#' @examples
#' contrasts <- read_contrasts(playbase::playbase_example('contrasts.csv'))
read_contrasts <- function(file) {
    df <- .read_as_matrix(file)
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
.read_as_matrix <- function(file) {
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