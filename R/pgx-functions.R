##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##



#' @title Create a Phenotype Matrix
#'
#' @description This function creates a phenotype matrix from a PGX object.
#'
#' @param pgx A list representing a PGX object containing the data to be analyzed.
#' @param phenotype A character string specifying the name of the phenotype to use for creating the matrix.
#'
#' @details This function takes a PGX object `pgx` and a character string `phenotype` as input and
#' creates a phenotype matrix using the specified phenotype. The phenotype matrix is created by extracting
#' the specified phenotype from the `samples` data frame in the `pgx` object and converting it into a model matrix.
#' The resulting phenotype matrix is returned as a numeric matrix.
#'
#' @return A numeric matrix representing the phenotype matrix created from the specified phenotype.
#'
#' @export
pgx.phenoMatrix <- function(pgx, phenotype) {
  y <- pgx$samples[, phenotype]
  mm <- t(stats::model.matrix(~ 0 + y))
  rownames(mm) <- sub("^y", "", rownames(mm))
  colnames(mm) <- rownames(pgx$samples)
  as.matrix(mm)
}


#' @title Compact Positions
#'
#' @description This function compacts a matrix of positions by removing white space.
#'
#' @param pos A numeric matrix of positions to be compacted.
#' @param d An optional numeric value specifying the maximum distance between adjacent positions after compaction.
#' The default value is 0.01.
#'
#' @details This function takes a numeric matrix of positions `pos` and an optional numeric value `d` as input and compacts the positions by removing white space.
#' For each column in `pos`, the function calculates the distances between adjacent positions and replaces any distance greater than `d` with `d`.
#' The resulting compacted positions are returned as a numeric matrix of the same size as `pos`.
#'
#' @return A numeric matrix of the same size as `pos`, containing the compacted positions.
#'
#' @export
pos.compact <- function(pos, d = 0.01) {
  ## make positions more dense removing white space
  for (i in 1:ncol(pos)) {
    x <- pos[, i]
    dr <- d * diff(range(x))
    names(x) <- 1:nrow(pos)
    ii <- order(x)
    x1 <- cumsum(c(x[ii[1]], pmin(diff(x[ii]), dr)))
    pos[, i] <- x1[order(as.integer(names(x1)))]
  }
  pos
}

#' Find non-overlapping boxes for points
#'
#' @param df Data frame containing point positions
#' @param xcol Name of column in \code{df} containing x coordinates
#' @param ycol Name of column in \code{df} containing y coordinates
#' @param box_padding_x Padding to add to box width
#' @param box_padding_y Padding to add to box height
#' @param point_padding_x Minimum spacing between points in x direction
#' @param point_padding_y Minimum spacing between points in y direction
#' @param xlim Limits for x axis
#' @param ylim Limits for y axis
#' @param force Repulsion force between boxes
#' @param maxiter Maximum number of iterations
#'
#' @return Data frame with updated non-overlapping box coordinates
#'
#' @details This function takes a data frame of point positions and calculates non-overlapping
#' box positions around each point using a physical simulation.
#'
#' It returns a data frame containing the updated box corner positions after repulsion.
#' The x1,y1,x2,y2 columns correspond to the box corners.
#'
#' @examples
#' \dontrun{
#' points <- data.frame(x = runif(10), y = runif(10))
#' boxes <- util.findboxes(points, "x", "y", 0.2, 0.2, 0.1, 0.1)
#' }

#' @export
util.findboxes <- function(df, xcol, ycol,
                           box_padding_x, box_padding_y,
                           point_padding_x, point_padding_y,
                           xlim, ylim,
                           force = 1e-7, maxiter = 20000) {
  # https://github.com/slowkow/ggrepel/issues/24

  # x and y posiitons as a dataframe
  posdf <- df[c(xcol, ycol)]

  # returnd a df where columns are points
  boxdf <- apply(posdf, 1, function(row) {
    xval <- row[xcol]
    yval <- row[ycol]
    return(c(
      xval - box_padding_x / 2,
      yval - box_padding_y / 2,
      xval + box_padding_x / 2,
      yval + box_padding_y / 2
    ))
  })
  # columns are x1,y1,x2,y2
  boxmatrix <- as.matrix(t(boxdf))

  moved <- ggrepel:::repel_boxes(
    data_points = as.matrix(posdf),
    point_padding_x = point_padding_x,
    point_padding_y = point_padding_y,
    boxes = boxmatrix,
    xlim = xlim,
    ylim = ylim,
    hjust = 0.5,
    vjust = 0.5,
    force = force,
    maxiter = maxiter
  )

  finaldf <- cbind(posdf, moved)
  names(finaldf) <- c("x1", "y1", "x2", "y2")
  return(finaldf)
}


#' @title Generate Star Symbols
#'
#' @description This function generates a string of star symbols.
#'
#' @param n A numeric value specifying the number of star symbols to generate.
#' @param pch An optional character string specifying the symbol to use for the stars.
#' The default value is a black star.
#'
#' @details This function takes a numeric value `n` and an optional character string `pch` as
#' input and generates a string of `n` star symbols using the specified `pch` symbol.
#' If `n` is 0, an empty string is returned.
#' The resulting string of star symbols is returned as a character vector of length 1.
#'
#' @return A character vector of length 1, containing the generated string of star symbols.
#'
#' @examples
#' \dontrun{
#' # example code
#' star.symbols(3)
#' star.symbols(5, pch = "*")
#' }
#' @export
star.symbols <- function(n, pch = "\u2605") {
  if (n == 0) {
    return("")
  }
  paste(rep(pch, n), collapse = "")
}


#' Search for file in system path
#'
#' @param path Character vector of directories to search
#' @param file Character string of filename to search for
#'
#' @return Full path of file if found, otherwise NULL
#'
#' @details Searches the directories in \code{path} to find the file specified by \code{file}.
#' Returns the full path if the file is found in one of the directories.
#' Returns NULL if the file is not found.
#'
#' @export
search_path <- function(paths, file) {
  dir <- paths[which(file.exists(file.path(paths, file)))]
  if (length(dir) == 0) {
    return(NULL)
  }
  file.path(dir[1], file)
}


#' @title Row Scale a Matrix
#'
#' @description This function scales the rows of a matrix by subtracting the row
#' means and dividing by the row standard deviations.
#'
#' @param x A numeric matrix to be scaled.
#'
#' @details This function takes a numeric matrix `x` as input and scales its rows by
#' subtracting the row means and dividing by the row standard deviations.
#' The resulting scaled matrix is returned.
#'
#' @return A numeric matrix of the same size as `x`, containing the scaled values.
#'
#' @examples
#' \dontrun{
#' # example code
#' x <- matrix(c(1, 2, 3, 4), nrow = 2)
#' rowscale(x)
#' }
#'
#' @export
rowscale <- function(x) {
  x <- x - Matrix::rowMeans(x, na.rm = TRUE)
  x / (1e-4 + sqrt(rowMeans(x**2, na.rm = TRUE)))
}


#' @title Wrap Strings to a Specified Width
#'
#' @description This function wraps a character vector of strings to a specified width.
#'
#' @param str A character vector of strings to be wrapped.
#' @param n A numeric value specifying the maximum width of each line in characters.
#'
#' @details This function takes a character vector of strings `str` and a numeric value `n` as input and wraps each string in `str` to a maximum width of `n` characters.
#' The resulting wrapped strings are returned as a character vector, where each element contains one or more lines separated by newline characters.
#'
#' @return A character vector of the same length as `str`, containing the wrapped strings.
#'
#' @examples
#' \dontrun{
#' # example code
#' str <- c("The quick brown fox jumps over the lazy dog")
#' strwrap2(str, 10)
#' }
#'
#' @export
strwrap2 <- function(str, n) {
  sapply(str, function(s) paste(base::strwrap(s, n), collapse = "\n"))
}


#' @title Add Opacity to Hexadecimal Colors
#'
#' @description This function adds opacity to a vector of hexadecimal color values.
#'
#' @param hexcol A character vector of hexadecimal color values.
#' @param opacity A numeric value specifying the opacity to be added, ranging from 0 (transparent) to 1 (opaque).
#'
#' @details This function takes a character vector of hexadecimal color values `hexcol` and
#' a numeric value `opacity` as input and adds the specified opacity to each color value.
#' The function uses the `toRGB` function from the `plotly` package to convert the hexadecimal
#' color values to RGBA format, then modifies the alpha channel according to the specified `opacity` value.
#' The resulting RGBA color values are returned as a character vector.
#'
#' @return A character vector of the same length as `hexcol`, containing the RGBA color values with added opacity.
#'
#' @examples
#' \dontrun{
#' # example code
#' hexcol <- c("#FF0000", "#00FF00", "#0000FF")
#' add_opacity(hexcol, 0.5)
#' }
#' @export
add_opacity <- function(hexcol, opacity) {
  col1 <- rep(NA, length(hexcol))
  ii <- which(!is.na(hexcol))
  rgba <- strsplit(gsub("rgba\\(|\\)", "", plotly::toRGB(hexcol[ii], opacity)), split = ",")
  rgba <- apply(do.call(rbind, rgba), 2, as.numeric)
  if (length(hexcol) == 1) rgba <- matrix(rgba, nrow = 1)
  col1[ii] <- grDevices::rgb(rgba[, 1] / 255, rgba[, 2] / 255, rgba[, 3] / 255, rgba[, 4])
  col1
}




#' @title Check for Required Fields in a PGX Object
#'
#' @description This function checks if a PGX object contains all required fields.
#'
#' @param pgx A list representing a PGX object to be checked.
#'
#' @details This function takes a list `pgx` representing a PGX object as input and checks if it contains all required fields.
#' The required fields are "counts", "samples", "genes", "model.parameters", "X", "gx.meta", "gset.meta", "gsetX", and "GMT".
#' If any of the required fields are not present in the `pgx` object, a warning message is printed.
#' The function returns a logical value indicating whether all required fields are present in the `pgx` object.
#'
#' @return A logical value indicating whether all required fields are present in the `pgx` object.
#'
#' @export
pgx.checkObject <- function(pgx) {
  must.have <- c(
    "counts", "samples", "genes", "model.parameters",
    "X", "gx.meta", "GMT"
  )
  not.present <- setdiff(must.have, names(pgx))
  if (length(not.present) > 0) {
    not.present <- paste(not.present, collapse = " ")
    message("[pgx.checkObject] WARNING!!! object does not have: ", not.present)
  }
  all(must.have %in% names(pgx))
}


#' Calculate group-wise summary statistics from a matrix
#'
#' @param X Numeric matrix with samples in columns.
#' @param group Factor indicating group membership for each column of X.
#' @param FUN Function to calculate row-wise statistics. Default is rowMeans.
#' @param dir Direction to apply FUN. 1=row-wise, 2=column-wise.
#'
#' @return Matrix of summary statistics, with rows corresponding to groups.
#'
#' @details This function calculates a row-wise or column-wise summary statistic
#'   for each group defined in the \code{group} vector.
#'
#'   It applies the function specified by \code{FUN} (default is rowMeans) to summarize
#'   either rows or columns of \code{X}, depending on the \code{dir} argument.
#'
#'   The result is a matrix with rows corresponding to the unique groups, and columns
#'   equal to the number of rows/columns summarized.
#'
#' @examples
#' \dontrun{
#' X <- matrix(rnorm(100 * 20), 100, 20)
#' groups <- sample(letters[1:5], 20, replace = TRUE)
#' means <- matGroupMeans(X, groups)
#' }
#' @export
matGroupMeans <- function(X, group, FUN = rowMeans, dir = 1, reorder = TRUE) {
  if (dir == 2) X <- t(X)
  mX <- do.call(cbind, tapply(
    1:ncol(X), group,
    function(i) FUN(X[, i, drop = FALSE], na.rm = TRUE)
  ))
  if (!reorder) mX <- mX[, unique(group), drop = FALSE]
  if (dir == 2) mX <- t(mX)
  mX
}

#' Calculate group-wise row means (like base::rowsum) from a
#' matrix. Faster than matGroupMeans.
#'
#' @export
rowmean <- function(X, group, reorder = TRUE) {
  sumX <- base::rowsum(X, group, na.rm = TRUE, reorder = reorder)
  nX <- base::rowsum(1 * (!is.na(X)), group, reorder = reorder)
  sumX / nX
}

#' @describeIn trimsame0 trimsame is a function that trims common prefixes and/or
#' suffixes from a character vector by applying trimsame0 forwards and/or backwards.
#' @export
trimsame <- function(s, split = " ", ends = TRUE, summarize = FALSE) {
  if (all(is.na(s)) || all(s == "")) {
    return(s)
  }
  if (ends) {
    return(trimsame.ends(s, split = split, summarize = summarize))
  }
  return(trimsame0(s, split = split, summarize = summarize))
}


#' @describeIn trimsame0 trimsame.ends is a function that trims common prefixes and
#' suffixes from a character vector by applying trimsame0 forwards and backwards.
#' @export
trimsame.ends <- function(s, split = " ", summarize = FALSE) {
  s1 <- trimsame0(s, split = split, summarize = summarize)
  s2 <- sapply(strsplit(s1, split = split), function(x) paste(rev(x), collapse = " "))
  s3 <- trimsame0(s2, split = split, summarize = summarize, rev = TRUE)
  s4 <- sapply(strsplit(s3, split = split), function(x) paste(rev(x), collapse = " "))
  s4
}


#' @title Trim Common Prefix from Strings
#'
#' @description This function trims the common prefix from a character vector of strings.
#'
#' @param s A character vector of strings to be trimmed.
#' @param split An optional character string specifying the delimiter used to split the strings into words.
#' The default value is a space (" ").
#' @param summarize An optional logical value indicating whether to summarize the common prefix using the first letter of each word.
#' The default value is `FALSE`.
#' @param rev An optional logical value indicating whether to reverse the order of the words in the summarized prefix.
#' The default value is `FALSE`.
#'
#' @details The function takes a character vector of strings `s` as input and searches for a common prefix among the strings.
#' The common prefix is defined as the longest sequence of characters that is shared by all strings and ends with the specified delimiter.
#' If a common prefix is found, it is trimmed from all strings.
#' If the `summarize` parameter is `TRUE`, the common prefix is summarized using the first letter of each word and appended to the beginning of each string.
#' If the `rev` parameter is `TRUE`, the order of the words in the summarized prefix is reversed.
#'
#' @return A character vector of the same length as `s`, containing the trimmed strings.
#'
#' @export
trimsame0 <- function(s, split = " ", summarize = FALSE, rev = FALSE) {
  if (all(is.na(s)) || all(s == "")) {
    return(s)
  }
  s <- strsplit(s, split)
  s.orig <- s

  ##
  i <- 1
  is.same <- FALSE
  while (i < 1000 && !is.same) {
    s1 <- sapply(s, "[", 1)
    slen <- sapply(s, length)
    ss <- setdiff(s1, NA)
    if (all(ss == ss[1])) {
      s <- lapply(s, "[", -1)
    } else {
      is.same <- FALSE
    }
    i <- i + 1
  }
  sapply(s, length)

  if (sapply(s, length) == 0 || all(s == "")) {
    sx <- sapply(s.orig, "[", 1)
    sx
    return(sx)
  }

  i <- 1
  is.same <- FALSE
  while (i < 1000 && !is.same) {
    slen <- sapply(s, length)
    s2 <- sapply(s, tail, 1)
    ss <- setdiff(s2, NA)
    if (all(ss == ss[1])) {
      s <- mapply(head, s, slen - 1)
    } else {
      is.same <- FALSE
    }
    i <- i + 1
  }

  if (sapply(s, length) == 0 || all(s == "")) {
    s <- sapply(s.orig, "[", 1)
  }
  s <- sapply(s, paste, collapse = split)
  s
}


#' @title Tag Duplicate Strings
#'
#' @description This function tags duplicate strings by appending spaces to them.
#'
#' @param s A character vector of strings to be tagged.
#'
#' @details The function takes a character vector of strings `s` as
#' input and searches for duplicate strings. For each set of duplicate strings,
#' the function appends a number of spaces to each string, such that the first
#' occurrence of the string has no spaces appended, the second occurrence has one space appended, and so on.
#' The resulting tagged strings are returned as a character vector.
#'
#' @return A character vector of the same length as `s`, containing the tagged strings.
#'
#' @export
tagDuplicates <- function(s) {
  ## Tag duplicate with blanks
  ##
  jj <- which(duplicated(s))
  t <- s[jj][1]
  for (t in unique(s[jj])) {
    ii <- which(s == t)
    spaces <- paste(rep(" ", length(ii)), collapse = "")
    blanks <- substring(spaces, 0, 0:(length(ii) - 1))

    s[ii] <- paste0(s[ii], blanks)
  }
  s <- gsub("[.]1$", "", s)
  s
}


#' @title Wrap gene symbols in hyperlinks
#'
#' @description This function wraps gene symbols in a character vector
#' with HTML hyperlinks to external databases.
#'
#' @param s Character vector of gene symbols.
#' @param gs Character vector of gene symbols with database prefixes.
#'
#' @details This function takes a character vector of gene symbols \code{s} and a corresponding
#' character vector \code{gs} containing gene symbols prefixed with database identifiers (e.g. ENSEMBL:).
#' It searches for matches in \code{gs} and replaces the plain symbols in \code{s} with HTML hyperlinks.
#'
#' The following database prefixes are recognized:
#' - ENSEMBL: - Links to Ensembl database
#' - UNIPROT: - Links to Uniprot database
#' - LNCRNADB: - Links to LncRNADB database
#' - BIOPLEX: - Links to BIOPLEX database
#'
#' @return Character vector \code{s} with gene symbols replaced by hyperlinks.
#'
#' @examples
#' \dontrun{
#' s <- c("TP53", "MYC")
#' gs <- c("ENSEMBL:ENSG00000171791", "UNIPROT:P01106")
#' s1 <- wrapHyperLink(s, gs)
#' print(s1)
#' }
#'
#' @export
wrapHyperLink <- function(s, gs) {
  gs <- as.character(gs)
  s1 <- s <- as.character(s)
  ## GEO/GSE accession
  jj <- grep("GSE[0-9]", gs, ignore.case = TRUE)
  if (length(jj)) {
    acc <- sub("[-_ ].*", "", gsub("^.*GSE", "GSE", gs[jj], ignore.case = TRUE))
    url <- paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", acc)
    s1[jj] <- paste0("<a href='", url, "' target='_blank'>", s[jj], "</a>")
  }

  ## GDS accession
  jj <- grep("GDS[0-9]", gs, ignore.case = TRUE)
  if (length(jj)) {
    acc <- sub("[-_ ].*", "", gsub("^.*GDS", "GDS", gs[jj], ignore.case = TRUE))
    url <- paste0("https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=", acc)
    s1[jj] <- paste0("<a href='", url, "' target='_blank'>", s[jj], "</a>")
  }

  ## Reactome accession
  jj <- grep("R-HSA-[0-9][0-9]", gs)
  if (length(jj)) {
    id <- sub("^.*R-HSA-", "R-HSA-", gs[jj])
    url <- paste0("https://reactome.org/content/detail/", id)
    s1[jj] <- paste0("<a href='", url, "' target='_blank'>", s[jj], "</a>")
  }

  ## Wikipathways accession
  jj <- grep("_WP[0-9][0-9]", gs)
  if (length(jj)) {
    id <- sub("^.*_WP", "WP", gs[jj])
    id <- stringr::str_extract(id, "WP\\d+")
    url <- paste0("https://www.wikipathways.org/index.php/Pathway:", id)
    s1[jj] <- paste0("<a href='", url, "' target='_blank'>", s[jj], "</a>")
  }

  ## MSigDB accession
  jj <- grep("^H:|^C[1-8]:|HALLMARK", gs)
  if (length(jj)) {
    gs1 <- sub("^.*:", "", gs[jj])
    url <- paste0("http://software.broadinstitute.org/gsea/msigdb/cards/", gs1)
    s1[jj] <- paste0("<a href='", url, "' target='_blank'>", s[jj], "</a>")
  }

  ## GO accession
  jj <- grep("\\(GO_.*\\)$", gs)
  if (length(jj)) {
    id <- gsub("^.*\\(GO_|\\)$", "", gs[jj])
    url <- paste0("http://amigo.geneontology.org/amigo/term/GO:", id)
    s1[jj] <- paste0("<a href='", url, "' target='_blank'>", s[jj], "</a>")
  }

  ## ELSEVIRE accession
  jj <- grep("ELSEVIER:", gs)
  if (length(jj)) {
    id <- sub("^.*:", "", gs[jj])
    base_url <- paste0("https://mammalcedfx.pathwaystudio.com/app/sd")
    url <- paste0(base_url, "?full=true&layout=flat&entitylist=", URLencode(id), "&separator=", ",")
    s1[jj] <- paste0("<a href='", url, "' target='_blank'>", s[jj], "</a>")
  }

  # BIOPLANET accession
  jj <- grep("BIOPLANET:", gs)
  if (length(jj)) {
    # bioplanet is down, placeholder for when it comes back alive
    url <- "https://tripod.nih.gov/bioplanet/"
    s1[jj] <- paste0("<a href='", url, "' target='_blank'>", s[jj], "</a>")
  }

  # LNCHUB accession
  jj <- grep("LNCHUB:", gs)
  if (length(jj)) {
    # bioplanet is down, placeholder for when it comes back alive
    id <- sub("^.*:", "", gs[jj])
    base_url <- paste0("https://maayanlab.cloud/lnchub/?lnc=")
    url <- paste0(base_url, URLencode(id))
    s1[jj] <- paste0("<a href='", url, "' target='_blank'>", s[jj], "</a>")
  }
  # BIOPLEX accession
  jj <- grep("BIOPLEX:", gs)
  if (length(jj)) {
    # bioplanet is down, placeholder for when it comes back alive
    base_url <- "https://bioplex.hms.harvard.edu/explorer/externalQuery.php?geneQuery="
    url <- sapply(gs[jj], function(x) {
      split_x <- strsplit(x, "_")[[1]]
      if (length(split_x) > 1) {
        url <- paste0(base_url, URLencode(split_x[[2]]))
        return(url)
      } else {
        return(base_url)
      }
    })

    s1[jj] <- paste0("<a href='", url, "' target='_blank'>", s[jj], "</a>")
  }
  return(s1)
}


#' @title Reverse Comparison Strings
#'
#' @description This function reverses the order of comparison strings in the format "A_vs_B".
#'
#' @param comp A character vector of comparison strings to be reversed.
#'
#' @details The function takes a character vector of comparison strings `comp` as
#' input and reverses the order of the strings in the format "A_vs_B" or "A_VS_B".
#' The function also preserves any prefix before a colon (":") or postfix after an at
#' symbol ("@") in the input strings.
#'
#' @return A character vector of the same length as `comp`, containing the reversed comparison strings.
#'
#' @examples
#' \dontrun{
#' # example code
#' comp <- c("A_vs_B", "C_VS_D:prefix@postfix")
#' reverse.AvsB(comp)
#' }
#' @export
reverse.AvsB <- function(comp) {
  reverse.AvsB.1 <- function(comp) {
    prefix <- postfix <- ""
    if (any(grepl("[:]", comp))) prefix <- sub(":.*", "", comp)
    if (any(grepl("[@]", comp))) postfix <- sub(".*@", "", comp)
    comp0 <- gsub(".:|@.*", "", comp)
    ab <- paste(rev(strsplit(comp0, split = "_vs_|_VS_")[[1]]), collapse = "_vs_")
    gsub("^:|@$", "", paste0(prefix, ":", ab, "@", postfix))
  }
  as.character(sapply(comp, reverse.AvsB.1))
}


#' Check if two groups represent a positive vs negative comparison
#'
#' @param pgx pgx object containing the data
#'
#' @return Logical indicating if the two groups represent a positive vs negative comparison
#'
#' @details This function examines two groups from the pgx object to determine if they represent a positive vs negative comparison.
#' It checks the mean expression of signature genes in each group to see if one group has overall higher expression.
#' It also looks for keywords like "neg", "untr", "wt" etc in the group names.
#' Based on these checks it returns a logical indicating if the groups represent a positive vs negative comparison.
#'
#' @export
is.POSvsNEG <- function(pgx) {
  ## Determines automagically from contrast matrix if notation is
  ## 'A_vs_B' or 'B_vs_A' (which group is positive in the contrast
  ## matrix). Too complicated... maybe we should just require one
  ## definition...
  ##
  ## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ## !!!!!!!!!!!! We should get rid of this... !!!!!!!!!!!!!
  ## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  cntrmat <- pgx$model.parameters$contr.matrix
  design <- pgx$model.parameters$design


  ## rely only on contrasts with '_vs_'
  cntrmat <- cntrmat[, grep("_vs_", colnames(cntrmat)), drop = FALSE]
  grp1 <- sapply(strsplit(colnames(cntrmat), split = "_vs_"), "[", 1)
  grp2 <- sapply(strsplit(colnames(cntrmat), split = "_vs_"), "[", 2)
  grp1 <- sub(".*[:]|@.*", "", grp1)
  grp2 <- sub(".*[:]|@.*", "", grp2)
  is.null(design)
  is.PosvsNeg1 <- NA
  grp1
  grp2
  if (FALSE && !is.null(design)) {
    is.pn <- rep(NA, length(grp1))
    i <- 1
    for (i in 1:length(grp1)) {
      j1 <- grep(grp1[i], rownames(cntrmat), fixed = TRUE)
      j2 <- grep(grp2[i], rownames(cntrmat), fixed = TRUE)
      grp1.sign <- mean(cntrmat[j1, i], na.rm = TRUE)
      grp2.sign <- mean(cntrmat[j2, i], na.rm = TRUE)
      grp1.sign
      grp2.sign
      if (!is.nan(grp1.sign) && !is.nan(grp2.sign)) {
        is.pn[i] <- (grp1.sign > grp2.sign)
      }
    }
    is.pn
    is.PosvsNeg1 <- mean(is.pn, na.rm = TRUE) > 0
  } else {
    ## This uses the experiment matrix (sample-based contrast) and
    ## the annotation to determine is A_vs_B or B_vs_A was
    ## intended.
    ##
    expmat <- pgx$model.parameters$exp.matrix
    is.pn <- rep(NA, length(grp1))
    i <- 1
    for (i in 1:length(grp1)) {
      a1 <- apply(pgx$samples, 1, function(a) mean(grepl(grp1[i], a)))
      a2 <- apply(pgx$samples, 1, function(a) mean(grepl(grp2[i], a)))
      j1 <- which(a1 > a2) ## samples with phenotype  more in grp1
      j2 <- which(a2 >= a1) ## samples with phenotype  more in grp2
      s1 <- s2 <- 0
      if (length(j1)) s1 <- rowMeans(expmat[j1, i, drop = FALSE] > 0, na.rm = TRUE)
      if (length(j2)) s2 <- rowMeans(expmat[j2, i, drop = FALSE] > 0, na.rm = TRUE)
      if (mean(s1) > mean(s2)) is.pn[i] <- TRUE
      if (mean(s2) > mean(s1)) is.pn[i] <- FALSE
    }

    is.PosvsNeg1 <- mean(is.pn, na.rm = TRUE) > 0
  }

  ## look for keywords
  grp1.neg2 <- mean(grepl("neg|untr|ref|wt|ctr|control", tolower(grp1)))
  grp2.neg2 <- mean(grepl("neg|untr|ref|wt|ctr|control", tolower(grp2)))

  is.PosvsNeg2 <- NA
  if (grp1.neg2 > 0 || grp2.neg2 > 0) {
    is.PosvsNeg2 <- (grp2.neg2 > grp1.neg2)
  }

  ok <- setdiff(c(is.PosvsNeg1, is.PosvsNeg2), NA)[1] ## priority to first test??
  if (is.na(ok) || length(ok) == 0) ok <- TRUE ## DEFAULT if not known !!!!!
  ok
}


#' @title Check if a Variable is Categorical
#'
#' @description This function checks if a variable is categorical based on its class,
#' number of unique values, and proportion of non-missing values.
#'
#' @param x A vector representing the variable to be checked.
#' @param max.ncat An optional numeric value specifying the maximum number of unique
#' values for a variable to be considered categorical.
#' If not provided, the default value is one less than the length of `x`.
#' @param min.ncat An optional numeric value specifying the minimum number of unique
#' values for a variable to be considered categorical. The default value is 2.
#'
#' @details The function takes a vector `x` as input and checks if it is a factor or character vector.
#' If `x` is not a factor or character vector, the function returns `FALSE`.
#' Otherwise, the function calculates the number of unique values in `x` (excluding missing values) and the proportion of non-missing values in `x`.
#' If the number of unique values is greater than 80% of the proportion of non-missing values, the function returns `FALSE`.
#' Otherwise, if the number of unique values is within the range specified by the `min.ncat` and `max.ncat` parameters, the function returns `TRUE`.
#' Otherwise, the function returns `FALSE`.
#'
#' @return A logical value indicating whether the input variable is categorical.
#'
#' @export
is.categorical <- function(x, max.ncat = NULL, min.ncat = 2) {
  max.ncat <- length(x) - 1
  is.factor <- any(class(x) %in% c("factor", "character"))
  n.unique <- length(unique(setdiff(x, NA)))
  n.notna <- length(x[!is.na(x)])
  is.id <- (n.unique > 0.8 * n.notna)
  is.factor2 <- (is.factor & !is.id & n.unique >= min.ncat & n.unique <= max.ncat)
  return(is.factor2)
}


#' Discretize phenotype matrix
#'
#' @param df Data frame containing phenotype data
#' @param min.ncat Minimum number of categories for discretization. Default is 2.
#' @param max.ncat Maximum number of categories for discretization. Default is 20.
#' @param remove.dup Logical indicating whether to remove duplicate factor levels. Default is FALSE.
#'
#' @return Discretized phenotype data frame
#'
#' @details This function discretizes the phenotype data frame \code{df} into categorical variables.
#' Numerical variables are converted into "high" and "low" values based on the median.
#' Categorical variables are kept as factors.
#'
#' The number of categories is controlled by \code{min.ncat} and \code{max.ncat} parameters.
#' Duplicate factor levels can be removed by setting \code{remove.dup = TRUE}.
#'
#' @export
pgx.discretizePhenotypeMatrix <- function(df, min.ncat = 2, max.ncat = 20, remove.dup = FALSE) {
  catpheno <- pgx.getCategoricalPhenotypes(
    df,
    max.ncat = max.ncat, min.ncat = min.ncat, remove.dup = remove.dup
  )
  numpheno <- pgx.getNumericalPhenotypes(df)

  numpheno <- setdiff(numpheno, catpheno) ## already in categories?
  df.num <- c()
  if (length(numpheno)) {
    df.num <- as.matrix(df[, numpheno, drop = FALSE])
    is.high <- t(t(df.num) > apply(df.num, 2, stats::median, na.rm = TRUE))
    df.num[is.high] <- "high"
    df.num[!is.high] <- "low"
  }
  df1 <- df[, 0]
  if (length(catpheno)) df1 <- cbind(df1, df[, catpheno, drop = FALSE])
  if (length(numpheno)) df1 <- cbind(df1, df.num)
  rownames(df1) <- rownames(df)
  df1
}


#' @title Get Numerical Phenotypes
#'
#' @description This function retrieves the column names of numerical
#' phenotypes from a data frame.
#'
#' @param df A data frame containing the data to be analyzed.
#'
#' @details The function takes a data frame `df` as input and searches for columns that
#' represent numerical phenotypes.
#' Numerical phenotypes are identified as columns that are numeric, have a majority of
#' unique values, and do not contain certain keywords in their names or values (e.g., "sample",
#' "id", "replicate", "patient", "donor", "individual", "year", "month", "day", "survival").
#' The names of the identified columns are returned as a character vector.
#'
#' @return A character vector representing the column names of the numerical phenotypes in the input data frame.
#'
#' @export
pgx.getNumericalPhenotypes <- function(df) {
  is.bad <- 0
  is.bad1 <- grepl("^sample$|[_.]id$|replic|rep|patient|donor|individ", tolower(colnames(df)))
  is.bad2 <- grepl("year|month|day|^efs|^dfs|surv|follow", tolower(colnames(df)))
  is.bad3 <- apply(df, 2, function(x) any(grepl("^sample|patient|replicate|donor|individ", x, ignore.case = TRUE)))
  is.bad <- (is.bad1 | is.bad2 | is.bad3)


  numratio <- apply(df, 2, function(x) length(unique(x))) / nrow(df)

  numpheno <- (apply(df, 2, is.num) & !is.bad & numratio > 0.5)

  names(which(numpheno == TRUE))
}


#' @title Get Categorical Phenotypes
#'
#' @description
#' Identifies categorical phenotype columns in a data frame.
#'
#' @param df Data frame containing phenotype data.
#' @param min.ncat Minimum number of categories required. Default 2.
#' @param max.ncat Maximum number of categories allowed. Default 20.
#' @param remove.dup Logical indicating whether to remove duplicate factor
#' levels. Default FALSE.
#'
#' @details
#' This function examines a phenotype data frame \code{df} and identifies
#' columns that represent categorical variables.
#' It first excludes any columns that appear to contain sample IDs.
#' It then selects columns that are factors, or character/numeric with
#' between \code{min.ncat} and \code{max.ncat} unique values.
#'
#' The function returns a character vector of column names representing putative categorical phenotypes.
#' Setting \code{remove.dup = TRUE} will remove any duplicate factor levels prior to counting number of categories.
#'
#' @return
#' A character vector of categorical phenotype column names
#'
#' @export
pgx.getCategoricalPhenotypes <- function(df, min.ncat = 2, max.ncat = 20, remove.dup = FALSE) {
  is.bad <- 0

  ## ... exclude sample IDs

  is.bad1 <- grepl("^sample$|[_.]id$|patient|donor|individ", tolower(colnames(df)))

  ## ... exclude numerical dates/age/year
  is.bad2 <- grepl("ratio|year|month|day|^age$|^efs|^dfs|surv|follow", tolower(colnames(df)))

  is.num <- sapply(df, class) == "numeric"
  is.bad2 <- (is.bad2 & is.num) ## no numeric

  ## ... exclude any sample ID coded in columns...
  is.bad3 <- apply(df, 2, function(x) {
    mean(grepl("^sample|patient|donor|individ",
      x[!is.na(x)],
      ignore.case = TRUE
    )) > 0.8
  })

  is.bad <- (is.bad2 | is.bad3)


  ## auto-determine which are factors
  is.factor <- apply(df, 2, is.categorical)

  n.unique <- apply(df, 2, function(x) length(unique(setdiff(x, c(NA, "NA", "")))))
  n.notna <- apply(df, 2, function(x) length(x[!is.na(x)]))
  is.id <- (n.unique > 0.9 * n.notna)

  is.factor2 <- (!is.bad & is.factor & !is.id & n.unique >= min.ncat & n.unique <= max.ncat)
  is.factor2

  ## take reduced matrix
  df1 <- df[, which(is.factor2), drop = FALSE]
  nlevel <- apply(df1, 2, function(x) length(unique(x)))
  nchars <- apply(df1, 2, function(x) max(nchar(iconv(x, "latin1", "ASCII", sub = ""))))
  df1 <- df1[, order(nlevel, -nchars), drop = FALSE]


  if (remove.dup && ncol(df1) > 1) {
    i <- 1
    j <- 2
    is.dup <- rep(FALSE, ncol(df1))
    for (i in 1:(ncol(df1) - 1)) {
      is.dup[i] <- FALSE
      for (j in (i + 1):ncol(df1)) {
        is.dup[i] <- is.dup[i] || all(rowSums(table(df1[, i], df1[, j]) != 0) == 1)
      }
    }
    df1 <- df1[, which(!is.dup), drop = FALSE]
  }
  return(colnames(df1))
}


## Determine if the gene type is mouse
#' @title Determine organism from count matrix
#'
#' @description Determines if count data is from human or mouse based on gene identifiers.
#'
#' @param pgx A pgx object with the pgx$organism information and the pgx$counts slot for the
#' guessing approach.
#' @param capitalise logical: by default FALSE. Parameter to capitalise the first letter of the
#' specie if mouse or human.
#' @details This function retreives the pgx$organism slot. If it is not found, then it examines
#' the gene identifiers in the row names of a count matrix to determine if the data is from human
#' or mouse (main organism supported in the old playbase version). It checks if the identifiers
#' match common patterns found in mouse genes, like "rik", "loc", "orf". If more than 20% match
#'  these mouse patterns, it assigns the organism as "mouse". Otherwise it assigns "human".
#'
#' The function calculates the fraction of row names that do NOT match the mouse gene patterns.
#' If this fraction is >0.8, it assigns "human". This relies on the assumption that human data
#' will tend to have more uppercase ENSEMBL identifiers.
#'
#' @return Character string the organism.
#' @export
pgx.getOrganism <- function(pgx, capitalise = FALSE) {
  pgx.counts <- pgx$counts
  if (!is.null(pgx$organism)) {
    org <- pgx$organism
  } else {
    rownames.counts <- grep("^rik|^loc|^orf", rownames(pgx$counts),
      value = TRUE,
      ignore.case = TRUE, invert = TRUE
    )
    cap.fraction <- mean(grepl("^[A-Z][a-z]+", rownames.counts), na.rm = TRUE)
    is.mouse <- (cap.fraction > 0.8)
    org <- ifelse(is.mouse, "mouse", "human")
  }

  if (capitalise && org %in% c("mouse", "human")) {
    if (org == "mouse") org <- "Mouse"
    if (org == "human") org <- "Human"
  }
  return(org)
}


#' @title Get Levels of Group Variables
#'
#' @description This function retrieves the levels of group variables in a data frame.
#'
#' @param Y A data frame containing the data to be analyzed.
#'
#' @details The function takes a data frame `Y` as input and searches for
#' columns that represent group variables. Group variables are identified as
#' columns that have a name that does not contain "title", "name", "sample", or
#' "patient" and have a majority of non-unique values.
#' Numeric columns are excluded from the search.
#' The levels of the identified group variables are then extracted and returned as a character vector.
#'
#' @return A character vector representing the levels of the group variables in the input data frame.
#'
#' @export
getLevels <- function(Y) {
  yy <- Y[, grep("title|name|sample|patient", colnames(Y), invert = TRUE), drop = FALSE] ## NEED RETHINK!!!!
  is.grpvar <- apply(yy, 2, function(y) max(table(y)) > 1)
  is.numeric <- apply(yy, 2, function(y) (length(table(y)) / length(y)) > 0.5)
  is.grpvar <- is.grpvar & !is.numeric
  yy <- yy[, is.grpvar, drop = FALSE]


  levels <- lapply(1:ncol(yy), function(i) unique(paste0(colnames(yy)[i], "=", yy[, i])))
  levels <- sort(unlist(levels))
  return(levels)
}


#' @title Select samples from selected levels
#'
#' @description Selects rows from a data frame corresponding to selected factor levels.
#'
#' @param Y Data frame containing sample data.
#' @param levels Character vector of selected factor level descriptions in format 'factor=level'.
#'
#' @details This function takes a data frame Y containing sample data, and a character vector
#' levels defining a subset of factor levels to select. It parses the level descriptions,
#' identifies matching samples in Y, and returns the rownames of selected samples.
#'
#' The levels should be strings in format 'factor=level' defining the
#' factor name and level value.  Samples in Y matching any of the
#' provided levels will be selected. Note: across factors samples as
#' intersected ('and' operator), within a factor samples are combined
#' ('or' operator). Not very logical but this is what people in practice want.
#'
#' @return A character vector of selected sample names.
#'
#' @export
selectSamplesFromSelectedLevels <- function(Y, levels) {
  if (is.null(levels) || all(levels == "")) {
    return(rownames(Y))
  }

  # fill ="" will (unfortunately) still return NA when level is "NA"... which crashes when phenotype is ""
  pheno <- data.table::tstrsplit(levels, "=", keep = 1) %>%
    unlist()
  ptype <- data.table::tstrsplit(levels, "=", keep = 2) %>%
    unlist()
  # force replace "NA" by NA
  ptype[ptype == "NA"] <- NA
  ##  sel <- rep(FALSE, nrow(Y))
  sel <- rep(TRUE, nrow(Y))
  for (ph in unique(pheno)) {
    # ph = pheno[1]
    k <- which(pheno == ph)
    ##    sel <- sel | (Y[, ph] %in% ptype[k])
    sel <- sel & (Y[, ph] %in% ptype[k])
  }

  return(rownames(Y)[which(sel)])
}


#' @title Get Gene Information
#'
#' @description This function retrieves gene information from the biomaRt package.
#'
#' @param eg A character vector of Entrez Gene IDs for which to retrieve information.
#' @param fields A character vector of fields to retrieve, with default values of
#' "symbol", "name", "alias", "map_location", and "summary".
#'
#' @details The function takes a character vector of Entrez Gene IDs `eg` and a
#' character vector of fields `fields` as input.
#' The function uses the `getGene` function from the `biomaRt` package to retrieve
#' the specified fields for each Entrez Gene ID.
#' The resulting information is returned as a named list, where each element
#' corresponds to one of the specified fields and contains a character vector of the retrieved values.
#'
#' @return A named list containing the retrieved gene information for each of the specified fields.
#'
#' @export
getMyGeneInfo <- function(eg, fields = c("symbol", "name", "alias", "map_location", "summary")) {
  info <- lapply(fields, function(f) biomaRt::getGene(eg, fields = f)[[1]])
  names(info) <- fields
  info <- lapply(info, function(x) ifelse(length(x) == 3, x[[3]], "(not available)"))
  info <- sapply(info, paste, collapse = ",")

  return(info)
}


#' Get human gene information from annotation packages
#'
#' @param eg Character vector of gene symbols
#' @param as.link Logical indicating whether to return symbols as hyperlinks. Default is TRUE.
#'
#' @return Named character vector with gene information.
#'
#' @details This function retrieves basic gene information for human genes from Bioconductor annotation packages.
#' It gets the symbol, name, map location, OMIM IDs, KEGG IDs and GO terms for each input gene symbol.
#'
#' If as.link is TRUE, the gene symbols are returned as hyperlinks to GeneCards and OMIM IDs are returned as hyperlinks.
#'
#' @examples
#' \dontrun{
#' genes <- c("TP53", "BRCA1", "ABCD1")
#' info <- getHSGeneInfo(genes)
#' }
#' @export
getHSGeneInfo <- function(eg, as.link = TRUE) {
  env.list <- c(
    "symbol" = org.Hs.eg.db::org.Hs.egSYMBOL,
    "name" = org.Hs.eg.db::org.Hs.egGENENAME,
    "map_location" = org.Hs.eg.db::org.Hs.egMAP,
    "OMIM" = org.Hs.eg.db::org.Hs.egOMIM,
    "KEGG" = org.Hs.eg.db::org.Hs.egPATH,
    "GO" = org.Hs.eg.db::org.Hs.egGO
  )

  info <- lapply(env.list, function(env) AnnotationDbi::mget(eg, envir = env, ifnotfound = NA)[[1]])
  names(info) <- names(env.list)
  gene.symbol <- toupper(AnnotationDbi::mget(as.character(eg),
    envir = org.Hs.eg.db::org.Hs.egSYMBOL
  ))[1]
  info[["symbol"]] <- gene.symbol

  ## create link to GeneCards
  if (as.link) {
    genecards.link <- "<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=GENE' target='_blank'>GENE</a>"
    info[["symbol"]] <- gsub("GENE", info[["symbol"]], genecards.link)
  }

  ## create link to OMIM
  if (as.link) {
    omim.link <- "<a href='https://www.omim.org/entry/OMIM' target='_blank'>OMIM</a>"
    info[["OMIM"]] <- sapply(info[["OMIM"]], function(x) gsub("OMIM", x, omim.link))
  }

  ## create link to KEGG
  kegg.link <- "<a href='https://www.genome.jp/kegg-bin/show_pathway?map=hsaKEGGID&show_description=show' target='_blank'>KEGGNAME (KEGGID)</a>"
  for (i in 1:length(info[["KEGG"]])) {
    kegg.id <- info[["KEGG"]][[i]]
    kegg.id <- setdiff(kegg.id, NA)
    if (length(kegg.id) > 0) {
      kegg.name <- AnnotationDbi::mget(kegg.id, envir = KEGG.db::KEGGPATHID2NAME, ifnotfound = NA)[[1]]
      if (!is.na(kegg.name) && as.link) {
        info[["KEGG"]][[i]] <- gsub("KEGGNAME", kegg.name, gsub("KEGGID", kegg.id, kegg.link))
      } else {
        info[["KEGG"]][[i]] <- kegg.name
      }
    }
  }

  ## create link to GO
  if (!is.na(info[["GO"]][1])) {
    go.evidence <- c("EXP", "IDA", "IPI", "IMP", "IGI", "IEP")
    amigo.link <- "<a href='http://amigo.geneontology.org/amigo/term/GOID' target='_blank'>GOTERM (GOID)</a>"
    sel <- which(sapply(info[["GO"]], "[[", 2) %in% go.evidence &
      sapply(info[["GO"]], "[[", 3) %in% c("BP"))
    sel
    info[["GO"]] <- info[["GO"]][sel]

    ## sometimes GO.db is broken...
    suppressWarnings(try.out <- try(AnnotationDbi::Term(AnnotationDbi::mget("GO:0000001",
      envir = GO.db::GOTERM,
      ifnotfound = NA
    )[[1]])))
    go.ok <- (class(try.out) != "try-error")
    if (go.ok && length(sel) > 0) {
      i <- 1
      for (i in 1:length(info[["GO"]])) {
        go_id <- info[["GO"]][[i]][[1]]
        go_term <- AnnotationDbi::Term(AnnotationDbi::mget(go_id, envir = GO.db::GOTERM, ifnotfound = NA)[[1]])
        if (as.link) {
          info[["GO"]][[i]] <- gsub("GOTERM", go_term, gsub("GOID", go_id, amigo.link))
        } else {
          info[["GO"]][[i]] <- go_term
        }
      }
    } else {
      info[["GO"]] <- NULL
    }
  }

  return(info)
}


#' @describeIn pgx.getGeneSetCollections Get the gene families.
#' @export
pgx.getGeneFamilies <- function(genes, min.size = 10, max.size = 500) {
  read.gmt <- function(gmt.file, dir = NULL, add.source = FALSE, nrows = -1) {
    f0 <- gmt.file
    if (strtrim(gmt.file, 1) == "/") dir <- NULL
    if (!is.null(dir)) f0 <- paste(sub("/$", "", dir), "/", gmt.file, sep = "")

    gmt <- utils::read.csv(f0, sep = "!", header = FALSE, comment.char = "#", nrows = nrows)[, 1]
    gmt <- as.character(gmt)

    gmt <- sapply(gmt, strsplit, split = "\t")
    names(gmt) <- NULL
    gmt.name <- sapply(gmt, "[", 1)
    gmt.source <- sapply(gmt, "[", 2)
    gmt.genes <- sapply(gmt, function(x) {
      if (length(x) < 3) {
        return("")
      }
      paste(x[3:length(x)], collapse = " ")
    })

    gset <- strsplit(gmt.genes, split = "[ \t]")
    gset <- lapply(gset, function(x) setdiff(x, c("", "NA", NA)))
    names(gset) <- gmt.name
    if (add.source) {
      names(gset) <- paste0(names(gset), " (", gmt.source, ")")
    }

    return(gset)
  }

  ## -----------------------------------------------------------------------------
  ## Gene families
  ## -----------------------------------------------------------------------------

  families <- list()
  families[["<all>"]] <- genes ## X is sorted

  gmt.kea <- read.gmt(playdata::get_file("kinase_substrates_kea.gmt"))
  gmt.chea <- read.gmt(playdata::get_file("tf_targets_chea.gmt"))
  families[["Kinases (KEA)"]] <- names(gmt.kea)
  families[["Transcription factors (ChEA)"]] <- names(gmt.chea)

  ## Read standard HGNC gene families (www.genefamilies.org)
  gmt.hgnc <- read.gmt(playdata::get_file("hgnc_genefamilies_EDITED.gmt"))

  gmt.hgnc.size <- sapply(gmt.hgnc, length)
  gmt.hgnc <- gmt.hgnc[which(gmt.hgnc.size >= 50 & gmt.hgnc.size <= 1000)]

  names(gmt.hgnc) <- paste0(names(gmt.hgnc), " (HGNC)")

  families <- c(families, gmt.hgnc)

  families[["Interleukins (IL)"]] <- genes[grep("^IL[1-9]", genes)]
  families[["Chemokines"]] <- genes[grep("CCL|CCR|CXCR|CXCL|XCL|CX3", genes)]
  families[["Ribosomal proteins"]] <- genes[grep("^RPS|^RPL", genes)]
  families[["Ribosomal (mitochondrial)"]] <- genes[grep("^MRPL|^MRPS", genes)]

  families[["G-protein family"]] <- genes[grep("^GN[ABG]|^GPBAR|^GPER|^FZD", genes)]
  families[["Heatshock proteins"]] <- genes[grep("^HSP", genes)]
  families[["Integrin family"]] <- genes[grep("^ITG", genes)]
  families[["MAPK family"]] <- genes[grep("^MAP[1-9]|^MAPK", genes)]

  families[["Myosins"]] <- genes[grep("^MYO[1-9]", genes)]
  families[["Protein phosphatase (PPP)"]] <- genes[grep("^PPP", genes)]
  families[["PTP family"]] <- genes[grep("^PTP", genes)]
  families[["Small nucleolar RNA"]] <- genes[grep("^SNOR", genes)]
  families[["Toll-like receptors"]] <- genes[grep("^TLR", genes)]

  families[["EIF factors"]] <- genes[grep("^EIF", genes)]
  families[["GPR proteins"]] <- genes[grep("^GPR", genes)]
  families[["CDCC proteins"]] <- genes[grep("^CDCC", genes)]
  families[["KIAA proteins"]] <- genes[grep("^KIAA", genes)]
  families[["TRIM proteins"]] <- genes[grep("^TRIM", genes)]
  families[["TNF proteins"]] <- genes[grep("^TNF", genes)]
  families[["SLC proteins"]] <- genes[grep("^SLC", genes)]
  families[["ZBTB proteins"]] <- genes[grep("^ZBTB", genes)]
  families[["TMEM family"]] <- genes[grep("^TMEM", genes)]
  families[["STAT family"]] <- genes[grep("^STAT", genes)]
  families[["DNA/RNA polymerases"]] <- genes[grep("^POL", genes)]
  families[["Proteasome"]] <- genes[grep("^PSM", genes)]
  families[["IFN/IFIT family"]] <- genes[grep("^IFN|^IFIT", genes)]
  families[["Nuclear receptors"]] <- genes[grep("^NR[0-9]|^RXR|^ESR|^PGR$|^AR$|^HNF4|^ROR|^PPAR|^THR|^VDR", genes)]
  families[["Cytochrome family"]] <- genes[grep("^CYP|^CYB|^CYC|^COX|^COA", genes)]
  families[["Micro RNA"]] <- genes[grep("^MIR", genes)]

  ## sort all
  families <- lapply(families, function(f) intersect(f, genes))

  ## ----------- filter on size
  nsize <- sapply(families, length)
  sel <- which(nsize >= min.size & nsize < max.size)
  sel <- unique(c(which(names(families) == "<all>"), sel))
  families <- families[sel]

  return(families)
}


#' @title Get Gene Set Collections
#'
#' @description
#' Extracts gene set collections from a list of gene sets based on size filters.
#'
#' @param gsets A list of gene sets
#' @param min.size The minimum gene set size to include. Default 10.
#' @param max.size The maximum gene set size to include. Default 500.
#'
#' @details This function takes a list of gene sets and extracts collections
#' of gene sets that fall within a specified size range, defined by min.size
#' and max.size parameters. It removes any gene sets that are smaller than
#' `min.size` or larger than `max.size`. This allows filtering gene sets into size-based
#' collections for different types of analyses.
#'
#' @return A list containing the extracted gene set collections.
#'
#' @export
pgx.getGeneSetCollections <- function(gsets = rownames(playdata::GSETxGENE)) {
  ## -----------------------------------------------------------------------------
  ## Gene set collections
  ## -----------------------------------------------------------------------------

  collections <- list()

  ## ----------- add main collections from gene set prefixes
  gsets.db <- sub("_.*", "", gsets)
  gsets.db <- sub(":.*", "", gsets.db)
  gsets.groups <- tapply(gsets, gsets.db, list)
  collections <- c(collections, gsets.groups)
  collections[["<all>"]] <- gsets
  return(collections)
}


## -----------------------------------------------------------------------------
## Generic module functions
## -----------------------------------------------------------------------------


#' Filter probes from gene expression data
#'
#' @param genes Gene expression data.frame with probes as row names.
#' @param gg Character vector of probes to filter for.
#'
#' @return Character vector of matching probe names.
#'
#' @details This function filters a gene expression data.frame to return only probes matching the input vector.
#' It checks the probe name, short probe name, and gene name for matches.
#'
#' @examples
#' \dontrun{
#' data <- utils::read.csv("expression.csv", row.names = 1)
#' probes <- c("FOO", "BAR")
#' matches <- filterProbes(data, probes)
#' }
#' @export
filterProbes <- function(genes, gg) {
  ## check probe name, short probe name or gene name for match
  p0 <- (toupper(sub(".*:", "", rownames(genes))) %in% toupper(gg))
  p1 <- (toupper(rownames(genes)) %in% toupper(gg))
  p2 <- (toupper(as.character(genes$gene_name)) %in% toupper(gg))
  if ("human_ortholog" %in% colnames(genes)) {
    p3 <- (toupper(as.character(genes$human_ortholog)) %in% toupper(gg))
  } else {
    p3 <- rep(FALSE, nrow(genes))
  }

  # Ensure all p* are valids
  p_list <- list(p0, p1, p2, p3)
  p_list <- p_list[sapply(p_list, length) > 0]

  # Combine list using OR operator
  jj <- which(Reduce("|", p_list))
  if (length(jj) == 0) {
    return(NULL)
  }
  return(rownames(genes)[jj])
}


#' Rename rownames of counts matrix by annotation table
#'
#' @param counts Numeric matrix of counts, with genes/probes as rownames.
#' @param annot_table Data frame with rownames matching counts and annotation columns.
#' @param new_id_col Column name in annot_table containing new identifiers. Default 'symbol'.
#'
#' @return Matrix with rownames changed to values from annot_table.
#' Duplicate new rownames are summed.
#'
#' @details Renames rownames of counts matrix using an annotation data frame.
#' Looks up the `new_id_col` in the annot_table and replaces counts rownames.
#' Handles special cases like missing values.
#' Sums duplicate rows after renaming.
#'
#' @export
rename_by <- function(counts, annot_table, new_id_col = "symbol") {
  symbol <- annot_table[rownames(counts), new_id_col]

  # Guard agaisn human_hommolog == NA
  if (all(is.na(symbol))) {
    symbol <- annot_table[rownames(counts), "symbol"]
  }

  # Sum columns of rows with the same gene symbol
  if (is.matrix(counts) | is.data.frame(counts)) {
    rownames(counts) <- symbol
    return(counts[!rownames(counts) %in% c("", "NA"), , drop = FALSE])
  } else {
    return(symbol)
  }
}


#' Compute feature scores
#'
#' @param X Numeric matrix of expression data, with genes in rows and samples in columns.
#' @param Y Data frame containing sample metadata.
#' @param features Named list of gene sets (character vectors of gene names) to compute scores for.
#'
#' @return A numeric matrix of feature scores, with gene sets in rows and samples in columns.
#'
#' @details This function computes a score for each specified gene set that
#' represents its average internal correlation within each sample group.
#'
#' For each sample group, it takes the top 100 most variable genes in the gene set.
#' It then calculates the correlation distance matrix between these genes within that sample group.
#' The score is the mean of the upper triangle of the distance matrix.
#'
#' Higher scores indicate the gene set has high internal correlation and is
#' coherently activated within that sample group.
#'
#' @examples
#' \dontrun{
#' # Generate random data
#' X <- matrix(rnorm(100 * 50), 100, 50)
#' Y <- data.frame(group = sample(letters[1:5], 50, replace = TRUE))
#' features <- list(
#'   feat1 = sample(rownames(X), 50),
#'   feat2 = sample(rownames(X), 50)
#' )
#'
#' scores <- computeFeatureScore(X, Y, features)
#' }
#' @export
computeFeatureScore <- function(X, Y, features) {
  sdx <- apply(X, 1, stats::sd)
  names(sdx) <- rownames(X)
  S <- matrix(NA, nrow = length(features), ncol = ncol(Y))
  rownames(S) <- names(features)
  colnames(S) <- colnames(Y)
  for (k in 1:ncol(Y)) {
    grp <- Y[colnames(X), k]
    grp <- as.character(grp)
    score <- rep(NA, length(features))
    names(score) <- names(features)
    i <- 1
    for (i in 1:length(features)) {
      pp <- features[[i]]

      pp <- Matrix::head(pp[order(-sdx[pp])], 100)
      mx <- t(apply(X[pp, ], 1, function(x) tapply(x, grp, mean)))

      D <- 1 - stats::cor(mx, use = "pairwise")
      diag(D) <- NA
      score[i] <- mean(D, na.rm = TRUE)
    }
    S[, k] <- score
  }
  if (is.null(S)) {
    return(NULL)
  }
  return(S)
}


## -----------------------------------------------------------------------------
## Generic helper functions
## -----------------------------------------------------------------------------


#' Check if a vector contains date values
#'
#' @param x A vector to check if it contains date values
#'
#' @return Logical. TRUE if the vector contains valid date values, FALSE otherwise.
#'
#' @details This function checks if a vector contains valid date values that can be converted to Date class.
#' It first replaces any "NA" values with actual NA values.
#' It then tries converting the vector to Date class using various date formats.
#' If none of the dates throw an error, it returns TRUE indicating the vector contains valid dates.
#'
#' @examples
#' \dontrun{
#' x <- c("2020-01-01", "01/15/2021", "foo")
#' is.Date(x) # TRUE
#' }
#' @export
is.Date <- function(x) {
  ## From https://stackoverflow.com/questions/18178451/is-there-a-way-to-check-if-a-column-is-a-date-in-r
  if (!all(is.na(as.Date(
    as.character(x),
    format = c("%d/%m/%Y", "%d-%m-%Y", "%Y/%m/%d", "%Y-%m-%d")
  )))) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#' @title Calculate group means
#'
#' @description
#' Calculates the column means within groups defined by a grouping variable.
#'
#' @param mat Numeric matrix with columns as samples.
#' @param group Grouping vector or factor.
#' @param FUN Function for aggregation. Default is mean.
#'
#' @return Matrix with group means.
#'
#' @details This function calculates the column means of \code{X} within groups
#'  defined by \code{y}. It calculates the mean for each column within each
#' group. The output is a matrix with rows corresponding to groups and columns
#' corresponding to samples.
#'
#' @examples
#' \dontrun{
#' mat <- matrix(rnorm(100), ncol = 10)
#' groups <- gl(2, 5)
#' means <- averageByGroup(mat, groups)
#' }
#'
#' @export
averageByGroup <- function(mat, group, FUN = mean) {
  out <- do.call(cbind, tapply(
    1:ncol(mat), group,
    function(i) rowMeans(mat[, i, drop = FALSE])
  ))
  return(out)
}


#' @title Make Acronym from String
#'
#' @description This function creates an acronym from a character string.
#'
#' @param x A character string from which to create an acronym.
#'
#' @details The function takes a character string `x` as input and splits it into words using `_`, `-`, and space as delimiters.
#' If the resulting list of words contains only one word, the first two characters of the word are returned.
#' Otherwise, the first character of each word is extracted, capitalized, and concatenated to form the acronym.
#'
#' @return A character string representing the acronym created from the input string.
#'
#' @examples
#' \dontrun{
#' # example code
#' makeAcronym("United Nations")
#' makeAcronym("World Health Organization")
#' }
#' @export
makeAcronym <- function(x) {
  xp <- strsplit(x, split = "[_ -]")
  sapply(xp, function(s) {
    if (length(s) == 1) {
      return(substring(s, 1, 2))
    }
    toupper(paste(substring(s, 1, 1), collapse = ""))
  })
}


#' @title Relevel Factor by First Occurrence
#'
#' @description This function relevels a factor by the order of first
#' occurrence of its levels.
#'
#' @param f A factor to be releveled.
#'
#' @details The function takes a factor `f` as input and relevels it by
#' the order of first occurrence of its levels in the data.
#' The new factor is created with the same levels as the original factor,
#' but with the levels reordered according to their first appearance in the data.
#'
#' @return A factor with the same levels as the input factor, but with the
#' levels reordered according to their first appearance in the data.
#'
#' @examples
#' \dontrun{
#' # example code
#' f <- factor(c("b", "a", "c", "b", "a"))
#' relevelFactorFirst(f)
#' }
#' @export
relevelFactorFirst <- function(f) {
  factor(f, levels = f[!duplicated(f)])
}


#' @title Extreme correlation
#'
#' @description Compute extreme (high and low) correlation between a
#' query signature and reference dataset
#'
#' @param query_sig Numeric vector containing the query signature values
#' for each gene
#' @param ref_set Numeric matrix containing the reference dataset
#' (samples in columns, genes in rows)
#' @param n Integer specifying number of top extreme correlations to return on each side. Default is 200.
#'
#' @details This function takes a query signature vector and a reference dataset matrix as input.
#' It calculates Pearson correlation between the query signature and each column of the reference dataset.
#' The function returns the most extreme (high and low) correlations, taking the top \code{n} correlations on each side.
#' This allows extracting the samples in the reference data that have the strongest positive and negative correlation to the query signature.
#'
#' @return A numeric vector containing the top \code{n} most extreme correlation coefficients (positive and negative).
#'
#' @examples
#' \dontrun{
#' # Generate random query signature
#' query_sig <- rnorm(1000)
#' names(query_sig) <- paste0("GENE", 1:1000)
#'
#' # Generate random reference dataset
#' ref_set <- matrix(rnorm(1000 * 50), ncol = 50)
#' rownames(ref_set) <- names(query_sig)
#'
#' # Get top 200 extreme correlations
#' extreme_cor <- extremeCorrelation(query_sig, ref_set, n = 200)
#' }
#' @export
extremeCorrelation <- function(query_sig, ref_set, n = 200) {
  gg <- intersect(rownames(ref_set), names(query_sig))
  if (n > 0) {
    gg <- gg[unique(c(Matrix::head(order(query_sig), n), utils::head(order(-query_sig), n)))]
  }
  rho <- stats::cor(ref_set[gg, ], query_sig[gg], use = "pairwise")
  rho <- rho[order(-rowMeans(rho**2, na.rm = TRUE)), , drop = FALSE]
  if (NCOL(rho) == 1) rho <- rho[, 1]
  return(rho)
}

#' Convert any gene alias to official gene symbol
#' @description
#' This function takes a character vector of gene aliases and converts them to HUGO gene symbols.
#' The organism can be specified as either human or mouse, or the function can attempt to determine
#' the organism based on the input gene aliases. Unmapped gene aliases can be kept as their original values
#' or replaced with NA.
#'
#' @param s A character vector of gene aliases to be converted to HUGO gene symbols.
#' @param org A character string specifying the organism, either "hs" for human or "mm" for mouse.
#' If not specified, the function will attempt to determine the organism based on the input gene aliases.
#' @param na.orig A logical value indicating whether to keep the original values of
#' gene aliases that could not be mapped to HUGO gene symbols. If TRUE, the original values will be kept,
#' otherwise they will be replaced with NA.
#'
#' @return A character vector of HUGO gene symbols corresponding to the input gene aliases.
#'
#' @examples
#' \dontrun{
#' symbols <- c("TP53", "P53", "Cd19", NA)
#' alias2hugo(symbols)
#' alias2hugo(symbols, org = "mm")
#' }
#' @export
alias2hugo <- function(s, org = NULL, na.orig = TRUE) {
  ## determine if we deal with human or mouse
  hs.symbol <- unlist(as.list(org.Hs.eg.db::org.Hs.egSYMBOL))
  mm.symbol <- unlist(as.list(org.Mm.eg.db::org.Mm.egSYMBOL))
  if (is.null(org)) {
    is.human <- mean(s %in% hs.symbol, na.rm = TRUE) >= mean(s %in% mm.symbol, na.rm = TRUE)
    org <- ifelse(is.human, "hs", "mm")
  }
  org

  nna <- which(!is.na(s) & s != "" & s != " ")
  s1 <- trimws(s[nna])
  hugo <- NULL
  if (org == "hs") {
    hugo <- suppressWarnings(
      limma::alias2SymbolTable(s1, species = "Hs")
    )
  } else if (org == "mm") {
    hugo <- suppressWarnings(
      limma::alias2SymbolTable(s1, species = "Mm")
    )
  } else {
    stop("[alias2hugo] invalid organism")
  }
  jj <- which(is.na(hugo))
  if (na.orig && length(jj)) hugo[jj] <- s1[jj]
  hugo0 <- rep(NA, length(s))
  hugo0[nna] <- hugo
  return(hugo0)
}

#' @describeIn breakstring2 Breaks a character string into substrings of a
#' specified length, separated by a specified line break character.
#' @param force logical not used
#' @export
breakstring <- function(s, n, nmax = 999, force = FALSE, brk = "\n") {
  if (is.na(s)) {
    return(NA)
  }
  s <- substring(as.character(s), 1, nmax)
  if (nchar(s) < n) {
    return(s)
  }
  b <- substring(s, 1, n)
  n1 <- n + 1
  for (i in 1:10) {
    if (n1 > nchar(s)) break
    b1 <- substring(s, n1, n1 + n)
    b <- paste0(b, brk, b1)
    n1 <- n1 + n + 1
  }
  return(b)
}


#' Break up a long string into multiple lines
#'
#' @param s The input string
#' @param n The maximum line length
#' @param brk The string to insert for line breaks. Default is "\\n".
#' @param nmax The maximum length to process. Longer strings are truncated. Default is 999.
#'
#' @return A string with line breaks inserted.
#'
#' @details This function takes a long string \code{s} and inserts line breaks \code{brk} every \code{n} characters.
#' It first truncates the input string to a maximum length \code{nmax} if longer.
#' The string is then split into words, and words are added to each line until it reaches length \code{n}.
#' A line break \code{brk} is inserted, and the process continues until all words have been processed.
#'
#' @export
breakstring2 <- function(s, n, brk = "\n", nmax = 999) {
  if (is.na(s)) {
    return(NA)
  }
  s <- substring(as.character(s), 1, nmax)
  if (is.na(s)) {
    return(NA)
  }
  if (nchar(s) < n) {
    return(s)
  }
  a <- ""
  words <- paste0(strsplit(s, split = " ")[[1]], " ")
  words
  it <- 0
  len1 <- sum(sapply(words, nchar))
  len1
  while (len1 > n && it < 100 && length(words) > 1) {
    len <- cumsum(sapply(words, nchar))
    len
    k <- which(len > n)[1]
    k
    if (k == 1) k <- 2
    a <- paste0(a, paste0(words[1:(k - 1)], collapse = ""), brk)
    a
    words <- words[k:length(words)]
    words
    it <- it + 1
    len1 <- sum(sapply(words, nchar))
    len1
  }
  a <- paste0(a, paste0(words, collapse = ""), brk)
  a <- sub(paste0(brk, "$"), "", gsub(paste0(" ", brk), brk, a))
  return(a)
}


#' @describeIn shortstring0 The shortstring function truncates a string to a specified length.
#' @export
shortstring <- function(s, n, dots = 1) {
  sapply(s, shortstring0, n = n, dots = dots)
}


#' @title Truncate a string to a maximum length
#'
#' @description
#' Truncates a string to a maximum length by replacing the middle part with ellipses.
#'
#' @param s The string to truncate.
#' @param n The maximum length of the truncated string.
#' @param dots The number of characters to devote to the ellipses.
#'
#' @details
#' This function truncates long strings to a specified maximum length
#' by replacing the middle part with ellipses.
#'
#' It first converts the input string to UTF-8 encoding and removes any HTML entities.
#' If the string is already shorter than the maximum length, it is returned unchanged.
#' Otherwise, it truncates the beginning and end of the string to fit within the max length when concatenated with the ellipses.
#'
#' The number of characters to use for the ellipses is set by the \code{dots}
#' parameter. By default, it uses one dot per 10 characters of max length.
#'
#' @return The truncated string.
#'
#' @examples
#' \dontrun{
#' shortstring0("This is a very long string", 20)
#' # Returns "This is a ...string"
#'
#' shortstring0("Short string", 20)
#' # Returns "Short string"
#' }
#' @export
shortstring0 <- function(s, n, dots = 1) {
  s0 <- iconv(as.character(s), to = "UTF-8")
  s0 <- gsub("[&].*[;]", "", s0) ## HTML special garbage...
  jj <- which(nchar(s0) > n)
  if (length(jj) == 0) {
    return(s0)
  }
  s <- s0[jj]
  if (dots < 1) {
    n1 <- ceiling(dots * n)
    s1 <- substring(s, 1, n1)
    s2 <- substring(s, nchar(s) - n + nchar(s1), nchar(s))
    aa <- paste0(s1, "...", s2)
  } else {
    aa <- paste0(substring(s, 1, n), "...")
  }
  s0[jj] <- aa
  s0
}


#' @title Sort Data Frame by P-value Column
#'
#' @description This function sorts a data frame by a column containing p-values.
#'
#' @param x A data frame representing the input data.
#' @param p.col An optional character string specifying the name of the column
#' containing p-values.
#'
#' @details The function takes a data frame `x` as input and searches
#' for a column containing p-values. If the `p.col` parameter is not provided,
#' the function searches for a column with a name containing "p.value", "p",
#' "p-val", or "pval" (case-insensitive). If multiple columns match, the first one is used.
#'
#' The data frame is then sorted in ascending order by the p-value column and returned.
#'
#' @return A data frame representing the input data sorted by the p-value column.
#'
#' @examples
#' \dontrun{
#' # example code
#' x <- data.frame(a = c(1, 2, 3), p.value = c(0.05, 0.01, 0.1))
#' psort(x)
#' }
#' @export
psort <- function(x, p.col = NULL) {
  j <- grep("p.value|^p$|p-val|pval", tolower(colnames(x)))[1]
  x[order(x[, j]), ]
}


#' @title Tidy Data Frame
#'
#' @description This function tidies a data frame by removing columns with all
#' missing values, converting character columns to numeric or factor columns,
#' and trimming white spaces.
#'
#' @param Y A data frame representing the input data.
#'
#' @details The function takes a data frame `Y` as input and performs the
#' following steps:
#' 1. Remove columns with all missing values.
#' 2. Replace character columns with only "NA" values with actual missing values (NA).
#' 3. Remove columns with all missing values again.
#' 4. Trim white spaces from the beginning and end of character columns.
#' 5. Convert character columns to numeric columns if possible.
#' 6. Convert character columns to factor columns if they have less than or equal to
#'  3 unique values or if their names contain certain keywords (e.g., "batch",
#' "replicate", "type", "cluster", "group").
#'
#' The resulting tidied data frame is returned.
#'
#' @return A data frame representing the tidied input data.
#'
#' @export
tidy.dataframe <- function(Y) {
  Y <- Y[, which(colMeans(is.na(Y)) < 1), drop = FALSE]
  Y <- apply(Y, 2, function(x) sub("^NA$", NA, x)) ## all characters
  Y <- Y[, which(colMeans(is.na(Y)) < 1), drop = FALSE]
  Y <- apply(Y, 2, function(x) gsub("^[ ]*|[ ]*$", "", x))
  suppressWarnings(num.Y <- apply(Y, 2, function(x) as.numeric(as.character(x))))
  is.numeric <- (0.8 * colMeans(is.na(num.Y)) <= colMeans(is.na(Y)))
  nlevel <- apply(Y, 2, function(x) length(unique(x)))
  is.factor <- (!is.numeric | (is.numeric & nlevel <= 3))
  is.factor <- (is.factor | grepl("batch|replicat|type|clust|group", colnames(Y)))
  new.Y <- data.frame(Y, check.names = FALSE)
  new.Y[, which(is.numeric)] <- num.Y[, which(is.numeric), drop = FALSE]
  for (i in which(is.numeric)) new.Y[[i]] <- num.Y[, i]
  for (i in which(is.factor)) new.Y[[i]] <- factor(as.character(new.Y[, i]))
  new.Y <- data.frame(new.Y, check.names = FALSE)
  return(new.Y)
}


#' @title Get Classes of Parameters
#'
#' @description This function returns the classes of the parameters in a data frame.
#'
#' @param A A data frame representing the input data.
#'
#' @details The function takes a data frame `A` as input and applies the `class`
#' function to each column of the data frame using the `sapply` function.
#' The resulting vector of classes is returned.
#'
#' @return A character vector representing the classes of the parameters in
#' the input data frame.
#'
#' @examples
#' \dontrun{
#' # example code
#' A <- data.frame(x = c(1, 2, 3), y = c("a", "b", "a"))
#' param.class(A)
#' }
#' @export
param.class <- function(A) sapply(tidy.dataframe(A), class)


#' @describeIn isanumber The function tests if the input \code{x} is of numeric type.
#' @export
is.num <- function(y) {
  suppressWarnings(numy <- as.numeric(as.character(y)))
  t1 <- !all(is.na(numy)) && is.numeric(numy)
  t2 <- length(unique(y)) > 0.33 * length(y)
  (t1 && t2)
}


#' Check if a vector contains numeric data
#'
#' @param x A vector to check if it contains numeric data
#' @param y A vector to check if it contains numeric data
#'
#' @return Logical. TRUE if the vector contains mostly numeric data, FALSE otherwise.
#'
#' @details This function checks if a vector contains mostly numeric data.
#' It first replaces any "NA" values with actual NA values.
#' It then removes any empty values and converts the vector to numeric.
#' If more than 50% of the resulting numeric vector is not NA, it returns TRUE.
#'
#' @examples
#' \dontrun{
#' x <- c("1", "2", "foo", "3")
#' isanumber(x) # FALSE
#'
#' x <- c("1", "2", NA, "3")
#' isanumber(x) # TRUE
#' }
#' @export
isanumber <- function(x) {
  x <- sub("NA", NA, x)
  x[which(x == "")] <- NA
  suppressWarnings(nx <- as.numeric(x[!is.na(x)]))
  (length(nx) > 0 && mean(!is.na(nx)) > 0.5)
}


#' @title Expand annotation matrix
#'
#' @description Expands a phenotype annotation matrix into dummy variables.
#'
#' @param A Data frame containing the annotation variables.
#'
#' @details This function takes an annotation data frame and expands any categorical variables into
#' dummy variables using model.matrix.
#'
#' For each column, it determines if the variable is numeric or categorical.
#' Numeric variables are ranked. Categorical variables are expanded into dummy variables.
#'
#' @return An expanded annotation matrix with dummy variables.
#'
#' @export
expandAnnotationMatrix <- function(A) {
  expandPhenoMatrix(A)
}


#' @title Expand phenotype matrix
#'
#' @description Expands a phenotype data matrix into dummy variables while optionally
#' dropping the reference level.
#'
#' @param pheno Data frame containing the phenotype variables.
#' @param drop.ref Logical indicating whether to drop the reference level for each factor.
#'
#' @details This function takes a phenotype data matrix and expands any categorical variables into
#' dummy variables
#'
#' For each column, it determines if the variable is numeric or categorical. Numeric variables are
#' dichotomized. Categorical variables are expanded into dummy variables using \code{model.matrix}.
#'
#' @return An expanded phenotype matrix with dummy variables suitable for regression modeling.
#' @export
expandPhenoMatrix <- function(M, drop.ref = TRUE, keep.numeric = FALSE, check = TRUE) {
  ## get expanded annotation matrix
  a1 <- tidy.dataframe(M)
  nlevel <- apply(a1, 2, function(x) length(setdiff(unique(x), NA)))
  nterms <- colSums(!is.na(a1))
  nratio <- nlevel / nterms
  if (inherits(a1, "data.frame")) {
    a1.typed <- utils::type.convert(a1, as.is = TRUE)
    y.class <- sapply(a1.typed, function(a) class(a)[1])
  } else {
    ## matrix??
    a1.typed <- utils::type.convert(a1, as.is = TRUE)
    y.class <- apply(a1.typed, 2, function(a) class(a)[1])
  }

  ## these integers are probably factors... (mostly...)
  is.fac <- rep(FALSE, ncol(a1))
  is.int <- (y.class == "integer")
  ii <- which(is.int)
  is.fac[ii] <- apply(a1[, ii, drop = FALSE], 2, function(x) {
    nlev <- length(unique(x[!is.na(x)]))
    max(x, na.rm = TRUE) %in% c(nlev, nlev - 1) ## boolean
  })
  is.fac2 <- (y.class == "integer" & nlevel <= 3 & nratio < 0.66)
  y.class[is.fac | is.fac2] <- "character"

  ## select allowed columns: numeric or with "sensible" levels
  y.isnum <- (y.class %in% c("numeric", "integer"))
  kk <- which(y.isnum | (!y.isnum & nlevel > 1 & nratio < 0.66))
  if (length(kk) == 0) {
    return(NULL)
  }
  a1 <- a1[, kk, drop = FALSE]
  a1.isnum <- y.isnum[kk]

  i <- 1
  m1 <- list()
  for (i in 1:ncol(a1)) {
    if (a1.isnum[i]) {
      suppressWarnings(x <- as.numeric(a1[, i]))
      if (keep.numeric) {
        m0 <- matrix(x, ncol = 1)
        colnames(m0) <- "#"
      } else {
        if (drop.ref) {
          m0 <- matrix((x > stats::median(x, na.rm = TRUE)), ncol = 1)
          colnames(m0) <- "high"
        } else {
          mx <- stats::median(x, na.rm = TRUE)
          m0 <- matrix(cbind(x <= mx, x > mx), ncol = 2)
          colnames(m0) <- c("low", "high")
        }
      }
    } else if (drop.ref && nlevel[i] == 2) {
      x <- as.character(a1[, i])
      x1 <- Matrix::tail(sort(x), 1)
      m0 <- matrix(x == x1, ncol = 1)
      colnames(m0) <- x1
    } else {
      x <- as.character(a1[, i])
      x[is.na(x) | x == "NA" | x == " "] <- "_"
      m0 <- stats::model.matrix(~ 0 + x)
      colnames(m0) <- sub("^x", "", colnames(m0))
    }
    rownames(m0) <- rownames(a1)
    ## remove "_"
    if ("_" %in% colnames(m0)) {
      m0 <- m0[, -which(colnames(m0) == "_")]
    }
    m1[[i]] <- m0
  }

  ## create level names
  names(m1) <- colnames(a1)
  for (i in 1:length(m1)) {
    colnames(m1[[i]]) <- paste0(names(m1)[i], "=", colnames(m1[[i]]))
  }
  m1 <- do.call(cbind, m1)
  colnames(m1) <- sub("=#", "", colnames(m1))
  rownames(m1) <- rownames(M)
  return(m1)
}


#' @title P-value for Pearson's Correlation Coefficient
#'
#' @description This function calculates the p-value for Pearson's correlation coefficient.
#'
#' @param x A numeric value representing Pearson's correlation coefficient.
#' @param n A numeric value representing the sample size.
#'
#' @details The p-value is calculated using the standard normal distribution,
#' where the test statistic is calculated as `z = x / sqrt((1 - x^2) / (n - 2))`.
#' The p-value represents the probability of observing a correlation coefficient as
#' extreme or more extreme than `x` given that the null hypothesis (no correlation) is true.
#'
#' @return A numeric value representing the p-value for Pearson's correlation coefficient.
#'
#' @examples
#' \dontrun{
#' # example code
#' cor.pvalue(0.8, 100)
#' }
#' @export
cor.pvalue <- function(x, n) 2 * stats::pnorm(-abs(x / ((1 - x**2) / (n - 2))**0.5))


#' @title Get gene sets from playbase data
#'
#' @description Retrieves gene sets from the playbase data package matching a pattern.
#' Allows optionally including additional custom gene sets.
#'
#' @param sets Character vector of gene set names to retrieve. Default NULL retrieves all sets.
#' @param pattern Pattern to match gene set names. Default NULL matches all sets.
#'
#' @details This function extracts gene sets from the playbase data package.
#' It returns gene sets matching the provided pattern.
#' If a sets parameter is provided, only those specific sets are returned.
#' By default it returns all gene sets in playbase.
#'
#' @return A named list containing the gene sets matching the criteria.
#'
#' @export
getGSETS_playbase <- function(gsets = NULL, pattern = NULL) {
  if (is.null(gsets)) gsets <- names(playdata::iGSETS)
  if (!is.null(pattern)) {
    gsets <- grep(pattern, gsets, value = TRUE)
  }
  gsets <- intersect(gsets, names(playdata::iGSETS))
  lapply(playdata::iGSETS[gsets], function(idx) playdata::GSET_GENES[idx])
}

#' Normalize Matrix by Row
#'
#' Normalizes a matrix by dividing each row by the sum of its elements.
#'
#' @param G The matrix to be normalized.
#'
#' @return The normalized matrix.
#'
#' @export
normalize_rows <- function(G) {
  # efficient normalization using linear algebra
  row_sums <- Matrix::rowSums(G)
  D <- Matrix::Diagonal(x = 1 / row_sums)
  G_scaled <- D %*% G
  rownames(G_scaled) <- rownames(G)
  colnames(G_scaled) <- colnames(G)
  return(G_scaled)
}

#' @export
normalize_cols <- function(G) {
  # efficient normalization using linear algebra
  col_sums <- Matrix::colSums(G)
  D <- Matrix::Diagonal(x = 1 / col_sums)
  G_scaled <- G %*% D
  rownames(G_scaled) <- rownames(G)
  colnames(G_scaled) <- colnames(G)
  return(G_scaled)
}

#' @export
make_unique <- function(s) {
  has.dup <- sum(duplicated(s)) > 0
  if(!has.dup) return(s)
  n = 1
  while(has.dup) {
    jj <- which(duplicated(s))
    s[jj] <- paste0(sub("[.][1-9]*","",s[jj]),'.',n)
    has.dup <- sum(duplicated(s)) > 0
    n <- n + 1
  }
  s
}



## =====================================================================================
## =========================== END OF FILE =============================================
## =====================================================================================
