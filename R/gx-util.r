##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


#' Impute missing values with row medians
#'
#' @param X Numeric matrix
#'
#' @return Matrix with missing values imputed
#'
#' @details
#' This function imputes missing values in a numeric matrix
#' by replacing them with the median value for each row.
#'
#' It first calculates the median of each row, ignoring NA values.
#' For rows where the median is NA, it takes the median of the row medians.
#' It then creates a copy of the input matrix and replaces all NA elements
#' with the corresponding row median value.
#'
#' @examples
#' \dontrun{
#' mat <- matrix(c(1:10, NA), nrow = 5)
#' imputed <- imputeMedian(mat)
#' }
#'
#' @export
imputeMedian <- function(X) {
  if (NCOL(X) == 1) {
    mx <- median(X, na.rm = TRUE)
  } else {
    mx <- apply(X, 1, median, na.rm = TRUE)
  }
  mx[is.na(mx)] <- median(mx, na.rm = TRUE)
  impX <- X
  impX[is.na(impX)] <- 0
  impX <- impX + is.na(X) * mx
  return(impX)
}

#' @title Calculate group means
#'
#' @description
#' Calculates the column means within groups defined by a grouping variable.
#'
#' @param X Numeric matrix with columns as samples.
#' @param y Grouping vector or factor.
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
averageByGroup <- function(X, y) {
  t(apply(X, 1, function(x) tapply(x, y, mean)))
}


#' @title Calculate geometric mean
#'
#' @description Calculates the geometric mean of a numeric vector.
#'
#' @param x A numeric vector.
#'
#' @return The geometric mean as a numeric value.
#'
#' @details The geometric mean is calculated by taking the nth root of the
#' product of n numbers, where n is the length of the input vector. It is
#' more appropriate than the arithmetic mean for positively skewed
#' distributions.
#'
#' @examples
#' \dontrun{
#' x <- c(1, 10, 100)
#' gmean(x)
#' }
#'
#' @export
gmean <- function(x) {
  ## geometric mean
  exp(mean(log(x + 1e-40)))
}




#' @title Histogram of gene expression values
#'
#' @description Generate a histogram of gene expression values.
#'
#' @param gx Gene expression matrix with genes in rows and samples in columns.
#' @param main Title for the plot.
#' @param ylim Limits for the y-axis.
#'
#' @details This function generates a histogram of the gene expression values in \code{gx}.
#' It first creates a histogram of all values using \code{\link[graphics]{hist}}.
#' It then overlays density histograms for each sample, with colors corresponding to column numbers.
#'
#' @return A histogram plot is generated, no value is returned.
#'
#' @examples
#' \dontrun{
#' gx <- matrix(rnorm(100 * 10), 100, 10)
#' gx.hist(gx)
#' }
gx.hist <- function(gx, main = "", ylim = NULL) {
  h0 <- hist(as.vector(gx),
    breaks = 120, main = main,
    col = "grey", freq = FALSE, ylim = ylim, xlab = "signal"
  )
  i <- 1
  for (i in 1:ncol(gx)) {
    h1 <- hist(gx[, i], breaks = h0$breaks, plot = FALSE)
    lines(h0$mids, h1$density, col = i + 1)
  }
}



#' @title Convert values to colors
#'
#' @description
#' Converts a numeric vector into a vector of colors based on the value.
#'
#' @param z A numeric vector of values to convert to colors.
#' @param zlim Limit values for the color scale. If missing, will use range of z.
#' @param col Colors to use. Passed to \code{\link[grDevices]{heat.colors}}.
#' @param breaks Break points for coloring. If missing, breaks are calculated automatically.
#'
#' @return A character vector of colors corresponding to the input values.
#'
#' @details
#' this function converts a vector of values("z") to a vector of color
#' levels. One must define the number of colors. The limits of the color
#' scale("zlim") or the break points for the color dplyr::changes("breaks") can
#' also be defined. when breaks and zlim are defined, breaks overrides zlim.
#' @seealso
#' Original post:
#'   \url{http://menugget.blogspot.ch/2011/09/converting-values-to-color-levels.html}
#'
#' @examples
#' \dontrun{
#' z <- rnorm(100)
#' cols <- val2col(z, zlim = c(-3, 3))
#' barplot(height = z, col = cols)
#' }
#' @export
val2col <- function(z, zlim, col = heat.colors(12), breaks) {
  if (!missing(breaks)) {
    if (length(breaks) != (length(col) + 1)) {
      stop("must have one more break than colour")
    }
  }
  if (missing(breaks) & !missing(zlim)) {
    zlim[2] <- zlim[2] + c(zlim[2] - zlim[1]) * (1E-3) # adds a bit to the range in both directions
    zlim[1] <- zlim[1] - c(zlim[2] - zlim[1]) * (1E-3)
    breaks <- seq(zlim[1], zlim[2], length.out = (length(col) + 1))
  }
  if (missing(breaks) & missing(zlim)) {
    zlim <- range(z, na.rm = TRUE)
    zlim[2] <- zlim[2] + c(zlim[2] - zlim[1]) * (1E-3) # adds a bit to the range in both directions
    zlim[1] <- zlim[1] - c(zlim[2] - zlim[1]) * (1E-3)
    breaks <- seq(zlim[1], zlim[2], length.out = (length(col) + 1))
  }
  CUT <- cut(z, breaks = breaks)
  colorlevels <- col[match(CUT, levels(CUT))] # assign colors to heights for each point
  return(colorlevels)
}


#' Convert gene symbols to official HUGO gene symbols
#'
#' @param genes Character vector of gene symbols to convert
#' @param remove.non.hugo Logical indicating whether to remove non-HUGO symbols. Default is TRUE.
#' @param silent Logical indicating whether to suppress messages about conversions. Default is FALSE.
#' @param take.only.first Logical indicating whether to take only first HUGO symbol for aliases. Default is FALSE.
#' @param split.char Character used to split multiple aliases. Default ";".
#' @param unknown Character string to use for unknown symbols. Default "unknown_gene".
#'
#' @return Character vector of official HUGO gene symbols.
#'
#' @details This function converts a character vector of gene symbols to official HUGO gene symbols.
#' It first removes any aliases by mapping to the primary HUGO symbol list.
#' For any remaining non-HUGO symbols, it attempts to find the official symbol by alias mapping.
#'
#' If take.only.first is TRUE, only the first HUGO symbol is taken when a gene maps to multiple aliases.
#' Multiple aliases are concatenated by split.char when take.only.first=FALSE.
#'
#' Non-HUGO symbols that cannot be converted are replaced with unknown by default.
#' Messages about conversions are printed unless silent=TRUE.
#'
#' @examples
#' \dontrun{
#' genes <- c("EGFR", "CDKN2A", "FOO", "BAR")
#' hugo_symbols <- symbol2hugo(genes)
#' }
#'
#' @export
symbol2hugo <- function(genes, remove.non.hugo = TRUE, silent = FALSE,
                        take.only.first = FALSE, split.char = ";", unknown = "unknown_gene") {
  HUGO.SYMBOLS <- unique(unlist(as.list(org.Hs.eg.db::org.Hs.egSYMBOL)))
  ss <- as.character(genes)
  ss <- gsub("Sep 0", "SEPT", ss) # typical XLS error
  ss[is.na(ss) | ss == ""] <- unknown
  ii <- which(!(ss %in% HUGO.SYMBOLS) & ss != unknown)
  length(ii)
  if (length(ii) == 0) {
    return(genes)
  }
  if (!silent) cat("trying to convert", length(ii), "aliases to HUGO\n")
  ss0 <- sapply(ss[ii], strsplit, split = split.char)
  ee0 <- lapply(ss0, function(s) unlist(mget(s, envir = org.Hs.eg.db::org.Hs.egALIAS2EG, ifnotfound = NA)))
  ee0 <- lapply(ee0, function(e) {
    e[is.na(e) | e == "" | is.nan(e)] <- unknown
    e
  })
  gg <- lapply(ee0, function(e) unlist(mget(e, envir = org.Hs.eg.db::org.Hs.egSYMBOL, ifnotfound = NA)))
  if (remove.non.hugo) {
    gg <- lapply(gg, intersect, HUGO.SYMBOLS)
  }
  gg.len <- lapply(gg, length)
  sum(gg.len > 1)
  if (sum(gg.len > 1) && !silent) {
    cat("warning:", sum(gg.len > 1), "entrezID have multiple symbols\n")
  }
  if (sum(gg.len > 1) && take.only.first) {
    gg <- sapply(gg, "[", 1)
  } else {
    gg <- sapply(gg, paste, collapse = split.char)
  }
  if (!silent) {
    cat("updating", length(gg), "deprecated symbols\n")
  }
  gg[is.na(gg) | gg == ""] <- unknown
  ss[ii] <- gg
  ss
}

