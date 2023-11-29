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
    mx <- stats::median(X, na.rm = TRUE)
  } else {
    mx <- apply(X, 1, stats::median, na.rm = TRUE)
  }
  mx[is.na(mx)] <- stats::median(mx, na.rm = TRUE)
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
  h0 <- graphics::hist(as.vector(gx),
    breaks = 120, main = main,
    col = "grey", freq = FALSE, ylim = ylim, xlab = "signal"
  )
  i <- 1
  for (i in 1:ncol(gx)) {
    h1 <- graphics::hist(gx[, i], breaks = h0$breaks, plot = FALSE)
    graphics::lines(h0$mids, h1$density, col = i + 1)
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