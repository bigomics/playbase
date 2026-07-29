##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##


#' @title Bar and bee swarm plot with significance stars
#'
#' @description
#' Creates a bar plot with overlaid bee swarm plot and significance stars.
#'
#' @param x Factor for x-axis groups.
#' @param y Numeric vector of values to plot.
#' @param first Subset of groups to plot first.
#' @param width Width of bars.
#' @param bar Logical to draw bar plot.
#' @param bee Logical to draw bee swarm plot.
#' @param sig.stars Logical to draw significance stars.
#' @param ymax Maximum y-axis value.
#' @param bee.cex Size of points in bee swarm plot.
#' @param max.stars Maximum number of stars to draw.
#' @param srt Rotation angle of x-axis labels.
#' @param xoff X-axis offset.
#' @param names.cex Size of x-axis labels.
#' @param names Logical to draw x-axis labels.
#' @param max.points Maximum number of points to draw.
#' @param col Color of points.
#' @param ... Additional graphics parameters.
#'
#' @details
#' This function creates a bar plot of \code{y} values grouped by \code{x}, with optional
#' bee swarm plot overlayed on top representing the distribution of values.
#' Significance stars highlighting differences between groups can also be added.
#'
#' @return None. Produces a bar/bee swarm plot.
#'
#' @export
gx.b3plot <- function(x, y, first = NULL,
                      width = 1, bar = TRUE, bee = TRUE, sig.stars = FALSE,
                      ymax = NULL, bee.cex = 0.3, max.stars = 5, srt = NULL, xoff = 0,
                      names.cex = 1, names = TRUE, max.points = 100, col = "grey80",
                      ...) {
  stats.segments <- function(x, y, xoffset = 0, lwd = 2) {
    bx <- graphics::boxplot(y ~ x, plot = FALSE)
    nx <- length(bx$n)
    x0 <- xoffset + (1:nx)
    #
    graphics::segments(x0 - 0.1, bx$conf[1, ], x0 + 0.1, bx$conf[1, ], lwd = lwd)
    graphics::segments(x0 - 0.1, bx$conf[2, ], x0 + 0.1, bx$conf[2, ], lwd = lwd)
    graphics::segments(x0, bx$conf[1, ], x0, bx$conf[2, ], lwd = lwd * 0.5)
  }
  ylevel <- levels(y)
  y <- as.character(y)
  if (any(is.na(y))) {
    y[is.na(y)] <- "NA"
    ylevel <- c(ylevel, "NA")
  }
  y <- factor(y, levels = ylevel, exclude = NULL)
  if (!is.null(first)) y <- stats::relevel(y, ref = first)
  mx <- tapply(x, y, stats::median, na.rm = TRUE)

  sig <- yc <- NULL
  if (sig.stars) {
    y.levels <- unique(y)
    yc <- utils::combn(y.levels, 2)
    pv <- rep(NA, ncol(yc))
    i <- 1
    for (i in 1:ncol(yc)) {
      grp <- yc[, i]
      pv[i] <- stats::t.test(x[which(y == grp[1])], x[which(y == grp[2])])$p.value
    }
    pv
    sig <- c("", "*", "**", "***")[1 + 1 * (pv < 0.05) + 1 * (pv < 0.01) + 1 * (pv < 0.001)]
    sig
    nthree <- sum(sig == "***")
    jj <- which(sig != "")
    jj <- jj[order(pv[jj])] ## only top 4 ??
    jj <- Matrix::head(jj, max(max.stars, nthree)) ## only top 5 ??
    j <- 1
    yc <- apply(yc, 2, as.integer)
    yc
    dd <- abs(as.vector(diff(yc)))
    jj <- jj[order(dd[jj])]
    sig <- sig[jj]
    yc <- yc[, jj, drop = FALSE]
  }

  dx <- max(x, na.rm = TRUE) * 0.11
  ylim <- c(xoff, max(x, na.rm = TRUE) * 1.3)
  if (!is.null(ymax)) ylim <- c(xoff, ymax)
  if (min(x, na.rm = TRUE) < 0) ylim <- c(1.3 * min(c(x, xoff)), max(x, na.rm = TRUE) * 1.3)
  if (sig.stars) {
    if (ncol(yc) > 8) dx <- dx / 5
    ylim[2] <- ylim[2] * 1.05 + (2 + NCOL(yc)) * dx
  }

  if (length(col) == 1) {
    col <- rep(col, length(mx))
    names(col) <- names(mx)
  } else {
    col <- col[match(names(mx), names(col))]
    col[is.na(col)] <- "grey90"
  }

  bx <- graphics::barplot(mx,
    width = 0.6666, space = 0.5, ylim = ylim, offset = xoff,
    names.arg = NA, col = col, ...
  )
  if (is.null(srt)) {
    nnchar <- sum(sapply(unique(y), nchar))
    srt <- ifelse(nnchar > 24, 30, 0)
  }
  pos <- ifelse(srt == 0, 1, 2)

  n <- length(unique(y))
  if (names == TRUE) {
    y0 <- min(ylim, na.rm = TRUE) - diff(ylim) * 0.05
    graphics::text(bx[, 1], y0, names(mx),
      cex = names.cex,
      srt = srt, adj = ifelse(srt == 0, 0.5, 0.965), xpd = TRUE,
      pos = pos, offset = 0
    )
  }
  if (bee) {
    jj <- 1:length(x)
    if (max.points > 0 && length(jj) > max.points) {
      jj <- unlist(tapply(jj, y, function(i) Matrix::head(sample(i), max.points)))
    }
    ## !!!!!!!!! NEED CHECK!! can be very slow if jj is large !!!!!!!!!!!
    beeswarm::beeswarm(x[jj] ~ y[jj], add = TRUE, at = 1:n - 0.33, pch = 19, cex = bee.cex, col = "grey20")
  }
  if (bar) stats.segments(y, x, xoffset = -0.333, lwd = 1.4)

  if (sig.stars) {
    i <- 1
    for (i in 1:NCOL(yc)) {
      grp <- yc[, i]
      xmax <- max(x, na.rm = TRUE) * 1.05 + dx * i
      j1 <- grp[1] - 0.4
      j2 <- grp[2] - 0.4
      graphics::segments(j1, xmax, j2, xmax, lwd = 0.5)
      if (ncol(yc) <= 8) {
        graphics::text((j1 + j2) / 2, xmax, labels = sig[i], pos = 1, offset = -0.33, adj = 0, cex = 1.4)
      }
    }
  }
}

#' @title Histogram of gene expression values
#'
#' @description
#' Generate a histogram of gene expression values.
#'
#' @param gx Gene expression matrix with genes in rows and samples in columns.
#' @param main Title for the plot.
#' @param ylim Limits for the y-axis.
#'
#' @details
#' This function generates a histogram of the gene expression values in \code{gx}.
#' It first creates a histogram of all values using \code{\link[graphics]{hist}}.
#' It then overlays density histograms for each sample, with colors corresponding to column numbers.
#'
#' @return
#' A histogram plot is generated, no value is returned.
#'
#' @examples
#' \dontrun{
#' gx <- matrix(rnorm(100 * 10), 100, 10)
#' gx.hist(gx)
#' }
#' @export
gx.hist <- function(gx, lines.col=NULL, ...) {
  h0 <- graphics::hist(
    as.vector(gx),
    col = "grey90", border = "grey80", freq = FALSE,
    ...
  )
  i <- 1
  if(!is.null(lines.col) && length(lines.col)==1) {
    lines.col <- rep(lines.col,ncol(gx))
  }

  for (i in 1:ncol(gx)) {
    h1 <- graphics::hist(gx[, i], breaks = h0$breaks, plot = FALSE)
    col1 <- ifelse(is.null(lines.col), i+1, lines.col[i])
    graphics::lines(h0$mids, h1$density, col = col1)
  }
}
