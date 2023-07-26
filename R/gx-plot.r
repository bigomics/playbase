##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


#' Barplot for grouped data
#'
#' @param x Numeric matrix with samples in columns
#' @param main Title for the plot
#' @param cex.main Size of main title
#' @param cex.names Size of sample names
#' @param cex.legend Size of legend text
#' @param srt Rotation angle of sample names
#' @param xlab X-axis label
#' @param ylab Y-axis label
#' @param group Vector or factor of group assignments
#' @param group.names Optional names for the groups
#' @param bar.names Optional names for each bar
#' @param voff Vertical offset between bars
#' @param legend Show legend for groups?
#'
#' @details
#' This function creates a barplot from a numeric matrix, with samples in columns.
#' Bars are grouped if a group vector is provided. Error bars show the standard error.
#' Many graphical parameters can be adjusted, including title, axis labels,
#' legend, and text sizes.
#'
#' @return None. Produces a barplot.
#'
#' @return
#' @export
#'
#' @examples
gx.barplot <- function(x, main = "", cex.main = 1.2, cex.names = 0.85,
                       cex.legend = 0.9, srt = 0, xlab = "", ylab = "",
                       group = NULL, group.names = NULL, bar.names = NULL,
                       voff = 2, legend = TRUE) {
  cpal <- rep(RColorBrewer::brewer.pal(12, "Paired"), 99)
  if (is.null(group)) group <- rep(1, ncol(x))
  group2 <- as.vector(sapply(group, rep, 2))
  col1 <- cpal[(group2 - 1) * 2 + rep(1:2, NCOL(x))]
  if (NCOL(x) == 1) col1 <- col1[2]
  ylim <- c(0, 1.15 * max(x, na.rm = TRUE))
  str.ht <- strheight("XXX", cex = cex.names)
  str.ht <- 0.04 * ylim[2]
  space <- c(0, 0.4)
  is.multi <- (NCOL(x) > 1)
  is.multi
  if (!is.multi) space <- c(0.2, 0)

  barplot(x,
    beside = TRUE, col = col1,
    main = main, cex.main = cex.main,
    xlab = xlab, ylab = ylab, ylim = ylim,
    las = 1, space = space,
    names.arg = rep("", length(x)),
    mgp = c(2, 0.9, 0)
  )

  if (is.multi) {
    sp0 <- nrow(x) + space[2]
    tx0 <- sp0 * (0:(ncol(x) - 1)) + 1
    tx1 <- tx0 + 0.5 * nrow(x) + space[1]
    tx <- as.vector(mapply(tx0, tx1, FUN = c))
    tx2 <- (tx0 + tx1) / 2
    if (is.null(bar.names)) {
      bar.names <- rep(rownames(x), NCOL(x))
      bar.names <- toupper(substring(bar.names, 1, 3))
    }
    mtext(bar.names, side = 1, line = 0, at = tx, adj = 0.5, cex = 0.5, srt = 0, xpd = TRUE)
    if (is.null(group.names)) {
      group.names <- colnames(x)
    }
    if (srt == 0) {
      mtext(group.names,
        side = 1, line = 1.5, at = tx2,
        adj = 0.5, cex = cex.names, srt = 0, xpd = TRUE
      )
    } else {
      text(tx2, -voff * str.ht, group.names,
        cex = cex.names,
        srt = srt, xpd = TRUE, pos = 2, offset = 0
      )
    }
  } else {
    text(1.2 * (1:length(x)) - 0.5, -voff * str.ht, names(x),
      cex = cex.names,
      srt = srt, xpd = TRUE, pos = 2, offset = 0 * str.ht
    )
  }

  if (legend) {
    nn <- names(x)
    if (is.multi) nn <- rownames(x)
    legend("topleft",
      legend = nn, fill = col1,
      xpd = TRUE, bty = "n", inset = c(-0.0, -0.04), cex = cex.legend,
      y.intersp = 0.8, x.intersp = 0.3
    )
  }
}


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
#'
#' @examples
gx.b3plot <- function(x, y, first = NULL,
                      width = 1, bar = TRUE, bee = TRUE, sig.stars = FALSE,
                      ymax = NULL, bee.cex = 0.3, max.stars = 5, srt = NULL, xoff = 0,
                      names.cex = 1, names = TRUE, max.points = 100, col = "grey80",
                      ...) {
  stats.segments <- function(x, y, xoffset = 0, lwd = 2) {
    bx <- boxplot(y ~ x, plot = FALSE)
    nx <- length(bx$n)
    x0 <- xoffset + (1:nx)
    #
    segments(x0 - 0.1, bx$conf[1, ], x0 + 0.1, bx$conf[1, ], lwd = lwd)
    segments(x0 - 0.1, bx$conf[2, ], x0 + 0.1, bx$conf[2, ], lwd = lwd)
    segments(x0, bx$conf[1, ], x0, bx$conf[2, ], lwd = lwd * 0.5)
  }
  ylevel <- levels(y)
  y <- as.character(y)
  if (any(is.na(y))) {
    y[is.na(y)] <- "NA"
    ylevel <- c(ylevel, "NA")
  }
  y <- factor(y, levels = ylevel, exclude = NULL)
  if (!is.null(first)) y <- relevel(y, ref = first)
  mx <- tapply(x, y, median, na.rm = TRUE)

  sig <- yc <- NULL
  if (sig.stars) {
    y.levels <- unique(y)
    yc <- combn(y.levels, 2)
    pv <- rep(NA, ncol(yc))
    i <- 1
    for (i in 1:ncol(yc)) {
      grp <- yc[, i]
      pv[i] <- t.test(x[which(y == grp[1])], x[which(y == grp[2])])$p.value
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
  ylim <- c(xoff, max(x) * 1.3)
  if (!is.null(ymax)) ylim <- c(xoff, ymax)
  if (min(x) < 0) ylim <- c(1.3 * min(c(x, xoff)), max(x) * 1.3)
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

  bx <- barplot(mx,
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
    y0 <- min(ylim) - diff(ylim) * 0.05
    text(bx[, 1], y0, names(mx),
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
      segments(j1, xmax, j2, xmax, lwd = 0.5)
      if (ncol(yc) <= 8) {
        text((j1 + j2) / 2, xmax, labels = sig[i], pos = 1, offset = -0.33, adj = 0, cex = 1.4)
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
