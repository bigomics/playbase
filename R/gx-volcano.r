##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


#' @title Volcano Plot
#'
#' @description
#' Generate a volcano plot to visualize significance versus magnitude of change.
#'
#' @param x Numeric vector of effect sizes (fold changes, differences, etc).
#' @param pv Numeric vector of p-values.
#' @param gene Character vector of gene IDs corresponding to x and pv.
#' @param ma_plot Logical, if TRUE creates an MA-plot instead.
#' @param ma Matrix of log-fold changes and mean average values for MA-plot.
#' @param p.sig Significance threshold p-value.
#' @param lfc Log2 fold change threshold.
#' @param render Plot rendering method (base, ggplot, scatterD3).
#' @param n Maximum number of points to plot.
#' @param highlight Vector of genes to highlight.
#' @param main Title for the plot.
#' @param cex Point size.
#' @param lab.cex Label size for highlighted genes.
#' @param nlab Number of top genes to label.
#' @param xlim Limits of the x-axis.
#' @param ylim Limits of the y-axis.
#' @param use.fdr Use FDR instead of p-values.
#' @param xlab Label for x-axis.
#' @param ylab Label for y-axis.
#' @param use.rpkm Use RPKM/FPKM values instead of raw counts.
#' @param maxchar Maximum characters for gene labels.
#' @param hi.col Color for highlighted points.
#' @param cex.main Size of main title.
#' @param axes Display axes?
#' @param cex.axis Size of axis labels.
#'
#' @details This function takes a vector of p-values and effect sizes (fold changes, differences, etc)
#' and generates a volcano plot to visualize significance versus magnitude of change. Points exceeding
#' both p-value and effect size thresholds are highlighted. The top most significant genes can be
#' labeled. Both raw and FDR corrected p-values can be used.
#'
#' @return A volcano plot is generated, no value is returned.
#'
#' @export
gx.volcanoPlot.XY <- function(x, pv, gene, ma_plot = FALSE, ma = NULL, p.sig = 0.05, lfc = 1,
                              render = "scatterD3", n = 1000, highlight = NULL, main = "",
                              cex = 1, lab.cex = 1, nlab = 10, xlim = NULL, ylim = NULL, use.fdr = FALSE,
                              xlab = "effect size", ylab = "significance (-log10p)",
                              use.rpkm = FALSE, maxchar = 40, hi.col = "#1e60bb",
                              cex.main = 1.2, axes = TRUE, cex.axis = 1) {
  jj <- which(!is.na(x) & !is.na(pv) & !is.na(gene))
  x <- x[jj]
  pv <- pv[jj]
  gene <- gene[jj]

  if (n > 0) {
    jj <- unique(c(1:100, utils::head(sample(1:length(x), replace = TRUE), n - 100)))
    if (!is.null(highlight)) {
      hgene <- intersect(highlight, gene)
      jj <- unique(c(match(hgene, gene), jj))
    }
  }

  pmin <- 1e-99
  if (!is.null(ylim)) pmin <- 10^(-max(ylim, na.rm = TRUE))
  if (is.null(ylim)) ylim <- c(0, max(-log10(pmin + pv)))
  y <- pmin(-log10(pv), max(ylim, na.rm = TRUE))
  nsig <- c("down" = sum(x <= -lfc & pv <= p.sig), "up" = sum(x >= lfc & pv <= p.sig))

  ## highlight significant
  impt <- function(g) {
    j <- match(g, gene)
    x1 <- scale(x, center = FALSE)[j]
    y1 <- scale(y, center = FALSE)[j]
    x <- sign(x1) * (x1**2 + 0.25 * y1**2)
    names(x) <- g
    x
  }

  ## highlight with color
  lab <- rep(NA, length(x))
  ii <- c()
  if (!is.null(highlight)) {
    ii <- match(intersect(highlight, gene), gene)
    lab[ii] <- gene[ii]
  } else {
    ii <- which(abs(x) >= lfc & pv <= p.sig)
  }
  klr <- rep(1, length(gene))
  if (length(ii) > 0) klr[ii] <- 2
  klr2 <- c("#BBBBBB", hi.col)[klr]


  ## add labels for some
  j0 <- Matrix::head(ii[order(impt(gene[ii]))], nlab)
  j1 <- Matrix::head(ii[order(-impt(gene[ii]))], nlab)
  jj <- unique(c(j0, j1))
  lab[jj] <- gene[jj]

  cex.wt <- 1
  if (is.null(xlim)) xlim <- c(-1.1, 1.1) * max(abs(x), na.rm = TRUE)
  dy <- 0.04 * (max(y, na.rm = TRUE) - min(y, na.rm = TRUE))
  if (is.null(ylim)) ylim <- c(min(y, na.rm = TRUE), max(y, na.rm = TRUE) + dy)
  if (is.infinite(ylim[2])) ylim[2] <- max(y, na.rm = TRUE) + dy
  lab <- substring(lab, 1, maxchar)

  if (ma_plot && !is.null(ma)) {
    y <- x
    x <- ma
    ylim <- 1.1 * xlim
    xlim <- NULL
  }

  gene.txt <- substring(gene, 1, 30) ## shortened for labels
  plt <- NULL

  if (render == "scatterD3") {
    tooltip_text <- paste(
      gene, "<br>x=", round(x, digits = 3),
      "<br>p=", round(pv, digits = 4)
    )
    jj <- order(klr)
    labjj <- NULL
    if (lab.cex > 0) labjj <- lab[jj]
    plt <- scatterD3::scatterD3(
      x = x[jj], y = y[jj],
      point_size = cex * 10,
      point_opacity = 0.66,
      xlab = xlab, ylab = ylab, col_var = klr[jj], legend_width = 0,
      lab = labjj, labels_size = lab.cex * 10,
      tooltip_text = tooltip_text[jj], caption = main,
      tooltips = TRUE, colors = c("2" = hi.col, "1" = "#BBBBBB")
    )
  } else if (render == "plotly") {
    gene2 <- paste0("  ", gene.txt, "  ")
    ann <- data.frame(gene = gene2, x = x, y = y)[jj, , drop = FALSE]
    ann.left <- ann[which(ann$x < 0), ]
    ann.right <- ann[which(ann$x >= 0), ]
    plt <- plotly::plot_ly(
      x = x, y = y, type = "scattergl", mode = "markers",
      marker = list(size = 5 * cex, color = "#BBBBBB"), hoverinfo = "text",
      text = tt
    ) %>%
      plotly::layout(
        xaxis = list(title = xlab, range = xlim),
        yaxis = list(title = ylab, range = ylim)
      )
    if (lab.cex > 0 && nrow(ann.left) > 0) {
      plt <- plt %>%
        plotly::add_annotations(
          x = ann.left$x, y = ann.left$y, text = ann.left$gene,
          xref = "x", yref = "y", xanchor = "right", showarrow = FALSE,
          font = list(size = 10 * lab.cex, color = hi.col)
        )
    }
    if (lab.cex > 0 && nrow(ann.right) > 0) {
      plt <- plt %>%
        plotly::add_annotations(
          x = ann.right$x, y = ann.right$y, text = ann.right$gene,
          xref = "x", yref = "y", xanchor = "left", showarrow = FALSE,
          font = list(size = 10 * lab.cex, color = hi.col)
        )
    }
    ## plt
  } else {
    if (!is.null(highlight)) cex.wt <- 0.5 * (1 + 1 * (gene %in% highlight))
    plot(
      x = x, y = y, pch = 19, cex = 0.4 * cex * cex.wt, xlim = xlim, ylim = ylim,
      col = klr2, xlab = xlab, ylab = ylab, main = main, cex.main = cex.main,
      axes = axes, cex.axis = cex.axis, cex.lab = cex.axis
    )
    psig <- -log10(p.sig)
    if (ma_plot) {
      graphics::abline(h = c(-1, 1) * lfc, v = 0, lty = 3, col = "grey50", lwd = 0.5)
      graphics::legend("topright",
        legend = paste(nsig["up"], "UP"), bty = "n",
        cex = 1, text.col = "grey50"
      )
      graphics::legend("bottomright",
        legend = paste(nsig["down"], "DOWN"), bty = "n",
        cex = 1, text.col = "grey50"
      )
    } else {
      graphics::abline(v = c(-1, 1) * lfc, h = c(-1, 1) * psig, lty = 3, col = "grey50", lwd = 0.5)
      graphics::legend("bottomright",
        legend = paste(nsig["up"], "UP"), bty = "n",
        cex = 1, text.col = "grey50"
      )
      graphics::legend("bottomleft",
        legend = paste(nsig["down"], "DOWN"), bty = "n",
        cex = 1, text.col = "grey50"
      )
    }
    if (length(jj) > 0 && lab.cex > 0) {
      graphics::text(
        x = x[jj], y = y[jj], labels = gene.txt[jj],
        pos = 3, cex = 0.65 * lab.cex, offset = 0.3, col = hi.col
      )
    }
  }
  plt
}

#' @title Volcano plot from LIMMA results
#'
#' @description Generate a volcano plot from a LIMMA differential
#' expression analysis result table.
#'
#' @param tab Data frame containing statistics from LIMMA analysis, with
#' columns for log2 fold change, p-values, and adjusted p-values.
#' @param render Plotting method (either base R graphics or scatterD3).
#' @param n Maximum number of points to plot.
#' @param highlight Vector of genes to highlight.
#' @param p.sig P-value cutoff for significance.
#' @param cex Point size.
#' @param lab.cex Label size for highlighted genes.
#' @param nlab Number of top genes to label.
#' @param xlim Limits of the x-axis.
#' @param ylim Limits of the y-axis.
#' @param use.fdr Use FDR adjusted p-values.
#' @param use.rpkm Use RPKM/TPM values for x-axis.
#' @param ma.plot Overlay MA plot lines?
#' @param cex.main Title size.
#' @param main Plot title.
#' @param cex.axis Axis label size.
#' @param axes Display axes?
#'
#' @details This function takes a LIMMA result table and generates a volcano plot
#' to visualize significance versus log2 fold change. Points exceeding the specified
#' p-value and fold change thresholds are highlighted. The top most significant
#' genes can be labeled. Both raw and FDR corrected p-values can be used.
#'
#' @return A volcano plot is generated, no value is returned.
#' @examples
#' \dontrun{
#' n <- 100
#' tab <- data.frame(
#'   logFC = rnorm(n),
#'   P.Value = runif(n),
#'   adj.P.Val = p.adjust(runif(n)),
#'   AveExpr = runif(n)
#' )
#' z <- gx.volcanoPlot.LIMMA(tab)
#' z
#' }
#' @export
#'
gx.volcanoPlot.LIMMA <- function(tab, render = "scatterD3", n = 1000, highlight = NULL, p.sig = 0.05,
                                 cex = 1, lab.cex = 1, nlab = 15, xlim = NULL, ylim = NULL, use.fdr = FALSE,
                                 use.rpkm = FALSE, ma.plot = FALSE, cex.main = 1.2,
                                 main = "", cex.axis = 1, axes = TRUE) {
  tab <- tab[order(tab$P.Value), ]
  gene <- as.character(tab[, grep("^gene$|^external_gene_name$", colnames(tab))])
  if (n > 0) {
    jj <- unique(c(1:100, utils::head(sample(1:nrow(tab), replace = TRUE), n - 100)))
    if (!is.null(highlight)) {
      hgene <- intersect(highlight, gene)
      jj <- unique(c(match(hgene, gene), jj))
    }
    tab <- tab[jj, ]
    gene <- gene[jj]
  }

  if (ma.plot) {
    if (use.rpkm && "AveRPKM" %in% colnames(tab)) {
      x <- tab$AveRPKM
      xlab <- "average RPKM  (log2)"
    } else {
      x <- tab$AveExpr
      xlab <- "average CPM  (log2)"
    }
    y <- tab$logFC
    ylab <- "fold change  (logFC)"
  } else {
    x <- tab$logFC
    if (use.fdr) {
      pmin <- 0.9 * min(tab$adj.P.Val[which(tab$adj.P.Val > 0)], na.rm = TRUE)
      y <- -log10(pmax(tab$adj.P.Val, pmin))
    } else {
      pmin <- 0.9 * min(tab$P.Value[which(tab$P.Value > 0)], na.rm = TRUE)
      y <- -log10(pmax(tab$P.Value, pmin))
    }
    xlab <- "effect size  (logFC)"
    ylab <- "significance  (-log10p)"
    if (use.fdr) ylab <- "significance  (-log10q)"
  }

  ## highlight significant
  impt <- function(g) {
    j <- match(g, gene)
    x1 <- tab$logFC[j]
    y1 <- -log10(tab$P.Value)[j]
    if (use.fdr) y1 <- -log10(tab$adj.P.Val)[j]
    x <- sign(x1) * (x1**2 + 0.3 * y1**2)
    names(x) <- g
    x
  }
  if (!is.null(highlight)) {
    gene1 <- intersect(highlight, gene)
    ii <- match(gene1, gene)
  } else {
    x1 <- tab$logFC
    y1 <- -log10(tab$P.Value)
    ii <- which(tab$P.Value < p.sig)
    if (use.fdr) {
      y1 <- -log10(tab$adj.P.Val)
      ii <- which(tab$adj.P.Val < p.sig)
    }
  }
  j0 <- Matrix::head(ii[order(impt(gene[ii]))], nlab)
  j1 <- Matrix::head(ii[order(-impt(gene[ii]))], nlab)
  jj <- unique(c(j0, j1))

  lab <- rep(NA, nrow(tab))
  lab[jj] <- gene[jj]
  klr <- rep(1, length(gene))
  klr[ii] <- 2
  klr2 <- c("#BBBBBB", "#1e60bb")[klr]

  if (is.null(xlim) && !ma.plot) {
    xlim <- c(-1.1, 1.1) * max(abs(x), na.rm = TRUE)
  }
  if (is.null(ylim) && ma.plot) {
    ylim <- c(-1.1, 1.1) * max(abs(y), na.rm = TRUE)
  }

  ## truncate outlier y values
  if (!is.null(ylim)) y <- pmin(y, max(ylim, na.rm = TRUE))

  tt <- paste(
    "gene=", gene, "<br>logFC=", round(tab$logFC, digits = 4),
    "<br>AveExpr=", round(tab$AveExpr, digits = 4),
    "<br>P.Value=", round(tab$P.Value, digits = 4),
    "<br>adj.P.Val=", round(tab$adj.P.Val, digits = 4)
  )

  if (render == "scatterD3") {
    scatterD3::scatterD3(
      x = x, y = y, point_size = cex * 10, point_opacity = 0.66,
      xlab = xlab, ylab = ylab, col_var = klr, legend_width = 0,
      lab = lab, labels_size = lab.cex * 8.5, tooltip_text = tt,
      tooltips = TRUE, colors = c("2" = "#1e60bb", "1" = "#BBBBBB"),
      caption = main
    )
  } else if (render == "plotly") {
    gene2 <- paste0("  ", gene, "  ")
    ann <- data.frame(gene = gene2, x = x, y = y)[jj, , drop = FALSE]
    ann.left <- ann[which(ann$x < 0), ]
    ann.right <- ann[which(ann$x >= 0), ]
    plotly::plot_ly(
      x = x, y = y, type = "scattergl", mode = "markers",
      marker = list(size = 5 * cex, color = "#BBBBBB"), hoverinfo = "text",
      text = tt
    ) %>%
      plotly::add_annotations(
        x = ann.left$x, y = ann.left$y, text = ann.left$gene,
        xref = "x", yref = "y", xanchor = "right", showarrow = FALSE,
        font = list(size = 10 * lab.cex, color = "#1e60bb")
      ) %>%
      plotly::add_annotations(
        x = ann.right$x, y = ann.right$y, text = ann.right$gene,
        xref = "x", yref = "y", xanchor = "left", showarrow = FALSE,
        font = list(size = 10 * lab.cex, color = "#1e60bb")
      ) %>%
      plotly::layout(
        xaxis = list(title = xlab, range = xlim),
        yaxis = list(title = ylab, range = ylim)
      )
  } else {
    if (!is.null(highlight)) cex.wt <- 0.5 * (1 + 1 * (gene %in% highlight))
    plot(
      x = x, y = y, pch = 19, cex = 0.4 * cex * cex.wt, xlim = xlim, ylim = ylim,
      col = klr2, xlab = xlab, ylab = ylab, main = main, cex.main = cex.main,
      cex.axis = cex.axis, axes = axes
    )
    psig <- -log10(p.sig)
    if (!ma.plot) {
      graphics::abline(v = c(-1, 1), h = c(-1, 1) * psig, lty = 3, col = "grey50", lwd = 0.5)
    }
    if (ma.plot) {
      graphics::abline(h = c(-1, 1), v = 0, lty = 3, col = "grey50", lwd = 0.5)
    }
    if (length(jj) > 0) {
      graphics::text(
        x = x[jj], y = y[jj], labels = gene[jj],
        pos = 3, cex = 0.65 * lab.cex, offset = 0.3, col = "#1e60bb"
      )
    }
  }
}
