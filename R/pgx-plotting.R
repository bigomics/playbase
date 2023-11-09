##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

########################################################################
## Plotting functions
########################################################################


#' @title heatmapWithAnnot
#'
#' @param F Numeric data matrix to visualize as a heatmap
#' @param anno.type Type of annotation to show - "boxplot" or "barplot"
#' @param bar.height Height of the annotation barplot in cm
#' @param map.height Height of the heatmap in cm
#' @param row_fontsize Fontsize for row labels
#' @param column_fontsize Fontsize for column labels
#' @param inset Inset margins for annotation
#' @param mar Heatmap margins
#' @param legend Show legend for annotation colors
#' @param ... Other arguments passed to heatmap functions
#'
#' @return A heatmap grob object
#'
#' @description Creates a heatmap visualization with side annotation bars.
#'
#' @details This function generates a heatmap from the input data matrix F.
#' Side annotation bars can be added, either as boxplots or barplots.
#' The annotation height, heatmap height, font sizes, margins, and other
#' heatmap parameters can be customized.
#'
#' @export
heatmapWithAnnot <- function(F, anno.type = c("boxplot", "barplot"),
                             bar.height = NULL, map.height = NULL,
                             row_fontsize = 9, column_fontsize = 9,
                             inset = c(-0.025, -0.1),
                             mar = c(0, 0, 0, 0), legend = TRUE,
                             ...) {
  col1 <- 1:ncol(F)
  anno.type <- anno.type[1]
  if (!is.null(map.height)) map.height <- grid::unit(map.height, "cm")
  if (!is.null(bar.height)) bar.height <- grid::unit(bar.height, "cm")

  if (anno.type == "barplot") {
    ha <- ComplexHeatmap::HeatmapAnnotation(
      up = ComplexHeatmap::anno_barplot(pmax(F, 0.01), gp = gpar(fill = col1)),
      down = ComplexHeatmap::anno_barplot(pmin(F, -0.01), gp = gpar(fill = col1)),
      annotation_height = c(0.5, 0.5) * bar.height,
      annotation_label = c("", "")
    )
  } else {
    ha <- ComplexHeatmap::HeatmapAnnotation(
      up = ComplexHeatmap::anno_boxplot(
        t(F),
        gp = gpar(fill = "grey80"),
        outline = TRUE, size = grid::unit(1, "mm"),
        box_width = 0.75
      ),
      annotation_height = bar.height,
      annotation_label = c("")
    )
  }

  ht <- ComplexHeatmap::Heatmap(
    t(F),
    name = "logFC",
    top_annotation = ha,
    row_names_gp = grid::gpar(fontsize = row_fontsize),
    column_names_gp = grid::gpar(fontsize = column_fontsize),
    height = map.height,
    heatmap_legend_param = list(
      legend_height = grid::unit(8, "mm")
    ),
    ...
  )

  ComplexHeatmap::draw(ht,
    padding = grid::unit(mar, "cm"),
    heatmap_legend_side = "left"
  )

  if (anno.type == "boxplot") {
    ComplexHeatmap::decorate_annotation("up", {
      grid::grid.text("logFC",
        x = grid::unit(0, "npc") - grid::unit(7, "mm"),
        y = grid::unit(0.5, "npc") - grid::unit(0, "mm"),
        gp = grid::gpar(fontsize = 10),
        rot = 90, hjust = "center"
      )
    })
  }
  if (anno.type == "barplot") {
    ComplexHeatmap::decorate_annotation("up", {
      grid::grid.text("cumulative logFC",
        x = grid::unit(0, "npc") - grid::unit(7, "mm"),
        y = grid::unit(0, "npc") - grid::unit(0, "mm"),
        gp = grid::gpar(fontsize = 10),
        rot = 90, hjust = "center"
      )
    })
  }
  if (anno.type == "barplot" && legend) {
    graphics::legend("topright", colnames(F),
      fill = col1,
      cex = 0.63, y.intersp = 0.8,
      inset = inset, xpd = TRUE
    )
  }
}


#' @title Repel overlapping text labels
#'
#' @param x Numeric vector of x coordinates for text labels
#' @param y Numeric vector of y coordinates for text labels
#' @param words Character vector of words/labels to plot
#' @param cex Size multiplier for text labels
#' @param rotate90 Rotate labels by 90 degrees
#' @param xlim x-axis limits
#' @param ylim y-axis limits
#' @param tstep Translation step size
#' @param rstep Rotation step size
#' @param maxiter Maximum number of iterations
#' @param ... Additional graphics parameters to graphics::text()
#'
#' @return The adjusted coordinates, rotated labels, and plotting parameters
#'
#' @description Repels overlapping text labels by iteratively adjusting positions.
#'
#' @details This function takes a set of text labels defined by x,y coordinates and words.
#' It iteratively adjusts the positions to repel overlapping labels. Rotation by 90 degrees
#' and size scaling can also be applied.
#'
#' The algorithm translates or rotates labels in small steps to minimize overlap.
#' The step sizes and maximum number of iterations can be controlled.
#'
#' Useful for scatterplots, word clouds, and other graphics where label overlap is a problem.
#'
#' @export
repelwords <- function(x, y, words, cex = 1, rotate90 = FALSE,
                       xlim = c(-Inf, Inf), ylim = c(-Inf, Inf),
                       tstep = 0.1, rstep = 0.1, maxiter = 2000, ...) {
  ## From wordcloud::wordlayout
  tails <- "g|j|p|q|y"
  n <- length(words)
  sdx <- stats::sd(x, na.rm = TRUE)
  sdy <- stats::sd(y, na.rm = TRUE)
  if (sdx == 0) {
    sdx <- 1
  }
  if (sdy == 0) {
    sdy <- 1
  }
  if (length(cex) == 1) {
    cex <- rep(cex, n)
  }
  if (length(rotate90) == 1) {
    rotate90 <- rep(rotate90, n)
  }
  boxes <- list()

  theta <- 0
  for (i in 1:length(words)) {
    rotWord <- rotate90[i]
    r <- 0
    x1 <- xo <- x[i]
    y1 <- yo <- y[i]


    wid <- graphics::strwidth(words[i], cex = cex[i])
    ht <- graphics::strheight(words[i], cex = cex[i])
    if (grepl(tails, words[i])) {
      ht <- ht + ht * 0.2
    }
    if (rotWord) {
      tmp <- ht
      ht <- wid
      wid <- tmp
    }
    isOverlapped <- TRUE
    iter <- 1

    while (isOverlapped && iter < maxiter) {
      if (!wordcloud:::is_overlap(
        x1 - 0.5 * wid, y1 - 0.5 * ht, wid,
        ht, boxes
      ) &&
        x1 - 0.5 * wid > xlim[1] &&
        y1 - 0.5 * ht > ylim[1] &&
        x1 + 0.5 * wid < xlim[2] &&
        y1 + 0.5 * ht < ylim[2]) {
        boxes[[length(boxes) + 1]] <-
          c(x1 - 0.5 * wid, y1 - 0.5 * ht, wid, ht)
        isOverlapped <- FALSE
      } else {
        theta <- theta + tstep
        r <- r + rstep * tstep / (2 * pi)

        x1 <- xo + sdx * r * cos(theta)
        y1 <- yo + sdy * r * sin(theta)
      }
      iter <- iter + 1
    }
    if (isOverlapped && iter == maxiter) {
      message("[repelwords] WARNING maximum iterations reached: iter = ", iter)
      boxes[[length(boxes) + 1]] <-
        c(x1 - 0.5 * wid, y1 - 0.5 * ht, wid, ht)
    } else {

    }
  }
  result <- do.call(rbind, boxes)

  colnames(result) <- c("x", "y", "width", "ht")
  rownames(result) <- words
  result
}

## =================================================================================
## PGX level plotting API
## =================================================================================




#' @title Scatter plot for PGX object
#'
#' @param pgx PGX object with results
#' @param pheno Phenotype data.frame or vector
#' @param gene Gene name to color by expression
#' @param contrast Contrast name to color by correlation
#' @param level Data level to use ("gene", "geneset")
#' @param plotlib Plotting library ("base", "plotly")
#' @param pos Positions for points (if level="geneset")
#' @param ... Other arguments passed to plotting functions
#'
#' @return A scatter plot grob or plotly object
#'
#' @description Generates a scatter plot from a PGX object.
#'
#' @details Samples or genes can be plotted, colored by phenotype,
#' gene expression, or correlation with a contrast. For gene sets,
#' positions need to be provided. Interactive plotly plots or
#' static ggplot2 plots can be generated.
#'
#' Useful for exploring sample relationships and associations.
#'
#' @export
pgx.scatterPlot <- function(pgx, pheno = NULL, gene = NULL,
                            contrast = NULL, level = "gene",
                            plotlib = "base", pos = NULL, ...) {
  ## Scatter t-SNE plot samples (or genes) colored on phenotype,
  ## gene expression, geneset expresssion or (correlation with)
  ## contrast.
  ##
  if (level == "geneset" && is.null(pos)) {
    stop("FATAL:: geneset scatter needs position vector")
  }

  if (level == "geneset") {
    X <- pgx$gsetX
  } else {
    X <- pgx$X
  }

  plt <- NULL

  if (!is.null(gene)) {
    if (is.null(pos) && level == "gene") pos <- pgx$tsne2d
    var <- X[gene, rownames(pos)]
    title <- gene
    plt <- pgx.scatterPlotXY(
      pos, var,
      plotlib = plotlib, #
      xlab = colnames(pos)[1], ylab = colnames(pos)[2],
      type = "numeric", ...
    )
  }
  if (!is.null(pheno)) {
    if (is.null(pos) && level == "gene") pos <- pgx$tsne2d
    var <- pgx$samples[rownames(pos), pheno]
    title <- pheno
    plt <- pgx.scatterPlotXY(
      pos, var,
      plotlib = plotlib, #
      xlab = colnames(pos)[1], ylab = colnames(pos)[2],
      ...
    )
  }
  if (!is.null(contrast)) {
    if (is.null(pos) && "cluster.genes" %in% names(pgx)) {
      pos <- pgx$cluster.genes$pos[[1]][, 1:2]
    }
    if (is.null(pos)) {
      stop("must supply positions for gene contrasts")
    }
    if (level == "gene") {
      var <- pgx$gx.meta$meta[[contrast]]$meta.fx
      names(var) <- rownames(pgx$gx.meta$meta[[contrast]])
      var <- var[rownames(pos)]
      tooltip <- probe2symbol(rownames(counts), ngs$genes)
    }
    if (level == "geneset") {
      var <- pgx$gset.meta$meta[[contrast]]$meta.fx
      names(var) <- rownames(pgx$gset.meta$meta[[contrast]])
      var <- var[rownames(pos)]
      tooltip <- rownames(pos)
    }

    plt <- pgx.scatterPlotXY(
      pos, var,
      plotlib = plotlib, #
      xlab = colnames(pos)[1], ylab = colnames(pos)[2],
      tooltip = tooltip,
      ...
    )
  }

  if (plotlib == "base") {
    return(NULL)
  }
  plt
}


#' @title Plot a Scatterplot Matrix
#'
#' @description
#' Creates a scatterplot matrix visualization comparing multiple datasets.
#'
#' @param F Numeric matrix with samples in columns.
#' @param F2 Optional second numeric matrix to compare against F.
#' @param hilight Vector of sample names to highlight.
#' @param cex Point size in scatterplots.
#' @param cex.axis Axis text size.
#' @param cex.space Spacing between scatterplots.
#'
#' @details
#' This function generates a scatterplot matrix comparing all columns of
#' the input matrix \code{F} against each other. If a second matrix \code{F2} is
#' provided, columns of \code{F} are compared against matching columns in \code{F2}.
#'
#' Samples specified in \code{hilight} are colored differently to stand out.
#' The \code{cex}, \code{cex.axis}, and \code{cex.space} parameters control graphical
#' element sizes.
#'
#' @return
#' A scatterplot matrix is generated comparing datasets. No value is returned.
#'
#' @export
plot_SPLOM <- function(F, F2 = NULL, hilight = NULL, cex = 0.5, cex.axis = 1, cex.space = 0.2) {
  if (is.null(F2)) F2 <- F
  symm <- all(colnames(F) == colnames(F2))
  gg <- intersect(rownames(F), rownames(F2))
  F <- F[gg, , drop = FALSE]
  F2 <- F2[gg, , drop = FALSE]

  x0 <- range(as.vector(apply(F, 2, stats::quantile, probs = c(0.001, 0.999))))
  x1 <- range(as.vector(apply(F2, 2, stats::quantile, probs = c(0.001, 0.999))))
  x0 <- range(as.vector(F))
  x1 <- range(as.vector(F2))
  x0 <- x0 + c(-1, 1) * diff(x0) * 0.05
  x1 <- x1 + c(-1, 1) * diff(x1) * 0.05

  nr <- 3
  nc <- 4
  nr <- ncol(F2)
  nc <- ncol(F)
  graphics::par(
    mfrow = c(nr, nc), mar = c(1, 1, 1, 1) * cex.space, oma = c(4, 4, 0, 0),
    bty = "o", xpd = FALSE
  )
  i <- 1
  for (j in ncol(F2):1) {
    for (i in 1:ncol(F)) {
      plot(F[, i], F2[, j],
        pch = 20, col = "grey70",
        cex = cex, cex.axis = 0.95 * cex.axis, cex.lab = cex.axis,
        xaxt = ifelse(j == 1, "s", "n"),
        yaxt = ifelse(i == 1, "s", "n"),
        xlab = "", ylab = "",
        xlim = x0, ylim = x1
      )
      graphics::abline(v = 0, h = 0, lty = 3, lwd = 0.6)
      if (j == 1) graphics::mtext(colnames(F)[i], 1, line = 2.7, cex = 0.9 * cex.axis)
      if (i == 1) graphics::mtext(colnames(F2)[j], 2, line = 2.7, cex = 0.9 * cex.axis)

      if (!is.null(hilight) && length(hilight) > 0) {
        hilight <- intersect(hilight, rownames(F))
        ii <- match(hilight, rownames(F))
        graphics::points(F[ii, i], F2[ii, j], pch = 20, col = "red3")
        graphics::text(F[ii, i], F2[ii, j], labels = hilight, cex = 0.85, pos = 3, col = "black")
      }

      ## write correlation value
      rho <- stats::cor(F[, i], F2[, j], use = "pairwise")
      rr <- paste("r =", round(rho, digits = 3))
      graphics::legend("topleft", legend = rr, bty = "n", cex = 1)
    }
  }
}


#' @title Create a sankey diagram from a list of matrices
#'
#' @param matlist A list of matrices to compare
#' @param contrast Optional contrast name
#'
#' @return A plotly sankey diagram object
#'
#' @description
#' Generates a sankey diagram comparing matrices in a list
#'
#' @details
#' This function takes a list of matrices as input matlist.
#' It standardizes each matrix, calculates cross-tables between matrices,
#' and generates a sankey diagram showing the relationships.
#'
#' The contrast parameter can optionally be used to weight the edges.
#'
#' Useful for visualizing relationships and flows between data matrices.
#'
#' @export
pgx.SankeyFromMatrixList.PLOTLY <- function(matlist, contrast = NULL) {
  ## ------------------------------------------------------------------------
  ## Prepare matrices
  ## ------------------------------------------------------------------------
  X <- list()
  for (i in 1:length(matlist)) {
    X[[i]] <- matlist[[i]] - rowMeans(matlist[[i]])
    X[[i]] <- X[[i]] / apply(X[[i]], 1, stats::sd)
  }

  ## Counts cross-table between matrices
  M <- list()
  i <- 1
  for (i in 1:(length(X) - 1)) {
    mm <- pmax(X[[i]], 0) %*% t(pmax(X[[i + 1]], 0))
    mm <- mm**4
    mm <- mm / mean(mm)

    M[[i]] <- mm
  }

  ## Correlation
  R <- list()
  for (i in 1:(length(X) - 1)) {
    r1 <- stats::cor(t(X[[i]]), t(X[[i + 1]]))
    R[[i]] <- pmax(r1, 0)
  }

  ## Edge value (i.e. capacity) : rho * contrast/FC
  cty.mode <- 1
  F <- R
  if (!is.null(contrast)) {
    fc <- lapply(X, function(m) stats::cor(t(m), contrast)[, 1])
    i <- 1
    for (i in 1:length(R)) {
      if (cty.mode == 1) node.wt <- outer(pmax(fc[[i]], 0), pmax(fc[[i + 1]], 0))
      if (cty.mode == 3) node.wt <- abs(outer(fc[[i]], fc[[i + 1]]))
      if (cty.mode == 2) node.wt <- pmax(outer(fc[[i]], fc[[i + 1]]), 0)
      ww <- R[[i]] * node.wt
      ww <- ww / max(ww) ## normalize??
      F[[i]] <- ww
    }
  }
  fill0 <- FALSE
  fill0 <- !is.null(contrast)
  pgx.SankeyFromMRF.PLOTLY(M = M, R = R, F = F, fill = fill0, labels = NULL)
}


#' @title Create a sankey diagram from multiple matrices
#'
#' @param M A list of matrices
#' @param R A list of correlation matrices
#' @param F A list of flow matrices
#' @param fill Logical, whether to color edges by flow
#' @param labels Optional edge labels
#'
#' @return A plotly sankey diagram object
#'
#' @description
#' Generates a sankey diagram by combining multiple matrices.
#'
#' @details
#' This function takes a list of matrices \code{M}, correlation matrices \code{R},
#' and flow matrices \code{F}. It converts each matrix pair \code{M[[i]]} and
#' \code{R[[i]]} into a graph, merges the graphs, and creates a sankey diagram
#' to show relationships and flows between matrices.
#'
#' If \code{fill=TRUE}, edge colors are determined by the flow values \code{F[[i]]}.
#' Optional \code{labels} can be provided to label the edges.
#'
#' Useful for visualizing relationships and flows between multiple data matrices.
#'
#' @export
pgx.SankeyFromMRF.PLOTLY <- function(M, R, F, fill = TRUE, labels = NULL) {
  rho2graph <- function(A, min.rho = 0) {
    idx <- which(A > min.rho, arr.ind = TRUE)
    ee <- cbind(rownames(A)[idx[, 1]], colnames(A)[idx[, 2]])
    gr <- igraph::graph_from_edgelist(ee, directed = TRUE)
    gr
  }

  gr.list <- list()
  i <- 1
  for (i in 1:length(M)) {
    gr <- rho2graph(M[[i]], min.rho = 0)
    ee <- igraph::get.edgelist(gr)
    pct <- M[[i]] / sum(M[[i]], na.rm = TRUE) * 100
    igraph::E(gr)$count <- pct[ee]
    igraph::E(gr)$rho <- R[[i]][ee]
    igraph::E(gr)$weight <- F[[i]][ee]
    gr.list[[i]] <- gr
  }

  ## merge graphs
  gr <- gr.list[[1]]
  if (length(gr.list) > 1) {
    i <- 2
    for (i in 2:length(gr.list)) {
      gr <- igraph::union(gr, gr.list[[i]])
      igraph::E(gr)$count <- rowSums(cbind(igraph::E(gr)$count_1, igraph::E(gr)$count_2), na.rm = TRUE)
      igraph::E(gr)$weight <- rowSums(cbind(igraph::E(gr)$weight_1, igraph::E(gr)$weight_2), na.rm = TRUE)
      igraph::E(gr)$rho <- rowSums(cbind(igraph::E(gr)$rho_1, igraph::E(gr)$rho_2), na.rm = TRUE)
    }
    gr <- igraph::delete_edge_attr(gr, "weight_1")
    gr <- igraph::delete_edge_attr(gr, "weight_2")
    gr <- igraph::delete_edge_attr(gr, "rho_1")
    gr <- igraph::delete_edge_attr(gr, "rho_2")
    gr <- igraph::delete_edge_attr(gr, "count_1")
    gr <- igraph::delete_edge_attr(gr, "count_2")
  }

  matnames <- c(list(rownames(M[[1]])), lapply(M, colnames))
  vlevel <- sapply(igraph::V(gr)$name, grep, matnames)

  ## create Sankey plot
  nv <- length(igraph::V(gr))
  col1 <- rep(RColorBrewer::brewer.pal(12, "Set3"), 100)[1:nv]
  ee <- igraph::get.edges(gr, igraph::E(gr)) - 1

  ee.label <- paste(
    "rho=", round(igraph::E(gr)$rho, 2),
    "<br>count=", round(igraph::E(gr)$count, 2), "%",
    "<br>weight=", round(igraph::E(gr)$weight, 2)
  )

  nodes <- data.frame(label = igraph::V(gr)$name, color = col1)
  nodes$info <- paste(igraph::V(gr)$name, unlist(labels)[igraph::V(gr)$name])
  nodes$label <- sub(".*[:]", "", nodes$label)
  nodes$x <- (vlevel - 1) / max(vlevel - 1)

  if (fill) {
    wt <- igraph::E(gr)$weight
    ev2 <- 0.05 + 0.55 * (wt / max(wt))
    col2 <- paste("rgba(80,80,120,", ev2, ")")
  } else {
    col2 <- paste("rgba(80,80,120,0.2)")
  }

  links <- data.frame(
    source = ee[, 1],
    target = ee[, 2],
    count = igraph::E(gr)$count,
    label = ee.label,
    weight = igraph::E(gr)$weight,
    color = col2
  )

  fig <- plotly::plot_ly(
    type = "sankey",
    domain = list(
      x =  c(0, 1),
      y =  c(0, 1)
    ),
    orientation = "h",
    valueformat = ".2f",
    valuesuffix = "",
    arrangement = "snap",
    node = list(
      label = nodes$label,
      x = nodes$x,
      y = 0.01 * (1:length(nodes$x)),
      color = nodes$color,
      pad = 15,
      thickness = 15,
      line = list(
        color = "black",
        width = 0.5
      )
    ),
    link = list(
      source = links$source,
      target = links$target,
      value =  links$count,
      label =  links$label,
      color =  links$color
    )
  )
  fig
}


#' @title Create a sankey diagram from phenotypes
#'
#' @param pgx PGX object with results
#' @param phenotypes Data frame of sample phenotypes
#' @param mat The results matrix to use, default is pgx$MAT
#' @param fill Whether to color edges by matrix values
#' @param nmin Minimum sample size for phenotype levels
#' @param title Plot title
#'
#' @return A plotly sankey diagram object
#'
#' @description
#' Generates a sankey diagram from sample phenotypes
#'
#' @details This function takes a PGX object and sample phenotype data frame.
#' It aggregates the results matrix mat by phenotypes, filters low count levels,
#' calculates a cross-table, and visualizes the relationships between phenotypes
#' as a sankey diagram.
#'
#' Useful for exploring sample relationships and phenotype flows.
#'
#' @export
pgx.SankeyFromPhenotypes.PLOTLY <- function(pgx, phenotypes, mat = NULL,
                                            fill = NULL, nmin = 1, title = "") {
  prefix1 <- "from"
  prefix2 <- "to"
  table2graph <- function(tab, prefix1, prefix2) {
    rownames(tab) <- paste0(prefix1, ":", rownames(tab))
    colnames(tab) <- paste0(prefix2, ":", colnames(tab))
    n <- nrow(tab) + ncol(tab)
    A <- matrix(0, nrow = n, ncol = n)
    rownames(A) <- colnames(A) <- c(rownames(tab), colnames(tab))
    A[rownames(tab), colnames(tab)] <- tab
    gr <- igraph::graph_from_adjacency_matrix(A, mode = "directed", diag = FALSE, weighted = TRUE)

    gr
  }

  if (is.null(mat)) {
    phenotypes <- intersect(phenotypes, colnames(pgx$samples))
    A <- pgx$samples[, phenotypes]
  }
  if (!is.null(mat)) {
    A <- mat
    phenotypes <- colnames(mat)
  }

  gr.list <- list()
  for (i in 1:(ncol(A) - 1)) {
    p1 <- A[, phenotypes[i]]
    p2 <- A[, phenotypes[i + 1]]
    P <- table(p1, p2)
    P[which(P < nmin, arr.ind = TRUE)] <- 0
    gr.list[[i]] <- table2graph(P, phenotypes[i], phenotypes[i + 1])
  }

  gr <- gr.list[[1]]
  if (length(gr.list) > 1) {
    i <- 2
    for (i in 2:length(gr.list)) {
      gr <- igraph::union(gr, gr.list[[i]])
      igraph::E(gr)$weight <- rowSums(cbind(igraph::E(gr)$weight_1, igraph::E(gr)$weight_2), na.rm = TRUE)
    }
    gr <- igraph::delete_edge_attr(gr, "weight_1")
    gr <- igraph::delete_edge_attr(gr, "weight_2")
  }

  ## create Sankey plot
  nv <- length(igraph::V(gr))
  col1 <- rep(RColorBrewer::brewer.pal(12, "Set3"), 100)[1:nv]
  ee <- igraph::ends(gr, igraph::E(gr))
  ee <- apply(ee, 2, function(v) match(v, igraph::V(gr)$name) - 1)

  nodes <- data.frame(label = igraph::V(gr)$name, color = col1)
  nodes$label <- sub(".*[:]", "", nodes$label)
  links <- data.frame(
    source = ee[, 1],
    target = ee[, 2],
    value = igraph::E(gr)$weight,
    label = ""
  )

  fig <- plotly::plot_ly(
    type = "sankey",
    domain = list(
      x =  c(0, 1),
      y =  c(0, 1)
    ),
    orientation = "h",
    valueformat = ".0f",
    valuesuffix = "cells",
    node = list(
      label = nodes$label,
      color = nodes$color,
      pad = 15,
      thickness = 15,
      line = list(
        color = "black",
        width = 0.5
      )
    ),
    link = list(
      source = links$source,
      target = links$target,
      value = links$value,
      label = links$label,
      color = "rgba(80,80,120,0.15)"
    )
  )
  fig
}


#' @describeIn pgx.SankeyFromPhenotypes.PLOTLY pgx.SankeyFromPhenotypes.GGPLOT (ggplot2 version)
#' @export
pgx.SankeyFromPhenotypes.GGPLOT <- function(pgx, phenotypes, mat = NULL, fill = NULL,
                                            sort = FALSE, nmin = 1, title = "") {
  if (is.null(mat)) {
    phenotypes <- intersect(phenotypes, colnames(pgx$samples))
    phenotypes <- Matrix::head(phenotypes, 5)
    pp <- intersect(c(phenotypes, fill), colnames(pgx$samples))
    A <- pgx$samples[, pp]
  }
  if (!is.null(mat)) {
    A <- mat
    phenotypes <- colnames(mat)
  }

  pt <- apply(A, 1, paste, collapse = ":")
  N <- table(pt)
  df <- do.call(rbind, strsplit(names(N), split = ":"))
  colnames(df) <- phenotypes
  zap.tiny <- (N < nmin)
  if (sum(zap.tiny) > 0) {
    n0 <- sum(N[zap.tiny])
    df <- df[!zap.tiny, ]
    N <- N[!zap.tiny]
  }
  df <- data.frame(df)
  df$Frequency <- N

  if (sort) {
    relevelBig2small <- function(idx) factor(idx, levels = names(sort(table(idx))))
    df[, phenotypes] <- data.frame(lapply(df[, phenotypes], relevelBig2small))
  }

  v <- sapply(phenotypes, function(s) {
    # The f sym likely comes from x package, but better confirm
    a <- sym(s)
    # The f enquo likely comes from rlang, but better confirm
    enquo(a)
  })
  if (length(phenotypes) == 2) {
    p <- ggplot2::ggplot(df, ggplot2::aes(y = Frequency, axis1 = !!v[[1]], axis2 = !!v[[2]]))
  }
  if (length(phenotypes) == 3) {
    p <- ggplot2::ggplot(df, ggplot2::aes(
      y = Frequency, axis1 = !!v[[1]], axis2 = !!v[[2]],
      axis3 = !!v[[3]]
    ))
  }
  if (length(phenotypes) == 4) {
    p <- ggplot2::ggplot(df, ggplot2::aes(
      y = Frequency, axis1 = !!v[[1]], axis2 = !!v[[2]],
      axis3 = !!v[[3]], axis4 = !!v[[4]]
    ))
  }
  if (length(phenotypes) == 5) {
    p <- ggplot2::ggplot(df, ggplot2::aes(
      y = Frequency, axis1 = !!v[[1]], axis2 = !!v[[2]],
      axis3 = !!v[[3]], axis4 = !!v[[4]], axis5 = !!v[[5]]
    ))
  }
  if (length(phenotypes) == 6) {
    p <- ggplot2::ggplot(df, ggplot2::aes(
      y = Frequency, axis1 = !!v[[1]], axis2 = !!v[[2]],
      axis3 = !!v[[3]], axis4 = !!v[[4]], axis5 = !!v[[5]],
      axis6 = !!v[[6]]
    ))
  }
  if (!is.null(fill)) {
    a <- ggplot2::sym(fill)
    vf <- ggplot2::enquo(a)
    p <- p + geom_alluvium(ggplot2::aes(fill = !!vf, color = !!vf), width = 1 / 12, alpha = 0.4)
  } else {
    p <- p + geom_alluvium(width = 1 / 12, alpha = 0.4)
  }
  p <- p +
    geom_stratum(width = 1 / 3, fill = "grey", color = "grey50") +
    ggplot2::geom_label(stat = "stratum", ggplot2::aes(label = ggplot2::after_stat(stratum))) +
    ggplot2::scale_x_discrete(limits = phenotypes, expand = c(.05, .05)) +
    ggplot2::scale_fill_brewer(type = "qual", palette = "Set1", direction = -1) +
    ggplot2::scale_color_brewer(type = "qual", palette = "Set1", direction = -1) +
    ggplot2::ggtitle(title) + ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 13, vjust = +7)
    )
  p
}




#' @title Plot contrasts from a PGX analysis
#'
#' @description Generate plots for contrasts from a pharmacogenomic (PGX) analysis.
#'
#' @param pgx A PGX object containing the analysis results.
#' @param contrast Character vector of contrast names to plot. Default plots all.
#' @param type Type of plot. Options are "scatter", "volcano", "MA", "UMAP".
#' @param set.par Logical to set graphical parameters before plotting. Default is TRUE.
#' @param par.sq Logical to set square plotting layout. Default is FALSE.
#' @param ... Additional parameters passed to plotting functions.
#'
#' @details This function generates plots to visualize contrasts from a PGX analysis.
#' The \code{pgx} object should contain results for the specified contrasts.
#'
#' The \code{type} parameter determines the type of plot. Options are:
#' - "scatter": Default scatter plot of log2 fold-changes vs. p-values.
#' - "volcano": Volcano plot of log2 fold-changes vs. -log10 p-values.
#' - "MA": Meta-analysis plots showing the contribution of each study.
#' - "UMAP": UMAP projection plot colored by log2 fold-change.
#'
#' If \code{set.par=TRUE}, graphical parameters are set before plotting.
#' If \code{par.sq=TRUE}, a square layout is used.
#' Additional parameters can be passed to the plotting functions via \code{...}.
#'
#' @return
#' A grid of plots visualizing the specified contrasts is generated.
#' The return value is a list of plot grob objects.
#'
#' @export
pgx.plotContrast <- function(pgx, contrast = NULL, type = "scatter",
                             set.par = TRUE, par.sq = FALSE, ...) {
  if (is.null(contrast)) {
    contrast <- colnames(pgx$model.parameters$exp.matrix)
  }

  if (set.par) {
    nc <- ceiling(sqrt(length(contrast)))
    nr <- ceiling(length(contrast) / nc)
    if (par.sq) nr <- nc
    graphics::par(mfrow = c(nr, nc))
  }

  plist <- list()
  for (i in 1:length(contrast)) {
    ct <- contrast[i]

    if (type == "volcano") {
      p <- pgx.Volcano(pgx, contrast = ct, ...)
    } else if (type == "MA") {
      p <- pgx.plotMA(pgx, contrast = ct, ...)
    } else if (type == "UMAP") {
      p <- pgx.plotGeneUMAP(pgx, contrast = ct, set.par = FALSE, ...)
    } else {
      ## scatter
      p <- pgx.contrastScatter(pgx, contrast = ct, ...)
    }
    plist[[i]] <- p
  }

  if (length(plist) == 1) plist <- plist[[1]]
  plist
}


#' @describeIn pgx.plotContrast Create a volcano plot from a PGX object
#'
#' @param pgx A PGX object containing differential expression results.
#' @param contrast The contrast name to extract results for.
#' @param level The data level to extract ("gene", "exon", etc). Default "gene".
#' @param methods The meta-analysis methods to include. Default "meta".
#' @param psig P-value cutoff for significance. Default 0.05.
#' @param fc Fold change cutoff. Default 1.
#' @param cex Point size. Default 1.
#' @param cex.lab Label size. Default 1.
#' @param ntop Number of top genes to highlight. Default 20.
#' @param p.min Minimum p-value for y-axis. Default NULL.
#' @param fc.max Maximum fold change for x-axis. Default NULL.
#' @param hilight Vector of genes to highlight. Default NULL.
#' @param cpal Color palette. Default c("grey60", "red3").
#' @param title Plot title. Default NULL.
#' @param plotlib Plotting library to use. Default "base".
#'
#' @export
pgx.Volcano <- function(pgx, contrast, level = "gene", methods = "meta",
                        psig = 0.05, fc = 1, cex = 1, cex.lab = 1, ntop = 20,
                        p.min = NULL, fc.max = NULL, hilight = NULL, #
                        cpal = c("grey60", "red3"), title = NULL,
                        plotlib = "base") {
  if (is.integer(contrast)) contrast <- names(pgx$gx.meta$meta)[contrast]
  res <- NULL
  if (level == "gene") {
    res <- pgx.getMetaMatrix(pgx, methods, level = "gene")
  } else if (level == "geneset") {
    res <- pgx.getMetaMatrix(pgx, methods, level = "geneset")
  } else {
    stop("FATAL:: invalid level=", level)
  }
  f <- res$fc[, contrast]
  q <- res$qv[, contrast]

  sig <- (q <= psig & abs(f) >= fc)
  xy <- cbind(fc = f, y = -log10(q))



  if (is.null(hilight)) {
    wt <- rowSums(scale(xy, center = FALSE)**2)
    hilight <- rownames(xy)[order(-wt)]
    hilight <- intersect(hilight, names(sig[sig == TRUE]))
  }
  hilight <- Matrix::head(hilight, ntop) ## label

  xlim <- ylim <- NULL
  if (!is.null(fc.max)) {
    xlim <- c(-1.1, 1.1) * fc.max
  }
  if (!is.null(p.min)) {
    ylim <- c(0, -log10(p.min))
  }

  if (is.null(title)) title <- contrast
  p <- pgx.scatterPlotXY(
    xy,
    var = sig, type = "factor", title = title,
    xlab = "differential expression (log2FC)",
    ylab = "significance (-log10q)",
    hilight = hilight, #
    cex = cex, cex.lab = cex.lab, cex.title = 1.0,
    xlim = xlim, ylim = ylim,
    legend = FALSE, col = cpal, opacity = 1,
    plotlib = plotlib
  )

  p
}


#' @describeIn pgx.plotContrast Meta-analysis plots for PGX results
#' @param pgx PGX object containing results
#' @param contrast Contrast name to plot
#' @param level Data level to extract ("gene", "exon", etc)
#' @param psig P-value cutoff for significance
#' @param fc Fold change cutoff
#' @param cex Point size
#' @param cex.lab Label size
#' @param hilight Vector of genes to highlight
#' @param ntop Number of top genes to label
#' @param plotlib Plotting library to use
#' @export
pgx.plotMA <- function(pgx, contrast, level = "gene", psig = 0.05, fc = 1,
                       cex = 1, cex.lab = 0.8, hilight = NULL, ntop = 20,
                       plotlib = "base") {
  if (is.integer(contrast)) contrast <- names(pgx$gx.meta$meta)[contrast]

  if (level == "gene") {
    f <- pgx$gx.meta$meta[[contrast]]$meta.fx
    q <- pgx$gx.meta$meta[[contrast]]$meta.q
    gg <- rownames(pgx$gx.meta$meta[[contrast]])
    names(f) <- names(q) <- gg
    m <- rowMeans(pgx$X[gg, ])
  } else if (level == "geneset") {
    f <- pgx$gset.meta$meta[[contrast]]$meta.fx
    q <- pgx$gset.meta$meta[[contrast]]$meta.q
    gg <- rownames(pgx$gset.meta$meta[[contrast]])
    names(f) <- names(q) <- gg
    m <- rowMeans(pgx$gsetX[gg, ])
  } else {
    stop("FATAL:: invalid level=", level)
  }

  sig <- (q <= psig & abs(f) >= fc)
  xy <- cbind(x = m, y = f)
  # <- gg
  cpal <- c("grey60", "red3")
  if (is.null(hilight)) {
    wt <- rowSums(scale(cbind(xy, m), center = FALSE)**2)
    hilight <- rownames(xy)[order(-wt)]
    hilight <- intersect(hilight, names(sig[sig == TRUE]))
  }
  hilight <- Matrix::head(hilight, ntop)

  p <- pgx.scatterPlotXY(
    xy,
    var = sig, type = "factor", title = contrast,
    xlab = "average expression  (log2)",
    ylab = "differential expression  (log2FC)",
    hilight = hilight, cex = 0.9 * cex,
    cex.lab = cex.lab, cex.title = 1.0,
    legend = FALSE, col = c("grey70", "red3"), opacity = 1,
    plotlib = plotlib
  )

  p
}


#' @describeIn pgx.plotContrast Scatter plot of contrast results
#' @param pgx PGX object with analysis results
#' @param contrast Contrast name to plot
#' @param hilight Vector of genes to highlight
#' @param cex Point size
#' @param cex.lab Label size
#' @param psig P-value cutoff for significance
#' @param fc Fold change cutoff
#' @param level Data level to extract ("gene", "exon", etc)
#' @param ntop Number of top genes to label
#' @param dir Direction to sort results (0 for p-value, 1 for logFC)
#' @param plotlib Plotting library ("base", "plotly")
#'
#' @export
pgx.contrastScatter <- function(pgx, contrast, hilight = NULL,
                                cex = 1, cex.lab = 0.8,
                                psig = 0.05, fc = 1, level = "gene",
                                ntop = 20, dir = 0, plotlib = "base") {
  if (is.numeric(contrast)) contrast <- names(pgx$gx.meta$meta)[contrast]
  exp.matrix <- pgx$model.parameters$exp.matrix
  ct <- exp.matrix[, contrast]
  ii <- which(ct < 0)
  jj <- which(ct > 0)
  if (level == "gene") {
    x0 <- rowMeans(pgx$X[, ii, drop = FALSE])
    x1 <- rowMeans(pgx$X[, jj, drop = FALSE])
    xy <- cbind(x0, x1)
    gg <- rownames(pgx$gx.meta$meta[[contrast]])
    fx <- pgx$gx.meta$meta[[contrast]]$meta.fx
    q <- pgx$gx.meta$meta[[contrast]]$meta.q
    names(fx) <- names(q) <- gg
  } else if (level == "gene") {
    x0 <- rowMeans(pgx$gsetX[, ii, drop = FALSE])
    x1 <- rowMeans(pgx$gsetX[, jj, drop = FALSE])
    xy <- cbind(x0, x1)
    gg <- rownames(pgx$gset.meta$meta[[contrast]])
    fx <- pgx$gset.meta$meta[[contrast]]$meta.fx
    q <- pgx$gset.meta$meta[[contrast]]$meta.q
    names(fx) <- names(q) <- gg
  } else {
    stop("FATAL:: invalid level=", level)
  }

  colnames(xy) <- c("", " ")
  grp1 <- gsub(".*[:]|_vs_.*", "", contrast)
  grp0 <- gsub(".*_vs_", "", contrast)
  xlab <- paste("expression", grp0, "  (logCPM)")
  ylab <- paste("expression", grp1, "  (logCPM)")

  if (is.null(hilight)) {
    top.gg <- c(
      Matrix::head(names(sort(fx)), ntop / 2),
      Matrix::tail(names(sort(fx)), ntop / 2)
    )
    if (dir > 0) {
      top.gg <- Matrix::tail(names(sort(fx)), ntop)
    }
    if (dir < 0) {
      top.gg <- Matrix::head(names(sort(fx)), ntop)
    }
    hilight <- top.gg
  }
  hilight <- Matrix::head(hilight, ntop)

  sig <- 1 * (q < psig & abs(fx) > fc)
  names(sig) <- gg


  tt <- contrast

  pgx.scatterPlotXY(
    xy,
    var = sig, type = "factor", title = tt,
    xlab = xlab, ylab = ylab,
    hilight = hilight, cex = 0.9 * cex,
    cex.lab = cex.lab, cex.title = 1.0,
    legend = FALSE, col = c("grey70", "red3"), opacity = 1,
    plotlib = plotlib
  )
}


#' @describeIn pgx.plotContrast Plot gene expression on UMAP projection
#'
#' @param pgx PGX object with results
#' @param contrast Contrast name to extract expression
#' @param value Gene name or ID to plot
#' @param pos Dataframe of UMAP coordinates
#' @param ntop Number of top genes to label
#' @param cex Point size
#' @param cex.lab Label size
#' @param hilight Vector of genes to highlight
#' @param title Plot title
#' @param zfix Fix color scale between genes
#' @param set.par Set graphical parameters before plotting
#' @param par.sq Use square plotting layout
#' @param level Data level to extract ("gene", "exon", etc)
#' @param plotlib Plotting library ("base", "ggplot")
#'
#' @export
pgx.plotGeneUMAP <- function(pgx, contrast = NULL, value = NULL,
                             pos = NULL, ntop = 20, cex = 1, cex.lab = 0.8,
                             hilight = NULL, title = NULL, zfix = FALSE,
                             set.par = TRUE, par.sq = FALSE,
                             level = "gene", plotlib = "ggplot") {
  if (!is.null(contrast)) {
    if (is.numeric(contrast)) contrast <- names(pgx$gx.meta$meta)[contrast]
    res <- NULL
    if (level == "gene") {
      res <- pgx.getMetaMatrix(pgx, level = "gene")
    } else if (level == "geneset") {
      res <- pgx.getMetaMatrix(pgx, level = "geneset")
    } else {
      stop("FATAL:: invalid level=", level)
    }
    F <- res$fc[, contrast, drop = FALSE]
  }

  if (!is.null(value)) {
    F <- cbind(value)
  }

  if (!is.null(pos)) {
    xy <- pos
  } else if ("cluster.genes" %in% names(pgx)) {
    xy <- pgx$cluster.genes$pos[["umap2d"]]
  } else {
    message("[WARNING] no cluster.genes in object. please supply positions")
    return(NULL)
  }

  F <- F[match(rownames(xy), rownames(F)), , drop = FALSE]

  rownames(F) <- rownames(xy)

  if (set.par) {
    nc <- ceiling(sqrt(ncol(F)))
    nr <- ceiling(ncol(F) / nc)
    if (par.sq) nr <- nc
    graphics::par(mfrow = c(nr, nc))
  }

  ## z-scale
  zlim <- NULL
  if (zfix) {
    zlim <- stats::quantile(F, probs = c(0.01, 0.99), na.rm = TRUE)
    zlim <- stats::quantile(F, probs = c(0.002, 0.998), na.rm = TRUE)
    zlim
  }

  plist <- list()
  i <- 1
  for (i in 1:ncol(F)) {
    title1 <- contrast[i]
    if (is.null(title)) title1 <- paste(title, contrast[i])

    f1 <- F[, i]

    hilight1 <- hilight
    if (is.null(hilight)) {
      hilight1 <- names(sort(-abs(f1)))
    }
    hilight1 <- Matrix::head(hilight1, ntop) ## label
    opacity <- ifelse(length(hilight1) > 0, 0.66, 1)

    p1 <- pgx.scatterPlotXY(
      xy,
      var = f1, type = "numeric",
      xlab = "UMAP-x  (genes)",
      ylab = "UMAP-y  (genes)",
      hilight = hilight1,
      zlim = zlim, zsym = TRUE, softmax = 1,
      cex = cex, cex.lab = cex.lab,
      title = title1, cex.title = 1.0,
      legend = TRUE,
      opacity = 0.5,
      plotlib = plotlib
    )

    plist[[i]] <- p1
  }
  if (plotlib == "base") {
    return()
  }
  if (length(plist) == 1) plist <- plist[[1]]
  return(plist)
}


#' @title Plot expression for a gene
#'
#' @param pgx PGX object with expression data
#' @param probe Gene name or ID to plot
#' @param comp Sample grouping variable
#' @param logscale Log2 transform expression values
#' @param level Data level to extract ("gene", "exon", etc)
#' @param grouped Group samples
#' @param srt Label rotation angle
#' @param cex Point size
#' @param collapse.others Collapse non-highlighted samples
#' @param showothers Show non-highlighted samples
#' @param max.points Maximum number of points to plot
#' @param group.names Group names for legend
#' @param main Plot title
#' @param xlab X axis label
#' @param ylab Y axis label
#' @param names Show sample names
#' @param plotly.annotations Plotly annotations
#' @param plotly.margin Plotly margins
#' @param plotlib Plotting library
#'
#' @return A gene expression plot
#'
#' @description
#' Generate a gene expression plot from a PGX object.
#'
#' @details
#' This function extracts the expression for a given \code{probe}
#' from the \code{pgx} object at the specified \code{level}.
#' It generates a dot plot with samples grouped on the x-axis.
#' Use \code{comp} to specify grouping.
#'
#' Options to collapse non-highlighted samples, log-transform values,
#' show sample names, set axis labels, etc.
#'
#' @export
pgx.plotExpression <- function(pgx, probe, comp, logscale = TRUE,
                               level = "gene", grouped = FALSE, srt = NULL, cex = 1,
                               collapse.others = TRUE, showothers = TRUE,
                               max.points = 200, group.names = NULL,
                               main = NULL, xlab = NULL, ylab = NULL, names = TRUE,
                               plotly.annotations = NULL,
                               plotly.margin = NULL,
                               plotlib = "base") {
  if (is.null(probe)) {
    return(NULL)
  }
  if (is.na(probe)) {
    return(NULL)
  }

  if (level == "geneset" && !probe %in% rownames(pgx$gsetX)) {
    graphics::frame() ## emtpy image
    return(NULL)
  }

  ## ------------- determine groups
  expmat <- pgx$model.parameters$exp.matrix
  cntrmat <- pgx$model.parameters$contr.matrix
  expmat <- expmat[rownames(pgx$samples), , drop = FALSE]

  if (inherits(comp, "numeric")) comp <- colnames(expmat)[comp]
  if (!is.null(group.names) && length(group.names) != 2) stop("group.names must be length=2")
  if (is.null(main)) main <- probe
  comp

  if (!is.null(group.names)) {
    message("[pgx.plotExpression] using group names from argument")
  }

  ## if a named contrast table is available it is safer
  if (is.null(group.names) && "contrasts" %in% names(pgx)) {
    message("[pgx.plotExpression] parsing group names from pgx$contrast labels")
    contr.labels <- pgx$contrasts[, comp]
    contr.idx <- expmat[, comp]
    group1 <- names(which.max(table(contr.labels[contr.idx > 0])))
    group0 <- names(which.max(table(contr.labels[contr.idx < 0])))
    group.names <- c(group0, group1)
  }
  group.names

  ## Otherwise we guess from the contrast title but this is dangerous
  if (is.null(group.names) && grepl("_vs_|_VS_", comp)) {
    message("[pgx.plotExpression] parsing group names contrast name (warning!)")
    comp1 <- sub(".*:", "", comp) ## remove prefix
    group.names <- strsplit(comp1, split = "_vs_|_VS_")[[1]]
    group.names <- rev(group.names) ## first is main group
    ## !!!!!!!! NEED RETHINK !!!!!!!! which one is reference
    ## class???  determine if notation is A_vs_B or B_vs_A (some
    ## people do different) and in sample-based contrasts we don't
    ## know group names.
    if (!is.POSvsNEG(pgx)) {
      ## A_vs_B or B_vs_A notation !!!
      message("[pgx.plotExpression] WARNING! A_vs_B reversed?")
      group.names <- rev(group.names) ## reversed!!
    }
    group.names <- gsub("@.*", "", group.names) ## strip postfix
  }
  group.names

  ## label rotation
  if (is.null(srt)) {
    srt <- 0
    if (max(nchar(group.names)) >= 7) srt <- 45
  }

  ## create groups
  ct <- expmat[, comp]
  names(ct) <- rownames(expmat)
  samples <- rownames(expmat)[which(ct != 0)]
  grp0.name <- grp1.name <- NULL
  if (!is.null(group.names)) {
    grp0.name <- group.names[1]
    grp1.name <- group.names[2]
  } else {
    grp1.name <- comp
    grp0.name <- "REF"
  }

  xgroup <- c("other", grp0.name, grp1.name)[1 + 1 * (ct < 0) + 2 * (ct > 0)]
  names(xgroup) <- rownames(pgx$samples)
  jj <- which(!(xgroup %in% xgroup[samples]))

  if (length(jj) > 0 && collapse.others) {
    xgroup <- as.character(xgroup)
    xgroup[jj] <- "other"
  }
  names(xgroup) <- rownames(expmat)

  if (inherits(xgroup, "character")) {
    xgroup <- as.character(xgroup)

    levels0 <- group.names
    if ("other" %in% xgroup) levels0 <- c(levels0, "other")
    xgroup <- factor(xgroup, levels = levels0)
  }

  ## ------------- set color of samples
  ngrp <- length(group.names)
  grp.klr <- rep(RColorBrewer::brewer.pal(12, "Paired"), 99)[1:ngrp]
  names(grp.klr) <- group.names

  if (any(grepl("other", xgroup))) {
    grp.klr <- c("other" = "#d9d9d9", grp.klr)
  }

  ## -------------- get expression value
  if (level == "geneset") {
    gx <- pgx$gsetX[probe, rownames(pgx$samples)]
  } else {
    gx <- pgx$X[rownames(pgx$X) == probe, rownames(pgx$samples)]
  }
  
  if (!logscale) {
    gx <- 2**(gx)
  }

  ## -------------- remove others
  if (showothers == FALSE && any(grepl("other", xgroup))) {
    jj <- grep("other", xgroup, invert = TRUE)
    xgroup <- factor(xgroup[jj], levels = setdiff(levels(xgroup), "other"))
    gx <- gx[jj]
  }

  ## -------------- plot grouped or ungrouped
  if (is.null(main)) main <- probe
  ## if(ncol(X) <= 20) {
  if (!grouped) {
    ## NOT GROUPED BARPLOTS
    nx <- length(gx)
    if (is.null(ylab)) {
      ylab <- "expression (log2CPM)"
      if (!logscale) ylab <- "expression (CPM)"
    }
    klr <- grp.klr[as.character(xgroup)]
    klr[is.na(klr)] <- "#e5e5e5"

    if (plotlib == "plotly") {
      fig <- playbase::pgx.barplot.PLOTLY(
        data = data.frame(
          gx = gx,
          xgroup = factor(names(gx), levels = names(gx))
        ),
        x = "xgroup",
        y = "gx",
        grouped = FALSE,
        title = main,
        yaxistitle = ylab,
        xaxistitle = xlab,
        annotations = plotly.annotations,
        margin = plotly.margin
      )
      if (nx > 20 || names == FALSE) { # remove x axis to labels if condition is met
        fig <- plotly::layout(p = fig, xaxis = list(showticklabels = FALSE))
      }
      return(fig)
    } else {
      ## plot using base graphics
      gx.min <- 0
      if (min(gx) < 0) gx.min <- min(gx)
      ylim <- c(gx.min, 1.3 * max(gx))
      bx <- graphics::barplot(gx[],
        col = klr[], ylim = ylim,
        ## offset = 0, ylim=c(gx.min,max(gx)),
        las = 3, ylab = ylab, names.arg = NA, border = NA
      )
    }
  } else {
    ## GROUPED PLOTS
    if (is.null(ylab)) {
      ylab <- "expression (log2CPM)"
      if (!logscale) ylab <- "expression (CPM)"
    }
    bee.cex <- c(0.3, 0.1, 0.05)[cut(length(gx), c(0, 100, 500, 99999))]

    xlevels <- levels(xgroup)
    grp.klr1 <- grp.klr[as.character(xlevels)]
    grp.klr1[is.na(grp.klr1)] <- "e5e5e5"
    names(grp.klr1) <- as.character(xlevels)


    if (plotlib == "plotly") {
      ## plot using plotly
      fig <- playbase::pgx.barplot.PLOTLY(
        data = data.frame(
          gx = gx,
          xgroup = xgroup
        ),
        x = "xgroup",
        y = "gx",
        grouped = TRUE,
        fillcolor = grp.klr1,
        title = main,
        yaxistitle = ylab,
        xaxistitle = xlab,
        annotations = plotly.annotations,
        margin = plotly.margin
      )
      return(fig)
    } else {
      ## plot using base graphics
      gx.b3plot(gx, xgroup, #
        first = xlevels[1],
        col = grp.klr1, ylab = ylab, bee.cex = bee.cex,
        max.points = max.points, xlab = xlab, names = names,
        ## sig.stars=TRUE, max.stars=5,
        las = 3, names.cex = cex, srt = srt
      )
      graphics::title(main, cex.main = 1.0, line = 0)
      return()
    }
  }
}





#' @title Visualize phenotype matrix as a heatmap
#'
#' @param annot Phenotype matrix with samples as rows and variables as columns.
#'
#' @return A heatmap visualization of the phenotype matrix.
#'
#' @description
#' Generates a heatmap to visualize a phenotype matrix.
#'
#' @details
#' This function takes a phenotype matrix \code{annot} with samples as rows and phenotype
#' variables as columns. It converts the columns to factors and visualizes the matrix
#' as a heatmap, with samples in rows and variables in columns.
#'
#' Useful for exploring sample metadata.
#'
#' @export
pgx.plotPhenotypeMatrix <- function(annot) {
  ## ------- set colors
  colors0 <- rep("Set2", ncol(annot))
  names(colors0) <- colnames(annot)

  grid_params <- iheatmapr::setup_colorbar_grid(
    nrows = 3,
    y_length = 0.32,
    y_spacing = 0.39,
    x_spacing = 0.17,
    x_start = 1.1,
    y_start = 0.96
  )

  ## maximize plot area
  mar <- list(l = 150, r = 0, b = 5, t = 0, pad = 3)


  plt <- NULL
  empty.X <- matrix(NA, nrow = 1, ncol = nrow(annot))
  colnames(empty.X) <- rownames(annot)
  plt <- iheatmapr::main_heatmap(
    empty.X,
    name = "phenotype data",
    show_colorbar = FALSE,
    colorbar_grid = grid_params,
    layout = list(margin = mar)
  )

  ## ------- annot need to be factor
  annotF <- data.frame(as.list(annot), stringsAsFactors = TRUE)
  rownames(annotF) <- rownames(annot)
  cvar <- which(sapply(annotF, is.factor))
  annotX <- annotF
  annotX[, cvar] <- data.frame(lapply(annotF[, cvar], as.integer))
  annotX[, ] <- data.frame(lapply(annotX[, ], as.numeric))
  hc <- fastcluster::hclust(stats::dist(annotX))
  plt <- plt %>% iheatmapr::add_col_dendro(hc, size = 20 * ncol(annot) / 6)

  col_annot_height <- 11
  plt <- plt %>%
    iheatmapr::add_col_annotation(
      annotation = annotF[, ],
      size = col_annot_height,
      buffer = 0.005, side = "bottom",
      colors = colors0
    )
  colcex <- 1
  if (nrow(annot) < 100 && colcex > 0) {
    plt <- plt %>% iheatmapr::add_col_labels(
      side = "bottom",
      size = 20 * ncol(annot) / 6,
      font = list(size = 11 * colcex)
    )
  }

  return(plt)
}


#' @title Visualize phenotype matrix as a heatmap
#'
#' @param annot Phenotype matrix with samples as rows and variables as columns.
#' @param annot.ht Height of annotation labels in mm. Default is 5.
#' @param cluster.samples Logical, cluster samples before plotting. Default is TRUE.
#'
#' @return A heatmap visualization of the phenotype matrix.
#'
#' @description
#' Generates a heatmap to visualize a phenotype matrix.
#'
#' @details
#' This function takes a phenotype matrix \code{annot} with samples as rows and phenotype
#' variables as columns. It converts the columns to factors and visualizes the matrix
#' as a heatmap, with samples in rows and variables in columns.
#'
#' Samples are clustered before plotting if \code{cluster.samples=TRUE}.
#' The annotation label height is set by \code{annot.ht}.
#'
#' Usef
#' @export
pgx.plotPhenotypeMatrix0 <- function(annot, annot.ht = 5, cluster.samples = TRUE) {
  cvar <- pgx.getCategoricalPhenotypes(
    annot,
    min.ncat = 2, max.ncat = 10, remove.dup = FALSE
  )
  fvar <- pgx.getNumericalPhenotypes(annot)
  annot.cvar <- annot[, cvar, drop = FALSE]
  annot.fvar <- annot[, fvar, drop = FALSE]
  annot.df <- cbind(annot.fvar, annot.cvar)
  colnames(annot.df) <- paste(colnames(annot.df), "        ")

  if (cluster.samples) {
    annotx <- expandAnnotationMatrix(annot.df)
    hc <- fastcluster::hclust(stats::dist(annotx)) ## cluster samples
    annot.df <- annot.df[hc$order, ]
  }

  npar <- apply(annot.df, 2, function(x) length(setdiff(unique(x), NA)))
  isnum <- c(rep(1, ncol(annot.fvar)), rep(0, ncol(annot.cvar)))
  is.binary <- apply(annot.df, 2, function(x) length(setdiff(unique(x), NA)) == 2)
  is.binary <- apply(annot.df, 2, function(x) all(x %in% c(0, 1, NA, TRUE, FALSE, "T", "F", "NA")))

  ## set colorscale for each annotation parameter
  ann.colors <- list()
  for (i in 1:length(npar)) {
    prm <- colnames(annot.df)[i]
    klrs <- rev(grDevices::grey.colors(npar[i], start = 0.4, end = 0.85)) ## continous scale
    if (npar[i] == 1) klrs <- "#d8d8d8"
    if (npar[i] > 3 && !isnum[i]) klrs <- rep(RColorBrewer::brewer.pal(8, "Set2"), 99)[1:npar[i]]


    names(klrs) <- sort(unique(annot.df[, i]))
    klrs <- klrs[!is.na(names(klrs))]
    ann.colors[[prm]] <- klrs
  }

  show.legend <- (!is.binary & npar <= 10)

  ha <- ComplexHeatmap::HeatmapAnnotation(
    df = annot.df,
    col = ann.colors,
    na_col = "#3b4252",
    simple_anno_size = grid::unit(annot.ht, "mm"), ## BioC 3.8!!
    show_legend = show.legend
  )

  show_colnames <- TRUE
  show_colnames <- (nrow(annot.df) < 100)
  nullmat <- matrix(0, 0, nrow(annot.df))
  colnames(nullmat) <- rownames(annot.df)
  colnames(nullmat) <- paste(colnames(nullmat), "     ")
  h <- ComplexHeatmap::Heatmap(nullmat,
    top_annotation = ha,
    show_column_names = show_colnames
  )

  ## some space between heatmap and annotation
  nullmat2 <- matrix(0, 0, nrow(annot.df) * 0.15)
  h2 <- ComplexHeatmap::Heatmap(nullmat2)

  ComplexHeatmap::draw(h + h2,
    padding = grid::unit(c(1, 10, 1, 10), "mm"),
    adjust_annotation_extension = TRUE
  )
}


#' @title Split heatmap from PGX object
#'
#' @param ngs PGX object containing expression data
#' @param splitx Sample groups to split heatmap columns
#' @param top.mode Method for selecting top rows ("specific", "sd", "pca")
#' @param annot.pheno Dataframe of sample annotations
#' @param row_annot_width Width of row annotation labels
#' @param scale Scaling method for expression data ("row.center", "row.zscore", etc)
#' @param ntop Number of top rows to include
#' @param colors Vector of colors to use for heatmap
#' @param rowcex Row text size scaling factor
#' @param colcex Column text size scaling factor
#'
#' @return A split heatmap grob object
#'
#' @description
#' Generates an interactive split heatmap visualization from a PGX object.
#'
#' @details
#' This function takes a PGX object and generates a split heatmap using the ComplexHeatmap package.
#' The heatmap can be split into groups by a sample metadata variable in \code{splitx}.
#'
#' The top \code{ntop} most variable rows are selected based on \code{top.mode}, which can be "specific", "sd", or "pca".
#' Additional sample annotations can be provided in \code{annot.pheno}.
#'
#' Various graphical parameters like colors, scaling, text sizes can be customized.
#' The output is a ComplexHeatmap grob object containing the split heatmap plot.
#'
#' @export
pgx.splitHeatmap <- function(ngs, splitx = NULL, top.mode = "specific",
                             annot.pheno = NULL, row_annot_width = 0.03,
                             scale = "row.center", ntop = 50, colors = NULL,
                             rowcex = rowcex, colcex = colcex) {
  X0 <- log2(1 + ngs$counts)
  X0[is.na(X0)] <- 0

  if (!is.null(splitx) && splitx[1] %in% colnames(ngs$samples)) {
    splitx <- ngs$samples[, splitx]
  } else {
    top.mode <- "sd"
    splitx <- NULL
  }

  if (top.mode == "pca") {
    cX <- X0 - rowMeans(X0, na.rm = TRUE)
    dr <- svd(cX, nu = 5)$u
    ntop1 <- ceiling(ntop / ncol(dr))
    jj <- as.vector(apply(dr, 2, function(x) Matrix::head(order(-abs(x)), ntop1)))

    X1 <- X0[jj, ]
    idx <- paste0("PC", as.vector(mapply(rep, 1:ncol(dr), ntop1)))
  } else if (top.mode == "specific") {
    grpX <- tapply(colnames(X0), splitx, function(k) rowMeans(X0[, k, drop = FALSE], na.rm = TRUE))
    grpX <- do.call(cbind, grpX)
    cat("dim(grpX)=", dim(grpX), "\n")
    cat("ntop=", ntop, "\n")
    ntop1 <- ceiling(ntop / ncol(grpX))
    grpX <- grpX - rowMeans(grpX, na.rm = TRUE) ## relative to others
    jj <- as.vector(apply(grpX, 2, function(x) Matrix::head(order(-x), ntop1)))
    X1 <- X0[jj, ]
    idx <- paste0("M", as.vector(mapply(rep, 1:ncol(grpX), ntop1)))
  } else {
    X1 <- Matrix::head(X0[order(-apply(X0, 1, stats::sd)), ], ntop)
    hc <- fastcluster::hclust(stats::as.dist(1 - stats::cor(t(X1), use = "pairwise")), method = "ward.D2")
    idx <- paste0("S", stats::cutree(hc, 5))
  }

  ## ----- Get valid phenotype variables
  if (is.null(annot.pheno)) {
    annot.pheno <- grep("^sample|id$|replicate|ratio|year|month|day",
      tolower(colnames(ngs$samples)),
      invert = TRUE
    )
    annot.pheno <- pgx.getCategoricalPhenotypes(ngs$samples, max.ncat = 12, min.ncat = 2)
  }
  annot.pheno <- intersect(annot.pheno, colnames(ngs$samples))
  Y <- ngs$samples[, annot.pheno, drop = FALSE]
  Y <- data.frame(apply(Y, 2, as.character))
  rownames(Y) <- rownames(ngs$samples)

  sampletips <- colnames(X1)
  genetips <- rownames(X1)

  ## ----- call plotting function
  plt <- pgx.splitHeatmapFromMatrix(
    X = X1, annot = Y,
    xtips = genetips, ytips = sampletips,
    idx = idx, splitx = splitx,
    row_annot_width = row_annot_width,
    scale = scale, colors = colors,
    rowcex = rowcex, colcex = colcex
  )
  return(plt)
}

## =================================================================================
## Lower level R level plotting functions
## =================================================================================




#' Convert a plotly plot to a ggplot object
#'
#' @param plot The plotly plot object to convert
#' @param width The width of the resulting ggplot, in pixels. Default NULL.
#' @param height The height of the resulting ggplot, in pixels. Default NULL.
#' @param scale The scaling factor to apply. Default 1.
#' @param hjust The horizontal adjustment to apply. Default 0.
#' @param vjust The vertical adjustment to apply. Default 0.
#'
#' @return A ggplot object representing the converted plotly plot.
#'
#' @details This function takes a plotly plot object and converts it to a ggplot object.
#' The plotly plot is rendered to a temporary PNG file using the orca() function. This PNG
#' is then read in and added as an annotation to a blank ggplot. Optional parameters allow
#' adjusting the width, height, scaling, and justification of the resulting ggplot.
#'
#' @export
plotly2ggplot <- function(plot, width = NULL, height = NULL, scale = 1, hjust = 0, vjust = 0) {
  tmpfile <- tempfile()
  tmpfile <- "tmp.png"
  unlink(tmpfile)
  unlink("tmp.png")
  unlink("tmp_1.png")
  try(plotly::orca(plot, file = tmpfile, width = width, height = height, scale = 2))
  img <- png::readPNG("tmp_1.png")
  img.grob <- grid::rasterGrob(img, interpolate = TRUE)
  ymin <- xmin <- 1 - scale
  xmax <- ymax <- scale
  message("converting to ggplot...")
  gg <- ggplot2::ggplot(data.frame(x = 0:1, y = 0:1), ggplot2::aes_(x = ~x, y = ~y)) +
    ggplot2::geom_blank() +
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(
      0,
      0
    )) +
    ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(
      0,
      0
    )) +
    ggplot2::annotation_custom(img.grob, xmin = xmin +
      hjust, xmax = xmax + hjust, ymin = ymin + vjust, ymax = ymax +
      vjust) +
    ggplot2::theme_void()
  return(gg)
}


#' @title GSEA Enrichment Plot
#'
#' @description Creates an enrichment plot visualization for gene set enrichment analysis.
#'
#' @param fc A numeric vector of gene-level stats (e.g. log2 fold changes) ordered by decreasing statistic.
#' @param gset A character vector of gene ids in the gene set.
#' @param cex The magnification to be used for points and text.
#' @param main The plot title.
#' @param xlab The x-axis label.
#' @param ticklen Length of the tick marks as a fraction of the height of a line of text.
#' @param ylab The y-axis label.
#' @param yth The length of the y-axis tick marks.
#' @param tooltips An optional vector of tooltips to display for each point.
#' @param cex.text The magnification to be used for point labels.
#' @param cex.title The magnification to be used for the title.
#' @param cbar.width The width of the color bar.
#'
#' @details This function takes a ranked list of gene-level stats and a gene set, and produces an enrichment plot visualization using ggplot2.
#' It computes a running enrichment score and plots it against the ranked list. The goal is to visualize if the gene set is enriched at either end of the ranked list.
#'
#' @return A ggplot object containing the enrichment plot.
#'
#' @export
gsea.enplotly <- function(fc, gset, cex = 1, main = NULL, xlab = NULL, ticklen = 0.25,
                          ylab = NULL, yth = 1, tooltips = NULL, cex.text = 1, cex.title = 1.4,
                          cbar.width = 32) {
  if (is.null(xlab)) {
    xlab <- "Rank in ordered dataset"
  }
  if (is.null(ylab)) {
    ylab <- "Ranked list metric"
  }
  if (is.null(main)) {
    main <- "Enrichment plot"
  }

  ## compute running metrix
  fc <- sort(fc, decreasing = TRUE)




  ## weighted cumulative random walk
  x0 <- 1 * (names(fc) %in% gset)
  x0 <- x0 * abs(fc)
  n0 <- sum(!(names(fc) %in% gset))
  n1 <- sum(names(fc) %in% gset)
  r0 <- cumsum(x0 == 0) / sum(x0 == 0)
  r1 <- cumsum(x0) / (1e-4 + sum(x0))

  rnk.trace <- (r1 - r0)
  rnk.trace <- rnk.trace / max(abs(rnk.trace)) * 0.8

  qq <- range(fc)
  y1 <- qq[2]
  y0 <- 0.8 * qq[1]
  dy <- ticklen * (y1 - y0)
  if (max(rnk.trace) >= abs(min(rnk.trace))) rnk.trace <- rnk.trace * abs(y1)
  if (max(rnk.trace) < abs(min(rnk.trace))) rnk.trace <- rnk.trace * abs(y0)

  cc <- sign(fc) * rank(abs(fc))
  df <- data.frame(x = rank(-fc), y = fc, trace = rnk.trace, cc = cc)

  ## downsample
  ii <- which(rownames(df) %in% gset)
  ii <- unique(c(ii, seq(1, nrow(df), round(nrow(df) / 1000))))
  df <- df[ii, ]
  df <- df[order(-df$y), ]

  cpal <- colorspace::diverge_hcl(64, c = 60, l = c(30, 100), power = 1)
  cpal <- colorspace::diverge_hcl(64)

  ## colorbar segments
  db <- nrow(df) / 11
  bb <- round(seq(1, nrow(df), db))
  cbar.x <- df$x[bb]
  cbar.xend <- df$x[c(bb[-1], nrow(df))]

  fx <- df$y[bb + round(0.5 * db)]

  r1 <- (fx / max(abs(fx), na.rm = TRUE))
  r1 <- abs(r1)**0.66 * sign(r1)
  suppressWarnings(
    cc <- gplots::colorpanel(21, "royalblue3", "grey90", "indianred3")
  )
  irnk <- 1 + round((length(cc) - 1) * (1 + r1) * 0.5)

  cbar <- data.frame(x = cbar.x, xend = cbar.xend, color = cc[irnk])

  df$text <- rownames(df)
  jj <- which(rownames(df) %in% gset)
  tooltips2 <- rownames(df)[jj]
  if (!is.null(tooltips)) {
    sel.tt <- match(tooltips2, names(tooltips))
    tooltips2 <- paste0("<b>", tooltips2, "</b><br>", tooltips[sel.tt])
  }

  ii <- seq(1, nrow(df), round(nrow(df) / 200))
  fig <- plotly::plot_ly() %>%
    plotly::add_trace(
      ## -------- grey ordered line
      x = ~ df$x[ii],
      y = ~ df$y[ii],
      type = "scatter",
      mode = "lines",
      fill = "tozeroy",
      fillcolor = "#BBBBBB",
      line = list(color = "#BBBBBB", width = 0),
      hoverinfo = "skip",
      mode = "none"
    ) %>%
    plotly::add_trace(
      ## -------- green score line
      x = ~ df$x, y = ~ df$trace,
      type = "scatter",
      mode = "lines",
      hoverinfo = "skip",
      line = list(
        color = "#00EE00",
        width = cex * 4
      )
    ) %>%
    plotly::add_trace(
      ## -------- black points of (geneset genes)
      x = ~ df$x[jj],
      y = ~ df$y[jj],
      type = "scatter",
      mode = "markers",
      marker = list(
        color = "#444444",
        size = cex * 4
      ),
      text = tooltips2,
      hoveron = "points",
      hoverinfo = "text"
    ) %>%
    plotly::add_segments(
      ## -------- black segments1 (geneset genes)
      x = df$x[jj],
      xend = df$x[jj],
      y = y0 - 0.98 * dy,
      yend = y0,
      type = "scatter",
      mode = "lines",
      line = list(color = "#444444", width = 1.5 * cex),
      text = rownames(df)[jj],
      hoveron = "points",
      hoverinfo = "text"
    )

  ## colorbar/color scale
  for (i in 1:nrow(cbar)) {
    fig <- fig %>%
      plotly::add_segments(
        x = cbar$x[i],
        xend = cbar$xend[i],
        y = y0 - 0.95 * dy,
        yend = y0 - 0.95 * dy,
        ## type = 'scatter', mode='lines',
        line = list(color = cbar$color[i], width = cbar.width)
      )
  }

  if (!is.null(yth) && yth >= 0) {
    ## static labels
    jj2 <- jj[which(abs(df$y[jj]) >= yth)]
    if (length(jj2)) {
      fig <- fig %>%
        plotly::add_annotations(
          x = df$x[jj2],
          y = df$y[jj2],
          text = rownames(df)[jj2],
          xref = "x",
          yref = "y",
          xanchor = "left",
          xshift = 4,
          yshift = 3,
          showarrow = FALSE
        )
    }
  }

  fig <- fig %>%
    plotly::layout(
      font = list(size = 12 * cex.text),
      title = list(text = main, y = 0.99, font = list(size = 12 * cex.title)),
      xaxis = list(title = xlab, gridwidth = 0.3),
      yaxis = list(title = ylab, gridwidth = 0.3, range = c(y0 - 1.1 * dy, y1))
    ) %>%
    plotly::config(toImageButtonOptions = list(format = "svg")) %>%
    plotly::hide_legend()

  fig
}

#' GSEA enrichment plot
#'
#' @param fc A numeric vector of gene-level stats (e.g. log2 fold changes) ordered by decreasing statistic.
#' @param gset A character vector of gene ids in the gene set.
#' @param cex Point expansion factor, passed to ggplot2.
#' @param main Plot title text.
#' @param xlab X axis label text.
#' @param ylab Y axis label text.
#'
#' @return A ggplot2 object containing the enrichment plot.
#'
#' @title Generate a GSEA enrichment plot with ggplot2
#'
#' @description This function takes a ranked list of gene-level stats and a gene set, and produces an enrichment plot visualization using ggplot2.
#'
#' @details It computes a running enrichment score based on the ranked gene list \code{fc} and gene set membership \code{gset}. The enrichment score shows if the gene set is enriched at either end of the ranked list. The plot is generated using ggplot2, with point size scaled by \code{cex}. Axis labels and plot title can be customized.
#'
#' @export
ggenplot <- function(fc, gset, cex = 1, main = NULL, xlab = NULL, ylab = NULL) {
  if (is.null(xlab)) {
    xlab <- "Rank in ordered dataset"
  }
  if (is.null(ylab)) {
    ylab <- "Ranked list metric"
  }
  if (is.null(main)) {
    main <- "Enrichment plot"
  }

  ## compute running metrix
  fc <- sort(fc, decreasing = TRUE)




  ## weighted cumulative random walk
  x0 <- 1 * (names(fc) %in% gset)
  x0 <- x0 * abs(fc)
  n0 <- sum(!(names(fc) %in% gset))
  n1 <- sum(names(fc) %in% gset)
  r0 <- cumsum(x0 == 0) / sum(x0 == 0)
  r1 <- cumsum(x0) / (1e-4 + sum(x0))

  rnk.trace <- (r1 - r0)
  rnk.trace <- rnk.trace / max(abs(rnk.trace)) * 0.8

  qq <- range(fc)
  y1 <- qq[2]
  y0 <- qq[1]
  dy <- 0.2 * (y1 - y0)
  if (max(rnk.trace) >= abs(min(rnk.trace))) rnk.trace <- rnk.trace * abs(y1)
  if (max(rnk.trace) < abs(min(rnk.trace))) rnk.trace <- rnk.trace * abs(y0)

  cc <- sign(fc) * rank(abs(fc))
  df <- data.frame(rank = rank(-fc), fc = fc, run = rnk.trace, cc)
  jj <- which(names(fc) %in% gset)


  cpal <- colorspace::diverge_hcl(64, c = 60, l = c(30, 100), power = 1)
  ii <- 1 + 32 + range(round(32 * (fc / max(abs(fc)))))
  cpal <- colorspace::diverge_hcl(65)[ii[1]:ii[2]]

  cex.title <- 1

  ggplot2::ggplot(data = df, ggplot2::aes(x = rank, y = fc, color = cc)) +
    ggplot2::geom_segment(ggplot2::aes(x = rank, y = 0, xend = rank, yend = fc), color = "grey70") +
    ggplot2::geom_line(ggplot2::aes(x = rank, y = run), color = "green", size = 0.8) +
    ggplot2::geom_segment(
      data = df[jj, ], color = "black", size = 0.3 * cex,
      ggplot2::aes(x = rank, y = (y0 - dy), xend = rank, yend = (y0))
    ) +
    ggplot2::geom_segment(ggplot2::aes(x = rank, y = y0 - 1 * dy, xend = rank, yend = y0 - 0.6 * dy)) +
    ggplot2::scale_color_gradientn(colors = cpal) +
    ggplot2::ylim(y0 - 1.1 * dy, y1) +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab) +
    ggplot2::ggtitle(main) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(size = 9 * cex.title),
      axis.text.x = ggplot2::element_text(size = 8, vjust = +5),
      axis.text.y = ggplot2::element_text(size = 8, hjust = 0),
      axis.title.x = ggplot2::element_text(size = 9, vjust = +5.5),
      axis.title.y = ggplot2::element_text(size = 9, hjust = +0.5)
    )
}




#' @title Scatter Plot with Filled Colors
#'
#' @description
#' Creates a scatter plot with points colored and sized based on input vectors.
#'
#' @param x Numeric vector for x-axis values.
#' @param y Numeric vector for y-axis values. Default is NULL.
#' @param col Vector of values used to determine point colors. Default is NULL.
#' @param shape Vector of values used to determine point shapes. Default is NULL.
#' @param main Title for the plot. Default is NULL.
#' @param cex Point expansion factor. Default is 1.
#' @param pch Point shape. Default is 20.
#' @param legend Show legend for colors? Default is TRUE.
#' @param xlab Label for x-axis. Default is NULL.
#' @param ylab Label for y-axis. Default is NULL.
#' @param legend.ysp Vertical space for legend. Default is 0.8.
#' @param cex.legend Legend text size. Default is 1.
#' @param barscale Color bar size scaling factor. Default is 0.3.
#' @param opacity Opacity of points. Default is 1.
#' @param gamma Gamma correction factor for colors. Default is 1.
#'
#' @details This function generates a scatter plot from x and y vectors. Points can be colored
#' and sized based on input col and shape vectors. A color bar legend is shown when col is provided.
#' Various graphical parameters like main title, axis labels, legend, and text sizes can be adjusted.
#'
#' @return A ggplot2 scatter plot object is returned.
#'
#' @export
plot_ggscatterFILL <- function(x, y = NULL, col = NULL, shape = NULL,
                               main = NULL, cex = 1,
                               pch = 20, legend = TRUE, xlab = NULL, ylab = NULL,
                               legend.ysp = 0.8, cex.legend = 1,
                               barscale = 0.3, opacity = 1, gamma = 1) {
  if (is.null(y) && NCOL(x) == 2) {
    xlab <- colnames(x)[1]
    ylab <- colnames(x)[2]
    y <- x[, 2]
    x <- x[, 1]
  }
  df <- data.frame(x = x, y = y)
  if (!is.null(col)) df$col <- col
  if (!is.null(shape)) df$shape <- shape

  p <- ggplot2::ggplot(df, ggplot2::aes(x, y, color = col, shape = shape)) +
    ggplot2::geom_point(shape = pch, alpha = opacity, size = 2.0 * cex) +
    ggplot2::ggtitle(main) +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab)
  p
  if (!is.null(col)) {
    cpal <- rev(RColorBrewer::brewer.pal(11, "RdYlBu"))

    zr <- range(col)
    zz <- round(c(zr[1], zr[2]), digits = 2)
    cgamma <- seq(0, 1, 1 / (length(cpal) - 1))**(1 / gamma)
    p <- p +
      ggplot2::guides(
        colour = ggplot2::guide_colourbar(
          barwidth = 1.5 * barscale, barheight = 8 * barscale
        ),
        shape = ggplot2::guide_legend(override.aes = list(size = 2.5 * cex.legend))
      ) +
      ggplot2::scale_color_gradientn(
        colors = cpal, breaks = zz, values = cgamma,
        labels = c(zz[1], zz[2])
      ) +
      ggplot2::expand_limits(color = zr + c(-0.01, 0.01))
  }
  if (legend) {
    p <- p + ggplot2::theme(
      legend.justification = c(0, 0),
      legend.position = c(0.01, 0.01)
    )
  } else {
    p <- p + ggplot2::theme(legend.position = "none")
  }
  p <- p + ggplot2::theme(legend.title = ggplot2::element_blank())
  return(p)
}


#' @title Scatterplot with ggplot2
#'
#' @param x Numeric vector or matrix. If matrix, first column is x values and second column is y values.
#' @param y Optional numeric vector of y values. Ignored if x is a matrix.
#' @param col Vector of colors for points.
#' @param main Title for the plot.
#' @param cex Point expansion factor.
#' @param col.scale Vector mapping colors to a sequential palette.
#' @param shape Vector of shapes for points.
#' @param pch Point shape to use if shape not provided.
#' @param legend Logical to show legend.
#' @param legend.ysp Vertical position of legend (0-1 scale).
#' @param cex.legend Legend text size scaling.
#' @param legend.pos Legend position ("right", "left", "top", "bottom").
#' @param xlab x-axis label.
#' @param ylab y-axis label.
#' @param base_size Base font size for plot text.
#'
#' @return A ggplot scatterplot object
#'
#' @description
#' Generates a scatterplot using ggplot2 with customized colors, shapes, legend, etc.
#'
#' @details
#' This function creates a scatterplot from x and y vectors or a matrix containing x and y in columns 1 and 2.
#' Colors can be specified via the col parameter and mapped to a color palette scale using col.scale.
#' Point shapes can be customized via the shape parameter.
#' The plot can be customized via parameters like legend position/text size, axis labels, and base font size.
#' The plot is rendered using ggplot2 and returned as a ggplot object.
#'
#' @examples
#' \dontrun{
#' x <- rnorm(100)
#' y <- rnorm(100)
#' plot_ggscatter(x, y)
#'
#' mat <- cbind(x, y)
#' plot_ggscatter(mat,
#'   legend.pos = "bottom",
#'   shape = c(1, 2), col.scale = 1:100
#' )
#' }

#' @export
plot_ggscatter <- function(x, y = NULL, col = NULL, main = NULL,
                           cex = 1, col.scale = NULL, shape = NULL, pch = 20,
                           legend = TRUE, legend.ysp = 0.8, cex.legend = 1,
                           legend.pos = "right",
                           xlab = NULL, ylab = NULL, base_size = 12) {
  if (is.null(y) && NCOL(x) == 2) {
    y <- x[, 2]
    x <- x[, 1]
  }
  if (is.null(col.scale)) col.scale <- rep(1:9, 99)
  legend.justification <- NULL
  legend.position <- legend.pos
  if (legend.position == "bottomright") {
    legend.justification <- c(1, 0)
    legend.position <- c(1.0, 0.0)
  } else if (legend.position == "topleft") {
    legend.justification <- c(0, 1)
    legend.position <- c(.0, 1.0)
  } else if (legend.position == "topright") {
    legend.justification <- c(1, 1)
    legend.position <- c(1.0, 1.0)
  } else if (legend.position == "bottomleft") {
    legend.justification <- c(0, 0)
    legend.position <- c(.0, 0.0)
  }

  df <- data.frame(x = x, y = y)
  if (!is.null(col)) df$col <- col
  if (!is.null(shape)) df$shape <- shape

  is.factor <- class(utils::type.convert(as.character(col), as.is = TRUE)) == "factor"
  if (is.factor) {
    p <- ggplot2::ggplot(df, ggplot2::aes(y = y, x = x, color = col, shape = shape)) +
      ggplot2::geom_point(size = 2.0 * cex) +
      ggplot2::scale_color_manual(values = col.scale, name = "") +
      ggplot2::ggtitle(main) +
      ggplot2::xlab(xlab) +
      ggplot2::ylab(ylab)
  } else {
    p <- ggplot2::ggplot(df, ggplot2::aes(y = y, x = x, color = col, shape = shape)) +
      ggplot2::geom_point(size = 2.0 * cex) +
      ggplot2::ggtitle(main) +
      ggplot2::xlab(xlab) +
      ggplot2::ylab(ylab)
    if (!is.null(col.scale)) {
      if (inherits(col, c("numeric", "integer"))) {
        p <- p + ggplot2::scale_color_gradient(low = col.scale[1], high = col.scale[2])
      } else {
        p <- p + ggplot2::scale_fill_gradient(low = col.scale[1], high = col.scale[2])
      }
    }
  }
  p

  if (legend) {
    p <- p + ggplot2::theme(
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 9 * cex.legend),
      legend.key.size = grid::unit(legend.ysp * 0.8, "lines"),
      legend.key = ggplot2::element_rect(color = "transparent", fill = scales::alpha("white", 0.5)),
      legend.justification = legend.justification,
      legend.position = legend.position,
      legend.margin = ggplot2::margin(1, 3, 2, 1),
      legend.box.just = "right",
      legend.box.background = ggplot2::element_rect(color = "#888888", size = 0.25),
      legend.box.margin = ggplot2::margin(1, 2, 1, 1)
    ) +
      ggplot2::guides(
        color = ggplot2::guide_legend(override.aes = list(size = 2.0 * cex.legend)),
        shape = ggplot2::guide_legend(override.aes = list(size = 2.0 * cex.legend))
      )
  } else {
    p <- p + ggplot2::theme(legend.position = "none")
  }

  return(p)
}


#' Violin plot using ggplot2
#'
#' @param x Factor variable for x-axis groups.
#' @param y Numeric vector of values to plot.
#' @param group Optional grouping variable.
#' @param main Plot title.
#' @param ylim Limits for y-axis. Calculated automatically if NULL.
#' @param add.dots Logical to add jittered data points. Default TRUE.
#' @param col Violin fill color. Default "#AAAAAA".
#' @param cex Point expansion factor. Default 1.
#' @param xlab Label for x-axis. Default "".
#' @param ylab Label for y-axis. Default "y".
#' @param srt Rotation angle for x-axis labels. Default 0.
#' @param pdodge Amount to dodge violins when group is used. Default 1.5.
#' @param n.dodge Number of levels to dodge. Default 1.
#' @param base_size Base font size. Default 13.
#'
#' @return A ggplot2 violin plot object.
#' @export
plot_ggviolin <- function(x, y, group = NULL, main = "", ylim = NULL, add.dots = TRUE,
                          col = "#AAAAAA", cex = 1, xlab = "", ylab = "y", srt = 0,
                          pdodge = 1.5, n.dodge = 1, base_size = 13) {
  df <- data.frame(y = y, x = x, group = "")
  if (!is.null(group)) {
    df$group <- group
    col <- RColorBrewer::brewer.pal(8, "Set3")
  }
  if (is.null(ylim)) ylim <- range(y)

  p <- ggplot2::ggplot(df, ggplot2::aes(y = y, x = x, fill = group)) +
    ggplot2::ggtitle(main) +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab) +
    ggplot2::ylim(ylim[1], ylim[2]) +
    ggplot2::scale_fill_manual(values = col) +
    ggplot2::scale_x_discrete(guide = guide_axis(angle = srt)) +
    ggplot2::geom_violin(trim = TRUE, position = ggplot2::position_dodge(pdodge)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = srt, vjust = 0)) +
    ggplot2::theme_minimal(base_size = base_size)
  if (is.null(group)) {
    p <- p + ggplot2::theme(legend.position = "none")
  }
  if (add.dots && is.null(group)) {
    p <- p +
      ggplot2::geom_jitter(shape = 20, size = 1.2 * cex, position = ggplot2::position_jitter(0.07))
  }
  p
}


#' @title Bar plot using ggplot2
#'
#' @description
#' Creates a bar plot visualization using ggplot2.
#'
#' @param mat Numeric data matrix with samples in rows and groups in columns.
#' @param xlab Character label for x-axis. Default is "x".
#' @param ylab Character label for y-axis. Default is "y".
#' @param srt Rotation angle for x-axis labels. Default is 0.
#' @param main Character title for the plot. Default is NULL.
#' @param las Style for axis labels. Default is NULL.
#' @param col Vector of colors to use for bars. Default is NULL.
#' @param beside Logical to draw bars beside each other. Default is FALSE.
#' @param legend.pos Numeric vector of legend position. Default is c(0.016, 1).
#' @param legend.cex Legend text size. Default is 1.
#' @param bar_width Width of bars. Default is 0.7.
#' @param base_size Base font size for plot. Default is 12.
#' @param group.name Label for legend groups. Default is "group".
#'
#' @details
#' This function takes a numeric matrix and generates a bar plot
#' visualization using ggplot2. Bars represent the groups, colored by row names.
#' Many graphical parameters like labels, colors, text sizes, etc. can be adjusted.
#'
#' @return
#' A ggplot bar plot object is returned.
#'
#' @export
plot_ggbarplot <- function(mat, xlab = "x", ylab = "y", srt = 0, main = NULL,
                           las = NULL, col = NULL, beside = FALSE,
                           legend.pos = c(0.016, 1), legend.cex = 1,
                           bar_width = 0.7, base_size = 12, group.name = "group") {
  if (NCOL(mat) == 1) mat <- rbind(mat)
  mat <- mat[nrow(mat):1, , drop = FALSE]
  df <- reshape2::melt(t(mat), value.name = "value")
  colnames(df)[1:2] <- c("x", "y")


  df$y <- factor(df$y, levels = rownames(mat))
  df$x <- factor(df$x, levels = colnames(mat))
  if (!is.null(las) && las == 3) srt <- 90

  cpal <- rev(grDevices::grey.colors(nrow(mat)))

  if (nrow(mat) == 1) cpal <- "grey70"
  if (!is.null(col)) cpal <- rep(col, 99)
  posmode <- ifelse(beside, "dodge", "stack")
  x <- y <- value <- NULL
  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = value, fill = y)) +
    ggplot2::geom_bar(
      stat = "identity", color = "black", size = 0.3,
      width = bar_width, position = posmode
    ) +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab) +
    ggplot2::labs(fill = group.name) +
    ggplot2::ggtitle(main) +
    ggplot2::scale_fill_manual(values = cpal) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(angle = srt)) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = srt, vjust = 0),
      axis.title.x = ggplot2::element_text(size = 10),
      axis.title.y = ggplot2::element_text(size = 10)
    )

  p <- p + ggplot2::theme(
    legend.title = ggplot2::element_blank(),
    legend.justification = legend.pos,
    legend.text = ggplot2::element_text(size = 9 * legend.cex),
    legend.position = legend.pos,
    legend.key.size = grid::unit(9 * legend.cex, "pt"),
    legend.key.height = grid::unit(7 * legend.cex, "pt")
  )

  if (nrow(mat) == 1) {
    p <- p + ggplot2::theme(legend.position = "none")
  }
  p
}


#' @title Violin plot with scatter dots
#'
#' @description Generate a violin plot visualization with embedded scatter dots
#'
#' @param x Factor for x-axis groups
#' @param y Numeric vector of values to plot
#' @param group Optional grouping variable
#' @param xlab Label for x-axis
#' @param ylab Label for y-axis
#' @param srt Rotation angle for x-axis labels
#' @param cex.lab Expansion factor for axis labels
#' @param cex.main Expansion factor for main title
#' @param jitter Jittering factor for scatter dots
#' @param vcol Violin fill color
#' @param plotlib Plotting library to use (base, ggplot, plotly)
#' @param maxbee Maximum number of dots to show in violins
#' @param ... Other arguments passed to plotting functions
#'
#' @return A violin plot
#'
#'
#' @details This function takes a factor \code{x} and numeric vector \code{y} and generates a violin plot.
#' An optional \code{group} variable can be used to produce grouped/facetted violins.
#' Scatter dots are added to each violin with jittering controlled by \code{jitter}.
#' The \code{maxbee} parameter limits the number of dots.
#' The plot can be rendered using base, ggplot2 or plotly graphics depending on \code{plotlib}.
#'
#' @export
pgx.violinPlot <- function(x, y, group = NULL, xlab = "", ylab = "",
                           srt = 0, cex.lab = 1, cex.main = 1.1,
                           jitter = 0.015, vcol = "grey85", main = NULL,
                           plotlib = "base", maxbee = NULL, ...) {
  fig <- NULL
  if (plotlib == "base") {
    if (is.null(maxbee)) maxbee <- 100 * length(setdiff(unique(y), NA))
    if (jitter > 0) x <- x + jitter * diff(range(x)) * stats::rnorm(length(x))


    vioplot::vioplot(x ~ y,
      col = vcol, #
      plotCentre = "line", xaxt = "n", xlab = xlab, ylab = ylab,
      ...
    )
    if (maxbee > 0) {
      ii <- Matrix::head(sample(length(x)), maxbee)
      beeswarm::beeswarm(x[ii] ~ y[ii], add = TRUE, pch = 20, cex = 0.6, col = "grey10")
    }
    yy <- sort(unique(y))
    graphics::text(
      x = 1:length(yy),
      y = graphics::par("usr")[3] - 0.03 * diff(range(x)),
      labels = yy,
      xpd = NA,
      srt = srt,
      adj = ifelse(srt == 0, 0.5, 0.965),
      cex = cex.lab
    )
  }
  if (plotlib == "ggplot") {
    d1 <- !is.null(maxbee) && maxbee > 0
    fig <- plot_ggviolin(y, x,
      add.dots = d1, xlab = xlab, ylab = ylab, srt = srt,
      main = main, cex = cex.main
    )
  }
  if (plotlib == "plotly") {
    df <- data.frame(x = x, y = y, group = 0)
    if (!is.null(group)) df$group <- group
    fig <- df %>%
      plotly::plot_ly(
        x = ~y,
        y = ~x,
        type = "violin",
        box = list(
          visible = TRUE
        ),
        meanline = list(
          visible = TRUE
        )
      ) %>%
      plotly::layout(
        font = list(size = 16 * cex.lab),
        xaxis = list(
          title = xlab
        ),
        yaxis = list(
          title = ylab,
          zeroline = FALSE
        )
      )
  }
  fig
}


#' @title Scatter Plot XY
#'
#' @description This function creates a scatter plot of two variables
#' using one of several available plotting libraries.
#'
#' @param ... Additional arguments passed to the selected plotting function.
#' @param plotlib A character string specifying the plotting library to use.
#'   Possible values are "base", "plotly", "ggplot", and "scatterD3".
#'
#' @details Depending on the value of the `plotlib` parameter, this function
#' calls one of several internal functions to create the scatter plot using the
#' specified plotting library. The available plotting libraries are:
#' - "base": The base R graphics system.
#' - "plotly": The plotly library for interactive plots.
#' - "ggplot": The ggplot2 library for static plots based on the Grammar of Graphics.
#' - "scatterD3": The scatterD3 library for interactive scatter plots using D3.js.
#'
#' @return A plot object created by the selected plotting library.
#'
#' @export
pgx.scatterPlotXY <- function(..., plotlib = "base") {
  if (plotlib == "plotly") {
    pgx.scatterPlotXY.PLOTLY(...)
  } else if (plotlib == "ggplot") {
    pgx.scatterPlotXY.GGPLOT(...)
  } else if (plotlib == "scatterD3") {
    pgx.scatterPlotXY.D3(...)
  } else {
    pgx.scatterPlotXY.BASE(...)
  }
}


#' Scatter plot with base graphics
#'
#' @param var Numeric vector of values to plot on x/y axes.
#' @param type Variable type, either "continuous" or "factor".
#' @param col Colors for points.
#' @param title Plot title.
#' @param zlim Value limits for coloring points by zvar.
#' @param zlog Log scale zvar colors?
#' @param zsym Symmetrize zvar colors?
#' @param softmax Softmax transform zvar colors?
#' @param pch Point shapes.
#' @param cex Point expansion factor.
#' @param cex.lab Label size.
#' @param cex.title Title size.
#' @param cex.legend Legend size.
#' @param zoom Zoom factor.
#' @param legend Show legend?
#' @param bty Legend box type.
#' @param legend.ysp Legend y-spacing.
#' @param legend.pos Legend position.
#' @param lab.pos Label positions, matrix with x, y.
#' @param repel Repel overlapping labels?
#' @param xlab,ylab Axis labels.
#' @param xlim,ylim Axis limits.
#' @param dlim Axis limit expansion fraction.
#' @param hilight2 Secondary highlight indices.
#' @param hilight.cex Size for highlight points.
#' @param lab.xpd Allow labels outside plot area?
#' @param hilight Indices to highlight.
#' @param hilight.col Color for highlighted points.
#' @param hilight.lwd Line width for highlighted points.
#' @param label.clusters Label clusters?
#' @param cex.clust Cluster label size.
#' @param tstep,rstep Label repel increments.
#' @param tooltip Tooltip text.
#' @param theme Theme parameters.
#' @param set.par Reset par?
#' @param axt,xaxs,yaxs Axis type, style, limits.
#' @param labels Point labels.
#' @param label.type Label type.
#' @param opacity Opacity for all elements.
#'
#' @return None. Plot is produced as a side-effect.
#'
#' @export
pgx.scatterPlotXY.BASE <- function(pos, var = NULL, type = NULL, col = NULL, title = "",
                                   zlim = NULL, zlog = FALSE, zsym = FALSE, softmax = FALSE, pch = 20,
                                   cex = NULL, cex.lab = 1, cex.title = 1.2, cex.legend = 1,
                                   zoom = 1, legend = TRUE, bty = "o", legend.ysp = 0.85,
                                   legend.pos = "bottomleft", lab.pos = NULL, repel = TRUE,
                                   xlab = NULL, ylab = NULL, xlim = NULL, ylim = NULL, dlim = 0.05,
                                   hilight2 = hilight, hilight.cex = NULL, lab.xpd = TRUE,
                                   hilight = NULL, hilight.col = NULL, hilight.lwd = 0.8,
                                   label.clusters = FALSE, cex.clust = 1.5,
                                   tstep = 0.1, rstep = 0.1,
                                   tooltip = NULL, theme = NULL, set.par = TRUE,
                                   axt = "s", xaxs = TRUE, yaxs = TRUE,
                                   labels = NULL, label.type = NULL, opacity = 1) {
  ## automatically set pointsize of dots
  if (is.null(cex)) {
    nr <- nrow(pos)
    i <- as.integer(cut(nr, breaks = c(0, 100, 500, 1000, 5000, Inf)))
    cex <- c(2, 1.4, 1, 0.7, 0.4)[i]
  }
  if (is.null(hilight.cex)) hilight.cex <- cex
  if (!is.null(var) && !is.null(ncol(var))) {
    var <- array(var[, 1], dimnames = list(rownames(var)))
  }
  if (!is.null(var) && length(var) == nrow(pos) && is.null(names(var))) {
    names(var) <- rownames(pos)
  }
  if (is.null(var)) {
    var <- as.character(rep("_", nrow(pos)))
    names(var) <- rownames(pos)
  }
  if (is.null(type)) {
    type <- c("numeric", "factor")[1 + class(var) %in% c("factor", "character")]
  }
  if (is.null(colnames(pos))) {
    colnames(pos) <- c("x", "y")
  }
  var <- var[match(rownames(pos), names(var))]

  ## normalize pos
  xlim0 <- range(pos[, 1])
  ylim0 <- range(pos[, 2])
  if (zoom != 1) {
    cx <- mean(range(pos[, 1]))
    cy <- mean(range(pos[, 2]))
    dx <- diff(range(pos[, 1]))
    dy <- diff(range(pos[, 2]))
    xlim0 <- cx + 0.5 * c(-1, 1.05) * dx / zoom
    ylim0 <- cy + 0.5 * c(-1, 1.05) * dy / zoom
  }

  if (length(dlim) == 1) dlim <- rep(dlim, 2)
  xlim0[1] <- xlim0[1] - dlim[1] * diff(xlim0)
  xlim0[2] <- xlim0[2] + dlim[1] * diff(xlim0)
  ylim0[1] <- ylim0[1] - dlim[2] * diff(ylim0)
  ylim0[2] <- ylim0[2] + dlim[2] * diff(ylim0)

  if (is.null(xlab)) xlab <- colnames(pos)[1]
  if (is.null(ylab)) ylab <- colnames(pos)[2]

  if (set.par) {
    par.save <- graphics::par()
    graphics::par(
      mar = c(2.6, 2.8, 1.9, 0.5), mgp = c(1.5, 0.4, 0),
      cex.axis = 1, cex.lab = 1, tcl = -0.3, las = 1
    )
  }

  if (is.null(labels)) {
    labels <- rownames(pos)
  }

  ## Plot the discrete variables
  if (type == "factor") {
    z1 <- factor(var)
    nz <- length(levels(z1))
    if (is.null(col) && nz > 2) {
      col1 <- c(
        RColorBrewer::brewer.pal(8, "Set1"),
        RColorBrewer::brewer.pal(8, "Set2"),
        RColorBrewer::brewer.pal(12, "Set3")
      ) # omics_pal_d("dark")(8))
    } else if (is.null(col) && nz == 2) {
      col1 <- rev(grDevices::grey.colors(2, end = 0.8))
      col1 <- c("#AAAAAA55", "#555555FF")
      col1 <- c("#00008855", "#AA0000FF")
      col1 <- c("#CCCCCC55", "#AA0000FF")
    } else if (is.null(col) && nz == 1) {
      col1 <- c("#22222255")
    } else {
      col1 <- col
    }
    col1 <- Matrix::head(rep(col1, 99), nz)
    pt.col <- col1[z1]
    pt.col[is.na(pt.col)] <- "#DDDDDD33"
    pt.col0 <- pt.col
    if (opacity < 1) {
      pt.col <- add_opacity(pt.col, opacity)
      col1 <- add_opacity(col1, opacity**0.33)
    }

    jj <- order(-table(pt.col)[pt.col]) ## plot less frequent points last...
    if (length(cex) > 1) cex <- cex[jj]
    plot(pos[jj, , drop = FALSE],
      col = pt.col[jj], pch = 20, cex = cex,
      xlim = xlim0, ylim = ylim0,
      xlab = xlab, ylab = ylab,
      axes = FALSE,
      bty = "n"
    )
    if (bty != "n") graphics::box(lwd = 0.8, bty = bty, col = "black")
    graphics::Axis(side = 1, labels = xaxs)
    graphics::Axis(side = 2, labels = yaxs)
    graphics::grid(lwd = 0.8)

    ## label cluster
    if (label.clusters) {
      mpos <- apply(pos, 2, function(x) tapply(x, z1, stats::median))
      mlab <- rownames(mpos)

      graphics::text(mpos, labels = mlab, cex = cex.clust, lheight = 0.8)
    }

    ## discrete parameter legend
    nlev <- length(levels(z1))
    if (legend.pos != "none" && legend && nlev < 30) {
      cex1 <- ifelse(length(levels(z1)) >= 8, 0.85, 1)
      cex1 <- ifelse(length(levels(z1)) >= 15, 0.75, cex1)
      graphics::legend(legend.pos, levels(z1),
        bty = bty, fill = col1,
        ## cex=cex1, y.intersp=0.75, x.intersp=0.7
        cex = 0.85 * cex1 * cex.legend,
        inset = c(0.01, 0.02),
        y.intersp = legend.ysp, x.intersp = 0.75
      )
    }
  }

  ## Plot continous variable
  if (type == "numeric") {
    z <- var

    if (is.null(zlim)) {
      ## global zlim
      zlim <- range(z, na.rm = TRUE)
    }

    if (zsym) {
      zlim <- c(-1, 1) * max(abs(zlim), na.rm = TRUE)
    }

    ## z1 is normalized [0;1] for coloring
    z1 <- (z - min(zlim)) / diff(zlim)
    z1 <- pmin(pmax(z1, 0), 1) ## clip
    if (softmax) {
      z1 <- 0.5 * (tanh(4 * (z1 - 0.5)) + 1)
    }

    ## -------------- set colors

    cpal <- colorspace::diverge_hcl(64, c = 60, l = c(30, 100), power = 1)
    cpal <- grDevices::colorRampPalette(c("#313695", "#FFFFDF", "#A50026"))(64)
    if (!is.null(col)) {
      cpal <- grDevices::colorRampPalette(col)(64)
    } else {
      cpal <- grDevices::colorRampPalette(c("#3136B5", "#FFFFDF", "#B50026"))(64)
    }
    cpal <- sapply(1:length(cpal), function(i) add_opacity(cpal[i], 0.2 + 0.8 * abs(i - 32.5) / 32))

    pt.col <- cpal[1 + ceiling(z1 * (length(cpal) - 1))]


    pt.col0 <- pt.col
    if (opacity < 1) {
      pt.col <- add_opacity(pt.col, opacity)
      cpal <- add_opacity(cpal, opacity**0.33)
    }
    pt.col[is.na(pt.col)] <- "#DDDDDD33" ## missing values color

    jj <- 1:nrow(pos)
    jj <- order(abs(z), na.last = FALSE) ## higher values last??
    if (length(cex) > 1) cex <- cex[jj]
    plot(pos[jj, ],
      col = pt.col[jj],
      pch = 20, cex = cex,
      xlim = xlim0, ylim = ylim0,
      xlab = xlab, ylab = ylab,
      axes = FALSE,
      bty = "n"
    )
    if (bty != "n") graphics::box(lwd = 0.8, bty = bty, col = "black")
    graphics::Axis(side = 1, labels = xaxs)
    graphics::Axis(side = 2, labels = yaxs)
    graphics::grid(lwd = 0.8)

    ## colorscale bar
    if (legend.pos != "none" && legend) {
      zr <- range(z, na.rm = TRUE)
      if (!is.null(zlim)) zr <- zlim
      if (zsym) zr <- c(-1, 1) * max(abs(zr), na.rm = TRUE)
      if (zlog) zr <- round(10**zr - 1) ## ???
      zr <- 0.01 * c(ceiling(100 * zr[1]), floor(100 * zr[2])) ## round
      graphics::legend(legend.pos,
        cex = 0.8 * cex.legend, #
        y.intersp = 0.18, x.intersp = 0.5, border = NA, bty = bty,
        inset = c(0.01, 0.02),
        fill = rev(cpal[seq(from = 1, to = length(cpal), length.out = 11)]),
        legend = c(zr[2], rep(NA, 9), zr[1]),
      )
    }
  }

  if (!is.null(hilight) && length(hilight) > 0) {
    jj <- which(rownames(pos) %in% hilight)
    if (length(jj)) {
      hcol1 <- hilight.col
      if (is.null(hcol1)) hcol1 <- pt.col0[jj]
      hcex <- hilight.cex
      if (length(hcex) > 1) hcex <- hcex[jj]
      graphics::points(pos[jj, , drop = FALSE], pch = 20, col = hcol1, cex = 1.05 * hcex)
      graphics::points(pos[jj, , drop = FALSE], pch = 1, lwd = hilight.lwd, cex = 0.85 * hcex)
    }
  }

  if (!is.null(hilight2) && length(hilight2) > 0) {
    jj <- which(rownames(pos) %in% hilight2)
    if (length(jj)) {
      if (length(cex.lab) == 1) cex.lab <- rep(cex.lab, nrow(pos))
      df <- data.frame(x = pos[jj, 1], y = pos[jj, 2], z = labels[jj], cex = cex.lab[jj])
      if (is.null(lab.pos)) {
        ## repelling text using wordcloud package
        df2 <- rbind(df, df)
        df2$z[1:nrow(df)] <- "x"
        xlim1 <- ylim1 <- c(-Inf, Inf)
        if (!lab.xpd) {
          xlim1 <- xlim0
          ylim1 <- ylim0
        }
        cex.lab2 <- c(rep(1, nrow(df)), df$cex)
        if (repel) {
          nc <- repelwords(df2$x, df2$y, df2$z,
            xlim = xlim1, ylim = ylim1,
            tstep = tstep, rstep = rstep,
            cex = cex.lab2
          )
          nc <- Matrix::tail(nc, nrow(df))
          lab.pos <- data.frame(x = nc[, 1] + .5 * nc[, 3], y = nc[, 2] + .5 * nc[, 4])
        } else {
          lab.pos <- data.frame(x = pos[jj, 1], y = pos[jj, 2])
        }
        rownames(lab.pos) <- rownames(pos)[jj]
      } else {
        lab.pos <- lab.pos[match(rownames(df), rownames(lab.pos)), ]
      }
      graphics::segments(df$x, df$y, lab.pos$x, lab.pos$y, col = "#222222AA", lwd = 0.85)
      graphics::text(lab.pos$x, lab.pos$y, labels = df$z, cex = 0.7 * df$cex)
    }
  }

  ## parameter name
  if (!is.null(title) && title != "") {
    graphics::mtext(title, 3, adj = 0, padj = -0.35, cex = 0.9 * cex.title)
  }

  if (set.par) {

  }
  out <- list()
  out$lab.pos <- lab.pos
  invisible(out)
}


#' Scatter plot with ggplot2
#'
#' @param pos Data frame containing x and y coordinates.
#' @param var Grouping variable.
#' @param type Variable type, either "continuous" or "factor".
#' @param col Colors.
#' @param cex Point expansion factor.
#' @param cex.lab Label size factor.
#' @param cex.title Title size factor.
#' @param cex.clust Cluster label size factor.
#' @param cex.legend Legend text size factor.
#' @param cex.axis Axis text size factor.
#' @param gridcolor Grid line color.
#' @param bgcolor Background color.
#' @param zoom Zoom factor.
#' @param legend Show legend? Logical.
#' @param bty Legend box type.
#' @param hilight Indices of points to highlight.
#' @param zlim Color scale limits.
#' @param zlog Log scale z axis? Logical.
#' @param softmax Apply softmax to z axis? Logical.
#' @param zsym Symmetrize z axis? Logical.
#' @param xlab,ylab Axis labels.
#' @param cmin,cmax Min and max color values.
#' @param xlim,ylim Axis limits.
#' @param hilight2 Secondary highlight indices.
#' @param hilight.col Color for highlighted points.
#' @param hilight.lwd Line width for highlighted points.
#' @param hilight.cex Size for highlighted points.
#' @param opacity Opacity factor, between 0 and 1.
#' @param label.clusters Logical, label clusters?
#' @param labels Point labels.
#' @param legend.ysp,legend.pos Legend position controls.
#' @param tooltip Tooltip text.
#' @param theme Theme parameters.
#' @param set.par Reset par? Logical.
#' @param label.type Label type, "text" or "box".
#' @param base_size Base point size.
#' @param title Plot title.
#' @param barscale Bar width scaling for box labels.
#' @param axis,box,guide Axis, box, and legend controls.
#'
#' @return A ggplot2 scatterplot object.
#'
#' @examples
#' \dontrun{
#' x <- stats::rnorm(100)
#' y <- stats::rnorm(100)
#' df <- data.frame(x = x, y = y)
#' p <- pgx.scatterPlotXY.GGPLOT(df)
#' }
#' @export
pgx.scatterPlotXY.GGPLOT <- function(pos, var = NULL, type = NULL, col = NULL, cex = NULL,
                                     cex.lab = 0.8, cex.title = 1.2, cex.clust = 1.5,
                                     cex.legend = 1, cex.axis = 1, gridcolor = NULL, bgcolor = NULL,
                                     zoom = 1, legend = TRUE, bty = "n", hilight = NULL,
                                     zlim = NULL, zlog = FALSE, softmax = FALSE, zsym = FALSE,
                                     xlab = NULL, ylab = NULL, cmin = 0, cmax = 1, xlim = NULL, ylim = NULL,
                                     hilight2 = hilight, hilight.col = "black",
                                     hilight.lwd = 0.8, hilight.cex = NULL,
                                     opacity = 1, label.clusters = FALSE, labels = NULL,
                                     legend.ysp = 0.85, legend.pos = "bottomleft",
                                     tooltip = NULL, theme = NULL, set.par = TRUE,
                                     label.type = c("text", "box"), base_size = 11,
                                     title = NULL, barscale = 0.8, axis = TRUE, box = TRUE, guide = "legend") {
  if (!is.null(var) && !is.null(ncol(var))) {
    var <- var[, 1]
  }
  if (!is.null(var) && length(var) == nrow(pos) && is.null(names(var))) {
    names(var) <- rownames(pos)
  }
  if (is.null(var)) {
    var <- rep("_", nrow(pos))
    names(var) <- rownames(pos)
  }
  var <- var[match(rownames(pos), names(var))]
  names(var) <- rownames(pos)

  if (is.null(type)) {
    type <- c("numeric", "factor")[1 + class(var) %in% c("factor", "character")]
  }
  label.type <- label.type[1]
  if (is.null(colnames(pos))) {
    colnames(pos) <- c("x", "y")
  }

  ## automatically set pointsize of dots
  if (is.null(cex)) {
    nr <- nrow(pos)
    i <- as.integer(cut(nr, breaks = c(0, 100, 500, 1000, 5000, Inf)))
    cex <- c(2, 1.4, 1, 0.7, 0.4)[i]
  }
  if (is.null(hilight.cex)) {
    hilight.cex <- cex
  }

  ## normalize pos
  if (is.null(xlim)) xlim <- range(pos[, 1])
  if (is.null(ylim)) ylim <- range(pos[, 2])
  if (zoom != 1) {
    cx <- mean(range(pos[, 1]))
    cy <- mean(range(pos[, 2]))
    dx <- diff(range(pos[, 1]))
    dy <- diff(range(pos[, 2]))
    xlim <- cx + 0.5 * c(-1, 1.05) * dx / zoom
    ylim <- cy + 0.5 * c(-1, 1.05) * dy / zoom
  }

  ax <- list(
    title = "",
    zeroline = FALSE,
    showline = FALSE,
    showticklabels = FALSE,
    showgrid = FALSE
  )

  legend.justification <- NULL
  legend.position <- legend.pos
  if (legend.pos == "bottomright") {
    legend.justification <- c(1, 0)
    legend.position <- c(1.0, 0.0)
  } else if (legend.pos == "topleft") {
    legend.justification <- c(0, 1)
    legend.position <- c(.0, 1.0)
  } else if (legend.pos == "topright") {
    legend.justification <- c(1, 1)
    legend.position <- c(1.0, 1.0)
  } else if (legend.pos == "bottomleft") {
    legend.justification <- c(0, 0)
    legend.position <- c(0.0, 0.0)
  }

  if (is.null(xlab)) xlab <- colnames(pos)[1]
  if (is.null(ylab)) ylab <- colnames(pos)[2]

  plt <- NULL
  ## Plot the discrete variables
  if (type == "factor") {
    z1 <- var
    if (!inherits(z1, "factor")) z1 <- factor(as.character(z1))
    nz <- length(levels(z1))
    col1 <- NULL
    if (!is.null(col)) {
      col1 <- col
    } else if (nz == 2) {
      col1 <- c("#eeeeee", "#f23451")
    } else if (nz == 1) {
      col1 <- "#3b4252"
    } else {
      col1 <- c(
        RColorBrewer::brewer.pal(8, "Set1"),
        RColorBrewer::brewer.pal(8, "Set2"),
        RColorBrewer::brewer.pal(12, "Set3")
      )
    }
    col1 <- Matrix::head(rep(col1, 99), nz)
    pt.col <- col1[z1]
    pt.col[is.na(pt.col)] <- "#d8d8d8"
    if (opacity < 1) {
      pt.col <- add_opacity(pt.col, opacity**0.33)
      col1 <- add_opacity(col1, opacity**0.33)
    }


    tt <- ifelse(!is.null(title), title, "var")
    tt <- "value"
    tooltip <- rownames(pos)
    tooltip <- paste(rownames(pos), "<br>", tt, "=", z1)
    label1 <- rownames(pos)
    if (!is.null(labels) && length(labels) == nrow(pos)) label1 <- labels

    df <- data.frame(
      x = pos[, 1], y = pos[, 2], name = rownames(pos),
      text = tooltip, label = label1
    )
    jj <- order(-table(pt.col)[pt.col]) ## plot less frequent points last...
    df <- df[jj, ]
    pt.col <- pt.col[jj]
    cex1 <- ifelse(length(cex) > 1, cex[jj], cex)
    x <- y <- NULL
    plt <- ggplot2::ggplot(df, ggplot2::aes(x, y), legend = legend) +
      ggplot2::geom_point(
        shape = 21,
        alpha = opacity,
        size = 1.8 * cex1,
        fill = pt.col,
        color = "#444444",
        stroke = 0.2
      ) +
      ggplot2::scale_color_manual(
        values = col1,
        name = title,
        na.value = "#DDDDDD44"
      )

    if (!is.null(theme)) {
      plt <- plt + theme
    } else {
      plt <- plt + ggplot2::theme_bw(
        base_size = base_size
      )
    }

    ## label cluster
    if (label.clusters) {
      name <- NULL
      ## Put a cluster label at the median position of the group
      mpos <- apply(pos, 2, function(x) tapply(x, z1, stats::median))
      mlab <- rownames(mpos)
      df1 <- data.frame(x = mpos[, 1], y = mpos[, 2], name = rownames(mpos))
      if (label.type == "text") labelFUN <- ggrepel::geom_text_repel
      if (label.type == "box") labelFUN <- ggrepel::geom_label_repel
      plt <- plt +
        labelFUN(
          data = df1,
          ggplot2::aes(x = x, y = y, label = name),
          size = 3.0 * cex.clust,
          color = "black",
          label.size = 0.10,
          max.overlaps = 99,
          fill = scales::alpha(c("white"), 0.7),
          segment.color = "grey30",
          segment.size = 0.2,
          box.padding = grid::unit(0.4, "lines"),
          point.padding = grid::unit(0.0, "lines")
        )
    }

    nlev <- length(levels(z1))
    if (legend && nlev <= 10) {
      plt <- plt +
        ggplot2::theme(
          legend.title = ggplot2::element_blank(),
          legend.text = ggplot2::element_text(size = 9 * cex.legend),
          legend.key.size = grid::unit(legend.ysp * 0.8 * cex.legend, "lines"),
          legend.key.height = grid::unit(legend.ysp * 0.8 * cex.legend, "lines"),
          legend.key = element_rect(color = "transparent", fill = scales::alpha("white", 0.0)),
          legend.justification = legend.justification,
          legend.position = legend.position,
          legend.background = ggplot2::element_rect(fill = scales::alpha("white", 0.5)),
          legend.margin = margin(0, 4, 4, 4),
          legend.box.just = "right",
          legend.box.background = ggplot2::element_rect(color = "#888888", size = 0.25),
          legend.box.margin = margin(0.8, 1, 1, 1)
        ) +
        ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 2.8 * cex.legend)))
    } else {
      plt <- plt + ggplot2::theme(legend.position = "none")
    }
  }

  ## Plot the continous variables
  if (type == "numeric") {
    z <- as.numeric(var)
    cpal <- rev(viridis::viridis(11))
    cpal <- rev(RColorBrewer::brewer.pal(11, "RdYlBu")) ## default
    if (!is.null(col)) {
      cpal <- col
    }
    if (opacity < 1) {
      cpal <- add_opacity(cpal, opacity**0.33)
    }


    tt <- ifelse(!is.null(title), title, "var")
    tt <- "value"
    tooltip <- paste(rownames(pos), "<br>", tt, "=", round(z, digits = 4))
    label1 <- rownames(pos)
    if (!is.null(labels) && length(labels) == nrow(pos)) label1 <- labels
    df <- data.frame(
      x = pos[, 1], y = pos[, 2], name = rownames(pos),
      variable = z, text = tooltip, label = label1
    )
    jj <- order(abs(z), na.last = FALSE)
    df <- df[jj, ] ## strongest last??
    cex1 <- ifelse(length(cex) > 1, cex[jj], cex)

    ## determine range for colorbar
    zr <- range(z)
    if (zsym && min(zr, na.rm = TRUE) < 0) zr <- c(-1, 1) * max(abs(zr), na.rm = TRUE)
    zz <- round(c(zr[1], zr[2]), digits = 2)
    variable <- NULL
    plt <- ggplot2::ggplot(df, ggplot2::aes(x, y, fill = variable)) +
      ggplot2::geom_point(
        shape = 21,
        alpha = opacity,
        size = 1.8 * cex1,
        color = "#444444",
        stroke = 0.2
      ) +
      ggplot2::scale_fill_gradientn(
        colors = cpal,
        limits = c(cmin, cmax),
        guide = guide,
        breaks = zz,
        labels = c(zz[1], zz[2]),
        na.value = "#DDDDDD44"
      ) +
      ggplot2::expand_limits(color = zr + c(-0.01, 0.01))


    ## colorscale bar
    if (legend) {
      xmax <- round(max(z, na.rm = TRUE), 2)
      xmin <- round(min(z, na.rm = TRUE), 2)
      plt <- plt +
        ggplot2::guides(
          colour = ggplot2::guide_colourbar(
            barwidth  = 0.5 * barscale,
            barheight = 2.2 * barscale
          )
        ) +
        ggplot2::theme(
          legend.title = ggplot2::element_blank(),
          legend.text = ggplot2::element_text(size = 9 * cex.legend),
          legend.justification = legend.justification,
          legend.position = legend.position,
          legend.background = ggplot2::element_blank(),
          legend.key = ggplot2::element_blank()
        )
    } else {
      plt <- plt + theme(legend.position = "none")
    }
  }

  if (!is.null(hilight)) {
    ## this hilights some points (with color and size)  at certain positions
    sel <- which(df$name %in% hilight)
    if (is.null(hilight.col)) hilight.col <- "transparent"
    plt <- plt +
      ggplot2::geom_point(
        data = df[sel, ],
        mapping = ggplot2::aes(x, y),
        size = 2.0 * hilight.cex,
        shape = 21,
        stroke = 0.5 * hilight.lwd,
        fill = hilight.col,
        color = "grey20"
      )
  }

  if (!is.null(hilight2)) {
    ## this put text labels at certain positions
    if (label.type == "text") labelFUN <- ggrepel::geom_text_repel
    if (label.type == "box") labelFUN <- ggrepel::geom_label_repel
    ## geom_text_repel(
    plt <- plt + labelFUN(
      data = subset(df, name %in% hilight2),
      ggplot2::aes(label = label),
      size = 3.0 * cex.lab,
      color = "black",
      label.size = 0.08,
      max.overlaps = 99,
      fill = scales::alpha(c("white"), 0.6),
      segment.color = "grey20",
      segment.size = 0.5,
      box.padding = grid::unit(0.25, "lines"),
      point.padding = grid::unit(0.2, "lines")
    )
  }


  if (!is.null(bgcolor)) {
    plt <- plt + ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = bgcolor)
    )
  }

  if (!is.null(gridcolor)) {
    plt <- plt + ggplot2::theme(
      panel.grid.major = ggplot2::element_line(
        size = 0.4,
        linetype = "solid",
        colour = gridcolor
      ),
      panel.grid.minor = ggplot2::element_line(
        size = 0.15,
        linetype = "solid",
        colour = gridcolor
      )
    )
  }

  if (box) {
    plt <- plt +
      ggplot2::theme(
        panel.border = ggplot2::element_rect(fill = NA, color = "grey20", size = 0.15)
      )
  }

  ## additional theme
  plt <- plt +
    ggplot2::xlim(xlim[1], xlim[2]) +
    ggplot2::ylim(ylim[1], ylim[2]) +
    ggplot2::xlab(xlab) + ggplot2::ylab(ylab) + ggplot2::ggtitle(title) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 22 * cex.title, hjust = 0, vjust = -1),
      axis.text.x = ggplot2::element_text(size = 18 * cex.axis),
      axis.text.y = ggplot2::element_text(size = 18 * cex.axis),
      axis.title.x = ggplot2::element_text(size = 22 * cex.axis, vjust = +2),
      axis.title.y = ggplot2::element_text(size = 22 * cex.axis, vjust = +0),
      plot.margin = ggplot2::margin(1, 1, 1, 1, "mm") ## ??
    )

  if (axis == FALSE) {
    plt <- plt +
      ggplot2::theme(
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank()
      )
  }

  plt
}


#' Interactive scatter plot using plotly
#'
#' @param pos Data frame containing x and y coordinates
#' @param var Grouping variable
#' @param type Variable type, either "continuous" or "factor"
#' @param col Colors
#' @param cex Point expansion factor
#' @param cex.lab Label size factor
#' @param cex.title Title size factor
#' @param cex.clust Cluster label size factor
#' @param cex.legend Legend text size factor
#' @param cex.axis Axis text size factor
#' @param xlab,ylab Axis labels
#' @param xlim,ylim Axis limits
#' @param axis Show axis? Logical
#' @param zoom Zoom factor
#' @param legend Show legend? Logical
#' @param bty Legend box type
#' @param hilight Indices of points to highlight
#' @param hilight2 Secondary highlight indices
#' @param hilight.col Color for highlighted points
#' @param hilight.cex Size for highlighted points
#' @param hilight.lwd Line width for highlighted points
#' @param zlim Color scale limits
#' @param zlog Log scale z axis? Logical
#' @param zsym Symmetrize z axis? Logical
#' @param softmax Apply softmax to z axis? Logical
#' @param opacity Opacity factor, between 0 and 1
#' @param bgcolor Background color
#' @param box Draw box around plot? Logical
#' @param label.clusters Logical, label clusters?
#' @param labels Point labels
#' @param label.type Label type, "text" or "box"
#' @param tooltip Tooltip text
#' @param theme Theme parameters
#' @param set.par Reset par? Logical
#' @param title Plot title
#' @param title.y Title position
#' @param gridcolor Grid line color
#' @param source Character label for plot source
#' @param key Color key specification
#' @param displayModeBar Show plotly modebar? Logical
#'
#' @return A plotly interactive scatterplot object
#'
#' @export
pgx.scatterPlotXY.PLOTLY <- function(pos,
                                     var = NULL, type = NULL, col = NULL,
                                     cex = NULL, cex.lab = 0.8, cex.title = 1.2,
                                     cex.clust = 1.5, cex.legend = 1, cex.axis = 1,
                                     xlab = NULL, ylab = NULL, xlim = NULL, ylim = NULL,
                                     axis = TRUE, zoom = 1, legend = TRUE, bty = "n",
                                     hilight = NULL, hilight2 = hilight, hilight.col = NULL,
                                     hilight.cex = NULL, hilight.lwd = 0.8,
                                     zlim = NULL, zlog = FALSE, zsym = FALSE, softmax = FALSE,
                                     opacity = 1, bgcolor = NULL, box = TRUE,
                                     label.clusters = FALSE, labels = NULL, label.type = NULL,
                                     tooltip = NULL, theme = NULL, set.par = TRUE,
                                     title = "", title.y = 1, gridcolor = NULL,
                                     source = NULL, key = NULL,
                                     displayModeBar = FALSE) {
  if (!is.null(var) && NCOL(var) > 1) {
    var <- stats::setNames(var[, 1], rownames(var))
  }
  if (!is.null(var) && length(var) == nrow(pos) && is.null(names(var))) {
    names(var) <- rownames(pos)
  }
  if (is.null(var)) {
    var <- rep("_", nrow(pos))
    names(var) <- rownames(pos)
  }

  var <- var[match(rownames(pos), names(var))]
  names(var) <- rownames(pos)

  if (is.null(type)) {
    type <- c("numeric", "factor")[1 + class(var) %in% c("factor", "character")]
  }

  ## automatically set pointsize of dots
  if (is.null(cex)) {
    nr <- nrow(pos)
    i <- as.integer(cut(nr, breaks = c(0, 100, 500, 1000, 5000, Inf)))
    cex <- c(2, 1.4, 1, 0.7, 0.4)[i]
  }
  if (is.null(colnames(pos))) {
    colnames(pos) <- c("x", "y")
  }
  if (is.null(hilight.cex)) hilight.cex <- cex

  ## normalize pos
  xlim0 <- range(pos[, 1])
  ylim0 <- range(pos[, 2])
  if (zoom != 1) {
    cx <- mean(range(pos[, 1]))
    cy <- mean(range(pos[, 2]))
    dx <- diff(range(pos[, 1]))
    dy <- diff(range(pos[, 2]))
    xlim0 <- cx + 0.5 * c(-1, 1.05) * dx / zoom
    ylim0 <- cy + 0.5 * c(-1, 1.05) * dy / zoom
  }
  if (is.null(xlab)) xlab <- colnames(pos)[1]
  if (is.null(ylab)) ylab <- colnames(pos)[2]

  if (type == "numeric") var <- round(var, digits = 4)
  hoverinfo <- "text"
  tooltip1 <- paste0(
    rownames(pos),
    "<br>value = ", var,
    "<br>x = ", round(pos[, 1], 2), "; y = ", round(pos[, 2], 2)
  )
  if (!is.null(tooltip) && length(tooltip) == length(var)) {
    tooltip1 <- paste0(tooltip1, "<br>", tooltip)
  }
  if (!is.null(tooltip) && tooltip == FALSE) {
    tooltip1 <- NA
    hoverinfo <- "none"
  }

  label1 <- rownames(pos)
  if (!is.null(labels) && labels[1] != FALSE) label1 <- labels
  label1 <- gsub("[)(]", "", label1) ## HTML labels do not like it...
  label1 <- gsub("[_]", " ", label1) ## HTML labels do not like it...
  label1 <- sapply(label1, function(s) paste(strwrap(s, 50), collapse = "<br>"))

  plt <- NULL
  ## Plot the discrete variables
  if (type == "factor") {
    z1 <- factor(var)
    nz <- length(levels(z1))
    cpal <- NULL
    if (!is.null(col)) {
      cpal <- col
    } else if (nz == 2) {
      cpal <- rev(grDevices::grey.colors(2, end = 0.8))
      cpal <- c("#AAAAAA55", "#555555FF")
      cpal <- c("#00008855", "#AA0000FF")
      cpal <- c("#CCCCCC55", "#AA0000FF")
    } else if (nz == 1) {
      cpal <- c("#22222255")
    } else {
      cpal <- c(
        RColorBrewer::brewer.pal(8, "Set1"),
        RColorBrewer::brewer.pal(8, "Set2"),
        RColorBrewer::brewer.pal(12, "Set3")
      )
    }
    cpal <- Matrix::head(rep(cpal, 99), nz)

    if (opacity < 1) {
      cpal <- add_opacity(cpal, opacity**0.33)
    }

    ## prepare dataframe
    df <- data.frame(
      x = pos[, 1],
      y = pos[, 2],
      name = rownames(pos),
      value = z1,
      size = 7 * cex,
      text = tooltip1,
      label = label1
    )

    ## plot less frequent points first... (NEED RETHINK)
    jj <- order(-table(z1)[z1])
    df <- df[jj, ]
  }

  ## Plot the continous variables
  if (type == "numeric") {
    z <- as.numeric(var)
    z1 <- NULL
    if (is.null(col)) {
      cpal <- rev(RColorBrewer::brewer.pal(11, "RdYlBu"))
    } else {
      cpal <- col
    }
    if (opacity < 1) {
      cpal <- add_opacity(cpal, opacity**0.33)
    }


    df <- data.frame(
      x = pos[, 1],
      y = pos[, 2],
      name = rownames(pos),
      value = z,
      size = 7 * cex,
      text = tooltip1,
      label = label1
    )
    rownames(df) <- rownames(pos)

    ## plot low values first
    df <- df[order(abs(z)), , drop = FALSE]

    cmin <- min(df$value, na.rm = TRUE)
    cmax <- max(df$value, na.rm = TRUE)
    if (zsym) {
      cmax <- max(abs(df$value), na.rm = TRUE)
      cmin <- -cmax
    }
  }


  ## ---------------- call PLOTLY -----------
  if (is.null(source)) source <- paste0(sample(LETTERS, 10), collapse = "")
  ## plt <- plotly::plot_ly(df,
  plt <- plotly::plot_ly(
    source = source,
    showlegend = FALSE
  )

  any(is.na(df$value))
  if (any(is.na(df$value))) {
    jj <- which(is.na(df$value))
    plt <- plt %>%
      plotly::add_markers(
        data = df[jj, , drop = FALSE],
        x = ~x, y = ~y,
        colors = cpal,
        text = ~text,
        hoverinfo = hoverinfo,
        marker = list(
          size = ~size,
          opacity = opacity,
          color = "#DDDDDD44"
        ),
        showlegend = FALSE,
        key = ~name,
        mode = "markers",
        type = "scattergl"
      )
  }

  jj <- which(!is.na(df$value))
  plt <- plt %>%
    plotly::add_markers(
      data = df[jj, , drop = FALSE],
      x = ~x, y = ~y,
      color = ~value,
      colors = cpal,
      text = ~text,
      hoverinfo = hoverinfo,
      marker = list(
        size = ~size,
        opacity = opacity,
        line = list(
          color = "#444444",
          width = 0.2
        )
      ),
      showlegend = FALSE,
      key = ~name,
      mode = "markers",
      type = "scattergl"
    )

  if (!is.null(hilight)) {
    jj <- which(rownames(df) %in% hilight)
    col1 <- "transparent"
    if (!is.null(hilight.col)) col1 <- hilight.col
    plt <- plt %>%
      ## add_trace(
      plotly::add_markers(
        data = df[jj, ],
        x = ~x, y = ~y,
        color = ~value, colors = cpal,
        color = NULL,
        key = ~name,
        mode = "markers", type = "scattergl", #
        text = ~text,
        hoverinfo = hoverinfo,
        marker = list(
          color = col1,
          size = 5 * hilight.cex,
          showlegend = FALSE,
          showscale = FALSE,
          line = list(
            color = "#000000",
            width = 1.0 * hilight.lwd
          )
        )
      )
  }

  if (!is.null(hilight2)) {
    jj <- which(rownames(df) %in% hilight2)
    plt <- plt %>%
      plotly::add_annotations(
        data = df[jj, , drop = FALSE],
        x = ~x,
        y = ~y,
        text = ~label,
        yanchor = "bottom",
        xanchor = "center", ## left,center,right
        showarrow = FALSE,
        showlegend = FALSE,
        font = list(size = 12 * cex.lab),
        xref = "x", yref = "y"
      )
  }
  ## label cluster
  if (label.clusters) {
    mpos <- apply(pos, 2, function(x) tapply(x, z1, stats::median))
    # If there is only one cluster
    if (length(unique(z1)) == 1) {
      mpos <- data.frame(mpos[1], mpos[2])
      mlab <- unique(z1)
    } else {
      mlab <- rownames(mpos)
    }

    plt <- plt %>%
      plotly::add_annotations(
        x = mpos[, 1], y = mpos[, 2],
        showarrow = FALSE, text = mlab,
        font = list(size = 15 * cex.clust),
        xref = "x", yref = "y"
      )
  }

  if (!legend) {
    plt <- plt %>%
      plotly::hide_colorbar() %>%
      plotly::hide_legend() %>%
      plotly::layout(showlegend = FALSE)
  }

  if (legend && type == "numeric") {
    plt <- plt %>%
      plotly::colorbar(
        limits = c(cmin, cmax),
        len = 0.15,
        thickness = 9,
        x = 0.01,
        y = 0.1,
        title = "",
        tickfont = list(size = 9)
      )
  }

  if (axis == FALSE) {
    no.axis <- list(
      title = "",
      zeroline = FALSE,
      showline = FALSE,
      showticklabels = FALSE
    )
    plt <- plt %>% plotly::layout(
      xaxis = no.axis,
      yaxis = no.axis
    )
  }

  if (!is.null(gridcolor)) {
    plt <- plt %>% plotly::layout(
      xaxis = list(gridcolor = gridcolor),
      yaxis = list(gridcolor = gridcolor)
    )
  }

  if (!is.null(bgcolor)) {
    plt <- plt %>% plotly::layout(
      plot_bgcolor = bgcolor
    )
  }

  if (box == TRUE) {
    plt <- plt %>%
      plotly::layout(
        shapes = list(
          list(
            type = "rect",
            line = list(
              color = "#666",
              width = 0.1
            ),
            xref = "paper",
            yref = "paper",
            y0 = -0.0, y1 = 1.0,
            x0 = -0.0, x1 = 1.0
          )
        )
      )
  }


  ## add legend and title
  plt <- plt %>%
    plotly::layout(
      showlegend = legend,
      xaxis = list(
        title = xlab,
        titlefont = list(size = 12 * cex.axis)
      ),
      yaxis = list(
        title = ylab,
        titlefont = list(size = 12 * cex.axis)
      )
    )

  if (!is.null(title) && title != "") {
    plt <- plt %>%
      plotly::layout(
        annotations = list(
          text = title,
          font = list(size = 14 * cex.title),
          xref = "paper",
          yref = "paper",
          yanchor = "bottom",
          xanchor = "left",
          align = "right",
          x = 0,
          y = title.y,
          showarrow = FALSE
        ),
        ## add top margin??
        margin = list(
          l = 10,
          r = 10,
          b = 10,
          t = 10
        )
      )
  } else {
    ## symmetric margin
    plt <- plt %>%
      plotly::layout(
        margin = list(l = 10, r = 10, b = 10, t = 10)
      )
  }
  plt
}


#' @title Scatter Plot XY using scatterD3
#'
#' @description This function creates a scatter plot of two variables
#' using the scatterD3 library for interactive scatter plots.
#'
#' @param pos A data frame or matrix containing the x and y coordinates of the points to be plotted.
#' @param var An optional numeric vector specifying the value of a third variable to be used for coloring the points.
#' @param type Not used.
#' @param col Not used.
#' @param cex A numeric value specifying the size of the points.
#' @param cex.lab Not used.
#' @param cex.title Not used.
#' @param cex.clust Not used.
#' @param cex.legend Not used.
#' @param zoom Not used.
#' @param legend Not used.
#' @param bty Not used.
#' @param hilight Not used.
#' @param zlim Not used.
#' @param zlog Not used.
#' @param softmax Not used.
#' @param xlab A character string specifying the label for the x-axis.
#' @param ylab A character string specifying the label for the y-axis.
#' @param hilight2 Not used.
#' @param opacity Not used.
#' @param label.clusters Not used.
#' @param labels Not used.
#' @param legend.ysp Not used.
#' @param legend.pos Not used.
#' @param tooltip Not used.
#' @param theme Not used.
#' @param set.par Not used.
#' @param title Not used.
#' @param barscale Not used.
#'
#' @export
pgx.scatterPlotXY.D3 <- function(pos, var = NULL, type = NULL, col = NULL, cex = 1,
                                 cex.lab = 0.8, cex.title = 1.2, cex.clust = 1.5, cex.legend = 1,
                                 zoom = 1, legend = TRUE, bty = "n", hilight = NULL,
                                 zlim = NULL, zlog = FALSE, softmax = FALSE,
                                 xlab = NULL, ylab = NULL, hilight2 = hilight,
                                 opacity = 1, label.clusters = FALSE, labels = NULL,
                                 legend.ysp = 0.85, legend.pos = "bottomleft",
                                 tooltip = NULL, theme = NULL, set.par = TRUE,
                                 title = NULL, barscale = 0.8) {
  if (is.null(colnames(pos))) {
    colnames(pos) <- c("x", "y")
  }
  if (is.null(xlab)) {
    xlab <- colnames(pos)[1]
  }
  if (is.null(ylab)) {
    ylab <- colnames(pos)[2]
  }
  x <- y <- z <- NULL # Init vars to prevent warning
  df <- data.frame(x = pos[, 1], y = pos[, 2], z = var, names = rownames(pos))
  if (!is.null(var)) {
    plt <- scatterD3::scatterD3(
      data = df, x = x, y = y,
      col_var = z,
      xlab = xlab, ylab = ylab,
      point_size = 32 * cex,
      legend_width = 70,
      col_lab = "value"
    )
  } else {
    plt <- scatterD3::scatterD3(
      data = df, x = x, y = y, #
      xlab = xlab, ylab = ylab,
      point_size = 32 * cex,
      col_lab = "value"
    )
  }
  plt
}




#' Stacked barplot
#'
#' @param x Numeric data matrix. Rows are stacked, columns are groups.
#' @param showlegend Logical indicating if legend should be displayed.
#' @param ylab Character string for y-axis label. Default NULL.
#' @param xlab Character string for x-axis label. Default NULL.
#' @param horiz Logical indicating if bars should be horizontal. Default FALSE.
#'
#' @return A plotly stacked barplot object.
#'
#' @details This function takes a numeric data matrix and generates an interactive
#' stacked barplot using plotly. Rows of the input matrix are stacked, and columns
#' represent groups. Set horiz=TRUE to generate horizontal stacked bars.
#'
#' @examples
#' \dontrun{
#' x <- matrix(rnorm(100), ncol = 5)
#' pgx.stackedBarplot(x)
#' }
#' @export
pgx.stackedBarplot <- function(x,
                               showlegend,
                               ylab = NULL,
                               xlab = NULL,
                               horiz = FALSE) {
  x_plot <- cbind(data.frame(groups = rownames(x)), x)
  x_plot <- data.table::melt(x_plot, id.vars = "groups", value.name = "Effect")

  if (horiz == FALSE) {
    x_plot$groups <- factor(x_plot$groups, levels = rownames(x))
  } else {
    c1 <- which(colnames(x_plot) == "variable")
    c2 <- which(colnames(x_plot) == "Effect")
    c3 <- which(colnames(x_plot) == "groups")
    colnames(x_plot)[c1] <- "Effect"
    colnames(x_plot)[c2] <- "groups"
    colnames(x_plot)[c3] <- "variable"
  }

  plotly::plot_ly(
    x_plot,
    x = ~groups,
    y = ~Effect,
    type = "bar",
    name = ~variable,
    color = ~variable
  ) %>%
    plotly::layout(
      showlegend = showlegend,
      barmode = "stack",
      yaxis = list(title = ylab),
      xaxis = list(title = xlab)
    )
}


#' @title Dark Mode Theme
#'
#' @description
#' Changes the theme of a plotly plot to a dark mode style.
#'
#' @param p A plotly plot object
#' @param dim The number of dimensions for the plot (2 or 3)
#'
#' @details
#' This function modifies a plotly plot to use a dark theme style.
#' It changes the colors of the background, axes, text, etc to darker colors.
#'
#' The \code{p} argument is a plotly plot object.
#' The \code{dim} argument specifies whether it is a 2D or 3D plot.
#'
#' It sets the plot and paper background colors to dark greys.
#' The font and axis colors are changed to light greys.
#' For 3D plots the z-axis is also styled.
#'
#' @return
#' The modified plotly plot object with a dark theme.
#'
#' @export
darkmode <- function(p, dim = 2) {
  font.par <- list(
    color = "#AAA"
  )
  axis.par <- list(
    color = "#AAA",
    linecolor = "#AAA",
    tickcolor = "#AAA",
    zerolinecolor = "#AAA"
  )
  p <- plotly::layout(p,
    plot_bgcolor = "rgb(10,20,40)",
    paper_bgcolor = "rgb(10,20,40)",
    xaxis = axis.par,
    yaxis = axis.par,
    title = font.par
  )
  if (dim == 3) {
    p <- plotly::layout(p,
      zaxis = axis.par
    )
  }
  return(p)
}

#' Interactive MA plot using plotly
#'
#' @param x Numeric vector of average expression values (log2 CPM)
#' @param y Numeric vector of effect sizes (log2 fold change)
#' @param names Vector of gene or feature names
#' @param source Character label for plot source. Default "plot1".
#' @param group.names Character vector of group names for legend.
#' @param xlab Label for x-axis. Default "average expression (log2.CPM)".
#' @param ylab Label for y-axis. Default "effect size (log2.FC)".
#' @param lfc Cutoff for absolute log fold change. Default 1.
#' @param psig Cutoff for adjusted p-value. Default 0.05.
#' @param showlegend Show legend? Logical. Default TRUE.
#' @param highlight Vector of genes to highlight. Default NULL.
#' @param marker.size Marker size. Default 5.
#' @param label Vector of labels for highlighted genes. Default NULL.
#' @param label.cex Text size for labels. Default 1.
#' @param marker.type Marker type (scatter, line, etc). Default "scatter".
#' @param displayModeBar Show plotly modebar? Logical. Default TRUE.
#'
#' @return A plotly interactive MA plot object.
#'
#' @export
plotlyMA <- function(x, y, names, source = "plot1",
                     group.names = c("group1", "group2"),
                     xlab = "average expression (log2.CPM)",
                     ylab = "effect size (log2.FC)",
                     lfc = 1, psig = 0.05, showlegend = TRUE, highlight = NULL,
                     marker.size = 5, label = NULL, label.cex = 1,
                     marker.type = "scatter", displayModeBar = TRUE) {
  if (is.null(highlight)) highlight <- names
  i0 <- which(!names %in% highlight)
  i1 <- which(names %in% highlight)

  p <- plotly::plot_ly(
    type = marker.type, mode = "markers"
  )

  p <- p %>%
    plotly::event_register("plotly_hover") %>%
    plotly::event_register("plotly_selected")

  if (length(i0)) {
    p <- p %>%
      plotly::add_trace(
        x = x[i0],
        y = y[i0],
        text = names[i0],
        marker = list(
          size = marker.size,
          color = "#ccc"
        ),
        showlegend = showlegend
      )
  }

  if (length(i1)) {
    p <- p %>%
      plotly::add_trace(
        x = x[i1],
        y = y[i1],
        text = names[i1],
        marker = list(
          size = marker.size,
          color = "#1f77b4"
        ),
        showlegend = showlegend
      )
  }

  if (!is.null(label) && length(label) > 0) {
    i2 <- which(names %in% label)
    p <- p %>%
      plotly::add_annotations(
        x = x[i2],
        y = y[i2],
        text = names[i2],
        font = list(
          size = 12 * label.cex,
          color = "#1f77b4"
        ),
        showarrow = FALSE,
        yanchor = "bottom",
        yshift = 2,
        textposition = "top"
      )
  }

  x1 <- 1.05 * max(x)
  yy <- 1.05 * max(abs(y))
  abline1 <- list(
    type = "line", y0 = -lfc, y1 = -lfc, x0 = 0, x1 = x1,
    line = list(dash = "dot", width = 1, color = "grey")
  )
  abline2 <- list(
    type = "line", y0 = +lfc, y1 = +lfc, x0 = 0, x1 = x1,
    line = list(dash = "dot", width = 1, color = "grey")
  )

  xrange <- c(0, 1) * max(abs(x)) * 1.05
  xrange <- range(x, na.rm = TRUE)
  yrange <- c(-1, 1) * max(abs(y), na.rm = TRUE) * 1.05
  xaxis <- list(title = xlab, range = xrange, gridwidth = 0.2, showgrid = FALSE)
  yaxis <- list(title = ylab, range = yrange, gridwidth = 0.2, showgrid = FALSE)

  p <- p %>%
    plotly::layout(
      shapes = list(abline1, abline2),
      xaxis = xaxis, yaxis = yaxis,
      showlegend = showlegend,
      hovermode = "closest",
      dragmode = "select"
    ) %>%
    plotly::config(displayModeBar = displayModeBar)

  xann <- c(1, 1) * 0.95
  yann <- c(0, 0.99)
  ann.text <- paste("UP in", group.names[c(2, 1)])

  p <- p %>%
    plotly::add_annotations(
      x = xann,
      y = yann,
      text = ann.text,
      font = list(size = 10),
      xanchor = c("right", "right"),
      align = c("right", "right"),
      showarrow = FALSE,
      xref = "paper",
      yref = "paper",
      borderpad = 3,
      bordercolor = "black",
      borderwidth = 0.6
    )

  p
}


#' Interactive volcano plot using plotly
#'
#' @param x Numeric vector of effect sizes (log fold changes)
#' @param y Numeric vector of p-values (-log10 transformed)
#' @param names Vector of gene or feature names
#' @param source Character label for plot source. Default "plot1".
#' @param group.names Character vector of group names for legend.
#' @param xlab Label for x-axis. Default "effect size (logFC)".
#' @param ylab Label for y-axis. Default "significance (-log10p)".
#' @param lfc Cutoff for absolute log fold change. Default 1.
#' @param psig Cutoff for adjusted p-value. Default 0.05.
#' @param showlegend Show legend? Logical. Default TRUE.
#' @param highlight Vector of genes to highlight. Default NULL.
#' @param marker.size Marker size. Default 5.
#' @param label Vector of labels for highlighted genes. Default NULL.
#' @param label.cex Text size for labels. Default 1.
#' @param marker.type Marker type (scatter, line, etc). Default "scatter".
#' @param displayModeBar Show plotly modebar? Logical. Default TRUE.
#'
#' @return A plotly interactive volcano plot object.
#'
#' @export
plotlyVolcano <- function(x, y, names, source = "plot1", group.names = c("group1", "group2"),
                          xlab = "effect size (logFC)", ylab = "significance (-log10p)",
                          lfc = 1, psig = 0.05, showlegend = TRUE, highlight = NULL,
                          marker.size = 5, label = NULL, label.cex = 1,
                          marker.type = "scatter", displayModeBar = TRUE) {
  if (is.null(highlight)) highlight <- names
  i0 <- which(!names %in% highlight)
  i1 <- which(names %in% highlight)
  p <- plotly::plot_ly(
    type = marker.type, mode = "markers"
    ## source=source, key=1:length(x)
  )
  p <- p %>%
    plotly::event_register("plotly_hover") %>%
    plotly::event_register("plotly_selected")
  if (length(i0)) {
    p <- p %>%
      plotly::add_trace(
        x = x[i0],
        y = y[i0],
        text = names[i0],
        marker = list(
          size = marker.size,
          color = "#ccc"
        ),
        showlegend = showlegend
      )
  }
  if (length(i1)) {
    p <- p %>%
      plotly::add_trace(
        x = x[i1],
        y = y[i1],
        text = names[i1],
        marker = list(
          size = 5,
          color = "#1f77b4"
        ),
        showlegend = showlegend
      )
  }
  if (!is.null(label) && length(label) > 0) {
    i2 <- which(names %in% label)
    p <- p %>%
      plotly::add_annotations(
        x = x[i2],
        y = y[i2],
        text = names[i2],
        font = list(
          size = 12 * label.cex,
          color = "#1f77b4"
        ),
        showarrow = FALSE,
        yanchor = "bottom",
        yshift = 2,
        textposition = "top"
      )
  }
  y0 <- -log10(psig)
  y1 <- 1.05 * max(y)
  xx <- 1.05 * max(abs(x))
  abline1 <- list(
    type = "line", x0 = -lfc, x1 = -lfc, y0 = 0, y1 = y1,
    line = list(dash = "dot", width = 1, color = "grey")
  )
  abline2 <- list(
    type = "line", x0 = +lfc, x1 = +lfc, y0 = 0, y1 = y1,
    line = list(dash = "dot", width = 1, color = "grey")
  )
  abline3 <- list(
    type = "line", x0 = -xx, x1 = +xx, y0 = y0, y1 = y0,
    line = list(dash = "dot", width = 1, color = "grey")
  )
  max.absx <- max(max(abs(x), na.rm = TRUE), lfc * 1.2)
  max.absy <- max(max(abs(y), na.rm = TRUE), y0 * 1.2)
  xrange <- c(-1, 1) * max.absx * 1.05
  if (min(x) >= 0) xrange <- c(0, 1) * max.absx * 1.05
  yrange <- c(0, 1) * max.absy * 1.05
  xaxis <- list(title = xlab, range = xrange, showgrid = FALSE) # titlefont = list(size = 12))
  yaxis <- list(title = list(text = ylab, standoff = 20L), range = yrange, showgrid = FALSE) # titlefont = list(size = 12))

  p <- p %>%
    plotly::layout(
      shapes = list(abline1, abline2, abline3),
      xaxis = xaxis, yaxis = yaxis,
      showlegend = showlegend,
      hovermode = "closest",
      dragmode = "select"
    ) %>%
    plotly_modal_default() %>%
    plotly::layout(
      margin = list(l = 0, b = 1, t = 10, r = 10),
      font = list(size = 12),
      legend = list(
        font = list(size = 12)
      )
    )
  p
}


#' Interactive cytoPlot using plotly
#'
#' @param pgx Expression data matrix
#' @param gene1 First gene for x-axis
#' @param gene2 Second gene for y-axis
#' @param samples Samples to include
#' @param nbinsx Number of bins for x-axis
#' @param nbinsy Number of bins for y-axis
#' @param lab.unit Units for axis labels. Default "(log2CPM)".
#' @param reversescale Reverse color scale? Default TRUE.
#' @param marker.size Size of scatter markers. Default 5.
#' @param contour.coloring Add density shading? "none", "color", or
#' "black". Default "none".
#' @param marker.color Color of markers. Default "black".
#' @param showgrid Show grid lines? Logical. Default TRUE.
#'
#' @return A plotly interactive cytoPlot object
#'
#' @details This function generates an interactive cytoPlot using the plotly package.
#' Expression values for two genes are binned and plotted with one gene on each axis.
#' Samples to include and graphical parameters like color scale, marker size, etc. can be adjusted.
#'
#' @export
plotlyCytoplot <- function(pgx,
                           gene1,
                           gene2,
                           samples,
                           nbinsx,
                           nbinsy,
                           lab.unit = "(log2CPM)",
                           reversescale = TRUE,
                           marker.size = 5,
                           contour.coloring = "none",
                           marker.color = "black",
                           showgrid = TRUE) {
  samples <- samples
  if (is.null(samples)) {
    samples <- colnames(pgx$X)
  }

  samples <- intersect(samples, colnames(pgx$X))

  x1 <- pgx$X[gene1, samples]
  x2 <- pgx$X[gene2, samples]

  names(x1) <- samples
  names(x2) <- samples

  m1 <- mean(x1)
  m2 <- mean(x2)

  ## select samples in different quadrants
  j1 <- length(samples[which(x1 < m1 & x2 > m2)])
  j2 <- length(samples[which(x1 > m1 & x2 < m2)])
  j3 <- length(samples[which(x1 > m1 & x2 > m2)])
  j4 <- length(samples[which(x1 < m1 & x2 < m2)])

  xlab1 <- paste(gene1, lab.unit, collapse = "  ")
  ylab1 <- paste(gene2, lab.unit, collapse = "  ")

  xaxis <- list(title = xlab1, range = range(x1), gridwidth = 0.2, showgrid = TRUE, showline = TRUE, autorange = TRUE)
  yaxis <- list(title = ylab1, range = range(x2), gridwidth = 0.2, showgrid = TRUE, showline = TRUE, autorange = TRUE)


  p <- plotly::plot_ly(
    x = x1,
    y = x2,
    marker = list(size = marker.size),
    showlegend = TRUE
  ) %>%
    plotly::add_trace(
      x = x1,
      y = x2,
      type = "histogram2dcontour",
      contours = list(coloring = contour.coloring),
      nbinsx = nbinsx,
      nbinsy = nbinsy
    ) %>%
    plotly::add_trace(
      x = x1,
      y = x2,
      text = names(x1),
      hoverinfo = "text",
      marker = list(size = marker.size, color = "black")
    ) %>%
    plotly::layout(
      shapes = list(list(
        type = "line",
        x0 = 0,
        x1 = 1,
        xref = "paper",
        y0 = mean(x1),
        y1 = mean(x1),
        line = list(color = "#bebebe", dash = "dot")
      ), list(
        type = "line",
        y0 = 0,
        y1 = 1,
        yref = "paper",
        x0 = mean(x2),
        x1 = mean(x2),
        line = list(color = "#bebebe", dash = "dot")
      )),
      xaxis = xaxis,
      yaxis = yaxis
    )


  N <- length(x1)

  quadrants <- c(j3, j1, j2, j4)

  positions <- matrix(c(1, 1, 0, 1, 1, 0, 0, 0), ncol = 2, byrow = TRUE)

  for (i in 1:4) {
    p <- p %>% plotly::add_annotations(
      x = positions[, 1][i], y = positions[, 2][i],
      text = paste(round(100 * quadrants[i] / N, 2), "%"),
      showarrow = FALSE,
      font = list(size = 15),
      xref = "paper",
      yref = "paper",
      showlegend = FALSE
    )
  }

  if (!is.null(pgx$deconv) && length(pgx$deconv) > 0) {
    inferred.celltype <- pgx$deconv[[1]][["meta"]]
    lab1 <- Matrix::head(names(sort(-Matrix::colSums(inferred.celltype[j1, , drop = FALSE]))), 3)
    pos1 <- apply(cbind(x1, x2)[j1, , drop = FALSE], 2, stats::median)
    p <- p %>% plotly::add_annotations(
      x = pos1[1], y = pos1[2],
      text = paste(lab1, collapse = "\n"),
      showarrow = FALSE,
      font = list(size = 15)
    )

    lab2 <- Matrix::head(names(sort(-Matrix::colSums(inferred.celltype[j2, , drop = FALSE]))), 3)
    pos2 <- apply(cbind(x1, x2)[j2, , drop = FALSE], 2, stats::median)
    p <- p %>% plotly::add_annotations(
      x = pos2[1], y = pos2[2],
      text = paste(lab2, collapse = "\n"),
      showarrow = FALSE,
      font = list(size = 15)
    )

    lab3 <- Matrix::head(names(sort(-Matrix::colSums(inferred.celltype[j3, , drop = FALSE]))), 3)
    pos3 <- apply(cbind(x1, x2)[j3, , drop = FALSE], 2, stats::median)
    p <- p %>% plotly::add_annotations(
      x = pos3[1], y = pos3[2],
      text = paste(lab3, collapse = "\n"),
      showarrow = FALSE,
      font = list(size = 15)
    )
  }
  return(p)
}


#' Correlation-based Clustering
#'
#' @title Correlation-based Clustering
#'
#' @description This function performs hierarchical clustering on a
#' dataset using a correlation-based distance metric.
#'
#' @param x A numeric matrix or data frame containing the data to be clustered.
#'
#' @details The function first calculates the pairwise correlation between
#' the columns of the input data using the `cor` function from the `stats`
#' package. The correlation matrix is then transformed into a distance matrix by
#' subtracting it from 1. The resulting distance matrix is used as input to the
#' `hclust` function from the `fastcluster` package, which performs hierarchical
#' clustering using Ward's method.
#'
#' @return An object of class `hclust` representing the hierarchical clustering
#' of the input data.
#'
#' @examples
#' \dontrun{
#' # Generate some example data
#' set.seed(123)
#' x <- matrix(rnorm(1000), ncol = 10)
#'
#' # Perform correlation-based clustering
#' hc <- corclust(x)
#'
#' # Plot the resulting dendrogram
#' plot(hc)
#' }
#'
#' @export
corclust <- function(x) {
  dd <- stats::as.dist(1 - stats::cor(t(x), use = "pairwise"))
  hc <- fastcluster::hclust(dd, method = "ward.D2")
  hc
}

## Override add_col_annotation to be able to suppress titles
##
##
# setMethod(add_col_annotation,

iheatmapr.add_col_annotation <- function(p,
                                         annotation,
                                         colors = NULL,
                                         side = c("top", "bottom"),
                                         size = 0.05,
                                         buffer = 0.015,
                                         inner_buffer = buffer / 2,
                                         show_title = TRUE,
                                         show_legend = TRUE,
                                         layout = list()) {
  side <- match.arg(side)
  # Convert to data.frame
  x <- as.data.frame(annotation)

  for (i in seq_len(ncol(x))) {
    if (is.character(x[, i]) || is.factor(x[, i]) || is.logical(x[, i])) {
      if (!is.null(colors) && colnames(x)[i] %in% names(colors)) {
        tmp_colors <- colors[[colnames(x)[i]]]
      } else {
        tmp_colors <- iheatmapr:::pick_discrete_colors(as.factor(x[, i]), p)
      }
      p <- iheatmapr::add_col_groups(p,
        x[, i],
        name = colnames(x)[i],
        title = colnames(x)[i],
        colors = tmp_colors,
        show_colorbar = show_legend,
        side = side,
        size = size,
        buffer = if (i == 1) {
          buffer
        } else {
          inner_buffer
        },
        layout = layout,
        show_title = show_title
      )
    } else if (is.numeric(x[, i])) {
      if (!is.null(colors) && colnames(x)[i] %in% names(colors)) {
        tmp_colors <- colors[[colnames(x)[i]]]
      } else {
        tmp_colors <- pick_continuous_colors(
          zmid = 0,
          zmin = min(x[, i]),
          zmax = max(x[, i]), p
        )
      }
      p <- iheatmapr::add_col_signal(p,
        x[, i],
        name = colnames(x)[i],
        colors = tmp_colors,
        side = side,
        size = size,
        buffer = if (i == 1) {
          buffer
        } else {
          inner_buffer
        },
        layout = layout,
        show_title = show_title
      )
    } else {
      stop("Input should be character, factor, logical, or numeric")
    }
  }
  return(p)
}


#' Split heatmap from matrix
#'
#' @param X Numeric data matrix
#' @param annot Data frame with row and column annotations
#' @param idx Rows to highlight
#' @param splitx Columns to split heatmap
#' @param xtips Custom column tooltips
#' @param ytips Custom row tooltips
#' @param row_clust Cluster rows? Default is TRUE.
#' @param row_annot_width Width for row annotations. Default is 0.03.
#' @param scale Scaling for data. Default is "row.center".
#' @param colors Vector of colors to use. Default is RColorBrewer Set1.
#' @param lmar Label margin parameter. Default is 60.
#' @param rowcex Row text size scaling. Default is 1.
#' @param colcex Column text size scaling. Default is 1.
#' @param show_legend Show color legend? Default is TRUE.
#'
#' @return A plotly split heatmap object
#'
#' @details This function generates an interactive split heatmap from a data matrix.
#' Rows and/or columns can be split according to the input annotations.
#' Data is scaled and clustered using hierarchical clustering by default.
#' Various graphical parameters like colors, labels, text sizes can be adjusted.
#'
#' @export
pgx.splitHeatmapFromMatrix <- function(X, annot, idx = NULL, splitx = NULL,
                                       xtips = NULL, ytips = NULL, row_clust = TRUE,
                                       row_annot_width = 0.03, scale = "row.center",
                                       colors = NULL, lmar = 60,
                                       rowcex = 1, colcex = 1, show_legend = TRUE) {
  ## constants
  col_annot_height <- 0.021
  if (!is.null(idx)) idx <- as.character(idx)
  if (!is.null(splitx)) splitx <- as.character(splitx)

  ## --------- defaults
  if (is.null(xtips)) xtips <- colnames(X)
  if (is.null(ytips)) ytips <- rownames(X)
  if (is.null(names(xtips))) names(xtips) <- colnames(X)
  if (is.null(names(ytips))) names(ytips) <- rownames(X)

  PLOTLY_COLORS <- c(
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"
  )

  ## --------- scaling
  if ("row.center" %in% scale) X <- X - rowMeans(X, na.rm = TRUE)
  if ("row" %in% scale) X <- t(scale(t(X)))
  if ("col.center" %in% scale) X <- t(t(X) - colMeans(X, na.rm = TRUE))
  if ("col" %in% scale) X <- scale(X)
  if ("row.bmc" %in% scale) {
    if (is.null(splitx)) {
      X <- X - rowMeans(X, na.rm = TRUE)
    } else {
      for (ux in unique(splitx)) {
        jj <- which(splitx == ux)
        X[, jj] <- X[, jj, drop = FALSE] - rowMeans(X[, jj, drop = FALSE], na.rm = TRUE)
      }
    }
  }

  ## ------ split Y-axis (genes) by factor
  hc.order <- function(x) {
    suppressWarnings(dd <- stats::as.dist(1 - stats::cor(t(x), use = "pairwise")))
    if (sum(is.na(dd))) dd[is.na(dd)] <- 1
    hc <- fastcluster::hclust(dd, method = "ward.D2")
    rownames(x)[hc$order]
  }
  if (!is.null(idx)) {
    if (row_clust) {
      kk <- tapply(rownames(X), idx, function(k) c(hc.order(X[k, ]), "   "))
    } else {
      kk <- tapply(rownames(X), idx, function(k) c(k, "   "))
    }
    idx <- tapply(idx, idx, function(k) c(k, NA))
    idx <- as.vector(unlist(idx))
    kk <- unlist(kk)
    kk <- kk[1:(length(kk) - 1)] ## remove trailing spacer
    idx <- idx[1:(length(idx) - 1)]
    X <- rbind(X, "   " = 0)[kk, ]

    ## invert
    X <- X[nrow(X):1, ]
    idx <- rev(idx)
  } else {
    if (row_clust) {
      kk <- hc.order(X[, ])
      X <- X[kk, ]
    }
  }

  ## ------ split X-axis by some group factor
  if (!is.null(splitx)) {
    xx <- tapply(colnames(X), splitx, function(i) X[, i, drop = FALSE])
  } else {
    xx <- list("Samples" = X)
  }

  ## ------- set colors
  if (!is.null(annot)) {
    colors0 <- rep("Set2", ncol(annot))
    names(colors0) <- colnames(annot)
    if (!is.null(colors) && any(names(colors) %in% names(colors0))) {
      for (v in intersect(names(colors), names(colors0))) colors0[[v]] <- colors[[v]]
    }

    ## ------- annot need to be factor
    annotF <- data.frame(as.list(annot), stringsAsFactors = TRUE)
    rownames(annotF) <- rownames(annot)
  } else {
    colors0 <- NULL
  }

  grid_params <- iheatmapr::setup_colorbar_grid(
    nrows = 5,
    y_length = 0.16,
    y_spacing = 0.18,
    x_spacing = 0.2,
    x_start = 1.1,
    y_start = 0.85
  )

  ## maximize plot area
  mar <- list(l = 50, r = 50, b = 100, t = 100, pad = 4)
  mar <- list(l = lmar, r = 0, b = 5, t = 0, pad = 3)
  ex <- ncol(X) / ncol(xx[[1]]) ## expansion factor
  hc <- NULL
  if (NCOL(xx[[1]]) > 1) {
    hc <- fastcluster::hclust(stats::as.dist(1 - stats::cor(xx[[1]], use = "pairwise")))
  }
  dd <- 0.08 * ex ## left dendrogram width

  hovertemplate <- "Row: %{y} <br>Column: %{x}<br>Value: %{z}<extra> </extra>"
  tooltip <- iheatmapr::setup_tooltip_options(
    prepend_row = "Row: ", prepend_col = "Column: ",
    prepend_value = "Value: "
  )

  x1 <- xx[[1]]


  plt <- iheatmapr::main_heatmap(
    x1,
    name = "expression",
    colorbar_grid = grid_params,
    x = xtips[colnames(x1)],
    y = ytips[rownames(x1)],
    tooltip = tooltip,
    layout = list(margin = mar)
  )

  if (!is.null(hc)) {
    plt <- plt %>% iheatmapr::add_col_dendro(hc, size = 0.06)
  }

  if (length(xx) > 1) {
    plt <- plt %>% iheatmapr::add_col_title(names(xx)[1], side = "top")
  }

  if (!is.null(annot)) {
    plt <- plt %>%
      iheatmapr.add_col_annotation(
        annotation = annotF[colnames(xx[[1]]), , drop = FALSE],
        size = col_annot_height,
        buffer = 0.005, side = "bottom",
        show_legend = show_legend, show_title = TRUE,
        colors = colors0
      )
  }

  if (ncol(X) < 100 && colcex > 0) {
    plt <- plt %>% iheatmapr::add_col_labels(
      side = "bottom",
      size = 0.15 * colcex,
      font = list(size = 11 * colcex)
    )
  }

  if (length(xx) > 1) {
    sizes <- sapply(xx, ncol) / ncol(xx[[1]])
    sizes
    i <- 5
    i <- 2
    for (i in 2:length(xx)) {
      x1 <- xx[[i]]


      hc <- NULL
      if (NCOL(x1) > 1) {
        hc <- fastcluster::hclust(stats::as.dist(1 - stats::cor(x1, use = "pairwise")))
      }

      plt <- plt %>%
        iheatmapr::add_main_heatmap(
          x1, ,
          name = "expression",
          x = xtips[colnames(x1)],
          y = ytips[rownames(x1)],
          tooltip = tooltip,
          size = sizes[i],
          buffer = 0.007 * ex
        )

      if (!is.null(hc)) {
        plt <- plt %>%
          iheatmapr::add_col_dendro(hc, size = 0.06)
        ## add_col_clustering() %>%
      }

      plt <- plt %>%
        iheatmapr::add_col_title(names(xx)[i], side = "top") %>%
        iheatmapr.add_col_annotation(
          annotation = data.frame(annotF[colnames(x1), , drop = FALSE]),
          size = col_annot_height,
          buffer = 0.005, side = "bottom",
          show_title = FALSE, show_legend = show_legend,
          colors = colors0
        )

      if (ncol(X) < 100 && colcex > 0) {
        plt <- plt %>%
          iheatmapr::add_col_labels(
            side = "bottom",
            size = 0.15 * colcex,
            font = list(size = 11 * colcex)
          )
      }
    }
  }

  ## ----------- row annotation (i.e. gene groups)
  if (!is.null(idx) && !is.null(row_annot_width) && row_annot_width > 0) {
    plt <- iheatmapr::add_row_annotation(
      plt, data.frame("gene.module" = as.factor(idx)),
      size = row_annot_width * ex,
      colors = list("gene.module" = c(omics_pal_d("muted_light")(nlevels(as.factor(idx))))),
      show_colorbar = show_legend
    )
  }

  ## ----------- add gene/geneset names
  if (rowcex > 0) {
    gnames <- rownames(X)
    gnames <- gsub("[&].*[;]", "", gnames) ## HTML special garbage...
    gnames <- gsub("^.*:", "", gnames) ## strip prefix
    gnames <- shortstring(gnames, 25) ## shorten
    gnames <- sub("   ", "-", gnames)

    maxlen <- max(sapply(gnames, nchar))
    w <- ifelse(maxlen >= 20, 0.45, 0.20)
    s1 <- ifelse(maxlen >= 20, 9, 11) * rowcex
    plt <- iheatmapr::add_row_labels(
      plt,
      side = "right", ticktext = gnames,
      size = w * ex, font = list(size = s1)
    )
  }

  return(plt)
}


#' Box plot using plotly
#'
#' @param data Data frame to plot.
#' @param x Column in data to use for x-axis grouping. Default NULL.
#' @param y Column in data to use for y-axis values. Default NULL.
#' @param title Plot title text. Default NULL.
#' @param color Box fill color. Default "#3181de".
#' @param fillcolor Box inside fill color. Default "#2fb5e3".
#' @param linecolor Box border color. Default "#3181de".
#' @param hoverinfo Determines hover label info. Default "y".
#' @param hoverformat Format for hover labels. Default ".2f".
#' @param yaxistitle Show y-axis title? Logical. Default FALSE.
#' @param xaxistitle Show x-axis title? Logical. Default FALSE.
#' @param font_family Font family for text. Default "Lato".
#' @param margin Plot margins. Default c(10,10,10,10).
#'
#' @return A plotly boxplot object.
#'
#' @export
pgx.boxplot.PLOTLY <- function(
    data,
    x = NULL,
    y = NULL,
    title = NULL,
    color = "#3181de",
    fillcolor = "#2fb5e3",
    linecolor = "#3181de",
    hoverinfo = "y",
    hoverformat = ".2f",
    yaxistitle = FALSE,
    xaxistitle = FALSE,
    font_family = "Lato",
    margin = list(l = 10, r = 10, b = 10, t = 10)) {
  plotly::plot_ly(
    data = data,
    x = ~ get(x),
    y = ~ get(y),
    type = "box",
    marker = list(
      color = color,
      fillcolor = fillcolor
    ),
    line = list(color = linecolor),
    hoverinfo = hoverinfo
  ) %>%
    plotly::layout(
      title = title,
      yaxis = list(
        title = yaxistitle,
        hoverformat = hoverformat
      ),
      xaxis = list(title = xaxistitle),
      font = list(family = font_family),
      margin = margin
    )
}


#' Bar plot using plotly
#'
#' @param data Data frame to plot
#' @param x Column in data to use for x-axis. Default NULL.
#' @param y Column in data to use for y-axis. Default NULL.
#' @param title Plot title. Default NULL.
#' @param color Bar color. Default "#3181de".
#' @param fillcolor Bar fill color. Default "#2fb5e3".
#' @param linecolor Bar border color. Default "#3181de".
#' @param titlecolor Title text color. Default "#1f77b4".
#' @param hoverinfo Determines hover label info. Default "y".
#' @param hoverformat Format for hover labels. Default ".2f".
#' @param yaxistitle Show y-axis title? Default FALSE.
#' @param xaxistitle Show x-axis title? Default FALSE.
#' @param xlen Length of x-axis. Default NULL for automatic.
#' @param yrange Limits for y-axis. Default NULL for automatic.
#' @param font_family Font family for text. Default "Lato".
#' @param margin Plot margins. Default c(0,0,0,0).
#' @param grouped Calculate mean and SD across groups? Default TRUE.
#' @param annotations Additional plot annotations. Default NULL.
#'
#' @return A plotly bar plot object
#'
#' @export
pgx.barplot.PLOTLY <- function(
    data,
    x = NULL,
    y = NULL,
    title = NULL,
    color = "#3181de",
    fillcolor = "#2fb5e3",
    linecolor = "#3181de",
    titlecolor = "#1f77b4",
    hoverinfo = "y",
    hoverformat = ".2f",
    yaxistitle = FALSE,
    xaxistitle = FALSE,
    xlen = NULL,
    yrange = NULL,
    font_family = "Lato",
    margin = list(l = 0, r = 0, b = 0, t = 0),
    grouped = TRUE, # true will calculate mean +/- (sd) across groups
    annotations = NULL) {
  if (is.null(x)) x <- 1
  if (is.null(y)) y <- 2

  # calculate error bars
  # calculate summary statistics for groups
  if (grouped) {
    data <- do.call(
      data.frame,
      stats::aggregate(
        data[[y]],
        list(data[[x]]),
        function(val) {
          c(mean = mean(val), sd = stats::sd(val))
        }
      )
    )
    colnames(data) <- c(x, y, "sd")
  }

  ngroups <- length(unique(data[[x]]))
  bargap <- ifelse(ngroups == 2, 0.5, NA)

  error_y <- NULL
  if (grouped) {
    error_y <- list(
      array = data[["sd"]],
      thickness = 1,
      color = "#000000"
    )
  }

  data[["short.x"]] <- data[[x]]
  if (!is.null(xlen)) {
    sx <- shortstring(data[[x]], xlen)
    i <- 1
    ## make unique: sometimes shortened names gets duplicated
    while (sum(duplicated(sx)) && i < 1000) {
      sx[which(duplicated(sx))] <- paste0(sx[which(duplicated(sx))], " ")
      i <- i + 1
    }
    data[["short.x"]] <- factor(sx, levels = sx)
  }

  p <- plotly::plot_ly(
    data = data,
    x = data[["short.x"]],
    hovertext = data[[x]] ## original long text
  ) %>%
    plotly::add_bars(
      y = data[[y]],
      error_y = error_y,
      marker = list(
        color = fillcolor
      ),
      line = list(
        color = linecolor
      ),
      textposition = "none",
      hoverinfo = hoverinfo,
      hovertemplate = paste0(
        "<b>%{hovertext}</b><br>",
        "%{yaxis.title.text}: %{y:", hoverformat, "}<br>",
        "<extra></extra>"
      )
    ) %>%
    plotly::layout(
      title = list(
        text = title,
        font = list(color = titlecolor)
      ),
      yaxis = list(
        title = yaxistitle,
        hoverformat = hoverformat,
        range = yrange
      ),
      xaxis = list(
        title = xaxistitle
      ),
      font = list(
        family = font_family
      ),
      margin = margin,
      bargap = bargap,
      annotations = annotations
    )

  return(p)
}
