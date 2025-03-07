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
    col = circlize::colorRamp2(colors = c(omics_colors("brand_blue"), omics_colors("grey"), omics_colors("red")), breaks = c(-10, 0, 10)),
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


#' @export
pgx.plotEnrichmentDotPlot <- function(pgx, contrast,
                                      ntop = 10, filter = NULL, dir = "both",
                                      main = "Enrichment Analysis") {
  gs <- pgx$gset.meta$meta[[contrast]]
  df <- data.frame(
    pathway = rownames(gs), fx = gs$meta.fx,
    pval = gs$meta.q, size = abs(gs$meta.fx)
  )
  df <- df[order(-df$fx), ]

  if (!is.null(filter)) {
    sel <- grep(filter, df$pathway, ignore.case = TRUE)
    df$pathway <- gsub(filter, "", df$pathway)
  }

  if (dir == "up") {
    sel <- head(sel, ntop)
  } else if (dir == "down") {
    sel <- tail(sel, ntop)
  } else {
    ## both
    sel <- c(head(sel, ntop / 2), tail(sel, ntop / 2))
  }
  df <- df[unique(sel), ]

  ## cleaning...
  df$pathway <- sub(".Homo.sapiens", "", df$pathway)
  df$pathway <- substring(df$pathway, 1, 60)
  df$pathway <- factor(df$pathway, levels = rev(df$pathway))


  ggplot2::ggplot(df, ggplot2::aes(x = fx, y = pathway)) +
    ggplot2::geom_point(ggplot2::aes(color = pval, size = size)) +
    ggplot2::scale_size_area(max_size = 12) +
    ggplot2::scale_color_gradient(low = "red", high = "blue") +
    ggplot2::labs(x = "Enrichment score", y = NULL, color = "P-value", size = "Score") +
    ggplot2::ggtitle(main) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = 30),
      axis.text = ggplot2::element_text(size = 30),
      title = ggplot2::element_text(size = 36),
      legend.title = ggplot2::element_text(size = 20),
      legend.text = ggplot2::element_text(size = 14)
    )
}

#' @export
pgx.dimPlot <- function(X, y, method = c("tsne", "pca", "umap"), nb = NULL, ...) {
  ## method=c('tsne','pca','umap','pacmap')
  jj <- head(order(-matrixStats::rowSds(X, na.rm = TRUE)), 1000)
  X1 <- X[jj, ]
  X1 <- X1 - rowMeans(X1, na.rm = TRUE)
  if (any(is.na(X1))) {
    X1 <- svdImpute2(X1)
  }
  if (ncol(X1) < 20) {
    X1 <- cbind(X1, X1, X1)
  }

  if (is.null(nb)) nb <- ceiling(min(15, dim(X) / 8))
  for (m in method) {
    if (m == "umap") pos <- try(uwot::umap(t(X1), n_neighbors = max(2, nb)))
    if (m == "tsne") {
      pos <- try(Rtsne::Rtsne(t(X1),
        perplexity = 2 * nb,
        check_duplicates = FALSE
      )$Y)
    }
    if (m == "pca") pos <- try(irlba::irlba(X1, nv = 2, nu = 0)$v)
    if (m == "pacmap") pos <- try(pacmap(t(X1)))
    if ("try-errror" %in% class(pos)) {
      pos <- matrix(0, nrow = ncol(X), ncol = 2)
      rownames(pos) <- colnames(X)
      pgx.scatterPlotXY(pos, var = y, title = m, ...)
    } else {
      pos <- pos[1:ncol(X), ]
      rownames(pos) <- colnames(X)
      pgx.scatterPlotXY(pos, var = y, title = m, ...)
    }
  }
}

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
pgx.scatterPlot <- function(pgx, samples = NULL, pheno = NULL,
                            gene = NULL, geneset = NULL,
                            contrast = NULL, level = "gene",
                            method = "tsne", pos = NULL,
                            plotlib = "base", ...) {
  ## Scatter t-SNE plot samples (or genes) colored on phenotype,
  ## gene expression, geneset expresssion or (correlation with)
  ## contrast.
  ##

  if (is.null(pos)) {
    pos <- pgx$tsne2d
    cpos <- pgx$cluster$pos
    if (method == "tsne" && "tsne2d" %in% names(cpos)) pos <- cpos[["tsne2d"]]
    if (method == "pca" && "tsne2d" %in% names(cpos)) pos <- cpos[["pca2d"]]
    if (method == "umap" && "tsne2d" %in% names(cpos)) pos <- cpos[["umap2d"]]
  }

  if (nrow(pos) != nrow(pgx$samples)) {
    stop("[pgx.scatterPlot] dimension mismatch of positions")
  }

  plt <- NULL
  vartype <- c("numeric", "factor")[1]
  if (!is.null(pheno)) {
    var <- pgx$samples[rownames(pos), pheno]
    title <- pheno
    vartype <- "factor"
  }
  if (!is.null(contrast)) {
    ct.matrix <- pgx$contrasts
    if (is.null(ct.matrix)) ct.matrix <- pgx$model.parameters$exp.matrix
    ct.matrix <- playbase::contrastAsLabels(ct.matrix)
    var <- ct.matrix[rownames(pos), contrast]
    title <- contrast
    vartype <- "factor"
  }
  if (!is.null(gene)) {
    var <- pgx$X[gene, rownames(pos)]
    title <- gene
    vartype <- "numeric"
  }
  if (!is.null(geneset)) {
    var <- pgx$gsetX[geneset, rownames(pos)]
    title <- gene
    vartype <- "numeric"
  }
  if (!is.null(samples)) {
    sel <- 1:nrow(pos)
    if (is.character(samples)) {
      sel <- which(rownames(pos) %in% samples)
    }
    if (is.integer(samples)) {
      sel <- samples
    }
    pos <- pos[sel, ]
    var <- var[sel]
  }
  plt <- pgx.scatterPlotXY(
    pos,
    var,
    plotlib = plotlib, #
    xlab = colnames(pos)[1],
    ylab = colnames(pos)[2],
    type = vartype,
    ...
  )

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

  x0 <- range(as.vector(apply(F, 2, stats::quantile, probs = c(0.001, 0.999))), na.rm = TRUE)
  x1 <- range(as.vector(apply(F2, 2, stats::quantile, probs = c(0.001, 0.999))), na.rm = TRUE)
  x0 <- range(as.vector(F), na.rm = TRUE)
  x1 <- range(as.vector(F2), na.rm = TRUE)
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
    ## X[[i]] <- X[[i]] / apply(X[[i]], 1, stats::sd, na.rm = TRUE)
    X[[i]] <- X[[i]] / matrixStats::rowSds(X[[i]], na.rm = TRUE)
  }

  ## Counts cross-table between matrices
  M <- list()
  i <- 1
  for (i in 1:(length(X) - 1)) {
    mm <- pmax(X[[i]], 0) %*% t(pmax(X[[i + 1]], 0))
    mm <- mm**4
    mm <- mm / mean(mm, na.rm = TRUE)

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
      ww <- ww / max(ww, na.rm = TRUE) ## normalize??
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
  nodes$x <- (vlevel - 1) / max(vlevel - 1, na.rm = TRUE)

  if (fill) {
    wt <- igraph::E(gr)$weight
    ev2 <- 0.05 + 0.55 * (wt / max(wt, na.rm = TRUE))
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
                        psig = 0.05, fc = 1, cex = 1, sig.par=c("q","p")[1],
                        p.min = NULL, fc.max = NULL, hilight = NULL,
                        label = NULL, cex.lab = 1, ntop = 20,
                        colors = c(
                          up = "#f23451", notsig = "#707070AA",
                          notsel = "#cccccc88", down = "#3181de"
                        ),
                        filt = NULL, datatype = NULL,
                        title = NULL, xlim = NULL, ylim = NULL,
                        xlab = "effect size (logFC)",
                        ylab = NULL,
                        repel = FALSE, labeltype = NULL,
                        set.par = TRUE, plotlib = "base", data = FALSE) {
  if (is.integer(contrast)) contrast <- names(pgx$gx.meta$meta)[contrast]
  res <- NULL
  if (level == "gene") {
    res <- pgx.getMetaMatrix(pgx, methods, level = "gene")
  } else if (level == "geneset") {
    res <- pgx.getMetaMatrix(pgx, methods, level = "geneset")
  } else {
    stop("FATAL:: invalid level=", level)
  }
  
  ## Label or hilight gene sets
  gsets <- colnames(pgx$GMT)
  if(!is.null(hilight) && any(hilight %in% gsets)) {
    G1 <- pgx$GMT[,which(gsets %in% hilight),drop=FALSE]
    hilight <- names(which(Matrix::rowSums(G1!=0)>0))
  }
  if(!is.null(label) && any(label %in% gsets)) {
    G1 <- pgx$GMT[,which(gsets %in% label),drop=FALSE]
    label <- names(which(Matrix::rowSums(G1!=0)>0))
  }

  f <- res$fc[, contrast]
  if(sig.par=="p" && "pv" %in% names(res)) {
    q <- res$pv[, contrast]  ## should be p-value!!
    if(is.null(ylab)) ylab <- "significance (-log10p)"
  } else {
    q <- res$qv[, contrast]
    if(is.null(ylab)) ylab <- "significance (-log10q)"    
  }
  
  if (!is.null(filt)) {
    sel <- grep(filt, names(f))
    f <- f[sel]
    q <- q[sel]
  }
  
  sig <- (q <= psig & abs(f) >= fc)
  df <- data.frame(fc = f, y = -log10(q), q = q, sig = sig)

  if(!is.null(datatype)) {
    dt <- sub(":.*","",rownames(df))
    df <- df[ dt %in% datatype,]
  }
      
  hilight0 <- hilight
  if (FALSE && is.null(hilight)) {
    hilight <- rownames(res$fc)
  }

  if (is.null(label)) {
    label <- rownames(df)
    if (!is.null(hilight)) label <- intersect(label, hilight)
    ## if(!is.null(sig)) label <- intersect(label, names(sig[sig == TRUE]))
  }

  if(!is.null(label))  label <- map_probes(pgx$genes, label)
  if(!is.null(hilight))  hilight <- map_probes(pgx$genes, hilight)  

  wt <- rowSums(scale(df[,c("fc","y")], center = FALSE)**2, na.rm = TRUE)
  label <- label[order(-wt[label])]
  label <- Matrix::head(label, ntop) ## label

  if(is.null(fc.max)) fc.max <- max(abs(df$fc),na.rm=TRUE)
  if (is.null(xlim) && !is.null(fc.max)) {
    xlim <- c(-1.1, 1.1) * fc.max
  }

  if(is.null(p.min)) p.min <- max(min(q,na.rm=TRUE),10^-99)
  if (is.null(ylim) && !is.null(p.min)) {
    ylim <- c(0, -log10(p.min)) * 1.05
  }

  if (data) {
    return(df)
  }

  if (is.null(title)) title <- contrast

  ## if(is.null(colors)) {
  ##   colors = c(
  ##     up = "#f23451",
  ##     notsig = "#707070AA",
  ##     notsel = "#cccccc44",
  ##     down = "#3181de"
  ##   )
  ## }

  label.names = rownames(df)
  if( !is.null(labeltype) && labeltype %in% colnames(pgx$genes)) {
    ii <- match(rownames(df), rownames(pgx$genes))
    label.names <- pgx$genes[ii, labeltype]
  }
  
  p <- NULL
  if (plotlib == "ggplot") {
    p <- ggVolcano(
      x = df$fc,
      y = df$y,
      names = rownames(df),
      title = title,
      label.names = label.names,
      facet = NULL,
      highlight = hilight,
      xlab = xlab,
      ylab = ylab,
      lfc = fc,
      psig = psig,
      showlegend = TRUE,
      marker.size = 2 * cex,
      marker.alpha = 0.7,
      label.cex = 5 * cex.lab,
      axis.text.size = 14,
      label = label,
      colors = colors
    )
  } else if (plotlib == "plotly") {
    p <- plotlyVolcano(
      x = df$fc,
      y = df$y,
      names = rownames(df),
      label.names = label.names,
      group.names = c("group1", "group2"),
      title = title,
      xlab = xlab,
      ylab = ylab,
      lfc = fc,
      psig = psig,
      showlegend = TRUE,
      highlight = hilight,
      marker.size = 8 * cex,
      label = label,
      label.cex = 1.2 * cex.lab,
      max.absy = NULL,
      color_up_down = TRUE,
      colors = colors,
      marker.type = "scatter",
      displayModeBar = TRUE,
      source = "plot1"
    )
  } else if (plotlib == "base") {
    if (is.null(hilight0)) hilight <- NULL
    plot_volcano(
      x = df$fc,
      y = df$q,
      names = rownames(df),
      label.names = label.names,
      highlight = hilight,
      label = label,
      ntop = 9999,  ## done earlier
      title = title,
      xlab = xlab,
      ylab = ylab,
      xlim = xlim,
      ylim = ylim,
      lfc = fc,
      psig = psig,
      sig = sig,
      showlegend = TRUE,
      cex = 1.2 * cex,
      cex.lab = 1.1 * cex.lab,
      colors = colors,
      set.par = set.par,      
      repel = repel
    )
  } else {
    p <- pgx.scatterPlotXY(
      df[, 1:2],
      var = df$sig,
      type = "factor",
      title = title,
      xlab = xlab,
      ylab = ylab,
      hilight = label, #
      cex = 1.3 * cex,
      cex.lab = 1.4 * cex.lab,
      cex.title = 1.0,
      xlim = xlim,
      ylim = ylim,
      legend = FALSE,
      col = colors,
      opacity = 1,
      set.par = set.par,
      plotlib = plotlib
    )
    if (plotlib == "base") {
      abline(v = 0, lty = 1, lwd = 0.5)
      abline(h = 0, lty = 1, lwd = 0.5)
      abline(v = c(-1, 1) * fc, lty = 2, lwd = 0.5)
      abline(h = -log10(psig), lty = 2, lwd = 0.5)
    }
  }

  p
}


plot_volcano <- function(x,
                         y,
                         names,
                         label.names = names,
                         facet = NULL,
                         highlight = NULL,
                         label = NULL,
                         ntop = 10,
                         xlab = "effect size (logFC)",
                         ylab = "significance (-log10p)",
                         lfc = 1,
                         psig = 0.05,
                         plim = 12,
                         sig = NULL,
                         showlegend = TRUE,
                         cex = 1,
                         cex.lab = 1,
                         xlim = NULL,
                         ylim = NULL,
                         marker.alpha = 0.7,
                         axis.text.size = 14,
                         title = "Volcano plot",
                         set.par = TRUE,                         
                         colors = c(
                           up = "#f23451",
                           notsig = "#707070AA",
                           notsel = "#cccccc88",
                           down = "#3181de"
                         ),
                         repel = TRUE,
                         ...) {
  if (is.null(sig)) {
    sig <- (y <= psig & abs(x) >= lfc)
  }
  if (FALSE && is.null(highlight)) {
    highlight <- rownames(df)
  }

  log10.y <- -log10(y + 10^(-plim))
  pos <- cbind(x, log10.y)
  if (!is.null(names)) rownames(pos) <- names

  if (is.null(label)) {
    wt <- rowSums(scale(pos, center = FALSE)**2, na.rm = TRUE)
    label <- rownames(pos)[order(-wt)]
    if (!is.null(highlight)) label <- intersect(label, highlight)
    ## if(!is.null(sig)) label <- intersect(label, names(sig[sig == TRUE]))
  }
  label <- head(label, ntop) ## label

  p <- pgx.scatterPlotXY.BASE(
    pos = pos,
    var = sig,
    type = "factor",
    title = title,
    labels = label.names,
    xlab = xlab,
    ylab = ylab,
    hilight = highlight, #
    hilight2 = label, #
    cex = 1.3 * cex,
    cex.lab = 1.4 * cex.lab,
    cex.title = 1.0,
    xlim = xlim,
    ylim = ylim,
    legend = FALSE,
    color_up_down = TRUE,
    col = colors[c("notsig", "up", "down")],
    na.color = colors["notsel"],
    opacity = 1,
    repel = repel,
    set.par = set.par,
    ...
  )

  abline(v = 0, lty = 1, lwd = 0.5)
  abline(h = 0, lty = 1, lwd = 0.5)
  abline(v = c(-1, 1) * lfc[1], lty = 2, lwd = 0.5)
  abline(h = -log10(psig[1]), lty = 2, lwd = 0.5)

  if (showlegend) {
    cc <- colors[c("notsig", "up", "down")]
    legend("bottomleft",
      legend = c("Not significant", "Up", "Down"),
      pch = 20, col = cc,
      pt.cex = 1.5 * cex, cex = 0.9, 
      bg = "#FFFFFF99", bty = "o", border = NA
    )
  }
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
ggVolcano <- function(x,
                      y,
                      names,
                      label.names = names,
                      facet = NULL,
                      highlight = NULL,
                      xlab = "effect size (logFC)",
                      ylab = "significance (-log10q)",
                      lfc = 1,
                      psig = 0.05,
                      ylim = NULL,
                      showlegend = TRUE,
                      marker.size = 2.5,
                      marker.alpha = 0.7,
                      label.cex = 4.5,
                      axis.text.size = 14,
                      label = NULL,
                      title = "Volcano plot",
                      colors = c(
                        up = "#f23451",
                        notsig = "#8F8F8F",
                        notsel = "#cccccc88",
                        down = "#3181de"
                      ),
                      girafe = FALSE) {
  if (is.null(highlight)) highlight <- names
  if (showlegend) {
    legend <- "right"
  } else {
    legend <- "none"
  }
  df <- data.frame(
    fc = x,
    y = y,
    name = names,
    tooltip = label.names
  )
  if (!is.null(facet)) {
    df$facet <- facet
  }
  df$category <- ifelse(
    df$y > -log10(psig) & df$fc > lfc, "Up",
    ifelse(df$y > -log10(psig) & df$fc < -lfc, "Down", "Not significant")
  )
  df$label <- ifelse(names %in% label | label.names %in% label, label.names, NA)
  df$category[!(names %in% highlight | label.names %in% highlight)] <- "Not selected"

  df$tooltip <- gsub("[\\'\\`-]", "", df$tooltip)
  df$name <- gsub("[\\'\\`-]", "", df$name)

  if (is.null(ylim)) ylim <- max(y, na.rm = TRUE) * 1.1

  plt <- ggplot2::ggplot(df, ggplot2::aes(x = fc, y = y)) +
    ggplot2::geom_point(
      ggplot2::aes(color = category),
      alpha = marker.alpha, size = marker.size
    ) +
    ggplot2::scale_color_manual(values = c(
      "Down" = colors[["down"]],
      "Up" = colors[["up"]],
      "Not significant" = colors[["notsig"]],
      "Not selected" = colors[["notsel"]]
    ))

  if (!girafe) {
    plt <- plt +
      ggplot2::geom_point(
        data = df[df$category == "Not selected", ],
        ggplot2::aes(x = fc, y = y),
        color = colors[["notsel"]],
        alpha = marker.alpha,
        size = marker.size
      ) +
      ggplot2::geom_point(
        data = df[df$category == "Down", ],
        ggplot2::aes(x = fc, y = y),
        color = colors[["down"]],
        alpha = 1,
        size = marker.size
      ) +
      ggplot2::geom_point(
        data = df[df$category == "Up", ],
        ggplot2::aes(x = fc, y = y),
        color = colors[["up"]],
        alpha = 1,
        size = marker.size
      ) +
      ggplot2::geom_point(
        data = df[df$category == "Not significant", ],
        ggplot2::aes(x = fc, y = y),
        color = colors[["notsig"]],
        alpha = 1,
        size = marker.size
      )
  }
  if (girafe) {
    plt <- plt +
      ggiraph::geom_point_interactive(
        data = df[df$category == "Not selected", ],
        ggplot2::aes(x = fc, y = y, tooltip = tooltip, data_id = name),
        color = colors[["notsel"]],
        alpha = marker.alpha,
        size = marker.size
      ) +
      ggiraph::geom_point_interactive(
        data = df[df$category == "Down", ],
        ggplot2::aes(x = fc, y = y, tooltip = tooltip, data_id = name),
        color = colors[["down"]],
        alpha = 1,
        size = marker.size
      ) +
      ggiraph::geom_point_interactive(
        data = df[df$category == "Up", ],
        ggplot2::aes(x = fc, y = y, tooltip = tooltip, data_id = name),
        color = colors[["up"]],
        alpha = 1,
        size = marker.size
      ) +
      ggiraph::geom_point_interactive(
        data = df[df$category == "Not significant", ],
        ggplot2::aes(x = fc, y = y, tooltip = tooltip, data_id = name),
        color = colors[["notsig"]],
        alpha = 1,
        size = marker.size
      )
  }

  plt <- plt +
    ggrepel::geom_label_repel(
      ggplot2::aes(label = label, color = category),
      size = label.cex, family = "lato", box.padding = 0.1, max.overlaps = 20, show.legend = FALSE
    )

  plt <- plt +
    ggplot2::geom_hline(yintercept = -log10(psig), linetype = "dashed", color = "gray") +
    ggplot2::geom_vline(xintercept = c(-lfc, lfc), linetype = "dashed", color = "gray") +
    ggplot2::geom_vline(xintercept = 0, linetype = "solid", color = "darkgrey") +
    ggplot2::scale_y_continuous(
      limits = c(0, ylim),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.1, 0))) +
    ggplot2::labs(x = xlab, y = ylab) +
    guides(colour = guide_legend(reverse = T)) +
    ggplot2::theme_minimal(base_size = 15) +
    ggplot2::ggtitle(title) +
    ggplot2::theme(
      legend.position = legend,
      legend.justification = "top",
      legend.title = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(color = "darkgrey"),
      axis.title = ggplot2::element_text(face = "plain", size = axis.text.size, family = "lato"),
      axis.text = ggplot2::element_text(family = "lato"),
      plot.margin = ggplot2::margin(l = 9, b = 3, t = 9, r = 9)
    )
  if (!is.null(facet)) {
    ncol_row <- ceiling(sqrt(length(unique(facet))))
    plt <- plt +
      ggplot2::facet_wrap(facet ~ ., nrow = ncol_row, ncol = ncol_row)
  }

  if (girafe) {
    plt <- ggiraph::girafe(ggobj = plt, pointsize = 12, width_svg = NULL, height_svg = NULL)
  }
  return(plt)
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
                       cex = 1, cex.lab = 0.8, hilight = NULL, label = NULL,
                       ntop = 20, plotlib = "base", data = FALSE) {
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
  df <- cbind(x = m, y = f)
  # <- gg
  cpal <- c("grey60", "red3")
  if (is.null(hilight)) {
    wt <- rowSums(scale(cbind(df, m), center = FALSE)**2, na.rm = TRUE)
    hilight <- rownames(df)[order(-wt)]
    hilight <- intersect(hilight, names(sig[sig == TRUE]))
  }
  hilight <- Matrix::head(hilight, ntop)

  if (data) {
    return(df)
  }

  p <- pgx.scatterPlotXY(
    df,
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
                                cex = 1, cex.lab = 0.8, label = NULL,
                                psig = 0.05, fc = 1, level = "gene",
                                ntop = 20, dir = 0, plotlib = "base",
                                data = FALSE) {
  if (is.numeric(contrast)) contrast <- names(pgx$gx.meta$meta)[contrast]
  exp.matrix <- pgx$model.parameters$exp.matrix
  ct <- exp.matrix[, contrast]
  ii <- which(ct < 0)
  jj <- which(ct > 0)
  if (level == "gene") {
    x0 <- rowMeans(pgx$X[, ii, drop = FALSE], na.rm = TRUE)
    x1 <- rowMeans(pgx$X[, jj, drop = FALSE], na.rm = TRUE)
    xy <- cbind(x0, x1)
    gg <- rownames(pgx$gx.meta$meta[[contrast]])
    fx <- pgx$gx.meta$meta[[contrast]]$meta.fx
    q <- pgx$gx.meta$meta[[contrast]]$meta.q
    names(fx) <- names(q) <- gg
  } else if (level == "gene") {
    x0 <- rowMeans(pgx$gsetX[, ii, drop = FALSE], na.rm = TRUE)
    x1 <- rowMeans(pgx$gsetX[, jj, drop = FALSE], na.rm = TRUE)
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

  if (is.null(hilight)) hilight <- names(fx)
  hilight <- intersect(hilight, names(fx))

  if (is.null(label)) {
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
    label <- top.gg
  }
  label <- Matrix::head(label, ntop)

  sig <- 1 * (q < psig & abs(fx) > fc)
  names(sig) <- gg

  tt <- contrast

  if (data) {
    return(xy)
  }

  pgx.scatterPlotXY(
    xy,
    var = sig, type = "factor", title = tt,
    xlab = xlab, ylab = ylab,
    hilight = hilight,
    hilight2 = label,
    cex = 0.9 * cex,
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
                             pos = NULL, ntop = 20, cex = 1,
                             cex.lab = 0.8, cex.legend = 1,
                             hilight = NULL, label = NULL,
                             title = NULL, zfix = FALSE,
                             set.par = TRUE, par.sq = FALSE,
                             level = "gene", plotlib = "ggplot",
                             data = FALSE, labeltype = "feature") {
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

  if (is.null(pos) && "cluster.genes" %in% names(pgx)) {
    pos <- pgx$cluster.genes$pos[["umap2d"]]
  } else if (is.null(pos) && !"cluster.genes" %in% names(pgx)) {
    message("[WARNING] no cluster.genes in object. please supply positions")
    return(NULL)
  }

  pos <- pos[which(rowMeans(is.na(pos)) == 0), , drop = FALSE]
  cm <- intersect(rownames(pos), rownames(F))
  F <- F[cm, , drop = FALSE]
  pos <- pos[cm, ]

  if (set.par) {
    nc <- ceiling(sqrt(ncol(F)))
    nr <- ceiling(ncol(F) / nc)
    if (par.sq) nr <- nc
    graphics::par(mfrow = c(nr, nc))
  }

  ## z-scale
  zlim <- NULL
  if (zfix) {
    zlim <- stats::quantile(F, probs = c(0.002, 0.998), na.rm = TRUE)
    zlim
  }

  # Rename labels
  if (!is.null(rownames(F))) {
    rownames(F) <- make.names(playbase::probe2symbol(rownames(F), pgx$genes, labeltype, fill_na = TRUE), unique = TRUE)
  }
  if (!is.null(rownames(pos))) {
    rownames(pos) <- make.names(playbase::probe2symbol(rownames(pos), pgx$genes, labeltype, fill_na = TRUE), unique = TRUE)
  }
  if (!is.null(label)) {
    label <- make.names(playbase::probe2symbol(label, pgx$genes, labeltype, fill_na = TRUE), unique = TRUE)
  }
  if (!is.null(hilight)) {
    hilight <- make.names(playbase::probe2symbol(hilight, pgx$genes, labeltype, fill_na = TRUE), unique = TRUE)
  }

  plist <- list()
  i <- 1
  for (i in 1:ncol(F)) {
    title1 <- contrast[i]
    if (is.null(title)) title1 <- paste(title, contrast[i])

    f1 <- F[, i]

    this.hilight <- hilight
    if (is.null(this.hilight)) {
      this.hilight <- names(sort(-abs(f1)))
    }
    this.hilight <- intersect(this.hilight, rownames(F))
    this.label <- label
    if (is.null(this.label)) {
      this.label <- Matrix::head(this.hilight, ntop) ## label
    }
    this.label <- intersect(this.label, rownames(F))
    opacity <- ifelse(length(this.hilight) > 0, 0.66, 1)

    if (data) {
      return(
        cbind(pos, f1)
      )
    }

    p1 <- pgx.scatterPlotXY(
      pos,
      var = f1,
      type = "numeric",
      xlab = "UMAP-x  (genes)",
      ylab = "UMAP-y  (genes)",
      hilight = this.hilight,
      hilight2 = this.label,
      hilight.lwd = 0.0,
      hilight2.lwd = 0.8,
      zlim = zlim,
      zsym = TRUE,
      softmax = 1,
      cex = cex,
      cex.lab = cex.lab,
      title = title1,
      cex.title = 1.0,
      cex.legend = cex.legend,
      legend = TRUE,
      opacity = 0.5,
      set.par = set.par,
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
                               plotly.margin = NULL, verbose = 0,
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
    if (verbose) message("[pgx.plotExpression] using group names from argument")
  }

  ## if a named contrast table is available it is safer
  if (is.null(group.names) && "contrasts" %in% names(pgx)) {
    if (verbose) message("[pgx.plotExpression] parsing group names from contrast labels")
    contr.labels <- pgx$contrasts[, comp]
    contr.idx <- expmat[, comp]
    group1 <- names(which.max(table(contr.labels[contr.idx > 0])))
    group0 <- names(which.max(table(contr.labels[contr.idx < 0])))
    group.names <- c(group0, group1)
  }
  group.names

  ## Otherwise we guess from the contrast title but this is dangerous
  if (is.null(group.names) && grepl("_vs_|_VS_", comp)) {
    if (verbose) message("[pgx.plotExpression] parsing group names contrast name (warning!)")
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
    if (max(nchar(group.names), na.rm = TRUE) >= 7) srt <- 45
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
    if ("other" %in% xgroup && !"other" %in% levels0) levels0 <- c(levels0, "other")
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
    gx <- pgx$X[which(rownames(pgx$X) == probe), rownames(pgx$samples)]
  }

  if (!logscale) {
    gx <- 2**(gx)
  }

  ## -------------- remove others
  if (showothers == FALSE && any(grepl("other", xgroup))) {
    jj <- !xgroup %in% "other"
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
      ylab <- "Expression (log2)"
      if (!logscale) ylab <- "Expression"
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
      if (min(gx, na.rm = TRUE) < 0) gx.min <- min(gx, na.rm = TRUE)
      ylim <- c(gx.min, 1.3 * max(gx, na.rm = TRUE))
      bx <- graphics::barplot(gx[],
        col = klr[], ylim = ylim,
        ## offset = 0, ylim=c(gx.min,max(gx)),
        las = 3, ylab = ylab, names.arg = NA, border = NA
      )
    }
  } else {
    ## GROUPED PLOTS
    if (is.null(ylab)) {
      ylab <- "Expression (log2)"
      if (!logscale) ylab <- "Expression"
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
      buffer = 0.005,
      side = "bottom",
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
    annot.df <- annot.df[hc$order, , drop = FALSE]
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
    ## X1 <- Matrix::head(X0[order(-apply(X0, 1, stats::sd, na.rm = TRUE)), ], ntop)
    X1 <- Matrix::head(X0[order(-matrixStats::rowSds(X0, na.rm = TRUE)), ], ntop)
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
                          cex.axis = 1, cbar.width = 32) {
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
  rnk.trace <- rnk.trace / max(abs(rnk.trace), na.rm = TRUE) * 0.8

  qq <- range(fc, na.rm = TRUE)
  y1 <- qq[2]
  y0 <- 0.8 * qq[1]
  dy <- ticklen * (y1 - y0)
  if (max(rnk.trace, na.rm = TRUE) >= abs(min(rnk.trace, na.rm = TRUE))) rnk.trace <- rnk.trace * abs(y1)
  if (max(rnk.trace, na.rm = TRUE) < abs(min(rnk.trace, na.rm = TRUE))) rnk.trace <- rnk.trace * abs(y0)

  cc <- sign(fc) * rank(abs(fc))
  df <- data.frame(x = rank(-fc), y = fc, trace = rnk.trace, cc = cc)

  ## downsample
  ii <- which(rownames(df) %in% gset)
  unknown_parameter <- round(nrow(df) / 1000)
  if (unknown_parameter == 0) unknown_parameter <- round(nrow(df) / 10)
  ii <- unique(c(ii, seq(1, nrow(df), unknown_parameter)))
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
    cc <- gplots::colorpanel(21, omics_colors("brand_blue"), omics_colors("grey"), omics_colors("red"))
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

  unknown_parameter <- round(nrow(df) / 200)

  if (unknown_parameter == 0) unknown_parameter <- round(nrow(df) / 5)

  ii <- seq(1, nrow(df), unknown_parameter)

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
      xaxis = list(
        title = xlab,
        gridwidth = 0.3,
        tickfont = list(size = 9 * cex.axis),
        titlefont = list(size = 11 * cex.axis)
      ),
      yaxis = list(
        title = ylab,
        gridwidth = 0.3,
        range = c(y0 - 1.1 * dy, y1),
        tickfont = list(size = 9 * cex.axis),
        titlefont = list(size = 11 * cex.axis)
      )
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
  rnk.trace <- rnk.trace / max(abs(rnk.trace), na.rm = TRUE) * 0.8

  qq <- range(fc, na.rm = TRUE)
  y1 <- qq[2]
  y0 <- qq[1]
  dy <- 0.2 * (y1 - y0)
  if (max(rnk.trace, na.rm = TRUE) >= abs(min(rnk.trace, na.rm = TRUE))) rnk.trace <- rnk.trace * abs(y1)
  if (max(rnk.trace, na.rm = TRUE) < abs(min(rnk.trace, na.rm = TRUE))) rnk.trace <- rnk.trace * abs(y0)

  cc <- sign(fc) * rank(abs(fc))
  df <- data.frame(rank = rank(-fc), fc = fc, run = rnk.trace, cc)
  jj <- which(names(fc) %in% gset)


  cpal <- colorspace::diverge_hcl(64, c = 60, l = c(30, 100), power = 1)
  ii <- 1 + 32 + range(round(32 * (fc / max(abs(fc), na.rm = TRUE))), na.rm = TRUE)
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

    zr <- range(col, na.rm = TRUE)
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
                           las = NULL, beside = FALSE,
                           legend.pos = c(0.016, 1), legend.cex = 1,
                           axis.cex = 1, label.cex = 1, horiz = FALSE,
                           cpal = "Blues",
                           bar_width = 0.7, base_size = 12, group.name = "group") {
  if (NCOL(mat) == 1) mat <- rbind(mat)
  mat <- mat[nrow(mat):1, , drop = FALSE]
  df <- reshape2::melt(t(mat), value.name = "value")
  colnames(df)[1:2] <- c("x", "y")

  df$y <- factor(df$y, levels = rownames(mat))
  df$x <- factor(df$x, levels = colnames(mat))
  if (!is.null(las) && las == 3) srt <- 90

  ##  colors <- rev(grDevices::grey.colors(nrow(mat)))
  colors <- rev(grDevices::hcl.colors(nrow(mat), cpal))
  if (nrow(mat) == 1) colors <- "grey70"

  posmode <- ifelse(beside, "dodge", "stack")
  x <- y <- value <- NULL
  if (horiz) {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = value, y = x, fill = y)) +
      ggplot2::geom_bar(
        stat = "identity", color = "black", size = 0.3,
        width = bar_width, position = posmode
      )
    p <- p + ggplot2::xlab(xlab) +
      ggplot2::ylab(ylab) +
      ggplot2::labs(fill = group.name) +
      ggplot2::ggtitle(main) +
      ggplot2::scale_fill_manual(values = colors) +
      ggplot2::theme_classic(base_size = base_size) +
      ggplot2::scale_y_discrete(guide = ggplot2::guide_axis(angle = srt)) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 0, vjust = 0, size = 12 * label.cex),
        axis.text.y = ggplot2::element_text(angle = 0, hjust = 0, size = 12 * label.cex),
        axis.title.x = ggplot2::element_text(size = 12 * axis.cex),
        axis.title.y = ggplot2::element_text(size = 12 * axis.cex)
      )
  } else {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = value, fill = y)) +
      ggplot2::geom_bar(
        stat = "identity", color = "black", size = 0.3,
        width = bar_width, position = posmode
      )
    p <- p + ggplot2::xlab(xlab) +
      ggplot2::ylab(ylab) +
      ggplot2::labs(fill = group.name) +
      ggplot2::ggtitle(main) +
      ggplot2::scale_fill_manual(values = colors) +
      ggplot2::theme_classic(base_size = base_size) +
      ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(angle = srt)) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = srt, vjust = 0, size = 10 * label.cex),
        axis.text.y = ggplot2::element_text(angle = 0, vjust = 0, size = 9 * label.cex),
        axis.title.x = ggplot2::element_text(size = 12 * axis.cex),
        axis.title.y = ggplot2::element_text(size = 12 * axis.cex)
      )
  }

  if (horiz) {
    legend.pos[1] <- 1 - legend.pos[1]
  }

  p <- p + ggplot2::theme(
    legend.title = ggplot2::element_blank(),
    legend.justification = legend.pos,
    legend.text = ggplot2::element_text(size = 9 * legend.cex),
    legend.position = legend.pos,
    legend.key.size = grid::unit(7 * legend.cex, "pt"),
    legend.key.height = grid::unit(5 * legend.cex, "pt")
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
    if (jitter > 0) x <- x + jitter * diff(range(x, na.rm = TRUE)) * stats::rnorm(length(x))


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
      y = graphics::par("usr")[3] - 0.03 * diff(range(x, na.rm = TRUE)),
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
                                   legend.pos = "bottomright", lab.pos = NULL, repel = TRUE,
                                   xlab = NULL, ylab = NULL, xlim = NULL, ylim = NULL, dlim = 0.02,
                                   hilight = NULL, hilight.cex = NULL, lab.xpd = TRUE,
                                   hilight.col = NULL, hilight.lwd = 0.2,
                                   hilight2 = hilight, hilight2.lwd = 0.6,
                                   label.clusters = FALSE, cex.clust = 1.5,
                                   tstep = 0.1, rstep = 0.1, na.color = "#AAAAAA44",
                                   tooltip = NULL, theme = NULL, set.par = TRUE,
                                   axt = "s", xaxs = TRUE, yaxs = TRUE, color_up_down = FALSE,
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
    type <- c("numeric", "factor")[1 + class(var) %in% c("factor", "character", "logical")]
  }
  if (is.null(colnames(pos))) {
    colnames(pos) <- c("x", "y")
  }
  var <- var[match(rownames(pos), names(var))]

  ## x/y limits
  if (!is.null(xlim)) {
    xlim0 <- xlim
  } else {
    xlim0 <- range(pos[, 1], na.rm = TRUE)
    if (zoom != 1) {
      cx <- mean(range(pos[, 1], na.rm = TRUE))
      dx <- diff(range(pos[, 1], na.rm = TRUE))
      xlim0 <- cx + 0.5 * c(-1, 1.05) * dx / zoom
    }
  }
  if (!is.null(ylim)) {
    ylim0 <- ylim
  } else {
    ylim0 <- range(pos[, 2], na.rm = TRUE)
    if (zoom != 1) {
      cy <- mean(range(pos[, 2], na.rm = TRUE))
      dy <- diff(range(pos[, 2], na.rm = TRUE))
      ylim0 <- cy + 0.5 * c(-1, 1.05) * dy / zoom
    }
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
      col <- c(
        RColorBrewer::brewer.pal(8, "Set1"),
        RColorBrewer::brewer.pal(8, "Set2"),
        RColorBrewer::brewer.pal(12, "Set3")
      ) # omics_pal_d("dark")(8))
    } else if (is.null(col) && nz == 2) {
      col <- rev(grDevices::grey.colors(2, end = 0.8))
      col <- c("#888888AA", "#AA0000FF", "#0000AAFF") ## grey/red
    } else if (is.null(col) && nz == 1) {
      col <- c("#22222255")
    }
    col1 <- head(rep(col, 99), nz)

    pt.col <- col1[z1]
    if (color_up_down && nz == 2 && length(col) == 3) {
      sel <- which(as.integer(z1) == 2 & pos[, 1] < 0)
      pt.col[sel] <- col[3]
    }

    ## NA points with light grey
    pt.col[which(is.na(pt.col))] <- na.color

    ## also dim not highlighted points
    if (!is.null(hilight) && length(hilight) > 0) {
      jj <- which(!rownames(pos) %in% hilight)
      if (length(jj)) pt.col[jj] <- na.color
    }

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
      graphics::legend(
        legend.pos,
        legend = levels(z1),
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
    z1 <- (z - min(zlim, na.rm = TRUE)) / diff(zlim)
    z1 <- pmin(pmax(z1, 0), 1) ## clip
    if (softmax) {
      z1 <- 0.5 * (tanh(4 * (z1 - 0.5)) + 1)
    }

    ## -------------- set colors

    if (!is.null(col)) {
      cpal <- grDevices::colorRampPalette(col)(64)
    } else {
      cpal <- grDevices::colorRampPalette(c(omics_colors("brand_blue"), omics_colors("grey"), omics_colors("red")))(64)
    }
    cpal <- sapply(1:length(cpal), function(i) add_opacity(cpal[i], 0.2 + 0.8 * abs(i - 32.5) / 32))

    pt.col <- cpal[1 + ceiling(z1 * (length(cpal) - 1))]


    pt.col0 <- pt.col
    if (opacity < 1) {
      pt.col <- add_opacity(pt.col, opacity)
      cpal <- add_opacity(cpal, opacity**0.33)
    }
    pt.col[is.na(pt.col)] <- na.color ## missing values color

    ## dim not highlighted points
    if (!is.null(hilight) && length(hilight) > 0) {
      jj <- which(!rownames(pos) %in% hilight)
      if (length(jj)) pt.col[jj] <- na.color
    }

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
        legend = c(zr[2], rep(NA, 9), zr[1])
      )
    }
  }

  ## hightlight points in hilight with a border
  hcex <- 0.9 * hilight.cex
  if (length(hcex) > 1) hcex <- hcex[jj]
  if (!is.null(hilight) && length(hilight) > 0) {
    jj <- which(rownames(pos) %in% hilight)
    if (length(jj)) {
      hcol1 <- hilight.col
      if (is.null(hcol1)) hcol1 <- pt.col0[jj]
      graphics::points(pos[jj, , drop = FALSE], pch = 20, col = hcol1, cex = 1.05 * hcex)
      graphics::points(pos[jj, , drop = FALSE], pch = 1, lwd = hilight.lwd, cex = 0.8 * hcex)
    }
  }

  ## labels points in hilight2 (should be renamed...IK
  if (!is.null(hilight2) && length(hilight2) > 0) {
    jj <- which(rownames(pos) %in% hilight2)
    if (length(jj)) {
      if (length(cex.lab) == 1) cex.lab <- rep(cex.lab, nrow(pos))
      df <- data.frame(x = pos[jj, 1], y = pos[jj, 2], z = labels[jj], cex = cex.lab[jj])
      if (repel) {
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
        graphics::points(df$x, df$y, pch = 1, lwd = hilight2.lwd, cex = 0.8 * hcex)
        graphics::segments(df$x, df$y, lab.pos$x, lab.pos$y, col = "#222222AA", lwd = 0.85)
        graphics::text(lab.pos$x, lab.pos$y, labels = df$z, cex = 0.7 * df$cex)
      } else {
        boxes <- sapply(0.8 * nchar(df$z), function(n) paste(rep("\U2588", n), collapse = ""))
        cex1 <- 0.7 * df$cex
        graphics::points(df$x, df$y, pch = 1, lwd = hilight2.lwd, cex = 0.8 * hcex)
        graphics::text(df$x, df$y, labels = boxes, col = "#FFFFFFAA", cex = 1.15 * cex1, pos = 3, offset = 0.5)
        graphics::text(df$x, df$y, labels = df$z, cex = cex1, pos = 3, offset = 0.45)
      }
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
                                     zlim = NULL, cmin = NULL, cmax = NULL,
                                     zlog = FALSE, softmax = FALSE, zsym = FALSE,
                                     xlab = NULL, ylab = NULL, xlim = NULL, ylim = NULL,
                                     hilight2 = hilight, hilight.col = NULL,
                                     hilight.lwd = 0.3, hilight2.lwd = 1, hilight.cex = NULL,
                                     na.color = "#AAAAAA55",
                                     opacity = 1, label.clusters = FALSE, labels = NULL,
                                     legend.ysp = 0.85, legend.pos = "bottomright",
                                     tooltip = NULL, theme = NULL, set.par = TRUE,
                                     label.type = c("text", "box"), base_size = 11,
                                     title = NULL, barscale = 0.8, axis = TRUE, box = TRUE,
                                     guide = "legend", girafe = FALSE) {
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
  if (is.null(zlim) && !is.null(cmin) && !is.null(cmax)) zlim <- c(cmin, cmax)
  
  ## automatically set pointsize of dots
  if (is.null(cex)) {
    nr <- nrow(pos)
    i <- as.integer(cut(nr, breaks = c(0, 100, 500, 1000, 5000, Inf)))
    cex <- c(2, 1.4, 1, 0.7, 0.4)[i]
  }
  if (is.null(hilight.cex)) {
    hilight.cex <- cex
  }

  if(0 && girafe) {
    cex.lab <- 0.9 * cex.lab
    cex.axis <- 0.9 * cex.axis
    cex <- 0.9 * cex
  }
  
  ## normalize pos
  if (is.null(xlim)) xlim <- range(pos[, 1], na.rm = TRUE)
  if (is.null(ylim)) ylim <- range(pos[, 2], na.rm = TRUE)
  if (zoom != 1) {
    cx <- mean(range(pos[, 1], na.rm = TRUE), na.rm = TRUE)
    cy <- mean(range(pos[, 2], na.rm = TRUE), na.rm = TRUE)
    dx <- diff(range(pos[, 1], na.rm = TRUE))
    dy <- diff(range(pos[, 2], na.rm = TRUE))
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

  if(girafe) {
    GEOM_POINT <- ggiraph::geom_point_interactive
  } else {
    GEOM_POINT <- ggplot2::geom_point
  }

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
      x = pos[, 1],
      y = pos[, 2],
      name = rownames(pos),
      text = tooltip,
      label = label1
    )
    jj <- order(-table(pt.col)[pt.col]) ## plot less frequent points last...
    df <- df[jj, ]
    pt.col <- pt.col[jj]
    cex1 <- ifelse(length(cex) > 1, cex[jj], cex)
    x <- y <- NULL    
    plt <- ggplot2::ggplot(
      data = df,
      ggplot2::aes(x = x, y = y, tooltip = text, data_id = name),
      legend = legend
    ) +
      GEOM_POINT(      
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
          legend.key = ggplot2::element_rect(color = "transparent", fill = scales::alpha("white", 0.0)),
          legend.justification = legend.justification,
          legend.position = legend.position,
          legend.background = ggplot2::element_rect(fill = scales::alpha("white", 0.5)),
          legend.margin = ggplot2::margin(0, 4, 4, 4),
          legend.box.just = "right",
          legend.box.background = ggplot2::element_rect(color = "#888888", size = 0.25),
          legend.box.margin = ggplot2::margin(0.8, 1, 1, 1)
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
    cpal <- c(omics_colors("brand_blue"), omics_colors("grey"), omics_colors("red")) # rev(RColorBrewer::brewer.pal(11, "RdYlBu")) ## default
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
      x = pos[, 1],
      y = pos[, 2],
      name = rownames(pos),
      variable = z,
      text = tooltip,
      label = label1
    )
    jj <- order(abs(z), na.last = FALSE)
    df <- df[jj, ] ## strongest last??
    cex1 <- ifelse(length(cex) > 1, cex[jj], cex)

    ## determine range for colorbar
    zr <- range(z, na.rm = TRUE)
    if (zsym && min(zr, na.rm = TRUE) < 0) zr <- c(-1, 1) * max(abs(zr), na.rm = TRUE)
    zz <- round(c(zr[1], zr[2]), digits = 2)
    ##    variable <- NULL
    if (is.null(zlim)) {
      cmin0 <- min(z, na.rm = TRUE)
      cmax0 <- max(z, na.rm = TRUE)
    } else {
      cmin0 <- zlim[1]
      cmax0 <- zlim[2]
    }

    plt <- ggplot2::ggplot(
      df,
      ggplot2::aes(x = x, y = y, fill = variable, tooltip = text, data_id = name),      
    ) +
      GEOM_POINT(
        shape = 21,
        alpha = opacity,
        size = 2.0 * cex1,
        color = "#444444",
        stroke = 0.2
      ) +
      ggplot2::scale_fill_gradientn(
        colors = cpal,
        limits = c(cmin0, cmax0),
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
  } ## end-of-is-numerice

  if (!is.null(hilight)) {
    ## this hilights some points (with color and size)  at certain positions
    sel <- which(df$name %in% hilight)
    if (is.null(hilight.col)) hilight.col <- "transparent"
    plt <- plt +
      GEOM_POINT(
        data = df[sel, ],
        mapping = ggplot2::aes(x, y),
        size = 2.1 * hilight.cex,
        shape = 21,
        stroke = 0.5 * hilight.lwd,
        fill = hilight.col,
        color = "grey20"
      )
  }

  if (!is.null(hilight2)) {
    ## this put text labels at certain positions with geom_text_repel(
    if (label.type == "text") labelFUN <- ggrepel::geom_text_repel
    if (label.type == "box") labelFUN <- ggrepel::geom_label_repel
    plt <- plt +
      GEOM_POINT(
        data = subset(df, name %in% hilight2),
        size = 2.1 * hilight.cex,
        shape = 21,
        stroke = 0.5 * hilight2.lwd,
        fill = hilight.col,
        color = "grey20"
      ) +
      labelFUN(
        data = subset(df, name %in% hilight2),
        ggplot2::aes(label = label),
        size = 5.0 * cex.lab,
        color = "black",
        max.overlaps = 99,
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
      axis.text.x = ggplot2::element_text(size = 12 * cex.axis),
      axis.text.y = ggplot2::element_text(size = 12 * cex.axis),
      axis.title.x = ggplot2::element_text(size = 18 * cex.axis, vjust = -3),
      axis.title.y = ggplot2::element_text(size = 18 * cex.axis, vjust = +5),
      plot.margin = ggplot2::margin(1, 1, 10, 10, "mm") ## ??
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

  if (girafe) {
    plt <- ggiraph::girafe(ggobj = plt, pointsize = 10, width_svg = NULL, height_svg = NULL)
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
                                     var = NULL, type = NULL, col = NULL, na.color = "#AAAAAA55",
                                     cex = NULL, cex.lab = 0.8, cex.title = 1.2,
                                     cex.clust = 1.5, cex.legend = 1, cex.axis = 1,
                                     xlab = NULL, ylab = NULL, xlim = NULL, ylim = NULL,
                                     axis = TRUE, zoom = 1, legend = TRUE, bty = "n",
                                     hilight = NULL, hilight2 = hilight, labels = rownames(pos),
                                     hilight.col = NULL, hilight.cex = NULL,
                                     hilight.lwd = 0.4, hilight2.lwd = 0.8,
                                     zlim = NULL, zlog = FALSE, zsym = FALSE, softmax = FALSE,
                                     opc.low = 1, opacity = 1, bgcolor = NULL, box = TRUE,
                                     label.clusters = FALSE, label.type = NULL,
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
  xlim0 <- range(pos[, 1], na.rm = TRUE)
  ylim0 <- range(pos[, 2], na.rm = TRUE)
  if (zoom != 1) {
    cx <- mean(range(pos[, 1], na.rm = TRUE), na.rm = TRUE)
    cy <- mean(range(pos[, 2], na.rm = TRUE), na.rm = TRUE)
    dx <- diff(range(pos[, 1], na.rm = TRUE))
    dy <- diff(range(pos[, 2], na.rm = TRUE))
    xlim0 <- cx + 0.5 * c(-1, 1.05) * dx / zoom
    ylim0 <- cy + 0.5 * c(-1, 1.05) * dy / zoom
  }
  if (is.null(xlab)) xlab <- colnames(pos)[1]
  if (is.null(ylab)) ylab <- colnames(pos)[2]

  if (type == "numeric") var <- round(var, digits = 4)
  hoverinfo <- "text"
  tooltip1 <- paste0(
    "<b>", rownames(pos), "</b>",
    "<br>value: <b>", var,
    "</b><br>x: <b>", round(pos[, 1], 2), "</b><br>y: <b>", round(pos[, 2], 2), "</b>"
  )
  if (!is.null(tooltip) && length(tooltip) == length(var)) {
    tooltip1 <- paste0(tooltip1, "<br>", tooltip)
  }
  if (!is.null(tooltip) && tooltip == FALSE) {
    tooltip1 <- NA
    hoverinfo <- "none"
  }

  if (is.null(labels)) {
    label1 <- rownames(pos)
  } else {
    label1 <- labels
  }

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
      size = 5 * cex,
      text = tooltip1,
      label = label1
    )

    ## plot less frequent points first... (NEED RETHINK)
    jj <- order(-table(z1)[z1])
    df <- df[jj, , drop = FALSE]
  }

  ## Plot the continous variables
  if (type == "numeric") {
    z <- as.numeric(var)
    z1 <- NULL
    if (is.null(col)) {
      cpal <- c(omics_colors("brand_blue"), omics_colors("grey"), omics_colors("red")) # rev(RColorBrewer::brewer.pal(11, "RdYlBu"))
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
      size = 5 * cex,
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

  ## plot NA values as light grey
  any(is.na(df$value))
  if (any(is.na(df$value))) {
    jj <- which(is.na(df$value))
    plt <- plt %>%
      plotly::add_markers(
        data = df[jj, , drop = FALSE],
        x = ~x,
        y = ~y,
        colors = cpal,
        text = ~text,
        hoverinfo = hoverinfo,
        marker = list(
          size = ~size,
          opacity = opacity,
          color = "#DDDDDD44",
          line = list(
            color = "#AAAAAA44",
            width = 0.2
          )
        ),
        showlegend = FALSE,
        key = ~name,
        mode = "markers",
        type = "scattergl"
      )
  }
  ## plot not missing values
  jj <- which(!is.na(df$value))
  pt.opacity <- 1
  if (!is.null(hilight)) {
    jj <- which(!is.na(df$value) & !rownames(df) %in% hilight)
    pt.opacity <- opc.low
  }
  plt <- plt %>%
    plotly::add_markers(
      data = df[jj, , drop = FALSE],
      x = ~x,
      y = ~y,
      color = ~value,
      colors = cpal,
      text = ~text,
      hoverinfo = hoverinfo,
      marker = list(
        size = ~size,
        opacity = opacity * pt.opacity,
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

  ## plot hilighted points
  if (!is.null(hilight)) {
    jj <- which(rownames(df) %in% hilight)
    col1 <- "transparent"
    if (!is.null(hilight.col)) col1 <- hilight.col
    plt <- plt %>%
      ## add_trace(
      plotly::add_markers(
        data = df[jj, ],
        x = ~x, y = ~y,
        color = ~value,
        colors = cpal,
        color = NULL,
        key = ~name,
        mode = "markers",
        type = "scattergl", #
        text = ~text,
        hoverinfo = hoverinfo,
        marker = list(
          ## color = col1,
          opacity = 1,
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

  ## plot hilighted points with label
  if (!is.null(hilight2)) {
    jj <- which(rownames(df) %in% hilight2)
    plt <- plt %>%
      plotly::add_markers(
        data = df[jj, ],
        x = ~x,
        y = ~y,
        color = ~value,
        colors = cpal,
        key = ~name,
        mode = "markers",
        type = "scattergl", #
        ## text = ~text,
        ## hoverinfo = hoverinfo,
        marker = list(
          ## color = col1,
          opacity = 1,
          size = 5 * hilight.cex,
          showlegend = FALSE,
          showscale = FALSE,
          line = list(
            color = "#000000",
            width = 1.0 * hilight2.lwd
          )
        )
      ) %>%
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

  ## cluster labels
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
        x = mpos[, 1],
        y = mpos[, 2],
        showarrow = FALSE,
        text = paste0("<b>", mlab, "</b>"),
        font = list(size = 16 * cex.clust),
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
pgx.scatterPlotXY.D3 <- function(pos, var = NULL, type = NULL, col = NULL, na.color = "#AAAAAA55",
                                 cex = 1, cex.lab = 0.8, cex.title = 1.2, cex.clust = 1.5, cex.legend = 1,
                                 zoom = 1, legend = TRUE, bty = "n", hilight = NULL,
                                 zlim = NULL, zlog = FALSE, softmax = FALSE,
                                 xlab = NULL, ylab = NULL, hilight2 = hilight,
                                 hilight.lwd = 0, hilight2.lwd = 0, zsym = TRUE,
                                 xlim = NULL, ylim = NULL,
                                 opacity = 1, label.clusters = FALSE, labels = NULL,
                                 legend.ysp = 0.85, legend.pos = "bottomright",
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
  df <- data.frame(
    x = pos[, 1],
    y = pos[, 2],
    z = var,
    names = rownames(pos)
  )
  if (!is.null(var)) {
    plt <- scatterD3::scatterD3(
      data = df,
      x = x,
      y = y,
      col_var = z,
      xlab = xlab,
      ylab = ylab,
      point_size = 32 * cex,
      legend_width = 70,
      col_lab = "value"
    )
  } else {
    plt <- scatterD3::scatterD3(
      data = df,
      x = x,
      y = y, #
      xlab = xlab,
      ylab = ylab,
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
    color = ~variable,
    hovertemplate = "<b>%{x}</b><br>%{fullData.name}: <b>%{y}</b><br><extra></extra>"
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
#' @param color_up_down Color up/down regulated features
#' @param marker.type Marker type (scatter, line, etc). Default "scatter".
#' @param displayModeBar Show plotly modebar? Logical. Default TRUE.
#'
#' @return A plotly interactive MA plot object.
#'
#' @export
plotlyMA <- function(x, y, names, label.names = names,
                     group.names = c("group1", "group2"),
                     xlab = "average expression (log2.CPM)",
                     ylab = "effect size (log2.FC)",
                     lfc = 1, psig = 0.05, showlegend = TRUE, highlight = NULL,
                     marker.size = 5, label = NULL, label.cex = 1,
                     color_up_down = TRUE,
                     colors = c(up = "#f23451", notsig = "#8F8F8F", down = "#3181de"),
                     marker.type = "scatter", source = "plot1",
                     displayModeBar = TRUE) {
  if (is.null(highlight)) highlight <- names
  i0 <- which(!names %in% highlight & !label.names %in% highlight)
  i1 <- which(names %in% highlight | label.names %in% highlight)

  p <- plotly::plot_ly(
    hovertemplate = "<b>%{text}</b><br><b>Average expression</b>: %{x:.2f}<br><b>Effect size</b>: %{y:.2f}<extra></extra>"
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
        type = marker.type,
        mode = "markers",
        marker = list(
          size = marker.size,
          color = "#ccc"
        ),
        showlegend = showlegend
      )
  }

  if (length(i1)) {
    if (color_up_down) {
      upreg <- y[i1] > 0
      dwreg <- y[i1] < 0
      p <- p %>%
        plotly::add_trace(
          x = x[i1][upreg],
          y = y[i1][upreg],
          text = names[i1][upreg],
          type = marker.type,
          mode = "markers",
          marker = list(
            size = marker.size,
            color = colors["up"]
          ),
          showlegend = showlegend
        )
      p <- p %>%
        plotly::add_trace(
          x = x[i1][dwreg],
          y = y[i1][dwreg],
          text = names[i1][dwreg],
          type = marker.type,
          mode = "markers",
          marker = list(
            size = marker.size,
            color = colors["down"]
          ),
          showlegend = showlegend
        )
    } else {
      p <- p %>%
        plotly::add_trace(
          x = x[i1],
          y = y[i1],
          text = names[i1],
          type = marker.type,
          mode = "markers",
          marker = list(
            size = marker.size,
            color = colors["notsig"]
          ),
          showlegend = showlegend
        )
    }
  }

  if (!is.null(label) && length(label) > 0) {
    i2 <- which(names %in% label | label.names %in% label)
    if (color_up_down) {
      upreg <- y[i2] > 0
      dwreg <- y[i2] < 0
      annot_text <- label.names[i2][upreg]
      if (length(annot_text) == 0) annot_text <- ""
      p <- p %>%
        plotly::add_annotations(
          x = x[i2][upreg],
          y = y[i2][upreg],
          text = annot_text,
          font = list(
            size = 12 * label.cex,
            color = colors["up"]
          ),
          showarrow = FALSE,
          yanchor = "bottom",
          yshift = 2,
          textposition = "top"
        )
      annot_text <- label.names[i2][dwreg]
      if (length(annot_text) == 0) annot_text <- ""
      p <- p %>%
        plotly::add_annotations(
          x = x[i2][dwreg],
          y = y[i2][dwreg],
          text = annot_text,
          font = list(
            size = 12 * label.cex,
            color = colors["down"]
          ),
          showarrow = FALSE,
          yanchor = "bottom",
          yshift = 2,
          textposition = "top"
        )
    } else {
      p <- p %>%
        plotly::add_annotations(
          x = x[i2],
          y = y[i2],
          text = label.names[i2],
          font = list(
            size = 12 * label.cex,
            color = colors["notsig"]
          ),
          showarrow = FALSE,
          yanchor = "bottom",
          yshift = 2,
          textposition = "top"
        )
    }
  }

  x1 <- 1.05 * max(x, na.rm = TRUE)
  yy <- 1.05 * max(abs(y), na.rm = TRUE)
  abline1 <- list(
    type = "line", y0 = -lfc, y1 = -lfc, x0 = 0, x1 = x1,
    line = list(dash = "dot", width = 1, color = "grey")
  )
  abline2 <- list(
    type = "line", y0 = +lfc, y1 = +lfc, x0 = 0, x1 = x1,
    line = list(dash = "dot", width = 1, color = "grey")
  )

  xrange <- c(0, 1) * max(abs(x), na.rm = TRUE) * 1.05
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
    ) %>%
    plotly_build_light(.)

  return(p)
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
#' @param color_up_down Color up/down regulated features
#' @param marker.type Marker type (scatter, line, etc). Default "scatter".
#' @param displayModeBar Show plotly modebar? Logical. Default TRUE.
#'
#' @return A plotly interactive volcano plot object.
#'
#' @export
plotlyVolcano <- function(x,
                          y,
                          names,
                          label.names = names,
                          group.names = c("group1", "group2"),
                          xlab = "effect size (logFC)",
                          ylab = "significance (-log10p)",
                          lfc = 1,
                          psig = 0.05,
                          showlegend = TRUE,
                          highlight = NULL,
                          marker.size = 5,
                          label = NULL,
                          label.cex = 1,
                          max.absy = NULL,
                          color_up_down = TRUE,
                          colors = c(
                            up = "#f23451", notsig = "#8F8F8F",
                            down = "#3181de", notsel = "#cccccc88"
                          ),
                          title = "Volcano plot",
                          marker.type = "scatter",
                          displayModeBar = TRUE,
                          source = "plot1") {
  if (is.null(highlight)) highlight <- names

  i0 <- which(!names %in% highlight & !label.names %in% highlight)
  i1 <- which(names %in% highlight | label.names %in% highlight)

  # Detect wich i1 genes are under the thresholds
  notsig.genes <- which(y <= -log10(psig) | abs(x) < lfc)
  ib <- intersect(notsig.genes, i1)
  i1 <- setdiff(i1, ib)

  p <- plotly::plot_ly(
    source = source,
    hovertemplate = "<b>%{text}</b><br><b>Effect size</b>: %{x:.2f}<br><b>Significance</b>: %{y:.2f}<extra></extra>"
  )
  p <- p %>%
    plotly::event_register("plotly_hover") %>%
    plotly::event_register("plotly_selected")

  ## not in highlight, not selected
  if (length(i0)) {
    p <- p %>%
      plotly::add_trace(
        x = x[i0],
        y = y[i0],
        text = names[i0],
        type = marker.type,
        mode = "markers",
        marker = list(
          size = marker.size,
          color = colors["notsel"]
        ),
        showlegend = FALSE
      )
  }


  ## highlighted, selected
  if (length(i1)) {
    if (color_up_down) {
      upreg <- x[i1] > 0
      dwreg <- x[i1] < 0
      p <- p %>%
        plotly::add_trace(
          x = x[i1][upreg],
          y = y[i1][upreg],
          text = names[i1][upreg],
          type = marker.type,
          mode = "markers",
          marker = list(
            size = marker.size,
            color = colors["up"]
          ),
          showlegend = showlegend,
          name = "Up"
        )
      p <- p %>%
        plotly::add_trace(
          x = x[i1][dwreg],
          y = y[i1][dwreg],
          text = names[i1][dwreg],
          type = marker.type,
          mode = "markers",
          marker = list(
            size = marker.size,
            color = colors["down"]
          ),
          showlegend = showlegend,
          name = "Down"
        )
    } else {
      p <- p %>%
        plotly::add_trace(
          x = x[i1],
          y = y[i1],
          text = names[i1],
          type = marker.type,
          mode = "markers",
          marker = list(
            size = marker.size,
            color = "black"
          ),
          showlegend = showlegend
        )
    }
  }

  ## highlighted but not significant
  if (length(ib)) {
    p <- p %>%
      plotly::add_trace(
        x = x[ib],
        y = y[ib],
        text = names[ib],
        type = marker.type,
        mode = "markers",
        marker = list(
          size = marker.size,
          color = colors["notsig"]
        ),
        showlegend = showlegend,
        name = "Not significant"
      )
  }

  ## label those points that are in 'label' vector. It should match
  ## either the names or label.names
  if (!is.null(label) && length(label) > 0) {
    i2 <- which(names %in% label | label.names %in% label)
    i2 <- i2[!i2 %in% ib]
    upreg <- x[i2] > 0
    dwreg <- x[i2] < 0
    if (TRUE && color_up_down) {
      annot_text <- label.names[i2][upreg]
      if (length(annot_text) == 0) annot_text <- ""
      p <- p %>%
        plotly::add_annotations(
          x = x[i2][upreg],
          y = y[i2][upreg],
          text = annot_text,
          font = list(
            size = 12 * label.cex,
            color = colors["up"]
          ),
          showarrow = FALSE,
          yanchor = "bottom",
          yshift = 2,
          textposition = "top"
        )
      annot_text <- label.names[i2][dwreg]
      if (length(annot_text) == 0) annot_text <- ""
      p <- p %>%
        plotly::add_annotations(
          x = x[i2][dwreg],
          y = y[i2][dwreg],
          text = annot_text,
          font = list(
            size = 12 * label.cex,
            color = colors["down"]
          ),
          showarrow = FALSE,
          yanchor = "bottom",
          yshift = 2,
          textposition = "top"
        )
    } else {
      p <- p %>%
        plotly::add_annotations(
          x = x[i2],
          y = y[i2],
          text = label.names[i2],
          font = list(
            size = 12 * label.cex,
            color = "black"
          ),
          showarrow = FALSE,
          yanchor = "bottom",
          yshift = 2,
          textposition = "top"
        )
    }
    idl <- which(names[ib] %in% label)
    if (length(idl)) {
      p <- p %>%
        plotly::add_annotations(
          x = x[ib][idl],
          y = y[ib][idl],
          text = label.names[ib][idl],
          font = list(
            size = 12 * label.cex,
            color = colors["notsig"]
          ),
          showarrow = FALSE,
          yanchor = "bottom",
          yshift = 2,
          textposition = "top"
        )
    }
  } ## end of if-label

  ## draw vertical logFC threshold lines at specified logFC
  y0 <- -log10(psig)
  y1 <- 0.95 * max(y, na.rm = TRUE)
  xx <- 1.05 * max(abs(x), na.rm = TRUE)
  abline1 <- list(
    type = "line", x0 = -lfc, x1 = -lfc, y0 = 0, y1 = y1,
    line = list(dash = "dot", width = 1, color = "grey")
  )
  abline2 <- list(
    type = "line", x0 = +lfc, x1 = +lfc, y0 = 0, y1 = y1,
    line = list(dash = "dot", width = 1, color = "grey")
  )
  ## draw horizontal significance lines at p=0.05
  abline3 <- list(
    type = "line", x0 = -xx, x1 = +xx, y0 = y0, y1 = y0,
    line = list(dash = "dot", width = 1, color = "grey")
  )
  if (lfc == 0) {
    significance.lines <- list(abline3)
  } else {
    significance.lines <- list(abline1, abline2, abline3)
  }

  ## calculate explicit ranges for x and y so we can add some padding.
  max.absx <- max(max(abs(x), na.rm = TRUE), lfc * 1.2, na.rm = TRUE)
  if (is.null(max.absy)) {
    max.absy <- max(max(abs(y), na.rm = TRUE), y0 * 1.2, na.rm = TRUE)
  }
  xrange <- c(-1, 1) * max.absx * 1.05
  if (min(x, na.rm = TRUE) >= 0) xrange <- c(0, 1) * max.absx * 1.05
  yrange <- c(0, 1) * max.absy * 1.05
  xaxis <- list(title = xlab, range = xrange, showgrid = FALSE)
  yaxis <- list(title = list(text = ylab, standoff = 20L), range = yrange, showgrid = FALSE)

  p <- p %>%
    plotly::layout(
      shapes = significance.lines,
      xaxis = xaxis,
      yaxis = yaxis,
      showlegend = showlegend,
      hovermode = "closest",
      dragmode = "select"
    ) %>%
    plotly::layout(
      margin = list(l = 0, b = 1, t = 10, r = 10),
      font = list(size = 12),
      legend = list(
        font = list(size = 12)
      )
    ) %>%
    plotly_build_light(.)
  return(p)
}


#' @export
plotlyVolcano_multi <- function(FC,
                                Q,
                                by_sig = TRUE,
                                fdr = 0.05,
                                lfc = 0,
                                label = NULL,
                                share_axis = FALSE,
                                title_y = "ylab",
                                title_x = "xlab",
                                cex = 3, ## marker size
                                label.cex = 1, ## label size
                                yrange = 0.5,
                                n_rows = 2,
                                margin_l = 45,
                                margin_b = 60,
                                title_y_offset = -0.025,
                                title_x_offset = -0.12,
                                interplot_margin = c(0.01, 0.0, 0.05, 0.1),
                                annotation_args = list(),
                                layout_args = list(),
                                highlight = NULL,
                                pval_cap = 1e-12,
                                ...) {
  dots <- list(...)

  ## Get tables and genes
  fc <- as.matrix(FC)
  qv <- as.matrix(Q)
  titles <- colnames(fc)

  all_genes <- rownames(FC)
  all_labels <- all_genes
  if ("names" %in% names(dots)) all_genes <- dots$names
  if ("label.names" %in% names(dots)) all_labels <- dots$label.names

  # Prepare collection list
  nplots <- min(24, length(titles))
  sub_plots <- vector("list", length = length(nplots))
  for (i in 1:nplots) {
    # Input variables
    fx <- fc[, i, drop = FALSE]
    qval <- qv[, i, drop = FALSE]
    title_i <- titles[i]

    # Set marker colour group
    hi.genes <- NULL
    if (by_sig) {
      is.sig <- (qval <= fdr & abs(fx) >= lfc)
      hi.genes <- all_genes[is.sig]
    }
    if (!is.null(highlight)) {
      hi.genes1 <- all_genes[all_genes %in% highlight | all_labels %in% highlight]
      if (!is.null(hi.genes)) {
        hi.genes <- intersect(hi.genes, hi.genes1)
      } else {
        hi.genes <- hi.genes1
      }
    }

    # Take -log and add 1e-12 to remove 0, and avoid Infs
    qval <- -log10(qval + pval_cap)

    # Set labels
    sub.label <- label
    if (!is.null(label) && is.list(label)) {
      sub.label <- label[[i]]
    }

    # Set title
    if (share_axis) {
      title_loc <- -log10(min(qv + pval_cap, na.rm = TRUE))
      title_loc <- title_loc + title_loc * (yrange / 20)
    } else {
      title_loc <- max(qval, na.rm = TRUE) + max(qval, na.rm = TRUE) * (yrange / 20)
    }
    # Call volcano plot
    sub_plots[[i]] <- plotlyVolcano(
      x = fx,
      y = qval,
      label = sub.label,
      marker.type = "scattergl",
      marker.size = cex,
      highlight = hi.genes,
      label.cex = label.cex,
      group.names = c("group1", "group0"),
      psig = fdr,
      lfc = lfc,
      max.absy = title_loc + title_loc * (yrange / 10),
      showlegend = FALSE,
      ...
      # Add plot title
    ) %>%
      plotly::add_annotations(
        text = paste("<b>", title_i, "</b>"),
        font = list(size = 15),
        showarrow = FALSE,
        bgcolor = "white",
        xanchor = "centre",
        yanchor = "bottom",
        x = 0,
        y = title_loc
      ) %>%
      plotly_build_light(.)
  }

  # Pass argument scale_per_plot to subplot
  shareY <- shareX <- ifelse(share_axis, TRUE, FALSE)

  # Arrange subplots
  if (is.null(n_rows)) {
    n_rows <- floor(sqrt(nplots))
  }
  if (nplots <= 2) {
    n_rows <- 1
  }
  suppressWarnings(
    all_plts <- plotly::subplot(sub_plots,
      nrows = n_rows,
      margin = interplot_margin,
      titleY = FALSE,
      titleX = FALSE,
      shareX = shareX,
      shareY = shareY
    ) %>%
      # Add common axis titles
      plotly::layout(
        annotations = modifyList(list(
          list(
            x = title_y_offset,
            y = 0.5,
            text = title_y,
            font = list(size = 13),
            textangle = 270,
            showarrow = FALSE,
            xref = "paper",
            yref = "paper"
          ),
          list(
            x = 0.5,
            y = title_x_offset,
            text = title_x,
            font = list(size = 13),
            showarrow = FALSE,
            xref = "paper",
            yref = "paper"
          )
        ), annotation_args),
        margin = modifyList(
          list(
            l = margin_l,
            b = margin_b
          ),
          layout_args
        )
      )
  )

  return(all_plts)
}


#' Build \code{plotly} data with low computation cost
#'
#' @description
#' Before illustrating data using \code{plotly}, it must be built
#' (\code{figure$x$data} are need to be made using \code{figure$x$attrs}).
#' However, because a lot of procedures are necessary,
#' the computation cost is relatively high.
#' With this function, the data is built in quite short time by omitting
#' several procedures for high-frequency data.
#' Note that this function is not universally applicable to all \code{plotly}
#' objects but made for high-frequency scatter data.
#' \code{plotly::plotly_build} function may return better results in
#' specific cases although it takes more time.
#' @param fig \code{plotly} object.
#' Note that \code{fig$x$attrs} is not \code{NULL} and
#' each \code{fig$x$attrs} element has an element named \code{x}.
#' This function generates \code{fig$x$data} using \code{fig$x$attrs}.
#' @param vars_hf Character, optional.
#' Variable names where high frequency data is included.
#' It must include \code{x}.
#' @return built \code{plotly} object
#' @examples
#' data(noise_fluct)
#' plotly_build_light(
#'   plotly::plot_ly(
#'     x = noise_fluct$time,
#'     y = noise_fluct$f500,
#'     name = "level",
#'     type = "scatter"
#'   )
#' )
#'
#' plotly_build_light(
#'   plotly::plot_ly(
#'     data = noise_fluct,
#'     x = ~time,
#'     y = ~f500,
#'     name = "level",
#'     type = "scatter"
#'   )
#' )
#' @export
plotly_build_light <- function(
    fig, vars_hf = c("x", "y", "text", "hovertext")) {
  # check_arguments
  stopifnot(inherits(fig, "plotly"))
  stopifnot(inherits(fig$x$attrs, "list"))
  stopifnot(!is.null(names(fig$x$attrs)))

  # just do plotly_build if the data is not large
  n_x <- fig$x$attrs %>%
    purrr::discard(~ is.null(.x$type) || is.na(.x$type)) %>%
    purrr::map_int(~ length(.x$x)) %>%
    max()

  if (n_x > 1e3) {
    # evaluate the trace, if necessary
    traces_div <- fig$x$attrs %>%
      purrr::discard(~ is.null(.x$type) || is.na(.x$type)) %>%
      purrr::imodify(
        function(trace, uid) {
          trace_eval <- purrr::modify_if(
            trace,
            lazyeval::is_formula,
            ~ lazyeval::f_eval(.x, plotly::plotly_data(fig, uid))
          )

          attrs_length <- purrr::map_int(trace_eval, length)

          vars_long <- names(trace_eval[attrs_length == attrs_length["x"]])
          vars_long <- intersect(vars_long, vars_hf)

          data_long <- trace_eval[vars_long] %>%
            data.table::setDT() %>%
            .[, lapply(.SD, list)]

          trace_data <- purrr::pmap(
            data_long,
            function(...) {
              c(trace[setdiff(names(trace), vars_long)], list(...))
            }
          )

          return(trace_data)
        }
      ) %>%
      unlist(recursive = FALSE)

    # replace attributes with the ones without high frequency data
    # then build it
    fig$x$attrs <- purrr::map(
      traces_div,
      ~ .x[setdiff(names(.x), vars_hf)]
    )
    fig_built <- plotly::plotly_build(fig)

    # directly input the high frequency data to the plotly data
    fig_built$x$data <- purrr::map2(
      fig_built$x$data, traces_div,
      ~ c(.x, .y[intersect(names(.y), vars_hf)])
    )
  } else {
    fig_built <- plotly::plotly_build(fig)
  }

  return(fig_built)
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
                           nbinsx = 10,
                           nbinsy = 10,
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

  m1 <- mean(x1, na.rm = TRUE)
  m2 <- mean(x2, na.rm = TRUE)

  ## select samples in different quadrants
  j1 <- length(samples[which(x1 < m1 & x2 > m2)]) ## top-left
  j2 <- length(samples[which(x1 > m1 & x2 < m2)]) ## bottom-right
  j3 <- length(samples[which(x1 > m1 & x2 > m2)]) ## top-right
  j4 <- length(samples[which(x1 < m1 & x2 < m2)]) ## bottom-left

  xlab1 <- paste(gene1, lab.unit, collapse = "  ")
  ylab1 <- paste(gene2, lab.unit, collapse = "  ")
  xaxis <- list(
    title = xlab1, range = range(x1, na.rm = TRUE), gridwidth = 0.2, showgrid = TRUE, showline = TRUE,
    zerolinewidth = 0, zerolinecolor = "#fff", autorange = TRUE
  )
  yaxis <- list(
    title = ylab1, range = range(x2, na.rm = TRUE), gridwidth = 0.2, showgrid = TRUE, showline = TRUE,
    zerolinewidth = 0, zerolinecolor = "#fff", autorange = TRUE
  )

  p <- plotly::plot_ly(
    x = x1,
    y = x2,
    text = names(x1),
    hoverinfo = "text",
    type = "scatter",
    mode = "markers",
    marker = list(size = marker.size, color = marker.color),
    showlegend = TRUE
  ) %>%
    plotly::add_trace(
      type = "histogram2dcontour",
      mode = NULL,
      contours = list(coloring = contour.coloring),
      nbinsx = nbinsx,
      nbinsy = nbinsy
    )

  p <- p %>%
    ## plotly::add_trace(
    ##      x = x1,
    ##      y = x2,
    ##      text = names(x1),
    ##      hoverinfo = "text",
    ##      type = "scatter",
    ##      mode = "markers",
    ##      marker = list(size = marker.size, color = "black")
    ## ) %>%
    plotly::layout(
      shapes = list(list(
        type = "line",
        x0 = 0,
        x1 = 1,
        xref = "paper",
        y0 = m2,
        y1 = m2,
        line = list(color = "#bebebe", dash = "dot")
      ), list(
        type = "line",
        y0 = 0,
        y1 = 1,
        yref = "paper",
        x0 = m1,
        x1 = m1,
        line = list(color = "#bebebe", dash = "dot")
      )),
      xaxis = xaxis,
      yaxis = yaxis
    )

  quadrants <- c(j3, j1, j2, j4)
  N <- sum(quadrants)
  positions <- matrix(c(0.98, 0.98, 0.02, 0.98, 0.98, 0.02, 0.02, 0.02), ncol = 2, byrow = TRUE)

  for (i in 1:4) {
    p <- p %>% plotly::add_annotations(
      x = positions[i, 1],
      y = positions[i, 2],
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
      x = pos1[1],
      y = pos1[2],
      text = paste(lab1, collapse = "\n"),
      showarrow = FALSE,
      font = list(size = 15)
    )

    lab2 <- Matrix::head(names(sort(-Matrix::colSums(inferred.celltype[j2, , drop = FALSE]))), 3)
    pos2 <- apply(cbind(x1, x2)[j2, , drop = FALSE], 2, stats::median)
    p <- p %>% plotly::add_annotations(
      x = pos2[1],
      y = pos2[2],
      text = paste(lab2, collapse = "\n"),
      showarrow = FALSE,
      font = list(size = 15)
    )

    lab3 <- Matrix::head(names(sort(-Matrix::colSums(inferred.celltype[j3, , drop = FALSE]))), 3)
    pos3 <- apply(cbind(x1, x2)[j3, , drop = FALSE], 2, stats::median)
    p <- p %>% plotly::add_annotations(
      x = pos3[1],
      y = pos3[2],
      text = paste(lab3, collapse = "\n"),
      showarrow = FALSE,
      font = list(size = 15)
    )

    ## bottom-left??
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
        tmp_colors <- iheatmapr:::pick_continuous_colors(
          zmid = 0,
          zmin = min(x[, i], na.rm = TRUE),
          zmax = max(x[, i], na.rm = TRUE), p
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


#' Split heatmap from matrix using iHeatmapR
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
#' @param colors Vector of colors to use for annotation. Default is RColorBrewer Set1.
#' @param lmar Label margin parameter. Default is 60.
#' @param rowcex Row text size scaling. Default is 1.
#' @param colcex Column text size scaling. Default is 1.
#' @param show_legend Show color legend? Default is TRUE.
#' @param return_x_matrix Return the X matrix of the plot along the plotly object (as a list)
#'
#' @return A plotly split heatmap object
#'
#' @details This function generates an interactive split heatmap from a data matrix.
#' Rows and/or columns can be split according to the input annotations.
#' Data is scaled and clustered using hierarchical clustering by default.
#' Various graphical parameters like colors, labels, text sizes can be adjusted.
#'
#' @export
pgx.splitHeatmapFromMatrix <- function(X, annot = NULL, idx = NULL, splitx = NULL,
                                       xtips = NULL, ytips = NULL, row_clust = TRUE,
                                       row_annot_width = 0.03, scale = "row.center",
                                       colors = NULL, lmar = 60, na_text = NULL,
                                       rowcex = 1, colcex = 1, show_legend = TRUE,
                                       return_x_matrix = FALSE) {
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
    if (nrow(x) == 1) {
      return(rownames(x))
    }
    suppressWarnings(dd <- stats::as.dist(1 - stats::cor(t(x), use = "pairwise")))
    if (sum(is.na(dd))) dd[is.na(dd) | is.nan(dd)] <- 1
    hc <- fastcluster::hclust(dd, method = "ward.D2")
    rownames(x)[hc$order]
  }
  if (!is.null(idx)) {
    if (row_clust) {
      kk <- tapply(rownames(X), idx, function(k) c(hc.order(X[k, , drop = FALSE]), "   "))
    } else {
      kk <- tapply(rownames(X), idx, function(k) c(k, "   "))
    }
    idx <- tapply(idx, idx, function(k) c(k, NA))
    idx <- as.vector(unlist(idx))
    kk <- unlist(kk)
    kk <- kk[1:(length(kk) - 1)] ## remove trailing spacer
    idx <- idx[1:(length(idx) - 1)]
    X <- rbind(X, "   " = NA)[kk, ]

    ## invert
    X <- X[nrow(X):1, ]
    idx <- rev(idx)
  } else {
    if (row_clust) {
      kk <- hc.order(X[, , drop = FALSE])
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
  zmax <- max(abs(X), na.rm = TRUE)

  plt <- iheatmapr::main_heatmap(
    x1,
    name = "expression",
    colorbar_grid = grid_params,
    x = xtips[colnames(x1)],
    y = ytips[rownames(x1)],
    colors = c(omics_colors("brand_blue"), omics_colors("grey"), omics_colors("red")),
    zmid = 0,
    zmin = -zmax,
    zmax = zmax,
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
          colors = c("royalblue3", "#EEEEE4", "indianred3"),
          zmid = 0,
          zmin = -zmax,
          zmax = zmax,
          tooltip = tooltip,
          size = sizes[i],
          buffer = 0.007 * ex
        )

      if (!is.null(hc)) {
        plt <- plt %>%
          iheatmapr::add_col_dendro(hc, size = 0.06)
        ## add_col_clustering() %>%
      }

      if (!is.null(annot)) {
        plt <- plt %>%
          iheatmapr::add_col_title(names(xx)[i], side = "top") %>%
          iheatmapr.add_col_annotation(
            annotation = data.frame(annotF[colnames(x1), , drop = FALSE]),
            size = col_annot_height,
            buffer = 0.005, side = "bottom",
            show_title = FALSE, show_legend = show_legend,
            colors = colors0
          )
      } else {
        plt <- plt %>% iheatmapr::add_col_title(names(xx)[i], side = "top")
      }

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
    cls <- list("gene.module" = c(omics_pal_d("muted_light")(nlevels(as.factor(idx)))))
    if (length(cls[[1]]) == 1) {
      cls[[1]] <- rep(cls[[1]], 2)
    }
    plt <- iheatmapr::add_row_annotation(
      plt, data.frame("gene.module" = as.factor(idx)),
      size = row_annot_width * ex,
      colors = cls,
      show_colorbar = show_legend
    )
  }

  ## ----------- add gene/geneset names
  if (rowcex > 0) {
    gnames <- rownames(X)
    gnames <- gsub("[&].*[;]", "", gnames) ## HTML special garbage...
    #gnames <- gsub("^.*:", "", gnames) ## strip prefix
    gnames <- shortstring(gnames, 40) ## shorten
    gnames <- sub("   ", "-", gnames) ## divider

    maxlen <- max(sapply(gnames, nchar), na.rm = TRUE)
    w <- ifelse(maxlen >= 20, 0.45, 0.20)
    s1 <- ifelse(maxlen >= 20, 9, 11) * rowcex
    plt <- iheatmapr::add_row_labels(
      plt,
      side = "right", ticktext = gnames,
      size = w * ex, font = list(size = s1)
    )
  }

  if (0) {
    ## ---------- Always convert to plotly??? (new)
    plt <- plt %>% iheatmapr::to_plotly_list()
    plt <- plotly::as_widget(plt)
  }

  if (return_x_matrix) {
    return(list(
      plt = plt,
      X = X
    ))
  } else {
    return(plt)
  }
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
    fillcolor = "#A6CEE3",
    linecolor = "#3181de",
    titlecolor = "#1f77b4",
    hoverinfo = "y",
    hoverformat = ".2f",
    yaxistitle = FALSE,
    xaxistitle = FALSE,
    xlen = NULL,
    yrange = NULL,
    barmode = "relative",
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
        data[y],
        list(data[[x]]),
        function(val) {
          c(mean = mean(val, na.rm = TRUE), sd = stats::sd(val, na.rm = TRUE))
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

  if (is.null(fillcolor)) {
    fillcolor <- RColorBrewer::brewer.pal(9, "Set1")
  }
  if (length(fillcolor) == 1) {
    fillcolor <- rep(fillcolor, length(y))
  }


  p <- plotly::plot_ly()
  show_legend <- ifelse(length(y) > 1, TRUE, FALSE)
  for (i in y) {
    p <- p %>% plotly::add_trace(
      type = "bar",
      x = data[[x]],
      y = data[, i],
      marker = list(color = fillcolor[which(i == y)]),
      name = gsub("y.", "", i),
      error_y = error_y,
      hovertext = data[[x]],
      textposition = "none",
      cliponaxis = FALSE,
      hoverinfo = hoverinfo,
      hovertemplate = paste0(
        "<b>%{hovertext}</b><br>",
        "<b>%{yaxis.title.text}:</b> %{y:", hoverformat, "}<br>",
        "<extra></extra>"
      ),
      showlegend = show_legend
    )
  }
  p <- p %>%
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

  if (length(y) > 1) {
    p <- p %>%
      plotly::layout(
        barmode = barmode,
        legend = list(orientation = "h", bgcolor = "transparent", y = 1.2)
      )
  }

  if (length(y) > 1) {
    p <- p %>%
      plotly::layout(
        barmode = barmode,
        legend = list(orientation = "h", bgcolor = "transparent", y = 1.2)
      )
  }

  return(p)
}

#' @export
pgx.plotActivation <- function(pgx,
                               features = NULL,
                               contrasts = NULL,
                               what = "geneset",
                               matrix = NULL,
                               plotlib = "base",
                               filter = NULL,
                               cexCol = 1.4,
                               cexRow = 1,
                               normalize = FALSE,
                               rotate = FALSE,
                               maxterm = 40,
                               maxfc = 100,
                               mar = c(15, 30),
                               tl.cex = 0.85,
                               row.nchar = 60,
                               showscale = TRUE,
                               cexBar = 0.66) {
  if (what == "geneset") {
    score <- pgx.getMetaMatrix(pgx, level = "geneset")$fc
  }
  if (what == "gene") {
    score <- pgx.getMetaMatrix(pgx, level = "gene")$fc
  }
  if (what == "drugs") {
    score <- pgx$drugs[[1]]$X
  }
  if (what == "matrix") {
    if (is.null(matrix)) stop("must provide matrix")
    score <- matrix
  }
  dim(score)

  ## avoid errors!!!
  score[is.na(score) | is.infinite(score)] <- 0

  if (!is.null(contrasts)) {
    score <- score[, contrasts, drop = FALSE]
  }
  if (!is.null(filter)) {
    for (f in filter) {
      sel <- grep(f, rownames(score))
      score <- score[sel, , drop = FALSE]
      rownames(score) <- sub(f, "", rownames(score), ignore.case = TRUE)
    }
  }

  if (!is.null(features)) {
    features <- intersect(features, rownames(score))
    score <- score[features, , drop = FALSE]
  }

  ## max number terms
  score <- score[order(-rowSums(score**2, na.rm = TRUE)), , drop = FALSE]

  ## max number terms
  score <- head(score, maxterm)

  ## max comparisons/FC
  score <- score[, head(order(-colSums(score**2, na.rm = TRUE)), maxfc), drop = FALSE]
  score <- score + 1e-3 * matrix(rnorm(length(score)), nrow(score), ncol(score))
  dim(score)

  ## normalize colums
  if (normalize) {
    ## column scale???
    score <- t(t(score) / (1e-8 + apply(abs(score), 2, max, na.rm = TRUE)))
  }

  if (NCOL(score) == 1) {
    score <- score[order(-score[, 1]), 1, drop = FALSE]
    score <- cbind(score, score)
    colnames(score)[2] <- ""
  } else {
    d1 <- as.dist(1 - cor(t(score), use = "pairwise"))
    d2 <- as.dist(1 - cor(score, use = "pairwise"))
    d1 <- dist(score)
    d2 <- dist(t(score))
    d1[is.na(d1)] <- 1
    d2[is.na(d2)] <- 1
    ii <- 1:nrow(score)
    jj <- 1:ncol(score)
    ii <- hclust(d1)$order
    jj <- hclust(d2)$order
    score <- score[ii, jj, drop = FALSE]
  }

  dim(score)
  colnames(score) <- substring(colnames(score), 1, 30)
  rownames(score) <- substring(rownames(score), 1, row.nchar)
  colnames(score) <- paste0(colnames(score), " ")

  if (rotate) score <- t(score)

  bluered.pal <- colorRamp(colors = c("royalblue3", "#ebeffa", "white", "#faeeee", "indianred3"))
  bluered.pal <- colorRamp(colors = c("royalblue3", "grey90", "indianred3"))
  bluered.pal <- colorRamp(playbase::omics_pal_c("blue_red_grey", reverse = TRUE)(30))

  score <- score[nrow(score):1, , drop = FALSE]
  x_axis <- colnames(score)
  y_axis <- rownames(score)

  fig <- NULL
  if (plotlib == "base") {
    gx.heatmap(
      score,
      dist.method = "euclidean", ## important
      scale = "none", ## important
      mar = mar,
      keysize = 0.4,
      key = FALSE,
      cexRow = cexRow,
      cexCol = cexCol,
      softmax = FALSE
    )
  }
  if (plotlib == "plotly") {
    fig <- plotly::plot_ly(
      x = x_axis,
      y = y_axis,
      z = score,
      zauto = FALSE,
      zmin = -max(abs(score)),
      zmax = +max(abs(score)),
      type = "heatmap",
      colors = bluered.pal,
      showscale = showscale
    ) %>%
      plotly::colorbar(
        len = 0.5 * cexBar,
        thickness = 30 * cexBar
      )
  }
  return(fig)
}


#' @export
visPrint <- function(visnet, file, width = 3000, height = 3000, delay = 0, zoom = 1) {
  is.pdf <- grepl("pdf$", file)
  if (is.pdf) {
    width <- width * 600
    height <- height * 600
  }
  vis2 <- htmlwidgets::createWidget(
    name = "visNetwork",
    x = visnet$x,
    width = width, height = height,
    package = "visNetwork"
  )
  tmp.html <- paste0(tempfile(), "-visnet.html")
  tmp.png <- paste0(tempfile(), "-webshot.png")
  visNetwork::visSave(vis2, file = tmp.html)
  webshot2::webshot(
    url = tmp.html,
    file = tmp.png,
    selector = "#htmlwidget_container",
    delay = delay,
    zoom = zoom,
    cliprect = "viewport",
    vwidth = width,
    vheight = height
  )
  if (is.pdf) {
    cmd <- paste("convert", tmp.png, "-density 600", file)
    system(cmd)
  } else {
    file.copy(tmp.png, file, overwrite = TRUE)
  }
  unlink(tmp.html)
}

#'
#' @export
ggLollipopPlot <- function( values, sizes=NULL, xlab="value",
                           cex.text=1) {
  df <- data.frame( x = names(values), y = values, size = 1 )
  if(!is.null(sizes)) df$size <- sizes
  ggplot2::ggplot(df, ggplot2::aes(x = reorder(x, +y), y = y)) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = reorder(x, +y),
        xend = reorder(x, +y),
        y = 0, yend = y),
      color = "gray",
      lwd = 1.3, show.legend = FALSE) +
    ggplot2::geom_point(
      mapping = ggplot2::aes(size = size),
      pch = 21, bg = 4, col = 1, show.legend = FALSE) +
    ggplot2::xlab("") +
    ggplot2::ylab(xlab) +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size=14*cex.text),
      axis.text.y = ggplot2::element_text(size=14*cex.text)
    )  
}

#' @export
plotlyLasagna <- function(df, znames=NULL) {

  zz <- sort(unique(df$z))
  min.x <- min(df$x, na.rm=TRUE)
  max.x <- max(df$x, na.rm=TRUE)
  min.y <- min(df$y, na.rm=TRUE)
  max.y <- max(df$y, na.rm=TRUE)
  
  fig <- plotly::plot_ly()
  z=1
  for(z in zz) {
    sel <- which( df$z == z )
    df1 <- df[sel,]
    mx <- c( min.x, min.x, max.x, max.x)
    my <- c( min.y, max.y, min.y, max.y)
    mz <- rep(z, 4)

    fig <- fig %>%
      plotly::add_mesh(x = mx, y = my, z = mz,
        color = "grey", opacity = 0.2)

    fig <- fig %>%
      plotly::add_trace(
        data = df1, x = ~x, y = ~y, z = z,
        mode = "markers", type = "scatter3d",
        showlegend = FALSE,
        marker = list(size = 3)
      )
    
    ztext <- z
    if(!is.null(znames) && (is.integer(z) || z %in% names(znames))) {
      ztext <- znames[z]
    }
    
    fig <- fig %>%
      plotly::add_trace(x = min.x, y = max.y, z = z, 
        type = "scatter3d", mode = "text",
        text = ztext,
        textfont = list(size = 20),
        showlegend = FALSE,
        inherit=FALSE)
    
  }
  fig
  fig <- fig %>%
    plotly::layout(
      scene = list(
        xaxis = list(title = ""),
        yaxis = list(title = ""),
        zaxis = list(title = "", showticklabels=FALSE, tickfont = list(size=0)),
        aspectmode = "cube",
        showlegend = FALSE
      )
    )
  fig
}


#' Plots time series facetted by module/colors. X is an expression
#' matrix. Time is a vector of times for each sample. 
#' 
#' @export
plotTimeSeries.modules <- function(time, xx, modules, main="",
                                   plottype="parcoord", legend=TRUE) {

  mX <- t(rowmean(scale(t(xx)), time))
  mx.time <- as.numeric(gsub("[_-].*|[^0-9]","",colnames(mX)))
  mx.time[is.na(mx.time)] <- 0
  mX <- mX[,order(mx.time)]

  dbg("[plotTimeSeries.modules] modules = ", head(modules))
  dbg("[plotTimeSeries.modules] dim.mX = ", dim(mX))  
  
  is.color <- mean(modules %in% WGCNA::standardColors(435)) > 0.8
  is.color
  if(is.numeric(modules)) {
    dbg("[plotTimeSeries.modules] modules is numeric")
    colors <- WGCNA::standardColors(435)[rank(modules)]
  } else if(is.color) {
    dbg("[plotTimeSeries.modules] modules is a color")    
    colors <- modules
  } else {
    dbg("[plotTimeSeries.modules] modules is factor")    
    colors <- WGCNA::standardColors(435)[as.integer(factor(modules))]
  }

  dbg("[plotTimeSeries.modules] colors = ", head(colors))    

  kpal <- sort(unique(as.character(colors)))
  kpal <- adjustcolor(kpal, alpha.f=0.33, 0.9, 0.9, 0.9)

  df <- data.frame( modules=modules, colors=colors,
    mX, check.names=FALSE)

  dbg("[plotTimeSeries.modules] plottype = ", plottype)
  
  if(plottype == "parcoord") {
    
    gg <- GGally::ggparcoord(
      data = df, 
      columns = c(3:ncol(df)), 
      groupColumn = 1,
      title = main)
    gg <- gg +
      ggplot2::geom_line(color = adjustcolor("grey",alpha.f=0.95)) 

  } else {

    mx.df <- reshape2::melt(df)
    mx.df$time <- as.vector(mapply(rep, mx.time, nrow(mX)))
    head(mx.df)
    gg <- ggplot2::ggplot(
      mx.df, ggplot2::aes(x=time, y=value, color=modules)) +
    ##  ggplot2::geom_line(color = adjustcolor("grey",alpha.f=0.95)) 
      ggplot2::geom_point(size=0.6)
    gg
  }

  gg <- gg + 
    ##ggplot2::geom_line(color = adjustcolor("grey",alpha.f=0.95)) +
    ## ggplot2::scale_color_manual(values=kpal) +
    ggplot2::xlab("time") + ggplot2::ylab("expression") +
    ggplot2::stat_summary(ggplot2::aes(group=colors),
      fun=mean, geom="smooth", se=FALSE, size=1.3) +
    ggplot2::facet_wrap(~modules) 
  gg  

  dbg("[plotTimeSeries.modules] legend = ", legend)  
  
  if(legend==FALSE) {
    gg <- gg + ggplot2::theme(legend.position="none")
  }
  return(gg)
}


#'
#' @export
plotTimeSeries.groups <- function(time, y, group=NULL, main="",
                                  lwd=3, xlab="time", time.factor=TRUE) {
  if(is.null(group)) group <- rep(1,length(time))
  groups <- sort(unique(group))
  ngroup <-length(groups)

  time[grepl("ctrl|co|control|bl",time,ignore.case=TRUE)] <- "0"
  ctime <- as.numeric(gsub("[-_].*|[^0-9]","",time))
  xlim <- range(ctime)
  if(time.factor) {
    ctime <- factor(ctime, levels=sort(unique(ctime)))
    ntime <- length(unique(ctime))
    xlim <- c(0,ntime)
  }
  ntime <- length(unique(ctime))
  
  i=1
  for(i in 1:ngroup) {
    j <- which(group == groups[i])
    x1 <- ctime[j]
    y1 <- y[j]
    klr <- 1+i
    klr1 <- adjustcolor(klr,alpha.f=0.5)
    klr2 <- adjustcolor(klr,red.f=0.9,green.f=0.9,blue.f=0.9) 
    x1pos <- sort(unique(x1))
    if(time.factor) x1pos <- 1:ntime
    boxplot( y1 ~ x1, col=klr1, boxwex=0.25,
      add=i>1, at = x1pos,
      xlab = xlab, ylab="expression",
      xlim=xlim, ylim=range(y))
    if(length(unique(x1))>1) {
      sy <- tapply(y1, x1, median, na.rm=TRUE)
      nk <- length(sy)*0.8
      sx <- as.numeric(names(sy))
      if(time.factor) sx <- factor(sx, levels=levels(ctime))
      lines(spline(sx,sy,n=nk),col=klr,lwd=3)
    }
    points(x1,y1,col=klr2,cex=1.8,pch=20) 
  }

  title(main=main)
  if(ngroup>1) {
    legend("bottomright",legend=groups,fill=2:99,y.intersp=0.8,cex=0.9)
  }
}


#'
#' @export
plotTimeSeries.modules_groups <- function(time, X, modules, group) {
  mx.list <- list()
  g=group[1]
  for(g in unique(group)) {
    ii <- which(group == g)
    mX <- t(rowmean(scale(t(X[,ii])), time[ii]))
    mx.list[[g]] <- data.frame(module=modules, group=g, mX, check.names=FALSE)
  }
  df <- do.call(rbind, mx.list)
  tt <- paste0(df$module,"_",df$group)
  df.avg <- rowmean(df[,3:ncol(df)], tt)
  df2 <- data.frame( do.call(rbind,strsplit(rownames(df.avg),split="_")),
                    df.avg, check.names=FALSE )
  colnames(df2)[1:2] <- c("module","group")

  mlabs <- sort(unique(df2$module))
  mlabs <- paste0(mlabs," (n=",table(df$module)[mlabs],")")
  names(mlabs) <- 1:length(mlabs)
  labeller <- ggplot2::labeller(module = mlabs)

  GGally::ggparcoord(
    data = df2,
    columns = c(3:ncol(df2)), 
    groupColumn = 2) +
    ggplot2::geom_line( ggplot2::aes(group=group), size=1) +
    ggplot2::labs(title="WGCNA modules", subtitle="Average expression in groups") + 
    ggplot2::facet_wrap(~module, labeller = labeller)

}
