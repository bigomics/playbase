##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


#' Create heatmap of top marker genes
#'
#' Identify top marker genes for groups and visualize in a heatmap.
#'
#' @param X Numeric matrix of gene expression values (genes in rows, samples in columns).
#' @param splitx Factor indicating groups for each sample.
#' @param n Number of top markers to identify per group.
#' @param ... Additional arguments passed to gx.splitmap().
#'
#' @return Heatmap created by gx.splitmap() using top markers.
#'
#' @details This function identifies the top differentially expressed genes (markers)
#' between groups defined by splitx. It calculates the mean expression of each gene in
#' each group, subtracts the overall mean, and finds the n genes with highest/lowest
#' mean expression per group. These top marker genes are extracted from X and visualized
#' in a heatmap using gx.splitmap().
#'
#' @seealso
#' \code{\link{gx.splitmap}} to see the heatmap creation function.
#'
#' @examples
#' \dontrun{
#' # TODO
#' }
#' @export
gx.markermap <- function(X, splitx, n = 5, ...) {
  F1 <- tapply(1:ncol(X), splitx, function(i) rowMeans(X[, i, drop = FALSE], na.rm = TRUE))
  F1 <- do.call(cbind, F1)
  F1 <- F1 - rowMeans(F1, na.rm = TRUE)
  topg <- apply(F1, 2, function(x) Matrix::head(order(-x), n))
  topg <- apply(topg, 2, function(i) rownames(F1)[i])
  gg <- as.vector(topg)
  gg.idx <- as.vector(sapply(colnames(topg), rep, n))
  gx.splitmap(X[gg, ], splitx = splitx, split = gg.idx, ...)
}


#' Create a heatmap of top genes for PCA components
#'
#' @param X Numeric matrix of gene expression values (genes in rows, samples in columns)
#' @param nv Number of principal components to use
#' @param ngenes Number of top genes to show for each component
#' @param ... Additional arguments passed to gx.splitmap()
#'
#' @return Heatmap created by gx.splitmap()
#'
#' @details This function calculates the top \code{ngenes} genes (by loading) for the first
#' \code{nv} principal components of \code{X}. It creates sample labels denoting the
#' component membership. The top genes are extracted and shown in a heatmap using
#' \code{gx.splitmap()}.
#'
#' @examples
#' \dontrun{
#' x <- matrix(rnorm(500 * 50), 500, 50)
#' gx.PCAheatmap(x, nv = 3, ngenes = 5)
#' }
#'
#' @export
gx.PCAheatmap <- function(X, nv = 5, ngenes = 10, ...) {
  if (inherits(X, "list") && "X" %in% names(X)) {
    X <- X$X - rowMeans(X$X, na.rm = TRUE)
  }
  X <- X - rowMeans(X, na.rm = TRUE)
  res <- irlba::irlba(X, nv = nv)
  gg <- apply(res$u, 2, function(x) Matrix::head(rownames(X)[order(-abs(x))], ngenes))
  sx <- paste0("PC.", as.vector(sapply(1:nv, rep, 10)))
  gx.splitmap(X[gg, ], split = sx, ...)
}


#' Visualize top genes for PCA components
#'
#' @param X Numeric matrix of gene expression data (genes in rows, samples in columns)
#' @param nv Number of principal components to extract
#' @param ngenes Number of top genes to visualize per component
#'
#' @details This function calculates the top \code{ngenes} genes (by loading) for the first
#' \code{nv} principal components of \code{X}. For each component, it extracts the top genes,
#' centers and scales their expression profiles, and visualizes them as a heatmap using
#' \code{gx.imagemap()}.
#'
#' @examples
#' \dontrun{
#' x <- matrix(rnorm(500 * 50), 500, 50)
#' gx.PCAcomponents(x, nv = 5, ngenes = 10)
#' }
#'
#' @export
gx.PCAcomponents <- function(X, nv = 20, ngenes) {
  if (inherits(X, "list") && "X" %in% names(X)) X <- X$X
  res <- irlba::irlba(X - rowMeans(X, na.rm = TRUE), nv = nv)
  for (i in 1:nv) {
    gg <- Matrix::head(rownames(X)[order(-abs(res$u[, i]))], ngenes)
    X1 <- X[gg, ]
    ## X1 <- (X1 - rowMeans(X1, na.rm = TRUE)) / (1e-4 + apply(X1, 1, stats::sd, na.rm = TRUE)) ## scale??
    X1 <- (X1 - rowMeans(X1, na.rm = TRUE)) / (1e-4 + matrixStats::rowSds(X1, na.rm = TRUE))
    colnames(X1) <- NULL
    gx.imagemap(X1, main = paste0("PC", i), cex = 0.8)
  }
}

#' Create an image heatmap of a data matrix
#'
#' @param X Numeric data matrix with genes/features as rows and samples/conditions as columns.
#' @param main Title for the heatmap.
#' @param cex Scaling factor for text size.
#' @param clust Logical for whether to cluster rows and columns.
#'
#' @details This function creates an image heatmap visualization of the input data matrix \code{X}.
#' If \code{clust=TRUE}, it will cluster rows and columns independently using hierarchical
#' clustering before plotting.
#'
#' It uses \code{\link[Matrix]{image}} to generate the heatmap, and adds axis labels, titles, etc.
#'
#' @examples
#' \dontrun{
#' data <- matrix(rnorm(200), 20, 10)
#' gx.imagemap(data, main = "My Heatmap")
#' }
#'
#' @export
gx.imagemap <- function(X, main = "", cex = 1, cex.main = 1.8,
                        clust = TRUE, col = NULL) {
  if (clust) {
    ii <- fastcluster::hclust(stats::dist(X))$order
    jj <- fastcluster::hclust(stats::dist(t(X)))$order
    X <- X[ii, jj]
  }
  if (is.null(col)) col <- heat.colors(64)
  graphics::image(
    1:ncol(X), 1:nrow(X), t(X),
    col = col,
    xaxt = "n", yaxt = "n",
    xlab = "", ylab = ""
  )
  if (cex > 0) {
    graphics::mtext(colnames(X),
      side = 1, at = 1:ncol(X), las = 3,
      cex = 0.8 * cex, line = 0.5
    )
    graphics::mtext(rownames(X),
      side = 4, at = 1:nrow(X), las = 1,
      cex = 0.8 * cex, line = 0.5
    )
  }
  if (cex < 0) {
    graphics::mtext("samples",
      side = 1, las = 1,
      cex = abs(cex), line = 0.65
    )
    graphics::mtext("genes",
      side = 4, las = 3,
      cex = abs(cex), line = 0.65
    )
  }
  graphics::title(main = main, cex.main = cex.main, line = 0.5)
}


#' Create a split heatmap of gene expression data using
#' ComplexHeatmap.
#'
#' Visualize a gene expression matrix as a heatmap, with samples split into groups.
#'
#' @param gx Gene expression matrix with genes in rows and samples in columns.
#' @param split Vector defining sample groups, or number of groups for clustering.
#' @param splitx External group assignments for samples.
#' @param clust.method Hierarchical clustering method for genes.
#' @param dist.method Distance metric for gene clustering.
#' @param col.dist.method Distance for sample clustering.
#' @param plot.method Heatmap plotting method (heatmap.2 or ComplexHeatmap).
#' @param scale Scaling method for the data matrix.
#' @param softmax value
#' @param order.groups value
#' @param symm.scale value
#' @param cluster_rows value
#' @param cluster_columns value
#' @param sort_columns value
#' @param row.annot Annotations for genes.
#' @param col.annot Annotations for samples.
#' @param annot.ht value
#' @param nmax Maximum number of genes to include.
#' @param cmax value
#' @param main Title for the heatmap.
#' @param verbose value
#' @param denoise value
#' @param cexRow value
#' @param cexCol value
#' @param mar value
#' @param rownames_width value
#' @param rowlab.maxlen value
#' @param title_cex value
#' @param column_title_rot value
#' @param column_names_rot value
#' @param show_legend value
#' @param show_key value
#' @param zlim value
#' @param show_rownames value
#' @param lab.len value
#' @param key.offset value
#' @param show_colnames value
#' @param use.nclust value
#' @param color_low value
#' @param color_mid value
#' @param color_high value
#'
#' @return A heatmap grob object
#'
#' @details This function subsets the input gene expression matrix to the top variable
#'  genes, clusters genes and samples, and visualizes the result as a heatmap with
#' samples split into groups defined by \code{split} or \code{splitx}.
#'
#' Many options are provided for clustering, scaling, annotations, and heatmap display
#'  parameters. By default it uses Euclidean distance and Ward clustering.
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' gx <- matrix(rnorm(250)^2, nrow = 25, ncol = 10)
#' rownames(gx) <- sample(LETTERS, 25)
#' colnames(gx) <- sample(letters, 10)
#' # Create a grouping variable
#'
#' p <- gx.splitmap(
#'   gx,
#'   scale = "row",
#'   cluster_rows = FALSE,
#'   show_key = FALSE,
#'   cluster_columns = FALSE
#' )
#' }
#' @export
gx.splitmap <- function(gx, split = 5, splitx = NULL,
                        clust.method = "ward.D2",
                        dist.method = "euclidean",
                        col.dist.method = "euclidean",
                        plot.method = "heatmap.2",
                        scale = "row",
                        softmax = 0,
                        order.groups = "clust",
                        symm.scale = FALSE,
                        cluster_rows = TRUE,
                        cluster_columns = TRUE,
                        sort_columns = NULL,
                        col.annot = NULL,
                        row.annot = NULL,
                        annot.ht = 3,
                        annot.cex = 1,
                        nmax = 1000,
                        cmax = NULL,
                        main = " ",
                        verbose = 1,
                        denoise = 0,
                        cexRow = 1,
                        cexCol = 1,
                        na_col = "green",
                        na_text = NULL,
                        mar = c(5, 5, 5, 5),
                        rowlab.maxlen = 20,
                        rownames_width = rowlab.maxlen + 10,
                        title_cex = 1.2,
                        column_title_rot = 0,
                        column_names_rot = 90,
                        show_legend = TRUE,
                        show_key = TRUE,
                        zlim = NULL,
                        show_rownames = nmax,
                        lab.len = 80,
                        key.offset = c(0.05, 0.99),
                        show_colnames = NULL,
                        use.nclust = FALSE,
                        data = FALSE,
                        color_low = "#3181de",
                        color_mid = "#eeeeee",
                        color_high = "#f23451") {
  ComplexHeatmap::ht_global_opt(fast_hclust = TRUE)
  graphics::par(xpd = FALSE)

  if (is.null(show_colnames)) {
    show_colnames <- ifelse(ncol(gx) > 100, FALSE, TRUE)
  }
  ## give unique name if duplicated
  ndup <- sum(duplicated(rownames(gx)))
  if (ndup > 0) {
    rownames(gx) <- tagDuplicates(rownames(gx))
    if (!is.null(row.annot)) rownames(row.annot) <- rownames(gx)
  }

  if (!is.null(split) && length(split) == 1 && split == 1) split <- NULL
  if (!is.null(splitx) && length(splitx) == 1 && splitx == 1) splitx <- NULL
  if (is.null(main)) main <- "  "
  if (length(mar) == 1) mar <- rep(mar[1], 4) ## old style
  if (length(mar) == 2) mar <- c(mar[1], 5, 5, mar[2]) ## old style

  cor.hclust <- function(x) {
    corx <- stats::cor(x, use = "pairwise")
    corx[is.na(corx)] <- 0
    fastcluster::hclust(stats::as.dist(1 - corx), method = "ward.D2")
  }

  ## -------------------------------------------------------------
  ## scaling options
  ## -------------------------------------------------------------
  if ("col" %in% scale || "both" %in% scale) {
    tgx <- t(gx) - colMeans(gx, na.rm = TRUE)
    gx <- t(tgx / (1e-4 + apply(gx, 2, stats::sd, na.rm = TRUE))) ## small EPS maintains SD order!
    remove(tgx)
  }
  if ("row" %in% scale || "both" %in% scale) {
    gx <- gx - rowMeans(gx, na.rm = TRUE)
    gx <- gx / (1e-4 + matrixStats::rowSds(gx, na.rm = TRUE)) ## small EPS maintains SD order!
  }
  if ("col.center" %in% scale) {
    gx <- t(t(gx) - colMeans(gx, na.rm = TRUE))
  }
  if ("row.center" %in% scale) {
    gx <- gx - rowMeans(gx, na.rm = TRUE)
  }
  if ("row.bmc" %in% scale && is.null(splitx)) {
    gx <- gx - rowMeans(gx, na.rm = TRUE)
  }
  if ("row.bmc" %in% scale && !is.null(splitx)) {
    if (inherits(splitx, "numeric") && length(splitx) == 1) {
      ## ii <- Matrix::head(order(-apply(gx, 1, stats::sd, na.rm = TRUE)), 1000) ## NEED RETHINK!
      ii <- Matrix::head(order(-matrixStats::rowSds(gx, na.rm = TRUE)), 1000)
      hc <- cor.hclust(gx[ii, ])
      splitx <- paste0("cluster", stats::cutree(hc, splitx))
      names(splitx) <- colnames(gx)
    }
    if (inherits(splitx, "character") && length(splitx) == 1 &&
      !is.null(col.annot) && splitx %in% colnames(col.annot)) {
      splitx <- as.character(col.annot[, splitx])
      names(splitx) <- colnames(gx)
    }
    for (ux in unique(splitx)) {
      jj <- which(splitx == ux)
      gx[, jj] <- gx[, jj, drop = FALSE] - rowMeans(gx[, jj, drop = FALSE], na.rm = TRUE)
    }
  }

  ## -------------------------------------------------------------
  ## Take top SD features
  ## -------------------------------------------------------------
  ##
  if (!is.null(splitx) && length(splitx) == ncol(gx)) {
    names(splitx) <- colnames(gx)
  }

  jj1 <- 1:nrow(gx)
  jj2 <- 1:ncol(gx)
  if (!is.null(nmax) && nmax < nrow(gx) && nmax > 0) {
    ## jj1 <- Matrix::head(order(-apply(gx, 1, stats::sd, na.rm = TRUE)), nmax)
    jj1 <- Matrix::head(order(-matrixStats::rowSds(gx, na.rm = TRUE)), nmax)
  }
  if (!is.null(cmax) && cmax < ncol(gx) && cmax > 0) {
    jj2 <- Matrix::head(order(-apply(gx, 2, stats::sd, na.rm = TRUE)), cmax)
  }
  gx <- gx[jj1, jj2]
  if (!is.null(row.annot)) {
    row.annot <- row.annot[jj1, , drop = FALSE]
  }
  if (!is.null(col.annot)) {
    col.annot <- col.annot[, which(colMeans(is.na(col.annot)) < 1), drop = FALSE]
    col.annot <- col.annot[jj2, , drop = FALSE]
  }

  ## --------------------------------------------
  ## split columns
  ## --------------------------------------------
  do.splitx <- !is.null(splitx)
  idx2 <- NULL
  if (do.splitx && inherits(splitx, "numeric") && length(splitx) == 1) {
    hc <- cor.hclust(gx)
    idx2 <- paste0("cluster", stats::cutree(hc, splitx))
  }
  if (do.splitx && inherits(splitx, "character") && length(splitx) == 1 &&
    !is.null(col.annot) && splitx %in% colnames(col.annot)) {
    idx2 <- as.character(col.annot[, splitx])
  }
  if (do.splitx && length(splitx) == ncol(gx)) {
    idx2 <- as.character(splitx[colnames(gx)])
  }
  if (!is.null(idx2)) {
    grp <- tapply(colnames(gx), idx2, list)
    ngrp <- length(grp)
  } else {
    grp <- list(colnames(gx))
    names(grp)[1] <- main ## title
    ngrp <- 1
  }
  remove(idx2)
  names(grp) <- paste0(names(grp), ifelse(names(grp) == "", " ", "")) ## avoid empty names

  ## --------------------------------------------
  ## split rows
  ## --------------------------------------------
  do.split <- !is.null(split)
  split.idx <- NULL
  if (do.split && inherits(split, "numeric") && length(split) == 1) {
    hc <- cor.hclust(t(gx))
    split.idx <- paste0("group", stats::cutree(hc, split))
  }
  if (do.split && inherits(split, "character") && length(split) == 1 &&
    !is.null(row.annot) && split %in% colnames(row.annot)) {
    split.idx <- row.annot[, split]
  }
  if (do.split && length(split) > 1) {
    split.idx <- split[jj1]
  }
  if (do.split && length(split) == 1 && split[1] == 1) {
    do.split <- FALSE
    split.idx <- NULL
  }

  if (!is.null(row.annot)) {
    cat(
      "[gx.splitmap] 7: duplicated annot.rownames =",
      sum(duplicated(rownames(row.annot))), "\n"
    )
  }

  if (!is.null(split.idx)) {
    row.annot <- cbind(row.annot, cluster = split.idx)
    row.annot <- data.frame(row.annot, check.names = FALSE, stringsAsFactors = FALSE)
    colnames(row.annot) <- tolower(colnames(row.annot))
    rownames(row.annot) <- rownames(gx)
  }

  ## --------------------------------------------
  ## PCA denoising
  ## --------------------------------------------
  if (denoise > 0) {
    sv <- irlba::irlba(gx, nv = denoise)
    dn <- dimnames(gx)
    gx <- sv$u %*% diag(sv$d) %*% t(sv$v)
    dimnames(gx) <- dn
  }

  ## --------------------------------------------
  ## column sorting (for what???)
  ## --------------------------------------------
  if (!cluster_columns && !is.null(sort_columns)) {
    y <- col.annot[, sort_columns, drop = FALSE]
    y1 <- as.integer(factor(y))
    jj <- fastcluster::hclust(stats::dist(t(gx)), method = "ward.D2")$order
    jj <- jj[order(1e6 * y1[jj] + 1:length(jj))]
    gx <- gx[, jj]
    col.annot <- col.annot[jj, ]
  }

  ## -------------------------------------------------------------------------------
  ## column  HeatmapAnnotation objects
  ## -------------------------------------------------------------------------------
  col.ha <- vector("list", ngrp)
  if (!is.null(col.annot)) {
    col.annot <- col.annot[, which(colMeans(is.na(col.annot)) < 1), drop = FALSE]
    col.pars <- lapply(data.frame(col.annot), function(x) sort(unique(x[!is.na(x)])))
    npar <- sapply(col.pars, length)
    par.type <- sapply(col.annot, class)
    col.colors <- list()
    i <- 1
    for (i in 1:length(npar)) {
      prm <- colnames(col.annot)[i]
      klrs <- rev(grDevices::grey.colors(npar[i], start = 0.4, end = 0.85))
      if (npar[i] == 1) klrs <- "#E6E6E6"
      names(klrs) <- col.pars[[i]]
      x <- col.annot[, i]
      is.number <- (par.type[i] %in% c("integer", "numeric"))
      if (npar[i] > 3 & !is.number) {
        klrs <- rep(RColorBrewer::brewer.pal(8, "Set2"), 99)[1:npar[i]]
        names(klrs) <- col.pars[[i]]
        klrs <- klrs[!is.na(names(klrs))]
      }
      if (any(is.na(x))) klrs <- c(klrs, "NA" = "grey90")
      if (is.number) {
        x <- as.numeric(as.character(x))
        rr <- range(x, na.rm = TRUE)
        rr <- rr + c(-1, 1) * 1e-8
        klrs <- circlize::colorRamp2(rr, c("grey85", "grey40"))
      }
      col.colors[[prm]] <- klrs
    }

    i <- 1
    for (i in 1:ngrp) {
      jj <- grp[[i]]
      ap <- list(
        title_gp = grid::gpar(fontsize = 3.3 * annot.ht * annot.cex),
        labels_gp = grid::gpar(fontsize = 2.9 * annot.ht * annot.cex),
        grid_width = grid::unit(1 * annot.ht, "mm"),
        grid_height = grid::unit(1 * annot.ht, "mm")
      )
      aa <- rep(list(ap), ncol(col.annot))
      names(aa) <- colnames(col.annot)
      col.ha[[i]] <- ComplexHeatmap::HeatmapAnnotation(
        df = col.annot[jj, , drop = FALSE],
        col = col.colors, na_col = "#FCFCFC",
        simple_anno_size = grid::unit(0.85 * annot.ht * annot.cex, "mm"), ## BioC 3.8!!
        show_annotation_name = (i == ngrp),
        show_legend = show_legend & (npar <= 20),
        annotation_name_gp = grid::gpar(fontsize = 2.9 * annot.ht * annot.cex),
        annotation_legend_param = aa
      )
    }
  }

  ## row annotation bars
  row.ha <- NULL
  if (!is.null(row.annot)) {
    npar <- apply(row.annot, 2, function(x) length(setdiff(x, c(NA, "NA"))))
    row.colors <- list()
    i <- 1
    for (i in 1:length(npar)) {
      prm <- colnames(row.annot)[i]
      x <- row.annot[, i]
      klrs <- rev(grDevices::grey.colors(npar[i], start = 0.3, end = 0.8))
      if (npar[i] == 1) klrs <- "#E6E6E6"
      if (npar[i] > 0) klrs <- omics_pal_d("muted_light")(npar)[1:npar[i]]
      names(klrs) <- sort(setdiff(unique(x), NA))
      if (any(is.na(x))) klrs <- c(klrs, "NA" = "grey90")
      row.colors[[prm]] <- klrs
    }

    row.ha <- ComplexHeatmap::rowAnnotation(
      df = row.annot,
      col = row.colors,
      show_annotation_name = show_colnames,
      show_legend = FALSE,
      annotation_name_gp = grid::gpar(fontsize = 3.3 * annot.ht * annot.cex),
      simple_anno_size = grid::unit(annot.ht, "mm"), ## BioC 3.8!!
      width = grid::unit(annot.ht * ncol(row.annot), "mm")
    )
  }

  ## -------------------------------------------------------------
  ## Plotting methods
  ## -------------------------------------------------------------

  if (softmax) {
    gx <- tanh(0.5 * gx / stats::sd(gx, na.rm = TRUE))
  }

  ## ------------- colorscale options
  col_scale <- NULL
  if (!is.null(zlim)) {
    if (length(zlim) == 3) {
      zz <- c(-zlim[1], zlim[2], zlim[3])
    } else {
      zmean <- mean(gx[, ], na.rm = TRUE)
      zz <- c(zlim[1], zmean, zlim[2])
      if (zlim[1] < 0) zz <- c(zlim[1], 0, zlim[2])
    }
    col_scale <- circlize::colorRamp2(zz, c(color_low, color_mid, color_high))
  } else if (symm.scale) {
    colmax <- 1
    colmax <- max(abs(gx[, ]), na.rm = TRUE)
    col_scale <- circlize::colorRamp2(
      c(-colmax, 0, colmax),
      c(color_low, color_mid, color_high)
    )
  } else {
    colmin <- min(gx[, ], na.rm = TRUE)
    colmax <- max(gx[, ], na.rm = TRUE)
    colmean <- mean(gx[, ], na.rm = TRUE)
    col_scale <- circlize::colorRamp2(
      c(colmin, colmean, colmax),
      c(color_low, color_mid, color_high)
    )
  }

  ## ------------- cluster blocks
  grp.order <- 1:ngrp
  if (!is.null(order.groups) && ngrp > 1 && order.groups[1] == "clust") {
    ## Reorder cluster indices based on similarity clustering
    mx <- do.call(cbind, lapply(grp, function(i) rowMeans(gx[, i, drop = FALSE], na.rm = TRUE)))
    mx <- t(scale(t(mx)))
    grp.order <- fastcluster::hclust(stats::dist(t(mx)))$order
  }
  if (!is.null(order.groups) && ngrp > 1 && length(order.groups) == ngrp) {
    grp.order <- match(order.groups, names(grp))
  }

  # Get plot data (for csv downloads)
  if (data) {
    return(gx)
  }

  ## ------------- draw heatmap

  ## global row clustering if no split
  if (cluster_rows && !do.split) {
    cluster_rows <- as.dendrogram(hclust(dist(gx)))
  }

  hmap <- NULL
  for (i in grp.order) {
    jj <- grp[[i]]

    coldistfun1 <- function(x) stats::dist(x)
    rowdistfun1 <- function(x, y) 1 - stats::cor(x, y)
    gx0 <- gx[, jj, drop = FALSE]
    grp.title <- names(grp)[i]

    gx0.colnames <- NULL
    if (show_colnames) {
      gx0.colnames <- ComplexHeatmap::HeatmapAnnotation(
        text = ComplexHeatmap::anno_text(colnames(gx0),
          rot = column_names_rot,
          gp = grid::gpar(fontsize = 11 * cexCol),
          location = grid::unit(1, "npc"), just = "right"
        ),
        annotation_height = ComplexHeatmap::max_text_width(
          colnames(gx0), grid::gpar(fontsize = 11 * cexCol)
        ) * sin((column_names_rot / 180) * pi)
      )
    }

    cell_fun <- NULL
    if (!is.null(na_text)) {
      cell_fun <- function(j, i, x, y, width, height, fill) {
        grid::grid.text(ifelse(is.na(gx0[i, j]), na_text, ""), x, y)
      }
    }

    hmap <- hmap + ComplexHeatmap::Heatmap(
      gx0,
      col = col_scale, ## from input
      cluster_rows = cluster_rows,
      cluster_columns = cluster_columns,
      clustering_distance_rows = dist.method,
      clustering_distance_columns = col.dist.method,
      clustering_method_rows = "ward.D2",
      clustering_method_columns = "ward.D2",
      show_heatmap_legend = FALSE,
      split = split.idx,
      row_title_gp = grid::gpar(fontsize = 11), ## cluster titles
      row_names_gp = grid::gpar(fontsize = 10 * cexRow),
      top_annotation = col.ha[[i]],
      column_title = grp.title,
      name = names(grp)[i],
      column_title_rot = column_title_rot,
      column_title_gp = grid::gpar(fontsize = 11 * title_cex, fontface = "plain"),
      show_row_names = FALSE,
      show_column_names = FALSE,
      bottom_annotation = gx0.colnames,
      column_names_gp = grid::gpar(fontsize = 11 * cexCol),
      na_col = na_col,
      cell_fun = cell_fun
    )
  }

  rownames.ha <- NULL
  if (FALSE && show_rownames < nrow(gx) && show_rownames > 0) {
    ## Show rownames with linked lines
    ## !!!!!!! BUGGY!!!!!!!!!!!!!
    if (1 && !is.null(split.idx) && length(unique(split.idx)) > 1) {
      nshow <- 10
      nshow <- show_rownames / length(unique(split.idx))
      nshow <- max(nshow, 10)
      subidx <- tapply(1:length(split.idx), split.idx, function(ii) {
        ## ii[utils::head(order(-apply(gx[ii, , drop = FALSE], 1, stats::sd, na.rm = TRUE)), nshow)]
        ii[utils::head(order(-matrixStats::rowSds(gx[ii, , drop = FALSE], na.rm = TRUE)), nshow)]
      })
      subset <- unlist(subidx)
    } else {
      ## subset <- Matrix::head(order(-apply(gx, 1, stats::sd, na.rm = TRUE)), show_rownames)
      subset <- Matrix::head(order(-matrixStats::rowSds(gx, na.rm = TRUE)), show_rownames)
    }
    lab <- rownames(gx)[subset]
    lab <- substring(lab, 1, lab.len)
    rownames.ha <- ComplexHeatmap::rowAnnotation(
      link = ComplexHeatmap::anno_mark(
        at = subset, labels = lab,
        link_width = grid::unit(0.8, "cm"),
        labels_gp = grid::gpar(fontsize = 10 * cexRow)
      ),
      width = grid::unit(0.8, "cm") + 0.8 * ComplexHeatmap::max_text_width(lab)
    )
  }

  if (1 && is.null(rownames.ha) && show_rownames > 0) {
    ## empty matrix just for rownames on the far right
    empty.mat <- matrix(nrow = nrow(gx), ncol = 0)
    rowlab <- substring(trimws(rownames(gx)), 1, rowlab.maxlen)
    rownames(empty.mat) <- rowlab
    if (!is.null(col.annot)) rowlab <- c(rowlab, colnames(col.annot))
    textwidth <- ComplexHeatmap::max_text_width(
      paste0(rowlab, "XXXXXXXXXX"),
      gp = grid::gpar(fontsize = 10 * cexRow)
    )
    textwidth <- grid::unit(rownames_width, "mm")
    rownames.ha <- ComplexHeatmap::Heatmap(
      empty.mat,
      row_names_max_width = textwidth,
      row_names_gp = grid::gpar(fontsize = 10 * cexRow),
      show_row_names = show_rownames > 0
    )
  }

  if (!is.null(row.ha)) {
    hmap <- hmap + row.ha
  }
  if (!is.null(rownames.ha) && show_rownames) {
    hmap <- hmap + rownames.ha
  }

  ComplexHeatmap::draw(hmap,
    annotation_legend_side = "right",
    padding = grid::unit(mar, "mm"), gap = grid::unit(1.0, "mm")
  )

  if (show_key && !is.null(col_scale)) {
    brk <- attr(col_scale, "breaks")
    lgd <- ComplexHeatmap::Legend(
      col_fun = col_scale, at = brk[c(1, 3)], labels = round(brk[c(1, 3)], 2),
      border = "transparent", labels_gp = grid::gpar(fontsize = 8),
      legend_width = grid::unit(0.08, "npc"), legend_height = grid::unit(0.01, "npc"),
      title = "\n", direction = "horizontal"
    )
    ComplexHeatmap::draw(lgd, x = grid::unit(key.offset[1], "npc"), y = grid::unit(key.offset[2], "npc"), just = c("left", "top"))
  }

  # Return TRUE so that unit test is possible
  return(TRUE)
}

#

#' Create a heatmap of a gene expression matrix
#'
#' @param gx Gene expression matrix with genes in rows and samples in columns
#' @param values value
#' @param clust.method Hierarchical clustering method for rows and columns
#' @param dist.method Distance metric for clustering genes
#' @param col.dist.method Distance metric for clustering samples
#' @param plot.method Heatmap plotting method (heatmap.2 or ComplexHeatmap)
#' @param col Colors to use for heatmap
#' @param softmax value
#' @param scale Scaling method for the data matrix
#' @param verbose value
#' @param symm value
#' @param col.annot Annotations for samples
#' @param row.annot Annotations for genes
#' @param annot.ht value
#' @param nmax Maximum number of genes to include
#' @param cmax value
#' @param show_colnames value
#' @param indent.names value
#' @param ... Additional arguments passed to heatmap plotting function
#'
#' @return A heatmap grob object
#'
#' @details This function subsets the input matrix to the top variable genes,
#' clusters genes and samples using the specified methods, and visualizes the
#' result as a heatmap using heatmap.2 or ComplexHeatmap.
#'
#' Gene and sample clustering is optional. The data can be scaled by rows or
#' columns before clustering. A color palette is used to map values to colors.
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' gx <- matrix(rnorm(250)^2, nrow = 25, ncol = 10)
#' rownames(gx) <- sample(LETTERS, 25)
#' colnames(gx) <- sample(letters, 10)
#' p <- gx.heatmap(gx)
#' }
#' @export
gx.heatmap <- function(gx, values = NULL,
                       clust.method = "ward.D2",
                       dist.method = "pearson",
                       col.dist.method = "euclidean",
                       plot.method = "heatmap.2",
                       col = grDevices::colorRampPalette(c("royalblue3", "grey90", "indianred3"))(64),
                       softmax = FALSE,
                       scale = "row", verbose = 1, symm = FALSE,
                       col.annot = NULL, row.annot = NULL, annot.ht = 1, annot.cex = 1,
                       nmax = 1000, cmax = NULL, show_colnames = TRUE,
                       indent.names = FALSE,
                       ...) {
  graphics::par(xpd = FALSE)

  if (verbose > 1) cat("input.dim.gx=", dim(gx), "\n")

  ## -------------------------------------------------------------
  ## scaling options
  ## -------------------------------------------------------------
  if (verbose > 1) cat("dim.gx=", dim(gx), "\n")
  sym0 <- FALSE
  if (sum(gx < 0, na.rm = TRUE) > 0 || scale %in% c("row", "col")) sym0 <- TRUE
  ## scaling options
  if ("col" %in% scale || "both" %in% scale) {
    tgx <- t(gx) - colMeans(gx, na.rm = TRUE)
    gx <- t(tgx / (1e-4 + apply(gx, 2, stats::sd, na.rm = TRUE))) ## small eps maintains SD order!
    remove(tgx)
  }
  if ("row" %in% scale || "both" %in% scale) {
    gx <- gx - rowMeans(gx, na.rm = TRUE)
    ## gx <- gx / (1e-4 + apply(gx, 1, stats::sd, na.rm = TRUE)) ## small eps maintains SD order!
    gx <- gx / (1e-4 + matrixStats::rowSds(gx, na.rm = TRUE))
  }
  if ("col.center" %in% scale) {
    gx <- t(t(gx) - colMeans(gx, na.rm = TRUE))
  }
  if ("row.center" %in% scale) {
    gx <- gx - rowMeans(gx, na.rm = TRUE)
  }

  ## -------------------------------------------------------------
  ## scaling options
  ## -------------------------------------------------------------
  jj1 <- 1:nrow(gx)
  jj2 <- 1:ncol(gx)
  if (!is.null(nmax) && nmax < nrow(gx) && nmax > 0) {
    ## jj1 <- Matrix::head(order(-apply(gx, 1, stats::sd, na.rm = TRUE)), nmax)
    jj1 <- Matrix::head(order(-matrixStats::rowSds(gx, na.rm = TRUE)), nmax)
  }
  if (!is.null(cmax) && cmax < ncol(gx) && cmax > 0) {
    jj2 <- Matrix::head(order(-apply(gx, 2, stats::sd, na.rm = TRUE)), cmax)
  }
  if (symm && ncol(gx) == nrow(gx)) {
    jj2 <- jj1
  }
  gx <- gx[jj1, jj2]
  if (!is.null(col.annot)) col.annot <- col.annot[jj2, , drop = FALSE]
  if (!is.null(row.annot)) row.annot <- row.annot[jj1, , drop = FALSE]

  ## ------------------------------------------------------------
  ## Distance calculation
  ## ------------------------------------------------------------
  fillNA <- function(x) {
    nx <- x
    nx[is.na(nx)] <- 0
    nx <- x + is.na(x) * rowMeans(x, na.rm = TRUE)
    nx
  }

  ## sample dimension (column)
  h1 <- NULL
  if (!is.null(clust.method)) {
    if (col.dist.method == "pearson") {
      suppressWarnings(cx <- stats::cor(gx, use = "pairwise.complete.obs"))
      cx[is.na(cx)] <- 0 ## really???
      d1 <- stats::as.dist(1 - cx)
    } else if (col.dist.method %in% c("multi.dist", "multidist")) {
      gx0 <- t(gx)
      if (any(is.na(gx0))) gx0 <- fillNA(gx0)
      d1 <- multi.dist(gx0)
    } else {
      ## euclidean
      gx0 <- t(gx)
      if (any(is.na(gx0))) gx0 <- fillNA(gx0)
      d1 <- stats::dist(gx0, method = col.dist.method)
    }
    h1 <- fastcluster::hclust(d1, method = clust.method)
  }

  ## gene dimension (rows)
  h2 <- NULL
  if (symm && ncol(gx) == nrow(gx)) {
    h2 <- h1
  } else if (!is.null(clust.method)) {
    if (dist.method == "pearson") {
      suppressWarnings(cx <- stats::cor(t(gx), use = "pairwise.complete.obs"))
      cx[is.na(cx)] <- 0
      d2 <- stats::as.dist(1 - cx)
    } else if (dist.method %in% c("multi.dist", "multidist")) {
      gx0 <- gx
      if (any(is.na(gx0))) gx0 <- fillNA(gx0)
      d2 <- multi.dist(gx0)
    } else {
      gx0 <- gx
      if (any(is.na(gx0))) gx0 <- fillNA(gx0)
      d2 <- stats::dist(gx0, method = dist.method)
    }
    h2 <- fastcluster::hclust(d2, method = clust.method)
  }
  dd <- c("both", "row", "column", "none")[1 + 1 * is.null(h1) + 2 * is.null(h2)]

  ## ------------------------------------------------------------
  ## annotation bars
  ## ------------------------------------------------------------
  ## ident names??
  if (indent.names > 0) {
    nn.sp <- sapply(floor((1:nrow(gx)) / indent.names), function(n) paste(rep(" ", n), collapse = ""))
    rownames(gx) <- paste(nn.sp, rownames(gx))
  }

  cc0 <- NULL
  cc1 <- NULL
  if (!is.null(col.annot)) {
    plot.method <- "heatmap.3"
    col.annot <- col.annot[, which(colMeans(is.na(col.annot)) < 1), drop = FALSE]
    col.annot <- as.data.frame(col.annot)
    aa <- col.annot
    is.num <- (sapply(aa, class) == "numeric")
    if (sum(is.num)) {
      for (j in which(is.num)) aa[, j] <- cut(aa[, j], breaks = 16)
    }
    ry <- sapply(aa, function(x) as.integer(as.factor(x)))
    ry <- t(ry)
    ry <- ry - apply(ry, 1, min, na.rm = TRUE)
    ry <- round((ry / (1e-8 + apply(ry, 1, max, na.rm = TRUE))) * 7)
    ry <- t(apply(ry, 1, function(x) as.integer(as.factor(x))))
    cc0 <- t(apply(ry, 1, function(x) rev(grDevices::grey.colors(max(x, na.rm = TRUE)))[x]))
    jj <- which(apply(cc0, 1, function(x) length(unique(x))) > 3 & !is.num)
    if (length(jj)) {
      klrs <- rep(RColorBrewer::brewer.pal(8, "Set2"), 99)
      klrs.mat <- matrix(klrs[ry[jj, , drop = FALSE] + 1], nrow = length(jj))
      cc0[jj, ] <- klrs.mat
    }
    cc0 <- t(cc0)
    colnames(cc0) <- colnames(col.annot)
    rownames(cc0) <- rownames(col.annot)
  }
  if (!is.null(row.annot)) {
    plot.method <- "heatmap.3"
    ry <- sapply(as.data.frame(row.annot), function(x) as.integer(as.factor(x)))
    ry <- t(ry)
    ry <- ry - apply(ry, 1, min, na.rm = TRUE)
    ry <- round((ry / (1e-8 + apply(ry, 1, max, na.rm = TRUE))) * 7)
    ry <- t(apply(ry, 1, function(x) as.integer(as.factor(x))))
    cc1 <- t(apply(ry, 1, function(x) rev(grDevices::grey.colors(max(x, na.rm = TRUE)))[x]))
    jj <- which(apply(cc1, 1, function(x) length(unique(x))) > 3)
    if (length(jj)) {
      klrs <- rep(RColorBrewer::brewer.pal(8, "Set2"), 99)
      cc1[jj, , drop = FALSE] <- matrix(klrs[ry[jj, , drop = FALSE] + 1], nrow = length(jj))
    }
    rownames(cc1) <- colnames(row.annot)
    colnames(cc1) <- rownames(row.annot)
  }

  ## -------------------------------------------------------------
  ## Plotting methods
  ## -------------------------------------------------------------
  if (!is.null(values)) gx <- values[rownames(gx), colnames(gx)]
  if (softmax) gx <- tanh(0.5 * gx / stats::sd(gx, na.rm = TRUE)) ## just for plotting...
  side.height <- 0.1 * annot.ht * NCOL(cc0)

  if (!show_colnames) {
    colnames(gx) <- rep("", ncol(gx))
  }

  if (0 && symm && !is.null(cc0)) {
    cc1 <- t(cc0)
  }
  if (!is.null(cc0)) cc0 <- cc0[colnames(gx), , drop = FALSE]
  if (!is.null(cc1)) cc1 <- cc1[, rownames(gx), drop = FALSE]

  ## draw heatmap
  if (plot.method == "heatmap.3" && !is.null(cc0) && !is.null(cc1)) {
    if (verbose > 1) cat("plotting with heatmap.3 + both ColSideColors\n")
    if (is.null(h1) && is.null(h2)) {
      heatmap.3(gx,
        Colv = NULL, Rowv = NULL,
        dendrogram = dd, col = col, scale = "none",
        symkey = sym0, symbreaks = sym0, trace = "none",
        side.height.fraction = side.height,
        ColSideColors = cc0,
        RowSideColors = cc1,
        ...
      )
    } else {
      heatmap.3(gx,
        Colv = stats::as.dendrogram(h1), Rowv = stats::as.dendrogram(h2),
        dendrogram = dd, col = col, scale = "none",
        symkey = sym0, symbreaks = sym0, trace = "none",
        side.height.fraction = side.height,
        ColSideColors = cc0,
        ## RowSideColors = matrix(cc0[,1],nrow=1) )
        RowSideColors = cc1,
        ...
      )
    }
  } else if (plot.method == "heatmap.3" && !is.null(cc0) && is.null(cc1)) {
    if (verbose > 1) cat("plotting with heatmap.3 + ColSideColors\n")
    if (is.null(h1) && is.null(h2)) {
      heatmap.3(gx,
        Colv = NULL, Rowv = NULL,
        dendrogram = dd, col = col, scale = "none",
        symkey = sym0, symbreaks = sym0, trace = "none",
        side.height.fraction = side.height,
        ColSideColors = cc0,
        ...
      )
    } else {
      heatmap.3(gx,
        Colv = stats::as.dendrogram(h1), Rowv = stats::as.dendrogram(h2),
        dendrogram = dd, col = col, scale = "none",
        symkey = sym0, symbreaks = sym0, trace = "none",
        side.height.fraction = side.height,
        ColSideColors = cc0,
        ...
      )
    }
  } else if (plot.method == "heatmap.3" && is.null(cc0) && !is.null(cc1)) {
    if (verbose > 1) cat("plotting with heatmap.3 + RowSideColors\n")
    if (is.null(h1) && is.null(h2)) {
      heatmap.3(gx,
        Colv = NULL, Rowv = NULL,
        dendrogram = dd, col = col, scale = "none",
        symkey = sym0, symbreaks = sym0, trace = "none",
        side.height.fraction = side.height,
        RowSideColors = cc1,
        ...
      )
    } else {
      heatmap.3(gx,
        Colv = stats::as.dendrogram(h1), Rowv = stats::as.dendrogram(h2),
        dendrogram = dd, col = col, scale = scale,
        symkey = sym0, symbreaks = sym0, trace = "none",
        side.height.fraction = side.height,
        RowSideColors = cc1,
        ...
      )
    }
  } else if (plot.method == "heatmap.3" && is.null(cc0) && is.null(cc1)) {
    if (verbose > 1) cat("plotting with heatmap.3 no ColSideColors\n")
    if (is.null(h1) && is.null(h2)) {
      heatmap.3(gx,
        Colv = NULL, Rowv = NULL,
        dendrogram = dd, col = col, scale = "none",
        symkey = sym0, symbreaks = sym0, trace = "none",
        side.height.fraction = side.height,
        ...
      )
    } else {
      heatmap.3(gx,
        Colv = stats::as.dendrogram(h1), Rowv = stats::as.dendrogram(h2),
        dendrogram = dd, col = col, scale = "none",
        symkey = sym0, symbreaks = sym0, trace = "none",
        side.height.fraction = side.height,
        ...
      )
    }
  } else {
    if (verbose > 1) cat("plotting with heatmap.2\n")
    if (is.null(h1) && is.null(h2)) {
      gplots::heatmap.2(gx,
        Colv = NULL, Rowv = NULL,
        dendrogram = dd, col = col, scale = "none",
        symkey = sym0, symbreaks = sym0, trace = "none", ...
      )
    } else {
      gplots::heatmap.2(gx,
        Colv = stats::as.dendrogram(h1), Rowv = stats::as.dendrogram(h2),
        dendrogram = dd, col = col, scale = "none",
        symkey = sym0, symbreaks = sym0, trace = "none", ...
      )
    }
  }

  res <- c()
  res$col.clust <- h1
  res$row.clust <- h2
  invisible(res)

  # Return res so that unit test is possible
  return(res)
}


#' Create an annotated clustered heatmap
#'
#' @param x Numeric data matrix to visualize
#' @param nc Number of sample clusters
#' @param nr Number of gene clusters
#' @param na value
#' @param q value
#' @param p value
#' @param method Distance metric for clustering (pearson, euclidean, etc)
#' @param nca value
#' @param col.annot Annotations for samples
#' @param row.annot Annotations for genes
#' @param plot value
#' @param na.fill value
#' @param nrlab value
#' @param nclab value
#' @param labrow value
#' @param labcol value
#' @param ... value
#'
#' @return A heatmap grob object if plot = FALSE; otherwise NULL
#'
#' @details This function clusters the rows and columns of a data matrix,
#' calculates cluster means, and visualizes the result as an annotated heatmap.
#'
#' It performs soft clustering, where each value is a mixture of the cluster mean
#' and original value. Side color bars are added from sample and gene annotations.
#'
#' @seealso
#' \code{\link[gplots]{heatmap.2}} for the underlying heatmap functions used
#'
#' @examples
#' \dontrun{
#' mat <- matrix(rnorm(200), 20, 10)
#' amat <- matrix(sample(letters, 200, replace = TRUE), 20, 10)
#' clustermap(mat, col.annot = amat)
#' }
#'
#' @export
clustermap <- function(x, nc = 6, nr = 6, na = 4, q = 0.80, p = 2,
                       method = "multidist",
                       nca = 24, col.annot = NULL, row.annot = NULL, plot = TRUE,
                       na.fill = FALSE, nrlab = nrow(x), nclab = ncol(x),
                       labrow = rownames(x), labcol = colnames(x),
                       ...) {
  ## non-linear transformation


  if (method == "pearson") {
    d1 <- stats::as.dist(stats::cor(t(x), use = "pairwise"))
    d2 <- stats::as.dist(stats::cor(x, use = "pairwise"))
    d1[is.na(d1)] <- mean(d1, na.rm = TRUE)
    d2[is.na(d2)] <- mean(d2, na.rm = TRUE)
    d1 <- (max(d1, na.rm = TRUE) - d1)**(p / 2)
    d2 <- (max(d2, na.rm = TRUE) - d2)**(p / 2)
    h1 <- fastcluster::hclust(d1, method = "ward.D2")
    h2 <- fastcluster::hclust(d2, method = "ward.D2")
  } else if (method == "minkowski") {
    d1 <- stats::dist(x, method = "minkowski", p = p)
    d2 <- stats::dist(t(x), method = "minkowski", p = p)
    d1[is.na(d1)] <- mean(d1, na.rm = TRUE)
    d2[is.na(d2)] <- mean(d2, na.rm = TRUE)
    h1 <- fastcluster::hclust(d1, method = "ward.D2")
    h2 <- fastcluster::hclust(d2, method = "ward.D2")
  } else if (method %in% c("multi.dist", "multidist")) {
    d1 <- multi.dist(x)
    d2 <- multi.dist(t(x))
    d1[is.na(d1)] <- mean(d1, na.rm = TRUE)
    d2[is.na(d2)] <- mean(d2, na.rm = TRUE)
    h1 <- fastcluster::hclust(d1, method = "ward.D2")
    h2 <- fastcluster::hclust(d2, method = "ward.D2")
  } else {
    stop("unknown distance method: ", method)
  }

  my.col <- grDevices::colorRampPalette(c("royalblue3", "grey90", "indianred3"))(64)
  c1 <- stats::cutree(h1, nr)
  c2 <- stats::cutree(h2, nc)
  kxmap <- function(x, c1, c2, q) {
    nr <- length(unique(c1))
    nc <- length(unique(c2))
    mx <- matrix(0, nr, nc)
    for (i in 1:nr) {
      for (j in 1:nc) {
        xx <- x[which(c1 == i), which(c2 == j)]
        mx[i, j] <- mean(xx, na.rm = TRUE)
      }
    }
    kx <- x * 0
    for (i in 1:nrow(kx)) {
      for (j in 1:ncol(kx)) {
        kx[i, j] <- mx[c1[i], c2[j]]
      }
    }
    x0 <- x
    if (na.fill == TRUE) {
      x0[is.na(x)] <- kx[is.na(x)] ## fill with cluster mean??
    }
    kx <- q * kx + (1 - q) * x0
    return(kx)
  }
  kx <- kxmap(x, c1, c2, q = q)

  ## order annotation
  order.annot <- function(A, kx, c1, n) {
    mkx <- apply(kx, 2, function(x) tapply(x, c1, mean))
    rho1 <- stats::cor(t(A), t(mkx), use = "pairwise")
    r <- apply(-abs(rho1), 2, order)
    jj <- unique(as.vector(t(r)))
    A <- A[jj, , drop = FALSE]
    A <- A[!duplicated(A), , drop = FALSE]
    Matrix::head(A, n)
  }

  ## create annotation side bars
  cc0 <- matrix("white", nrow = ncol(kx), ncol = 1)
  cc1 <- matrix("white", nrow = 1, ncol = nrow(kx))
  rownames(cc0) <- 1:ncol(kx)
  colnames(cc0) <- " "
  rownames(cc1) <- " "
  colnames(cc1) <- 1:nrow(kx)
  if (!is.null(col.annot)) {
    col.annot0 <- col.annot
    col.annot <- order.annot(col.annot, kx, c1, n = nca)
    ax <- kxmap(col.annot, 1:nrow(col.annot), c2, q = q)
    ## col.annot <- (col.annot - rowMeans(col.annot, na.rm = TRUE)) /
    ##  (1e-8 + apply(col.annot, 1, stats::sd, na.rm = TRUE))
    col.annot <- (col.annot - rowMeans(col.annot, na.rm = TRUE)) /
      (1e-8 + matrixStats::rowSds(col.annot, na.rm = TRUE))
    d3 <- stats::dist(col.annot)
    d3 <- stats::as.dist(1 - stats::cor(t(col.annot), use = "pairwise"))
    d3[is.na(d3)] <- mean(d3, na.rm = TRUE) ## impute missing...
    h3 <- fastcluster::hclust(d3, method = "ward.D2")
    na0 <- min(na, dim(col.annot))
    col.annot <- kxmap(col.annot, stats::cutree(h3, na0), c2, q = q)
    col.annot <- col.annot[h3$order, , drop = FALSE]
    ry <- t(apply(col.annot, 1, rank, na.last = "keep"))
    ry <- ry - apply(ry, 1, min, na.rm = TRUE)
    ry <- (ry / (apply(ry, 1, max, na.rm = TRUE) + 1e-8)) * 31
    cc0 <- matrix(rev(grDevices::grey.colors(32))[ry + 1], nrow = nrow(ry), ncol = ncol(ry))
    colnames(cc0) <- colnames(col.annot)
    rownames(cc0) <- rownames(col.annot)
    cc0 <- t(cc0)
  }

  if (plot == TRUE) {
    j1 <- nrow(x) - which(diff(c1[h1$order]) != 0)
    j2 <- which(diff(c2[h2$order]) != 0)
    nh <- max(ncol(cc0), nrow(cc1))**0.5
    k1 <- ((1:nrow(x) %% max(1, floor(nrow(x) / nrlab))) == 0)
    k2 <- ((1:ncol(x) %% max(1, floor(ncol(x) / nclab))) == 0)
    if (sum(!k1) > 0) {
      labrow[h1$order[!k1]] <- ""
      labrow[h1$order[k1]] <- paste(labrow[h1$order[k1]], "+ . . .")
    }
    if (sum(!k2) > 0) {
      labcol[h2$order[!k2]] <- ""
      labcol[h2$order[k2]] <- paste(". . . +", labcol[h2$order[k2]])
    }
    heatmap.3(kx,
      Colv = stats::as.dendrogram(h2), Rowv = stats::as.dendrogram(h1),
      col = my.col, labRow = labrow, labCol = labcol,
      ColSideColors = cc0, ...
    )
  }

  ## return order clustermap
  kx <- kx[h1$order, h2$order]
  invisible(kx)
}


#' Create a heatmap from downsampled data
#'
#' @param x Numeric data matrix
#' @param m Number of rows to downsample to
#' @param n Number of columns to downsample to
#' @param ... Additional arguments passed to gx.heatmap()
#'
#' @return A heatmap grob object
#'
#' @details This function downsamples the rows and columns of a large data matrix
#' to reduce its size, then visualizes the downsampled matrix as a heatmap.
#'
#' It downsamples by clustering rows and columns hierarchically, taking the
#' cluster means, and selecting representative rows/columns. The downsampled
#' matrix is plotted with gx.heatmap().
#'
#' @return
#' A heatmap grob object from gx.heatmap()
#'
#' @examples
#' \dontrun{
#' mat <- matrix(rnorm(10000), 1000, 10)
#' frozenmap(mat, m = 100, n = 5)
#' }
#'
#' @export
frozenmap <- function(x, m = 8, n = 8, ...) {
  x.downsample <- function(x, n = 100, k = 1, dist.method = "pearson") {
    if (n >= nrow(x)) {
      return(x)
    }
    if (dist.method == "multi.dist") {
      d <- multi.dist(x)
    } else if (dist.method == "pearson") {
      d <- stats::as.dist(1 - stats::cor(t(x), use = "pairwise"))
      d[is.na(d)] <- mean(d, na.rm = TRUE)
    } else {
      d <- stats::dist(x)
    }
    h <- fastcluster::hclust(d)
    idx <- stats::cutree(h, k = n)
    ordx <- -apply(h$merge, 1, min)
    ordx <- ordx[ordx > 0]
    ## sdx <- apply(x, 1, stats::sd, na.rm = TRUE)
    sdx <- matrixStats::rowSds(x, na.rm = TRUE)
    midx <- function(j) j[which.min(match(j, ordx))]
    maxsd <- function(j) j[which.max(sdx[j])]
    maxsd2 <- function(j, k) Matrix::head(j[order(-sdx[j])], k)
    ii <- sapply(1:n, function(i) which(idx == i))
    jj <- lapply(1:n, function(i) maxsd2(which(idx == i), k))
    nn <- (table(idx) > sapply(jj, length))
    hx <- t(sapply(jj, function(j) colMeans(x[j, , drop = FALSE], na.rm = TRUE)))
    hx[is.nan(hx)] <- NA
    rn <- sapply(jj, function(j) paste(rownames(x)[j], collapse = " | "))
    rownames(hx) <- paste(rn, c("", "+ ...")[1 + nn])
    hx
  }

  cx <- x.downsample(x, n = m)
  cx <- t(x.downsample(t(cx), n = n))
  res <- gx.heatmap(cx, ...)
  cx <- cx[res$row.clust$order, res$col.clust$order]
  invisible(cx)
}


#' Calculate multiple distance metrics between rows
#'
#' @param x Numeric matrix with samples in rows
#' @param p Power parameter for Minkowski distance
#' @param method Vector of distance metrics to use (pearson, euclidean, manhattan)
#'
#' @return Distance matrix
#'
#' @description This function calculates multiple distance metrics between the rows of a matrix.
#' It returns a combined distance matrix by taking the element-wise minimum of the individual distances.
#'
#' @details For each distance metric specified in \code{method}, it calculates the pairwise distance between rows of \code{x}.
#' It also calculates distances on a scaled version of \code{x} and on the ranked values.
#' The resulting distance matrices are scaled to [0,1] and combined by taking the minimum of each element across matrices.
#'
#'
#' @return
#' A combined distance matrix
#'
#' @examples
#' \dontrun{
#'
#' data <- matrix(rnorm(100), 10, 10)
#' rownames(mat) <- replicate(10, paste(sample(LETTERS, 4, replace = TRUE), collapse = ""))
#' D <- multi.dist(data)
#' }
#'
#' @export
multi.dist <- function(x, p = 4, method = c("pearson", "euclidean", "manhattan")) {
  mm <- matrix(Inf, nrow(x), nrow(x))
  sx <- t(scale(t(x)))
  rx <- t(apply(x, 1, rank, na.last = "keep"))
  xx <- list(x, sx, rx) ## parameter transform
  for (this.x in xx) {
    for (m in method) {
      if (m == "pearson") {
        d1 <- stats::as.dist(1 - stats::cor(t(this.x), use = "pairwise"))
      } else {
        d1 <- stats::dist(this.x, method = m, p = p)
      }
      m1 <- matrix(rank(as.vector(as.matrix(d1)), na.last = "keep"), nrow = nrow(x))
      m1 <- m1 / max(m1, na.rm = TRUE)
      diag(m1) <- 0
      mm <- pmin(mm, m1, na.rm = TRUE)
    }
  }
  colnames(mm) <- rownames(x)
  rownames(mm) <- rownames(x)
  stats::as.dist(mm)
}


##########################################################################
## from http://www.biostars.org/p/18211/
## https://gist.github.com/nachocab/3853004
##########################################################################

# CODE (see https://gist.github.com/nachocab/3853004)


#' Draw a heatmap with dendrograms and annotations
#'
#' @param x Numeric matrix of values to visualize
#' @param Rowv Draw row dendrogram?
#' @param Colv Draw column dendrogram?
#' @param distfun Distance function for dendrogram
#' @param hclustfun Hierarchical clustering function
#' @param dendrogram Draw dendrograms on rows, columns, both or neither
#' @param symm value
#' @param scale Scale rows, columns or neither before clustering
#' @param na.rm value
#' @param revC value
#' @param add.expr value
#' @param breaks Breakpoints for color scale
#' @param symbreaks value
#' @param col value
#' @param colsep value
#' @param rowsep value
#' @param sepcolor value
#' @param sepwidth value
#' @param cellnote value
#' @param notecex value
#' @param notecol value
#' @param na.color value
#' @param trace Draw trace lines?
#' @param tracecol value
#' @param hline value
#' @param vline value
#' @param linecol value
#' @param margins value
#' @param ColSideColors value
#' @param RowSideColors value
#' @param side.height.fraction value
#' @param cexRow value
#' @param cexCol value
#' @param labRow Labels for rows
#' @param labCol Labels for columns
#' @param key value
#' @param keysize value
#' @param density.info value
#' @param denscol value
#' @param symkey value
#' @param densadj value
#' @param main Title for plot
#' @param xlab value
#' @param ylab value
#' @param lmat value
#' @param lhei value
#' @param lwid value
#' @param NumColSideColors value
#' @param NumRowSideColors value
#' @param KeyValueName value
#' @param ... Additional arguments passed to image() and densCols()
#'
#' @return None. Produces a heatmap plot.
#'
#' @details This function creates an annotated heatmap with optional
#' row and column dendrograms. Scaling and clustering can be applied
#' to rows or columns before plotting. Many options are provided for
#' customizing the color scale, labels, titles, etc.
#'
#' @seealso \code{\link[stats]{hclust}} for clustering methods.
#'
#' @examples
#' \dontrun{
#' mat <- matrix(rnorm(200), 20, 10)
#' heatmap.3(mat, Colv = FALSE, scale = "column", dendrogram = "row")
#' }
#'
#' @export
heatmap.3 <- function(x,
                      Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
                      distfun = stats::dist,
                      hclustfun = stats::hclust,
                      dendrogram = c("both", "row", "column", "none"),
                      symm = FALSE,
                      scale = c("none", "row", "column"),
                      na.rm = TRUE,
                      revC = identical(Colv, "Rowv"),
                      add.expr,
                      breaks,
                      symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",
                      col = "heat.colors",
                      colsep,
                      rowsep,
                      sepcolor = "white",
                      sepwidth = c(0.05, 0.05),
                      cellnote,
                      notecex = 1,
                      notecol = "cyan",
                      na.color = graphics::par("bg"),
                      trace = c("none", "column", "row", "both"),
                      tracecol = "cyan",
                      hline = stats::median(breaks),
                      vline = stats::median(breaks),
                      linecol = tracecol,
                      margins = c(5, 5),
                      ColSideColors,
                      RowSideColors,
                      side.height.fraction = 0.3,
                      cexRow = 0.2 + 1 / log10(nr),
                      cexCol = 0.2 + 1 / log10(nc),
                      labRow = NULL,
                      labCol = NULL,
                      key = TRUE,
                      keysize = 1.5,
                      density.info = c("none", "histogram", "density"),
                      denscol = tracecol,
                      symkey = max(x < 0, na.rm = TRUE) || symbreaks,
                      densadj = 0.25,
                      main = NULL,
                      xlab = NULL,
                      ylab = NULL,
                      lmat = NULL,
                      lhei = NULL,
                      lwid = NULL,
                      NumColSideColors = 1,
                      NumRowSideColors = 1,
                      KeyValueName = "Value", ...) {
  invalid <- function(x) {
    if (missing(x) || is.null(x) || length(x) == 0) {
      return(TRUE)
    }
    if (is.list(x)) {
      return(all(sapply(x, invalid)))
    } else if (is.vector(x)) {
      return(all(is.na(x)))
    } else {
      return(FALSE)
    }
  }

  x <- as.matrix(x)
  scale01 <- function(x, low = min(x, na.rm = TRUE), high = max(x, na.rm = TRUE)) {
    x <- (x - low) / (high - low)
    x
  }
  retval <- list()
  scale <- if (symm && missing(scale)) {
    "none"
  } else {
    match.arg(scale)
  }
  dendrogram <- match.arg(dendrogram)
  trace <- match.arg(trace)
  density.info <- match.arg(density.info)
  if (length(col) == 1 && is.character(col)) {
    col <- get(col, mode = "function")
  }
  if (!missing(breaks) && (scale != "none")) {
    warning(
      "Using scale=\"row\" or scale=\"column\" when breaks are",
      "specified can produce unpredictable results.", "Please consider using only one or the other."
    )
  }
  if (is.null(Rowv) || all(is.na(Rowv))) {
    Rowv <- FALSE
  }
  if (is.null(Colv) || all(is.na(Colv))) {
    Colv <- FALSE
  } else if (Colv[1] == "Rowv" && !all(isTRUE(Rowv))) {
    Colv <- FALSE
  }
  if (length(di <- dim(x)) != 2 || !is.numeric(x)) {
    stop("`x' must be a numeric matrix")
  }
  nr <- di[1]
  nc <- di[2]
  if (nr <= 1 || nc <= 1) {
    stop("`x' must have at least 2 rows and 2 columns")
  }
  if (!is.numeric(margins) || length(margins) != 2) {
    stop("`margins' must be a numeric vector of length 2")
  }
  if (missing(cellnote)) {
    cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
  }
  if (!inherits(Rowv, "dendrogram")) {
    if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
      c("both", "row"))) {
      if (is.logical(Colv) && (Colv)) {
        dendrogram <- "column"
      } else {
        dendrogram <- "none"
      }
      warning(
        "Discrepancy: Rowv is FALSE, while dendrogram is `",
        dendrogram, "'. Omitting row dendogram."
      )
    }
  }
  if (!inherits(Colv, "dendrogram")) {
    if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
      c("both", "column"))) {
      if (is.logical(Rowv) && (Rowv)) {
        dendrogram <- "row"
      } else {
        dendrogram <- "none"
      }
      warning(
        "Discrepancy: Colv is FALSE, while dendrogram is `",
        dendrogram, "'. Omitting column dendogram."
      )
    }
  }
  if (inherits(Rowv, "dendrogram")) {
    ddr <- Rowv
    rowInd <- stats::order.dendrogram(ddr)
  } else if (is.integer(Rowv)) {
    hcr <- hclustfun(distfun(x))
    ddr <- stats::as.dendrogram(hcr)
    ddr <- stats::reorder(ddr, Rowv)
    rowInd <- stats::order.dendrogram(ddr)
    if (nr != length(rowInd)) {
      stop("row dendrogram ordering gave index of wrong length")
    }
  } else if (isTRUE(Rowv)) {
    Rowv <- rowMeans(x, na.rm = na.rm)
    hcr <- hclustfun(distfun(x))
    ddr <- stats::as.dendrogram(hcr)
    ddr <- stats::reorder(ddr, Rowv)
    rowInd <- stats::order.dendrogram(ddr)
    if (nr != length(rowInd)) {
      stop("row dendrogram ordering gave index of wrong length")
    }
  } else {
    rowInd <- nr:1
  }
  if (inherits(Colv, "dendrogram")) {
    ddc <- Colv
    colInd <- stats::order.dendrogram(ddc)
  } else if (identical(Colv, "Rowv")) {
    if (nr != nc) {
      stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
    }
    if (exists("ddr")) {
      ddc <- ddr
      colInd <- stats::order.dendrogram(ddc)
    } else {
      colInd <- rowInd
    }
  } else if (is.integer(Colv)) {
    hcc <- hclustfun(distfun(if (symm) {
      x
    } else {
      t(x)
    }))
    ddc <- stats::as.dendrogram(hcc)
    ddc <- stats::reorder(ddc, Colv)
    colInd <- stats::order.dendrogram(ddc)
    if (nc != length(colInd)) {
      stop("column dendrogram ordering gave index of wrong length")
    }
  } else if (isTRUE(Colv)) {
    Colv <- colMeans(x, na.rm = na.rm)
    hcc <- hclustfun(distfun(if (symm) {
      x
    } else {
      t(x)
    }))
    ddc <- stats::as.dendrogram(hcc)
    ddc <- stats::reorder(ddc, Colv)
    colInd <- stats::order.dendrogram(ddc)
    if (nc != length(colInd)) {
      stop("column dendrogram ordering gave index of wrong length")
    }
  } else {
    colInd <- 1:nc
  }
  retval$rowInd <- rowInd
  retval$colInd <- colInd
  retval$call <- match.call()
  x <- x[rowInd, colInd]
  x.unscaled <- x
  cellnote <- cellnote[rowInd, colInd]
  if (is.null(labRow)) {
    labRow <- if (is.null(rownames(x))) {
      (1:nr)[rowInd]
    } else {
      rownames(x)
    }
  } else {
    labRow <- labRow[rowInd]
  }
  if (is.null(labCol)) {
    labCol <- if (is.null(colnames(x))) {
      (1:nc)[colInd]
    } else {
      colnames(x)
    }
  } else {
    labCol <- labCol[colInd]
  }
  if (scale == "row") {
    retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
    x <- sweep(x, 1, rm)
    ## retval$rowSDs <- sx <- apply(x, 1, stats::sd, na.rm = na.rm)
    retval$rowSDs <- sx <- matrixStats::rowSds(x, na.rm = na.rm)
    x <- sweep(x, 1, sx, "/")
  } else if (scale == "column") {
    retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
    x <- sweep(x, 2, rm)
    retval$colSDs <- sx <- apply(x, 2, stats::sd, na.rm = na.rm)
    x <- sweep(x, 2, sx, "/")
  }
  if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
    if (missing(col) || is.function(col)) {
      breaks <- 16
    } else {
      breaks <- length(col) + 1
    }
  }
  if (length(breaks) == 1) {
    if (!symbreaks) {
      breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
        length = breaks
      )
    } else {
      extreme <- max(abs(x), na.rm = TRUE)
      breaks <- seq(-extreme, extreme, length = breaks)
    }
  }
  nbr <- length(breaks)
  ncol <- length(breaks) - 1
  if (inherits(col, "function")) {
    col <- col(ncol)
  }
  min.breaks <- min(breaks, na.rm = TRUE)
  max.breaks <- max(breaks, na.rm = TRUE)
  x[x < min.breaks] <- min.breaks
  x[x > max.breaks] <- max.breaks
  if (missing(lhei) || is.null(lhei)) {
    lhei <- c(keysize, 4)
  }
  if (missing(lwid) || is.null(lwid)) {
    lwid <- c(keysize, 4)
  }
  if (missing(lmat) || is.null(lmat)) {
    lmat <- rbind(4:3, 2:1)

    if (!missing(ColSideColors)) {
      if (!is.character(ColSideColors) || nrow(ColSideColors) != nc) {
        stop("'ColSideColors' must be a matrix of nrow(x) rows")
      }
      lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)

      lhei <- c(lhei[1], side.height.fraction * NumColSideColors, lhei[2])
    }

    if (!missing(RowSideColors)) {
      if (!is.character(RowSideColors) || ncol(RowSideColors) != nr) {
        stop("'RowSideColors' must be a matrix of ncol(x) columns")
      }
      lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[, 2] + 1)

      lwid <- c(lwid[1], side.height.fraction * NumRowSideColors, lwid[2])
    }
    lmat[is.na(lmat)] <- 0
  }

  if (length(lhei) != nrow(lmat)) {
    stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
  }
  if (length(lwid) != ncol(lmat)) {
    stop("lwid must have length = ncol(lmat) =", ncol(lmat))
  }
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op))

  plotly::layout(lmat, widths = lwid, heights = lhei, respect = FALSE)

  if (!missing(RowSideColors)) {
    if (!is.matrix(RowSideColors)) {
      graphics::par(mar = c(margins[1], 0, 0, 0.5))
      Matrix::image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
    } else {
      graphics::par(mar = c(margins[1], 0, 0, 0.5))
      rsc <- t(RowSideColors[, rowInd, drop = F])
      rsc.colors <- matrix()
      rsc.names <- names(table(rsc))
      rsc.i <- 1
      for (rsc.name in rsc.names) {
        rsc.colors[rsc.i] <- rsc.name
        rsc[rsc == rsc.name] <- rsc.i
        rsc.i <- rsc.i + 1
      }
      rsc <- matrix(as.numeric(rsc), nrow = dim(rsc)[1])
      Matrix::image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
      if (length(colnames(RowSideColors)) > 0) {
        graphics::axis(1, 0:(dim(rsc)[2] - 1) / (dim(rsc)[2] - 1), colnames(RowSideColors), las = 2, tick = FALSE)
      }
    }
  }

  if (!missing(ColSideColors)) {
    if (!is.matrix(ColSideColors)) {
      graphics::par(mar = c(0.5, 0, 0, margins[2]))
      Matrix::image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    } else {
      graphics::par(mar = c(0.5, 0, 0, margins[2]))
      csc <- ColSideColors[colInd, , drop = F]
      csc.colors <- matrix()
      csc.names <- names(table(csc))
      csc.i <- 1
      for (csc.name in csc.names) {
        csc.colors[csc.i] <- csc.name
        csc[csc == csc.name] <- csc.i
        csc.i <- csc.i + 1
      }
      csc <- matrix(as.numeric(csc), nrow = dim(csc)[1])
      Matrix::image(csc, col = as.vector(csc.colors), axes = FALSE)
      if (length(colnames(ColSideColors)) > 0) {
        graphics::axis(2, 0:(dim(csc)[2] - 1) / max(1, (dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
      }
    }
  }

  graphics::par(mar = c(margins[1], 0, 0, margins[2]))
  x <- t(x)
  cellnote <- t(cellnote)
  if (revC) {
    iy <- nr:1
    if (exists("ddr")) {
      ddr <- rev(ddr)
    }
    x <- x[, iy]
    cellnote <- cellnote[, iy]
  } else {
    iy <- 1:nr
  }
  Matrix::image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
  retval$carpet <- x
  if (exists("ddr")) {
    retval$rowDendrogram <- ddr
  }
  if (exists("ddc")) {
    retval$colDendrogram <- ddc
  }
  retval$breaks <- breaks
  retval$col <- col
  if (!invalid(na.color) & any(is.na(x))) { # load
    mmat <- ifelse(is.na(x), 1, NA)
    Matrix::image(1:nc, 1:nr, mmat,
      axes = FALSE, xlab = "", ylab = "",
      col = na.color, add = TRUE
    )
  }
  graphics::axis(1, 1:nc,
    labels = labCol, las = 2, line = -0.2, tick = 0,
    cex.axis = cexCol
  )
  if (!is.null(xlab)) {
    graphics::mtext(xlab, side = 1, line = margins[1] - 1.25)
  }
  graphics::axis(4, iy,
    labels = labRow, las = 2, line = -0.2, tick = 0,
    cex.axis = cexRow
  )
  if (!is.null(ylab)) {
    graphics::mtext(ylab, side = 4, line = margins[2] - 1.25)
  }
  if (!missing(add.expr)) {
    eval(substitute(add.expr))
  }
  if (!missing(colsep)) {
    for (csep in colsep) graphics::rect(xleft = csep + 0.5, ybottom = rep(0, length(csep)), xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  }
  if (!missing(rowsep)) {
    for (rsep in rowsep) graphics::rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  }
  min.scale <- min(breaks, na.rm = TRUE)
  max.scale <- max(breaks, na.rm = TRUE)
  x.scaled <- scale01(t(x), min.scale, max.scale)
  if (trace %in% c("both", "column")) {
    retval$vline <- vline
    vline.vals <- scale01(vline, min.scale, max.scale)
    for (i in colInd) {
      if (!is.null(vline)) {
        graphics::abline(
          v = i - 0.5 + vline.vals, col = linecol,
          lty = 2
        )
      }
      xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
      xv <- c(xv[1], xv)
      yv <- 1:length(xv) - 0.5
      graphics::lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (trace %in% c("both", "row")) {
    retval$hline <- hline
    hline.vals <- scale01(hline, min.scale, max.scale)
    for (i in rowInd) {
      if (!is.null(hline)) {
        graphics::abline(h = i + hline, col = linecol, lty = 2)
      }
      yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
      yv <- rev(c(yv[1], yv))
      xv <- length(yv):1 - 0.5
      graphics::lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (!missing(cellnote)) {
    graphics::text(
      x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
      col = notecol, cex = notecex
    )
  }
  graphics::par(mar = c(margins[1], 0, 0, 0))
  if (dendrogram %in% c("both", "row")) {
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  } else {
    graphics::plot.new()
  }
  graphics::par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
  if (dendrogram %in% c("both", "column")) {
    plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
  } else {
    graphics::plot.new()
  }
  if (!is.null(main)) {
    graphics::title(main, cex.main = 1.5 * op[["cex.main"]])
  }
  if (key) {
    graphics::par(mar = c(5, 4, 2, 1), cex = 0.75)
    tmpbreaks <- breaks
    if (symkey) {
      max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
      min.raw <- -max.raw
      tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
      tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
    } else {
      min.raw <- min(x, na.rm = TRUE)
      max.raw <- max(x, na.rm = TRUE)
    }

    z <- seq(min.raw, max.raw, length = length(col))
    Matrix::image(
      z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
      xaxt = "n", yaxt = "n"
    )
    graphics::par(usr = c(0, 1, 0, 1))
    lv <- pretty(breaks)
    xv <- scale01(as.numeric(lv), min.raw, max.raw)
    graphics::axis(1, at = xv, labels = lv)
    if (scale == "row") {
      graphics::mtext(side = 1, "Row Z-Score", line = 2)
    } else if (scale == "column") {
      graphics::mtext(side = 1, "Column Z-Score", line = 2)
    } else {
      graphics::mtext(side = 1, KeyValueName, line = 2)
    }
    if (density.info == "density") {
      dens <- stats::density(x, adjust = densadj, na.rm = TRUE)
      omit <- dens$x < min(breaks, na.rm = TRUE) | dens$x > max(breaks, na.rm = TRUE)
      dens$x <- dens$x[-omit]
      dens$y <- dens$y[-omit]
      dens$x <- scale01(dens$x, min.raw, max.raw)
      graphics::lines(dens$x, dens$y / max(dens$y, na.rm = TRUE) * 0.95,
        col = denscol,
        lwd = 1
      )
      graphics::axis(2, at = pretty(dens$y) / max(dens$y, na.rm = TRUE) * 0.95, pretty(dens$y))
      graphics::title("Color Key\nand Density Plot")
      graphics::par(cex = 0.5)
      graphics::mtext(side = 2, "Density", line = 2)
    } else if (density.info == "histogram") {
      h <- graphics::hist(x, plot = FALSE, breaks = breaks)
      hx <- scale01(breaks, min.raw, max.raw)
      hy <- c(h$counts, h$counts[length(h$counts)])
      graphics::lines(hx, hy / max(hy, na.rm = TRUE) * 0.95,
        lwd = 1, type = "s",
        col = denscol
      )
      graphics::axis(2, at = pretty(hy) / max(hy, na.rm = TRUE) * 0.95, pretty(hy))
      graphics::title("Color Key\nand Histogram")
      graphics::par(cex = 0.5)
      graphics::mtext(side = 2, "Count", line = 2)
    } else {
      graphics::title("Color Key")
    }
  } else {
    graphics::plot.new()
  }
  retval$colorTable <- data.frame(
    low = retval$breaks[-length(retval$breaks)],
    high = retval$breaks[-1], color = retval$col
  )
  invisible(retval)
}
