##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


#' @title Create data for hive plot visualization
#'
#' @description
#' Generates data frames required to visualize associations between omics datasets
#' using a hive plot.
#'
#' @param res Data frame of correlation statistics between features and phenotype.
#' @param rho.min Minimum correlation threshold.
#' @param cxi Node circle size parameter.
#' @param use.alpha Use transparency based on correlation strength.
#' @param ntop Number of top features to include.
#' @param ew Edge width parameter.
#' @param nv Number of vertices for node polygons.
#' @param dx.dir X axis direction.
#' @param rx Hive plot x-radius.
#'
#' @details
#' This function takes a data frame of correlation statistics between features in
#' gene expression, copy number, methylation datasets vs a phenotype. It filters,
#' orders, and processes the data into formats needed for generating a hive plot using the
#' function.
#'
#' The main outputs are two nested data frames defining the nodes and edges for the plot.
#' Visual parameters like node size, color, transparency, etc. are set based on the input
#' correlation data.
#'
#' @return
#' A list with two data frames 'nodes' and 'edges' containing node and edge properties
#' for generating a hive plot.
#' @export
omx.makeHivePlotData_ <- function(res, rho.min = 0.15, cxi = 0.11, use.alpha = TRUE,
                                  ntop = 1000, ew = 4, nv = 10, dx.dir = c(-1, 1), rx = 2.2) {
  ## omics score
  res <- res[order(-abs(res[, "score"])), ]
  if (ntop > 0) res <- Matrix::head(res, ntop)

  ## stuff...
  dx <- res[, "gx.rho"]
  dc <- res[, "cn.rho"]
  dm <- res[, "me.rho"]
  n <- nrow(res)
  names(dx) <- paste("x", rownames(res), sep = "")
  names(dc) <- paste("c", rownames(res), sep = "")
  names(dm) <- paste("m", rownames(res), sep = "")
  dx[is.na(dx)] <- 0
  dc[is.na(dc)] <- 0
  dm[is.na(dm)] <- 0

  ## define node properties
  dx0 <- dc0 <- dm0 <- rho.min ## threshold for coloring
  dx0 <- dc0 <- dm0 <- 0 ## threshold for coloring
  kr0 <- c("grey40", "green4", "red")[1 + 1 * (dx < -dx0) + 2 * (dx > dx0)]
  n0 <- data.frame(
    id = 1:n, lab = names(dx), axis = 1,
    radius = 2 * (rank(dx, na.last = "keep") / n - 0.5), value = dx,
    size = abs(dx) / max(abs(dx), na.rm = TRUE), color = kr0
  )
  kr1 <- c("grey40", "blue", "orange2")[1 + 1 * (dc < -dc0) + 2 * (dc > dc0)]
  n1 <- data.frame(
    id = (n + 1):(2 * n), lab = names(dc), axis = 2,
    radius = 2 * (rank(dc, na.last = "keep") / n - 0.5), value = dc,
    size = abs(dc) / max(abs(dc), na.rm = TRUE), color = kr1
  )
  kr2 <- c("grey40", "purple", "gold")[1 + 1 * (dm < -dm0) + 2 * (dm > dm0)]
  n2 <- data.frame(
    id = (2 * n + 1):(3 * n), lab = names(dm), axis = 3,
    radius = 2 * (rank(dm, na.last = "keep") / n - 0.5), value = dm,
    size = abs(dm) / max(abs(dm), na.rm = TRUE), color = kr2
  )
  nn <- rbind(n0, n1, n2)
  nn$size <- (0.10 + 0.60 * nn$size * 2) * 2
  #
  nn$lab <- as.character(nn$lab)
  nn$axis <- as.integer(nn$axis)
  nn$color <- as.character(nn$color)

  ## define edge properties
  ee1 <- data.frame(n1 = names(dx), n2 = names(dc), cor = res[, "cn.tau"])
  ee2 <- data.frame(n1 = names(dx), n2 = names(dm), cor = res[, "me.tau"])
  #
  ee <- rbind(ee1, ee2)
  ee <- ee[which(abs(ee$cor) > rho.min), ]
  ee$id1 <- nn$id[match(ee$n1, nn$lab)]
  ee$id2 <- nn$id[match(ee$n2, nn$lab)]
  ee$weight <- (abs(ee$cor) - rho.min) / (1 - rho.min)
  #
  if (use.alpha) {
    alpha <- abs(ee$weight * nn[ee$id1, "radius"] * nn[ee$id2, "radius"])
    alpha <- alpha - min(alpha)
    #
    alpha <- pmax((alpha / max(alpha))**1.0, 0.20)
  } else {
    alpha <- 1
  }
  ee$weight <- 0.0 + ew * (ee$weight)**1.5 ## add minimum and gamma
  ee$alpha <- alpha
  #
  ee$color <- c("steelblue", "darkorange")[1 + 1 * (ee$cor > 0)]
  cc2 <- cbind(t(grDevices::col2rgb(ee$color) / 255), alpha)
  ee$color <- apply(cc2[, ], 1, function(x) grDevices::rgb(x[1], x[2], x[3], x[4]))
  ee <- ee[!duplicated(ee), ]
  ee <- ee[which(as.character(ee$n1) != as.character(ee$n2)), ]

  ## node annotation matrix
  rr0 <- res[which(dx > 0 & abs(res[, "score"]) > 0), ] ## upregulated
  rr1 <- res[which(dx < 0 & abs(res[, "score"]) > 0), ] ## downregulated
  rr0 <- Matrix::head(rr0, nv)
  rr1 <- Matrix::head(rr1, nv)
  rr0 <- rr0[order(-rr0[, "gx.rho"]), ]
  rr1 <- rr1[order(+rr1[, "gx.rho"]), ]
  yy0 <- seq(1, 1 - nv * cxi, -cxi)[1:nrow(rr0)]
  yy1 <- seq(-1, -1 + nv * cxi, cxi)[1:nrow(rr1)] + 0.5
  rownames(rr0) <- paste("x", rownames(rr0), sep = "")
  rownames(rr1) <- paste("x", rownames(rr1), sep = "")
  yy0 <- yy0 - nn[rownames(rr0), ]$radius + 0.2
  yy1 <- yy1 - nn[rownames(rr1), ]$radius - 0.2 ## offset from top
  M0 <- data.frame(
    node.lab = rownames(rr0),
    node.text = sub("^x", "", rownames(rr0)),
    angle = atan(yy0 / rx) / pi * 180, radius = sqrt(rx**2 + yy0**2),
    offset = 0, hjust = 0, vjust = 0
  )
  M1 <- data.frame(
    node.lab = rownames(rr1),
    node.text = sub("^x", "", rownames(rr1)),
    angle = 180 - atan(yy1 / rx) / pi * 180, radius = sqrt(rx**2 + yy1**2),
    offset = 0, hjust = 1, vjust = 0
  )
  mm <- c()
  if (1 %in% dx.dir) mm <- rbind(mm, M0)
  if (-1 %in% dx.dir) mm <- rbind(mm, M1)
  rownames(res) <- sub("^x", "", rownames(res))

  ## HPD object
  hpd <- c()
  hpd$nodes <- nn
  hpd$edges <- ee
  hpd$type <- "2D"
  hpd$desc <- paste("3 axes --", nrow(nn), "nodes --", nrow(ee), "edges")
  #
  hpd$axis.cols <- rep("grey40", 3)
  hpd$score <- res
  hpd$nodes.ann <- mm
  class(hpd) <- "HivePlotData"
  return(hpd)
}


#' #' @title Create a hive plot visualization for an NGS object
#' #'
#' #' @description
#' #' Creates a hive plot to visualize associations between omics datasets
#' #' in an NGS object.
#' #'
#' #' @param ngs An NGS object containing multiple omics datasets.
#' #' @param pheno The phenotype data column name.
#' #' @param level The omics data level (1, 2, 3).
#' #' @param ntop Number of top features to include.
#' #' @param main Plot title.
#' #' @param axis.lab Axis labels c("GX", "CN", "ME").
#' #' @param bkgnd Background color.
#' #' @param rx Hive plot x-radius.
#' #' @param cex Text size factor.
#' #' @param slen Feature name substring length.
#' #'
#' #' @details
#' #' This function takes an NGS object containing gene expression (GX), copy number (CN),
#' #' and methylation (ME) data. It calculates correlation statistics between features in
#' #' each dataset vs. the phenotype. It selects the top \code{ntop} features and creates
#' #' a hive plot visualization, with hubs representing GX/CN/ME and edges representing
#' #' correlations between features across datasets related to the phenotype.
#' #'
#' #' The plot is customizable via parameters like \code{main}, \code{bkgnd}, \code{cex}, etc.
#' #'
#' #' @return
#' #' A hive plot grob object is returned, no value.
#' #'
#' #' @export
#' ngs.hiveplot <- function(ngs, pheno, level = 1, ntop = 400, main = "", axis.lab = c("GX", "CN", "ME"),
#'                          bkgnd = "white", rx = 2.2, cex = 1, slen = -1) {
#'   res <- omx$zstats[[level]]
#'   res <- cbind(res, omx$stats[[pheno]][[level]])
#'   res <- res[!grepl("^PHENO", rownames(res)), ]
#'   rownames(res) <- gsub("^.*:", "", rownames(res)) ## strip prefix
#'   if (slen > 0) rownames(res) <- substr(rownames(res), 1, slen)
#'   ann.cex <- cex * 0.7
#'
#'   hpd <- omx.makeHivePlotData_(res,
#'     rho.min = 0.15, ntop = ntop, rx = rx,
#'     cxi = (0.05 + ann.cex / 7)
#'   )
#'   write.csv(hpd$nodes.ann, file = "/tmp/annode.csv", row.names = FALSE)
#'   grid::grid.newpage()
#'   axlab.col <- ifelse(bkgnd == "black", "grey90", "grey15")
#'   HiveR::plotHive(hpd,
#'     np = FALSE, ch = 1.4, bkgnd = bkgnd,
#'     axLabs = c(paste0(main, "\n", axis.lab[1]), axis.lab[2], axis.lab[3]),
#'     axLab.pos = c(0.4, 0.33, 0.33),
#'     axLab.gpar = grid::gpar(col = axlab.col, cex = cex * 1.3),
#'     anNodes = "/tmp/annode.csv",
#'     anNode.gpar = grid::gpar(cex = 0.7 * cex, col = axlab.col, lwd = 0.50)
#'   )
#' }


#' @title Create Hive Plot for mixed network
#'
#' @param res Input data
#' @param ngs NGS object containing expression data
#' @param ct Contrast name or column in \code{ngs$samples} to use for node coloring
#' @param showloops Whether to show looping edges
#' @param numlab Number of node labels to show
#' @param cex Label size scaling factor
#'
#' @return Hive plot object
#'
#' @description Creates a Hive Plot to visualize a mixed gene-gene set network.
#'
#' @details  It extracts the network graph and formats it into a Hive Plot layout.
#'
#' Node importance scores are extracted from the NGS object based on the specified \code{ct} contrast.
#' These scores are used to determine node sizes in the plot.
#'
#' The number of node labels can be reduced by setting \code{numlab} to avoid overplotting.
#'
#' @export
mixHivePlot <- function(res, ngs, ct, showloops = FALSE, numlab = 6, cex = 1) {
  cat("<mixHivePlot> called\n")
  if (is.null(showloops)) showloops <- FALSE

  gr <- res$graph

  ## -------------------------------------------------------------
  ## Prepare the HivePlot data
  ## -------------------------------------------------------------
  df <- data.frame(res$edges)
  hpd <- edge2HPD(df, axis.cols = rep("grey", 3))

  hpd$edges$looping <- res$edges$looping
  hpd$edges$weight <- res$edges$importance
  loop.nodes <- unique(
    c(
      as.character(res$edges$from)[res$edges$looping],
      as.character(res$edges$to)[res$edges$looping]
    )
  )
  hpd$nodes$looping <- hpd$nodes$lab %in% loop.nodes
  hpd <- mineHPD(hpd, option = "rad <- tot.edge.count")
  hpd$nodes$degree <- hpd$nodes$radius
  hpd$edges$from <- hpd$nodes$lab[hpd$edges$id1]
  hpd$edges$to <- hpd$nodes$lab[hpd$edges$id2]

  if (is.null(ct) || is.null(ngs)) {
    fx <- tapply(igraph::V(gr)$importance, sub(".*:", "", igraph::V(gr)$name), max)
  } else if (ct %in% names(ngs$gx.meta$meta)) {
    ## use fold change as radial layout
    fc <- ngs$gx.meta$meta[[ct]]$meta.fx
    names(fc) <- rownames(ngs$gx.meta$meta[[ct]])
    gs <- ngs$gset.meta$meta[[ct]]$meta.fx
    names(gs) <- rownames(ngs$gset.meta$meta[[ct]])
    fc <- fc / max(abs(fc), na.rm = TRUE)
    gs <- gs / max(abs(gs), na.rm = TRUE)
    fx <- c(fc, gs)
  } else if (ct %in% colnames(ngs$samples)) {
    group <- ngs$samples[, ct]

    design <- stats::model.matrix(~ 0 + group)
    colnames(design) <- sub("group", "", colnames(design))
    fit <- limma::eBayes(limma::lmFit(ngs$X, design))
    stat <- limma::topTable(fit, number = Inf)
    fx <- stat$F
    names(fx) <- rownames(ngs$X)
  } else {
    stop("FATAL:: mixHivePlot: unknown contrast/conditions=", ct, "\n")
  }

  fx <- fx / max(abs(fx), na.rm = TRUE)
  g <- sub("[1-9]:", "", hpd$nodes$lab)
  hpd$nodes$radius <- rank(fx[g], na.last = "keep")
  hpd$nodes$radius <- 100 * hpd$nodes$radius / max(hpd$nodes$radius, na.rm = TRUE)

  maxgrp <- unlist(lapply(res$W, function(w) max.col(w)))

  names(maxgrp) <- as.vector(sapply(res$W, rownames))

  ## use importance as node size
  importance <- igraph::V(gr)$importance
  names(importance) <- igraph::V(gr)$name

  hpd$nodes$size <- abs(importance[hpd$nodes$lab])
  hpd$nodes$size <- 1.6 * (hpd$nodes$size / max(hpd$nodes$size))**0.5
  hpd$nodes$axis <- as.integer(sub(":.*", "", hpd$nodes$lab))
  hpd$nodes$color <- c("red3", "blue2")[1 + 1 * (fx[g] > 0)]
  wt1 <- hpd$edges$weight ## edge.importance
  wt1 <- rank(abs(wt1), na.last = "keep") * sign(wt1)
  hpd$edges$weight <- 3 * abs(wt1 / max(abs(wt1)))**2
  hpd$edges$color <- psych::alpha("grey70", 0.3)

  jj <- which(hpd$edges$looping)
  if (showloops && length(jj)) {
    hpd$edges[jj, ]
    hpd$edges$color <- psych::alpha("grey70", 0.2)
    hpd$edges$color[jj] <- psych::alpha("red3", 0.3)
  }

  axis.names <- names(res$X)
  makeAcronym <- function(x) {
    x <- gsub("[)(]", "", x)
    sapply(strsplit(x, split = "[_ -]"), function(s) {
      if (length(s) == 1) {
        return(substring(s, 1, 2))
      }
      toupper(paste(substring(s, 1, 1), collapse = ""))
    })
  }
  axis.names <- sapply(names(res$X), makeAcronym) ## see pgx-functions

  ## -------------------------------------------------------------
  ## Finally do the plotting
  ## -------------------------------------------------------------
  hpd$nodes$size <- cex * hpd$nodes$size
  hpd$edges$weight <- cex * hpd$edges$weight
  mr <- max(hpd$nodes$radius)
  plotHive(hpd,
    ch = 5, bkgnd = "white",
    axLabs = axis.names,
    axLab.pos = c(1, 1.2, 1.2) * 0.15 * mr, #

    axLab.gpar = grid::gpar(
      col = "black", fontsize = 18 * cex,
      lwd = 4, fontface = "bold"
    )
  )

  tt <- paste("edge.width = edge.importance",
    "node.size = variable importance",
    "axis = fold-change or F-stat",
    sep = "\n"
  )
  grid::grid.text(tt,
    x = 0.2 * mr, y = -1.0 * mr, default.units = "native",
    just = "left", gp = grid::gpar(fontsize = 9, col = "black")
  )

  ## axis 1
  rot.xy <- function(x, y, deg) {
    a <- deg * pi / 180
    rx <- cos(a) * x - sin(a) * y
    ry <- sin(a) * x + cos(a) * y
    cbind(rx, ry)
  }
  rot <- c(0, 120, 240)
  mr <- max(hpd$nodes$radius)
  yoff <- c(0, -0, +0)

  k <- 1
  for (k in 1:3) {
    kk <- which(hpd$nodes$axis == k)
    if (showloops) {
      kk <- which(hpd$nodes$axis == k & hpd$nodes$looping)
    }
    jj <- Matrix::head(kk[order(-hpd$nodes$size[kk])], numlab) ## number of labels
    rr <- hpd$nodes$radius
    rx <- rot.xy(0, rr[jj] + 5, rot[k])


    lab <- sub(".*:", "", hpd$nodes$lab[jj])
    ##    pt <- maptools::pointLabel(rx[, 1], rx[, 2], labels = lab, cex = cex * 2, doPlot = FALSE)
    ##    px <- cbind(pt$x, pt$y)
    px <- cbind(rx[, 1], rx[, 2])
    px[, 1] <- px[, 1] + 4
    grid::grid.text(lab,
      ## x = 10 + rx[,1], y = rx[,2],
      x = px[, 1], y = px[, 2],
      default.units = "native", just = "left",
      gp = grid::gpar(fontsize = 12 * cex, col = "black")
    )

    grid::grid.segments(rx[, 1], rx[, 2], px[, 1] - 1, px[, 2],
      default.units = "native"
    )
  }
}
