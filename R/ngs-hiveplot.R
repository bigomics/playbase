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
