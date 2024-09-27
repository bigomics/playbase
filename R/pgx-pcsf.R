##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

#' @title Compute PCSF solution from pgx object
#'
#' @description
#' Compute PCSF solution from pgx object for given contrast.
#'
#' @export
pgx.computePCSF <- function(pgx, contrast, level = "gene",
                            ntop = 250, ncomp = 3, beta = 1,
                            use.corweight = TRUE, rm.negedge = TRUE,
                            dir = "both") {
  F <- pgx.getMetaMatrix(pgx, level = level)$fc
  if (is.null(contrast)) {
    fx <- rowMeans(F**2, na.rm = TRUE)**0.5
  } else {
    if (!contrast %in% colnames(F)) {
      stop("[pgx.computePCSF] invalid contrast")
    }
    fx <- F[, contrast]
    names(fx) <- rownames(F)
  }
  if (level == "geneset") {
    fx <- fx[grep("^GOBP:", names(fx))] ## really?
    ## names(fx) <- gsub(".*:| \\(.*", "", names(fx))
  }

  if (dir == "both") {
    sel1 <- head(order(fx), ntop / 2)
    sel2 <- head(order(-fx), ntop / 2)
    sel <- c(sel1, sel2)
    fx <- fx[sel]
  }
  if (dir == "up") {
    sel <- head(order(-fx), ntop)
    fx <- fx[sel]
  }
  if (dir == "down") {
    sel <- head(order(fx), ntop)
    fx <- fx[sel]
  }

  fx[is.na(fx)] <- 0

  ## first pass
  nodes <- names(fx)
  get_edges <- function(nodes, use.corweight) {
    if (level == "gene") {
      data(STRING, package = "PCSF")
      sel <- (STRING$from %in% nodes & STRING$to %in% nodes)
      ee <- STRING[sel, ]
      if (use.corweight) {
        R <- cor(t(pgx$X[nodes, ]))
        if (rm.negedge) R[which(R < 0)] <- NA
        wt <- (1 - R[cbind(ee$from, ee$to)])
        ee$cost <- ee$cost * wt
      }
    }
    if (level == "geneset") {
      G <- pgx$GMT[, nodes]
      H <- cor(as.matrix(G))
      diag(H) <- 0
      jj <- which(abs(H) > 0.5, arr.ind = TRUE)
      ee <- data.frame(
        from = rownames(H)[jj[, 1]],
        to = rownames(H)[jj[, 2]],
        cost = (1 - H[jj])
      )
      if (use.corweight) {
        gg <- unique(c(ee$from, ee$to))
        R <- cor(t(pgx$gsetX[gg, ]))
        if (rm.negedge) R[which(R < 0)] <- NA
        wt <- (1 - R[cbind(ee$from, ee$to)])
        ee$cost <- ee$cost * wt
      }
    }
    ee
  }

  ee <- get_edges(names(fx), use.corweight = use.corweight)
  ee <- ee[!is.na(ee$cost), ]

  suppressMessages(suppressWarnings(
    ppi <- PCSF::construct_interactome(ee)
  ))
  prize1 <- abs(fx[igraph::V(ppi)$name])

  suppressMessages(suppressWarnings(
    pcsf <- try(PCSF::PCSF(ppi,
      terminals = prize1, w = 2,
      b = beta, mu = 5e-04, verbose = 0
    ))
  ))

  if (inherits(pcsf, "try-error")) {
    message("[pgx.computePCSF] WARNING: No solution. Please adjust beta or mu.")
    return(NULL)
  }

  igraph::V(pcsf)$foldchange <- fx[igraph::V(pcsf)$name]

  ## remove small clusters...
  cmp <- igraph::components(pcsf)
  if (ncomp < 1) {
    sel.kk <- which(cmp$csize > ncomp * max(cmp$csize, na.rm = TRUE))
  } else {
    sel.kk <- head(order(-cmp$csize), ncomp)
  }
  pcsf <- igraph::subgraph(pcsf, cmp$membership %in% sel.kk)
  class(pcsf) <- c("PCSF", "igraph")

  return(pcsf)
}


#' @title Plot PCSF graph from results object
#'
#' @description Plot PCSF graph (visnetwork or igraph) from result of
#'   pgx.computePCSF
#'
#' @export
plotPCSF <- function(pcsf,
                     highlightby = c("centrality", "prize")[1],
                     plotlib = c("visnet", "igraph")[1],
                     layout = "layout_with_kk", physics = TRUE,
                     node_cex = 30, label_cex = 30, nlabel = -1) {
  ## set node size
  fx <- igraph::V(pcsf)$foldchange
  wt <- abs(fx / mean(abs(fx), na.rm = TRUE))**0.7
  node_cex1 <- node_cex * pmax(wt, 1)
  node_cex1
  label_cex1 <- label_cex + 1e-8 * abs(fx)

  ## set colors as UP/DOWN
  igraph::V(pcsf)$type <- c("down", "up")[1 + 1 * (sign(fx) > 0)]

  ## set label size
  if (highlightby == "centrality") {
    wt <- igraph::E(pcsf)$weight
    ewt <- 1.0 / (0.01 * mean(wt, na.rm = TRUE) + wt) ##
    bc <- igraph::page_rank(pcsf, weights = ewt)$vector
    label_cex1 <- label_cex * (1 + 3 * (bc / max(bc, na.rm = TRUE))**3)
  }
  if (highlightby == "prize") {
    fx1 <- abs(fx)
    label_cex1 <- label_cex * (1 + 3 * (fx1 / max(fx1, na.rm = TRUE))**3)
  }

  igraph::V(pcsf)$label <- igraph::V(pcsf)$name
  if (nlabel > 0) {
    top.cex <- head(order(-label_cex1), nlabel)
    bottom.cex <- setdiff(1:length(label_cex1), top.cex)
    igraph::V(pcsf)$label[bottom.cex] <- ""
  }

  out <- NULL
  if (plotlib == "visnet") {
    library(igraph)
    bluered <- playdata::BLUERED(7)
    out <- visplot.PCSF(
      pcsf,
      style = 1,
      node_size = node_cex1,
      node_label_cex = label_cex1,
      invert.weight = TRUE,
      edge_width = 5,
      Steiner_node_color = "lightblue",
      Terminal_node_color = "lightgreen",
      extra_node_colors = list(
        "down" = bluered[2],
        "up" = bluered[6]
      ),
      edge_color = "lightgrey",
      width = "100%", height = 900,
      layout = layout,
      physics = physics
    )
  }

  if (plotlib == "igraph") {
    plotPCSF.IGRAPH(pcsf, fx0 = NULL, label.cex = label_cex1)
    out <- NULL
  }

  return(out)
}


#' @title Plot interactive visnetwork graph from PCSF results object
#'
#' @description Plot interactive visnetwork PCSF graph from result of
#'   pgx.computePCSF. To be called from plotPCSF.
#'
#' @return visNetwork plot object
#' @export
visplot.PCSF <- function(
    net, style = 0, edge_width = 5, node_size = 40, node_label_cex = 30,
    Steiner_node_color = "lightblue", Terminal_node_color = "lightgreen",
    Terminal_node_legend = "Terminal", Steiner_node_legend = "Steiner",
    layout = "layout_with_kk", physics = TRUE, layoutMatrix = NULL,
    width = 1800, height = 1800, invert.weight = FALSE,
    extra_node_colors = list(), edge_color = "lightgrey", ...) {
  if (missing(net)) {
    stop("Need to specify the subnetwork obtained from the PCSF algorithm.")
  }
  if (class(net)[1] != "PCSF" || class(net)[2] != "igraph") {
    stop("The subnetwork must be a \"PCSF\" object derived from an \"igraph\" class.")
  }
  if (any(edge_width < 1)) {
    stop("The edge_width must be greater than 1.")
  }
  if (any(node_size < 10)) {
    stop("The node_size must be greater than 10.")
  }
  prize <- abs(igraph::V(net)$prize)
  min1 <- 10
  max1 <- node_size
  r1 <- max1 - min1
  min2 <- min(prize, na.rm = TRUE)
  max2 <- max(prize, na.rm = TRUE)
  r2 <- max2 - min2
  adjusted_prize <- r1 * (prize - min2) / r2 + min1

  weight <- igraph::E(net)$weight
  if (invert.weight) weight <- 1 / (weight + 1e-10)
  min1 <- 1
  max1 <- edge_width
  r1 <- max1 - min1
  min2 <- min(weight, na.rm = TRUE)
  max2 <- max(weight, na.rm = TRUE)
  r2 <- max2 - min2
  adjusted_weight <- r1 * (weight - min2) / r2 + min1

  ## prepare nodes/edge dataframes
  nodes <- data.frame(
    id = 1:length(igraph::V(net)),
    name = igraph::V(net)$name
  )
  nodes$label <- igraph::V(net)$label
  nodes$group <- igraph::V(net)$type
  nodes$size <- adjusted_prize
  nodes$title <- nodes$name
  nodes$label.cex <- node_label_cex
  nodes$font.size <- node_label_cex

  edges <- data.frame(igraph::ends(net, es = igraph::E(net)), adjusted_weight)
  names(edges) <- c("from", "to", "value")
  edges$from <- match(edges$from, nodes$name)
  edges$to <- match(edges$to, nodes$name)

  visNet <- visNetwork::visNetwork(
    nodes, edges,
    width = width,
    height = height, ...
  ) %>%
    visNetwork::visNodes(shadow = list(enabled = TRUE, size = 12)) %>%
    visNetwork::visEdges(
      color = edge_color,
      scaling = list(min = 6, max = edge_width * 6)
    )

  if (length(extra_node_colors) > 0) {
    for (i in 1:length(extra_node_colors)) {
      en <- names(extra_node_colors)[i]
      visNet <- visNet %>% visNetwork::visGroups(
        groupname = en,
        color = list(
          background = extra_node_colors[[en]],
          border = "lightgrey"
        )
      )
      #      leg.groups[[i]] <- list(label = en, shape = "triangle",
      #        size = 13, color = list(background = extra_node_colors[[en]],
      #          border = "grey"), label.cex = 0.8)
    }
  }

  if (layout != "hierarchical") {
    visNet <- visNet %>%
      visNetwork::visIgraphLayout(
        layout = layout,
        physics = physics,
        layoutMatrix = layoutMatrix
      )
  }
  if (layout == "hierarchical") {
    visNet <- visNet %>% visNetwork::visHierarchicalLayout(direction = "UD")
  }

  visNet <- visNet %>%
    visNetwork::visPhysics(enabled = physics) %>%
    visNetwork::visOptions(highlightNearest = list(enabled = TRUE, degree = 2)) %>%
    visNetwork::visLegend(width = 0.05, useGroups = TRUE, position = "right")

  visNet
}

#' @title Plot static igraph network from PCSF results object
#'
#' @description Plot static igraph PCSF graph from result of
#'   pgx.computePCSF. To be called from plotPCSF.
#'
#' @return NULL
#' @param x
#'
plotPCSF.IGRAPH <- function(net, fx0 = NULL, label.cex = 1) {
  if (is.null(fx0)) fx0 <- igraph::V(net)$prize
  ## fx0 <- tanh(1.3 * fx0)
  fx0 <- fx0 / max(abs(fx0), na.rm = TRUE)
  label.cex <- label.cex / mean(label.cex, na.rm = TRUE)
  vertex.label.cex <- 0.8 * label.cex**0.80

  vertex.size <- 1.2 * igraph::V(net)$prize**0.9
  vertex.color <- "lightblue"
  cpal <- colorRampPalette(c("blue2", "grey90", "red3"))(33)
  vertex.color <- cpal[1 + 16 + round(16 * fx0)]
  edge.width <- (1 - igraph::E(net)$weight / max(igraph::E(net)$weight, na.rm = TRUE))

  pos <- igraph::layout_with_graphopt(
    net,
    niter = 5000,
    charge = 0.00001,
    mass = 30,
    spring.length = 1,
    spring.constant = 1
  )

  vv <- (vertex.size / mean(vertex.size, na.rm = TRUE))**1
  plot(
    net,
    vertex.size = vertex.size,
    vertex.color = vertex.color,
    vertex.label.cex = vertex.label.cex,
    vertex.label.dist = 0.3 + 0.7 * vv,
    vertex.label.degree = -0 * pi,
    vertex.label.family = "sans",
    edge.width = 5 * edge.width,
    layout = pos
  )
}

#' @export
pgx.getPCSFcentrality <- function(pgx, contrast, pcsf = NULL, plot = TRUE, n = 10) {
  if (is.null(pcsf)) {
    pcsf <- pgx.computePCSF(pgx, contrast, level = "gene", ntop = 250, ncomp = 3)
  }

  ## centrality
  ewt <- 1.0 / igraph::E(pcsf)$weight
  cc <- igraph::page_rank(pcsf, weights = ewt)$vector
  cc <- cc / mean(cc, na.rm = TRUE) ## normalize
  tail(sort(cc), 20)
  top.cc <- sort(cc, decreasing = TRUE)
  M <- pgx$gx.meta$meta[[contrast]]
  sel <- intersect(names(top.cc), rownames(M))
  sel <- head(sel, n)
  fc <- M[sel, c("meta.fx")]
  M <- data.frame(centrality = top.cc[sel], logFC = fc)
  M <- round(M, digits = 4)
  aa <- pgx$genes[rownames(M), c("gene_name", "gene_title")]
  aa <- cbind(aa, M)
  aa$gene_title <- substring(aa$gene_title, 1, 50)
  rownames(aa) <- NULL

  if (plot) {
    ## Plot your table with table Grob in the library(gridExtra)
    tab <- gridExtra::tableGrob(aa, rows = NULL)
    gridExtra::grid.arrange(tab)
  }

  return(aa)
}
