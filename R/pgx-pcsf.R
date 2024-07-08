##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


#' @export
pgx.computePCSF <- function(pgx, contrast, level = "gene", 
                            ntop = 250, ncomp = 3 ) {

  F <- pgx.getMetaMatrix(pgx, level = level)$fc
  if (is.null(contrast)) {    
    fx <- rowMeans(F**2)**0.5
  } else {
    if(!contrast %in% colnames(F)) {
      stop("[pgx.computePCSF] invalid contrast")
    }
    fx <- F[,contrast]
    names(fx) <- rownames(F)
  }
  if (level == "geneset") {
    fx <- fx[grep("^GOBP:", names(fx))]  ## really?
    ## names(fx) <- gsub(".*:| \\(.*", "", names(fx))
  }

  sel1 <- head(order(fx), ntop)
  sel2 <- head(order(-fx), ntop)
  sel <- c(sel1, sel2)
  fx <- fx[sel]
  fx[is.na(fx)] <- 0

  ## first pass
  get_edges <- function(nodes) {
    if (level == "gene") {
      genes <- names(fx)
      data(STRING, package = "PCSF")
      sel <- (STRING$from %in% nodes & STRING$to %in% nodes)
      ee <- STRING[sel, ]
    }
    if (level == "geneset") {
      G <- pgx$GMT[, nodes]
      R <- cor(as.matrix(G))
      diag(R) <- 0
      jj <- which(abs(R) > 0.5, arr.ind = TRUE)
      ee <- data.frame(
        from = rownames(R)[jj[, 1]],
        to = rownames(R)[jj[, 2]],
        cost = 1 - R[jj]
      )
    }
    ee
  }

  ee <- get_edges(names(fx))
  suppressMessages( suppressWarnings(
    ppi <- PCSF::construct_interactome(ee)
  ))
  prize1 <- abs(fx[igraph::V(ppi)$name])
  suppressMessages( suppressWarnings(
    pcsf <- PCSF::PCSF(ppi, terminals = prize1, w = 2, b = exp(.01), verbose = 0)
  ))

  ## remove small clusters...
  cmp <- igraph::components(pcsf)
  if (ncomp < 1) {
    sel.kk <- which(cmp$csize > ncomp * max(cmp$csize))
  } else {
    sel.kk <- head(order(-cmp$csize), ncomp)
  }
  pcsf <- igraph::subgraph(pcsf, cmp$membership %in% sel.kk)
  class(pcsf) <- c("PCSF", "igraph")
  
  return(pcsf)
}

#' @export
plotPCSF <- function(pcsf,
                     highlightby = c("centrality","prize")[1],
                     plot = c("visnet","igraph"),
                     node_cex = 30, label_cex = 30)
{
  ## set node size
  fx <- igraph::V(pcsf)$prize
  wt <- abs(fx / mean(abs(fx)))**0.8
  node_cex1 <- node_cex * pmax(wt, 1)
  node_cex1

  ## set colors as UP/DOWN
  igraph::V(pcsf)$type <- c("down", "up")[1 + 1 * (sign(fx) > 0)]

  do.physics <- TRUE
  layout <- "hierarchical"
  layout <- "layout_with_kk"

  ## set label size
  if (highlightby == "centrality") {
    wt <- igraph::E(pcsf)$weight
    ewt <- 1.0 / (0.01 * mean(wt) + wt) ##
    bc <- igraph::page_rank(pcsf, weights = ewt)$vector
    label_cex1 <- label_cex * (1 + 3 * (bc / max(bc))**3)
  }
  if (highlightby == "prize") {
    fx1 <- abs(fx)
    label_cex1 <- label_cex * (1 + 3 * (fx1 / max(fx1))**3)
  }

  ## ## set name
  ## if (level == "geneset") {
  ##   igraph::V(pcsf)$name <- gsub(".*:| \\(.*", "", igraph::V(pcsf)$name)
  ## }

  out <- NULL  
  if (plot == "visnet") {
    library(igraph)
    out <- visplot.PCSF(
      pcsf,
      style = 1,
      node_size = node_cex1,
      node_label_cex = label_cex1,
      invert.weight = TRUE,
      edge_width = 4,
      Steiner_node_color = "lightblue",
      Terminal_node_color = "lightgreen",
      extra_node_colors = list("down" = "blue", "up" = "red"),
      width = "100%", height = 900,
      layout = "layout_with_kk",
      physics = TRUE
    )
  }

  if (plot == "igraph") {
    plotPCSF.IGRAPH(pcsf, fx0 = NULL, label.cex = label_cex1)
    out <- NULL
  }

  return(out)
}


#' @param x
#'
#' @return
#' @export
visplot.PCSF <- function(
    net, style = 0, edge_width = 5, node_size = 40, node_label_cex = 30,
    Steiner_node_color = "lightblue", Terminal_node_color = "lightgreen",
    Terminal_node_legend = "Terminal", Steiner_node_legend = "Steiner",
    layout = "layout_with_fr", physics = TRUE, layoutMatrix = NULL,
    width = 1800, height = 1800, invert.weight = FALSE,
    extra_node_colors = list(), ...) {

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
  min2 <- min(prize)
  max2 <- max(prize)
  r2 <- max2 - min2
  adjusted_prize <- r1 * (prize - min2) / r2 + min1

  weight <- igraph::E(net)$weight
  if (invert.weight) weight <- 1 / (weight + 1e-10)
  min1 <- 1
  max1 <- edge_width
  r1 <- max1 - min1
  min2 <- min(weight)
  max2 <- max(weight)
  r2 <- max2 - min2
  adjusted_weight <- r1 * (weight - min2) / r2 + min1

  ## prepare nodes/edge dataframes
  nodes <- data.frame(1:length(igraph::V(net)), igraph::V(net)$name)
  names(nodes) <- c("id", "name")
  nodes$group <- igraph::V(net)$type
  nodes$size <- adjusted_prize
  nodes$title <- nodes$name
  nodes$label <- nodes$name
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
    visNetwork::visEdges(scaling = list(min = 6, max = edge_width * 6))

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

#' @param x
#'
#' @return
plotPCSF.IGRAPH <- function(net, fx0 = NULL, label.cex = 1) {

  if (is.null(fx0)) fx0 <- igraph::V(net)$prize
  ## fx0 <- tanh(1.3 * fx0)
  fx0 <- fx0 / max(abs(fx0),na.rm=TRUE)
  label.cex <- label.cex / mean(label.cex, na.rm = TRUE)
  vertex.label.cex <- 0.8 * label.cex**0.80

  vertex.size <- 1.2 * igraph::V(net)$prize**0.9
  vertex.color <- "lightblue"
  cpal <- colorRampPalette(c("blue2", "grey90", "red3"))(33)
  vertex.color <- cpal[1 + 16 + round(16 * fx0)]
  edge.width <- (1 - igraph::E(net)$weight / max(igraph::E(net)$weight))

  pos <- igraph::layout_with_graphopt(
    net,
    niter = 5000,
    charge = 0.00001,
    mass = 30,
    spring.length = 1,
    spring.constant = 1
  )

  vv <- (vertex.size / mean(vertex.size))**1
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
pgx.getPCSFcentrality <- function(pgx, contrast, pcsf = NULL, plot = TRUE, n=10) {

  if(is.null(pcsf)) {
    pcsf <- pgx.computePCSF(pgx, contrast, level = "gene", ntop = 250, ncomp = 3)
  }

  ## centrality
  ewt <- 1.0 / igraph::E(pcsf)$weight
  cc <- igraph::page_rank(pcsf, weights = ewt)$vector
  tail(sort(cc),20)
  top.cc <- sort(cc, decreasing = TRUE)  
  M <- pgx$gx.meta$meta[[contrast]]
  sel <- intersect(names(top.cc), rownames(M))
  sel <- head(sel, n)
  fc <- M[sel,c("meta.fx")]
  M <- data.frame(logFC = fc, centrality = top.cc[sel])
  M <- format(M, digits=3)
  aa <- pgx$genes[rownames(M),c("gene_name","gene_title")]
  aa <- cbind(aa, M)
  aa$gene_title <- substring(aa$gene_title, 1, 50)
  rownames(aa) <- NULL

  if(plot) {
    ##Plot your table with table Grob in the library(gridExtra)
    tab <- gridExtra::tableGrob(aa, rows=NULL)
    gridExtra::grid.arrange(tab)
  }
  
  return(aa)
}
