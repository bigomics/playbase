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
pgx.computePCSF <- function(pgx, contrast, datatypes = NULL,
                            as_prize = "fc", use_rank = FALSE,
                            use_rank2 = FALSE,                            
                            ntop = 250, ncomp = 3, beta = 1,
                            rm.negedge = TRUE, highcor = 0.8, 
                            dir = "both", level="gene",
                            ppi = c("STRING", "GRAPHITE")
                            ) {

  if(level!="gene") {
    message("[pgx.computePCSF] WARNING. 'level' is deprecated. Please use pgx.computePCSF_gset for geneset level PCSF.")
  }
  
  ## Compute correlation and foldchange (wrt pheno)
  Y <- pgx$model.parameters$exp.matrix[,]
  Y[which(Y==0)] <- NA
  R <- cor(Matrix::t(pgx$X), Y, use="pairwise")
  F <- pgx.getMetaMatrix(pgx, level="gene")$fc
  X <- pgx$X
  
  ## determine prize
  if(as_prize == "rho") {
    message("using rho as prize")
    Z <- R
  } else if(as_prize %in% c("rho.fc","fc.rho")) {
    message("using rho.fc as prize")
    Z <- sqrt(abs(R*F))*(sign(R)==sign(F))*sign(F)
  } else {
    ## foldchange as prize
    message("using fc as prize")    
    Z <- F
  }

  zx=NULL
  if (is.null(contrast)) {
    zx <- sqrt(rowMeans(Z**2, na.rm = TRUE))
  } else if (contrast %in% colnames(Z)) {
    zx <- Z[,contrast]
  } else {
    stop("[pgx.computePCSF] invalid contrast")
  }
  
  if(use_rank) {
    zx <- signedRank(zx)
  }
  
  ## align just to be sure
  gg <- intersect(rownames(X), names(zx))
  X <- X[gg,,drop=FALSE ]
  F <- F[gg,,drop=FALSE ]
  R <- R[gg,,drop=FALSE ]
  zx <- zx[gg]
  labels <- gg
  names(labels) <- gg
  
  ## filter on datatype or geneset collection
  if(!is.null(datatypes)) {
    dt <- mofa.get_prefix(rownames(X))
    datatypes <- intersect(datatypes, unique(dt))
    if(length(datatypes)==0) {
      message("[pgx.computePCSF2] ERROR. no valid datatypes")
      return(NULL)
    }
  }

  ## take top features by datatype
  if(!is.null(datatypes)) {
    dt <- mofa.get_prefix(names(zx))    
    ii <- tapply(zx, dt, function(x) head(names(sort(-abs(x))),2*ntop))
    ii <- unlist(ii[datatypes])
  } else {
    ii <- head(names(sort(-abs(zx))),2*ntop)
  }
  X <- X[ii,,drop=FALSE ]
  F <- F[ii,,drop=FALSE ]
  R <- R[ii,,drop=FALSE ]
  zx <- zx[ii]
  labels <- labels[ii]

  if(use_rank2) {
    zx <- signedRank(zx)
  }
  
  ##-----------------------------------
  ## Setup PPI
  ##-----------------------------------
  ppix <- c()
  if(is.null(ppi)) {
    message("skipping PPI edges")
  }
  if(!is.null(ppi)) {
    ppi.ismatrix <- is.data.frame(ppi) || is.matrix(ppi)
    ppi.ismatrix
    if (ppi.ismatrix) {
      ppix <- ppi
    } else if(all(is.character(ppi))) {
      if ("STRING" %in% ppi) {
        message("adding STRING PPI...")
        ppix <- rbind(ppix, PCSF::STRING)
      }
      if ("GRAPHITE" %in% ppi) {
        message("adding GRAPHITE PPI...")
        ppix <- rbind(ppix, playdata::GRAPHITE_PPI)
      }
    } else {
      stop("invalid PPI")
    }
    
    ## For now, we can strip prefixes in PPI
    ppix[, 1:2] <- apply(ppix[, 1:2], 2, function(x) sub(".*:", "", x))
    ppix.features <- unique(c(ppix[,1],ppix[,2]))

    ## convert PPI human symbol to organism symbol
    annot <- pgx$genes
    if(!is.null(datatypes) && !is.null(annot$data_type)) {
      sel <- which(annot$data_type %in% datatypes)
      annot <- annot[sel,,drop=FALSE]
    }
    nn.match <- apply(annot, 2, function(a) sum(ppix.features %in% a))
    nn.match
    match.col <- names(which.max(nn.match))
    ii <- match(ppix[,1], annot[,match.col])
    jj <- match(ppix[,2], annot[,match.col])
    ppix[,1] <- rownames(annot)[ii]
    ppix[,2] <- rownames(annot)[jj]
    kk <- which(!is.na(ppix[,1]) & !is.na(ppix[,2]))
    ppix <- ppix[kk,,drop=FALSE]
    dim(ppix)

    ## take all edges. disregard cost.
    ppix$cost <- 0.01 ### ???
  }
  
  ## Add high-correlation edges
  if( !is.null(highcor) && highcor <= 1) {
    corX <- cor(t(X), use = "pairwise")
    diag(corX) <- 0
    hce <- which( abs(corX) > highcor, arr.ind=TRUE)
    if(nrow(hce)>0) {
      message("adding ",nrow(hce)," very-high-correlation edges rho>", highcor)
      hce.ppi <- data.frame(
        from = rownames(corX)[hce[,1]],
        to = rownames(corX)[hce[,2]],
        cost = 0.1
      )
      ppix <- rbind(ppix, hce.ppi)
    }
  }
  dim(ppix)
  
  ## compute PCSF
  pcsf <- solvePCSF(
    X,
    fx = zx,  
    ppi = ppix,
    labels = labels,
    ntop = ntop,
    ncomp = ncomp,
    beta = beta,
    rm.negedge = rm.negedge,
    dir = dir
  )
  
  ## Reinstate foldchange to actual foldchange
  ## save fc & rho
  vname <- igraph::V(pcsf)$name
  if(!is.null(contrast)) {
    igraph::V(pcsf)$rho <- R[vname,contrast]
    igraph::V(pcsf)$foldchange <- F[vname,contrast]
  } else {
    rx <- sqrt(rowMeans(R[,]**2, na.rm=TRUE))
    fx <- sqrt(rowMeans(F[,]**2, na.rm=TRUE))
    igraph::V(pcsf)$rho <- rx[vname]
    igraph::V(pcsf)$foldchange <- fx[vname]
  }
  
  return(pcsf)
}


#' @title Compute PCSF solution from pgx object for GENESET level
#'
#' @export
pgx.computePCSF_gset <- function(pgx, contrast,
                                 as_prize = "fc", use_rank = FALSE,
                                 use_rank2 = FALSE,
                                 gmt.rho = 0.8, gset.filter = NULL,
                                 highcor = 0.9, ntop = 500,
                                 ncomp = 10, beta = 1, 
                                 rm.negedge = TRUE,
                                 dir = "both") {

  ##gmt.rho=0.8;gset.filter=NULL;highcor=0.9;ntop=250;ncomp=3;beta=1;rm.negedge=TRUE;dir="both"
  
  ## Compute correlation and foldchange (wrt pheno)
  Y <- pgx$model.parameters$exp.matrix[,]
  Y[which(Y==0)] <- NA
  R <- cor(Matrix::t(pgx$gsetX), Y, use="pairwise")
  F <- pgx.getMetaMatrix(pgx, level="geneset")$fc
  X <- pgx$gsetX
  
  ## determine prize
  if(as_prize == "rho") {
    message("using rho as prize")
    Z <- R
  } else if(as_prize == "rho.fc") {
    message("using rho.fc as prize")
    Z <- sqrt(abs(R*F))*(sign(R)==sign(F))*sign(F)
  } else {
    ## foldchange as prize
    message("using fc as prize")    
    Z <- F
  }

  zx=NULL
  if (is.null(contrast)) {
    zx <- sqrt(rowMeans(Z**2, na.rm = TRUE))
  } else if (contrast %in% colnames(Z)) {
    zx <- Z[,contrast]
  } else {
    stop("[pgx.computePCSF] invalid contrast")
  }
  
  if(use_rank) {
    zx <- signedRank(zx)
  }

  ## align just to be sure
  gg <- intersect(rownames(X), names(zx))
  X <- X[gg,,drop=FALSE ]
  F <- F[gg,,drop=FALSE ]
  R <- R[gg,,drop=FALSE ]
  zx <- zx[gg]
  
  ## filter on datatype or geneset collection
  if(!is.null(gset.filter)) {
    sel <- grep(gset.filter, names(zx))
    zx <- zx[sel]
  }  
  ## take 2*ntop 
  ii <- head(names(sort(-abs(zx))),2*ntop)
  X <- X[ii,,drop=FALSE ]
  F <- F[ii,,drop=FALSE ]
  R <- R[ii,,drop=FALSE ]
  zx <- zx[ii]

  if(use_rank2) {
    zx <- signedRank(zx)
  }
  
  ## Construct base PPI using GMT overlap
  aa <- intersect(names(zx), colnames(pgx$GMT))
  G <- (pgx$GMT[, aa] != 0)
  H <- qlcMatrix::corSparse(G)
  dim(H)
  diag(H) <- 0
  jj <- which(H > gmt.rho, arr.ind = TRUE)
  ppix <- c()
  if(length(jj)) {
    message("adding ",length(jj)," GMT edges at gmt.rho>", gmt.rho)
    ppix <- data.frame(
      from = colnames(G)[jj[, 1]],
      to = colnames(G)[jj[, 2]],
      cost = pmax(1 - H[jj], 0)
    )
  }
  
  ## Add high-correlation edges
  dim(X)
  corX <- cor(t(X), use = "pairwise")
  diag(corX) <- 0
  hce <- which(abs(corX) > highcor, arr.ind=TRUE)
  if(nrow(hce) > 0) {
    message("adding ",nrow(hce)," very-high-correlation edges rho>", highcor)
    ppi.hce <- data.frame(
      from = rownames(corX)[hce[,1]],
      to = rownames(corX)[hce[,2]],
      cost = 0.1
    )
    ppix <- rbind(ppix, ppi.hce)
  }
  message("PPI number of edges: n=", nrow(ppix))
  
  ## compute PCSF
  pcsf <- solvePCSF(
    X,
    fx = zx,
    ppi = ppix,
    ntop = ntop,
    ncomp = ncomp,
    beta = beta,
    rm.negedge = rm.negedge,
    dir = dir
  )
  
  ## Reinstate foldchange to actual foldchange
  ## save fc & rho
  vname <- igraph::V(pcsf)$name
  if(!is.null(contrast)) {
    igraph::V(pcsf)$rho <- R[vname,contrast]
    igraph::V(pcsf)$foldchange <- F[vname,contrast]
  } else {
    rx <- sqrt(rowMeans(R[,]**2, na.rm=TRUE))
    fx <- sqrt(rowMeans(F[,]**2, na.rm=TRUE))
    igraph::V(pcsf)$rho <- rx[vname]
    igraph::V(pcsf)$foldchange <- fx[vname]
  }

  return(pcsf)
}


solvePCSF <- function(X, fx, ppi, labels = NULL, ntop = 250, ncomp = 3,
                      beta = 1, rm.negedge = TRUE, dir = "both") {
  if (!is.null(labels)) names(labels) <- rownames(X)

  ppi.genes <- unique(c(ppi$from, ppi$to))  
  ppi.ratio <- mean(names(fx) %in% ppi.genes)
  ppi.num <- sum(names(fx) %in% ppi.genes)
  if(ppi.ratio < 0.10 || ppi.num < 10) {
    if(ppi.ratio < 0.10) message("[solvePCSF] WARNING: less than 10% genes are in PPI. ")
    if(ppi.num < 10) message("[solvePCSF] WARNING: less than 10 genes are in PPI. ")
    message("[solvePCSF] ERROR: Exiting. ")
    return(NULL)    
  }
  fx <- fx[which(names(fx) %in% ppi.genes)]
  
  if (dir == "both" && length(fx)>ntop) {
    sel1 <- head(order(fx), ntop / 2)
    sel2 <- head(order(-fx), ntop / 2)
    sel <- unique(c(sel1, sel2))
    fx <- fx[sel]
  }
  if (dir == "up" && length(fx)>ntop) {
    sel <- head(order(-fx), ntop)
    fx <- fx[sel]
  }
  if (dir == "down" && length(fx)>ntop) {
    sel <- head(order(fx), ntop)
    fx <- fx[sel]
  }
  fx[is.na(fx)] <- 0

  ## first pass
  nodes <- names(fx)
  get_edges <- function(nodes) {
    sel <- (ppi$from %in% nodes & ppi$to %in% nodes)
    ee <- ppi[which(sel), ]
    selx <- rownames(X) %in% union(ee$from, ee$to)
    table(selx)
    R <- cor(t(X[selx, , drop = FALSE]), use = "pairwise")
    if (rm.negedge) R[which(R < 0)] <- NA
    wt <- (1 - R[cbind(ee$from, ee$to)])
    ee$cost <- ee$cost * wt
    ee
  }
  ee <- get_edges(names(fx))
  ee <- ee[!is.na(ee$cost), ]
  
  # 1. Construct PPI interactome
  ##suppressMessages(suppressWarnings(xppi <- PCSF::construct_interactome(ee)))
  xppi <- igraph::graph_from_edgelist(as.matrix(ee[,1:2]), directed=FALSE)
  igraph::E(xppi)$weight <- as.numeric(ee[,3])
  xppi <- igraph::simplify(xppi, edge.attr.comb = "min") ## important!!!

  # 2. Solve the PCSF
  prize1 <- abs(fx[igraph::V(xppi)$name])
  suppressMessages(suppressWarnings(
    pcsf <- try(PCSF::PCSF(xppi,
      terminals = prize1, w = 2,
      b = beta, mu = 5e-04, verbose = 0
    ))
  ))

  if (inherits(pcsf, "try-error")) {
    message("[pgx.solvePCSF] WARNING: No solution. Please adjust beta or mu.")
    return(NULL)
  }

  # set foldchange
  igraph::V(pcsf)$foldchange <- fx[igraph::V(pcsf)$name]
  if (!is.null(labels)) {
    igraph::V(pcsf)$label <- labels[igraph::V(pcsf)$name]
  }

  ## determine node datatype (SYMBOL or CHEBI)
  is.multiomics <- all(grepl(":",igraph::V(pcsf)$name))
  if(is.multiomics) {
    dtype <- sub(":.*","",igraph::V(pcsf)$name)
    igraph::V(pcsf)$type <- dtype
  } else {
    igraph::V(pcsf)$type <- ""
  }

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
                     sizeby = "prize",
                     colorby = "foldchange",
                     highlightby = c("centrality", "prize")[1],
                     labels = NULL,
                     plotlib = c("visnet", "igraph")[1],
                     layout = c("layout_with_kk","hierarchical")[1],
                     layoutMatrix = NULL,
                     asp = 1,
                     physics = TRUE,
                     node_cex = 1,
                     node_alpha = 0.8,
                     node_gamma = 1,
                     label_cex = 1,
                     label_gamma = 2,                     
                     nlabel = -1,
                     border_width = 1.5,
                     edge_cex = 1,
                     edge_length = NULL,
                     cut.clusters = FALSE,
                     cut.layout = "kk",
                     cut.resolution = 0.05,
                     as_grid = TRUE,
                     nlargest = -1
                     ) {


  ## determine node size vector
  fx <- NULL
  if(is.character(sizeby) && sizeby %in% names(igraph::vertex_attr(pcsf))) {
    fx <- igraph::vertex_attr(pcsf, sizeby)
  } else if(is.numeric(sizeby)) {
    fx <- sizeby
  } else {
    fx <- igraph::V(pcsf)$foldchange
  }
  names(fx) <- igraph::V(pcsf)$name
  wt <- abs(fx / mean(abs(fx), na.rm = TRUE))
  node_cex1 <- node_cex * (0.1 + wt)**node_gamma
  
  ## determine color vector
  color.fx <- NULL
  if(is.character(colorby) && colorby %in% names(igraph::vertex_attr(pcsf))) {
    color.fx <- igraph::vertex_attr(pcsf, colorby)
  } else if(is.numeric(colorby)) {
    color.fx <- colorby
  } else {
    color.fx <- igraph::V(pcsf)$foldchange
  }
  names(color.fx) <- igraph::V(pcsf)$name

  ## set colors for UP/DOWN. For VisNetwork they are set in
  ## 'extra_node_colors' in the visPlot.
  up_down <- c("down", "up")[1 + 1 * (sign(color.fx) > 0)]
  dtype <- igraph::V(pcsf)$type
  igraph::V(pcsf)$group <- paste0(dtype, "_", up_down)
  bluered <- playdata::BLUERED(7)
  bluered <- adjustcolor(bluered, alpha.f = node_alpha)
  groups <- sort(setdiff(igraph::V(pcsf)$group, c(NA)))
  ngroups <- length(groups)
  extra_node_colors <- c(bluered[2], bluered[6])[1 + grepl("_up$", groups)]
  names(extra_node_colors) <- groups

  ## set shape for different groups
  shapes <- head(rep(c(
    "dot", "triangle", "star", "diamond", "square",
    "triangleDown", "hexagon"
  ), 99), ngroups)
  levels <- names(sort(-table(dtype)))
  gtype <- gsub("_down|_up","",groups)
  extra_node_shapes <- shapes[factor(gtype, levels=levels)]
  names(extra_node_shapes) <- groups

  ## set bordercolors for different groups
  border_colors <- c("lightgrey",rep(c("yellow","magenta","cyan","green"),99))
  extra_node_borders <- border_colors[factor(gtype, levels=levels)]
  names(extra_node_borders) <- groups
  
  ## set label and label sizes
  if (!is.null(labels)) {
    if(!is.null(names(labels))) labels <- labels[igraph::V(pcsf)$name]
    igraph::V(pcsf)$label <- labels
  } else {
    igraph::V(pcsf)$label <- igraph::V(pcsf)$name
  }
  label_cex1 <- label_cex + 1e-8 * abs(fx)  
  if (highlightby == "centrality") {
    ewt <- igraph::E(pcsf)$weight
    wt <- 1.0 / (0.01 * mean(ewt, na.rm = TRUE) + ewt) ##
    bc <- (igraph::page_rank(pcsf, weights = wt)$vector)**0.8
    label_cex1 <- label_cex * (0.2 + 3 * (bc / max(bc, na.rm = TRUE))**label_gamma)
  }
  if (highlightby == "prize") {
    fx1 <- abs(fx)
    label_cex1 <- label_cex * (0.2 + 3 * (fx1 / max(fx1, na.rm = TRUE))**label_gamma)
  }
  if (highlightby %in% c("centrality.prize","prize.centrality")) {
    ewt <- igraph::E(pcsf)$weight
    wt <- 1.0 / (0.01 * mean(ewt, na.rm = TRUE) + ewt) ##
    bc <- (igraph::page_rank(pcsf, weights = wt)$vector)**0.8
    bc <- bc * abs(fx)
    label_cex1 <- label_cex * (0.2 + 3 * (bc / max(bc, na.rm = TRUE))**label_gamma)
  }
  names(label_cex1) <- igraph::V(pcsf)$name

  if(is.null(layoutMatrix)) {
    res <- pcsf.cut_and_relayout(
      pcsf,
      cut = cut.clusters,
      ncomp = nlargest,
      component.wise = cut.clusters,
      cluster.method = "leiden",
      leiden.resolution = cut.resolution,
      layout = cut.layout,
      as_grid = as_grid
    )
    pcsf <- res$graph
    layoutMatrix <- res$layout
  } 

  ## take n-largest components
  if(0 && nlargest > 0) {
    cc <- igraph::components(pcsf)$membership
    fsq <- tapply( igraph::V(pcsf)$foldchange**2, cc, sum)    
    top.comp <- head(names(sort(fsq,decreasing=TRUE)),nlargest)
    sel <- which(cc %in% top.comp)
    pcsf <- igraph::subgraph(pcsf, sel)
  }   
  
  ## this equalizes label size per component
  cc <- igraph::components(pcsf)$membership
  comp <- unique(cc)
  for(k in comp) {
    ii <- which(cc == k)
    #f <- mean(label_cex1) / mean(label_cex1[ii])
    f <-  median(label_cex1) / median(label_cex1[ii])
    label_cex1[ii] <- label_cex1[ii] * f
  }

  ## this limits number of labels per component
  if (nlabel == 0) {
    igraph::V(pcsf)$label <- ""
  } else if (nlabel > 0) {    
    for(k in comp) {
      ii <- which(cc == k)
      top.cex <- ii[head(order(-label_cex1[ii]), nlabel)]
      jj <- setdiff(ii, top.cex)
      igraph::V(pcsf)$label[jj] <- ""
    }
  }

  ## align
  vv <- igraph::V(pcsf)$name
  if(!is.null(layoutMatrix)) layoutMatrix <- layoutMatrix[vv,,drop=FALSE]
  fx <- fx[vv]
  color.fx <- color.fx[vv]
  node_cex1 <- node_cex1[vv]
  label_cex1 <- label_cex1[vv]
  
  out <- NULL
  if (plotlib == "visnet") {
    library(igraph)  ## ??
    class(pcsf) <- c("PCSF", "igraph")    
    if(!is.null(layoutMatrix)) {
      layoutMatrix[,2] <- -layoutMatrix[,2]
    }
    out <- visplot.PCSF(
      pcsf,
      style = 1,
      #node_size = 10 + 2 * node_cex1,
      node_size = pmax(4*node_cex1, 10),      
      node_label_cex = 40 * label_cex1,
      invert.weight = TRUE,
      edge_width = 5 * edge_cex,
      edge_length = edge_length,
      border_width = border_width,
      Steiner_node_color = "lightblue",
      Terminal_node_color = "lightgreen",
      extra_node_colors = extra_node_colors,
      extra_node_shapes = extra_node_shapes,
      extra_node_borders = extra_node_borders,
      edge_color = "lightgrey",
      width = "100%",
      height = 900,
      layout = layout,
      layoutMatrix = layoutMatrix,
      physics = physics
    )

  }

  if (plotlib == "igraph") {
    plotPCSF.IGRAPH(
      pcsf,
      fx = fx,  ## size
      color.fx = color.fx,
      label.fx = label_cex1,
      label.cex = 1.3 * label_cex,
      vertex.cex = node_cex,
      vertex.alpha = node_alpha,
      edge.cex = 1 * edge_cex,
      layoutMatrix = layoutMatrix,
      asp = asp
    )
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
    border_width = 1, edge_length = 10, edge_color = "lightgrey", 
    width = 1800, height = 1800, invert.weight = FALSE,
    extra_node_colors = c(), extra_node_shapes = c(), extra_node_borders = c(),
    ...) {
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
  nodes$group <- igraph::V(net)$group
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
    visNetwork::visNodes(
      shadow = list(enabled = TRUE, size = 12),
      borderWidth = border_width
    ) %>%
    visNetwork::visEdges(
      length = edge_length,
      color = edge_color,
      scaling = list(min = 6, max = edge_width * 6)
    )

  if (length(extra_node_colors) > 0) {
    num_extra_nodes <- length(extra_node_colors)
    if(length(extra_node_borders) == 0) {
      extra_node_borders <- rep("lightgrey",num_extra_nodes)
      names(extra_node_borders) <- names(extra_node_colors)
    }
    for (i in 1:length(extra_node_colors)) {
      en <- names(extra_node_colors)[i]
      visNet <- visNet %>% visNetwork::visGroups(
        groupname = en,
        color = list(
          background = as.character(extra_node_colors[en]),          
          border = as.character(extra_node_borders[en])
        ),
        shape = as.character(extra_node_shapes[en])
      )
      #      leg.groups[[i]] <- list(label = en, shape = "triangle",
      #        size = 13, color = list(background = extra_node_colors[[en]],
      #          border = "grey"), label.cex = 0.8)
    }
  }

  if (layout != "hierarchical") {
    if(!is.null(layoutMatrix)) {
      layout <- "layout.norm"
      #physics <- FALSE
    }
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

  #  visNet <- visNet %>%
  #    visNetwork::visPhysics(enabled = physics) %>%
  #    visNetwork::visOptions(highlightNearest = list(enabled = TRUE, degree = 2)) %>%
  #    visNetwork::visLegend(width = 0.05, useGroups = TRUE, position = "right")

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
plotPCSF.IGRAPH <- function(net, fx = NULL, vertex.cex = 1,
                            color.fx = fx, vertex.alpha = 1,
                            label.fx = fx, label.cex = 4,
                            edge.cex = 1, 
                            layoutMatrix = NULL, ...) {

  if (is.null(fx)) fx <- igraph::V(net)$prize
  fx <- fx / max(abs(fx), na.rm = TRUE)
  vertex.size <- 5 * vertex.cex * (0.1 + abs(fx))**0.5

  color.fx <- color.fx / max(abs(color.fx),na.rm=TRUE)
  cpal <- colorRampPalette(c("blue2", "grey90", "red3"))(33)
  cpal <- adjustcolor(cpal, alpha.f = vertex.alpha)
  vertex.color <- cpal[1 + 16 + round(16 * color.fx)]

  label.fx <- label.fx / max(label.fx, na.rm = TRUE)
  vertex.label.cex <- label.cex * (0.2 + label.fx**2)
  if(label.cex==0) igraph::V(net)$label <- ''  
  ## somehow it is faster to bin label.cex
  vertex.label.cex <- round(10*vertex.label.cex)*0.1

  if(!is.null(layoutMatrix)) {
    if(is.null(rownames(layoutMatrix))) {
      stop("layoutMatrix must have rownames")
    }
  }

  if(is.null(layoutMatrix)) {
    g <- tapply( igraph::V(net), igraph::components(net)$membership,
      function(v) igraph::subgraph(net,v))
    g <- igraph::disjoint_union(g)
    pos <- igraph::layout_(g, igraph::with_kk(weights=NA), igraph::component_wise())
    rownames(pos) <- igraph::V(g)$name
    layoutMatrix <- pos[igraph::V(net)$name,,drop=FALSE]
  }

  ## convert cost to width: zero cost -> max width
  edge.width <- (1 - igraph::E(net)$weight / max(igraph::E(net)$weight, na.rm = TRUE))

  vnames <- igraph::V(net)$name
  if(!all(vnames %in% rownames(layoutMatrix))) {
    stop("ERROR: all nodes not in rownames layoutMatrix")
  }
  layoutMatrix <- layoutMatrix[vnames,,drop=FALSE]
    
  vv <- (vertex.size / mean(vertex.size, na.rm = TRUE))**1
  plot(
    net,
    vertex.size = vertex.size,
    vertex.color = vertex.color,
    vertex.frame.color = "grey50",
    vertex.frame.width = 0.3,    
    vertex.label.cex = vertex.label.cex,
    vertex.label.dist = 0.3 + 0.7 * vv,
    vertex.label.degree = -0 * pi,
    vertex.label.family = "sans",
    edge.width = 2 * edge.cex * edge.width,
    layout = layoutMatrix,
    ...
  )
  
}

#' @export
pgx.getPCSFcentrality <- function(pgx, contrast, pcsf, level="gene", 
                                  plot = TRUE, n = 10) {
  if (is.null(pcsf)) {
    stop("ERROR. must provide pcsf object")
  }

  ## centrality
  ewt <- igraph::E(pcsf)$weight
  wt <- 1.0 / (0.01 * mean(ewt, na.rm = TRUE) + ewt) ##
  cc <- igraph::page_rank(pcsf, weights = wt)$vector
  cc <- cc / (1e-8 + mean(cc, na.rm = TRUE)) ## normalize
  fc <- igraph::V(pcsf)$foldchange

  M <- data.frame(centrality = cc, logFC = fc)
  M <- round(M, digits = 2)
  rownames(M) <- igraph::V(pcsf)$name
  
  aa <- M[,0]
  if(level == "gene") {
    ft <- rownames(M)
    match.sum <- apply(pgx$genes, 2, function(a) sum(ft %in% a))
    if(all(match.sum==0)) {
      ft <- mofa.strip_prefix(rownames(M))
      match.sum <- apply(pgx$genes, 2, function(a) sum(ft %in% a))
    }
    ii <- match(ft, pgx$genes[,which.max(match.sum)])  
    aa <- pgx$genes[ii, c("feature", "symbol", "human_ortholog", "gene_title")]
  }
  if(level == "geneset") {
    aa <- data.frame( geneset = rownames(M))
  }
  aa <- cbind(aa, M)
  rownames(aa) <- rownames(M)
  aa <- aa[order(-aa$centrality), ]

  if (plot) {
    ## Plot your table with table Grob in the library(gridExtra)
    tab <- gridExtra::tableGrob(aa, rows = NULL)
    gridExtra::grid.arrange(tab)
  }

  return(aa)
}


#' Cut graph into components using Louvain clustering and calculate
#' merged layout.
#'
#' @export
pcsf.cut_and_relayout <- function(net, cut = TRUE,
                                  ncomp=-1, component.wise = TRUE,
                                  cluster.method = c("louvain","leiden")[1],
                                  leiden.resolution = 0.1,
                                  layout = c("kk","tree","circle","graphopt")[1],
                                  as_grid = TRUE) {

  if(cut == FALSE) {
    g <- net
  } else {
    ## we need component-wise layout
    mwt <- mean(igraph::E(net)$weight,na.rm=TRUE)
    ewt <- 1/(1e-3*mwt + igraph::E(net)$weight)
    if(cluster.method == "leiden") {
      clust <- igraph::cluster_leiden(
        net, weights=ewt, resolution = leiden.resolution)
    } else {
      clust <- igraph::cluster_louvain(net, weights=ewt)
    }
    ee <- igraph::E(net)[igraph::crossing(clust, net)]
    g <- igraph::delete_edges(net, ee)
    
    ## calculate component-wise layout
    g <- tapply( igraph::V(g), igraph::components(g)$membership,
      function(v) igraph::subgraph(g,v))
    
    if(ncomp>0) {
      fx <- igraph::V(net)$prize
      names(fx) <- igraph::V(net)$name
      gx <- sapply(g, function(x) sum(fx[igraph::V(x)$name]**2))
      g <- head( g[order(-gx)], ncomp )
    }
    
    ## create union and layout
    g <- igraph::disjoint_union(g)
  }
  
  layout <- head(layout,1)
  all.layouts = c("tree","circle","fr","kk","gem","graphopt")
  if(!layout %in% all.layouts) {
    layout <- "kk"
  }
  
  FUN <- switch(
    layout,
    tree = igraph::as_tree(),
    circle = igraph::in_circle(),
    fr = igraph::with_fr(weights=NA),
    kk = igraph::with_kk(weights=NA),
    mds = igraph::with_mds(),
    gem = igraph::with_gem(),    
    graphopt = igraph::with_graphopt(),
    igraph::with_kk(weights=NA)
  )

  if(component.wise) {
    pos <- igraph::layout_(g, FUN, igraph::component_wise())
  } else {
    pos <- igraph::layout_(g, FUN)
  }
  rownames(pos) <- igraph::V(g)$name

  ## Layout components on a grid
  if(as_grid) {
    mm <- igraph::components(g)$membership
    mm.order <- order(table(mm),decreasing=TRUE)
    ncomp <- length(unique(mm))
    nx <- ceiling(sqrt(ncomp))
    grid.pos <- t(sapply(0:(ncomp-1), function(i) c(i%%nx, i%/%nx)))
    if(layout=="tree") {
      dx <- tapply(pos[,1],mm,function(x) diff(range(x)))
      dy <- tapply(pos[,2],mm,function(y) diff(range(y)))    
      ydepth <- tapply(pos[,2], mm, function(y) length(unique(y)))
      ydepth
      k=1
      for(k in 1:length(unique(mm))) {
        ii <- which(mm == k)
        new_center <- grid.pos[k,,drop=FALSE]
        pos.ii <- pos[ii,,drop=FALSE]
        pos.ii[,2] <- as.integer(factor(rank(pos.ii[,2])))
        pos[ii,] <- pos.ii
      }
    }
    
    dx <- tapply(1:nrow(pos), mm, function(i) diff(range(pos[i,1])))
    dy <- tapply(1:nrow(pos), mm, function(i) diff(range(pos[i,2])))    
    row.dy <- tapply(dy, grid.pos[,2], max)
    row.dy    
    row.dy <- row.dy + 0.1*max(row.dy)  ## add gutter      
    dx <- 1.05 * max(dx)
    cx <- dx * grid.pos[,1]
    cy <- c(0,cumsum(row.dy))[1 + grid.pos[,2]]

    grid.pos <- cbind(cx, cy)
    grid.pos[,2] <- max(grid.pos[,2]) - grid.pos[,2]    
    k=1
    for(k in 1:length(unique(mm))) {
      ii <- which(mm == k)
      new_center <- grid.pos[k,]
      if(length(ii)==1) {
        pos[ii,] <- new_center
      } else {
        pos.mid <- colMeans(pos[ii,,drop=FALSE])      
        pos[ii,] <- t(t(pos[ii,]) - pos.mid + new_center)
      }
    }
  }
  
  list(graph=g, layout=pos)
}
