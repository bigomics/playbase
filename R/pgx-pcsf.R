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
                            rm.negedge = TRUE, use.corweight = TRUE,
                            dir = "both", ppi = c("STRING", "GRAPHITE"),
                            gset.rho = 0.8, gset.filter = NULL,
                            as.name = NULL) {

  ## get foldchange matrix
  F <- pgx.getMetaMatrix(pgx, level = level)$fc
  if (is.null(contrast)) {
    fx <- rowMeans(F**2, na.rm = TRUE)**0.5
  } else {
    if (!contrast %in% colnames(F)) {
      stop("[pgx.computePCSF] invalid contrast")
    }
    fx <- F[, contrast]
  }

  ## for multi-omics we must re-normalize by type otherwise FC can be
  ## quite different between datatypes.
  has.colons <- all(grepl(":", names(fx)))
  is.multiomics <- (pgx$datatype == "multi-omics") || has.colons
  is.multiomics
  if (is.multiomics) {
    dbg("[pgx.computePCSF] normalizing multi-omics logFC...")
    fx <- normalize_multifc(fx, by = "mad")
  }

  ## We need to harmonize organism symbol with the PPI database that
  ## is based on human ortholog. What node symbol/names do we use in
  ## the PCSF??
  if (level == "gene") {
    fx <- rename_by(fx, pgx$genes, "symbol") 
  }

  if (level == "geneset") {
    if (!is.null(gset.filter)) {
      fx <- fx[grep(gset.filter, names(fx))]
    }
  }

  if (level == "gene") {
    if (is.data.frame(ppi) || is.matrix(ppi)) {
      PPI <- ppi
    } else if (all(is.character(ppi))) {
      PPI <- c()
      if ("STRING" %in% ppi) {
        message("adding STRING PPI...")
        data(STRING, package = "PCSF")
        PPI <- rbind(PPI, STRING)
      }
      if ("GRAPHITE" %in% ppi) {
        message("adding GRAPHITE PPI...")
        PPI <- rbind(PPI, playdata::GRAPHITE_PPI)
      }
    } else {
      stop("invalid PPI")
    }
    PPI$cost <- 0.10 ### ???
  }

  if (level == "geneset") {
    aa <- intersect(names(fx), colnames(pgx$GMT))
    G <- pgx$GMT[, aa]
    ## H <- cor(as.matrix(G))
    H <- qlcMatrix::corSparse(G)
    dim(H)
    diag(H) <- 0
    jj <- which(H > gset.rho, arr.ind = TRUE)
    PPI <- data.frame(
      from = colnames(G)[jj[, 1]],
      to = colnames(G)[jj[, 2]],
      cost = pmax(1 - H[jj], 0)
    )
  }

  ## data matrix
  if (level == "gene") {
    X <- rename_by(pgx$X, pgx$genes, "symbol") ## safe
  }
  if (level == "geneset") {
    X <- pgx$gsetX
  }
  
  ## just to be sure
  gg <- intersect(rownames(X), names(fx))
  X <- X[gg, ]
  fx <- fx[gg]
  labels <- gg
  
  ## If multi-omics for metabolite set name as labels
  if (!is.null(as.name) && length(as.name) && as.name[1] != FALSE) {
    dt <- pgx$genes$data_type
    if (is.null(dt)) dt <- "gx"
    if (is.logical(as.name)) as.name <- unique(dt)
    if (any(dt %in% as.name)) {
      match.col <- which.max(apply(pgx$genes, 2, function(a) sum(gg %in% a)))
      ii <- match(labels, pgx$genes[, match.col])
      mx.name <- pgx$genes[ii, "gene_title"]
      dt <- pgx$genes[ii, "data_type"]
      jj <- which(dt %in% as.name)
      if (length(jj)) {
        labels[jj] <- paste0("[", mx.name[jj], "]")
      }
    }
  }

  ## all X and fx in symbol (HGNC or ChEBI) so we can strip of prefix
  ## in PPI
  PPI[, 1:2] <- apply(PPI[, 1:2], 2, function(x) sub(".*:", "", x))

  ## convert PPI human symbol to organism symbol
  ii <- match(PPI[,1], pgx$genes$human_ortholog)
  jj <- match(PPI[,2], pgx$genes$human_ortholog)
  PPI[,1] <- pgx$genes$symbol[ii]
  PPI[,2] <- pgx$genes$symbol[jj]
  kk <- which(!is.na(PPI[,1]) & !is.na(PPI[,2]))
  PPI <- PPI[kk,]

  ## compute PCSF
  pcsf <- solvePCSF(
    X, fx,
    ppi = PPI,
    labels = labels,
    ntop = ntop,
    ncomp = ncomp,
    beta = beta,
    rm.negedge = rm.negedge,
    dir = dir
  )

  return(pcsf)
}


#' @title Compute PCSF solution from pgx object for MULTI-OMICS
#'   (datatype prefixed) data.
#'
#' @export
pgx.computePCSF_multiomics <- function(pgx, contrast,
                                       datatypes = NULL,
                                       ntop = 250, ncomp = 3, beta = 1,
                                       rm.negedge = TRUE, highcor = 0.8, 
                                       dir = "both", ppi = c("STRING", "GRAPHITE"),
                                       as.name = NULL) {

  # For multi-omics we use correlation instead of fold-change as prize.
  y <- pgx$model.parameters$exp.matrix[,contrast]
  ii <- which(!is.na(y) & y!=0)
  rho <- cor(Matrix::t(pgx$X[,ii]), y[ii])[,1]  
  fc <- pgx.getMetaMatrix(pgx)$fc[,contrast]
  X <- pgx$X

  ## filter on datatype or geneset collection
  dt <- mofa.get_prefix(names(rho))
  table(dt)
  if(is.null(datatypes)) {
    datatypes <- unique(dt)
  }
  datatypes <- intersect(datatypes, unique(dt))

  if(length(datatypes)==0) {
    message("[pgx.computePCSF2] ERROR. no valid datatypes")
    return(NULL)
  }
  
  ## take 2*ntop 
  ii <- tapply(rho, dt, function(x) head(names(sort(-abs(x))),2*ntop))
  ii <- unlist(ii[datatypes])
  X  <- X[ii,,drop=FALSE]
  rho  <- rho[ii]
  fc  <- fc[ii]
  dim(X)
  
  ## We need to harmonize organism symbol with the PPI database that
  ## is based on human ortholog. What node symbol/names do we use in
  ## the PCSF??
  rho <- rename_by2(rho, pgx$genes, "symbol", keep.prefix=TRUE) 
  X   <- rename_by2(X,  pgx$genes, "symbol", keep.prefix=TRUE)
  fc  <- rename_by2(fc,  pgx$genes, "symbol", keep.prefix=TRUE) 
  
  ## align everything just to be sure
  gg <- intersect(rownames(X), names(rho))
  length(gg)
  X <- X[gg, ]
  rho <- rho[gg]
  fc <- fc[gg]
  labels <- gg
    
  ## Convert feature labels to gene_title if requested.
  if (!is.null(as.name) && length(as.name) && as.name[1] != FALSE) {
    dt <- pgx$genes$data_type
    if (is.null(dt)) dt <- "gx"
    if (is.logical(as.name)) as.name <- unique(dt)
    if (any(dt %in% as.name)) {
      labels0 <- mofa.strip_prefix(labels)
      ii <- match(labels0, pgx$genes[, "symbol"])
      mx.name <- pgx$genes$gene_title[ii]
      dt <- pgx$genes[ii, "data_type"]
      jj <- which(dt %in% as.name)
      if (length(jj)>0) {
        labels[jj] <- paste0(dt,":[", mx.name, "]")[jj]
        ##labels[jj] <- paste0("[", mx.name, "]")[jj]
      }
    }
  }

  ## Setup PPI
  ppi.ismatrix <- is.data.frame(ppi) || is.matrix(ppi)
  ppi.ismatrix
  if (ppi.ismatrix) {
    ppix <- ppi
  } else if(all(is.character(ppi))) {
    ppix <- c()
    if ("STRING" %in% ppi) {
      message("adding STRING PPI...")
      data(STRING, package = "PCSF")
      ppix <- rbind(ppix, STRING)
    }
    if ("GRAPHITE" %in% ppi) {
      message("adding GRAPHITE PPI...")
      ppix <- rbind(ppix, playdata::GRAPHITE_PPI)
    }
  } else {
    stop("invalid PPI")
  }
  ppix$cost <- 0.10 ### ???
  
  ## convert PPI human symbol to organism symbol
  ppix[, 1:2] <- apply(ppix[, 1:2], 2, function(x) sub(".*:", "", x))
  sel <- pgx$genes$data_type %in% c("mx","px")
  annot <- pgx$genes[sel,]
  ii <- match(ppix[,1], annot$human_ortholog)
  jj <- match(ppix[,2], annot$human_ortholog)
  dt.symbol <- paste0(annot$data_type,":",annot$symbol)
  ppix[,1] <- dt.symbol[ii]
  ppix[,2] <- dt.symbol[jj]
  kk <- which(!is.na(ppix[,1]) & !is.na(ppix[,2]))
  ppix <- ppix[kk,]
  dim(ppix)
  head(ppix)
  
  ## Add high-correlation edges
  if( !is.null(highcor) ) {
    dim(X)
    R <- cor(t(X), use = "pairwise")
    diag(R) <- 0
    hce <- which( abs(R) > highcor, arr.ind=TRUE)
    if(nrow(hce)>0) {
      message("adding ",nrow(hce)," very-high-correlation edges rho>", highcor)
      hce.ppi <- data.frame(
        from = rownames(R)[hce[,1]],
        to = rownames(R)[hce[,2]],
        cost = 0.1
      )
      ppix <- rbind(ppix, hce.ppi)
    }
  }
  dim(ppix)
  
  ## compute PCSF
  pcsf <- solvePCSF(
    X,
    fx = rho,
    ppi = ppix,
    labels = labels,
    ntop = ntop,
    ncomp = ncomp,
    beta = beta,
    rm.negedge = rm.negedge,
    dir = dir
  )
  
  ## Reinstate foldchange to actual foldchange
  V(pcsf)$foldchange <- fc[V(pcsf)$name]
  ##V(pcsf)$foldchange <- rho[V(pcsf)$name]
  
  return(pcsf)
}


#' @title Compute PCSF solution from pgx object for GENESET level
#'
#' @export
pgx.computePCSF_gset <- function(pgx, contrast,
                                 gmt.rho = 0.8, gset.filter = NULL,
                                 highcor = 0.9, 
                                 ntop = 250, ncomp = 3, beta = 1,
                                 rm.negedge = TRUE,
                                 dir = "both") {

  ## geneset parameters
  y <- pgx$model.parameters$exp.matrix[,contrast]
  ii <- which(!is.na(y) & y!=0)
  rho <- cor(Matrix::t(pgx$gsetX[,ii]), y[ii])[,1]  
  fc <- pgx.getMetaMatrix(pgx, level="geneset")$fc[,contrast]
  X <- pgx$gsetX

  ## filter on datatype or geneset collection
  if(!is.null(gset.filter)) {
    sel <- grep(gset.filter, names(rho))
    rho <- rho[sel]
    fc <- fc[sel]
    X <- X[sel,,drop=FALSE]
  }  

  ## take 2*ntop 
  ii <- head(names(sort(-abs(fc))),2*ntop)
  X  <- X[ii,,drop=FALSE]
  rho  <- rho[ii]
  fc  <- fc[ii]
  dim(X)
    
  ## align everything just to be sure
  labels <- rownames(X)

  ## Construct base PPI using GMT overlap
  aa <- intersect(names(fc), colnames(pgx$GMT))
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
  ## R <- qlcMatrix::corSparse(t(X))
  R <- cor(t(X), use = "pairwise")
  diag(R) <- 0
  hce <- which(abs(R) > highcor, arr.ind=TRUE)
  if(nrow(hce) > 0) {
    message("adding ",nrow(hce)," very-high-correlation edges rho>", highcor)
    ppi.hce <- data.frame(
      from = rownames(R)[hce[,1]],
      to = rownames(R)[hce[,2]],
      cost = 0.1
    )
    ppix <- rbind(ppix, ppi.hce)
  }
  message("PPI number of edges: n=", nrow(ppix))

  fx <- rho
  ##fx <- fc
  
  ## compute PCSF
  pcsf <- solvePCSF(
    X,
    fx = fc,
    ppi = ppix,
    labels = labels,
    ntop = ntop,
    ncomp = ncomp,
    beta = beta,
    rm.negedge = rm.negedge,
    dir = dir
  )
  
  ## Reinstate foldchange to actual foldchange
  V(pcsf)$foldchange <- fc[V(pcsf)$name]
  V(pcsf)$rho <- rho[V(pcsf)$name]
  
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
  dim(ee)

  suppressMessages(suppressWarnings(
    xppi <- PCSF::construct_interactome(ee)
  ))
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
                     highlightby = c("centrality", "prize")[1],
                     plotlib = c("visnet", "igraph")[1],
                     layout = "layout_with_kk", physics = TRUE,
                     node_cex = 30,
                     label_cex = 30, nlabel = -1, lab_exp = 3,
                     border_width = 1.5, cut.clusters = FALSE,
                     nlargest = -1, edge_width = 5) {
  ## set node size
  fx <- igraph::V(pcsf)$foldchange
  wt <- abs(fx / mean(abs(fx), na.rm = TRUE))**0.7
  node_cex1 <- node_cex * pmax(wt, 1)
  label_cex1 <- label_cex + 1e-8 * abs(fx)

  ## set colors as UP/DOWN
  up_down <- c("down", "up")[1 + 1 * (sign(fx) > 0)]
  dtype <- igraph::V(pcsf)$type
  igraph::V(pcsf)$group <- paste0(dtype, "_", up_down)
  bluered <- playdata::BLUERED(7)
  groups <- sort(setdiff(igraph::V(pcsf)$group, c(NA)))
  ngroups <- length(groups)
  extra_node_colors <- c(bluered[2], bluered[6])[1 + grepl("_up$", groups)]
  names(extra_node_colors) <- groups

  shapes <- head(rep(c(
    "dot", "triangle", "star", "diamond", "square",
    "triangleDown", "hexagon"
  ), 99), ngroups)
  levels <- names(sort(-table(dtype)))
  gtype <- gsub("_down|_up","",groups)
  extra_node_shapes <- shapes[factor(gtype, levels=levels)]
  names(extra_node_shapes) <- groups

  border_colors <- c("lightgrey",rep(c("yellow","magenta","cyan","green"),99))
  extra_node_borders <- border_colors[factor(gtype, levels=levels)]
  names(extra_node_borders) <- groups
  
  ## set label size
  if (highlightby == "centrality") {
    ewt <- igraph::E(pcsf)$weight
    wt <- 1.0 / (0.01 * mean(ewt, na.rm = TRUE) + ewt) ##
    bc <- igraph::page_rank(pcsf, weights = wt)$vector
    label_cex1 <- label_cex * (1 + 3 * (bc / max(bc, na.rm = TRUE))**lab_exp)
  }
  if (highlightby == "prize") {
    fx1 <- abs(fx)
    label_cex1 <- label_cex * (1 + 3 * (fx1 / max(fx1, na.rm = TRUE))**lab_exp)
  }
  if (highlightby %in% c("centrality.prize","prize.centrality")) {
    ewt <- igraph::E(pcsf)$weight
    wt <- 1.0 / (0.01 * mean(ewt, na.rm = TRUE) + ewt) ##
    bc <- igraph::page_rank(pcsf, weights = wt)$vector
    bc <- bc * abs(fx)
    label_cex1 <- label_cex * (1 + 3 * (bc / max(bc, na.rm = TRUE))**lab_exp)
  }

  if (is.null(igraph::V(pcsf)$label)) {
    igraph::V(pcsf)$label <- igraph::V(pcsf)$name
  }

  if(cut.clusters) {
    ewt <- 1/(1e-8 + igraph::E(pcsf)$weight)**2
    clust <- igraph::cluster_louvain(pcsf, weights=ewt)
    clust <- igraph::cluster_fast_greedy(pcsf, weights=ewt)
    table(clust$membership)
    ee <- igraph::E(pcsf)[igraph::crossing(clust, pcsf)]
    pcsf <- igraph::delete_edges(pcsf, ee)
  }

  ## take n-largest components
  if(nlargest > 0) {
    wc <- igraph::components(pcsf)$membership
    top.comp <- head(names(sort(table(wc),decreasing=TRUE)),nlargest)
    sel <- which(wc %in% top.comp)
    pcsf <- igraph::subgraph(pcsf, sel)
    label_cex1 <- label_cex1[sel]
    node_cex1 <- node_cex1[sel]    
  }   
  
  ## this equalizes label size per component
  wc <- igraph::components(pcsf)$membership
  comp <- unique(wc)
  for(k in comp) {
    ii <- which(wc == k)
    #f <- max(label_cex1) / max(label_cex1[ii])
    f <- mean(label_cex1) / mean(label_cex1[ii])    
    label_cex1[ii] <- label_cex1[ii] * f
  }

  ## this limits number of labels per component
  if (nlabel > 0) {
    for(k in comp) {
      ii <- which(wc == k)
      top.cex <- ii[head(order(-label_cex1[ii]), nlabel)]
      jj <- setdiff(ii, top.cex)
      igraph::V(pcsf)$label[jj] <- ""
    }
  }

  out <- NULL
  if (plotlib == "visnet") {
    library(igraph)
    class(pcsf) <- c("PCSF", "igraph")    
    out <- visplot.PCSF(
      pcsf,
      style = 1,
      node_size = node_cex1,
      node_label_cex = label_cex1,
      invert.weight = TRUE,
      edge_width = edge_width,
      border_width = border_width,
      Steiner_node_color = "lightblue",
      Terminal_node_color = "lightgreen",
      extra_node_colors = extra_node_colors,
      extra_node_shapes = extra_node_shapes,
      extra_node_borders = extra_node_borders,
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
    border_width = 1,
    width = 1800, height = 1800, invert.weight = FALSE,
    extra_node_colors = c(), extra_node_shapes = c(), extra_node_borders = c(),
    edge_color = "lightgrey", ...) {
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
  ewt <- igraph::E(pcsf)$weight
  wt <- 1.0 / (0.01 * mean(ewt, na.rm = TRUE) + ewt) ##
  cc <- igraph::page_rank(pcsf, weights = wt)$vector
  cc <- cc / (1e-8 + mean(cc, na.rm = TRUE)) ## normalize
  fc <- igraph::V(pcsf)$foldchange

  M <- data.frame(centrality = cc, logFC = fc)
  M <- round(M, digits = 2)
  match.col <- which.max(apply(pgx$genes, 2, function(a) sum(rownames(M) %in% a)))
  ii <- match(rownames(M), pgx$genes[,match.col])  
  aa <- pgx$genes[ii, c("feature", "symbol", "human_ortholog", "gene_title")]
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
