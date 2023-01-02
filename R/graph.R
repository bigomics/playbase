# hclustGraph
hclust_graph <- function(g, k = NULL, mc.cores = 2) {
  ## Hierarchical clustering of graph using iterative Louvain
  ## clustering on different levels. If k=NULL iterates until
  ## convergences.
  idx <- rep(1, length(igraph::V(g)))
  K <- c()
  maxiter <- 100
  if (!is.null(k)) maxiter <- k
  iter <- 1
  ok <- 1
  idx.len <- -1
  while (iter <= maxiter && ok) {
    old.len <- idx.len
    newidx0 <- newidx <- idx
    i <- idx[1]
    if (mc.cores > 1 && length(unique(idx)) > 1) {
      idx.list <- tapply(1:length(idx), idx, list)
      system.time(newidx0 <- parallel::mclapply(idx.list, function(ii) {
        subg <- igraph::induced_subgraph(g, ii)
        subi <- igraph::cluster_louvain(subg)$membership
        return(subi)
      }, mc.cores = mc.cores))
      newidx0 <- lapply(1:length(newidx0), function(i) paste0(i, "-", newidx0[[i]]))
      newidx0 <- as.vector(unlist(newidx0))
      newidx <- rep(NA, length(idx))
      newidx[as.vector(unlist(idx.list))] <- newidx0
    } else {
      for (i in unique(idx)) {
        ii <- which(idx == i)
        subg <- igraph::induced_subgraph(g, ii)
        subi <- igraph::cluster_louvain(subg)$membership
        newidx[ii] <- paste(i, subi, sep = "-")
      }
    }
    vv <- names(sort(table(newidx), decreasing = TRUE))
    idx <- as.integer(factor(newidx, levels = vv))
    K <- cbind(K, idx)
    idx.len <- length(table(idx))
    ok <- (idx.len > old.len)
    iter <- iter + 1
  }
  if (NCOL(K) == 1) K <- matrix(K, ncol = 1)
  rownames(K) <- igraph::V(g)$name
  if (!ok && is.null(k)) K <- K[, 1:(ncol(K) - 1), drop = FALSE]
  colnames(K) <- NULL
  return(K)
}

# pgx.computeCoreGOgraph
compute_core_go_graph <- function(ngs, fdr = 0.05) {
  require(igraph)

  ## test if there are GO terms
  mx <- ngs$gset.meta$meta[[1]]
  jj <- grep("^GO", rownames(mx))
  has.go <- (length(jj) > 10)
  if (!has.go) {
    return(NULL)
  }

  comparisons <- names(ngs$gset.meta$meta)
  subgraphs <- list()
  for (i in 1:length(comparisons)) {
    subgraphs[[i]] <- get_sig_go(
      ngs,
      comparison = comparisons[i],
      methods = NULL, ## should be actual selected methods!!
      fdr = fdr, nterms = 200, ntop = 20
    )
  }

  sub2 <- igraph::graph.union(subgraphs, byname = TRUE)
  A <- data.frame(igraph::vertex.attributes(sub2))
  rownames(A) <- A$name

  go_graph <- get_go_graph()
  go_graph <- igraph::induced.subgraph(go_graph, igraph::V(sub2)$name)
  A <- A[igraph::V(go_graph)$name, ]
  Q <- S <- V <- c()
  j1 <- grep("^score", colnames(A))
  j2 <- grep("^foldchange", colnames(A))
  j3 <- grep("^qvalue", colnames(A))
  for (i in 1:length(comparisons)) {
    S <- cbind(S, as.numeric(A[, j1[i]]))
    V <- cbind(V, as.numeric(A[, j2[i]]))
    Q <- cbind(Q, as.numeric(A[, j3[i]]))
  }
  rownames(Q) <- rownames(V) <- rownames(S) <- rownames(A)
  colnames(Q) <- colnames(V) <- colnames(S) <- comparisons

  ## can we match the GO terms with our gsetX values??
  go.sets <- grep("^GO.*\\(GO_", rownames(ngs$gsetX), value = TRUE)
  go.id <- gsub("GO_", "GO:", gsub(".*\\(|\\)", "", go.sets))
  matched.gset <- go.sets[match(igraph::V(go_graph)$name, go.id)]
  names(matched.gset) <- igraph::V(go_graph)$name

  ## compute FR layout
  layoutFR <- igraph::layout_with_fr(go_graph)
  layoutFR[, 1] <- 1.33 * layoutFR[, 1]
  rownames(layoutFR) <- igraph::V(go_graph)$name
  go_graph$layout <- layoutFR
  igraph::V(go_graph)$label <- igraph::V(go_graph)$Term
  res <- list(
    graph = go_graph, pathscore = S, foldchange = V,
    qvalue = Q, match = matched.gset
  )
  return(res)
}

# pgx.getSigGO
get_sig_go <- function(ngs, comparison, methods = NULL,
                       fdr = 0.20, nterms = 500, ntop = 100) {
  mx <- ngs$gset.meta$meta[[comparison]]
  jj <- grep("^GO", rownames(mx))
  length(jj)
  if (length(jj) == 0) {
    return(NULL)
  }
  mx <- mx[jj, , drop = FALSE]

  ## All methods????
  if (is.null(methods)) {
    methods <- colnames(unclass(mx$p))
  }
  methods

  ## recalculate meta values
  pv <- unclass(mx$p)[, methods, drop = FALSE]
  qv <- unclass(mx$q)[, methods, drop = FALSE]
  fc <- unclass(mx$fc)[, methods, drop = FALSE]
  pv[is.na(pv)] <- 0.999
  qv[is.na(qv)] <- 0.999
  fc[is.na(fc)] <- 0
  score <- fc * (-log10(qv))
  if (NCOL(pv) > 1) {
    ss.rank <- function(x) scale(sign(x) * rank(abs(x), na.last = "keep"), center = FALSE)
    fc <- rowMeans(scale(fc, center = FALSE), na.rm = TRUE)
    pv <- apply(pv, 1, max, na.rm = TRUE)
    qv <- apply(qv, 1, max, na.rm = TRUE)
    score <- rowMeans(apply(score, 2, ss.rank), na.rm = TRUE)
  }

  vinfo <- data.frame(geneset = rownames(mx), score = score, fc = fc, pv = pv, qv = qv)
  colnames(vinfo) <- c("geneset", "score", "fc", "pv", "qv") ## need
  rownames(vinfo) <- rownames(mx)
  remove(fc)

  terms <- AnnotationDbi::toTable(GO.db::GOTERM)[, 2:5]
  colnames(terms)[1] <- "go_id"
  terms <- terms[!duplicated(terms[, 1]), ]
  rownames(terms) <- terms[, 1]
  has.goid <- all(grepl(")$", rownames(vinfo)))
  if (has.goid) {
    ## rownames have GO ID at the end
    go_id <- gsub(".*\\(|\\)$", "", rownames(vinfo))
    go_id <- gsub("GO_|GO", "", go_id)
    go_id <- paste0("GO:", go_id)
    rownames(vinfo) <- go_id
    vinfo <- cbind(vinfo, terms[match(go_id, terms$go_id), , drop = FALSE])
  } else {
    ## rownames have no GO ID (downloaded from from MSigDB)
    vv <- sub("GO_|GOCC_|GOBP_|GOMF_", "", vinfo$geneset)
    idx <- match(vv, gsub("[ ]", "_", toupper(terms$Term)))
    jj <- which(!is.na(idx))
    vinfo <- cbind(vinfo[jj, , drop = FALSE], terms[idx[jj], , drop = FALSE])
    rownames(vinfo) <- vinfo$go_id
  }
  vinfo <- vinfo[which(!is.na(vinfo$go_id)), , drop = FALSE]

  ## Get full GO graph and assign node prizes
  go_graph <- get_go_graph()
  igraph::V(go_graph)$value <- rep(0, length(igraph::V(go_graph)))
  igraph::V(go_graph)[vinfo$go_id]$foldchange <- vinfo$fc
  igraph::V(go_graph)[vinfo$go_id]$qvalue <- vinfo$qv

  ## !!!!!!!!!!!!!!!!!!! THIS DEFINES THE SCORE !!!!!!!!!!!!!!!!!
  ## Value = "q-weighted fold-change"
  igraph::V(go_graph)[vinfo$go_id]$value <- vinfo$fc * (1 - vinfo$qv)**1

  get.vpath <- function(v1) {
    igraph::shortest_paths(go_graph, v1, "all")$vpath[[1]]
  }
  get.pathscore <- function(v1) {
    sp <- igraph::shortest_paths(go_graph, v1, "all")$vpath[[1]]
    sum((igraph::V(go_graph)[sp]$value))
  }

  sig.terms10 <- Matrix::head(rownames(vinfo)[order(vinfo$qv)], 10)
  sig.terms <- rownames(vinfo)[which(vinfo$qv <= fdr)]
  sig.terms <- unique(c(sig.terms, sig.terms10))
  sig.terms <- Matrix::head(sig.terms[order(vinfo[sig.terms, "qv"])], nterms) ## maximum number

  pathscore <- sapply(sig.terms, get.pathscore) ## SLOW!!!
  top.terms <- Matrix::head(sig.terms[order(-abs(pathscore))], ntop)

  ## total subgraph
  vv <- unique(unlist(sapply(top.terms, get.vpath)))
  vv <- igraph::V(go_graph)$name[vv]
  sub1 <- igraph::induced_subgraph(go_graph, vv)
  score1 <- (pathscore / max(abs(pathscore), na.rm = TRUE))[igraph::V(sub1)$name]
  igraph::V(sub1)$score <- pathscore[igraph::V(sub1)$name]
  igraph::V(sub1)$color <- gplots::bluered(32)[16 + round(15 * score1)]
  igraph::V(sub1)$label <- vinfo[igraph::V(sub1)$name, "Term"]

  return(sub1)
}

# getGOgraph
get_go_graph <- function() {
  terms <- AnnotationDbi::toTable(GO.db::GOTERM)[, 2:5]
  terms <- terms[!duplicated(terms[, 1]), ]
  rownames(terms) <- terms[, 1]

  BP <- AnnotationDbi::toTable(GO.db::GOBPPARENTS)
  MF <- AnnotationDbi::toTable(GO.db::GOMFPARENTS)
  CC <- AnnotationDbi::toTable(GO.db::GOCCPARENTS)
  bp.terms <- unique(c(BP[, 1], BP[, 2]))
  all.parents <- rbind(BP, MF, CC)
  go_graph <- igraph::graph_from_data_frame(all.parents, vertices = terms)
  return(go_graph)
}

# pgx.createOmicsGraph
create_omics_graph <- function(ngs, do.intersect = TRUE) {
  ## ======================================================================
  ## Create a graph object by merging nodes into
  ## clusters of genes/genesets.
  ##
  ##
  ## make bipartite igraph object
  ## ======================================================================


  ## ----------------------------------------------------------------------
  ## Read in gene/geneset graph structure
  ## ----------------------------------------------------------------------
  ## gr <- readRDS(file.path(FILES,"pgx-graph-geneXgset-XL.rds"))
  gr <- readRDS(file.path(FILES, "pgx-graph-geneXgset-XL-snn20.rds"))

  ## ----------------------------------------------------------------------
  ## Create large data matrix (includes all levels)
  ## ----------------------------------------------------------------------
  xx1 <- ngs$X
  gene <- toupper(ngs$genes[rownames(xx1), "gene_name"])
  rownames(xx1) <- paste0("{gene}", gene)
  xx2 <- ngs$gsetX
  rownames(xx2) <- paste0("{geneset}", rownames(ngs$gsetX))
  xx <- rbind(xx1, xx2)
  xx <- t(scale(t(xx))) ## scale??, then use innerproduct/cosine distance
  remove(xx1, xx2)

  ## ----------------------------------------------------------------------
  ## Prepare fold-change matrix
  ## ----------------------------------------------------------------------
  names(ngs$gx.meta$meta)
  FF <- sapply(ngs$gx.meta$meta, function(x) unclass(x$fc)[, "trend.limma"])
  FF <- FF / max(abs(FF), na.rm = TRUE)
  S <- sapply(ngs$gset.meta$meta, function(x) unclass(x$fc)[, "gsva"])
  S <- S / max(abs(S), na.rm = TRUE)
  rownames(S) <- paste0("{geneset}", rownames(S))
  fgene <- toupper(ngs$genes[rownames(FF), "gene_name"])
  rownames(FF) <- paste0("{gene}", fgene)

  kk <- intersect(colnames(FF), colnames(S))
  fc <- rbind(FF[, kk, drop = FALSE], S[, kk, drop = FALSE])
  remove(FF)
  remove(S)

  sel <- intersect(rownames(xx), rownames(fc))
  sel <- sort(intersect(sel, igraph::V(gr)$name))
  gr1 <- igraph::induced_subgraph(gr, sel)
  xx <- xx[igraph::V(gr1)$name, , drop = FALSE]
  fc <- fc[igraph::V(gr1)$name, , drop = FALSE]

  ## save the matched foldchange matrix
  gr1$foldchange <- fc
  gr1$scaled.data <- xx

  ## ----------------------------------------------------------------------
  ## should we recreate the SNNgraph in the intralayers????
  ## ----------------------------------------------------------------------
  ## this connect all points with at least 3 neighbours
  sel1 <- which(igraph::V(gr1)$level == "gene")
  pos1 <- gr$layout[igraph::V(gr1)[sel1]$name, ]
  r1 <- scran::buildSNNGraph(t(pos1), k = 3)
  igraph::V(r1)$name <- igraph::V(gr1)[sel1]$name

  sel2 <- which(igraph::V(gr1)$level == "geneset")
  pos2 <- gr$layout[igraph::V(gr1)[sel2]$name, ]
  r2 <- scran::buildSNNGraph(t(pos2), k = 3)
  igraph::V(r2)$name <- igraph::V(gr1)[sel2]$name
  new.gr <- igraph::graph.union(gr1, r1)
  new.gr <- igraph::graph.union(new.gr, r2)
  gr1 <- new.gr

  ## get rid of all weight attributes
  attr <- igraph::edge_attr_names(gr1)
  attr <- attr[grep("weight|rho", attr)]
  if (length(attr)) {
    for (i in 1:length(attr)) gr1 <- igraph::remove.edge.attribute(gr1, attr[i])
  }

  ## ----------------------------------------------------------------------
  ##  set correlation on edges
  ## ----------------------------------------------------------------------
  ee <- igraph::get.edges(gr1, igraph::E(gr1))
  ee.rho <- rep(NA, nrow(ee))
  bs <- 100000
  nblock <- ceiling(nrow(ee) / bs)
  for (k in 1:nblock) {
    i0 <- (k - 1) * bs + 1
    i1 <- min(i0 + bs - 1, nrow(ee))
    ii <- i0:i1
    ## fast method to compute rowwise-correlation (xx should be row-scaled)
    ee.rho[ii] <- rowMeans(xx[ee[ii, 1], ] * xx[ee[ii, 2], ])
  }
  igraph::E(gr1)$weight <- ee.rho ## replace weight with correlation

  ## ----------------------------------------------------------------------
  ## cluster graph
  ## ----------------------------------------------------------------------
  idx <- igraph::cluster_louvain(gr1, weights = abs(igraph::E(gr1)$weight))$membership
  igraph::V(gr1)$cluster <- idx
  return(gr1)
}

# pgx.computePathscores
compute_path_scores <- function(graph, strict.pos = TRUE) {
  ## add source/sink
  graph <- add_source_sink(graph)

  ## calculate weights for this particular contrast
  FF <- graph$foldchange
  P <- matrix(NA, nrow = length(igraph::V(graph)), ncol = ncol(FF))
  rownames(P) <- igraph::V(graph)$name
  colnames(P) <- colnames(FF)
  for (i in 1:ncol(FF)) {
    fc <- graph$foldchange[, i]
    ee <- igraph::get.edges(graph, igraph::E(graph))

    if (strict.pos) {
      f1 <- pmax(fc[ee[, 1]], 0) ## strictly positive
      f2 <- pmax(fc[ee[, 2]], 0)
      node.values <- sqrt(abs(f1 * f2)) ## strictly positive
      edge.rho <- pmax(igraph::E(graph)$weight, 0)
    } else {
      f1 <- fc[ee[, 1]] ## strictly positive
      f2 <- fc[ee[, 2]]
      node.values <- sqrt(abs(f1 * f2)) ## strictly positive
      edge.rho <- abs(igraph::E(graph)$weight)
    }
    score <- node.values * edge.rho ## always positive
    weights0 <- -log(pmax(score / max(score), 1e-8))
    length(weights0)
    summary(weights0)

    ## ----------------------------------------------------------
    ## Compute pathscore (solve all distances, 3-point SP)
    ## ----------------------------------------------------------
    dist.source <- igraph::distances(graph, v = "SOURCE", weights = weights0)
    dist.sink <- igraph::distances(graph, v = "SINK", weights = weights0)
    w1 <- rep(1, length(igraph::E(graph)))
    path.score <- exp(-(dist.source + dist.sink))
    names(path.score) <- igraph::V(graph)$name
    P[, i] <- path.score
  }

  return(P)
}

# pgx._addSourceSink
add_source_sink <- function(gr) {
  min.level <- min(gr$layout[igraph::V(gr)$name, 3])
  max.level <- max(gr$layout[igraph::V(gr)$name, 3])

  gr <- igraph::add_vertices(gr, 2, name = c("SOURCE", "SINK"))
  ss.layout <- rbind("SOURCE" = c(0, 0, -999), "SINK" = c(0, 0, 999))
  gr$layout <- rbind(gr$layout, ss.layout)
  nfc <- ncol(gr$foldchange)
  ss.foldchange <- rbind("SOURCE" = rep(1, nfc), "SINK" = rep(1, nfc))
  gr$foldchange <- rbind(gr$foldchange, ss.foldchange)
  gr$scaled.data <- rbind(gr$scaled.data,
    "SOURCE" = rep(NA, ncol(gr$scaled.data)),
    "SINK" = rep(NA, ncol(gr$scaled.data))
  )

  level <- gr$layout[igraph::V(gr)$name, 3]
  i1 <- which(level == min.level)
  i2 <- which(level == max.level)
  ee1 <- data.frame(from = "SOURCE", to = igraph::V(gr)$name[i1])
  ee2 <- data.frame(from = igraph::V(gr)$name[i2], to = "SINK")
  ee <- rbind(ee1, ee2)
  gr <- igraph::add_edges(gr, as.vector(t(ee)), weight = 1)
  return(gr)
}


# pgx.reduceOmicsGraph
reduce_omics_graph <- function(ngs) {
  ## ======================================================================
  ## Create a 'reduced' representation graph by merging nodes into
  ## clusters of genes/genesets.
  ##
  ##
  ## make bipartite igraph object
  ## ======================================================================

  ## get full omics graph
  gr <- ngs$omicsnet
  if (is.null(gr)) {
    stop("FATAL ERROR:: no omicsnet in ngs object. first run pgx.createOmicsGraph().")
  }

  ## ------------------------------------------------------------
  ## conform features
  ## ------------------------------------------------------------
  v1 <- which(igraph::V(gr)$level == "gene")
  v2 <- which(igraph::V(gr)$level == "geneset")
  g1 <- igraph::induced_subgraph(gr, v1)
  g2 <- igraph::induced_subgraph(gr, v2)
  g1 <- igraph::delete_edge_attr(g1, "weight")
  g2 <- igraph::delete_edge_attr(g2, "weight")
  h1 <- hclust_graph(g1)
  h2 <- hclust_graph(g2)
  apply(h1, 2, function(x) length(table(x)))
  apply(h2, 2, function(x) length(table(x)))
  hc1 <- paste0("{gene}cluster", h1[, ncol(h1)])
  hc2 <- paste0("{geneset}cluster", h2[, ncol(h2)])
  names(hc1) <- igraph::V(g1)$name
  names(hc2) <- igraph::V(g2)$name

  ## ------------------------------------------------------------
  ## create reduction transform matrix
  ## ------------------------------------------------------------
  idx0 <- c(h1[, 1], h2[, 1])[igraph::V(gr)$name]
  idx <- c(hc1, hc2)[igraph::V(gr)$name]
  R <- t(model.matrix(~ 0 + idx))
  R <- R / rowSums(R)
  colnames(R) <- igraph::V(gr)$name
  rownames(R) <- sub("^idx", "", rownames(R))

  ## reduce all variable: weights (correlation)
  rA <- (R %*% gr[, ]) %*% t(R) ## weights
  rA <- (rA + t(rA)) / 2

  ## reduced data matrix
  rX <- R %*% gr$scaled.data

  ## reduced fold-change
  rF <- R %*% gr$foldchange

  ## ------------------------------------------------------------
  ## keep group indices
  ## ------------------------------------------------------------
  grp.members <- tapply(names(idx), idx, list)
  grp.members0 <- grp.members
  grp.members <- lapply(grp.members, function(s) gsub("^\\{.*\\}", "", s))
  grp.label <- lapply(grp.members, function(s) paste(s, collapse = "\n"))

  ## ------------------------------------------------------------
  ## create reduced combined graph
  ## ------------------------------------------------------------
  gr1 <- igraph::graph_from_adjacency_matrix(as.matrix(rA),
    mode = "undirected",
    weighted = TRUE, diag = FALSE
  )

  ## rpos <- (R %*% gr$layout)
  ee <- igraph::get.edges(gr1, igraph::E(gr1))
  lev <- igraph::V(gr1)$level
  ee.type <- c("inter", "intra")[1 + 1 * (lev[ee[, 1]] == (lev[ee[, 2]]))]
  gr1 <- igraph::delete_edges(gr1, which(abs(igraph::E(gr1)$weight) < 0.01 & ee.type == "inter"))

  ## ------------------------------------------------------------
  ## compute new layout positions???
  ## ------------------------------------------------------------
  rpos <- lapply(grp.members0, function(m) colMeans(gr$layout[m, , drop = FALSE]))
  rpos <- do.call(rbind, rpos)
  gr1$layout <- scale(rpos)

  ## should we recreate the SNNgraph in the intralayers????
  ## this connect all points with at least 3 neighbours
  vtype <- gsub("\\}.*|^\\{", "", igraph::V(gr1)$name)
  sel1 <- which(vtype == "gene")
  pos1 <- gr1$layout[igraph::V(gr1)[sel1]$name, ]
  r1 <- scran::buildSNNGraph(t(pos1), k = 3)
  igraph::V(r1)$name <- rownames(pos1)

  sel2 <- which(vtype == "geneset")
  pos2 <- gr1$layout[igraph::V(gr1)[sel2]$name, ]
  r2 <- scran::buildSNNGraph(t(pos2), k = 3)
  igraph::V(r2)$name <- rownames(pos2)
  new.gr <- igraph::graph.union(r1, r2)
  new.gr <- igraph::remove.edge.attribute(new.gr, "weight_1")
  new.gr <- igraph::remove.edge.attribute(new.gr, "weight_2")
  new.gr <- igraph::graph.union(gr1, new.gr)
  jj <- which(is.na(igraph::E(new.gr)$weight))
  igraph::E(new.gr)$weight[jj] <- 0 ## actually we should recompute this....
  gr1 <- new.gr

  ## ------------------------------------------------------------
  ## add some graph data
  ## ------------------------------------------------------------
  igraph::V(gr1)$label <- grp.label
  igraph::V(gr1)$cluster <- tapply(idx0, idx, median) ## level 1 cluster index
  igraph::V(gr1)$level <- gsub("\\}.*|^\\{", "", igraph::V(gr1)$name)

  gr1$foldchange <- rF
  gr1$members <- grp.members
  gr1$scaled.data <- rX

  return(gr1)
}
