
#pos.compact
pos_compact <- function(pos, d = 0.01) {
  ## make positions more dense removing white space
  for (i in 1:ncol(pos)) {
    x <- pos[, i]
    dr <- d * diff(range(x))
    names(x) <- 1:nrow(pos)
    ii <- order(x)
    x1 <- cumsum(c(x[ii[1]], pmin(diff(x[ii]), dr)))
    pos[, i] <- x1[order(as.integer(names(x1)))]
  }
  pos
}

#pgx.clusterGenes
cluster_genes <- function(pgx, methods = c("pca", "tsne", "umap"), dims = c(2, 3),
                             reduce.pca = 50, perplexity = 30, level = "gene",
                             rank.tf = FALSE, center.rows = TRUE, scale.rows = FALSE,
                             X = NULL, umap.pkg = "uwot") {
  if (!is.null(X)) {
    invisible()
  } else if (!is.null(pgx$X) && level == "gene") {
    X <- pgx$X
  } else if (!is.null(pgx$gsetX) && level == "geneset") {
    X <- pgx$gsetX
  } else {
    return(pgx)
  }

  if (center.rows) {
    X <- X - rowMeans(X)
  }
  if (scale.rows) {
    X <- X / (1e-6 + apply(X, 1, sd))
  }
  if (rank.tf) {
    X <- scale(apply(X, 2, rank))
  } ## works nicely

  clust <- cluster_big_matrix(
    t(X),
    methods = methods,
    dims = dims,
    perplexity = perplexity,
    center.features = FALSE,
    scale.features = FALSE,
    reduce.sd = -1,
    reduce.pca = reduce.pca,
    find.clusters = FALSE,
    umap.pkg = umap.pkg
  )

  clust.index <- NULL
  clust.index <- clust$membership
  clust$membership <- NULL

  ## remove empty space in tSNE/UMAP
  ii <- grep("tsne|umap", names(clust))
  if (length(ii) > 0) {
    clust[ii] <- lapply(clust[ii], pos_compact) ## make more compact
  }

  ## put in slot 'gene cluster'
  if (level == "gene") {
    pgx$cluster.genes <- NULL
    pgx$cluster.genes$pos <- clust
    pgx$cluster.genes$index <- clust.index
  }
  if (level == "geneset") {
    pgx$cluster.gsets <- NULL
    pgx$cluster.gsets$pos <- clust
    pgx$cluster.gsets$index <- clust.index
  }

  pgx
}

#pgx.findLouvainClusters
find_louvain_clusters <- function(X, graph.method = "dist", level = 1, prefix = "c",
                                    gamma = 1, small.zero = 0.01) {
  ## find clusters from t-SNE positions
  idx <- NULL

  if (graph.method == "dist") {
    dist <- as.dist(dist(scale(X))) ## use 3D distance??
    gr <- igraph::graph_from_adjacency_matrix(1.0 / dist**gamma, diag = FALSE, mode = "undirected")
  } else if (graph.method == "snn") {
    suppressMessages(suppressWarnings(gr <- scran::buildSNNGraph(t(X), d = 50)))
  } else {
    stop("FATAL: unknown graph method ", graph.method)
  }

  ## should we iteratively cluster (louvain)???
  hc <- hclust_graph(gr, k = level) ##
  idx <- hc[, min(level, ncol(hc))]

  if (!is.null(idx) && small.zero > 0) {
    ## ------------ zap small clusters to "0"
    min.size <- pmax(3, small.zero * length(idx))
    small.clusters <- names(which(table(idx) < min.size))
    idx[which(idx %in% small.clusters)] <- "0"
  }

  ## rename levels with largest cluster first
  idx <- factor(idx, levels = names(sort(-table(idx))))
  levels(idx) <- paste0(prefix, 1:length(levels(idx)))
  idx <- as.character(idx)
  return(idx)
}

#pgx.clusterBigMatrix
cluster_big_matrix <- function(X, methods = c("pca", "tsne", "umap"), dims = c(2, 3),
                                 perplexity = 30, reduce.sd = 1000, reduce.pca = 50,
                                 center.features = TRUE, scale.features = FALSE,
                                 find.clusters = FALSE, svd.gamma = 1, umap.pkg = "uwot") {
  methods <- intersect(methods, c("pca", "tsne", "umap"))
  if (length(methods) == 0) methods <- "pca"

  ## Reduce dimensions by SD
  dimx <- dim(X) ## original dimensions
  namesx <- colnames(X)
  if (reduce.sd > 0 && nrow(X) > reduce.sd) {
    sdx <- apply(X, 1, sd, na.rm = TRUE)
    is.constant <- all(abs(sdx - mean(sdx, na.rm = TRUE)) < 1e-8)
    if (is.constant) {
      warning("SD is constant. Skipping SD reduction...\n")
    } else {
      X <- X[head(order(-sdx), reduce.sd), ]
    }
  }

  ## scale and augment if few samples
  ## X <- scale(X) ## columm normalization??
  if (center.features) {
    X <- X - rowMeans(X, na.rm = TRUE) ## do??
  }
  if (scale.features) {
    X <- X / apply(X, 1, sd, na.rm = TRUE)
  }

  ## impute on row median
  if (any(is.na(X))) {
    X <- imputeMedian(X)
  }

  if (ncol(X) <= 6) X <- cbind(X, X, X, X, X, X)
  if (nrow(X) <= 3) X <- rbind(X, X, X, X)

  ## add small variation...
  X <- X + 1e-3 * matrix(rnorm(length(X)), nrow(X), ncol(X))

  ## Further pre-reduce dimensions using SVD
  res.svd <- NULL
  if (reduce.pca > 0) {
    reduce.pca <- max(3, min(c(reduce.pca, dim(X) - 1)))
    cnx <- colnames(X)
    suppressMessages(suppressWarnings(
      res.svd <- irlba::irlba(X, nv = reduce.pca)
    ))
    X <- t(res.svd$v) * res.svd$d**svd.gamma ## really weight with D??
    colnames(X) <- cnx
  }

  all.pos <- list()

  if ("pca" %in% methods && 2 %in% dims) {
    if (is.null(res.svd)) {
      suppressMessages(suppressWarnings(
        res.svd <- irlba::irlba(X, nv = 3)
      ))
    }
    pos <- res.svd$v[, 1:2]
    rownames(pos) <- colnames(X)
    pos <- pos[1:dimx[2], ] ## if augmented
    colnames(pos) <- paste0("PC-", c("x", "y"))
    all.pos[["pca2d"]] <- pos
  }

  if ("pca" %in% methods && 3 %in% dims) {
    if (is.null(res.svd)) {
      suppressMessages(suppressWarnings(
        res.svd <- irlba::irlba(X, nv = 3)
      ))
    }
    pos <- res.svd$v[, 1:3]
    rownames(pos) <- colnames(X)
    pos <- pos[1:dimx[2], ] ## if augmented
    colnames(pos) <- paste0("PC-", c("x", "y", "z"))
    all.pos[["pca3d"]] <- pos
  }

  if ("tsne" %in% methods && 2 %in% dims) {
    perplexity <- pmax(min(ncol(X) / 4, perplexity), 2)
    pos <- Rtsne::Rtsne(t(X),
      dims = 2,
      is_distance = FALSE, check_duplicates = FALSE,
      perplexity = perplexity, num_threads = 0
    )$Y
    rownames(pos) <- colnames(X)
    pos <- pos[1:dimx[2], ] ## if augmented
    colnames(pos) <- paste0("tSNE-", c("x", "y"))
    all.pos[["tsne2d"]] <- pos
  }

  if ("tsne" %in% methods && 3 %in% dims) {
    perplexity <- pmax(min(dimx[2] / 4, perplexity), 2)
    pos <- Rtsne::Rtsne(t(X[, ]),
      dims = 3,
      is_distance = FALSE, check_duplicates = FALSE,
      perplexity = perplexity, num_threads = 0
    )$Y
    rownames(pos) <- colnames(X)
    pos <- pos[1:dimx[2], ] ## if augmented
    colnames(pos) <- paste0("tSNE-", c("x", "y", "z"))
    all.pos[["tsne3d"]] <- pos
  }

  if ("umap" %in% methods && 2 %in% dims) {
    if (umap.pkg == "uwot") {
      nb <- ceiling(pmax(min(dimx[2] / 4, perplexity), 2))
      pos <- uwot::umap(t(X[, ]),
        n_components = 2,
        n_neighbors = nb,
        local_connectivity = ceiling(nb / 15),
        min_dist = 0.1
      )
    } else {
      custom.config <- umap.defaults
      custom.config$n_components <- 2
      custom.config$n_neighbors <- pmax(min(dimx[2] / 4, perplexity), 2)
      pos <- umap::umap(t(X[, ]), custom.config)$layout
    }
    rownames(pos) <- colnames(X)
    pos <- pos[1:dimx[2], ] ## if augmented
    colnames(pos) <- paste0("UMAP-", c("x", "y"))
    all.pos[["umap2d"]] <- pos
  }

  if ("umap" %in% methods && 3 %in% dims) {
    if (umap.pkg == "uwot") {
      nb <- ceiling(pmax(min(dimx[2] / 4, perplexity), 2))
      pos <- uwot::umap(t(X[, ]),
        n_components = 3,
        n_neighbors = nb,
        local_connectivity = ceiling(nb / 15),
        min_dist = 0.1
      )
    } else {
      custom.config <- umap.defaults
      custom.config$n_components <- 3
      custom.config$n_neighbors <- pmax(min(dimx[2] / 4, perplexity), 2)
      pos <- umap::umap(t(X[, ]), custom.config)$layout
    }
    rownames(pos) <- colnames(X)
    pos <- pos[1:dimx[2], ] ## if augmented
    colnames(pos) <- paste0("UMAP-", c("x", "y", "z"))
    all.pos[["umap3d"]] <- pos
  }

  all.pos$membership <- NULL
  if (find.clusters) {
    stop("*** DEPRECATED *** please call find.clusters seperately")
  }

  return(all.pos)
}
