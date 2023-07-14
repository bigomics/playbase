##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##
##


#' Title
#'
#' @param pgx value
#' @param methods value
#' @param dims value
#' @param reduce.pca value
#' @param perplexity value
#' @param level value
#' @param rank.tf value
#' @param center.rows value
#' @param scale.rows value
#' @param X value
#' @param umap.pkg value
#'
#' @return
#' @export
#'
#' @examples
pgx.clusterGenes <- function(pgx, methods = c("pca", "tsne", "umap"), dims = c(2, 3),
                             reduce.pca = 50, perplexity = 30, level = "gene",
                             rank.tf = FALSE, center.rows = TRUE, scale.rows = FALSE,
                             X = NULL, umap.pkg = "uwot") {
  if (!is.null(X)) {
    message("using provided X matrix...")
  } else if (!is.null(pgx$X) && level == "gene") {
    message("using expression gene X matrix...")
    X <- pgx$X
  } else if (!is.null(pgx$gsetX) && level == "geneset") {
    message("using expression geneset X matrix...")
    X <- pgx$gsetX
  } else {
    message("WARNING:: could not find matrix X")
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
  }

  clust <- pgx.clusterBigMatrix(
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
    clust[ii] <- lapply(clust[ii], pos.compact)
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

  message("[pgx.clusterGenes] done!")
  pgx
}



#' Title
#'
#' @param pgx value
#' @param methods value
#' @param dims value
#' @param reduce.sd value
#' @param reduce.pca value
#' @param perplexity value
#' @param center.rows value
#' @param scale.rows value
#' @param X value
#' @param umap.pkg value
#' @param replace.orig value
#'
#' @return
#' @export
#'
#' @examples
pgx.clusterSamples2 <- function(pgx, methods = c("pca", "tsne", "umap"), dims = c(2, 3),
                                reduce.sd = 1000, reduce.pca = 50, perplexity = 30,
                                center.rows = TRUE, scale.rows = FALSE,
                                X = NULL, umap.pkg = "uwot", replace.orig = TRUE) {
  if (!is.null(X)) {
    message("using provided X matrix...")
  } else if (!is.null(pgx$X)) {
    message("using pgx$X matrix...")
    X <- pgx$X
  } else {
    message("using logCPM(pgx$counts)...")
    X <- logCPM(pgx$counts, total = NULL)
  }

  clust.pos <- pgx.clusterBigMatrix(
    X,
    methods = methods,
    dims = dims,
    perplexity = perplexity,
    center.features = center.rows,
    scale.features = scale.rows,
    reduce.sd = reduce.sd,
    reduce.pca = reduce.pca,
    find.clusters = FALSE,
    umap.pkg = umap.pkg
  )
  
  clust.index <- clust.pos$membership
  clust.pos$membership <- NULL

  if (replace.orig) {
    message("[pgx.clusterSamples2] update tsne2d/tsne3d and 'cluster' pheno...")
    pgx$samples$cluster <- clust.index
    pgx$tsne2d <- clust.pos[["tsne2d"]]
    pgx$tsne3d <- clust.pos[["tsne3d"]]
  } else {
    message("[pgx.clusterSamples2] skipping tsne2d/tsne3d update...")
  }

  pgx$cluster <- NULL
  pgx$cluster$pos <- clust.pos
  pgx$cluster$index <- clust.index

  message("[pgx.clusterSamples2] done!")
  pgx
}




#' Title
#'
#' @param pgx value
#' @param X value
#' @param skipifexists value
#' @param perplexity value
#' @param ntop value
#' @param npca value
#' @param prefix value
#' @param kclust value
#' @param dims value
#' @param find.clusters value
#' @param clust.detect value
#' @param row.center value
#' @param row.scale value
#' @param method value
#'
#' @return
#' @export
#'
#' @examples
pgx.clusterSamples <- function(pgx, X = NULL, skipifexists = FALSE, perplexity = 30,
                               ntop = 1000, npca = 50, prefix = "C", kclust = 1,
                               dims = c(2, 3), find.clusters = TRUE,
                               clust.detect = c("louvain", "hclust"),
                               row.center = TRUE, row.scale = FALSE,
                               method = c("tsne", "umap", "pca")) {
  clust.detect <- clust.detect[1]
  if (!is.null(X)) {
    message("using provided X matrix...")
  } else if (is.null(X) && !is.null(pgx$X)) {
    message("using pgx$X matrix...")
    X <- pgx$X
  } else if (is.null(X) && is.null(pgx$X)) {
    message("using pgx$counts matrix...")
    X <- logCPM(pgx$counts, total = NULL, prior = 1)
  } else {
    stop("[pgx.clusterSamples] FATAL ERROR")
  }

  res <- NULL
  res <- pgx.clusterMatrix(
    X,
    perplexity = perplexity, dims = dims,
    ntop = ntop, npca = npca, prefix = prefix,
    kclust = kclust, find.clusters = find.clusters,
    clust.detect = clust.detect,
    row.center = row.center, row.scale = row.scale,
    method = method
  )

  if (!is.null(res$pos2d)) pgx$tsne2d <- res$pos2d
  if (!is.null(res$pos3d)) pgx$tsne3d <- res$pos3d
  if (!is.null(res$idx)) {
    if (class(pgx$samples) == "data.frame") {
      pgx$samples$cluster <- as.character(res$idx)
    } else {
      ## matrix??
      if ("cluster" %in% colnames(pgx$samples)) {
        pgx$samples[, "cluster"] <- as.character(res$idx)
      } else {
        pgx$samples <- cbind(pgx$samples, cluster = as.character(res$idx))
      }
    }
  }

  return(pgx)
}


#' Title
#'
#' @param X value
#' @param method value
#' @param top.sd value
#' @param npca value
#'
#' @return
#' @export
#'
#' @examples
pgx.FindClusters <- function(X, method = c("kmeans", "hclust", "louvain", "meta"),
                             top.sd = 1000, npca = 50) {
  message("[FindClusters] called...")

  km.sizes <- c(2, 3, 4, 5, 7, 10, 15, 20, 25, 50, 100)
  km.sizes <- km.sizes[km.sizes < ncol(X)]
  km.sizes
  if (length(method) == 1 && method[1] == "meta") {
    method <- c("kmeans", "hclust", "louvain", "meta")
  }

  ## reduce dimensions
  X <- Matrix::head(X[order(apply(X, 1, sd)), ], top.sd)
  X <- t(scale(t(X))) ## scale features??
  if (nrow(X) > npca) {
    npca <- min(npca, dim(X) - 1)
    suppressMessages(suppressWarnings(
      out <- irlba::irlba(X, nv = npca)
    ))
    X <- t(out$v)
  }

  index <- list()

  ## perform K-means
  if ("kmeans" %in% method) {
    message("perform K-means...")
    km <- lapply(km.sizes, function(k) kmeans(t(X), k, iter.max = 10))
    km.idx <- do.call(cbind, lapply(km, function(r) r$cluster))
    colnames(km.idx) <- paste0("kmeans.", km.sizes)
    index[["kmeans"]] <- km.idx
  }

  ## perform hclust (on positions)
  if ("hclust" %in% method) {
    message("perform hclust...")
    hc <- fastcluster::hclust(dist(t(X)), method = "ward.D")
    hc.idx <- lapply(km.sizes, function(k) cutree(hc, k))
    hc.idx <- do.call(cbind, hc.idx)
    colnames(hc.idx) <- paste0("hclust.", km.sizes)
    index[["hclust"]] <- hc.idx
  }

  ## perform Louvain clustering
  if ("louvain" %in% method) {
    message("perform Louvain clustering...")
    gr <- scran::buildSNNGraph(X)
    gr.idx <- hclustGraph(gr, k = 4) ## iterative cluster until level3
    rownames(gr.idx) <- colnames(X)
    nc <- apply(gr.idx, 2, function(x) length(unique(x)))
    colnames(gr.idx) <- paste0("louvain.", nc)
    index[["louvain"]] <- gr.idx
  }

  ## find meta-index
  if ("meta" %in% method && length(index) > 1) {
    message("perform meta clustering...")
    K <- do.call(cbind, index)
    k.rows <- split(K, row(K))
    d1 <- outer(k.rows, k.rows, Vectorize(function(x, y) sum(x != y)))
    rownames(d1) <- colnames(d1) <- rownames(K)
    hc <- fastcluster::hclust(as.dist(d1))
    meta.idx <- do.call(cbind, lapply(km.sizes, function(k) cutree(hc, k)))
    colnames(meta.idx) <- paste0("meta.", km.sizes)
    rownames(meta.idx) <- rownames(K)
    index[["meta"]] <- meta.idx
  }

  ## sort cluster index from big to small clusters
  relevelBig2small <- function(idx) as.integer(factor(idx, levels = names(sort(-table(idx)))))
  for (i in 1:length(index)) {
    index[[i]] <- apply(index[[i]], 2, relevelBig2small)
  }

  return(index)
}


#' Title
#'
#' @param X value
#' @param methods value
#' @param dims value
#' @param perplexity value
#' @param reduce.sd value
#' @param reduce.pca value
#' @param center.features value
#' @param scale.features value
#' @param find.clusters value
#' @param svd.gamma value
#' @param umap.pkg value
#'
#' @return
#' @export
#'
#' @examples
pgx.clusterBigMatrix <- function(X, methods = c("pca", "tsne", "umap"), dims = c(2, 3),
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
      message("WARNING:: SD is constant. Skipping SD reduction...\n")
    } else {
      message("Reducing to ", reduce.sd, " max SD features...\n")
      X <- X[head(order(-sdx), reduce.sd), ]
    }
  }

  ## scale and augment if few samples
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
    reduce.pca
    message("Reducing to ", reduce.pca, " PCA dimenstions...\n")
    cnx <- colnames(X)
    suppressMessages(suppressWarnings(
      res.svd <- irlba::irlba(X, nv = reduce.pca)
    ))
    X <- t(res.svd$v) * res.svd$d**svd.gamma ## really weight with D??
    colnames(X) <- cnx
  }

  all.pos <- list()

  if ("pca" %in% methods && 2 %in% dims) {
    message("calculating PCA 2D/3D...\n")

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
    message("calculating t-SNE 2D...\n")

    perplexity <- pmax(min(ncol(X) / 4, perplexity), 2)
    perplexity
    pos <- Rtsne::Rtsne(t(X),
      dims = 2,
      ## pca = TRUE, partial_pca = TRUE,
      is_distance = FALSE, check_duplicates = FALSE,
      perplexity = perplexity, num_threads = 0
    )$Y
    rownames(pos) <- colnames(X)
    pos <- pos[1:dimx[2], ] ## if augmented
    colnames(pos) <- paste0("tSNE-", c("x", "y"))
    all.pos[["tsne2d"]] <- pos
  }

  if ("tsne" %in% methods && 3 %in% dims) {
    message("calculating t-SNE 3D...\n")

    perplexity <- pmax(min(dimx[2] / 4, perplexity), 2)
    perplexity
    pos <- Rtsne::Rtsne(t(X[, ]),
      dims = 3,
      ## pca = TRUE, partial_pca = TRUE,
      is_distance = FALSE, check_duplicates = FALSE,
      perplexity = perplexity, num_threads = 0
    )$Y
    rownames(pos) <- colnames(X)
    pos <- pos[1:dimx[2], ] ## if augmented
    colnames(pos) <- paste0("tSNE-", c("x", "y", "z"))
    all.pos[["tsne3d"]] <- pos
  }

  if ("umap" %in% methods && 2 %in% dims) {
    message("calculating UMAP 2D...\n")
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
    message("calculating UMAP 3D...\n")
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
    message("*** DEPRECATED *** please call seperately")
    message("calculating Louvain memberships (from reduced X)...")
    idx <- pgx.findLouvainClusters(t(X), level = 1, prefix = "c", small.zero = 0.01)
    all.pos$membership <- idx[1:dimx[2]]
  }

  return(all.pos)
}


#' Title
#'
#' @param X value
#' @param perplexity value
#' @param dims value
#' @param ntop value
#' @param npca value
#' @param prefix value
#' @param row.center value
#' @param row.scale value
#' @param find.clusters value
#' @param kclust value
#' @param clust.detect value
#' @param method value
#'
#' @return
#' @export
#'
#' @examples
pgx.clusterMatrix <- function(X, perplexity = 30, dims = c(2, 3),
                              ntop = 1000, npca = 50, prefix = "c",
                              row.center = TRUE, row.scale = FALSE,
                              find.clusters = TRUE, kclust = 1,
                              clust.detect = c("louvain", "hclust"),
                              method = c("tsne", "umap", "pca")) {
  method <- method[1]
  clust.detect <- clust.detect[1]
  X <- Matrix::head(X[order(-apply(X, 1, sd)), ], ntop)
  if (row.center) X <- X - rowMeans(X, na.rm = TRUE)
  if (row.scale) X <- (X / apply(X, 1, sd, na.rm = TRUE))

  ## some randomization is sometimes necessary if the data is 'too
  ## clean' and clusters become lines..


  ## ------------ find t-SNE clusters
  max.perplexity <- max(1, round((ncol(X) - 1) / 4))
  if (is.null(perplexity)) {
    perplexity <- max.perplexity
  }
  if (perplexity > max.perplexity) {
    message("[pgx.clusterMatrix] perplexity too large, decreasing perplexity to ", max.perplexity)
    perplexity <- max.perplexity
  }
  perplexity

  if (npca > 0) {
    npca <- min(npca, dim(X) - 1)
    message("performing tSNE on reduced PCA k=", npca)
    suppressMessages(suppressWarnings(
      svd <- irlba::irlba(X, nv = npca)
    ))
    sv <- svd$v %*% diag(svd$d[1:ncol(svd$v)])
    rownames(sv) <- colnames(X)
    X <- t(sv)
  }

  pos2 <- pos3 <- NULL
  if (method == "umap") {
    if (2 %in% dims) {
      pos2 <- uwot::umap(
        t(X),
        n_components = 2,
        metric = "euclidean",
        n_neighbors = perplexity,
        local_connectivity = ceiling(perplexity / 15),
        min_dist = 0.1
      )
      colnames(pos2) <- c("umap_1", "umap_2")
    }
    if (3 %in% dims) {
      pos3 <- uwot::umap(
        t(X),
        n_components = 3,
        metric = "euclidean",
        n_neighbors = perplexity,
        local_connectivity = ceiling(perplexity / 15),
        min_dist = 0.1
      )
      colnames(pos3) <- c("umap_1", "umap_2", "umap_3")
    }
  } else if (method == "tsne") {
    if (2 %in% dims) {
      pos2 <- Rtsne::Rtsne(t(X),
        dim = 2, perplexity = perplexity,
        ## pca = TRUE, partial_pca = TRUE,
        check_duplicates = FALSE, num_threads = 0
      )$Y
      colnames(pos2) <- c("tsne_1", "tsne_2")
    }
    if (3 %in% dims) {
      pos3 <- Rtsne::Rtsne(t(X),
        dim = 3, perplexity = perplexity,
        ## pca = TRUE, partial_pca = TRUE,
        check_duplicates = FALSE, num_threads = 0
      )$Y
      colnames(pos3) <- c("tsne_1", "tsne_2", "tsne_3")
    }
  } else if (method == "pca") {
    suppressMessages(suppressWarnings(
      svd <- irlba::irlba(X, nv = 3)
    ))
    if (2 %in% dims) {
      pos2 <- svd$v[, 1:2]
      colnames(pos2) <- c("pca_1", "pca_2")
    }
    if (3 %in% dims) {
      pos3 <- svd$v[, 1:3]
      colnames(pos3) <- c("pca_1", "pca_2", "pca_3")
    }
  }
  if (!is.null(pos2)) {
    rownames(pos2) <- colnames(X)
  }
  if (!is.null(pos3)) {
    rownames(pos3) <- colnames(X)
  }

  if (!is.null(pos2)) pos <- pos2
  if (!is.null(pos3)) pos <- pos3
  idx <- pgx.findLouvainClusters(pos, level = 1, prefix = "c", small.zero = 0.01)

  res <- list(pos2d = pos2, pos3d = pos3, idx = idx)
  return(res)
}


#' Title
#'
#' @param X value
#' @param prefix value
#' @param level value
#' @param gamma value
#' @param small.zero value
#'
#' @return
#' @export
#'
#' @examples
pgx.findLouvainClusters.SNN <- function(X, prefix = "c", level = 1, gamma = 1, small.zero = 0.01) {
  ## find clusters using graph clustering method
  message("perform Louvain clustering...")


  if (level == 1) {
    suppressMessages(suppressWarnings(
      gr <- scran::buildSNNGraph(t(X), d = 50)
    ))
  } else {
    ## finer clusters
    suppressMessages(suppressWarnings(
      gr <- scran::buildSNNGraph(t(X), d = 50, k = 2)
    ))
  }

  idx <- igraph::cluster_louvain(gr)$membership
  idx <- paste0(prefix, idx)

  if (!is.null(idx) && small.zero > 0) {
    ## ------------ zap small clusters to "0"
    sort(table(idx))
    min.size <- pmax(3, small.zero * length(idx))
    min.size
    small.clusters <- names(which(table(idx) < min.size))
    idx[which(idx %in% small.clusters)] <- "0"
    sort(table(idx))
  }

  ## rename levels with largest cluster first
  idx <- factor(idx, levels = names(sort(-table(idx))))
  levels(idx) <- paste0(prefix, 1:length(levels(idx)))
  idx <- as.character(idx)
  message("Found ", length(unique(idx)), " clusters...")
  return(idx)
}


#' Title
#'
#' @param X value
#' @param graph.method value
#' @param level value
#' @param prefix value
#' @param gamma value
#' @param small.zero value
#'
#' @return
#' @export
#'
#' @examples
pgx.findLouvainClusters <- function(X, graph.method = "dist", level = 1, prefix = "c",
                                    gamma = 1, small.zero = 0.01) {
  ## find clusters from t-SNE positions
  idx <- NULL
  message("Finding clusters using Louvain...\n")


  if (graph.method == "dist") {
    dist <- as.dist(dist(scale(X)))
    gr <- igraph::graph_from_adjacency_matrix(1.0 / dist**gamma, diag = FALSE, mode = "undirected")
  } else if (graph.method == "snn") {
    suppressMessages(suppressWarnings(gr <- scran::buildSNNGraph(t(X), d = 50)))
  } else {
    stop("FATAL: unknown graph method ", graph.method)
  }

  ## should we iteratively cluster (louvain)???
  hc <- hclustGraph(gr, k = level)
  idx <- hc[, min(level, ncol(hc))]

  if (!is.null(idx) && small.zero > 0) {
    ## ------------ zap small clusters to "0"
    sort(table(idx))
    min.size <- pmax(3, small.zero * length(idx))
    min.size
    small.clusters <- names(which(table(idx) < min.size))
    idx[which(idx %in% small.clusters)] <- "0"
    sort(table(idx))
  }

  ## rename levels with largest cluster first
  idx <- factor(idx, levels = names(sort(-table(idx))))
  levels(idx) <- paste0(prefix, 1:length(levels(idx)))
  idx <- as.character(idx)
  message("Found ", length(unique(idx)), " clusters...")
  return(idx)
}
