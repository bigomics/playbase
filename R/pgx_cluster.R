pos.compact <- function(pos, d=0.01) {
  ## make positions more dense removing white space
  for(i in 1:ncol(pos)) {
    x=pos[,i]
    dr = d*diff(range(x))
    names(x)=1:nrow(pos)
    ii <- order(x)
    x1 = cumsum(c(x[ii[1]],pmin(diff(x[ii]),dr)))
    pos[,i] = x1[order(as.integer(names(x1)))]
  }
  pos
}

pgx.clusterGenes <- function(pgx, methods=c("pca","tsne","umap"), dims=c(2,3),
                             reduce.pca=50, perplexity=30, level='gene',
                             rank.tf=FALSE, center.rows=TRUE, scale.rows=FALSE,
                             X=NULL, umap.pkg="uwot" )
{
  if(!is.null(X)) {
    message("using provided X matrix...")
  } else if(!is.null(pgx$X) && level=='gene') {
    message("using expression gene X matrix...")
    X <- pgx$X
  } else if(!is.null(pgx$gsetX) && level=='geneset') {
    message("using expression geneset X matrix...")
    X <- pgx$gsetX
  } else {
    message("WARNING:: could not find matrix X")
    return(pgx)
  }
  dim(X)
  ## X <- limma::normalizeQuantiles(X)  ##??
  if(center.rows)
    X <- X  - rowMeans(X)
  if(scale.rows)
    X <- X / (1e-6+apply(X,1,sd))
  if(rank.tf)
    X <- scale(apply(X,2,rank))  ## works nicely

  clust <- pgx.clusterBigMatrix(
    t(X), methods = methods,
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
  names(clust)
  ##clust.index <- paste0("c",clust.pos$membership)
  clust.index <- clust$membership
  clust$membership <- NULL

  if(1) {
    ## remove empty space in tSNE/UMAP
    ii <- grep("tsne|umap",names(clust))
    if(length(ii)>0) {
      clust[ii] <- lapply(clust[ii], pos.compact)  ## make more compact
    }
  }

  ## put in slot 'gene cluster'
  if(level=='gene') {
    pgx$cluster.genes <- NULL
    pgx$cluster.genes$pos  <- clust
    pgx$cluster.genes$index <- clust.index
  }
  if(level=='geneset') {
    pgx$cluster.gsets <- NULL
    pgx$cluster.gsets$pos  <- clust
    pgx$cluster.gsets$index <- clust.index
  }

  message("[pgx.clusterGenes] done!")
  pgx
}

pgx.findLouvainClusters <- function(X, graph.method='dist', level=1, prefix='c',
                                    gamma=1, small.zero=0.01)
{

  ## find clusters from t-SNE positions
  idx = NULL
  message("Finding clusters using Louvain...\n")


  if(graph.method=='dist') {
    dist = as.dist(dist(scale(X))) ## use 3D distance??
    ##dist = dist + 0.1*mean(dist)
    gr = igraph::graph_from_adjacency_matrix(1.0/dist**gamma, diag=FALSE, mode="undirected")
  } else if(graph.method=='snn') {
    suppressMessages( suppressWarnings( gr <- scran::buildSNNGraph(t(X), d=50) ))
  } else {
    stop("FATAL: unknown graph method ",graph.method)
  }

  ## should we iteratively cluster (louvain)???
  hc <- hclustGraph(gr,k=level)  ##
  dim(hc)
  idx <- hc[,min(level,ncol(hc))]
  ##idx <- igraph::cluster_louvain(gr)$membership
  table(idx)

  if(!is.null(idx) && small.zero>0) {
    ## ------------ zap small clusters to "0"
    sort(table(idx))
    min.size <- pmax(3, small.zero*length(idx))
    min.size
    small.clusters <- names(which(table(idx) < min.size))
    idx[which(idx %in% small.clusters)] <- "0"
    sort(table(idx))
  }

  ## rename levels with largest cluster first
  idx <- factor(idx, levels=names(sort(-table(idx))))
  levels(idx) <- paste0(prefix,1:length(levels(idx)))
  table(idx)
  idx <- as.character(idx)
  message("Found ",length(unique(idx))," clusters...")
  return(idx)
}

pgx.clusterBigMatrix <- function(X, methods=c("pca","tsne","umap"), dims=c(2,3),
                                 perplexity=30, reduce.sd = 1000, reduce.pca = 50,
                                 center.features=TRUE, scale.features=FALSE,
                                 find.clusters=FALSE, svd.gamma=1, umap.pkg="uwot")
{




  methods <- intersect(methods, c("pca","tsne","umap"))
  if(length(methods)==0) methods <- "pca"

  ## Reduce dimensions by SD
  dimx = dim(X)  ## original dimensions
  namesx = colnames(X)
  if(reduce.sd>0 && nrow(X)>reduce.sd) {
    sdx <- apply(X,1,sd,na.rm=TRUE)
    is.constant <- all( abs(sdx-mean(sdx,na.rm=TRUE)) < 1e-8 )
    if(is.constant) {
      message("WARNING:: SD is constant. Skipping SD reduction...\n")
    } else {
      message("Reducing to ",reduce.sd," max SD features...\n")
      X <- X[head(order(-sdx),reduce.sd),]
    }
  }

  ## scale and augment if few samples
  ## X <- scale(X) ## columm normalization??
  if(center.features) {
    X <- X - rowMeans(X,na.rm=TRUE) ## do??
  }
  if(scale.features) {
    X <- X / apply(X,1,sd,na.rm=TRUE)
  }

  ## impute on row median
  if(any(is.na(X))) {
    X <- imputeMedian(X)
  }

  if(ncol(X)<=6) X <- cbind(X,X,X,X,X,X)
  dim(X)
  #dbg("[pgx.clusterBigMatrix] 2a: dim(X)=",dim(X),"\n")

  if(nrow(X)<=3) X <- rbind(X,X,X,X)
  dim(X)
  #dbg("[pgx.clusterBigMatrix] 2b: dim(X)=",dim(X),"\n")

  ## add small variation...
  X <- X + 1e-3*matrix(rnorm(length(X)),nrow(X),ncol(X))

  ## Further pre-reduce dimensions using SVD
  res.svd <- NULL
  if(reduce.pca>0) {
    reduce.pca <- max(3,min(c(reduce.pca,dim(X)-1)))
    reduce.pca
    message("Reducing to ",reduce.pca," PCA dimenstions...\n")
    cnx = colnames(X)
    suppressMessages( suppressWarnings(
      res.svd <- irlba::irlba(X, nv=reduce.pca)
    ))
    X <- t(res.svd$v) * res.svd$d**svd.gamma  ## really weight with D??
    colnames(X) <- cnx
  }
  dim(X)

  #dbg("[pgx.clusterBigMatrix] 3: dim(X)=",dim(X),"\n")

  all.pos <- list()

  if("pca" %in% methods && 2 %in% dims) {
    message("calculating PCA 2D/3D...\n")

    if(is.null(res.svd)) {
      suppressMessages( suppressWarnings(
        res.svd <- irlba::irlba(X, nv=3)
      ))
    }
    pos <- res.svd$v[,1:2]
    rownames(pos) <- colnames(X)
    pos <- pos[1:dimx[2],]  ## if augmented
    colnames(pos) <- paste0("PC-",c("x","y"))
    all.pos[["pca2d"]] <- pos
  }

  if("pca" %in% methods && 3 %in% dims) {
    if(is.null(res.svd)) {
      suppressMessages( suppressWarnings(
        res.svd <- irlba::irlba(X, nv=3)
      ))
    }
    pos <- res.svd$v[,1:3]
    rownames(pos) <- colnames(X)
    pos <- pos[1:dimx[2],]  ## if augmented
    colnames(pos) <- paste0("PC-",c("x","y","z"))
    all.pos[["pca3d"]] <- pos
  }

  if("tsne" %in% methods && 2 %in% dims) {
    message("calculating t-SNE 2D...\n")

    perplexity <- pmax(min(ncol(X)/4,perplexity),2)
    perplexity
    pos <- Rtsne::Rtsne( t(X), dims = 2,
                         ## pca = TRUE, partial_pca = TRUE,
                         is_distance = FALSE, check_duplicates=FALSE,
                         perplexity = perplexity, num_threads=0)$Y
    rownames(pos) <- colnames(X)
    pos <- pos[1:dimx[2],]  ## if augmented
    colnames(pos) <- paste0("tSNE-",c("x","y"))
    all.pos[["tsne2d"]] <- pos
  }

  if("tsne" %in% methods && 3 %in% dims) {
    message("calculating t-SNE 3D...\n")

    perplexity <- pmax(min(dimx[2]/4,perplexity),2)
    perplexity
    pos <- Rtsne::Rtsne( t(X[,]), dims = 3,
                         ## pca = TRUE, partial_pca = TRUE,
                         is_distance = FALSE, check_duplicates=FALSE,
                         perplexity = perplexity, num_threads=0)$Y
    rownames(pos) <- colnames(X)
    pos <- pos[1:dimx[2],]  ## if augmented
    colnames(pos) <- paste0("tSNE-",c("x","y","z"))
    all.pos[["tsne3d"]] <- pos
  }

  if("umap" %in% methods && 2 %in% dims) {
    message("calculating UMAP 2D...\n")
    if(umap.pkg=="uwot") {

      nb = ceiling(pmax(min(dimx[2]/4,perplexity),2))
      pos <- uwot::umap( t(X[,]),
                         n_components = 2,
                         n_neighbors = nb,
                         local_connectivity = ceiling(nb/15),
                         min_dist = 0.1
      )
    } else {

      custom.config <- umap.defaults
      custom.config$n_components = 2
      custom.config$n_neighbors = pmax(min(dimx[2]/4,perplexity),2)
      pos <- umap::umap( t(X[,]), custom.config )$layout
    }
    rownames(pos) <- colnames(X)
    pos <- pos[1:dimx[2],]  ## if augmented
    colnames(pos) <- paste0("UMAP-",c("x","y"))
    all.pos[["umap2d"]] <- pos
  }

  if("umap" %in% methods && 3 %in% dims) {
    message("calculating UMAP 3D...\n")
    if(umap.pkg=="uwot") {

      nb = ceiling(pmax(min(dimx[2]/4,perplexity),2))
      pos <- uwot::umap( t(X[,]),
                         n_components = 3,
                         n_neighbors = nb,
                         local_connectivity = ceiling(nb/15),
                         min_dist = 0.1
      )

    } else {

      custom.config <- umap.defaults
      custom.config$n_components = 3
      custom.config$n_neighbors = pmax(min(dimx[2]/4,perplexity),2)
      pos <- umap::umap( t(X[,]), custom.config )$layout
    }
    dim(pos)
    rownames(pos) <- colnames(X)
    pos <- pos[1:dimx[2],]  ## if augmented
    colnames(pos) <- paste0("UMAP-",c("x","y","z"))
    all.pos[["umap3d"]] <- pos
  }

  length(all.pos)
  ##all.pos <- lapply(all.pos, pos.compact)  ## make more compact

  all.pos$membership <- NULL
  if(find.clusters) {
    message("*** DEPRECATED *** please call seperately")
    message("calculating Louvain memberships (from reduced X)...")
    ##X = X - rowMeans(X)
    idx <- pgx.findLouvainClusters(t(X), level=1, prefix='c', small.zero=0.01)
    table(idx)
    all.pos$membership <- idx[1:dimx[2]]
  }

  return(all.pos)
}