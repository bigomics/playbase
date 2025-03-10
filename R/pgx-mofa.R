#' Main function to call MOFA analysis on pgx object.
#'
#' @export
pgx.compute_mofa <- function(pgx, kernel = "MOFA", numfactors = 8,
                             ntop = 2000, gset.ntop = 2000,
                             add_gsets = FALSE,
                             factorizations = TRUE) {
  if (0) {
    numfactors <- 8
    ntop <- 2000
    gset.ntop <- 2000
    add_gsets <- FALSE
  }

  has.prefix <- (mean(grepl(":", rownames(pgx$X))) > 0.8)
  is.multiomics <- (pgx$datatype == "multi-omics" && has.prefix)

  xdata <- NULL
  if (is.multiomics) {
    message("[compute.multi_omics] splitting by omics-type ")
    xdata <- mofa.split_data(pgx$X)
  } else {
    message("[compute.multi_omics] splitting by gene role")
    xdata <- mofa.splitByGeneRole(pgx$X)
  }
  names(xdata)

  numfactors <- min(numfactors, min(dim(xdata[[1]])))
  numfactors

  ## add geneset layers if asked
  if (add_gsets) {
    message("[compute.multi_omics] adding geneset matrix")
    xdata <- c(xdata, list(gset = pgx$gsetX))
  }
  lapply(xdata, dim)

  ## cleanup samples
  samples <- pgx$samples
  samples$cluster <- NULL
  samples$sample <- NULL
  samples <- samples[, grep("^[.].*|^cluster", colnames(samples), invert = TRUE), drop = FALSE]
  discretized_samples <- pgx.discretizePhenotypeMatrix(samples)

  for (i in 1:length(xdata)) {
    d <- rownames(xdata[[i]])
    ## rownames(xdata[[i]]) <- stringi::stri_trans_general(d, "latin-ascii")
    rownames(xdata[[i]]) <- iconv(d, to = "ascii//TRANSLIT")
  }

  ## MOFA computation
  message("computing MOFA ...")
  dbg("[pgx.compute_mofa] numfactors =", numfactors)

  ## samples=discretized_samples;contrasts=pgx$contrasts;annot=pgx$genes;GMT=pgx$GMT;kernel="mofa";ntop=2000;numfactors=8
  mofa <- mofa.compute(
    xdata,
    samples = discretized_samples,
    contrasts = pgx$contrasts,
    annot = pgx$genes,
    pheno = NULL,
    GMT = pgx$GMT,
    kernel = kernel,
    scale_views = TRUE,
    ntop = ntop,
    gset.ntop = gset.ntop,
    max_iter = 200,
    numfactors = numfactors
  )


  if (factorizations) {
    all.kernels <- c(
      "mofa", "pca", "nmf", "nmf2", "mcia", "diablo", ## "wgcna",
      "rgcca", "rgcca.rgcca", "rgcca.rgccda", "rgcca.mcoa"
    )
    # all.kernels <- c("mofa","pca","nmf","nmf2","mcia","wgcna","diablo","rgcca"
    mofa$factorizations <- mofa.compute_factorizations(
      xdata,
      discretized_samples,
      pgx$contrasts,
      numfactors = numfactors,
      kernels = all.kernels
    )
  }

  ## LASAGNA computation
  message("computing LASAGNA model...")
  las.data <- list(
    X = xdata,
    samples = pgx$samples,
    contrasts = pgx$contrasts
  )
  lasagna <- lasagna.create_model(
    las.data,
    pheno = "contrasts", ntop = ntop, nc = 20,
    annot = pgx$genes, use.graphite = FALSE
  )

  ## pre-compute cluster positions
  xx <- mofa.split_data(lasagna$X)
  xx <- xx[names(xdata)] ## remove SOURCE/SINK

  message("computing cluster positions (samples)...")
  mofa$posx <- mofa.compute_clusters(xx, along = "samples")
  message("computing cluster positions (features)...")
  mofa$posf <- mofa.compute_clusters(xx, along = "features")
  mofa$lasagna <- lasagna
  mofa$snf <- playbase::snf.cluster(mofa$xx, pheno = NULL, plot = FALSE)

  ## compute multi-gsea enrichment
  F <- pgx.getMetaMatrix(pgx)$fc
  mofa$mgsea <- mgsea.compute_enrichment(
    F,
    G = pgx$GMT,
    annot = pgx$genes,
    filter = NULL,
    ntop = gset.ntop
  )
  return(mofa)
}



#' @export
mofa.compute <- function(xdata,
                         samples,
                         contrasts,
                         pheno = NULL,
                         kernel = "mofa",
                         factorname = NULL,
                         GMT = NULL,
                         annot = NULL,
                         scale_views = TRUE,
                         ntop = 2000,
                         gset.ntop = 2000,
                         max_iter = 1000,
                         gpu_mode = FALSE,
                         compute.enrichment = TRUE,
                         compute.graphs = TRUE,
                         numfactors = 8,
                         numfeatures = 100) {
  if (!is.null(ntop) && ntop > 0) {
    info("[mofa.compute] reducing data blocks: ntop = ", ntop)
    ## fast prioritization of features using SD or rho (correlation
    ## with contrasts).
    i <- 1
    num.contrasts <- sign(makeContrastsFromLabelMatrix(contrasts))
    for (i in 1:length(xdata)) {
      d <- xdata[[i]]
      rho <- cor(t(d), num.contrasts, use = "pairwise")
      rho[is.na(rho)] <- 0
      rho2 <- rowMeans(rho**2, na.rm = TRUE) * rowMeans(!is.na(d))
      d <- head(d[order(-rho2), , drop = FALSE], ntop)
      xdata[[i]] <- d
    }
  }
  lapply(xdata, dim)

  ## impute missing?
  nna <- sapply(xdata, function(x) sum(is.na(x)))
  nna
  if (any(nna > 0)) {
    message("[mofa.compute] warning: imputing missing values in X.")
    xdata[which(nna > 0)] <- lapply(xdata[which(nna > 0)], svdImpute2)
    ## X <- svdImpute2(X)
  }

  ## no dups
  xdata <- lapply(xdata, function(x) x[!duplicated(rownames(x)), ])

  ## Scale datatypes blocks? This is generally a good thing to do
  ## because the different datatypes may have different SD
  ## distributions.
  if (scale_views) {
    message("[mofa.compute] scaling blocks")
    xdata <- lapply(xdata, function(d) d - rowMeans(d, na.rm = TRUE))
    dt.sd <- sapply(xdata, function(d) {
      mean(matrixStats::rowSds(d, na.rm = TRUE),
        na.rm = TRUE
      )
    })
    xdata <- lapply(names(xdata), function(d) xdata[[d]] / (1e-4 + dt.sd[d]))
    names(xdata) <- names(dt.sd)
  }

  ## To be sure translate all feature names to ASCII (IK: should be
  ## done in creating genesets...)
  for (i in 1:length(xdata)) {
    d <- rownames(xdata[[i]])
    rownames(xdata[[i]]) <- iconv(d, to = "ascii//TRANSLIT")
  }

  ## add prefix???
  xdata <- mofa.prefix(xdata)
  X <- do.call(rbind, xdata)
  X <- X - rowMeans(X, na.rm = TRUE)
  dim(X)

  ## create phenotype levels from contrasts if missing. Note: this
  ## takes all contrasts in account, not a single comparison
  if (is.null(pheno)) {
    message("[mofa.compute] creating pheno from contrasts")
    ct <- contrasts
    ct[is.na(ct)] <- "NA"
    pheno <- apply(ct, 1, paste, collapse = "_")
    table(pheno)
  }

  model <- NULL
  kernel <- tolower(kernel)

  numfactors <- min(numfactors, ncol(X) - 1)

  if (kernel == "mofa") {
    message("[mofa.compute] computing MOFA factorization")
    obj <- MOFA2::create_mofa(xdata, groups = NULL)
    ## 'group' is used by MOFA
    samples1 <- samples
    samples1$group <- NULL
    MOFA2::samples_metadata(obj) <- data.frame(
      sample = rownames(samples1),
      samples1, check.names = FALSE
    )

    data_opts <- MOFA2::get_default_data_options(obj)
    data_opts$scale_views <- FALSE
    data_opts

    model_opts <- MOFA2::get_default_model_options(obj)
    model_opts$num_factors <- numfactors
    ## model_opts$likelihoods <- c("gaussian","gaussian","gaussian","bernoulli")

    train_opts <- MOFA2::get_default_training_options(obj)
    train_opts$gpu_mode <- gpu_mode
    train_opts$maxiter <- max_iter
    ## train_opts$gpu_mode <- FALSE

    obj <- MOFA2::prepare_mofa(
      object = obj,
      data_options = data_opts,
      model_options = model_opts,
      training_options = train_opts
    )

    ## where to save??? Do we need to save???
    outfile <- file.path(tempdir(), "mofa-model.hdf5")
    suppressMessages(suppressWarnings(
      model <- MOFA2::run_mofa(obj,
        outfile = outfile, save_data = TRUE,
        use_basilisk = TRUE
      )
    ))

    model <- MOFA2::load_model(outfile, remove_inactive_factors = FALSE)
    model <- MOFA2::impute(model)
    factors <- MOFA2::get_factors(model, factors = "all")
    weights <- MOFA2::get_weights(model, views = "all", factors = "all")
    W <- do.call(rbind, weights)
    F <- do.call(rbind, factors)
  } else if (tolower(kernel) %in% c("pca", "svd")) {
    message("[mofa.compute] computing PCA factorization")
    model <- irlba::irlba(X, nv = numfactors, maxit = max_iter, work = 100)
    W <- model$u
    F <- model$v %*% diag(model$d)
    rownames(W) <- rownames(X)
    rownames(F) <- colnames(X)
    colnames(W) <- paste0("PC", 1:ncol(W))
    colnames(F) <- paste0("PC", 1:ncol(F))
  } else if (tolower(kernel) %in% c("nmf", "nmf2")) {
    message("[mofa.compute] computing NMF factorization")
    if (kernel == "nmf2") {
      model <- mofa.intNMF(xdata,
        k = numfactors, method = "RcppML",
        separate.sign = TRUE,
        shift = FALSE, scale = TRUE,
        exponentiate = FALSE
      )
    } else {
      model <- mofa.intNMF(xdata,
        k = numfactors, method = "RcppML",
        separate.sign = FALSE,
        shift = TRUE, scale = TRUE,
        exponentiate = FALSE
      )
    }
    W <- do.call(rbind, model$W)
    F <- t(model$H)
    rownames(W) <- rownames(X)
    rownames(F) <- colnames(X)
    colnames(W) <- paste0("NMF", 1:ncol(W))
    colnames(F) <- paste0("NMF", 1:ncol(F))
  } else if (kernel %in% c(
    "diablo", "splsda", "rgcca", "sgcca", "sgccda",
    "rgcca.rgcca", "rgcca.rgccda", "rgcca.mcoa"
  )) {
    message(paste("[mofa.compute] computing", toupper(kernel), "factorization"))
    dtx <- sub(":.*", "", rownames(X))
    xx <- tapply(1:nrow(X), dtx, function(i) X[i, , drop = FALSE])
    if (is.null(pheno)) {
      stop("pheno must be provided for diablo")
    }
    if (any(is.na(pheno))) {
      message("[DIABLO] warning: missing values in phenotype. removing samples.")
      sel <- which(!is.na(pheno))
      pheno <- pheno[sel]
      xx <- lapply(xx, function(x) x[, sel, drop = FALSE])
      X <- X[, sel, drop = FALSE]
      samples <- samples[sel, ]
    }
    message("[mofa.compute] numfeatures = ", numfeatures)
    mix <- mixomics.compute_diablo(
      xx,
      y = pheno, method = kernel, ncomp = numfactors, nfeat = numfeatures
    )
    W <- mix$loadings
    F <- mix$scores
    model <- mix$res
  } else if (kernel == "mcia") {
    message("[mofa.compute] computing MCIA factorization")
    data_blocks <- lapply(xdata, t)
    data_blocks_mae <- nipalsMCIA::simple_mae(
      data_blocks,
      row_format = "sample",
      colData = samples
    )
    model <- nipalsMCIA::nipals_multiblock(
      data_blocks_mae,
      num_PCs = numfactors,
      plot = "none"
    )
    F <- model@global_scores
    W <- model@global_loadings
    dim(F)
    colnames(W) <- paste0("Factor", 1:ncol(W))
    colnames(F) <- paste0("Factor", 1:ncol(F))
  } else if (kernel == "wgcna") {
    message("[mofa.compute] computing WGCNA factorization")
    model <- wgcna.compute(
      X, samples,
      minmodsize = 10,
      power = 6,
      cutheight = 0.15,
      deepsplit = 2,
      networktype = "signed",
      tomtype = "signed",
      ngenes = -1
    )

    F <- as.matrix(model$net$MEs) ## factors (coefficients)
    ## we replace loading with feature-module (aka module-membership)
    ## correlation. the origonal loading (W) is not continuous.
    ## W <- model$W          ## loadings (but not continuous)
    W <- cor(model$datExpr, F)
  } else {
    message("[mofa.compute] FATAL invalid kernel:", kernel)
    return(NULL)
  }

  if (!is.null(factorname)) {
    colnames(W) <- paste0(factorname, 1:ncol(W))
    colnames(F) <- paste0(factorname, 1:ncol(F))
  }

  W <- W[rownames(X), , drop = FALSE]
  F <- F[colnames(X), , drop = FALSE]
  dt <- sub(":.*", "", rownames(W))
  ww <- tapply(1:nrow(W), dt, function(i) W[i, , drop = FALSE])
  xx <- tapply(1:nrow(X), dt, function(i) X[i, , drop = FALSE])
  ww <- ww[names(xdata)]
  xx <- xx[names(xdata)]

  ww_bysymbol <- ww
  if (!is.null(annot)) {
    ww_bysymbol <- mofa.prefix(ww)
    ww_bysymbol <- lapply(ww_bysymbol, function(w) rename_by2(w, annot, "symbol"))
  }
  ww <- mofa.strip_prefix(ww)
  xx <- mofa.strip_prefix(xx)

  ## variance explained
  w.ssq <- sapply(ww, function(w) colSums(w**2))
  w.ssq <- w.ssq * sqrt(colSums(F**2)) ## check!!
  V <- t(w.ssq / sum(w.ssq))
  rownames(V) <- names(xx)
  colnames(V) <- colnames(W)

  ## Covariate x Factor correlation
  Y <- samples
  Y$sample <- NULL
  Y$group <- NULL
  Y <- expandPhenoMatrix(Y, drop.ref = FALSE)
  Y <- as.matrix(Y)
  Z <- cor(Y, F, use = "pairwise")
  Z[is.na(Z)] <- 0
  colnames(Z) <- colnames(W)

  gsea <- NULL
  if (compute.enrichment) {
    message("computing factor enrichment...")
    symbolW <- do.call(rbind, mofa.prefix(ww_bysymbol))

    if (mean(rownames(symbolW) %in% rownames(GMT)) < 0.10) {
      message("WARNING! less than 10% of your features map to symbols. Please check or add a gene annotation table.")
    }

    dbg("[mofa.compute] dim.GMT = ", dim(GMT))
    dbg("[mofa.compute] dim.symbolW = ", dim(symbolW))
    dbg("[mofa.compute] rownames.symbolW = ", head(rownames(symbolW)))
    dbg("[mofa.compute] rownames.GMT = ", head(rownames(GMT)))

    gsea <- mofa.compute_enrichment(symbolW, G = GMT, ntop = gset.ntop)
  }

  ## create graphs
  graphs <- NULL
  if (compute.graphs) {
    graphs <- mofa.factor_graphs(F, W, X, y = pheno, n = numfeatures)
  }
  num.contrasts <- sign(makeContrastsFromLabelMatrix(contrasts))

  res <- list(
    model = model,
    samples = samples,
    contrasts = contrasts,
    pheno = pheno,
    K = num.contrasts,
    F = F,
    W = W,
    V = V,
    X = X,
    Z = Z,
    xx = xx,
    ww = ww,
    Y = Y,
    gsea = gsea,
    GMT = GMT,
    graphs = graphs
  )

  return(res)
}


#'
#'
#' @export
mofa.add_genesets <- function(xdata, GMT = NULL, datatypes = NULL) {
  names(xdata)
  if (is.null(datatypes)) {
    datatypes <- intersect(c("gx", "px", "mx"), names(xdata))
  }
  if (length(datatypes) == 0) {
    message("WARNING. No valid datatypes")
    return(xdata)
  }

  if (is.null(GMT)) {
    message("combining GSETxGENE and PATHBANK genesets...")
    stdGMT <- playdata::GSETxGENE
    colnames(stdGMT) <- paste0("SYMBOL:", colnames(stdGMT))
    pathbankGMT <- playdata::MSETxMETABOLITE
    GMT <- merge_sparse_matrix(pathbankGMT, stdGMT[, ])
    dim(GMT)
  } else {
    GMT <- Matrix::t(GMT) ## genesets on rows
  }

  if (is.null(GMT)) {
    message("ERROR: could not determing GMT. Please provide.")
    return(xdata)
  }

  ## convert non-ascii characters....
  ## rownames(GMT) <- iconv(rownames(GMT), "latin1", "ASCII", sub="")
  rownames(GMT) <- stringi::stri_trans_general(rownames(GMT), "latin-ascii")
  colnames(GMT) <- mofa.strip_prefix(colnames(GMT))

  gsetX <- list()
  dt <- datatypes[1]
  for (dt in datatypes) {
    rX <- xdata[[dt]] - rowMeans(xdata[[dt]], na.rm = TRUE)
    rownames(rX) <- mofa.strip_prefix(rownames(rX))
    sel <- intersect(rownames(rX), colnames(GMT))

    G <- GMT[, sel, drop = FALSE]
    rX <- rX[sel, , drop = FALSE]
    rho <- gset.rankcor(rX, G)$rho
    ii <- as.vector(which(rowMeans(is.na(rho)) < 1))
    rho <- rho[ii, , drop = FALSE]
    gsetX[[dt]] <- rho
  }
  lapply(gsetX, dim)

  ## make sure they align. Only intersection????
  names(gsetX) <- paste0("gset.", names(gsetX))
  gsetX
}


#' Compute clusters for MOFA factors.
#'
#' @export
mofa.compute_clusters <- function(xx, matF = NULL, along = "samples") {
  if (along == "samples") {
    ## cluster samples
    sxx <- lapply(xx, function(d) {
      d <- t(scale(t(d), scale = FALSE))
      d[is.na(d)] <- 0
      d
    })
    if (!is.null(matF)) sxx$MOFA <- t(matF)
    sxx <- lapply(sxx, function(d) d[matrixStats::rowSds(d, na.rm = TRUE) > 0.01, ])
    nb <- max(round(min(15, ncol(xx[[1]]) / 4)), 2)
    posx <- lapply(sxx, function(x) uwot::umap(t(x), n_neighbors = nb))
    for (i in 1:length(posx)) rownames(posx[[i]]) <- colnames(sxx[[i]])
  }

  if (along == "features") {
    ## cluster features
    sxx <- lapply(xx, function(d) {
      d <- t(scale(t(d), scale = FALSE))
      d[is.na(d)] <- 0
      d
    })
    posx <- list()
    for (i in 1:length(sxx)) {
      nb <- max(min(15, nrow(sxx[[i]]) / 4), 2)
      posx[[i]] <- uwot::umap(sxx[[i]], n_neighbors = nb)
      rownames(posx[[i]]) <- rownames(sxx[[i]])
    }
    names(posx) <- names(sxx)
  }

  posx <- lapply(posx, playbase::uscale)
  return(posx)
}


#' Compute enrichment for MOFA factors.
#'
#' @export
mofa.compute_enrichment <- function(W, G = NULL, filter = NULL, ntop = 1000) {
  if (is.null(G)) {
    info("[pgx.add_GMT] Adding metabolomics genesets")
    G1 <- Matrix::t(playdata::MSETxMETABOLITE)
    info("[pgx.add_GMT] Adding transcriptomics/proteomics genesets")
    G2 <- Matrix::t(playdata::GSETxGENE)
    rownames(G2) <- paste0("SYMBOL:", rownames(G2))
    G <- merge_sparse_matrix(G1, G2)
  } else {
    has.colons <- all(grepl(":", rownames(G)))
    has.colons
    if (!has.colons) {
      ## defaults prefix. If starts with letter SYMBOL, if all numbers CHEBI
      rownames(G) <- mofa.strip_prefix(rownames(G))
      ii <- grep("^[0-9]+$", rownames(G), ignore.case = TRUE)
      if (length(ii) > 0) rownames(G)[ii] <- paste0("CHEBI:", rownames(G)[ii])
      jj <- setdiff(1:nrow(G), ii)
      if (length(jj) > 0) rownames(G)[jj] <- paste0("SYMBOL:", rownames(G)[jj])
    }
  }

  dbg("[mofa.compute_enrichment] dim.G = ", dim(G))

  ## filter gene sets
  if (!is.null(filter)) {
    sel <- grep(filter, colnames(G))
    G <- G[, sel]
  }

  ## detect datatypes. associate with symbols
  dtypes <- unique(sub(":.*", "", rownames(W)))
  dbg("[mofa.compute_enrichment] dtypes = ", dtypes)
  dtype2gtype <- c(
    "gx" = "SYMBOL", "px" = "SYMBOL",
    "me" = "SYMBOL", "mx" = "CHEBI"
  )
  dtypes <- intersect(dtypes, names(dtype2gtype))
  dtypes

  ## build full gene set sparse matrix
  GMT <- G[0, ]
  dt <- dtypes[1]
  for (dt in dtypes) {
    gt <- dtype2gtype[dt]
    ii <- grep(gt, rownames(G))
    G1 <- G[ii, , drop = FALSE]
    rownames(G1) <- sub(gt, dt, rownames(G1))
    GMT <- rbind(GMT, G1)
  }
  table(sub(":.*", "", rownames(GMT)))

  dbg("[mofa.compute_enrichment] dim.GMT = ", dim(GMT))
  dbg("[mofa.compute_enrichment] rownames.GMT = ", head(rownames(GMT)))
  dbg("[mofa.compute_enrichment] rownames.W = ", head(rownames(W)))

  pp <- intersect(rownames(W), rownames(GMT))
  W <- W[pp, ]
  sel <- which(Matrix::colSums(GMT[pp, ] != 0) >= 3)

  dbg("[mofa.compute_enrichment] lenght(pp) = ", length(pp))
  dbg("[mofa.compute_enrichment] lenght(sel) = ", length(sel))

  GMT1 <- GMT[pp, sel]
  dim(GMT1)

  dbg("[mofa.compute_enrichment] dim.GMT1 = ", dim(GMT1))

  ## Perform geneset enrichment with fast rank-correlation. We could
  ## do with fGSEA instead but it is much slower.
  info("[pgx.add_GMT] Performing rankcor")
  normW <- apply(W, 2, function(x) normalize_multirank(x))
  normW[is.na(normW)] <- mean(normW, na.rm = TRUE)
  names(normW) <- rownames(W)
  rnk <- gset.rankcor(normW, GMT1, compute.p = TRUE)
  topsel <- apply(abs(rnk$rho), 2, function(r) head(order(-r), ntop))
  topsel <- head(unique(as.vector(t(topsel))), ntop)
  length(topsel)
  topsets <- rownames(rnk$rho)[topsel]

  ## prioritize top genesets to accelerate GSEA
  rnk <- lapply(rnk, function(x) x[topsets, ])
  topGMT <- GMT1[, topsets]
  gmt <- mat2gmt(topGMT)
  length(gmt)
  info("[pgx.add_GMT] Performing fGSEA...")
  gsea <- list()
  k <- colnames(W)[1]
  for (k in colnames(W)) {
    w <- W[, k]
    w <- w + 1e-4 * rnorm(length(w))
    w <- normalize_multirank(w)
    enr <- fgsea::fgsea(gmt, stats = w)
    enr <- data.frame(enr)
    rownames(enr) <- enr$pathway
    gsea[[k]] <- enr
  }
  gsea <- lapply(gsea, function(x) x[rownames(gsea[[1]]), ])

  dim(topGMT)
  dtype <- sub(":.*", "", rownames(topGMT))
  size <- tapply(1:nrow(topGMT), dtype, function(ii) Matrix::colSums(topGMT[ii, ] != 0))
  size <- do.call(cbind, size)
  size <- size[rownames(gsea[[1]]), ]
  for (i in 1:length(gsea)) {
    gsea[[i]]$size <- size
    gsea[[i]] <- gsea[[i]][, c("pathway", "NES", "pval", "padj", "size", "leadingEdge")]
    gsea[[i]] <- data.frame(lapply(gsea[[i]], as.matrix), check.names = FALSE)
  }

  ## extract rho, pval and qval
  rho <- sapply(gsea, function(x) x$NES)
  pval <- sapply(gsea, function(x) x$pval)
  padj <- sapply(gsea, function(x) x$padj)
  rownames(rho) <- rownames(gsea[[1]])
  rownames(pval) <- rownames(gsea[[1]])
  rownames(padj) <- rownames(gsea[[1]])

  list(table = gsea, NES = rho, pval = pval, padj = padj)
}


#'
#'
#' @export
pgx.compute_mofa_importance <- function(pgx, pheno,
                                        numfactors = 8, use.sdwt = TRUE,
                                        kernels = c(
                                          "mofa", "pca", "nmf", "nmf2", "mcia",
                                          "wgcna", "diablo", "rgcca",
                                          "rgcca.rgcca", "rgcca.rgccda",
                                          "rgcca.mcoa"
                                        )) {
  if (!"factorizations" %in% names(pgx$mofa)) {
    message("WARNING: no factorizatons in pgx. Computing factorizations...")
    xdata <- mofa.split_data(pgx$X)
    names(xdata)
    samples <- pgx$samples
    contrasts <- pgx$contrasts
    meta.res <- mofa.compute_factorizations(
      xdata,
      samples,
      contrasts,
      numfactors = numfactors,
      kernels = kernels
    )
  } else {
    meta.res <- pgx$mofa$factorizations
  }
  sdx <- NULL
  if (use.sdwt) {
    sdx <- matrixStats::rowSds(pgx$X) ## important
    sdx <- sdx**0.5 ## moderate weight
  }
  res <- mofa.compute_meta_importance(
    meta.res, pheno,
    sd.wt = sdx, plot = FALSE
  )

  return(res)
}


#'
#'
#' @export
mofa.compute_factorizations <- function(xdata,
                                        samples,
                                        contrasts,
                                        numfactors = 8,
                                        kernels =
                                          c(
                                            "mofa", "pca", "nmf", "nmf2",
                                            "mcia", "wgcna", "diablo", "rgcca",
                                            "rgcca.rgcca", "rgcca.rgccda",
                                            "rgcca.mcoa"
                                          )) {
  meta.res <- list()
  for (kernel in kernels) {
    cat("[mofa.compute_factorizations] computing for kernel", toupper(kernel), "\n")
    meta.res[[kernel]] <- mofa.compute(
      xdata,
      samples,
      contrasts = contrasts,
      GMT = NULL,
      annot = NULL,
      pheno = NULL,
      kernel = kernel,
      scale_views = TRUE,
      gpu_mode = FALSE,
      ntop = 2000,
      max_iter = 200,
      numfactors = numfactors,
      numfeatures = 30,
      compute.enrichment = FALSE,
      compute.graphs = FALSE
    )
  }
  meta.res <- lapply(meta.res, function(r) r[c("W", "F", "Z")])
  return(meta.res)
}


#'
#'
#' @export
mofa.plot_meta_importance <- function(meta.res, pheno, sd.wt = NULL,
                                      power = 6, plot = FALSE) {
  mofa.compute_meta_importance(
    meta.res = meta.res,
    pheno = pheno,
    sd.wt = sd.wt,
    power = power,
    plot = TRUE
  )
}


#'
#'
#' @export
mofa.compute_meta_importance <- function(meta.res, pheno, sd.wt = NULL,
                                         power = 6, plot = FALSE) {
  # pheno="condition=Her2"
  if (is.null(pheno)) {
    message("ERROR: must provide pheno")
    return(NULL)
  }

  zpheno <- rownames(meta.res[[1]]$Z)
  zpheno

  if (0) {
    meta.res[[1]]$Z
    pheno <- zpheno[2]
    pheno <- zpheno[1:2]
    pheno
  }

  if (!grepl("=", pheno)) {
    zpheno <- rownames(meta.res[[1]]$Z)
    pheno <- grep(pheno, zpheno, value = TRUE)
  }

  pheno <- intersect(pheno, zpheno)
  if (length(pheno) == 0) {
    message("ERROR: could not match pheno")
    return(NULL)
  }

  if (length(pheno) == 1) {
    zz <- lapply(meta.res, function(r) r$Z[pheno[1], ])
    topk <- sapply(zz, function(z) names(which.max(abs(z))))
    topk
    topsign <- sapply(names(zz), function(i) sign(zz[[i]][topk[i]]))
    names(topsign) <- names(zz)
    topsign
    ww <- lapply(meta.res, function(res) res$W)
    B <- sapply(names(topk), function(i) topsign[i] * ww[[i]][, topk[i]])
    if (!is.null(sd.wt)) B <- B * sd.wt[rownames(B)]
  } else {
    bb <- list()
    i <- 1
    for (i in 1:length(meta.res)) {
      Z <- meta.res[[i]]$Z[pheno, ]
      Z <- tanh(Z / sd(Z))
      W <- meta.res[[i]]$W
      if (!is.null(sd.wt)) W <- W * sd.wt[rownames(W)]
      K <- (W %*% t(Z))
      ## this is heuristic to equalize contribution of columns
      R <- apply(K, 2, function(x) sign(x) * rank(abs(x)) / length(x))
      wt <- apply(abs(R), 1, max)
      bb[[i]] <- wt
    }
    names(bb) <- names(meta.res)
    B <- do.call(cbind, bb)
  }

  ## Convert to weighted elevated rank. Highest rank is more
  ## important.
  R <- (apply(B, 2, rank) / nrow(B))**power
  R <- R[order(-rowMeans(R)), ]
  head(R)

  ## return R matrix if not plotting
  if (plot == FALSE) {
    return(R)
  }

  par(mfrow = c(1, 1), mar = c(12, 4.5, 2, 0))
  barplot(t(head(R, 60)),
    las = 3,
    ylab = "cumulative ranking", main = "variable importance"
  )
  legend("topright",
    legend = rev(colnames(R)),
    fill = rev(grey.colors(ncol(R))),
    cex = 0.8, y.intersp = 0.8, bg = "white",
    inset = c(0.035, 0)
  )
}


#' @export
mofa.plot_all_factortraits <- function(meta.res) {
  nr <- ceiling(sqrt(length(meta.res)))
  nc <- ceiling(length(meta.res) / nr)
  par(mfrow = c(nr, nc), mar = c(6, 1, 2, 10))
  for (i in 1:length(meta.res)) {
    gx.imagemap(abs(meta.res[[i]]$Z), main = names(meta.res)[i])
  }
}

## -------------------------------------------------------------
## ------------------ HELPER FUNCTIONS -------------------------
## -------------------------------------------------------------

#' @export
mofa.topsd <- function(xdata, ntop) {
  lapply(xdata, function(x) {
    sdx <- matrixStats::rowSds(x, na.rm = TRUE)
    head(x[order(-sdx), ], ntop)
  })
}

#' @export
mofa.scale_views <- function(xdata) {
  xdata <- lapply(xdata, function(d) d - rowMeans(d, na.rm = TRUE))
  dt.sd <- sapply(xdata, function(d) {
    mean(matrixStats::rowSds(d, na.rm = TRUE),
      na.rm = TRUE
    )
  })
  dt.sd
  xdata <- lapply(names(xdata), function(d) xdata[[d]] / (1e-9 + dt.sd[d]))
  names(xdata) <- names(dt.sd)
  xdata
}

#' @export
mofa.split_data <- function(X, keep.prefix = FALSE) {
  if (!all(grepl("[:]|SOURCE|SINK", rownames(X)))) {
    ## if no prefix, then single-omics. Assume gene expression.
    rownames(X) <- paste0("gx:", rownames(X))
  }
  dtype <- sub(":.*", "", rownames(X))
  xx <- tapply(1:nrow(X), dtype, function(i) X[i, , drop = FALSE])
  if (!keep.prefix) {
    xx <- mofa.strip_prefix(xx)
  }
  xx
}

#' @export
mofa.merge_data <- function(xx) {
  do.call(rbind, mofa.prefix(xx))
}

#' @export
mofa.prefix <- function(xx) {
  xx <- mofa.strip_prefix(xx)
  i <- 1
  for (i in 1:length(xx)) {
    dt <- paste0(names(xx)[i], ":")
    if (is.null(dim(xx[[i]]))) {
      names(xx[[i]]) <- paste0(dt, names(xx[[i]]))
    } else {
      rownames(xx[[i]]) <- paste0(dt, rownames(xx[[i]]))
    }
  }
  xx
}

#' @export
mofa.get_prefix <- function(x) {
  ifelse(grepl(":", x), sub(":.*", "", x), "")
}

#' @export
mofa.strip_prefix <- function(xx) {
  if (class(xx) == "character") {
    xx <- sub("[A-Za-z]+:", "", xx)
    return(xx)
  }
  if (class(xx) == "matrix") {
    rownames(xx) <- sub("[A-Za-z]+:", "", rownames(xx))
    return(xx)
  }
  if (class(xx) %in% c("list", "array") || is.list(xx)) {
    i <- 1
    for (i in 1:length(xx)) {
      dt <- paste0("^", names(xx)[i], ":")
      if (is.null(dim(xx[[i]]))) {
        names(xx[[i]]) <- sub(dt, "", names(xx[[i]]))
      } else {
        rownames(xx[[i]]) <- sub(dt, "", rownames(xx[[i]]))
      }
    }
    return(xx)
  }
  xx
}

#' @export
mofa.augment <- function(xx, n, z = 1) {
  sdxx <- sapply(xx, function(x) mean(matrixStats::rowSds(x)))
  ## enlarge
  for (i in 1:length(xx)) {
    rn <- rownames(xx[[i]])
    cn <- colnames(xx[[i]])
    xx[[i]] <- matrix(xx[[i]], nrow(xx[[i]]), n * ncol(xx[[i]]))
    rownames(xx[[i]]) <- rn
    colnames(xx[[i]]) <- rep(cn, n)
  }
  ## add noise
  for (i in 1:length(xx)) {
    rxx <- matrix(rnorm(length(xx[[i]])), nrow(xx[[i]]), ncol(xx[[i]]))
    xx[[i]] <- xx[[i]] + z * sdxx[[i]] * rxx
  }
  xx
}

#' @export
mofa.correlate_covariate <- function(res, Y) {
  K <- cov(res$F, Y)
  R <- cov(t(res$W), K)
  rr <- lapply(res$ww, function(n) R[rownames(n), , drop = FALSE])
  rr <- mofa.strip_prefix(rr)
  list(
    R = R,
    rr = rr
  )
}


## -------------------------------------------------------------
## -----------------PLOTTING FUNCTIONS -------------------------
## -------------------------------------------------------------

#' @export
mofa.plot_weights <- function(weights, k, ntop = 10, cex.names = 0.9,
                              maxchar = 999) {
  wt <- lapply(weights, function(w) w[, k])
  ## plot the weights for each datatype (Factor1)
  ## par(mfrow=c(2,3), mar=c(4,30,4,2))
  v <- names(wt)[1]
  for (v in names(wt)) {
    w1 <- wt[[v]]
    names(w1) <- gsub("[a-zA-Z]+:| \\(SMP.*", "", names(w1))
    w1 <- w1[order(-abs(w1))]
    w1 <- w1[!duplicated(names(w1))]
    names(w1) <- stringr::str_trunc(names(w1), maxchar)
    xlab <- paste(k, " weight")
    if (!is.null(ntop) && ntop > 0 && ntop < length(w1)) {
      w1 <- w1[w1 != 0]
      i0 <- head(order(w1), ntop)
      i1 <- tail(order(w1), ntop)
      if (length(i0) == 0) i1 <- tail(order(w1), 2 * ntop)
      if (length(i1) == 0) i0 <- head(order(w1), 2 * ntop)
      ii <- unique(c(i0, i1))
      barplot(w1[ii],
        horiz = TRUE, las = 1,
        xlab = xlab, cex.names = cex.names
      )
    } else {
      names(w1) <- NULL
      barplot(sort(w1),
        border = NA, col = "grey50",
        horiz = TRUE, las = 1, xlab = xlab, space = 0
      )
      abline(v = 0)
      mtext(paste0("sorted features (N=", length(w1), ")"), side = 2)
    }
    title(toupper(v), cex.main = 1.1)
  }
}

#' @export
mofa.plot_heatmap <- function(mofa,
                              gene_table = NULL,
                              k = NULL, ntop = 50,
                              features = NULL,
                              type = c("heatmap", "splitmap")[1],
                              annot = c("scores", "pheno")[1],
                              split = TRUE, maxchar = 999,
                              show_types = NULL,
                              mar = c(8, 12), annot.ht = 1,
                              show_legend = TRUE,
                              show_colnames = TRUE,
                              cexRow = 0.9) {
  if (!is.null(k)) {
    xx <- mofa.prefix(mofa$ww)
  } else {
    xx <- mofa.prefix(mofa$xx)
  }

  if (!is.null(show_types)) {
    xx <- xx[intersect(show_types, names(xx))]
  }

  if (!is.null(features)) {
    features2 <- mofa.strip_prefix(features)
    for (i in 1:length(xx)) {
      namesx <- mofa.strip_prefix(rownames(xx[[i]]))
      ii <- which(namesx %in% features2)
      xx[[i]] <- xx[[i]][ii, , drop = FALSE]
    }
  }

  ntop1 <- ntop / length(xx)

  if (!is.null(k)) {
    top <- lapply(xx, function(w) {
      w1 <- w[, k]
      names(w1) <- rownames(w)
      w1 <- w1[which(w1 != 0 & !is.na(w1))]
      head(names(sort(-abs(w1))), ntop1)
    })
  } else {
    top <- lapply(xx, function(x) {
      sdx <- matrixStats::rowSds(x, na.rm = TRUE)
      head(names(sort(-sdx)), ntop1)
    })
  }
  top <- unlist(top)
  topX <- mofa$X[top, , drop = FALSE]
  rownames(topX) <- stringr::str_trunc(rownames(topX), maxchar)

  if (annot == "scores") {
    aa <- data.frame(mofa$F)
  }
  if (annot == "pheno") {
    aa <- mofa$samples
  }
  rownames(aa) <- colnames(topX)

  if (!is.null(gene_table)) {
    topX <- rename_by2(topX, gene_table, "symbol", keep.prefix = TRUE)
  }

  if (type == "heatmap") {
    par(mar = c(0, 0, 0, 0))
    gx.heatmap(topX,
      mar = c(8, 10),
      key = FALSE,
      keysize = 1,
      cexRow = cexRow,
      col.annot = aa,
      show_colnames = show_colnames,
      annot.ht = 0.88 * annot.ht
    )
  } else if (type == "splitmap") {
    dx <- sub(":.*", "", rownames(topX))
    if (!split) dx <- NULL

    gx.splitmap(topX,
      split = dx,
      na_col = "white",
      softmax = TRUE,
      rowlab.maxlen = 80,
      rownames_width = 140,
      show_legend = show_legend,
      show_key = FALSE,
      col.annot = aa,
      show_colnames = show_colnames,
      annot.ht = 3.6 * annot.ht
    )
  } else if (type == "graph") {

  } else {
    message("unknown graph type")
  }
}

#'
#'
#'
#' @export
mofa.plot_factor_trait_correlation <- function(mofa,
                                               main = "Factor-Trait correlation",
                                               type = c("wgcna", "splitmap"),
                                               par = TRUE, cex_text = NULL,
                                               collapse = FALSE,
                                               cluster = TRUE,
                                               features = NULL,
                                               ...) {
  ## type = c("wgcna","splitmap")
  type <- type[1]

  ## Covariate x Factor correlation
  Z <- mofa$Z

  if (!is.null(features)) {
    sel <- intersect(features, rownames(mofa$W))
    Z <- rbind(Z, mofa$W[sel, , drop = FALSE])
  }

  if (collapse) {
    rownames(Z) <- sub("=.*", "", rownames(Z))
    Z <- rowmean(Z**2)**0.5
  }

  if (nrow(Z) == 1) {
    Z <- rbind(Z, " " = Z[1, ]) ## hack for single row...
  }

  if (cluster) {
    ii <- hclust(dist(Z))$order
    jj <- hclust(dist(t(Z)))$order
    Z <- Z[ii, jj]
  }

  if (type == "splitmap") {
    gx.splitmap(
      t(Z),
      nmax = 50, scale = "none", main = main,
      col.annot = mofa$samples, split = 1,
      show_legend = FALSE, show_key = FALSE
    )
  }

  if (type == "wgcna") {
    ## par(mfrow=c(1,1))
    if (par) par(mar = c(6, 5, 2, 1))
    ftext <- round(t(Z), digits = 2)
    if (is.null(cex_text)) {
      cex.text <- min(0.8, max(0.3, 6 / ncol(Z)))
    } else {
      cex.text <- cex_text
    }

    WGCNA::labeledHeatmap(
      Matrix = t(Z),
      xLabels = rownames(Z),
      yLabels = colnames(Z),
      textMatrix = ftext,
      cex.text = cex.text,
      colorLabels = TRUE,
      colors = WGCNA::blueWhiteRed(50),
      setStdMargins = FALSE,
      zlim = c(-1, 1),
      main = main
    )
  }
}


#' @export
mofa.plot_factor_corheatmap <- function(mofa,
                                        main = "Factor correlation heatmap",
                                        marx = 1,
                                        ...) {
  R <- cor(mofa$F, use = "pairwise")
  lablength <- max(nchar(rownames(R)))
  gx.heatmap(
    R,
    nmax = 50, scale = "none", main = main, sym = TRUE,
    key = FALSE, keysize = 0.8, mar = c(0.8, 1) * (2 + marx * lablength),
    ...
  )
}

#' @export
mofa.plot_factor_boxplots <- function(mofa, k = 1, pheno = NULL,
                                      by = "condition", par = TRUE) {
  samples <- mofa$samples

  if (by == "condition") {
    ## boxplots
    if (par) par(mfrow = c(2, 2), mar = c(5, 5, 3, 2))
    j <- 1
    if (is.integer(k)) k <- colnames(mofa$F)[k]
    for (j in 1:ncol(samples)) {
      cond <- colnames(samples)[j]
      y <- factor(samples[, j])
      f1 <- mofa$F[, k]
      tt <- paste(k, "by", colnames(samples)[j])
      boxplot(f1 ~ y,
        main = tt, xlab = cond,
        ylab = paste(fn, "score")
      )
    }
  }

  if (by == "factor") {
    k <- 1
    if (!is.null(pheno) && is.integer(pheno)) pheno <- colnames(res$samples)[pheno]
    if (par) par(mfrow = c(3, 4))
    for (k in 1:min(12, ncol(mofa$F))) {
      y <- factor(samples[, pheno])
      f1 <- mofa$F[, k]
      tt <- colnames(mofa$F)[k]
      boxplot(f1 ~ y, main = tt, xlab = pheno, ylab = "Score")
    }
  }
}


mofa.combine_layers <- function(xx, weights = 1) {
  k <- 1
  if (length(weights) == 1) weights <- rep(weights, length(xx))
  for (i in 1:length(xx)) {
    rownames(xx[[i]]) <- mofa.strip_prefix(rownames(xx[[i]])) ## strip
  }
  kk <- Reduce(intersect, lapply(xx, colnames))
  gg <- Reduce(intersect, lapply(xx, rownames))
  xx <- lapply(xx, function(x) x[gg, kk])
  cx.pos <- xx[[1]] * 0 + 1
  cx.neg <- xx[[1]] * 0 + 1
  for (i in 1:length(xx)) {
    xi <- xx[[i]] - rowMeans(xx[[i]], na.rm = TRUE)
    cx.pos <- cx.pos * pmax(xi, 0) * weights[i]
    cx.neg <- cx.neg * pmax(-xi, 0) * weights[i]
  }
  cx <- cx.pos - cx.neg
  cx ## combined
}


#' @export
mofa.plot_enrichment <- function(gsea, type = "barplot",
                                 remove.dup = TRUE, ntop = 20,
                                 filter = NULL, select = NULL,
                                 strip.names = FALSE, par = NULL,
                                 sort = TRUE,
                                 title = "multi-omics enrichment") {
  ## k=1;ntop=20
  S <- gsea
  if (nrow(S) == 0) {
    return(NULL)
  }

  if (!is.null(select) && length(select) > 0) {
    S <- S[select, , drop = FALSE]
  }

  if (!is.null(filter)) {
    S <- S[grep(filter, rownames(S)), , drop = FALSE]
  }

  if (remove.dup) {
    sname <- gsub("[a-zA-Z]+:| \\(.*", "", rownames(S))
    S <- S[which(!duplicated(sname)), , drop = FALSE]
  }
  S <- head(S, ntop)
  topS <- S$NES
  names(topS) <- S$pathway
  plot <- NULL
  if (strip.names) {
    names(topS) <- gsub(" \\(.*", "", names(topS))
  }
  if (sort) topS <- sort(topS)

  if (type == "barplot") {
    if (is.null(par)) {
      par(mfrow = c(1, 2), mar = c(4, 10, 4, 2))
      plot.new()
    } else {
      par(par)
    }
    barplot(rev(topS), horiz = TRUE, las = 1, xlab = "enrichment (NES)")
    title(title, outer = TRUE, cex.main = 1.4, line = -1.2)
  }

  if (type == "lollipop") {
    ## need double lolliplot lines for multi-omics
    plot <- ggLollipopPlot(topS, xlab = "enrichment (NES)")
  }
  plot
}


mofa.plotVar <- function(mofa, comp = 1:2, style = "correlation",
                         textlabel = TRUE) {
  if (style == "correlation") {
    xdata <- mofa$ww
    dtypes <- names(mofa$ww)
    for (i in 1:length(xdata)) {
      x <- mofa$xx[[i]]
      w <- mofa$ww[[i]][, comp]
      f <- mofa$F[, comp]
      sel <- which(rowSums(abs(w)) > 0) ## only non-zero variables
      rho <- cor(t(x[sel, , drop = FALSE]), f)
      if (i == 1) {
        plot(rho[, 1], rho[, 2],
          pch = 19, col = 1 + i,
          xlab = colnames(rho)[1],
          ylab = colnames(rho)[2],
          asp = 1, xlim = c(-1, 1), ylim = c(-1, 1)
        )
      } else {
        points(rho[, 1], rho[, 2], pch = 19, col = 1 + i)
      }
      if (textlabel) {
        text(rho[, 1], rho[, 2], rownames(rho), pos = 3)
      }
    }
    title("Correlation Circle Plot")
    abline(v = 0, h = 0, lty = 2)
    plotrix::draw.circle(0, 0, radius = c(0.5, 1))
    plotrix::draw.circle(0, 0, radius = 0.5)
    legend("bottomright",
      legend = dtypes,
      fill = 2:100, cex = 0.9, y.intersp = 0.85
    )
  }

  if (style == "loading") {
    dtypes <- names(mofa$ww)
    ww <- lapply(mofa$ww, function(w) w[, comp, drop = FALSE])
    maxw <- max(sapply(ww, function(w) max(abs(w), na.rm = TRUE)))
    for (i in 1:length(dtypes)) {
      w <- ww[[i]][, comp]
      sel <- which(rowSums(abs(w)) > 0) ## only non-zero variables
      if (i == 1) {
        plot(w[, 1], w[, 2],
          pch = 19, col = 1 + i,
          xlab = colnames(w)[1],
          ylab = colnames(w)[2],
          asp = 1, xlim = c(-1, 1) * maxw,
          ylim = c(-1, 1) * maxw
        )
      } else {
        points(w[, 1], w[, 2], pch = 19, col = 1 + i)
      }
      text(w[, 1], w[, 2], rownames(w), pos = 3)
    }
    title("Variable weights (Factor loading)")
    abline(v = 0, h = 0, lty = 2)
    legend("bottomright",
      legend = dtypes,
      fill = 2:100, cex = 0.9, y.intersp = 0.85
    )
  }
}

#' @export
mofa.factor_graphs <- function(F, W, X, y, n = 100, ewidth = 1, vsize = 1) {
  ## n=100;ewidth=1;vsize=1

  ## cluster-reduced graph
  create_graph <- function(gx, y, min.cor) {
    rho <- cor(t(gx), use = "pairwise")
    gr <- igraph::graph_from_adjacency_matrix(
      rho,
      weighted = TRUE, diag = FALSE, mode = "undirected"
    )
    ew <- igraph::E(gr)$weight
    ew <- (ew - min(ew, na.rm = TRUE)) /
      (max(ew, na.rm = TRUE) - min(ew, na.rm = TRUE))
    igraph::E(gr)$weight <- ew
    gr <- igraph::delete_edges(gr, edges = which(ew < min.cor))
    gr <- igraph::simplify(gr)
    sdx <- matrixStats::rowSds(gx)
    igraph::V(gr)$size <- 20 * vsize * (sdx / max(sdx))
    if (!is.null(y)) {
      y <- factor(as.character(y))
      ii <- which(!is.na(y))
      if (length(unique(y[ii])) == 2) {
        val <- cor(t(gx[, ii, drop = FALSE]), as.integer(y[ii]))[, 1]
      } else if (length(unique(y[ii])) > 2) {
        res <- gx.limmaF(gx[, ii, drop = FALSE], y[ii], fdr = 1, lfc = 0)
        val <- res$logFC
      } else {
        val <- NULL
      }
      # igraph::V(gr)$size <- 20*vsize*(abs(val)/max(abs(val)))
      if (!is.null(val)) igraph::V(gr)$color <- colorscale(val, gamma = 2)
    }
    igraph::E(gr)$width <- 3 * ewidth * abs(igraph::E(gr)$weight)**2
    gr
  }

  meGraph <- create_graph(t(F), y, min.cor = 0.33)
  me.size <- rowSums(abs(t(W)) > 0.5 * apply(abs(W), 2, max))
  igraph::V(meGraph)$size <- 25 * vsize * (me.size / max(me.size))**0.5

  topfeatures <- function(x, n) {
    names(x) <- rownames(W)
    x <- sort(x[x != 0])
    unique(c(names(tail(x, n)), names(head(x, n))))
  }
  topff <- apply(W, 2, function(x) topfeatures(x, n = n), simplify = FALSE)
  subgraphs <- list()
  k <- 1
  for (k in 1:length(topff)) {
    sel <- topff[[k]]
    if (length(sel) > 1) {
      gx <- X[sel, , drop = FALSE]
      subgraphs[[k]] <- create_graph(gx, y, min.cor = 0.33)
    } else {
      subgraphs[[k]] <- NULL
    }
  }
  names(subgraphs) <- colnames(W)
  list(
    factors = meGraph,
    features = subgraphs
  )
}

#' @export
mofa.plot_centrality <- function(res, k, show_types = NULL, transpose = FALSE,
                                 main = "centrality vs. factor weight",
                                 justdata = FALSE) {
  ## compute centrality
  gr <- res$graph$features[[k]]
  gg <- igraph::V(gr)$name
  gr <- igraph::mst(gr, weight = 1 / igraph::E(gr)$weight)
  ctx <- igraph::page_rank(gr)$vector
  names(ctx) <- gg

  rx <- res$W[, k]
  names(rx) <- rownames(res$W)

  gg <- intersect(names(ctx), names(rx))
  if (!is.null(show_types)) {
    dt <- sub(":.*", "", gg)
    gg <- gg[which(dt %in% show_types)]
  }
  rx <- rx[gg]
  ry <- ctx[gg]

  if (justdata) {
    df <- data.frame(feature = names(rx), weight = rx, centrality = ry)
    return(df)
  }

  dx <- 0.06 * c(-1, 1) * diff(range(rx))
  if (transpose) {
    plot(ry, rx,
      pch = 20, cex = 1,
      xlim = c(0, 1.05 * max(ry)),
      ylim = range(rx) + dx,
      xlab = "centrality  (page_rank)",
      ylab = paste0("logFC  (", test.type, ")")
    )
  } else {
    plot(rx, ry,
      pch = 20, cex = 1,
      ylim = c(0, 1.05 * max(ry)),
      xlim = range(rx) + dx,
      ylab = "centrality  (page_rank)",
      xlab = "factor weight"
    )
  }
  title(main)
  abline(h = 0, v = 0, lty = 2)
  nrx <- rx / mean(abs(rx))
  nry <- ry / mean(abs(ry))
  ii <- head(order(-(nrx**2 + nry**2)), 40)
  text(rx[ii], ry[ii], gg[ii], cex = 0.85, pos = 3, offset = 0.3)
}

#' @export
mofa.plot_module <- function(graph, mst = TRUE, nlabel = 10, rm.single = FALSE,
                             highlightby = "centrality.prize", cex = 1,
                             physics = TRUE,
                             plotlib = c("igraph", "visnet")) {
  if (rm.single) {
    vdeg <- igraph::degree(graph)
    graph <- igraph::subgraph(graph, vids = which(vdeg > 0))
  }
  if (mst) {
    wt <- igraph::E(graph)$weight
    graph <- igraph::mst(graph, weights = 1 / wt)
  }

  igraph::V(graph)$prize <- igraph::V(graph)$size
  igraph::V(graph)$type <- 1
  igraph::V(graph)$foldchange <- 1
  ew <- igraph::E(graph)$width
  igraph::E(graph)$width <- 10 * playbase::uscale(ew)**2

  if (plotlib == "igraph") {
    plot(graph)
  }
  if (plotlib == "visnet") {
    class(graph) <- c("PCSF", "igraph")
    igraph::E(graph)$weight <- 1 / igraph::E(graph)$weight
    plotPCSF(
      graph,
      highlightby = highlightby,
      nlabel = nlabel,
      edge_width = 5 * cex,
      node_cex = max(10, 30 * cex),
      label_cex = 25 * cex,
      layout = "layout_with_kk",
      physics = physics,
      plotlib = "visnet"
    )
  }
}

#' @export
mofa.predict <- function(mofa, newdata) {
  if (!is.null(mofa$model) && "block.splsda" %in% class(mofa$model)) {
    message("This seems to be a MixOmics Diablo model. Please use mixomics.predict() for better performance.")
  }

  if (is.list(newdata)) {
    newdata <- mofa.scale_views(newdata)
    newdata <- mofa.prefix(newdata)
    newdata <- do.call(rbind, newdata)
  }

  W <- mofa$W
  F <- mofa$F

  kk <- intersect(rownames(newdata), rownames(W))
  rho1 <- cor(newdata[kk, ], W[kk, ], use = "pairwise")
  rho2 <- cor(t(rho1), t(F), use = "pairwise")
  max.nb <- max.col(rho2)
  predicted <- mofa$samples[max.nb, ]
  rownames(predicted) <- colnames(newdata)
  colnames(predicted) <- paste0("predicted.", colnames(predicted))
  predicted
}


#'
#'
#' @export
mofa.splitByGeneRole <- function(X) {
  G <- playdata::GSETxGENE
  lig <- grep("^LIG", rownames(G), value = TRUE)
  lig.targets <- names(which(Matrix::colSums(G[lig, ]) > 0))
  lig <- unique(gsub("^[a-zA_Z]+:|[ ].*", "", lig))
  kin <- grep("^KINASE_", rownames(G), value = TRUE)
  kin.targets <- names(which(Matrix::colSums(G[kin, ]) > 0))
  kin <- unique(gsub("^[a-zA_Z]+:|[ ].*", "", kin))
  tf <- grep("^TF_", rownames(G), value = TRUE)
  tf.targets <- names(which(Matrix::colSums(G[tf, ]) > 0))
  tf <- unique(gsub("^[a-zA_Z]+:|[ ].*", "", tf))
  mir <- grep("^MIR_", rownames(G), value = TRUE)
  mir.targets <- names(which(Matrix::colSums(G[mir, ]) > 0))
  mir <- unique(gsub("^[a-zA_Z]+:|[ ].*", "", mir))
  mir <- toupper(sub("-", "", gsub("^hsa-|^mmu-|-5p$|-3p$", "", mir)))

  dim(X)
  gg <- rownames(X)[which(!rownames(X) %in% c(lig, tf, kin, mir))]
  data <- list(
    kin = X[rownames(X) %in% kin, ],
    ## tf = X[rownames(X) %in% tf,],
    tf = X[rownames(X) %in% intersect(kin.targets, tf), ],
    ## mir = X[rownames(X) %in% mir,],
    gx = X[gg, , drop = FALSE]
  )
  return(data)
}


#' @export
mofa.exampledata <- function(dataset = "geiger", ntop = 2000,
                             scale.views = TRUE, omx.tissue = NULL) {
  data <- NULL
  if (dataset == "geiger") {
    dir <- "~/Playground/opg-exampledata/metabolomics-kegg"
    counts <- read_counts(file.path(dir, "multiomics-counts.csv"))
    samples <- read_samples(file.path(dir, "multiomics-samples.csv"))
    ## X <- logCPM(counts)
    mindet <- min(counts[counts > 0])
    X <- log2(counts + mindet)
    X <- X + 1e-3 * matrix(rnorm(length(X)), nrow(X), ncol(X))
    ## X <- t(scale(t(X)))
    data <- list(
      px = X[grep("px:", rownames(X)), , drop = FALSE],
      mx = X[grep("mx:", rownames(X)), , drop = FALSE]
    )
    samples <- samples
    contrasts <- as.matrix(samples[, "activated", drop = FALSE])
    colnames(contrasts)[1] <- "act_vs_notact"
    rownames(contrasts) <- rownames(samples)
  }

  if (dataset == "brca") {
    library(mixOmics)
    data(breast.TCGA)
    data <- list(
      mir = t(breast.TCGA$data.train$mirna),
      gx = t(breast.TCGA$data.train$mrna),
      px = t(breast.TCGA$data.train$protein)
    )
    for (i in 1:3) rownames(data[[i]]) <- paste0(names(data)[i], ":", rownames(data[[i]]))
    samples <- data.frame(
      sample = colnames(data[[1]]),
      condition = breast.TCGA$data.train$subtype
    )
    rownames(samples) <- colnames(data[[1]])
    c1 <- c2 <- as.character(samples$condition)
    c1[c1 == "Her2"] <- NA
    c2[c2 == "LumA"] <- NA
    contrasts <- cbind("LumA_vs_Basal" = c1, "Her2_vs_Basal" = c2)
    rownames(contrasts) <- rownames(samples)
  }

  if (dataset == "cll") {
    library(MOFAdata)
    utils::data("CLL_data")
    data <- CLL_data
    lapply(data, dim)
    rownames(data$mRNA) <- convert_probetype("Human", rownames(data$mRNA), "SYMBOL")
    names(data) <- c("dr", "me", "gx", "mu")
    kk <- sort(colnames(data[[1]]))
    for (i in 1:length(data)) {
      data[[i]] <- data[[i]][, kk]
      data[[i]][is.nan(data[[i]])] <- NA
      ## collapse duplicates
      data[[i]] <- rowmean(data[[i]], rownames(data[[i]]))
      data[[i]][is.nan(data[[i]])] <- NA
      rownames(data[[i]]) <- paste0(names(data)[i], ":", rownames(data[[i]]))
    }
    data <- data[c("mu", "me", "gx", "dr")]
    sapply(data, function(x) mean(apply(x, 1, sd, na.rm = TRUE)))
    ## data$Drugs <- 0.01*data$Drugs
    CLL_metadata <- data.table::fread("ftp://ftp.ebi.ac.uk/pub/databases/mofa/cll_vignette/sample_metadata.txt")
    samples <- data.frame(CLL_metadata, row.names = CLL_metadata[["sample"]])
    samples <- samples[kk, ]
    c1 <- c("no", "yes")[1 + samples$IGHV]
    c2 <- c("no", "yes")[1 + samples$trisomy12]
    contrasts <- cbind("IGHV:yes_vs_no" = c1, "trisomy12:yes_vs_no" = c2)
    rownames(contrasts) <- rownames(samples)
  }

  if (dataset == "omx") {
    load("~/Playground/omxp5/omxdata/drugs/ctrpv2-bytissue-omx.rda", verbose = TRUE)
    names(omx)
    data <- omx$level[[1]]$mat
    data <- data[c("mt", "cn", "me", "gx")]
    names(data)

    lapply(data, dim)
    names(data)
    kk <- Reduce(intersect, lapply(data, colnames))
    data <- lapply(data, function(x) x[, kk])
    lapply(data, dim)
    names(data)
    for (i in 1:length(data)) {
      rownames(data[[i]]) <- paste0(names(data)[i], ":", rownames(data[[i]]))
    }
    dim(omx$pheno)
    samples <- omx$pheno[kk, , drop = FALSE]

    if (!is.null(omx.tissue)) {
      sel <- grep(omx.tissue, samples$tissue)
      samples <- samples[sel, ]
      data <- lapply(data, function(d) d[, rownames(samples)])
    }
    samples$tissue <- NULL

    dd <- colnames(samples)
    if (!is.null(omx.tissue)) dd <- sort(grep(omx.tissue, dd, value = 1))
    dd <- grep("BRD", dd, value = 1, invert = TRUE)
    samples <- samples[, dd]
    colnames(samples) <- sub("@", "_", colnames(samples))

    contrasts <- samples
    contrasts <- apply(contrasts, 2, function(x) c("R", "S")[1 + x])
    colnames(contrasts) <- paste0(colnames(contrasts), ":S_vs_R")
  }

  if (dataset == "tf") {
    pgx <- playdata::GEIGER_PGX
    X <- pgx$X
    data <- mofa.splitByGeneRole(pgx$X)
    samples <- pgx$samples
    contrasts <- pgx$contrasts
  }

  if (is.null(data)) stop("could not get dataset")

  ## strip prefix, convert to ASCII
  for (i in 1:length(data)) {
    rownames(data[[i]]) <- sub("^[a-zA_Z]+:", "", rownames(data[[i]]))
    rownames(data[[i]]) <- iconv(rownames(data[[i]]), "latin1", "ASCII", sub = "")
  }

  ## filter out and reduce
  data <- lapply(data, function(d) d[rowMeans(is.na(d)) < 0.5, ])
  data <- lapply(data, function(d) d[rowMeans(d**2, na.rm = TRUE) > 0, ])
  if (ntop > 0) {
    data <- lapply(data, function(x) head(x[order(-apply(x, 1, sd)), ], ntop))
  }

  list(X = data, samples = samples, contrasts = contrasts)
}


#' Normalize a rank vector for each datatype to [-1;1]. This is useful
#' for doing integrated GSEA on multi-omics data.
#'
#' @export
normalize_multirank <- function(rnk) {
  dtypes <- sub(":.*", "", names(rnk))
  nrnk <- rnk
  for (dt in unique(dtypes)) {
    ii <- which(dtypes == dt)
    nrnk[ii] <- ifelse(rnk[ii] < 0, rnk[ii] / max(-rnk[ii], na.rm = TRUE),
      rnk[ii] / max(rnk[ii], na.rm = TRUE)
    )
  }
  nrnk
}

#' Normalize a rank vector for each datatype to [-1;1]. This is useful
#' for doing integrated GSEA on multi-omics data.
#'
#' @export
normalize_multifc <- function(fc, by = c("sd", "mad")[1]) {
  dtypes <- sub(":.*", "", names(fc))
  nfc <- fc
  for (dt in unique(dtypes)) {
    ii <- which(dtypes == dt)
    if (by == "mad") {
      sdx <- mad(fc[ii], na.rm = TRUE)
    } else {
      sdx <- sd(fc[ii], na.rm = TRUE)
    }
    nfc[ii] <- fc[ii] / sdx
  }
  nfc
}


## ======================================================================
## ======================== MULTI-GSEA ==================================
## ======================================================================


#' Compute multi-enrichment for contrasts
#'
#' @export
mgsea.compute_enrichment <- function(F, annot, filter = NULL,
                                     G = NULL, ntop = 1000) {
  ## compute enrichment for all contrasts
  message("computing phenotype enrichment...")

  ## ww.types=filter=G=ntop=NULL
  is.multi <- all(grepl(":", rownames(F)))
  if (inherits(F, "matrix") && is.multi) F <- mofa.split_data(F)
  if (inherits(F, "matrix") && !is.multi) F <- list(gx = F)
  names(F)
  allowed_types <- c("mx", "px", "gx", "me")
  F <- F[names(F) %in% allowed_types]

  ## make sure F is symbol. Enrichment need rownames as symbol
  ## (as in GMT)
  if (!is.null(annot)) {
    F <- mofa.prefix(F)
    F <- lapply(F, function(f) rename_by2(f, annot, "symbol"))
  } else {
    ## should check if by symbol!!!!
  }

  ## determine type
  F.types <- ifelse(grepl("^mx|^metabo", names(F), ignore.case = TRUE),
    "ChEBI", "SYMBOL"
  )
  F.types

  ## If available please use internal pgx$GMT.
  has.mx <- any(F.types == "ChEBI")
  has.px <- any(F.types == "SYMBOL")
  if (is.null(G)) {
    info("[mgsea.compute_enrichment] WARNING. No GMT matrix provided. Not recommended. Please use GMT from pgx object.")
    G <- NULL
    ## add metabolomic gene sets
    if (has.mx) {
      info("[mgsea.compute_enrichment] Adding metabolomics genesets")
      G <- Matrix::t(playdata::MSETxMETABOLITE)
      rownames(G) <- sub("^[a-zA_Z]+:", "", rownames(G)) ## TEMPORARY!!!
    }
    ## add SYMBOL (classic) gene sets
    if (has.px) {
      info("[mgsea.compute_enrichment] Adding transcriptomics/proteomics genesets")
      if (is.null(G)) {
        G <- Matrix::t(playdata::GSETxGENE)
      } else {
        G <- merge_sparse_matrix(G, Matrix::t(playdata::GSETxGENE))
      }
    }
  }
  dim(G)

  ## filter gene sets
  if (!is.null(filter)) {
    sel <- grep(filter, colnames(G))
    G <- G[, sel]
  }

  ## Perform geneset enrichment with fast rank-correlation. We could
  ## do with fGSEA instead but it is much slower.
  rownames(G) <- mofa.strip_prefix(rownames(G))
  F <- mofa.strip_prefix(F)
  wnames <- unlist(lapply(F, rownames))
  pp <- intersect(wnames, rownames(G))
  G <- G[pp, , drop = FALSE]
  G <- G[, Matrix::colSums(G != 0) >= 3, drop = FALSE]
  dim(G)
  gg <- rownames(G)
  gsea <- list()
  i <- 1
  for (i in 1:length(F)) {
    w1 <- F[[i]]
    nc <- sum(rownames(w1) %in% gg)
    if (nc > 0) {
      pp <- intersect(rownames(w1), rownames(G))
      f1 <- gset.rankcor(w1[pp, , drop = FALSE], G[pp, ], compute.p = TRUE)
      dt <- names(F)[i]
      if (!all(is.na(f1$rho)) && !all(f1$rho == 0)) {
        gsea[[dt]] <- f1
      }
    }
  }
  names(gsea)

  ## extract rho, pval and qval
  rho <- lapply(gsea, function(x) x$rho)
  p.value <- lapply(gsea, function(x) x$p.value)

  ## correct missing values
  for (i in 1:length(gsea)) {
    rho[[i]][which(is.na(rho[[i]]))] <- 0
    p.value[[i]][which(is.na(p.value[[i]]))] <- 0.9999
  }

  ## As multi-omics score we take the absolute geometric average
  nt <- length(rho)
  multi.rho <- log(pmax(abs(rho[[1]]), 0.01))
  if (nt > 1) {
    for (i in 2:nt) {
      multi.rho <- multi.rho + log(pmax(abs(rho[[i]]), 0.01))
    }
  }
  multi.rho <- exp(multi.rho / nt)

  ## As multi-omics p.value we use Stouffer's method
  P <- sapply(p.value, as.vector)
  P[is.na(P)] <- 0.999
  P <- pmin(pmax(P, 1e-99), 0.999999)
  if (ncol(P) > 1) {
    multi.p <- apply(P, 1, function(x) metap::sumz(x)$p) ## slow....
  } else {
    multi.p <- P[, 1]
  }
  nr <- nrow(p.value[[1]])
  nc <- ncol(p.value[[1]])
  dimnames <- dimnames(p.value[[1]])
  multi.p <- matrix(multi.p, nrow = nr, ncol = nc, dimnames = dimnames)
  multi.q <- apply(multi.p, 2, p.adjust)
  dim(multi.p)

  ## determine is sign/direction is all same
  S <- sapply(rho, as.vector)
  multi.sign <- matrix(1, nr, nc)
  if (ncol(S) > 1) {
    multi.sign <- matrixStats::rowSds(sign(S), na.rm = TRUE) < 1e-8
    multi.sign <- matrix(multi.sign,
      nrow = nr, ncol = nc,
      dimnames = dimnames
    )
  }

  ## define a integrated score
  multi.score <- multi.rho * -log(multi.p) * multi.sign

  ## add some stats
  gsets <- rownames(rho[[1]])
  dt <- sapply(F, function(f) 1 * (rownames(G) %in% rownames(f)))
  dt <- colnames(dt)[max.col(dt)]
  sizes <- tapply(1:nrow(G), dt, function(i) Matrix::colSums(G[i, , drop = FALSE] != 0))
  sizes <- do.call(cbind, sizes)

  ## Create nice dataframe
  enr <- list()
  for (i in 1:ncol(multi.rho)) {
    df.rho <- sapply(rho, function(x) x[, i])
    df.pval <- sapply(p.value, function(x) x[, i])
    df.multi <- data.frame(
      score = multi.score[, i],
      rho = multi.rho[, i],
      sign = multi.sign[, i],
      p = multi.p[, i],
      q = multi.q[, i]
    )
    colnames(df.multi) <- paste0("multi.", colnames(df.multi))
    colnames(df.rho) <- paste0("rho.", colnames(df.rho))
    colnames(df.pval) <- paste0("pval.", colnames(df.pval))
    colnames(sizes) <- paste0("num.", colnames(sizes))
    df <- data.frame(
      df.multi,
      sizes,
      df.rho,
      df.pval
    )
    # df <- df[!is.na(df$multi.rho),,drop=FALSE]
    ct <- colnames(multi.rho)[i]
    enr[[ct]] <- df
  }

  if (!is.null(ntop) && ntop > 0) {
    i <- 1
    for (i in 1:length(enr)) {
      mscore <- enr[[i]]$multi.score
      rho.col <- grep("^rho[.]", colnames(enr[[i]]))[1]
      msign <- sign(enr[[i]][, rho.col])
      ii <- order(-mscore * msign)
      enr[[i]] <- head(enr[[i]][ii, ], ntop)
    }
  }

  return(enr)
}


#'
#' @export
mgsea.plot_scatter <- function(gsea, type1, type2, size.par = "p",
                               main = NULL, hilight = NULL) {
  ## gsea=res$gsea;type1="px";type2="mx"

  ## using pre-computed GSEA
  R <- gsea
  R <- R[, grep("^rho", colnames(R)), drop = FALSE]
  colnames(R) <- sub("^rho.", "", colnames(R))
  s1 <- R[, type1]
  s2 <- R[, type2]
  if (size.par == "q") {
    qcomb <- gsea$multi.q ## meta q-value
  } else {
    qcomb <- gsea$multi.p ## meta q-value
  }

  s1 <- s1 + 2e-3 * rnorm(length(s1))
  s2 <- s2 + 2e-3 * rnorm(length(s2))

  ## point size
  cex1 <- 0.1 + (1 - qcomb)**2
  col0 <- ifelse(is.null(hilight), "black", "grey60")
  plot(s1, s2,
    cex = 1.4 * cex1, col = col0,
    xlab = paste(type1, "enrichment (rho)"),
    ylab = paste(type2, "enrichment (rho)")
  )
  if (is.null(title)) {
    title <- paste("multiGSEA:", toupper(type1), "vs.", toupper(type2))
  } else {
    title <- main
  }
  title(title, cex.main = 1.3)

  abline(h = 0, v = 0, lty = 3)
  kk <- rownames(gsea)
  if (!is.null(hilight) && length(hilight) > 0) {
    if (length(intersect(hilight, kk))) {
      sel <- intersect(hilight, kk)
    } else {
      sel <- grep(hilight, kk, value = TRUE)
    }
    if (length(sel)) {
      sel <- match(sel, kk)
      points(s1[sel], s2[sel],
        pch = 20, col = "red2",
        lwd = 2, cex = 1.5 * cex1[sel]
      )
    }
  }
}


#' @export
mgsea.plot_barplot <- function(gsea,
                               type = c("barplot", "lollipop")[1],
                               remove.dup = TRUE,
                               ntop = 20,
                               filter = NULL,
                               select = NULL,
                               strip.names = FALSE,
                               par = NULL,
                               title = "multi-omics enrichment") {
  ## k=1;ntop=20
  S <- gsea
  if (nrow(S) == 0) {
    return(NULL)
  }

  if (!is.null(select) && length(select) > 0) {
    S <- S[select, , drop = FALSE]
  }

  if (!is.null(filter)) {
    S <- S[grep(filter, rownames(S)), , drop = FALSE]
  }

  if (remove.dup) {
    sname <- gsub("^[a-zA_Z]+:| \\(.*", "", rownames(S))
    S <- S[which(!duplicated(sname)), , drop = FALSE]
  }

  if ("multi.score" %in% names(S)) {
    S <- S[order(-abs(S$multi.score)), , drop = FALSE]
    topS <- abs(S[, grep("^rho", colnames(S)), drop = FALSE])
  } else if ("multi" %in% names(S)) {
    S <- S[order(-abs(S$multi$score)), , drop = FALSE]
    topS <- abs(S[["rho"]])
  }
  topS <- head(topS, ntop)
  topS <- topS[nrow(topS):1, , drop = FALSE]
  plot <- NULL
  if (strip.names) {
    rownames(topS) <- gsub(" \\(.*", "", rownames(topS))
  }

  if (type == "barplot") {
    if (is.null(par)) {
      par(mfrow = c(1, 2), mar = c(4, 10, 4, 2))
      plot.new()
    } else {
      par(par)
    }
    barplot(t(topS), horiz = TRUE, las = 1, xlab = "cumulative abs.enrichment (rho)")
    cc <- grey.colors(ncol(topS))
    legend("bottomright",
      legend = colnames(topS),
      fill = cc, inset = c(0.0, 0.02), bg = "white",
      y.intersp = 0.8, cex = 0.85
    )
    title(title, outer = TRUE, cex.main = 1.4, line = -1.2)
  }

  if (type == "lollipop") {
    ## need double lolliplot lines for multi-omics
    top.nes <- rowSums(topS)
    names(top.nes) <- rownames(topS)
    plot <- ggLollipopPlot(top.nes, xlab = "enrichment (NES)")
  }
  plot
}


#' Plots from results from compute.mofa()
#'
#' @export
mofa.plot_biplot <- function(r, pheno, nfeat = 5, comp = c(1, 2),
                             arrow = TRUE, cex = 1) {
  names(r)
  pheno <- factor(pheno)
  col <- 1 + as.integer(pheno)
  ff <- scale(r$F[, comp], center = FALSE)
  plot(ff, col = col, pch = 19, cex = 0.8 * cex)
  tt <- paste0(pheno, "_", rownames(r$F))
  tt <- rownames(r$F)
  text(ff, tt, col = col, cex = 0.8 * cex, pos = 3)
  legend("topright",
    legend = levels(pheno), fill = 2:7,
    cex = 0.8, y.intersp = 0.7
  )

  if (arrow) {
    R <- cor(t(r$W), t(r$Z))
    dim(R)
    head(R)
    top <- tail(names(sort(abs(R[, 2]))), nfeat)
    R[top, ]
    ww <- r$W[, comp]
    ww <- scale(ww, center = FALSE)
    ww[top, ]
    aa <- 1.0 * abs(ww[top, ])**0.5 * sign(ww[top, ])
    arrows(0, 0, aa[, 1], aa[, 2], length = 0.1, lwd = 1.6)
    text(1.1 * aa[, 1], 1.1 * aa[, 2], labels = rownames(aa), cex = 0.8)
  }
}


## ======================================================================
## ========================== MIXOMICS ==================================
## ======================================================================

## mixOmics related compute and plot functions

#' @export
mixomics.compute_diablo <- function(xdata, y, ncomp = 5, nfeat = 20,
                                    method = c(
                                      "diablo", "rgcca", "sgcca",
                                      "sgccda", "rgcca.rgcca",
                                      "rgcca.rgccda", "rgcca.mcoa"
                                    )) {
  xx <- xdata
  xx <- mofa.prefix(xx)
  xx <- lapply(xx, t)
  method <- tolower(method[1])

  ## number of features must be smaller than smallest dim
  minx <- min(sapply(xx, function(x) min(dim(x))))
  nfeat <- min(nfeat, minx)

  ## full design as default
  design <- (1 - diag(length(xx)))
  design

  dbg("[mixomics.compute_diablo] nfeat = ", nfeat)
  keepX <- lapply(names(xx), function(dt) rep(nfeat, ncomp))
  names(keepX) <- names(xx)

  if (method %in% c("splsda", "diablo")) {
    res <- mixOmics::block.splsda(
      X = xx,
      Y = y,
      keepX = keepX,
      ncomp = ncomp,
      design = design
    )
  } else if (method == "rgcca") {
    res <- mixOmics::wrapper.rgcca(
      X = xx,
      design = design,
      ncomp = ncomp,
      keepX = keepX
    )
  } else if (method == "rgcca.rgcca") {
    res <- RGCCA::rgcca(
      blocks = xx,
      method = "rgcca",
      connection = design, ## full connect
      tau = rep(1, length(xx)),
      ncomp = rep(ncomp, length(xx)),
      superblock = TRUE,
      verbose = FALSE
    )
    ## need res$loadings, res$X
    for (i in 1:length(res$a)) colnames(res$a[[i]]) <- paste0("comp", 1:ncomp)
    res$loadings <- res$a[names(xx)]
  } else if (method == "rgcca.rgccda") {
    xx2 <- c(response = data.frame(y), xx)
    res <- RGCCA::rgcca(
      blocks = xx2,
      response = 1,
      method = "rgcca",
      connection = design, ## full connect
      tau = rep(1, length(xx2)),
      ncomp = ncomp,
      superblock = TRUE,
      verbose = FALSE
    )
    ## need res$loadings, res$X
    for (i in 1:length(res$a)) colnames(res$a[[i]]) <- paste0("comp", 1:ncomp)
    res$loadings <- res$a[names(xx)]
  } else if (method == "rgcca.mcoa") {
    res <- RGCCA::rgcca(
      blocks = xx,
      method = "mcoa",
      connection = design, ## full connect
      tau = rep(1, length(xx)),
      ncomp = ncomp,
      superblock = TRUE,
      verbose = FALSE
    )
    ## need res$loadings, res$X
    for (i in 1:length(res$a)) colnames(res$a[[i]]) <- paste0("comp", 1:ncomp)
    res$loadings <- res$a[names(xx)]
  } else if (method == "sgcca") {
    res <- mixOmics::wrapper.sgcca(
      X = xx,
      design = design,
      ncomp = ncomp,
      keepX = keepX
    )
  } else if (method == "sgccda") {
    res <- mixOmics::wrapper.sgccda(
      X = xx,
      Y = y,
      design = design,
      ncomp = ncomp,
      keepX = keepX
    )
  } else {
    stop("[mixomics.compute_diablo] FATAL: invalid method", method)
  }

  ## extract weight/loading matrix
  dt <- names(xx)
  ww <- res$loadings[dt]
  W <- do.call(rbind, ww)

  ## extract factor matrix. NEED RETHINK!!!! instead of cor(X,W) we
  ## should actually solve W.F=X for F.
  ## X <- do.call(rbind, lapply(res$X,t))
  X <- t(do.call(cbind, xx))
  F <- cor(X, W) ## factors on phenotype???

  list(
    res = res,
    scores = F,
    loadings = W
  )
}

#' @export
mixomics.plot <- function(model, style, comp = 1, cutoff = 0.7, ...) {
  dt <- names(model$X)
  nt <- length(dt)
  block.colors <- rainbow(nt)

  if (style == "correlation") {
    mixOmics::plotDiablo(model, ncomp = 1)
  }

  if (style == "indiv") {
    mixOmics::plotIndiv(model,
      ind.names = FALSE, legend = TRUE,
      title = "DIABLO Sample Plots"
    )
  }

  if (style == "arrow") {
    mixOmics::plotArrow(model,
      ind.names = FALSE, legend = TRUE,
      title = "DIABLO"
    )
  }

  if (style == "var") {
    mixOmics::plotVar(model,
      var.names = FALSE,
      style = "graphics", legend = TRUE,
      pch = c(16, 17, 15), cex = c(1, 1, 1) * 1.2,
      col = block.colors
    )
  }

  if (style == "circos") {
    mixOmics::circosPlot(model,
      cutoff = cutoff, line = TRUE,
      color.blocks = block.colors,
      color.cor = c("black", "red"),
      size.labels = 1.5
    )
  }

  if (style == "network") {
    mixOmics::network(model,
      blocks = dt,
      color.node = block.colors,
      size.node = 0.04,
      cutoff = cutoff
    )
  }

  if (style == "loadings") {
    plotLoadings(model, comp = comp, contrib = "max", method = "median")
  }

  if (style == "cim") {
    cimDiablo(model,
      mar = c(8, 20), transpose = TRUE,
      legend.position = "bottomright", size.legend = 1
    )
  }
}

#' @export
mixomics.predict <- function(model, newdata) {
  newdata <- mofa.scale_views(newdata)
  newdata <- lapply(mofa.prefix(newdata), t)
  for (k in 1:length(newdata)) {
    dtype <- names(newdata)[k]
    if (dtype %in% names(model$X)) {
      mat <- newdata[[k]]
      jj <- match(colnames(model$X[[dtype]]), colnames(mat))
      newdata[[k]] <- mat[, jj]
    } else {
      newdata[[k]] <- NULL
    }
  }
  names(newdata)
  pred <- mixOmics:::predict.block.spls(model, newdata = newdata)
  predicted <- pred$WeightedVote$centroids.dist[, 2]
  predicted
}


## ===================================================================
## ======================= SNF =======================================
## ===================================================================

#' @export
snf.cluster <- function(xx, pheno = NULL, plot = TRUE) {
  has.na <- any(sapply(xx, function(x) sum(is.na(x)) > 0))
  has.inf <- any(sapply(xx, function(x) sum(is.infinite(x)) > 0))
  has.missing <- has.na || has.inf
  message("clusterSNF: has.missing = ", has.missing)
  if (has.missing) {
    xx <- lapply(xx, function(x) svdImpute2(x, infinite.na = TRUE))
  }

  Data <- lapply(xx, t)
  Dist <- list()
  dt <- names(xx)[1]
  for (dt in names(Data)) {
    x <- Data[[dt]]
    x <- SNFtool::standardNormalization(x) ## cannot handle NA
    Dist[[dt]] <- SNFtool::dist2(x, x)^(1 / 2)
  }

  ## next, construct similarity graphs
  ## First, set all the parameters:
  ## K = max(min(ncol(xx[[1]])/4, 15),2);
  K <- max(min(ncol(xx[[1]]), 10), 2) # number of neighbors, usually (10~30)
  alpha <- 0.5 # hyperparameter, usually (0.3~0.8)

  message("clusterSNF: K = ", K)
  message("clusterSNF: alpha = ", alpha)
  message("clusterSNF: length(Dist) = ", length(Dist))
  message("clusterSNF: names(Dist) = ", names(Dist))
  message("clusterSNF: dim(Dist[[1]]) = ", paste(dim(Dist[[1]]), collapse = "x"))

  Wlist <- lapply(Dist, function(d) SNFtool::affinityMatrix(d, K = K, alpha))

  ## next, we fuse all the graphs
  ## then the overall matrix can be computed by similarity network fusion(SNF):
  W <- SNFtool::SNF(Wlist, K = K, t = 100)

  ## With this unified graph W of size n x n,
  ## you can do either spectral clustering or Kernel NMF.
  ## If you need help with further clustering, please let us know.

  Data2 <- c(Dist, list(SNF = W))
  lapply(Data2, dim)
  ## posx <- lapply( Data2, function(x) Rtsne::Rtsne(t(x), is_distance=TRUE, perplexity=15)$Y)
  k <- max(min(ncol(xx[[1]]) / 4, 15), 1) # number of neighbors
  posx <- lapply(Data2, function(x) {
    Rtsne::Rtsne(t(x), perplexity = k, check_duplicates = FALSE)$Y
  })
  ## posx <- lapply( Data2, function(x) uwot::umap(t(x), n_neighbors=15))
  for (i in 1:length(posx)) rownames(posx[[i]]) <- colnames(xx[[1]])
  rownames(W) <- colnames(W) <- rownames(posx[[1]])

  ## -------------- graph & louvain ----------------
  gr <- igraph::graph_from_adjacency_matrix(
    W,
    mode = "undirected", diag = FALSE, weighted = TRUE
  )
  cl <- igraph::cluster_louvain(gr)$membership
  cl <- paste0("SNF", cl)
  table(cl)

  list(
    affinityMatrix = Wlist,
    W = W,
    posx = posx,
    pheno = pheno,
    graph = gr,
    cluster = cl
  )
}


#' @export
snf.plot_affinity <- function(snf, k = 1, par = TRUE) {
  i <- 1
  mat <- c(snf$affinityMatrix, list(SNF = snf$W))
  n <- ceiling(sqrt(length(mat)))
  if (par) par(mfrow = c(n, n), mar = c(6, 1, 2, 7))
  for (i in 1:length(mat)) {
    x <- mat[[i]]
    if (k != 1) x <- sign(x) * x^k
    diag(x) <- 0
    if (nrow(x) <= 20) {
      gx.imagemap(x)
    } else {
      rownames(x) <- NULL
      colnames(x) <- NULL
      gx.imagemap(x)
      mtext("samples", side = 1, line = 0.2)
      mtext("samples", side = 4, line = 0.2)
    }
    title(names(mat)[i], cex.main = 1.6)
  }
}

#' @export
snf.plot_clustering <- function(snf, y = NULL, par = TRUE, ...) {
  i <- 1
  n <- ceiling(sqrt(length(snf$posx)))
  if (par) par(mfrow = c(n, n), mar = c(6, 1, 2, 7))
  cc <- "grey30"
  if (!is.null(y)) cc <- factor(y)
  for (i in 1:length(snf$posx)) {
    plot(snf$posx[[i]],
      col = cc, pch = 19, cex = 0.8,
      xlab = "TSNE-x", ylab = "TSNE-y", ...
    )
    title(names(snf$posx)[i], cex.main = 1.6)
  }
}

#' @export
snf.plot_graph <- function(snf, min.rho = 0.5, plot = TRUE,
                           val = NULL, label = NULL,
                           vlabcex = 1, q.edge = 0.85,
                           ...) {
  A <- snf$affinityMatrix
  W <- snf$W
  diag(W) <- 0
  W <- W / max(W)
  for (i in 1:length(A)) {
    diag(A[[i]]) <- 0
    A[[i]] <- A[[i]] / max(A[[i]])
    ##  A[[i]] <- A[[i]] * (A[[i]] > 0.5)
  }

  rownames(W) <- colnames(W) <- rownames(snf$posx[[1]])
  if (!is.null(label)) {
    names(label) <- rownames(W)
  }
  ee <- which(W > min.rho, arr.ind = TRUE)
  # ee <- which(W > 0, arr.ind=TRUE)
  ee <- ee[which(ee[, 2] > ee[, 1]), ]
  af <- unlist(lapply(A, function(a) a[ee]))
  af <- af / max(af, na.rm = TRUE)
  tp <- as.vector(sapply(names(A), function(a) rep(a, nrow(ee))))
  df <- data.frame(ee, type = tp, affinity = af)
  head(df)
  df[, 1] <- rownames(W)[df[, 1]]
  df[, 2] <- rownames(W)[df[, 2]]
  gr <- igraph::graph_from_edgelist(
    as.matrix(df[, 1:2]),
    directed = FALSE
  )

  colors <- paste0(rainbow(length(A)), "55")
  colors <- head(rep(paste0(palette()[-1], "66"), 99), length(A))

  igraph::V(gr)$name
  igraph::E(gr)$type <- df$type
  igraph::E(gr)$weight <- (1 * df$affinity)**1.1
  igraph::E(gr)$color <- colors[factor(df$type)]
  igraph::V(gr)$label <- igraph::V(gr)$name
  if (!is.null(label)) {
    igraph::V(gr)$label <- label[igraph::V(gr)$name]
  }

  length(igraph::E(gr))
  ewt <- igraph::E(gr)$weight
  summary(ewt)
  qq <- quantile(ewt, probs = q.edge)[1] ## edge threshold
  gr <- igraph::delete_edges(gr, igraph::E(gr)[which(ewt < qq)])

  vcolors <- rev(rainbow(8))
  if (is.null(val)) {
    val <- factor(snf$cluster)
    names(val) <- igraph::V(snf$graph)$name
    val <- val[igraph::V(gr)$name]
  }
  igraph::V(gr)$color <- vcolors[factor(val)]

  if (plot) {
    par(mfrow = c(1, 1), mar = c(1, 1, 1, 1))
    dtypes <- unique(factor(df$type))
    plot(gr,
      vertex.size = 0,
      vertex.label = igraph::V(gr)$label,
      vertex.label.family = "sans",
      vertex.label.cex = 1.15 * vlabcex,
      edge.width = 6 * igraph::E(gr)$weight,
      edge.color = igraph::E(gr)$color,
      ...
    )
    legend("bottomright",
      legend = dtypes, fill = colors,
      cex = 0.9, y.intersp = 0.85
    )
  } else {
    return(gr)
  }
}

#' @export
snf.heatmap <- function(snf, X, samples, nmax = 50, do.split = TRUE, legend = TRUE) {
  rX <- X[order(-apply(X, 1, sd)), ]
  dt <- sub(":.*", "", rownames(X))
  n <- nmax / length(unique(dt))
  sel <- tapply(rownames(X), dt, function(a) head(a, n))
  sel <- unlist(sel)
  dt <- NULL
  if (do.split) dt <- sub(":.*", "", sel)
  gx.splitmap(rX[sel, ],
    splitx = snf$cluster,
    split = dt,
    nmax = 160, cexRow = 0.8,
    col.annot = samples,
    dist.method = "euclidean",
    show_legend = legend,
    softmax = TRUE, show_key = FALSE
  )
}

## ======================================================================
## ====================== LASAGNA =======================================
## ======================================================================


#' @export
plot_lasagna <- function(posx, vars = NULL, num_edges = 20) {
  for (i in 1:length(posx)) {
    a <- names(posx)[i]
    if (!all(grepl(a, names(vars[[i]])))) {
      rownames(posx[[i]]) <- paste0(a, ":", rownames(posx[[i]]))
      if (!is.null(vars) && is.list(vars)) {
        names(vars[[i]]) <- sub(paste0(a, ":"), "", names(vars[[i]]))
        names(vars[[i]]) <- paste0(a, ":", names(vars[[i]]))
      }
    }
  }

  ## Sample planes across datatyps
  colx <- "black"

  ## feature maps across datatypes
  x <- data.frame()
  for (i in names(posx)) {
    x <- rbind(x, data.frame(posx[[i]], type = i))
  }
  x$type <- factor(x$type, levels = names(posx))

  if (!is.null(vars)) {
    colorsx <- gplots::colorpanel(255, low = "blue3", mid = "grey80", high = "red3")
    if (!is.list(vars) && is.vector(vars)) {
      vars <- lapply(posx, function(p) vars[rownames(p)])
    }
    vars <- lapply(vars, function(x) (x / max(abs(x))))
    colx <- unlist(lapply(vars, function(v) colorsx[128 + ceiling(v * 127)]))
    dt <- names(posx)
    dt
    i <- 1
    ee <- c()
    for (i in 1:(length(dt) - 1)) {
      v1 <- vars[[i]][rownames(posx[[i]])]
      v2 <- vars[[i + 1]][rownames(posx[[i + 1]])]
      ii <- head(names(sort(-v1)), num_edges)
      jj <- head(names(sort(-v2)), num_edges)
      ee1 <- cbind(ii, jj)
      ee <- rbind(ee, ee1)
      ii <- head(names(sort(v1)), num_edges)
      jj <- head(names(sort(v2)), num_edges)
      ee2 <- cbind(ii, jj)
      ee <- rbind(ee, ee2)
    }
  }

  segment_mat <- cbind(
    match(ee[, 1], rownames(x)),
    match(ee[, 2], rownames(x))
  )
  segment_mat
  nseg <- nrow(segment_mat)

  grimon::grimon(
    x = x,
    label = names(posx),
    col = colx,
    format = "long",
    segment_mat = segment_mat[, ],
    segment_col = colx[segment_mat[, 1]],
    point_size = 4,
    z_interval = 2,
    optimize_coordinates = TRUE,
    maxiter = 1e3,
    score_function = "angle",
    # return_coordinates = TRUE,
    # plot_2d_panels = FALSE,
    segment_alpha = 0.8
  )
}

#'
#'
#' @export
plotly_lasagna <- function(posx, vars = NULL, num_edges = 20) {
  for (i in 1:length(posx)) {
    a <- names(posx)[i]
    if (!all(grepl(a, names(vars[[i]])))) {
      rownames(posx[[i]]) <- paste0(a, ":", rownames(posx[[i]]))
      if (!is.null(vars)) {
        names(vars[[i]]) <- paste0(a, ":", names(vars[[i]]))
      }
    }
  }

  ## scale each layer
  posx <- lapply(posx, uscale)

  ## feature maps across datatypes
  x <- data.frame()
  for (i in names(posx)) {
    x <- rbind(x, data.frame(posx[[i]], type = i))
  }
  x$type <- factor(x$type, levels = names(posx))
  colnames(x) <- c("x", "y", "z")

  znames <- c(
    "mx" = "Metabolomics",
    "gx" = "Transcriptomics",
    "mir" = "micro-RNA",
    "px" = "Proteomics"
  )

  plotlyLasagna(x, znames = znames)
}


#'
#' @export
lasagna.create_model <- function(data, pheno = "pheno", ntop = 1000, nc = -1,
                                 annot = NULL, use.gmt = TRUE, use.graphite = TRUE) {
  names(data)
  if (pheno == "pheno") {
    Y <- expandPhenoMatrix(data$samples, drop.ref = FALSE)
  } else {
    Y <- makeContrastsFromLabelMatrix(data$contrasts)
    revY <- -Y
    colnames(revY) <- reverse.AvsB(colnames(Y))
    Y <- cbind(Y, revY)
  }
  data$X$PHENO <- t(Y)
  names(data)

  xx <- data$X
  if (!is.null(ntop) && ntop > 0) {
    xx <- lapply(xx, function(x) head(x[order(-apply(x, 1, sd)), , drop = FALSE], ntop))
  }
  xx <- mofa.prefix(xx)
  X <- do.call(rbind, xx)
  remove(xx)

  ## add SOURCE/SINK
  X <- rbind(X, SOURCE = 1, SINK = 1)
  suppressWarnings(R <- cor(t(X), use = "pairwise"))
  R[is.na(R)] <- 0.01

  ## mask for metabolics PPI
  if (use.graphite) {
    xtypes <- setdiff(names(data$X), "PHENO")
    xtypes
    has.mx <- ("mx" %in% xtypes)
    has.px <- any(c("gx", "px") %in% xtypes)
    if (has.mx && has.px) {
      message("using GRAPHITE PPI structure")
      load("~/Playground/public-db/pathbank.org/GRAPHITE_PPI.rda", verbose = TRUE)
      head(GRAPHITE_PPI)
      GRAPHITE_PPI[, 1] <- paste0("px:", GRAPHITE_PPI[, 1])
      GRAPHITE_PPI[, 2] <- paste0("px:", GRAPHITE_PPI[, 2])
      GRAPHITE_PPI[, 1] <- sub("px:CHEBI:", "mx:", GRAPHITE_PPI[, 1])
      GRAPHITE_PPI[, 2] <- sub("px:CHEBI:", "mx:", GRAPHITE_PPI[, 2])
      gr <- igraph::graph_from_edgelist(as.matrix(GRAPHITE_PPI[, 1:2]), directed = "FALSE")
      A <- igraph::as_adjacency_matrix(igraph::simplify(gr))
      dt <- sub(":.*", "", rownames(R))
      i1 <- which(dt %in% c("gx", "px"))
      i2 <- which(dt == "mx")
      a1 <- match(rownames(R)[i1], rownames(A))
      a2 <- match(rownames(R)[i2], rownames(A))
      i1 <- i1[!is.na(a1)]
      a1 <- a1[!is.na(a1)]
      i2 <- i2[!is.na(a2)]
      a2 <- a2[!is.na(a2)]
      R[i1, i2] <- R[i1, i2] * as.matrix(A[a1, a2])
      R[i2, i1] <- R[i2, i1] * as.matrix(A[a2, a1])
    }
  }

  ## mask for GSETS/pathways ???
  if (use.gmt) {

  }

  ## mask for layer interaction
  dt <- sub(":.*", "", rownames(R))
  table(dt)
  names(data$X)
  layers <- c("SOURCE", names(data$X), "SINK")
  layer_mask <- as.matrix(R) * 0
  for (i in 1:(length(layers) - 1)) {
    ii <- which(dt == layers[i])
    jj <- which(dt == layers[i + 1])
    layer_mask[ii, jj] <- 1
    layer_mask[jj, ii] <- 1
  }
  R <- R * layer_mask

  ## Reduce connections
  if (!is.null(nc) && nc > 0) {
    message(paste("reducing edges to maximum", nc, "connections"))
    ## NEED CHECK!!! SOMETHING WRONG.
    xtypes <- setdiff(names(data$X), "PHENO")
    reduce_mask <- matrix(1, nrow(R), ncol(R))
    i <- 1
    for (i in 1:(length(xtypes) - 1)) {
      ii <- which(dt == xtypes[i])
      jj <- which(dt == xtypes[i + 1])
      R1 <- R[ii, jj, drop = FALSE]
      rii <- apply(abs(R1), 1, function(r) tail(sort(r), nc)[1])
      rjj <- apply(abs(R1), 2, function(r) tail(sort(r), nc)[1])
      rr <- abs(R1) >= rii | t(t(abs(R1)) >= rjj)
      reduce_mask[ii, jj] <- rr
      reduce_mask[jj, ii] <- t(rr)
    }
    R <- R * reduce_mask
  }

  ## create graph
  gr <- igraph::graph_from_adjacency_matrix(
    R,
    diag = FALSE, weighted = TRUE, mode = "undirected"
  )

  ii <- grep("^mx:", igraph::V(gr)$name)
  if (length(ii) && !is.null(annot)) {
    message(paste("translating metabolite ID to names"))
    mx.id <- sub("mx:", "", igraph::V(gr)$name[ii])
    mm <- match(igraph::V(gr)$name[ii], rownames(annot))
    mx.annot <- annot[mm, ]
    mx.id <- mx.annot$symbol
    mx.names <- paste0("mx:", mx.annot$gene_title, " (", mx.id, ")")
    igraph::V(gr)$name[ii] <- mx.names
    rownames(X)[ii] <- mx.names
  }

  ## R1 <- igraph::as_adjacency_matrix(gr, attr='weight')
  list(
    graph = gr,
    X = X,
    Y = Y,
    layers = layers
  )
}


#'
#' @export
lasagna.solve_SP <- function(obj, pheno, rm.neg = TRUE, vtop = 200) {
  ## condition with phenotype correlation
  PHENO <- paste0("PHENO:", pheno)
  if (is.character(pheno) && PHENO[1] %in% rownames(obj$X)) {
    pheno <- obj$X[PHENO, ]
  }
  suppressWarnings(rho <- cor(t(obj$X), pheno)[, 1])
  rho[is.na(rho)] <- 0.1 ## ???
  names(rho) <- rownames(obj$X)

  tail(sort(rho))
  head(sort(rho))
  F <- outer(rho, rho)
  F <- sqrt(abs(F))
  dim(F)

  ## weight edges with foldchange on nodes
  gr <- obj$graph
  R <- igraph::as_adjacency_matrix(gr, attr = "weight")
  F <- F[rownames(R), colnames(R)]
  R <- R * F
  gr <- igraph::graph_from_adjacency_matrix(
    R,
    diag = FALSE, weighted = TRUE, mode = "undirected"
  )
  gr <- igraph::simplify(gr)

  if (rm.neg) {
    message("removing negative edges")
    gr <- igraph::delete_edges(gr, edges = which(igraph::E(gr)$weight <= 0))
  }

  ## compute all(?) shortest paths
  wt <- igraph::E(gr)$weight
  max.wt <- max(wt, na.rm = TRUE)
  ## wt <- (max.wt - wt)**1
  7 ## wt <- (max.wt / pmax(wt, 1e-2) - 1)**1
  wt <- -log(wt + 1e-4)

  ## resctrict search to top logFC nodes
  v.nodes <- igraph::V(gr)
  if (vtop > 0) {
    v.names <- igraph::V(gr)$name
    v.types <- sub(":.*", "", v.names)
    vtop.nodes <- tapply(v.names, v.types, function(v) {
      head(names(sort(-abs(rho[v]))), vtop)
    })
    v.nodes <- unlist(vtop.nodes)
  }

  suppressWarnings({
    s1 <- igraph::shortest_paths(gr,
      from = "SOURCE", to = v.nodes, weights = wt,
      output = "both"
    )
    s2 <- igraph::shortest_paths(gr,
      from = "SINK", to = v.nodes, weights = wt,
      output = "both"
    )
  })

  p1 <- sapply(s1$epath, function(e) sum(wt[e]))
  p2 <- sapply(s2$epath, function(e) sum(wt[e]))
  sp.score <- p1 + p2
  names(sp.score) <- v.nodes

  v1 <- lapply(s1$vpath, function(i) igraph::V(gr)$name[i])
  v2 <- lapply(s2$vpath, function(i) rev(igraph::V(gr)$name[i]))
  vv <- mapply(union, v1, v2, SIMPLIFY = FALSE)
  names(vv) <- v.nodes

  ## check paths: this correct the shortest path
  splen <- sapply(vv, length)
  table(splen)
  splen0 <- min(splen[splen > 0])
  ii <- which(splen > splen0)
  if (length(ii)) {
    rmdupv <- function(v) v[!duplicated(sub(":.*", "", v))]
    vv[ii] <- lapply(vv[ii], function(v) rmdupv(v))
  }
  names(sp.score) <- v.nodes

  vv <- lapply(vv, function(v) setdiff(v, c("SOURCE", "SINK")))
  vv <- vv[!duplicated(names(vv))]
  vv <- vv[sapply(vv, length) > 0]
  spath <- sapply(vv, paste, collapse = "->")
  V <- data.frame(do.call(rbind, vv))

  ##  V$path <- spath
  rownames(V) <- make_unique(names(vv))
  vtype <- sub(":.*", "", V[1, ])
  colnames(V) <- vtype
  head(V)

  S <- apply(V, 2, function(v) rho[v])
  rownames(S) <- rownames(V)
  sp.score <- sp.score[rownames(S)]

  min.sp <- min(sp.score, na.rm = TRUE)
  sp.score <- exp(-2 * sp.score / min.sp)

  rho <- rho[rownames(S)]
  sp.score <- sp.score * sign(rho)

  sp <- data.frame(
    score = sp.score,
    rho = S,
    path = spath
  )
  rownames(sp) <- rownames(V)

  sp <- sp[grep("SOURCE|SINK", rownames(sp), invert = TRUE), , drop = FALSE]
  sp <- sp[order(-sp$score), ]
  sp
}

#'
#' @export
lasagna.plot_SP <- function(sp, ntop = 200, hilight = NULL, labcex = 1,
                            colorby = "pheno", plotlib = "ggplot") {
  if (!is.null(ntop) && ntop > 0) {
    sp <- sp[order(-sp$score), ]
    sp <- head(sp, ntop)
  }
  fig <- NULL

  if (plotlib == "plotly") {
    dimensions <- list()
    this.layer <- "mx"
    rho.layers <- grep("rho", colnames(sp), value = TRUE)
    for (var in rho.layers) {
      dd <- list(
        range = c(-1, 1),
        label = var,
        values = sp[, var]
      )
      dimensions <- c(dimensions, list(dd))
    }

    if (colorby == "pheno") {
      line <- list(
        color = ~rho.PHENO,
        colorscale = "RdBu",
        showscale = TRUE, cmin = -1, cmax = 1
      )
    }
    if (colorby == "score") {
      line <- list(
        color = ~score,
        colorscale = "RdBu",
        showscale = TRUE,
        cmin = min(sp$score, na.rm = TRUE),
        cmax = min(sp$score, na.rm = TRUE)
      )
    }

    fig <- sp %>% plotly::plot_ly(
      type = "parcoords",
      line = line,
      dimensions = dimensions
    )
  }

  if (plotlib == "ggplot") {
    ## data <- iris
    rho.cols <- grep("rho", colnames(sp))
    rho.cols
    S <- sp[, rho.cols]
    sp$alpha <- 0.2
    if (!is.null(hilight)) {
      sp$alpha <- c(0.04, 1)[1 + 1 * (rownames(sp) %in% hilight)]
    }

    fig <- GGally::ggparcoord(
      sp,
      scale = "globalminmax",
      columns = rho.cols,
      groupColumn = "score",
      # order = "anyClass",
      showPoints = TRUE,
      title = "LASAGNA Parallel Coordinates Plot"
    ) + geom_line(
      linewidth = 1,
      aes(linewidth = alpha, alpha = alpha)
    ) +
      scale_colour_gradient2() +
      scale_alpha(guide = "none") +
      labs(x = "", y = "correlation   (rho)")
    fig

    if (!is.null(hilight)) {
      i <- 10
      for (hi in hilight) {
        i <- match(hi, rownames(S))
        tx <- colnames(S)
        ty <- as.matrix(S)[i, ]
        pp <- strsplit(sp$path[i], split = "->")[[1]]
        cc <- c("white", "yellow")[1 + 1 * (pp == hi)]
        tt <- sub(":.*", "", pp)
        pp <- pp[which(!duplicated(tt))]
        pp <- gsub("^[a-zA_Z]+:|\\(.*", "", pp)
        pp <- stringr::str_wrap(pp, 15)
        fig <- fig + ggplot2::annotate("label",
          x = tx, y = ty, label = pp,
          fill = cc, size = (14 * labcex / ggplot2::.pt)
        )
      }
    }
  }

  fig
}


## ======================================================================
## ======================================================================
## ======================================================================

#'
#' @export
nonnegative_transform <- function(x,
                                  separate.sign = FALSE,
                                  shift = TRUE, scale = TRUE,
                                  exponentiate = FALSE) {
  ## break positive and negative separately
  if (separate.sign) {
    x <- (x - rowMeans(x)) ## / matrixStats::rowSds(xx)
    x <- rbind(pmax(x, 0), pmax(-x, 0))
  }

  ## (see IntNMF example) Make all data positive by shifting to positive direction.
  ## Also rescale the datasets so that they are comparable.
  if (shift) {
    # if (!all(dat1>=0)) dat1 <- pmax(dat1 + abs(min(dat1)), .Machine$double.eps)
    # dat1 <- dat1 / max(dat1)
    qmin <- quantile(x, probs = c(0.01))[1]
    if (qmin < 0) x <- pmax(x + abs(qmin), .Machine$double.eps)
  }

  ## (see IntNMF example) Make all data positive by shifting to positive direction.
  ## Also rescale the datasets so that they are comparable.
  if (scale) {
    qmax <- quantile(x, probs = c(0.99))[1]
    x <- x / qmax
  }

  ## do exponential transform. As NMF does multiplicative optimization,
  ## exponential values should work better
  if (exponentiate) {
    x <- pmax(exp(x - min(x)) - 1, 0)
  }

  x
}


#' Integrated NMF.
#'
#' @export
mofa.intNMF <- function(datasets, k = NULL, method = "RcppML",
                        opt.method = "RcppML", maxiter = 200,
                        separate.sign = FALSE,
                        shift = TRUE, scale = TRUE,
                        exponentiate = FALSE) {
  ## method=opt.method=opt.method="RcppML"
  # shift.pos=TRUE;scale=TRUE;separate.sign=TRUE;take.exponent=FALSE

  method <- tolower(method)
  opt.method <- tolower(opt.method)

  ## transform to non-negative
  datasets <- lapply(datasets, function(x) {
    nonnegative_transform(
      x,
      separate.sign = separate.sign,
      shift = shift, scale = scale,
      exponentiate = exponentiate
    )
  })


  lapply(datasets, min)
  lapply(datasets, dim)

  # Find optimum number of clusters for the data
  cophcor <- NULL
  if (is.null(k)) {
    message("[mofa.intNMF] optimizing k... ")
    if (opt.method == "intnmf") {
      tdata <- lapply(datasets, t)
      opt.k <- IntNMF::nmf.opt.k(
        dat = tdata,
        n.runs = 5,
        n.fold = 5,
        k.range = 2:7,
        result = TRUE,
        make.plot = FALSE,
        progress = FALSE
      )
      cophcor <- rowMeans(opt.k)
      k <- names(which.max(cophcor))
      k <- as.integer(sub("^k", "", k))
    } else if (opt.method == "rcppml") {
      cophenetic_coeff <- function(W, xdist) {
        cl <- paste0("C", max.col(W))
        K <- 1 - cor(t(model.matrix(~ 0 + cl)))
        cor(as.vector(K), as.vector(as.matrix(xdist)))
      }
      X <- do.call(rbind, datasets)
      xdist <- dist(t(X))
      cophcor <- c()
      k <- 2
      for (k in 2:20) {
        model <- RcppML::nmf(t(X), k = k, verbose = 0)
        id <- paste0("k=", k)
        cophcor[id] <- cophenetic_coeff(model$w, xdist)
      }
      k <- names(which.max(cophcor))
      k <- as.integer(sub("k=", "", k))
    } else {
      stop("[mofa.intNMF] invalid method", method)
    }
    message("[mofa.intNMF] optimal k = ", k)
  }


  # Find clustering assignment for the samples
  res <- list()
  if (method == "intnmf") {
    fit <- nmf.mnnals(
      dat = lapply(datasets, t),
      k = k,
      maxiter = maxiter,
      st.count = 20,
      n.ini = 15,
      ini.nndsvd = TRUE,
      seed = TRUE
    )
    names(fit)
    res$H <- t(fit$W)
    res$W <- lapply(fit$H, t)
    names(res$W) <- names(datasets)
    res$clusters <- fit$clusters
  } else if (method == "rcppml") {
    np <- sapply(datasets, nrow)
    X <- do.call(rbind, datasets)
    fit <- RcppML::nmf(X, k = k, maxit = maxiter, verbose = 0)
    wt <- rowSums(fit$h**2)**0.5
    res$H <- fit$h / wt
    res$W <- fit$w %*% diag(fit$d * wt)
    res$W <- tapply(1:nrow(res$W), rep(1:length(np), np), function(i) res$W[i, ])
  } else {
    stop("[mofa.intNMF] invalid method", method)
  }

  if (separate.sign) {
    for (i in 1:length(datasets)) {
      n <- nrow(res$W[[i]]) / 2
      wpos <- res$W[[i]][1:n, , drop = FALSE]
      wneg <- res$W[[i]][(n + 1):(2 * n), , drop = FALSE]
      res$W[[i]] <- wpos - wneg
    }
  }

  for (i in 1:length(datasets)) {
    n <- nrow(res$W[[i]])
    rownames(res$W[[i]]) <- colnames(datasets[[i]])[1:n]
    colnames(res$W[[i]]) <- paste0("NMF", 1:k)
  }
  rownames(res$H) <- paste0("NMF", 1:k)
  colnames(res$H) <- colnames(datasets[[1]])
  names(res$W) <- names(datasets)

  res$clusters <- max.col(t(res$H))
  res$cophcor <- cophcor

  return(res)
}


## ======================================================================
## ======================================================================
## ======================================================================
