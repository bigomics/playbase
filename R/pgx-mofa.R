##reticulate::use_miniconda("r-reticulate")

#' 
#'
#' @export
pgx.compute_mofa <- function(pgx, kernel="MOFA", numfactors=10,
                             ntop=1000, add_gsets=TRUE ) {

  has.prefix <- (mean(grepl(":",rownames(pgx$X))) > 0.8)
  is.multiomics <- ( pgx$datatype == "multi-omics" && has.prefix)
  xdata <- NULL
  if(is.multiomics) {
    message("[compute.multi_omics] splitting by omics-type ")
    xdata <- mofa.split_data(pgx$X)
  } else {
    message("[compute.multi_omics] splitting by gene role")
    xdata <- mofa.splitByGeneRole(pgx$X)
  }

  numfactors <- min( numfactors, min(dim(xdata[[1]])) )
  numfactors
  
  ## add geneset layers if asked
  if(add_gsets) {
    gsetX <- mofa.add_genesets(xdata, GMT=pgx$GMT)
    xdata <- c( xdata, gsetX)
  }
  lapply(xdata, dim)
  
  ## samples = pgx$samples
  ## contrasts = pgx$contrasts    
  ## pheno = NULL
  ## kernel = kernel,
  ## scale_views = TRUE 
  ## ntop = 1000
  ## max_iter = 200
  ## num_factors = numfactors

  for(i in 1:length(xdata)) {
    d <- rownames(xdata[[i]])
    ##rownames(xdata[[i]]) <- stringi::stri_trans_general(d, "latin-ascii")
    rownames(xdata[[i]]) <- iconv(d,to="ascii//TRANSLIT")
  }
  
  ## MOFA computation
  message("computing MOFA ...")              
  res <- mofa.compute(
    xdata,
    samples = pgx$samples,
    contrasts = pgx$contrasts,    
    pheno = NULL,
    GMT = pgx$GMT,
    kernel = kernel,
    scale_views = TRUE, 
    ntop = ntop,
    max_iter = 200,
    num_factors = numfactors) 

  names(res)
  
  ## LASAGNA computation
  message("computing LASAGNA model...")            
  las.data <- list(
    X = xdata,
    samples = pgx$samples,
    contrasts = pgx$contrasts
  )
  lasagna <- lasagna.create_model(
    las.data, pheno="contrasts", ntop=1000, nc=20,
    use.graphite = FALSE)
  
  ## pre-compute cluster positions
  xx <- mofa.split_data(lasagna$X)
  xx <- xx[names(xdata)]
   
  message("computing cluster positions (samples)...")
  res$posx <- mofa.compute_clusters(xx, along="samples")
  
  message("computing cluster positions (features)...")
  res$posf <- mofa.compute_clusters(xx, along="features")

  res$lasagna <- lasagna
  
  return(res)
}


#' @export
mofa.compute <- function(xdata,
                         samples,
                         contrasts,
                         kernel = "mofa",
                         pheno = NULL,
                         GMT = NULL,
                         scale_views = TRUE,
                         ntop = 2000,
                         max_iter = 1000,
                         gpu_mode = FALSE,
                         compute.clusters = TRUE,
                         num_factors = 10) {
  

  if(!is.null(ntop) && ntop>0) {
    dbg("[mofa.compute] reducing data blocks: ntop = ", ntop)
    ## fast prioritization of features using SD or rho (correlation
    ## with contrasts).
    i = 1
    num.contrasts <- sign(makeContrastsFromLabelMatrix(contrasts))
    for(i in 1:length(xdata)) {
      d <- xdata[[i]]      
      rho <- cor(t(d), num.contrasts, use="pairwise")
      rho[is.na(rho)] <- 0
      rho2 <- rowMeans(rho**2,na.rm=TRUE) * rowMeans(!is.na(d))
      d <- head(d[order(-rho2),,drop=FALSE],ntop)
      xdata[[i]] <- d
    }
  }
  lapply(xdata, dim)
  
  ## impute missing?
  nna <- sapply(xdata, function(x) sum(is.na(x)))
  nna
  if(any(nna>0)) {
    message("[mofa.compute] warning: imputing missing values in X.")
    xdata[which(nna>0)] <- lapply( xdata[which(nna>0)], svdImpute2)
    ## X <- svdImpute2(X)
  }
  
  ## scale datatypes??
  if(scale_views) {
    message("[mofa.compute] scaling blocks")
    xdata <- lapply(xdata, function(d) d - rowMeans(d,na.rm=TRUE))
    dt.sd <- sapply(xdata, function(d) mean(matrixStats::rowSds(d,na.rm=TRUE),
                                            na.rm=TRUE))
    dt.sd
    xdata <- lapply(names(xdata), function(d) xdata[[d]] / (1e-4 + dt.sd[d]) )
    names(xdata) <- names(dt.sd)
  }

  ## translate to ASCII (IK: should be done in creating genesets...)
  for(i in 1:length(xdata)) {
    d <- rownames(xdata[[i]])
    ##rownames(xdata[[i]]) <- stringi::stri_trans_general(d, "latin-ascii")
    rownames(xdata[[i]]) <- iconv(d,to="ascii//TRANSLIT")
  }

  ## add prefix???
  xdata <- mofa.prefix(xdata)
  X <- do.call( rbind, xdata )
  X <- X - rowMeans(X,na.rm=TRUE)
  dim(X)

  ## create pheno from contrasts if missing
  if(is.null(pheno)) {
    message("[mofa.compute] creating pheno from contrasts")
    ct <- contrasts
    ct[is.na(ct)] <- "NA"
    pheno <- apply(ct,1,paste,collapse="_")
    table(pheno)
  }
  
  model <- NULL
  kernel <- tolower(kernel)

  if( kernel == "mofa") {
    ##num_factors = 10; max_iter = 1000
    xdata <- mofa.prefix(xdata)  ## MOFA needs prefix
    obj <- MOFA2::create_mofa(xdata, groups=NULL)
    MOFA2::samples_metadata(obj) <- data.frame( sample=rownames(samples), samples )
    
    data_opts <- MOFA2::get_default_data_options(obj)
    data_opts$scale_views <- FALSE
    data_opts
    
    model_opts <- MOFA2::get_default_model_options(obj)
    model_opts$num_factors <- num_factors
    ##model_opts$likelihoods <- c("gaussian","gaussian","gaussian","bernoulli")
    
    train_opts <- MOFA2::get_default_training_options(obj)
    train_opts$gpu_mode <- gpu_mode
    train_opts$maxiter <- max_iter
    ##train_opts$gpu_mode <- FALSE  
  
    obj <- MOFA2::prepare_mofa(
      object = obj,
      data_options = data_opts,
      model_options = model_opts,
      training_options = train_opts
    )
    
    ## where to save??? Do we need to save???
    outfile = file.path(tempdir(),"mofa-model.hdf5")
    ##outfile = file.path(getwd(),"mofa-model.hdf5")
    ##  suppressMessages(suppressWarnings(
    model <- MOFA2::run_mofa(obj, outfile=outfile, save_data=TRUE,
                             use_basilisk=FALSE)

    ##model <- load_model(outfile, remove_inactive_factors = FALSE)
    model <- MOFA2::impute(model)
    factors <- MOFA2::get_factors(model, factors = "all")
    weights <- MOFA2::get_weights(model, views = "all", factors = "all")
    W <- do.call(rbind, weights)
    F <- do.call(rbind, factors)

  } else if( tolower(kernel) %in% c("pca","svd")) {
    model <- irlba::irlba( X, nv=num_factors, maxit=max_iter, work=100)
    W <- model$u 
    F <- model$v %*% diag(model$d)
    rownames(W) <- rownames(X)
    rownames(F) <- colnames(X)
    colnames(W) <- paste0("PC",1:ncol(W))
    colnames(F) <- paste0("PC",1:ncol(F))
  } else if( kernel == "diablo") {
    dtx <- sub(":.*","",rownames(X))
    xx <- tapply( 1:nrow(X), dtx, function(i) X[i,,drop=FALSE])
    if(is.null(pheno)) {
      stop("pheno must be provided for diablo")
    }
    if(any(is.na(pheno))) {
      message("[DIABLO] warning: missing values in phenotype. removing samples.")
      sel <- which(!is.na(pheno))
      pheno <- pheno[sel]
      xx <- lapply(xx, function(x) x[,sel,drop=FALSE])
      X <- X[,sel,drop=FALSE]
      samples <- samples[sel,]
    }
    mix <- mixomics.compute_diablo(xx, y=pheno, ncomp=num_factors) 
    W <- mix$loadings
    F <- mix$scores
    model <- mix
  } else if( kernel == "mcia") {
    data_blocks <- lapply(xdata, t)
    data_blocks_mae <- nipalsMCIA::simple_mae(
      data_blocks,
      row_format = "sample",
      colData = samples )
    mcia <- nipalsMCIA::nipals_multiblock(
      data_blocks_mae,
      num_PCs = num_factors,
      plot="none")
    F <- mcia@global_scores
    W <- mcia@global_loadings
    dim(F)
    colnames(W) <- paste0("factor",1:ncol(W))
    colnames(F) <- paste0("factor",1:ncol(F))
    model <- mcia
  } else if( kernel == "wmfcna") {
    model <- wmfcna.compute(
      X, y = pheno, k = num_factors, power = 1, vex = 1,
      ewidth = 1, min.cor = 0.7)
    W <- model$W
    F <- t(model$meExpr)
  } else {
    message("[mofa.compute] invalid kernel")
    return(NULL)
  }
  
  W <- W[rownames(X),]
  F <- F[colnames(X),]  
  dt <- sub(":.*","",rownames(X))
  ww <- tapply(1:nrow(W),dt,function(i) W[i,,drop=FALSE] )
  xx <- tapply(1:nrow(X),dt,function(i) X[i,,drop=FALSE] )
  ww <- ww[names(xdata)]
  xx <- xx[names(xdata)]
  ww <- mofa.strip_prefix(ww)
  xx <- mofa.strip_prefix(xx)  
  
  ## variance explained
  w.ssq <- sapply(ww, function(w) colSums(w**2))
  w.ssq <- w.ssq * sqrt(colSums(F**2)) ## check!!
  V <- t(w.ssq / sum(w.ssq))
  rownames(V) <- names(xx)
  colnames(V) <- colnames(W)

  ## Covariate x Factor correlation
  M <- samples
  M$sample <- NULL
  M$group <- NULL
  M <- expandPhenoMatrix(M, drop.ref=FALSE)
  M <- as.matrix(M)
  Z <- cor(M, F, use="pairwise")
  Z[is.na(Z)] <- 0
  colnames(Z) <- colnames(W)

  message("computing factor enrichment...")
  gsea <- mofa.compute_enrichment(
    ww, G=GMT, filter=NULL, ntop=NULL) 

  ## compute enrichment for all contrasts
  message("computing phenotype enrichment...")
  fc.gsea <- NULL
  ff <- NULL
  if(!is.null(contrasts)) {
    contrasts <- contrasts[rownames(samples),,drop=FALSE]

    mfc <- list()
    ct=colnames(contrasts)[1]
    for(ct in colnames(contrasts)) {
      y <- contrasts[,ct]
      dx <- lapply(xx, function(x) gx.limma(x, y, lfc=0, fdr=1))
      mfc[[ct]] <- lapply( dx, function(d) matrix(d$logFC, ncol=1,
        dimnames=list(rownames(d),ct)))
    }
    
    ff <- mfc[[1]]
    if(length(mfc)>1) {
      for(i in 2:length(mfc)) {
        for(j in 1:length(ff)) ff[[j]] <- cbind(ff[[j]], mfc[[i]][[j]])
      }
    }
    fc.gsea <- mofa.compute_enrichment(ff, G=GMT)
  }
  
  ## create graphs
  if( kernel == "wmfcna") {
    graphs <- list(
      factors = model$graph,
      features = model$subgraphs
    )
  } else {
    graphs <- mofa.factor_graphs(F, W, X, y=pheno, n=20)
  }

  num.contrasts <- sign(makeContrastsFromLabelMatrix(contrasts))
  
  res <- list(
    model = model,
    samples = samples,
    contrasts = contrasts,
    num.contrasts = num.contrasts,
    pheno = pheno,
    F = F,
    W = W,
    V = V,
    X = X,
    Z = Z,
    xx = xx,
    ww = ww,
    fc = ff,
    M = M,
    gsea = gsea,
    fc.gsea = fc.gsea,
    graphs = graphs
  )
  
  res
}

#' @export
mofa.scale_views <- function(xdata) {
  xdata <- lapply(xdata, function(d) d - rowMeans(d,na.rm=TRUE))
  dt.sd <- sapply(xdata, function(d) mean(matrixStats::rowSds(d,na.rm=TRUE),
                                          na.rm=TRUE))
  dt.sd
  xdata <- lapply(names(xdata), function(d) xdata[[d]] / (1e-9+dt.sd[d]) )
  names(xdata) <- names(dt.sd)
  xdata
}

#' @export
mofa.add_genesets <- function(xdata, GMT=NULL, datatypes=NULL) {
  names(xdata)
  if(is.null(datatypes)) {
    datatypes <- intersect(c("gx","px","mx"), names(xdata))
  }
  if(length(datatypes)==0) {
    stop("must provide datatypes")
  }

  if(is.null(GMT)) {
    ## if("mx" %in% datatypes) {
    ##   load("~/Playground/public-db/pathbank.org/PATHBANK-matrix.rda",verbose=TRUE)
    ##   message("[mofa.add_genesets] adding PATHBANK")
    ##   GMT <- PATHBANK
    ## } else {
    ##   message("[mofa.add_genesets] adding GSETxGENE")
    ##   stdGMT <- playdata::GSETxGENE
    ##   colnames(stdGMT) <- paste0("SYMBOL:",colnames(stdGMT))
    ## }
    message("combining GSETxGENE and PATHBANK genesets...")
    stdGMT <- playdata::GSETxGENE
    colnames(stdGMT) <- paste0("SYMBOL:",colnames(stdGMT))
    pathbankGMT <- playdata::MSETxMETABOLITE
    GMT <- merge_sparse_matrix( pathbankGMT, stdGMT[,])
    dim(GMT)
  } else {
    GMT <- Matrix::t(GMT)   ## genesets on rows
  }
  
  if(is.null(GMT)) {
    message("ERROR: could not determing GMT. Please provide.")
    return(NULL)
  }

  ## convert non-ascii characters....
  ##rownames(GMT) <- iconv(rownames(GMT), "latin1", "ASCII", sub="")
  rownames(GMT) <- stringi::stri_trans_general(rownames(GMT), "latin-ascii")
  colnames(GMT) <- sub(".*:","",colnames(GMT))
  
  gsetX <- list()
  dt=datatypes[1]
  for(dt in datatypes) {
    rX <- xdata[[dt]] - rowMeans(xdata[[dt]], na.rm=TRUE)
    rownames(rX) <- sub(".*:","",rownames(rX))
    sel <- intersect(rownames(rX), colnames(GMT))
    G <- GMT[,sel,drop=FALSE]
    rX <- rX[sel,,drop=FALSE]
    rho <- gset.rankcor(rX, G)$rho
    ii <- as.vector(which(rowMeans(is.na(rho)) < 1))    
    rho <- rho[ii,,drop=FALSE]
    gsetX[[dt]] <- rho
  }
  lapply(gsetX,dim)
  
  ## make sure they align. Only intersection????
#  kk <- Reduce(intersect, lapply(gsetX, rownames))
#  gsetX <- lapply( gsetX, function(m) m[kk,])
  
  names(gsetX) <- paste0("gset.",names(gsetX))  
  gsetX
}

#' Compute clusters for MOFA factors.
#' 
#' @export
mofa.compute_clusters <- function(xx, matF=NULL, along="samples") {

  if(along == "samples") {
    ## cluster samples
    sxx <- lapply(xx, function(d) {d=t(scale(t(d),scale=FALSE));d[is.na(d)]=0;d})
    if(!is.null(matF)) sxx$MOFA <- t(matF)
    sxx <- lapply(sxx, function(d) d[matrixStats::rowSds(d,na.rm=TRUE)>0.01,])
    nb <- round(min(15,ncol(xx[[1]])/4))
    posx <- lapply(sxx, function(x) uwot::umap(t(x),n_neighbors=nb))
    for(i in 1:length(posx)) rownames(posx[[i]]) <- colnames(sxx[[i]])
  }

  if(along == "features") {
    ## cluster features
    sxx <- lapply(xx, function(d) {d=t(scale(t(d),scale=FALSE));d[is.na(d)]=0;d})
    posx <- list()
    for(i in 1:length(sxx)) {
      nb <- min(15,nrow(sxx[[i]])/4)
      posx[[i]] <- uwot::umap(sxx[[i]],n_neighbors=nb)
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
mofa.compute_enrichment <- function(ww, ww.types=NULL, filter=NULL,
                                    G=NULL,  ntop=NULL) {

  ##ww.types=filter=G=ntop=NULL

  ## no genesets
  ww <- ww[grep("^gset",names(ww),invert=TRUE)]
  names(ww)
  
  ## determine type
  if(is.null(ww.types)) {
    ww.types <- ifelse( grepl("^mx|^metabo",names(ww),ignore.case=TRUE),
                       "ChEBI", "SYMBOL" )
  }
  ww.types

  ## USE internal GMT????
  has.mx <- any(ww.types == "ChEBI")
  has.px <- any(ww.types == "SYMBOL") 
  if(is.null(G)) {
    G <- NULL
    ## add metabolomic gene sets   
    if (has.mx) {
      info("[pgx.add_GMT] Adding metabolomics genesets")
      G <- Matrix::t(playdata::MSETxMETABOLITE)
      rownames(G) <- sub(".*:","",rownames(G)) ## TEMPORARY!!!
    }
    ## add SYMBOL (classic) gene sets 
    if (has.px) {
      info("[pgx.add_GMT] Adding transcriptomics/proteomics genesets")
      if(is.null(G)) {
        G <- Matrix::t(playdata::GSETxGENE)
      } else {
        G <- merge_sparse_matrix(G, Matrix::t(playdata::GSETxGENE))      }
    }
  }
  rownames(G) <- sub(".*:","",rownames(G))  
  dim(G)

  ww <- mofa.strip_prefix(ww)
  wnames <- unlist(lapply(ww, rownames))
  pp <- intersect(wnames,rownames(G))  
  G <- G[pp,]
  G <- G[, Matrix::colSums(G!=0)>=3]
  dim(G)
  gg <- rownames(G)
  
  gsea <- list()
  i=1
  for(i in 1:length(ww)) {
    w1 <- ww[[i]]
    nc <- sum(rownames(w1) %in% gg)
    if(nc>0) {
      f1 <- gset.rankcor( w1, G, compute.p=TRUE)
      dt <- names(ww)[i]
      gsea[[dt]] <- f1
    }
  }
  names(gsea)

  rho <- lapply( gsea, function(x) x$rho )
  p.value <- lapply( gsea, function(x) x$p.value )
  q.value <- lapply( gsea, function(x) x$q.value )      

  scores <- rho
  nt <- length(scores)
  for(i in 1:nt) {
    s <- (scores[[i]] * (1 - p.value[[i]]))  ## p or q??
    s[is.na(s)] <- 0  ## really??
    scores[[i]] <- s**2
  }
  ##multi.score <- sqrt( Reduce("+", x=scores) / nt )
  multi.score <- abs(Reduce("*", x=scores))^(1/(2*nt))  
  table(is.na(multi.score))
  
  enr <- list()
  for(i in 1:ncol(multi.score)) {
    df <- data.frame(
      multi.score = multi.score[,i],
      score = I(sapply(scores, function(x) x[,i])),
      rho  = I(sapply(rho, function(x) x[,i])),
      pval = I(sapply(p.value, function(x) x[,i])),
      qval = I(sapply(q.value, function(x) x[,i]))
    )
    df <- df[!is.na(df$multi.score),,drop=FALSE]
    ct <- colnames(multi.score)[i]
    enr[[ct]] <- df
  }
  lapply(enr, dim)
  
  return(enr)
}


#' @export
mofa.split_data <- function(X) {
  dtype <- sub(":.*","",rownames(X))
  xx <- tapply(1:nrow(X), dtype, function(i) X[i,,drop=FALSE] )
  xx <- mofa.strip_prefix(xx)
  xx
}

#' @export
mofa.merge_data <- function(xx) {
  do.call(rbind, mofa.prefix(xx))
}

#' @export
mofa.prefix <- function(xx) {
  xx <- mofa.strip_prefix(xx)
  i=1
  for(i in 1:length(xx)) {
    dt <- paste0(names(xx)[i],":")
    if(is.null(dim(xx[[i]]))) {
      names(xx[[i]]) <- paste0(dt,names(xx[[i]]))
    } else {
      rownames(xx[[i]]) <- paste0(dt,rownames(xx[[i]]))
    }
  }
  xx
}

mofa.strip_prefix <- function(xx) {
  i=1
  for(i in 1:length(xx)) {
    dt <- paste0("^",names(xx)[i],":")
    if(is.null(dim(xx[[i]]))) {
      names(xx[[i]]) <- sub(dt,"",names(xx[[i]]))
    } else {
      rownames(xx[[i]]) <- sub(dt,"",rownames(xx[[i]]))
    }
  }
  xx
}

#' @export
mofa.correlate_covariate <- function(res, Y) {  
  K  <- cov(res$F, Y)
  R <- cov( t(res$W), K)
  rr <- lapply(res$ww, function(n) R[rownames(n),,drop=FALSE])
  rr <- mofa.strip_prefix(rr)
  list(
    R = R,
    rr = rr
  )
}

#' @export
mofa.plot_weights <- function(weights, k, ntop=10, cex.names=0.9,
                              maxchar=999) {
  wt <- lapply( weights, function(w) w[,k])
  ## plot the weights for each datatype (Factor1)
  ##par(mfrow=c(2,3), mar=c(4,30,4,2))
  v = names(wt)[1]
  for(v in names(wt)) {
    w1 <- wt[[v]]
    names(w1) <- gsub(".*:| \\(SMP.*","",names(w1))
    w1 <- w1[order(-abs(w1))]
    w1 <- w1[!duplicated(names(w1))]
    names(w1) <- stringr::str_trunc(names(w1),maxchar)
    xlab <- paste(k," weight")
    if(!is.null(ntop) && ntop>0 && ntop < length(w1)) {
      w1 <- w1[w1!=0]
      i0 <- head(order(w1),ntop)
      i1 <- tail(order(w1),ntop)
      if(length(i0)==0) i1 <- tail(order(w1),2*ntop)
      if(length(i1)==0) i0 <- head(order(w1),2*ntop)    
      ii <- unique(c(i0,i1))
      barplot( w1[ii], horiz=TRUE, las=1,
              xlab=xlab, cex.names=cex.names)
    } else {
      names(w1) <- NULL
      barplot( sort(w1), border = NA, col="grey50",
        horiz=TRUE, las=1, xlab=xlab, space=0)
      abline(v=0)
      mtext(paste0("sorted features (N=",length(w1),")"), side=2)
    }
    title(toupper(v), cex.main=1.2)
  }
}

#' @export
mofa.plot_loading_heatmap <- function(mofa, k=NULL, ntop=50,
                                      type = c("heatmap","splitmap")[1],
                                      annot = c("scores","pheno")[1],
                                      split = TRUE, maxchar=999,
                                      show_types = NULL,
                                      ... ) {

  if(!is.null(k)) {
    xx <- mofa.prefix(mofa$ww)
  } else {
    xx <- mofa.prefix(mofa$xx)
  }
  
  if(!is.null(show_types)) {
    xx <- xx[intersect(show_types,names(xx))]
  }

  ntop1 <- ntop / length(xx)
  if(!is.null(k)) {
    top <- lapply(xx, function(w) {
      w1 <- w[,k];
      names(w1) <- rownames(w);
      w1 <- w1[which(w1!=0 & !is.na(w1))];
      head(names(sort(-abs(w1))),ntop1)
    })
  } else {
    top <- lapply(xx, function(x) {
      sdx <- matrixStats::rowSds(x,na.rm=TRUE)
      head(names(sort(-sdx)),ntop1)
    })
  }
  top <- unlist(top)
  topX <- mofa$X[top,,drop=FALSE]
  rownames(topX) <- stringr::str_trunc(rownames(topX),maxchar)

  if(annot=="scores") {
    aa <- data.frame(mofa$F)
  }
  if(annot=="pheno") {
    aa <- mofa$samples
  }
  rownames(aa) <- colnames(topX)

  
  if(type == "heatmap") {
    par(mar=c(0,0,0,0))
    gx.heatmap( topX, mar=c(8,10), key=FALSE, keysize=1,
               ##cexRow = 1, col.annot = aa )
               cexRow = 1, col.annot = aa, ... )    
  } else if(type == "splitmap") {
    dx <- sub(":.*","",rownames(topX))
    if(!split) dx <- NULL
    gx.splitmap( topX, 
                split=dx, na_col = "white", softmax=TRUE,
                rowlab.maxlen = 80, rownames_width = 140,
                show_legend=FALSE, show_key = FALSE,
                col.annot = aa, ... )
  } else if(type == "graph") {
    
  } else {
    message("unknown graph type")
  }

  
}

#' @export
mofa.plot_factor_trait <- function(mofa, Y=NULL,
                                   main="Factor-Trait heatmap",
                                   type = c("wgcna","splitmap"),
                                   par = TRUE, cex_text=NULL,
                                   cluster = TRUE,
                                   ...) {

  ##type = c("wgcna","splitmap")
  type <- type[1]
  F <- mofa$F
  if(is.null(Y)) {
    Y <- expandPhenoMatrix(mofa$samples, drop.ref=FALSE)
  }
  R <- cor(Y, F, use="pairwise")

  if(cluster) {
    ii <- hclust(dist(R))$order
    jj <- hclust(dist(t(R)))$order
    R <- R[ii,jj]
  }
  
  if(type=="splitmap") {
    gx.splitmap(
      t(R), nmax=50, scale='none', main=main,
      col.annot = mofa$samples, split=1,
      show_legend=FALSE, show_key=FALSE,
      ...)
  }
  
  if(type=="wgcna") {
    ##par(mfrow=c(1,1))
    if(par) par(mar=c(6,5,2,1))    
    ftext <- round(t(R), digits=2)
    if(is.null(cex_text)) {
      cex.text = min(0.8, max(0.3, 6/ncol(R)))
    } else {
      cex.text = cex_text
    }
    
    WGCNA::labeledHeatmap(
      Matrix = t(R),
      xLabels = rownames(R), 
      yLabels = colnames(R), 
      textMatrix = ftext,
      cex.text = cex.text,
      # cex.lab = 1, 
      #  ySymbols = colnames(res$F),
      colorLabels = TRUE, 
      colors = WGCNA::blueWhiteRed(50), 
      setStdMargins = FALSE, 
      zlim = c(-1,1),
      main = main,
      #...
    )
  }  
}

#' @export
mofa.plot_factor_correlation <- function(mofa, 
                                         main="Factor correlation heatmap",
                                         marx = 1,
                                         ...) {
  R <- cor(mofa$F, use="pairwise")
  lablength <- max(nchar(rownames(R)))
  gx.heatmap(
    R, nmax=50, scale='none', main=main, sym=TRUE,
    key=FALSE, keysize=0.8, mar=c(1,1)*(lablength)*marx,
    ...)
}

#' @export
mofa.plot_factor_boxplots <- function(mofa, k=1, pheno=NULL,
                                      by="condition", par=TRUE) {
  
  if(by == "condition") {
    ## boxplots
    if(par) par(mfrow=c(2,2), mar=c(5,5,3,2))
    for(j in 1:ncol(samples)) {
      cond <- colnames(samples)[j]
      y <- factor(samples[,j])
      f1 <- mofa$F[,k]
      fn <- colnames(mofa$F)[k]
      tt <- paste(fn,"by", colnames(samples)[j])
      boxplot( f1 ~ y, main = tt, xlab = cond,
              ylab = paste(fn,"score"))
    }
  }

  if(by == "factor") {
    k=1
    if(!is.null(pheno) && is.integer(pheno)) pheno <- colnames(res$samples)[pheno]
    if(par) par(mfrow=c(3,4))
    for(k in 1:min(12,ncol(mofa$F))) {
      y <- factor(samples[,pheno])
      f1 <- mofa$F[,k]
      tt <- paste0("Factor",k)
      boxplot( f1 ~ y, main=tt,
              xlab = pheno, ylab="Score")
    }
  }
  
}


#' Covariate x Factor correlation plot. Shows how factors are
#' correlated with response covariates.
#'
#' @export
mofa.plot_covariate_correlation <- function(mofa, collapse=FALSE,
                                            mar = c(5,10), nmax = 60, 
                                            ...) {

  ## Covariate x Factor correlation
  Z <- mofa$Z
  if(collapse) {
    rownames(Z) <- sub("=.*","",rownames(Z))
    Z <- rowmean(Z**2)**0.5
  }

  if(nrow(Z)==1)  {
    Z <- rbind(Z," "=Z[1,])  ## hack for single row...
  }

  gx.heatmap(
    Z, scale='none', mar = mar, nmax = nmax, 
    key=FALSE, keysize=0.8, cexCol=1, cexRow=1,
    ... )
}


mofa.combine_layers <- function(xx, weights=1) {
  k=1
  if(length(weights)==1) weights <- rep(weights,length(xx))  
  for(i in 1:length(xx)) {
    rownames(xx[[i]]) <- sub(".*:","",rownames(xx[[i]])) ## strip
  }
  kk <- Reduce(intersect, lapply(xx, colnames))
  gg <- Reduce(intersect, lapply(xx, rownames))    
  xx <- lapply(xx, function(x) x[gg,kk])
  cx.pos <- xx[[1]]*0 + 1
  cx.neg <- xx[[1]]*0 + 1
  for(i in 1:length(xx)) {
    xi <- xx[[i]] - rowMeans(xx[[i]], na.rm=TRUE)
    cx.pos <- cx.pos * pmax(xi,0) * weights[i]
    cx.neg <- cx.neg * pmax(-xi,0) * weights[i]
  }
  cx <- cx.pos - cx.neg
  cx  ## combined
}


#' @export
mofa.plot_enrichment <- function(gsea, type="barplot",
                                 remove.dup = TRUE, ntop=20,
                                 filter = NULL, select = NULL,
                                 strip.names = FALSE, par = NULL,
                                 title = "multi-omics enrichment"){
  ##k=1;ntop=20
  S <- gsea
  if(nrow(S)==0) return(NULL)

  if(!is.null(select) && length(select)>0) {
    S <- S[select,,drop=FALSE]
  } else {
    S <- S[order(-S$multi.score),,drop=FALSE]
  }
  
  if(!is.null(filter)) {
    S <- S[grep(filter,rownames(S)),,drop=FALSE]
  }
  
  if(remove.dup) {
    sname <- gsub(".*:| \\(.*","",rownames(S))
    S <- S[which(!duplicated(sname)),,drop=FALSE]
  }

  S <- abs(S$score)
  topS <- head(S, ntop)
  topS <- topS[nrow(topS):1,,drop=FALSE]
  plot <- NULL
  if(strip.names) {
    rownames(topS) <- gsub(" \\(.*","",rownames(topS))
  }
  
  if(type=='barplot') {
    if(is.null(par)) {
      par(mfrow=c(1,2), mar=c(4,10,4,2))
      plot.new()
    } else {
      par(par)
    }
    barplot( t(topS), horiz=TRUE, las=1, xlab="cumulative abs.enrichment")
    cc <- grey.colors(ncol(topS))
    legend("bottomright", legend=colnames(topS),
           fill=cc, inset=c(0.0,0.02), bg='white',
           y.intersp = 0.8, cex=0.85)
    title( title, outer=TRUE, cex.main=1.4, line=-1.2)
  }

  if(type=='lollipop') {
    ## need double lolliplot lines for multi-omics
    top.nes <- topS[,1]
    names(top.nes) <- rownames(topS)
    plot <- ggLollipopPlot( top.nes, xlab="enrichment (NES)" )
  }
  plot
}



#' 
#' @export
mofa.plot_multigsea <- function(gsea, type1, type2, k=1,
                                main=NULL, hilight=NULL) {
  ##gsea=res$gsea;type1="px";type2="mx"
  ## using pre-computed GSEA

  s1 <- gsea[[k]]$score[,type1]
  s2 <- gsea[[k]]$score[,type2]
  qq <- gsea[[k]]$pval[,c(type1, type2)]
  qq[is.na(qq)] <- 1
  ##qcomb <- apply(qq,1, function(pp) metap::meanz(pp)$p)
  qcomb <- apply(qq,1,max)
  s1 <- s1 + 1e-2*rnorm(length(s1))
  s2 <- s2 + 1e-2*rnorm(length(s2))
  cex1 <- (1-qcomb)**2
  col1 <- ifelse( is.null(hilight), "black", "grey60")
  plot( s1, s2, cex = 1.5*cex1, col=col1,
       xlab = paste(type1, "enrichment  [NES*(1-p)]"),
       ylab = paste(type2, "enrichment  [NES*(1-p)]")
  )
  if(is.null(title)) {
    title <- paste("multiGSEA:", toupper(type1), "vs.",
      toupper(type2) )
  } else {
    title <- main
  }
  title(title, cex.main=1.3)
  
  abline(h=0, v=0, lty=3)
  kk <- rownames(gsea[[k]])
  if(!is.null(hilight) && length(hilight)>0 ) {
    if( length(intersect(hilight,kk))) {
      sel <- intersect(hilight,kk)
    } else {
      sel <- grep(hilight, kk, value=TRUE)
    }
    if(length(sel)) {
      sel <- match(sel, kk)
      points( s1[sel], s2[sel], pch=20, col="red2",
             lwd=2, cex=1.4*cex1[sel] )
    }
  }
}

mofa.plotVar <- function(mofa, comp=1:2, style="correlation",
                         textlabel=TRUE) {

  if(style == "correlation") {
    xdata <- mofa$ww
    dtypes <- names(mofa$ww)
    for(i in 1:length(xdata)) {
      x <- mofa$xx[[i]]
      w <- mofa$ww[[i]][,comp]
      f <- mofa$F[,comp]
      sel <- which(rowSums(abs(w))>0)  ## only non-zero variables
      rho <- cor( t(x[sel,]), f)
      if(i == 1) {
        plot( rho[,1], rho[,2], pch=19, col=1+i,
             xlab = colnames(rho)[1],
             ylab = colnames(rho)[2],           
             asp=1, xlim=c(-1,1), ylim=c(-1,1))
      } else {
        points( rho[,1], rho[,2], pch=19, col=1+i )
      }
      if(textlabel) {
        text( rho[,1], rho[,2], rownames(rho), pos=3)
      }
    }
    title("Correlation Circle Plot")
    abline(v=0,h=0, lty=2)
    plotrix::draw.circle(0,0,radius = c(0.5,1))
    plotrix::draw.circle(0,0,radius=0.5)
    legend("bottomright", legend=dtypes,
           fill=2:100, cex=0.9, y.intersp=0.85)
  }

  if(style == "loading") {
    dtypes <- names(mofa$ww)
    ww <- lapply(mofa$ww, function(w) w[,comp,drop=FALSE])
    maxw <- max(sapply(ww, function(w) max(abs(w),na.rm=TRUE)))
    for(i in 1:length(dtypes)) {
      w <- ww[[i]][,comp]
      sel <- which(rowSums(abs(w))>0)  ## only non-zero variables
      if(i == 1) {
        plot( w[,1], w[,2], pch=19, col=1+i,
             xlab = colnames(w)[1],
             ylab = colnames(w)[2],           
             asp=1, xlim=c(-1,1)*maxw,
             ylim=c(-1,1)*maxw )
      } else {
        points( w[,1], w[,2], pch=19, col=1+i )
      }
      text( w[,1], w[,2], rownames(w), pos=3)      
    }
    title("Variable weights (Factor loading)")
    abline(v=0,h=0, lty=2)
    legend("bottomright", legend=dtypes,
           fill=2:100, cex=0.9, y.intersp=0.85)
  }
  
}

#' @export
mofa.factor_graphs <- function(F, W, X, y, n=10,
                               ewidth=1, vsize=1) {
  
  ## cluster-reduced graph
  create_graph <- function(gx, y, min.cor) {
    rho <- cor(t(gx),use="pairwise")
    gr <- igraph::graph_from_adjacency_matrix(
      rho, weighted=TRUE, diag=FALSE, mode = "undirected")
    ew <- igraph::E(gr)$weight
    ew <- (ew - min(ew,na.rm=TRUE)) /
      (max(ew,na.rm=TRUE)-min(ew,na.rm=TRUE))
    igraph::E(gr)$weight <- ew
    gr <- igraph::delete_edges(gr, edges = which(ew < min.cor))
    gr <- igraph::simplify(gr)
    sdx <- matrixStats::rowSds(gx)
    igraph::V(gr)$size <- 20*vsize*(sdx/max(sdx))
    if(!is.null(y)) {
      y <- factor(as.character(y))
      ii <- which(!is.na(y))
      if(length(unique(y[ii]))==2) {
        val <- cor( t(gx[,ii,drop=FALSE]), as.integer(y[ii]))[,1]
      } else if(length(unique(y[ii])) > 2) {
        res <- gx.limmaF( gx[,ii,drop=FALSE], y[ii], fdr=1, lfc=0)
        val <- res$logFC
      } else {
        val <- NULL
      }
      #igraph::V(gr)$size <- 20*vsize*(abs(val)/max(abs(val)))
      if(!is.null(val)) igraph::V(gr)$color <- colorscale(val, gamma=2)
    }
    igraph::E(gr)$width <- 3 * ewidth * abs(igraph::E(gr)$weight)**2
    gr
  }
  
  meGraph <- create_graph(t(F), y, min.cor=0.33 )
  me.size <- rowSums(abs(t(W)) > 0.5*apply(abs(W),2,max))  
  igraph::V(meGraph)$size <-  25*vsize*(me.size/max(me.size))**0.5

  topfeatures <- function(x,n) {
    names(x)<-rownames(W); x=sort(x[x!=0]);
    unique(c(names(tail(x,n)), names(head(x,n))))
  }
  topff <- apply(W, 2, function(x) topfeatures(x,n=n), simplify=FALSE)  
  subgraphs <- list()
  k=1
  for(k in 1:length(topff)) {
    sel <- topff[[k]]
    if(length(sel)>1) {
      gx <- X[sel,,drop=FALSE]
      subgraphs[[k]] <- create_graph(gx, y, min.cor=0.33 )
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
mofa.plot_centrality <- function(res, k, y, show_types=NULL,
                                 main="centrality vs. foldchange") {
  
  gr <- res$graph$features[[k]]
  gg <- igraph::V(gr)$name
  gr <- igraph::mst(gr, weight=1/igraph::E(gr)$weight)
  ctx <- igraph::page_rank(gr)$vector
  names(ctx) <- gg

  gg <- intersect(gg, rownames(res$X))
  xx <- res$X[gg,]
  
  if(!is.null(show_types) && length(show_types)) {
    dt <- sub(":.*","",rownames(xx))
    dbg("[mofa.plot_centrality_vs] 1b: unique(dt) = ",unique(dt) )
    sel <- which(dt %in% show_types)
    dbg("[mofa.plot_centrality_vs] 1b: lenght(sel) = ",length(sel) )    
    xx <- xx[sel,, drop=FALSE]
  }
  
  nlev <- length(unique(y[!is.na(y)]))
  test.type <- "" 
  if(nlev==2) {
    lm <- gx.limma(xx, y, fdr=1,lfc=0)
    test.type <- "limma" 
  } else { 
    lm <- gx.limmaF(xx, y, fdr=1,lfc=0)
    test.type <- "limmaF" 
  }
  
  rx <- lm$logFC
  names(rx) <- rownames(lm)
  gg <- rownames(xx)
  rx <- rx[gg]
  ry <- ctx[gg]
  dx <- 0.06*c(-1,1) * diff(range(rx))
  plot(rx, ry,
       pch = 20, cex=1,
       ylim = c(0,1.05*max(ry)),
       xlim = range(rx) + dx,
       ylab = "centrality  (page_rank)",
       xlab = paste0("logFC  (",test.type,")"))
  title(main)
  abline(h=0, v=0, lty=2)
  text(rx, ry, gg, pos=3, offset=0.3)
}

#' @export
mofa.plot_module <- function(graph, mst=TRUE, nlabel=10, rm.single=FALSE,
                             highlightby = "centrality.prize", cex=1,
                             physics = TRUE,
                             plotlib=c("igraph","visnet")) {
  if(rm.single) {
    vdeg <- igraph::degree(graph)
    graph <- igraph::subgraph(graph, vids=which(vdeg>0))
  }
  if(mst) {
    wt <- igraph::E(graph)$weight
    graph <- igraph::mst(graph, weights=1/wt)
  }

  igraph::V(graph)$prize <- igraph::V(graph)$size
  igraph::V(graph)$type <- 1
  igraph::V(graph)$foldchange <- 1
  ew <- igraph::E(graph)$width
  igraph::E(graph)$width <- 10 * playbase::uscale(ew)**2

  if(plotlib=="igraph") {
    plot(graph)
  }
  if(plotlib=="visnet") {
    
    class(graph) <- c("PCSF", "igraph")
    igraph::E(graph)$weight <- 1 / igraph::E(graph)$weight
    plotPCSF(
      graph,
      highlightby = highlightby,
      nlabel = nlabel,
      edge_width = 5*cex,
      node_cex = max(10,30*cex),
      label_cex = 25*cex,
      layout = "layout_with_kk",
      physics = physics,
      plotlib = "visnet"
    )

  }
  
}

#' @export
mofa.predict <- function(mofa, newdata) {

  if(!is.null(mofa$model) && "block.splsda" %in% class(mofa$model)) {
    message("This seems to be a MixOmics Diablo model. Please use mixomics.predict() for better performance.")
  }

  if(is.list(newdata)) {
    newdata <- mofa.scale_views(newdata)    
    newdata <- mofa.prefix(newdata)
    newdata <- do.call(rbind, newdata)
  }

  W <- mofa$W
  F <- mofa$F  

  kk <- intersect(rownames(newdata),rownames(W))
  rho1 <- cor( newdata[kk,], W[kk,], use="pairwise" )
  rho2 <- cor( t(rho1), t(F), use="pairwise"  )
  max.nb <- max.col(rho2)
  predicted <- mofa$samples[max.nb,]
  rownames(predicted) <- colnames(newdata)
  colnames(predicted) <- paste0("predicted.",colnames(predicted))
  predicted
}


#' 
#'
#' @export
mofa.splitByGeneRole <- function(X) {
  G <- playdata::GSETxGENE
  lig <- grep("^LIG", rownames(G), value=TRUE)
  lig.targets <- names(which(Matrix::colSums(G[lig,])>0))
  lig <- unique(gsub(".*:|[ ].*","",lig))    
  kin <- grep("^KINASE_", rownames(G), value=TRUE)
  kin.targets <- names(which(Matrix::colSums(G[kin,])>0))
  kin <- unique(gsub(".*:|[ ].*","",kin))
  tf <- grep("^TF_", rownames(G), value=TRUE)
  tf.targets <- names(which(Matrix::colSums(G[tf,])>0))
  tf <- unique(gsub(".*:|[ ].*","",tf))
  mir <- grep("^MIR_", rownames(G), value=TRUE)
  mir.targets <- names(which(Matrix::colSums(G[mir,])>0))
  mir <- unique(gsub(".*:|[ ].*","",mir))
  mir <- toupper(sub("-","",gsub("^hsa-|^mmu-|-5p$|-3p$","",mir)))
  
  dim(X)
  gg <- rownames(X)[which(!rownames(X) %in% c(lig,tf,kin,mir))]
  data <- list(
    kin = X[rownames(X) %in% kin,],
    ##tf = X[rownames(X) %in% tf,],
    tf = X[rownames(X) %in% intersect(kin.targets,tf),],            
    ## mir = X[rownames(X) %in% mir,],
    gx = X[gg,]
  )
  return(data)
}


#' @export
mofa.exampledata <- function(dataset="geiger", ntop=2000,
                             scale.views=TRUE, omx.tissue=NULL) {

  data <- NULL
  if(dataset == "geiger") {
    dir = "~/Playground/opg-exampledata/metabolomics-kegg"
    counts  <- read_counts(file.path(dir, "multiomics-counts.csv"))
    samples <- read_samples(file.path(dir, "multiomics-samples.csv"))
    X <- logCPM(counts)
    X <- X + 1e-3*matrix(rnorm(length(X)),nrow(X),ncol(X))
    X <- head( X[order(-apply(X,1,sd)),], 1000)
    ##X <- t(scale(t(X)))
    dim(X)
    data <- list(
      px = X[grep("px:",rownames(X)),],
      mx = X[grep("mx:",rownames(X)),]
    )
    samples <- samples
    contrasts <- as.matrix(samples[,"activated",drop=FALSE])
    colnames(contrasts)[1] = "act_vs_notact"
    rownames(contrasts) = rownames(samples)
  }

  if(dataset == "brca") {
    library(mixOmics)
    data(breast.TCGA)
    data <- list(
      mir = t(breast.TCGA$data.train$mirna),
      gx = t(breast.TCGA$data.train$mrna),
      px = t(breast.TCGA$data.train$protein)
    )
    for(i in 1:3) rownames(data[[i]]) <- paste0(names(data)[i],":",rownames(data[[i]]))
    samples <- data.frame( sample = colnames(data[[1]]),
      condition = breast.TCGA$data.train$subtype )
    rownames(samples) <- colnames(data[[1]])
    c1=c2=as.character(samples$condition)
    c1[c1=="Her2"] <- NA
    c2[c2=="LumA"] <- NA
    contrasts <- cbind("LumA_vs_Basal"=c1, "Her2_vs_Basal"=c2)
    rownames(contrasts) <- rownames(samples)
  }

  if(dataset == "cll") {
    library(MOFAdata)
    utils::data("CLL_data")
    data <- CLL_data
    lapply(data,dim)
    rownames(data$mRNA) <- convert_probetype("Human",rownames(data$mRNA),"SYMBOL")    
    names(data) <- c("dr","me","gx","mu")
    kk <- sort(colnames(data[[1]]))
    for(i in 1:length(data)) {
      data[[i]] <- data[[i]][,kk]
      data[[i]][is.nan(data[[i]])] <- NA
      ## collapse duplicates
      data[[i]] <- rowmean(data[[i]], rownames(data[[i]]))
      data[[i]][is.nan(data[[i]])] <- NA
      rownames(data[[i]]) <- paste0(names(data)[i],":",rownames(data[[i]]))
    }
    data <- data[c("mu","me","gx","dr")]
    sapply(data, function(x) mean(apply(x,1,sd,na.rm=TRUE)))
    ##data$Drugs <- 0.01*data$Drugs
    CLL_metadata <- data.table::fread("ftp://ftp.ebi.ac.uk/pub/databases/mofa/cll_vignette/sample_metadata.txt")
    samples <- data.frame(CLL_metadata, row.names= CLL_metadata[["sample"]])
    samples <- samples[kk,]
    c1 <- c("no","yes")[1+samples$IGHV]
    c2 <- c("no","yes")[1+samples$trisomy12]
    contrasts <- cbind("IGHV:yes_vs_no"=c1, "trisomy12:yes_vs_no"=c2)
    rownames(contrasts) = rownames(samples)
  }

  if(dataset == "omx") {    
    load("~/Playground/omxp5/omxdata/drugs/ctrpv2-bytissue-omx.rda",verbose=TRUE)
    names(omx)
    data <- omx$level[[1]]$mat
    data <- data[c("mt","cn","me","gx")]
    names(data)
    
    lapply(data, dim)
    names(data)
    kk <- Reduce(intersect, lapply(data, colnames))
    data <- lapply(data, function(x) x[,kk])
    lapply(data, dim)
    names(data)
    for(i in 1:length(data)) {
      rownames(data[[i]]) <- paste0(names(data)[i],":",rownames(data[[i]]))
    }
    dim(omx$pheno)
    samples <- omx$pheno[kk,]

    if(!is.null(omx.tissue)) {
      sel <- grep(omx.tissue,samples$tissue)
      samples <- samples[sel,]
      data <- lapply(data, function(d) d[,rownames(samples)])
    }
    samples$tissue <- NULL
    
    dd <- colnames(samples)
    if(!is.null(omx.tissue)) dd <- sort(grep(omx.tissue,dd,value=1))
    dd <- grep("BRD",dd,value=1,invert=TRUE)
    samples <- samples[,dd]
    colnames(samples) <- sub("@","_",colnames(samples))
    sel <- head(order(-colSums(!is.na(samples))),10)
    contrasts <- samples[,sel]
    colnames(contrasts) <- paste0(colnames(contrasts),":1_vs_0")
  }

  if(dataset == "tf") {
    pgx <- playdata::GEIGER_PGX
    X <- pgx$X
    data <- mofa.splitByGeneRole(pgx$X) 
    samples <- pgx$samples
    contrasts <- pgx$contrasts
  }

  if(dataset == "tf") {
    pgx <- playdata::GEIGER_PGX
    X <- pgx$X
    data <- mofa.splitByGeneRole(pgx$X) 
    samples <- pgx$samples
    contrasts <- pgx$contrasts
  }

  if(is.null(data)) stop("could not get dataset")
  dbg("[mofa.exampledata] names(data) = ", names(data))
  
  ## strip prefix, convert to ASCII
  for(i in 1:length(data)) {
    rownames(data[[i]]) <- sub(".*:","",rownames(data[[i]]))
    rownames(data[[i]]) <- iconv(rownames(data[[i]]), "latin1", "ASCII", sub="")
  }

  ## filter out and reduce
  data <- lapply(data, function(d) d[rowMeans(is.na(d))<0.5,])
  data <- lapply(data, function(d) d[rowMeans(d**2,na.rm=TRUE)>0,])
  
  list( X = data, samples = samples, contrasts=contrasts )
}


##======================================================================
##========================== MIXOMICS ==================================
##======================================================================

## mixOmics related functions


#' @export
mixomics.compute_diablo <- function(xdata, y, ncomp=5) {

  # for square matrix filled with 0.1s
  design = matrix(0.1, ncol = length(xdata), nrow = length(xdata), 
                  dimnames = list(names(xdata), names(xdata)))
  diag(design) = 0 # set diagonal to 0s
  
  xx <- xdata
  xx <- mofa.prefix(xx)
  xx <- lapply( xx, t)
  Y <- y

  ## ncomp <- 5  
  nfeat <- 20

  keepX = lapply(names(xx), function(dt) rep(nfeat,ncomp))
  names(keepX) <- names(xx) 

  model <- mixOmics::block.splsda(
    X = xx,
    Y = Y,
    keepX = keepX,
    ncomp = ncomp,
    design = design) 
  #model = mixOmics::mixOmics(X, Y.factor, ncomp = 2, keepX = c(20, 20))

  dt <- names(xx)
  F <- model$variates$Y
  ww <- model$loadings[dt]
  W <- do.call(rbind, ww)

  list(
    model = model,
    scores = F,
    loadings = W
  )
}

#' @export
mixomics.plot <- function(model, style, comp=1, cutoff = 0.7, ...) {
  
  dt <- names(model$X)
  nt <- length(dt)
  block.colors <- rainbow(nt)

  if(style == "correlation") {
    mixOmics::plotDiablo(model, ncomp = 1)
  }

  if(style == "indiv") {        
    mixOmics::plotIndiv(model, ind.names = FALSE, legend = TRUE, 
                        title = 'DIABLO Sample Plots')
  }

  if(style == "arrow") {      
    mixOmics::plotArrow(model, ind.names = FALSE, legend = TRUE, 
                        title = 'DIABLO')
  }

  if(style == "var") {    
    mixOmics::plotVar(model, var.names = FALSE, 
                      style = 'graphics', legend = TRUE,
                      pch = c(16, 17, 15), cex = c(1,1,1)*1.2, 
                      col = block.colors)
  }

  if(style == "circos") {    
    mixOmics::circosPlot(model, cutoff = cutoff, line = TRUE,
                         color.blocks= block.colors,
                         color.cor = c("black","red"),
                         size.labels = 1.5)
  }

  if(style == "network") {  
    mixOmics::network(model, blocks = dt,
                      color.node = block.colors,
                      size.node = 0.04,
                      cutoff = cutoff)
  }

  if(style == "loadings") {
    plotLoadings(model, comp = comp, contrib = 'max', method = 'median')
  }
  
  if(style == "cim") {
    cimDiablo(model, mar=c(8,20), transpose=TRUE,
              legend.position='bottomright', size.legend=1)
  }
    
}

#' @export
mixomics.predict <- function(model, newdata) {
  newdata <- mofa.scale_views(newdata)
  newdata <- lapply(mofa.prefix(newdata),t)
  for(k in 1:length(newdata)) {
    dtype <- names(newdata)[k]
    if(dtype %in% names(model$X)) {
      mat <- newdata[[k]]
      jj <- match( colnames(model$X[[dtype]]), colnames(mat))
      newdata[[k]] <- mat[,jj]
    } else {
      newdata[[k]] <- NULL
    }
  }
  names(newdata)  
  pred <- mixOmics:::predict.block.spls(model, newdata = newdata)    
  predicted = pred$WeightedVote$centroids.dist[,2]
  predicted
}


#' @describeIn mixHivePlot function generates a hive plot visualization of variable loadings
#' from a lmer model result object.
#' @export
mixPlotLoadings.DEPRECATED <- function(res, showloops = FALSE, cex = 1) {
  cat("<mixPlotLoadings> called\n")
  levels <- levels(res$Y)
  ny <- length(levels)
  klrpal <- c("blue2", "orange2")
  klrpal <- rep(RColorBrewer::brewer.pal(n = 8, "Set2"), 10)[1:ny]

  names(klrpal) <- levels

  plotly::layout(matrix(1:6, 1, 6), widths = c(1, 0.5, 1, 0.5, 1, 0.5))
  k <- 1
  for (k in 1:3) {
    W <- res$W[[k]]
    graphics::par(mar = c(5, 8 * cex, 4, 0), mgp = c(2.2, 0.8, 0))
    graphics::barplot(t(W),
      horiz = TRUE, las = 1,
      border = NA, col = klrpal,
      names.arg = sub(".*:", "", rownames(W)),
      xlim = c(0, 1.1) * max(rowSums(W, na.rm = TRUE)),
      cex.axis = 1 * cex, cex.names = 1.1 * cex,
      cex.lab = 1 * cex, xlab = "importance"
    )
    graphics::title(names(res$loadings)[k],
      cex.main = 1.3 * cex,
      adj = 0.33, xpd = NA
    )
    graphics::legend("topright",
      legend = names(klrpal),
      cex = 1.1 * cex, pch = 15, col = klrpal, #
      y.intersp = 0.85, inset = c(0.15, 0.03)
    )

    if (k < 99) {
      ## add correlation lines
      graphics::par(mar = c(5, 0, 4, 0))

      plot(0,
        type = "n", xlim = c(0, 1), ylim = c(0, nrow(W)),
        xaxt = "n", yaxt = "n", bty = "n", xlab = ""
      )

      g1 <- rownames(res$W[[k]])
      g2 <- rownames(res$W[[ifelse(k < 3, k + 1, 1)]])

      sel <- which(res$edges[, "from"] %in% c(g1, g2) &
        res$edges[, "to"] %in% c(g1, g2))
      sel
      ee <- res$edges[sel, ]
      ii <- apply(ee[, 1:2], 1, function(e) which(e %in% g1))
      jj <- apply(ee[, 1:2], 1, function(e) which(e %in% g2))
      ee$from <- res$edges[sel, ][cbind(1:nrow(ee), ii)]
      ee$to <- res$edges[sel, ][cbind(1:nrow(ee), jj)]


      lwd <- ee$importance

      lwd <- rank(abs(lwd), na.last = "keep")**1.5
      lwd <- 3.0 * cex * (lwd / max(lwd))
      lty <- 1 + 1 * (sign(ee$rho) < 0)
      xy <- cbind(match(ee$from, g1), match(ee$to, g2))
      xy[, 2] <- (xy[, 2] - 0.5) / length(g2) * length(g1)
      klr <- rep(psych::alpha("grey70", 0.3), nrow(ee))
      if (showloops) {
        klr <- rep(psych::alpha("grey70", 0.2), nrow(ee))
        klr[which(ee$looping)] <- psych::alpha("red3", 0.3)
      }
      graphics::segments(0, xy[, 1] - 0.5, 1, xy[, 2], lwd = lwd, col = klr, lty = lty)
      rr <- paste(round(range(abs(ee$rho)), 2), collapse = ",")

      graphics::title(sub = paste0("[", rr, "]"), line = -1.2, cex.sub = cex)
    }
  }
}



##===================================================================
##======================= SNF =======================================
##===================================================================

#' @export
snf.cluster <- function( xx, pheno=NULL, plot=TRUE ) {

  has.missing <- sapply(xx, function(x) sum(is.na(x))>0)
  has.missing
  xx <- lapply( xx, function(x) svdImpute2(x))
  
  Data <- lapply(xx, t)
  Dist <- list()
  dt=names(xx)[1]
  for(dt in names(Data)) {
    x <- Data[[dt]]
    x <- SNFtool::standardNormalization(x) ## cannot handle NA
    Dist[[dt]] <- SNFtool::dist2(x,x)^(1/2)
  }
  
  ## next, construct similarity graphs
  ## First, set all the parameters:
  K = min(ncol(xx[[1]])/4, 15);         # number of neighbors, usually (10~30)
  alpha = 0.5;                          # hyperparameter, usually (0.3~0.8)
  message("clusterSNF: K = ", K)
  message("clusterSNF: alpha = ", alpha)  
  Wlist = lapply( Dist, function(d) SNFtool::affinityMatrix(d, K=K, alpha))
    
  ## next, we fuse all the graphs
  ## then the overall matrix can be computed by similarity network fusion(SNF):
  W = SNFtool::SNF( Wlist, K=K, t=100)
  
  ## With this unified graph W of size n x n, 
  ## you can do either spectral clustering or Kernel NMF. 
  ## If you need help with further clustering, please let us know. 
  
  Data2 <- c(Dist, list(SNF=W))
  lapply(Data2,dim)
  ##posx <- lapply( Data2, function(x) Rtsne::Rtsne(t(x), is_distance=TRUE, perplexity=15)$Y)
  posx <- lapply( Data2, function(x)
    Rtsne::Rtsne(t(x), perplexity=K, check_duplicates=FALSE)$Y)
  ##posx <- lapply( Data2, function(x) uwot::umap(t(x), n_neighbors=15))
  for(i in 1:length(posx)) rownames(posx[[i]]) <- colnames(xx[[1]])
  rownames(W) = colnames(W) = rownames(posx[[1]])

  ## -------------- graph & louvain ----------------
  gr <- igraph::graph_from_adjacency_matrix(
    W, mode="undirected", diag=FALSE, weighted=TRUE)
  cl <- igraph::cluster_louvain(gr)$membership
  cl <- paste0("SNF",cl)
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
snf.plot_affinity <- function(snf, k=1, par=TRUE) {
  i=1
  mat <- c( snf$affinityMatrix, list(SNF=snf$W))
  n <- ceiling(sqrt(length(mat)))
  if(par) par(mfrow=c(n,n), mar=c(6,1,2,7))
  for(i in 1:length(mat)) {
    x <- mat[[i]]
    if(k!=1) x <- sign(x) * x^k
    diag(x) <- 0
    if(nrow(x) <= 20) { 
      gx.imagemap( x )
    } else {
      rownames(x) <- NULL
      colnames(x) <- NULL
      gx.imagemap( x )
      mtext("samples", side=1, line=0.2)
      mtext("samples", side=4, line=0.2)
    }
    title( names(mat)[i], cex.main=1.6 )
  }
}

#' @export
snf.plot_clustering <- function(snf, y=NULL, par=TRUE, ...) {
  i=1
  n <- ceiling(sqrt(length(snf$posx)))
  if(par) par(mfrow=c(n,n), mar=c(6,1,2,7))
  cc <- "grey30"
  if(!is.null(y)) cc <- factor(y)
  for(i in 1:length(snf$posx)) {
    plot( snf$posx[[i]], col=cc, pch=19, cex=0.8,
         xlab = "TSNE-x", ylab = "TSNE-y", ... )
    title( names(snf$posx)[i], cex.main=1.6 )
  }
}

#' @export
snf.plot_graph <- function(snf, min.rho=0.5, plot=TRUE,
                           val=NULL, label=NULL,
                           vlabcex = 1,
                           ...) {
  
  A <- snf$affinityMatrix
  W <- snf$W
  diag(W) <- 0
  W <- W / max(W)
  for(i in 1:length(A)) {
    diag( A[[i]] ) <- 0
    A[[i]] <- A[[i]] / max(A[[i]])
    ##  A[[i]] <- A[[i]] * (A[[i]] > 0.5)
  }
  
  rownames(W) = colnames(W) = rownames(snf$posx[[1]])
  if(!is.null(label)) {
    names(label) <- rownames(W)
  }
  ee <- which(W > min.rho, arr.ind=TRUE)
  ee <- ee[ which(ee[,2] > ee[,1]), ]
  af <- unlist(lapply(A, function(a) a[ee] ))
  af <- af / max(af,na.rm=TRUE)
  tp <- as.vector(sapply(names(A), function(a) rep(a, nrow(ee))))
  df <- data.frame(ee, type=tp, affinity=af)
  head(df)
  df[,1] <- rownames(W)[df[,1]]
  df[,2] <- rownames(W)[df[,2]]
  gr <- igraph::graph_from_edgelist(
    as.matrix(df[,1:2]), directed=FALSE)
  
  colors <- paste0( rainbow(length(A)), "55")
  colors <- head(rep(paste0( palette()[-1], "66"),99),length(A))

  igraph::V(gr)$name
  igraph::E(gr)$type <- df$type
  igraph::E(gr)$weight <- (1*df$affinity)**1.1
  igraph::E(gr)$color <- colors[factor(df$type)]
  igraph::V(gr)$label <- igraph::V(gr)$name
  if(!is.null(label)) {
    igraph::V(gr)$label <- label[igraph::V(gr)$name]
  }

  vcolors <- rev(rainbow(8))
  if(is.null(val)) {
    val <- factor(snf$cluster)
    names(val) <- igraph::V(snf$graph)$name
    val <- val[igraph::V(gr)$name]
  }
  igraph::V(gr)$color <- vcolors[factor(val)]
    
  if(plot) {
    par(mfrow=c(1,1), mar=c(1,1,1,1))
    dtypes <- unique(factor(df$type))
    plot(gr,
         vertex.size = 0,         
         vertex.label = igraph::V(gr)$label,
         vertex.label.family = "sans",
         vertex.label.cex = 1.2*vlabcex,
         # vertex.label.color = igraph::V(gr)$color,
         edge.width = 6 * igraph::E(gr)$weight,
         edge.color = igraph::E(gr)$color,    
         ...
         )
    legend("bottomright", legend=dtypes, fill=colors,
           cex=0.9, y.intersp=0.85 )
  } else {
    return(gr)
  }
}

#' @export
snf.heatmap <- function(snf, X, samples, nmax=50) {

  rX <- X[order(-apply(X,1,sd)),]
  dt <- sub(":.*","",rownames(X))
  n <- nmax / length(unique(dt))
  sel <- tapply(rownames(X),dt,function(a) head(a,n))
  sel <- unlist(sel)
  dt <- sub(":.*","",sel)
  gx.splitmap( rX[sel,], splitx=snf$cluster,
              split = dt, nmax=160, cexRow=0.8, 
              col.annot = samples,
              dist.method = "euclidean",
              softmax=TRUE, show_key=FALSE)

}

##======================================================================
##====================== LASAGNA =======================================
##======================================================================


#' @export
plot_lasagna <- function(posx, vars=NULL, num_edges=20) {
    
  for(i in 1:length(posx)) {
    a <- names(posx)[i]
    if(!all(grepl(a,names(vars[[i]]))) ) {
      rownames(posx[[i]]) <- paste0(a,":",rownames(posx[[i]]))
      if(!is.null(vars) && is.list(vars)) {
        names(vars[[i]]) <- sub(paste0(a,":"),"",names(vars[[i]]))
        names(vars[[i]]) <- paste0(a,":",names(vars[[i]]))        
      }
    }
  }

  ## Sample planes across datatyps
  colx = 'black'

  ## feature maps across datatypes
  x <- data.frame()
  for(i in names(posx)) {
    x <- rbind(x, data.frame( posx[[i]], type=i))
  }
  x$type <- factor(x$type, levels=names(posx))
    
  if(!is.null(vars)) {
    colorsx <- gplots::colorpanel(255,low="blue3",mid="grey80",high="red3")
    if(!is.list(vars) && is.vector(vars)) {
      vars <- lapply( posx, function(p) vars[rownames(p)] )
    }
    vars <- lapply(vars, function(x) (x / max(abs(x))))
    colx <- unlist(lapply(vars, function(v) colorsx[128+ceiling(v*127)]))
    dt <- names(posx)
    dt
    i=1
    ee <- c()
    for(i in 1:(length(dt)-1)) {
      v1 <- vars[[i]][rownames(posx[[i]])]
      v2 <- vars[[i+1]][rownames(posx[[i+1]])]
      ii <- head(names(sort(-v1)),num_edges)
      jj <- head(names(sort(-v2)),num_edges)
      ee1 <- cbind(ii,jj)
      ee <- rbind(ee, ee1)
      ii <- head(names(sort(v1)),num_edges)
      jj <- head(names(sort(v2)),num_edges)
      ee2 <- cbind(ii,jj)
      ee <- rbind(ee, ee2)
    }
  }

  segment_mat <- cbind(
    match(ee[,1],rownames(x)),
    match(ee[,2],rownames(x))
  )
  segment_mat
  nseg <- nrow(segment_mat)

  grimon::grimon(
    x = x,
    label = names(posx),
    col = colx,
    format = "long",
    segment_mat = segment_mat[,],
    segment_col = colx[segment_mat[,1]],
    point_size = 4,
    z_interval = 2,
    optimize_coordinates = TRUE,
    maxiter = 1e3,
    score_function = "angle",
    #return_coordinates = TRUE,
    #plot_2d_panels = FALSE,
    segment_alpha = 0.8)

}

#'
#' 
#' @export
plotly_lasagna <- function(posx, vars=NULL, num_edges=20) {
  
  for(i in 1:length(posx)) {
    a <- names(posx)[i]
    if(!all(grepl(a,names(vars[[i]]))) ) {
      rownames(posx[[i]]) <- paste0(a,":",rownames(posx[[i]]))
      if(!is.null(vars)) {
        names(vars[[i]]) <- paste0(a,":",names(vars[[i]]))
      }
    }
  }

  ## scale each layer
  posx <- lapply( posx, uscale )
  
  ## feature maps across datatypes
  x <- data.frame()
  for(i in names(posx)) {
    x <- rbind(x, data.frame( posx[[i]], type=i))
  }
  x$type <- factor(x$type, levels=names(posx))
  colnames(x) <- c("x","y","z")

  znames = c("mx" = "Metabolomics",
             "gx" = "Transcriptomics",
             "mir" = "micro-RNA",                           
             "px" = "Proteomics")
  
  plotlyLasagna(x, znames = znames)
  
}


#' 
#' @export
lasagna.create_model <- function(data, pheno="pheno", ntop=1000, nc=-1,
                                 use.gmt=TRUE, use.graphite=TRUE) {

  names(data)
  if(pheno == "pheno") {
    Y <- expandPhenoMatrix(data$samples, drop.ref=FALSE)
  } else {
    Y <- makeContrastsFromLabelMatrix(data$contrasts)
    revY <- -Y
    colnames(revY) <- reverse.AvsB(colnames(Y))
    Y <- cbind(Y, revY)
  }
  data$X$PHENO <- t(Y)
  names(data)

  xx <- data$X
  if(!is.null(ntop) && ntop>0) {
    xx <- lapply(xx, function(x) head(x[order(-apply(x,1,sd)),,drop=FALSE],ntop))
  }
  xx <- mofa.prefix(xx)
  X <- do.call(rbind, xx)
  remove(xx)
  
  ## add SOURCE/SINK
  X <- rbind(X, SOURCE=1, SINK=1)
  suppressWarnings( R <- cor(t(X), use="pairwise") )
  R[is.na(R)] <- 0.01
  
  ## mask for metabolics PPI
  if(use.graphite) {
    xtypes <- setdiff(names(data$X),"PHENO")
    xtypes
    has.mx <- ("mx" %in% xtypes)
    has.px <- any(c("gx","px") %in% xtypes)
    if(has.mx && has.px) {
      message("using GRAPHITE PPI structure")
      load("~/Playground/public-db/pathbank.org/GRAPHITE_PPI.rda",verbose=TRUE)
      head(GRAPHITE_PPI)
      GRAPHITE_PPI[,1] <- paste0("px:",GRAPHITE_PPI[,1])
      GRAPHITE_PPI[,2] <- paste0("px:",GRAPHITE_PPI[,2])
      GRAPHITE_PPI[,1] <- sub("px:CHEBI:","mx:",GRAPHITE_PPI[,1])
      GRAPHITE_PPI[,2] <- sub("px:CHEBI:","mx:",GRAPHITE_PPI[,2])
      gr <- igraph::graph_from_edgelist( as.matrix(GRAPHITE_PPI[,1:2]), directed="FALSE")
      A <- igraph::as_adjacency_matrix(igraph::simplify(gr))
      dt <- sub(":.*","",rownames(R))
      i1 <- which( dt %in% c('gx','px'))
      i2 <- which( dt == 'mx')
      a1 <- match( rownames(R)[i1], rownames(A))
      a2 <- match( rownames(R)[i2], rownames(A))
      i1 <- i1[!is.na(a1)]
      a1 <- a1[!is.na(a1)]
      i2 <- i2[!is.na(a2)]
      a2 <- a2[!is.na(a2)]
      R[i1, i2] <- R[i1, i2] * as.matrix(A[a1, a2])
      R[i2, i1] <- R[i2, i1] * as.matrix(A[a2, a1])
    }
  }

  ## mask for GSETS/pathways ???
  if(use.gmt) {
    
  }
  
  ## mask for layer interaction
  dt <- sub(":.*","",rownames(R))
  table(dt)
  names(data$X)  
  layers <- c("SOURCE",names(data$X),"SINK")
  layer_mask <- as.matrix(R) * 0
  for(i in 1:(length(layers)-1)) {
    ii <- which( dt == layers[i] )
    jj <- which( dt == layers[i+1] )
    layer_mask[ii,jj] <- 1
    layer_mask[jj,ii] <- 1
  }
  R <- R * layer_mask 

  ## Reduce connections
  if(!is.null(nc) && nc > 0) {
    message(paste("reducing edges to maximum", nc, "connections"))
    ## NEED CHECK!!! SOMETHING WRONG.
    xtypes <- setdiff(names(data$X),"PHENO")
    reduce_mask <- matrix(1, nrow(R), ncol(R))
    i=1
    for(i in 1:(length(xtypes)-1)) {
      ii <- which( dt == xtypes[i] )
      jj <- which( dt == xtypes[i+1] )
      R1 <- R[ii,jj,drop=FALSE]
      rii <- apply(abs(R1), 1, function(r) tail(sort(r),nc)[1])
      rjj <- apply(abs(R1), 2, function(r) tail(sort(r),nc)[1])    
      rr <- abs(R1) >= rii | t( t(abs(R1)) >= rjj)
      reduce_mask[ii,jj] <- rr
      reduce_mask[jj,ii] <- t(rr)
    }
    R <- R * reduce_mask 
  }

  ## create graph
  gr <- igraph::graph_from_adjacency_matrix(
    R, diag=FALSE, weighted=TRUE, mode="undirected")

  ii <- grep("^mx:", igraph::V(gr)$name)
  if(length(ii)) {
    message(paste("translating metabolite ID to names"))
    mx.id <- sub("mx:","", igraph::V(gr)$name[ii])
    mx.annot <- getMetaboliteAnnotation(mx.id)
    mx.names <- paste0("mx:",mx.annot$gene_title," (",mx.id,")")
    igraph::V(gr)$name[ii] <- mx.names
    rownames(X)[ii] <- mx.names    
  }

  ##R1 <- igraph::as_adjacency_matrix(gr, attr='weight')
  list(
    graph = gr,
    X = X,
    Y = Y,
    layers = layers
  )
}


#' 
#' @export
lasagna.solve_SP <- function(obj, pheno, rm.neg=TRUE, vtop=200) {
  
  ## condition with phenotype correlation
  PHENO = paste0("PHENO:",pheno)  
  if( is.character(pheno) && PHENO[1] %in% rownames(obj$X)) {
    pheno <- obj$X[PHENO,]
  }
  suppressWarnings( rho <- cor(t(obj$X), pheno)[,1] )
  rho[is.na(rho)] <- 0.1  ## ???
  names(rho) <- rownames(obj$X)

  tail(sort(rho))
  head(sort(rho))  
  F <- outer(rho, rho)
  F <- sqrt(abs(F))
  dim(F)
  
  ## weight edges with foldchange on nodes
  gr <- obj$graph
  R <- igraph::as_adjacency_matrix(gr, attr='weight')
  F <- F[rownames(R),colnames(R)]
  R <- R * F
  gr <- igraph::graph_from_adjacency_matrix(
    R, diag=FALSE, weighted=TRUE, mode="undirected")
  gr <- igraph::simplify(gr)

  if(rm.neg) {
    message("removing negative edges")
    gr <- igraph::delete_edges(gr, edges = which(igraph::E(gr)$weight <= 0))
  }
  
  ## compute all(?) shortest paths
  wt <- igraph::E(gr)$weight
  max.wt <- max(wt,na.rm=TRUE)
  ##wt <- (max.wt - wt)**1
  ##wt <- (max.wt / pmax(wt, 1e-2) - 1)**1
  wt <- -log(wt + 1e-4)

  ## resctrict search to top logFC nodes
  v.nodes <- igraph::V(gr)
  if(vtop > 0) {
    v.names <- igraph::V(gr)$name  
    v.types <- sub(":.*","",v.names)
    vtop.nodes <- tapply( v.names, v.types, function(v)
      head(names(sort(-abs(rho[v]))), vtop))
    v.nodes <- unlist(vtop.nodes)
  }

  suppressWarnings({
    s1 <- igraph::shortest_paths(gr, from="SOURCE", to=v.nodes, weights = wt,
                                 output='both')
    s2 <- igraph::shortest_paths(gr, from="SINK", to=v.nodes, weights = wt,
                                 output='both')
  })
  
  p1 <- sapply(s1$epath, function(e) sum(wt[e]))
  p2 <- sapply(s2$epath, function(e) sum(wt[e]))    
  sp.score <- p1 + p2
  names(sp.score) <- v.nodes
  
  v1 <- lapply(s1$vpath, function(i) igraph::V(gr)$name[i])
  v2 <- lapply(s2$vpath, function(i) rev(igraph::V(gr)$name[i]))
  vv <- mapply(union, v1, v2, SIMPLIFY=FALSE)
  names(vv) <- v.nodes

  ## check paths: this correct the shortest path
  splen <- sapply( vv, length )
  table(splen)
  splen0 <- min(splen[splen>0])
  ii <- which(splen > splen0)
  if(length(ii)) {
    rmdupv <- function(v) v[!duplicated(sub(":.*","",v))]
    vv[ii] <- lapply(vv[ii], function(v) rmdupv(v))
  }
  names(sp.score) <- v.nodes

  vv <- lapply(vv, function(v) setdiff(v, c("SOURCE","SINK")))  
  table(sapply(vv, length))
  vv <- vv[sapply(vv, length)>0]
  spath <- sapply(vv, paste, collapse="->")
  V <- data.frame(do.call(rbind, vv))

  ##  V$path <- spath
  rownames(V) <- names(vv)
  vtype <- sub(":.*","",V[1,])
  colnames(V) <- vtype
  head(V)
  
  S <- apply(V, 2, function(v) rho[v])
  rownames(S) <- rownames(V)
  sp.score <- sp.score[rownames(S)]
  
  min.sp <- min(sp.score, na.rm=TRUE)
  sp.score <- exp(-2*sp.score/min.sp)

  rho <- rho[rownames(S)]
  sp.score <- sp.score * sign(rho)
  
  sp <- data.frame(
    score = sp.score,
    rho = S,
    path = spath
  )
  rownames(sp) <- rownames(V)
  
  sp <- sp[grep("SOURCE|SINK",rownames(sp),invert=TRUE),,drop=FALSE]
  sp <- sp[order(-sp$score),]  
  sp
}

#' 
#' @export
lasagna.plot_SP <- function(sp, ntop=200, hilight=NULL, labcex=1,
                            colorby="pheno", plotlib="ggplot") {

  if(!is.null(ntop) && ntop>0) {
    sp <- sp[order(-sp$score),]
    sp <- head(sp, ntop)
  }
  fig <- NULL
  
  if(plotlib=="plotly") {
    dimensions <- list()
    this.layer <- "mx"
    rho.layers <- grep("rho",colnames(sp),value=TRUE)
    for(var in rho.layers) {
      dd <- list(range = c(-1,1),
                 label = var,
                 values = sp[,var])
      dimensions <- c(dimensions, list(dd))
    }
    
    if(colorby=="pheno") {
      line = list(
        color = ~rho.PHENO,
        colorscale = "RdBu",
        showscale = TRUE, cmin = -1, cmax=1)
    }
    if(colorby=="score") {
      line = list(
        color = ~score,
        colorscale = "RdBu",
        showscale = TRUE,
        cmin = min(sp$score,na.rm=TRUE),
        cmax = min(sp$score,na.rm=TRUE)
      )
    }
    
    fig <- sp %>% plotly::plot_ly(
      type = 'parcoords',
      line = line,
      dimensions = dimensions
    )
  }

  if(plotlib=="ggplot") {

    ##data <- iris
    rho.cols <- grep("rho",colnames(sp))
    rho.cols
    S <- sp[,rho.cols]
    sp$alpha <- 0.2
    if(!is.null(hilight)) {
      sp$alpha <- c(0.04, 1)[ 1 + 1*(rownames(sp) %in% hilight) ]
    }
    
    fig <- GGally::ggparcoord(
      sp,
      scale = "globalminmax",
      columns = rho.cols,
      groupColumn = "score",
      #order = "anyClass",
      showPoints = TRUE, 
      title = "LASAGNA Parallel Coordinates Plot"
    ) + geom_line(
      linewidth = 1,
      aes(linewidth = alpha, alpha = alpha) ) +
      scale_colour_gradient2() +
      scale_alpha(guide = 'none') +
      labs( x = "", y = "correlation   (rho)")
    fig

    if(!is.null(hilight)) {
      i = 10
      for(hi in hilight) {
        i <- match(hi, rownames(S))
        tx <- colnames(S)
        ty <- as.matrix(S)[i,]
        pp <- strsplit( sp$path[i], split="->")[[1]]
        cc <- c("white","yellow")[1 + 1*(pp == hi)]
        tt <- sub(":.*","",pp)
        pp <- pp[which(!duplicated(tt))]
        pp <- gsub(".*:|\\(.*","",pp)
        pp <- stringr::str_wrap(pp,15)
        fig <- fig + ggplot2::annotate("label", x=tx, y=ty, label=pp,
          fill=cc, size = (14*labcex/ ggplot2::.pt))
      }
    }
  }
  
  fig
}



##======================================================================
##======================================================================
##======================================================================


