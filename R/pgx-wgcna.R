##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##
##


#' @title WGCNA network construction and module detection
#'
#' @param pgx PGX object containing gene expression data
#' @param minmodsize Minimum module size cutoff
#' @param power Soft thresholding power for network construction
#' @param cutheight Cut height for module dendrogram
#' @param deepsplit Number of splits for module dendrogram
#' @param ngenes Number of genes to use (most variable)
#'
#' @return List containing WGCNA network and module results
#'
#' @description Constructs a weighted gene coexpression network and detects
#' modules using WGCNA on a PGX object.
#'
#' @details This function takes a PGX object containing gene expression data.
#' It constructs a weighted gene coexpression network using the WGCNA package.
#' Soft thresholding power is set by \code{power}.
#'
#' Modules are detected by cutting the module dendrogram at \code{cutheight} and
#' with \code{deepsplit} splits. Only the \code{ngenes} most variable genes are used.
#'
#' The output is a list containing the WGCNA network object and module results,
#' including module assignments, colors, and summary statistics.
#'
#' @export
pgx.wgcna <- function(
    pgx,
    minmodsize = 20,
    power = 12,
    cutheight = 0.15,
    deepsplit = 2,
    minKME = 0.3,
    networktype = "signed",
    tomtype = "signed",
    numericlabels = FALSE,
    ngenes = 4000,
    maxBlockSize = 5000,
    gset.filter = "PATHWAY|HALLMARK|^GO|^C[1-9]"
    ) {

  ##minmodsize=10;power=NULL;cutheight=0.15;deepsplit=2;ngenes=4000;networktype="signed";tomtype="signed";numericlabels=FALSE;ngenes=2000;gset.filter=NULL;minKME=0.8;maxBlockSize=5000

  samples <- pgx$samples
  ## no dot pheno
  samples <- samples[, grep("^[.]", colnames(samples), invert = TRUE), drop = FALSE]
  X <- pgx$X
  
  ## WGCNA does not work very well with scRNAseq due to sparsity.  To
  ## bypass the issue, hdWGCNA computes metacells.  Here we compute
  ## metacells too using Supercell,
  is.singlecell <- !is.null(pgx$datatype) && tolower(pgx$datatype) %in% c("scrnaseq","scrna-seq")
  is.toobig <- ncol(X) > 2000
  if (is.singlecell && is.toobig) {
    message("[pgx.wgcna] WGCNA: scRNAseq. >2K cells. Computing supercells with SuperCell.")
    counts <- pmax(2**X - 1,0)
    group <- samples[, "celltype"]
    nb <- round(ncol(counts)/1000)
    message("[pgx.wgcna] running SuperCell. nb = ", nb)    
    sc <- pgx.supercell(counts, samples, group = group, gamma = nb)
    message("[pgx.wgcna] SuperCell done: ", ncol(counts), " ->", ncol(sc$counts))
    message("[pgx.wgcna] Normalizing supercell matrix (logCPM)")
    X <- as.matrix(logCPM(sc$counts, total = 1e4, prior = 1))
    samples <- sc$meta
    remove(counts, group, sc)
    gc()
  }

  if(!is.null(pgx$datatype) &&  pgx$datatype == "multi-omics") {
    message("[pgx.wgcna] Performing multi-omics ComBat on datatype.")
    X <- playbase::normalizeMultiOmics(X, method = "combat")
  }
  
  message("[pgx.wgcna] start wgcna.compute...")
  wgcna <- wgcna.compute(
    X = X,
    samples = samples,
    minmodsize = minmodsize, # default: min(20,...)
    power = power, # default: 12 (for signed)
    cutheight = cutheight, # default: 0.15
    deepsplit = deepsplit, # default: 2
    minKME = minKME, # default: 0.30
    networktype = networktype, # default: unsigned (but signed is better...)
    tomtype = tomtype, # default: signed
    numericlabels = numericlabels,
    maxBlockSize = maxBlockSize,
    ngenes = ngenes
  )
  message("[pgx.wgcna] finished computing wgcna.compute!")
  
  ##---------------------------------------------------
  ## compute dimensionality reductions using TOM matrix
  ##---------------------------------------------------
  wTOM <- NULL
  if("TOM" %in% names(wgcna)) wTOM <- wgcna$TOM
  if(is.null(wTOM) && "svTOM" %in% names(wgcna)) wTOM  <- wgcna$svTOM %*% t(wgcna$svTOM)
  dissTOM <- 1 - wTOM
  rownames(dissTOM) <- colnames(dissTOM) <- colnames(wgcna$datExpr)
  clust <- pgx.clusterBigMatrix(dissTOM, methods = c("umap", "tsne", "pca"), dims = c(2))
  if ("cluster.genes" %in% names(pgx)) {
    posx <- pgx$cluster.genes$pos[["umap2d"]]
    clust[["umap2d"]] <- posx[colnames(wgcna$datExpr), ]
  }
  remove(dissTOM)
  
  ## ----------------------------------------------------
  ## Do geneset analysis
  ## ----------------------------------------------------
  res.gse <- wgcna.compute_enrichment(
    wgcna, pgx,
    method = c("fisher", "gsetcor", "mmcor"),
    ntop = 1000,
    filter = gset.filter
  )

  ## add to results object
  wgcna$gse <- res.gse$gsets
  wgcna$enriched_genes <- res.gse$genes
  wgcna$clust <- clust
  wgcna$networktype <- networktype
  wgcna$tomtype <- tomtype

  return(wgcna)
}

#' @export
wgcna.compute <- function(X,
                          samples,
                          minmodsize = 20,
                          power = 12,
                          cutheight = 0.15,
                          deepsplit = 2,
                          minKME = 0.3,
                          networktype = "signed",
                          tomtype = "signed",
                          reassignThreshold = 1e-6,
                          ngenes = 2000,
                          numericlabels = FALSE,
                          maxBlockSize = 8000,
                          merge.dendro = TRUE,
                          prefix = "ME",
                          sv.dim = 80,
                          verbose = 0
                          ) {

  ##minmodsize=20;power=12;cutheight=0.15;deepsplit=2;ngenes=2000;networktype="signed";tomtype="signed";numericlabels=FALSE;prefix="ME";minKME=0.3;reassignThreshold=1e-6;maxBlockSize=5000

  require(WGCNA)
  
  nmissing <- sum(is.na(X))
  if (nmissing > 0) {
    message("Found ", nmissing, " missing values in X. Imputing prior to WGCNA.")
    X <- svdImpute2(X)
  }

  X <- X[!duplicated(rownames(X)), ]

  if (ngenes > 0 && nrow(X) > ngenes) {
    sdx <- matrixStats::rowSds(X, na.rm = TRUE)
    X <- X[sdx > 0.1 * mean(sdx, na.rm = TRUE), ] ## filter low SD??
    X <- X[order(-matrixStats::rowSds(X, na.rm = TRUE)), ]
    X <- utils::head(X, ngenes)
  }

  message("[wgcna.compute] dim(X) = ", paste(dim(X), collapse = " x "))
  message("[wgcna.compute] dim(samples) = ", paste(dim(samples), collapse = " x "))

  datExpr <- t(X)

  if (is.null(power)) {
    ## Estimate best power
    message("[wgcna.compute] estimating optimal power...")
    powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
    sft <- WGCNA::pickSoftThreshold(
      datExpr,
      powerVector = powers,
      networkType = networktype,
      verbose = verbose
    )
    power <- sft$powerEstimate
    if (is.na(power)) power <- 20
  }

  ## adapt for small datasets (also done in WGCNA package)
  minmodsize <- min(minmodsize, ncol(datExpr) / 2)

  message("[wgcna.compute] minmodsize = ", minmodsize)
  message("[wgcna.compute] number of features = ", nrow(X))
  message("[wgcna.compute] minKME = ", minKME)
  message("[wgcna.compute] power = ", power)
  message("[wgcna.compute] mergeCutHeight = ", cutheight)

  message("[wgcna.compute] computing blockwise Modules...")
  ##WGCNA::enableWGCNAThreads()
  cor <- WGCNA::cor ## needed...
  net <- WGCNA::blockwiseModules(
    datExpr,
    power = power,
    networkType = networktype,
    TOMType = tomtype,
    minModuleSize = minmodsize,
    reassignThreshold = reassignThreshold,
    mergeCutHeight = cutheight,
    minKMEtoStay = minKME,
    numericLabels = numericlabels, ## numeric or 'color' labels
    deepSplit = deepsplit,
    maxBlockSize = maxBlockSize,
    verbose = verbose
  )
  cor <- stats::cor
  message("[wgcna.compute] finished blockwise Modules...")
    
  ##  prefix="ME"
  table(net$colors)
  names(net$MEs) <- sub("^ME", prefix, names(net$MEs))
  net$labels <- paste0(prefix, net$colors)

  ## clean up traits matrix
  datTraits <- data.frame(samples, check.names = FALSE)
  isdate <- apply(datTraits, 2, is.Date)
  datTraits <- datTraits[, !isdate, drop = FALSE]

  ## Expand multi-class discrete phenotypes into binary vectors
  datTraits <- utils::type.convert(datTraits, as.is = TRUE)
  tr.class <- sapply(datTraits, class)
  sel1 <- which(tr.class %in% c("factor", "character"))
  sel2 <- which(tr.class %in% c("integer", "numeric"))
  tr1 <- datTraits[, 0]
  tr2 <- datTraits[, 0]
  if (length(sel1)) {
    tr1 <- expandPhenoMatrix(datTraits[, sel1, drop = FALSE], drop.ref = FALSE)
    if (is.null(tr1)) tr1 <- datTraits[, 0]
  }
  if (length(sel2)) {
    ## keeping numeric phenotypes
    tr2 <- datTraits[, sel2, drop = FALSE]
  }
  datTraits <- cbind(tr1, tr2)

  ## list of genes in modules
  me.genes <- tapply(names(net$colors), net$colors, list)
  names(me.genes) <- paste0(prefix, names(me.genes))
  me.genes <- me.genes[names(net$MEs)]

  ## get colors of eigengene modules
  color1 <- wgcna.labels2colors(net$colors)
  me.colors <- color1[!duplicated(color1)]
  me.labels <- net$labels[!duplicated(color1)]
  names(me.colors) <- me.labels
  me.colors <- me.colors[names(net$MEs)]

  ## compute clustering based on TOM matrix
  message("[wgcna.compute] computing TOM matrix...")
  TOM <- NULL
  TOM <- WGCNA::TOMsimilarityFromExpr(
    datExpr,
    power = power,
    TOMType = tomtype,
    networkType = networktype,
    verbose = verbose
  )

  ## instead of the huge TOM matrix we save a smaller SVD.
  message("[wgcna.compute] reducing TOM. sv.dim = ", sv.dim)
  rownames(TOM) <- colnames(TOM) <- colnames(datExpr)
  sv <- irlba::irlba(TOM, nv = min(sv.dim,dim(X)))
  svTOM <- sv$v %*% diag(sqrt(sv$d))
  rownames(svTOM) <- colnames(datExpr)
  
  ## compute module eigenvectors (loading matrix)
  message("[wgcna.compute] computing module eigenvectors...")
  MVs <- list()
  for (clr in unique(net$colors)) {
    ii <- which(net$colors == clr)
    mX <- datExpr[, ii]
    mX <- scale(mX) ## NOTE: seems WGCNA is using full scaling
    sv1 <- irlba::irlba(mX, nv = 1, nu = 1)
    sv <- rep(0, ncol(datExpr))
    sv[ii] <- sv1$v[, 1] * sv1$d[1]
    MVs[[paste0(prefix, clr)]] <- sv
  }
  MVs <- as.matrix(do.call(cbind, MVs[names(net$MEs)]))
  rownames(MVs) <- colnames(datExpr)

  ## compute gene statistics
  message("[wgcna.compute] computing gene statistics...")
  stats <- wgcna.compute_geneStats(net, datExpr, datTraits, TOM=TOM) 
  
  ## module-traits matrix
  message("[wgcna.compute] computing module-traits matrix...")
  modTraits <- cor(net$MEs, datTraits, use="pairwise")
  
  ## merge dendrograms
  if(merge.dendro) {
    message("[wgcna.compute] merge_block_dendrograms...")  
    merged <- try(wgcna.merge_block_dendrograms(net, TOM))
    if(!inherits(merged,"try-error")) {
      net$merged_dendro <- merged
    } else {
      net$merged_dendro <- NULL
    }
  }
  remove(TOM)
  
  ## construct results object
  results <- list(
      datExpr = datExpr,
      datTraits = datTraits,
      ## TOM = TOM,  ## this can be BIG!!!
      svTOM = svTOM,  ## smaller singular vectors
      net = net,
      power = power,
      me.genes = me.genes,
      me.colors = me.colors,
      W = MVs,
      modTraits = modTraits,
      stats = stats
  )

  return(results)
}

wgcna.compute_geneStats <- function(net, datExpr, datTraits, TOM) {
  ## Define numbers of genes and samples
  nGenes <- ncol(datExpr)
  nSamples <- nrow(datExpr)

  ## Recalculate MEs with color labels
  moduleTraitCor <- cor(net$MEs, datTraits, use = "pairwise.complete")
  moduleTraitPvalue <- WGCNA::corPvalueStudent(moduleTraitCor, nSamples)

  ## Module membership correlation (with p-values)
  moduleMembership <- cor(datExpr, net$MEs, use = "p")
  MMPvalue <- WGCNA::corPvalueStudent(as.matrix(moduleMembership), nSamples)

  ## Gene-trait significance (trait correlation) (with p-values)
  traitSignificance <- cor(datExpr, datTraits, use = "p")
  GSPvalue <- WGCNA::corPvalueStudent(as.matrix(traitSignificance), nSamples)

  ## Fold-change
  is.binary <- apply(datTraits, 2, function(a) length(unique(a[!is.na(a)])) == 2)
  Y <- datTraits[, which(is.binary)]
  lm <- lapply(Y, function(y) {
    gx.limma(t(datExpr), y,
      lfc = 0, fdr = 1,
      sort.by = "none", verbose = 0
    )
  })
  foldChange <- sapply(lm, function(m) m$logFC)
  foldChangePvalue <- sapply(lm, function(m) m$P.Value)
  rownames(foldChange) <- rownames(lm[[1]])
  rownames(foldChangePvalue) <- rownames(lm[[1]])

  ## Gene Centrality. Compute centrality of gene in Module subgraph
  ## using TOM matrix.
  adj <- TOM
  diag(adj) <- NA
  adj[which(abs(adj) < 0.01)] <- 0
  gr <- igraph::graph_from_adjacency_matrix(
    adj,
    mode = "undirected", weighted = TRUE, diag = FALSE
  )
  geneCentrality <- rep(NA, nrow(adj))
  names(geneCentrality) <- rownames(adj)
  me.genes <- tapply(names(net$colors), net$colors, list)
  gg <- me.genes[[1]]
  for (gg in me.genes) {
    gr1 <- igraph::subgraph(gr, gg)
    ct <- igraph::page_rank(gr1, weights = NULL)$vector
    ct <- ct / mean(ct, na.rm = TRUE)
    geneCentrality[gg] <- ct
  }

  ## propotion variance explained (PVE) per module
  propVarExplained <- WGCNA::propVarExplained(
    datExpr,
    colors = net$colors, MEs = net$MEs
  )

  stats <- list(
    moduleTraitCor = moduleTraitCor,
    moduleTraitPvalue = moduleTraitPvalue,
    moduleMembership = moduleMembership,
    MMPvalue = MMPvalue,
    traitSignificance = traitSignificance,
    GSPvalue = GSPvalue,
    foldChange = foldChange,
    foldChangePvalue = foldChangePvalue,
    geneCentrality = geneCentrality
  )

  return(stats)
}


## ----------------------------------------------------
## Perform geneset analysis on modules
## ----------------------------------------------------
wgcna.compute_enrichment <- function(wgcna, pgx,
                                     method = c("fisher","gsetcor","mmcor","wcor"),
                                     ntop = 200,
                                     annot = NULL, GMT = NULL,
                                     filter = "^PATHWAY|^HALLMARK|^GO|^C[1-9]") {
  ## collapse features to symbol
  geneX <- t(wgcna$datExpr)
  symbol.col <- NULL
  ##annot=NULL;GMT=NULL
  if(!is.null(pgx)) {
    GMT <- pgx$GMT
    gsetX <- pgx$gsetX
    annot <- pgx$genes
    if(!is.null(annot)) {
      symbol.col <- intersect(c("symbol","gene_name"),colnames(annot))[1]
    }
  }

  if (is.null(GMT)) {
    stop("FATAL. must supply GMT matrix")
    return(NULL)
  }
  if (is.null(gsetX)) {
    stop("FATAL. must supply GMT matrix")
    return(NULL)
  }
  if(!is.null(annot)) {
    geneX <- rename_by(geneX, annot, symbol.col)  
  } else {
    rownames(geneX) <- gsub(".*:|[_ -].*","",rownames(geneX))
  }

  bg <- rownames(geneX)
  bg <- intersect(bg, rownames(GMT))
  if(length(bg)==0) {
    message("FATAL. no overlapping genes")
    return(NULL)
  }
  G1 <- GMT[bg,, drop=FALSE]
  if(!is.null(filter)) {
    sel <- grep(filter,colnames(G1))
    if(length(sel)) G1 <- G1[,sel, drop=FALSE]
  }
  G1 <- G1[, which(Matrix::colSums(G1 != 0) >= 4), drop = FALSE]
  gmt <- mat2gmt(G1)

  W <- as.matrix(wgcna$W)
  if (!is.null(annot)) {
    W <- rename_by(W, annot, symbol.col)
  }

  rho.list <- list()
  pval.list <- list()

  ## Perform fisher-test on ME genes
  if ("fisher" %in% method) {
    message("[wgcna.compute_enrichment] calculating Fisher tests...")
    i <- 1
    rho <- NULL
    pval <- NULL
    me.genes <- wgcna$me.genes
    for (i in 1:length(me.genes)) {
      gg <- me.genes[[i]]
      if (!is.null(annot)) {
        gg <- probe2symbol(me.genes[[i]], annot, query = symbol.col)
      }
      ## perform Fisher test
      rr <- try(gset.fisher(gg, gmt, background = bg, fdr = 1, min.genes = 4, verbose = 0))
      if (!"try-error" %in% class(rr)) {
        rr <- rr[match(names(gmt), rownames(rr)), ]
        if (is.null(rho)) {
          rho <- cbind(rr$odd.ratio)
          pval <- cbind(rr$p.value)
          colnames(rho) <- colnames(pval) <- names(me.genes)[i]
        } else {
          rho <- cbind(rho, rr$odd.ratio)
          pval <- cbind(pval, rr$p.value)
          me <- names(me.genes)[i]
          colnames(rho)[ncol(rho)] <- me
          colnames(pval)[ncol(pval)] <- me
        }
        rownames(rho) <- names(gmt)
        rownames(pval) <- names(gmt)
      }
    }
    ## handle infinite
    rho[is.infinite(rho)] <- 2 * max(rho, na.rm = TRUE)
    rho[is.na(rho)] <- 0
    pval[is.na(pval)] <- 1
    rho.list[["fisher"]] <- rho
    pval.list[["fisher"]] <- pval
  }

  ## Here we correlate geneset score (averageCLR) with the module
  ## eigengene (ME). This should select genesets correlated with the
  ## ME.
  if ("gsetcor" %in% method) {
    message("[wgcna.compute_enrichment] calculating geneset correlation...")
    gsetX <- gsetX[colnames(G1), ]
    ME <- as.matrix(wgcna$net$MEs)
    if (ncol(gsetX) > nrow(ME)) {
      kk <- sample(colnames(gsetX), nrow(ME))
      gsetX <- gsetX[, kk, drop = FALSE]
    }
    rc.rho <- cor(t(gsetX), ME, use = "pairwise")
    rc.pvalue <- cor.pvalue(rc.rho, n = nrow(ME))
    rho.list[["gsetcor"]] <- rc.rho
    pval.list[["gsetcor"]] <- rc.pvalue
  }

  ## Here we correlate genes with the module eigengene (ME) then do
  ## a gset.rankcor() on the ME correlation.
  if ("mmcor" %in% method) {
    message("[wgcna.compute_enrichment] calculating MM correlation...")
    mm <- cor(t(geneX), wgcna$net$MEs, use = "pairwise")
    rc <- gset.rankcor(mm, G1, compute.p = TRUE) ## NEEDS CHECK!!!
    rho.list[["mmcor"]] <- rc$rho
    pval.list[["mmcor"]] <- rc$p.value
  }

  ## Here we correlate genes with the module eigengene (ME) then do
  ## a gset.rankcor() on the ME correlation.
  if ("wcor" %in% method) {
    message("[wgcna.compute_enrichment] calculating W correlation...")
    rc <- gset.rankcor(W, G1, compute.p = TRUE) ## NEEDS CHECK!!!
    rho.list[["wcor"]] <- rc$rho
    pval.list[["wcor"]] <- rc$p.value
  }

  ## Compute meta rank and pval
  meta.p <- Reduce(pmax, pval.list) ## NEED RETHINK!!!
  meta.q <- apply(meta.p, 2, p.adjust, method = "fdr")
  rnk.list <- lapply(rho.list, function(x) apply(x, 2, rank) / nrow(x))
  meta.rnk <- Reduce("+", rnk.list) / length(rnk.list)

  ## create dataframe by module
  message("[wgcna.compute_enrichment] creating dataframes...")
  gse.list <- list()
  i <- 1
  for (i in 1:ncol(meta.p)) {
    k <- colnames(meta.p)[i]
    pv <- sapply(pval.list, function(x) x[, i])
    colnames(pv) <- paste0("p.", colnames(pv))
    df <- data.frame(
      module = k,
      geneset = rownames(meta.p),
      score = meta.rnk[, i],
      p.value = meta.p[, i],
      q.value = meta.q[, i],
      pv
    )
    df <- df[order(-df$score), ]
    if (!is.null(ntop) && ntop > 0) df <- head(df, n = ntop)
    gse.list[[k]] <- df
  }

  ## add genes
  gse.genes <- list()
  k <- names(gse.list)[1]
  for (k in names(gse.list)) {
    gset <- rownames(gse.list[[k]])
    gg <- wgcna$me.genes[[k]]
    if (!is.null(annot)) {
      gg <- probe2symbol(gg, annot, query = symbol.col)
    } else {
      gg <- gsub(".*:","",gg)
    }
    set.genes <- lapply(gmt[gset], function(s) intersect(s, gg))
    n0 <- sapply(gmt[gset], length)
    n1 <- sapply(set.genes, length)
    gse.genes[[k]] <- sort(table(unlist(set.genes)), decreasing = TRUE)
    set.genes <- sapply(set.genes, function(g) paste(sort(g), collapse = "|"))
    gse.list[[k]]$genes <- set.genes
    gse.list[[k]]$overlap <- paste0(n1, "/", n0)
  }

  ## compute gene enrichment
  gene.list <- list()
  ## message("[wgcna.compute_enrichment] computing gene enrichment...")
  ## k = names(gse.list)[1]
  ## for(k in names(gse.list)) {
  ##   me <- gse.list[[k]]
  ##   rnk <- me$score
  ##   names(rnk) <- rownames(me)
  ##   me.genes <- strsplit(me$genes,split='\\|')
  ##   num.genes <- sapply(me.genes,length)
  ##   if(any(num.genes>0)) {
  ##     ii <- which(num.genes>0)
  ##     me.names <- unlist(mapply(rep, rownames(me)[ii], num.genes[ii]))
  ##     me.genes <- me.genes[ii]
  ##     me.gmt <- tapply(me.names, unlist(me.genes), list)
  ##     res <- fgsea::fgsea(me.gmt, rnk)
  ##     res <- res[order(-res$NES),]
  ##     gene.list[[k]] <- head(res,100)
  ##   }
  ## }

  res <- list(
    gsets = gse.list,
    genes = gene.list
  )

  return(res)
}


#' @export
wgcna.runConsensusWGCNA <- function(exprList,
                                    datTraits,
                                    power = 12,
                                    minKME = 0.8,
                                    cutheight = 0.15,
                                    deepsplit = 2,
                                    maxBlockSize = 5000
                                    ) {

  ##exprList <- list(Set1 = V1, Set2 = V2)
  multiExpr = WGCNA::list2multiData(lapply(exprList,Matrix::t))

  # The variable exprSize contains useful information about the sizes of all
  # of the data sets now we run automatic module detection procedure
  net.list <- list()
  for (k in names(multiExpr)) {
    message("[wgcna.runConsensusWGCNA] computing WGCNA for ", k)
    net.list[[k]] <- WGCNA::blockwiseModules(
      datExpr = multiExpr[[k]]$data,
      power = power,
      networkType = "signed",
      TOMType = "signed",
      minModuleSize = 20,
      deepSplit = deepsplit,
      mergeCutHeight = cutheight, 
      numericLabels = FALSE,
      minKMEtoStay = minKME,
      saveTOMs = FALSE,
      maxBlockSize = maxBlockSize,      
      verbose = 0
    )
  }

  # now we run automatic consensus module detection
  message("[wgcna.runConsensusWGCNA] computing consensus modules...")
  net = WGCNA::blockwiseConsensusModules(
    multiExpr[],
    power = power, 
    networkType = "signed",
    TOMType = "signed",
    minModuleSize = 20,
    deepSplit = deepsplit,
    mergeCutHeight = cutheight, 
    numericLabels = FALSE,
    minKMEtoStay = minKME,
    maxBlockSize = maxBlockSize,
    saveTOMs = FALSE,
    verbose = 1
  )

  ## create and match colors
  colors <- sapply(net.list, function(net) net$colors)
  c0 <- net$colors
  matched.colors <- apply(colors, 2, function(k) WGCNA::matchLabels(k,c0))
  colors <- cbind(Consensus=c0, matched.colors)
  
  ## add labels to dendrogram
  for(i in 1:length(net$dendrograms)) {
    ii <- which(net$goodGenes & net$blocks==i)
    xnames <- names(net$colors)
    net$dendrograms[[i]]$labels <- xnames[ii]
  }
  
  ## merge dendrograms
  message("[wgcna.compute] merge_block_dendrograms...")  
  multiX <- Matrix::t(do.call(rbind,lapply(exprList,function(x)scale(t(x)))))
  merged <- try(wgcna.merge_block_dendrograms(net, multiX))
  if(!inherits(merged,"try-error")) {
    net$merged_dendro <- merged
  } else {
    net$merged_dendro <- NULL
  }
  
  ## create module-trait matrices for each set
  message("[wgcna.compute] computing module-traits matrices...")
  Z.list <- list()
  k=1
  for(k in names(net$multiME)) {
    M <- (net$multiME[[k]][[1]])
    Z.list[[k]] <- cor(M, datTraits[rownames(M),], use="pairwise")
  }
  z.check <- sapply(Z.list, function(z) colSums(z != 0, na.rm = TRUE) > 0)
  sel <- names(which(rowMeans(z.check) == 1))
  sel
  Z.list <- lapply(Z.list, function(z) z[, sel, drop = FALSE])
  Z.list

  ## create consensus module-trait matrix
  message("[wgcna.compute] computing consensus module-traits matrix...")
  mdim <- sapply(exprList, ncol)
  pv <- mapply(function(z, n) corPvalueStudent(z, n), Z.list, mdim, SIMPLIFY = FALSE)
  all.sig <- Reduce("*", lapply(pv, function(p) 1 * (p < 0.05)))
  all.sig
  all.pos <- Reduce("*", lapply(Z.list, function(z) sign(z) == 1))
  all.neg <- Reduce("*", lapply(Z.list, function(z) sign(z) == -1))
  disconcordant <- !(all.pos | all.neg)
  disconcordant
  zsign <- sign(Reduce("+", lapply(Z.list, sign)))
  consZ <- Reduce(pmin, lapply(Z.list, abs)) * zsign
  consZ[disconcordant] <- NA
  ## consZ[!all.sig] <- NA
  ydim <- sapply(exprList, ncol)

  res <- list(
    net = net,
    modTraits = consZ,
    dendro = net$merged_dendro,    
    colors = colors,
    zlist = Z.list,
    ydim = ydim
  )

  res
}


#'
#'
#' @export
wgcna.getGeneStats <- function(wgcna, module, trait, plot = TRUE,
                               showallmodules = TRUE, col = NULL,
                               main = NULL) {
  
  p1 <- c("moduleMembership", "MMPvalue")
  p2 <- c("traitSignificance", "GSPvalue", "foldChange", "foldChangePvalue")
  p3 <- c("geneCentrality")
  nrow <- ncol(wgcna$datExpr)
  df <- wgcna$stats[[p1[1]]][, 0] ## empty data.frame
  rownames(df) <- colnames(wgcna$datExpr)

  mm.cols <- colnames(wgcna$stats[[p1[1]]])
  if (!is.null(module) && module %in% mm.cols) {
    A1 <- sapply(wgcna$stats[p1], function(x) x[, module])
    df <- cbind(df, A1)
  }

  tt.cols <- colnames(wgcna$stats[[p2[1]]])  
  if (!is.null(trait) && trait %in% tt.cols) {
    A2 <- sapply(wgcna$stats[p2], function(x) x[, trait])
    df <- cbind(df, A2)
  }
  
  A3 <- wgcna$stats[[p3]]
  df <- cbind(df, centrality = A3)
  
  sel <- c("moduleMembership", "traitSignificance", "foldChange", "centrality")
  sel <- intersect(sel, colnames(df))
  df1 <- pmax(as.matrix(df[, sel]), 1e-8)
  score <- exp(rowMeans(log(df1)))
  
  labels <- wgcna$net$labels
  df <- data.frame(module = labels, score = score, df)
  df <- df[order(-df$score), ]
  if (!is.null(module) && !showallmodules) {
    sel <- which(df$module == module)
    df <- df[sel, , drop = FALSE]
  }
  
  if (plot) {
    cols <- c("moduleMembership", "traitSignificance", "foldChange", "centrality")
    cols <- intersect(cols, colnames(df))
    df1 <- df[, cols]
    col1 <- wgcna.labels2colors(wgcna$net$colors[rownames(df1)])
    if (!is.null(col)) col1 <- col
    pairs(df1, col = col1)
    if (is.null(main)) {
      main <- paste("Gene significance for module", module, "and trait", trait)
    }
    title(main, line = 3, cex.main = 1.15)
  }

  df
}

wgcna.merge_block_dendrograms <- function(net, X, method=1) {
  ## This function is fragile: it can give a C stack limit error. In
  ## case that happens you can increase the stack limit by running on
  ## the cmd line: >>> ulimit -s unlimited
  ## 
  hc = net$dendrograms
  
  ## merge block dendrogram
  mx <- list()
  for(b in 1:length(net$dendrograms)) {
    ii <- which(net$goodGenes & net$blocks==b)
    mx[[b]] <- colMeans(X[ii,])
    hc[[b]]$labels <- rownames(X)[ii]
  }

  if(length(hc)==1) {
    return(hc[[1]])
  }
  
  M <- do.call(rbind, mx)
  hclust_p = hclust(dist(M),method="average")
  dend_p = as.dendrogram(hclust_p)
  dend.list = lapply(hc, as.dendrogram)
  
  if(method==1) {
    merged = ComplexHeatmap::merge_dendrogram(dend_p, dend.list)
  } else {
    mrg <- hclust_p$merge    
    merged_branch <- list()
    k=1
    for(k in 1:nrow(mrg)) {
      i <- mrg[k,1]
      j <- mrg[k,2]
      if(i<0) d1 <- dend.list[[-i]]
      if(i>0) d1 <- merged_branch[[i]]      
      if(j<0) d2 <- dend.list[[-j]]
      if(j>0) d2 <- merged_branch[[j]]
      merged_branch[[k]] <- merge(d1,d2)  ## actual merge
    }
    merged <- merged_branch[[k]]
  }
  
  merged_hclust <- as.hclust(merged)
  merged_hclust
}

##=========================================================================
## PLOTTING FUNCTIONS
## =========================================================================

#'
#'
#' @export
wgcna.plotTOM <- function(wgcna, justdata = FALSE, block=NULL,
                          legend=TRUE, downsample=NULL) {

  datExpr <- wgcna$datExpr
  wTOM <- NULL  
  if(!is.null(wgcna$TOM)) {
    wTOM <- wgcna$TOM
  }
  ## if SV of TOM is stored, reconstruct TOM
  if(is.null(wTOM) && !is.null(wgcna$svTOM)) {
    wTOM  <- wgcna$svTOM %*% t(wgcna$svTOM)
  }
  if(is.null(wTOM)) {
    message("[wgcna.plotTOM] ERROR. no TOM matrix")
    return(NULL)
  }
  dissTOM <- 1 - wTOM
  rownames(dissTOM) <- colnames(dissTOM) <- colnames(datExpr)
  
  ## clustering wgcnaults
  moduleColors <- NULL
  if(is.null(block) && "merged_dendro" %in% names(wgcna$net)) {
    geneTree <- wgcna$net$merged_dendro
    gg <- geneTree$labels    
    if(length(gg)>0) {
      dissTOM <- dissTOM[gg,gg]
      moduleColors <- wgcna.labels2colors(wgcna$net$colors[gg])
    }
  }

  if(is.null(moduleColors)) {
    if(is.null(block)) block <- 1
    geneTree <- wgcna$net$dendrograms[[block]]
    ii <- which(wgcna$net$blocks == block & wgcna$net$goodGenes==TRUE)
    gg <- names(wgcna$net$color)[ii]
    dissTOM <- dissTOM[gg,gg]
    moduleColors <- wgcna.labels2colors(wgcna$net$colors[gg])
  }

  if (justdata) {
    return(dissTOM)
  }
 
  if(!is.null(downsample) && ncol(dissTOM)>downsample) {
    ii <- seq(1,ncol(dissTOM),length.out=downsample)
    dissTOM <- dissTOM[ii,ii]
    moduleColors <- moduleColors[ii]
    geneTree <- fastcluster::hclust(as.dist(dissTOM),
                                    method = "average")
  }
  
  if(legend) {
    par(oma = c(1,0,0,0), mar=c(0,0,0,0))
    plotly::layout(
      matrix(c(
        0, 0, 5, 0,
        0, 0, 2, 0,
        4, 1, 3, 6
      ), nr = 3, byrow = TRUE),
      widths = c(2.3, 0.5, 10, 3),
      heights = c(2.3, 0.5, 10)
    )
    WGCNA::TOMplot(
      dissTOM^7,
      geneTree,
      moduleColors,
      setLayout = FALSE,
      main = NULL
    )
    
    ## add color legend
    frame()
    legend(
      -0.1, 1,
      fill = wgcna$me.colors,
      legend = names(wgcna$me.colors),
      cex = 1.2, bty = "n", x.intersp = 0.5
    )
  } else {
    WGCNA::TOMplot(
      dissTOM^7,
      geneTree,
      moduleColors,
      setLayout = TRUE,
      main = NULL
    )
  }
  
}

#'
#'
#' @export
wgcna.plotDendroAndColors <- function(wgcna, main=NULL, block=NULL, extra.colors=NULL) {

  if(length(wgcna$net$dendrograms)>1) {
    message("warning: this wgcna has multiple blocks")
  }

  if(is.null(block) && "merged_dendro" %in% names(wgcna$net)) {
    geneTree <- wgcna$net$merged_dendro
  } else {
    if(is.null(block)) block <- 1
    geneTree <- wgcna$net$dendrograms[[block]]
  }
  colors <- wgcna.labels2colors(wgcna$net$colors)
  gg <- geneTree$labels
  if(is.null(gg) && !is.null(block)) {
    ii <- which(wgcna$net$blocks == block & wgcna$net$goodGenes==TRUE)
    gg <- names(wgcna$net$color)[ii]
  }
  if(is.null(gg) && is.null(block)) {
    ii <- which(wgcna$net$goodGenes==TRUE)
    gg <- names(wgcna$net$color)[ii]
  }
  colors <- colors[gg]
  groupLabels <- "Module colors"
  if(!is.null(extra.colors)) {
    jj <- match(gg, rownames(extra.colors))
    colors <- cbind(colors, extra.colors[jj,])
    groupLabels <- c("Module colors",colnames(extra.colors))
  }
  
  if(is.null(main)) main <- "Gene dendrogram and module colors"
  ## Plot the dendrogram and the module colors underneath
  WGCNA::plotDendroAndColors(
    dendro = geneTree,
    colors = colors,
    groupLabels = groupLabels,
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = FALSE, guideHang = 0.05,
    marAll = c(0.2, 5, 1, 0.2),
    main = main
  )
}


#' @export
wgcna.labels2colors <- function(colors, ...) {
  if (all(is.numeric(colors))) {
    colors <- WGCNA::labels2colors(colors, ...)
    return(colors)
  }
  stdColors <- c("grey", WGCNA::standardColors())
  if (all(colors %in% stdColors)) {
    return(colors)
  }
  icolors <- as.integer(factor(as.character(colors)))
  colors <- WGCNA::standardColors()[icolors]
  return(colors)
}


#' @export
wgcna.plotModuleTraitHeatmap <- function(wgcna, setpar = TRUE, cluster = FALSE,
                                         main = NULL, justdata = FALSE,
                                         pstar = TRUE) {
  ## Define numbers of genes and samples
  nGenes <- ncol(wgcna$datExpr)
  nSamples <- nrow(wgcna$datExpr)
  if ("stats" %in% names(wgcna)) {
    moduleTraitCor <- wgcna$stats$moduleTraitCor
    moduleTraitPvalue <- wgcna$stats$moduleTraitPvalue
  } else {
    moduleTraitCor <- cor(wgcna$net$MEs, wgcna$datTraits, use = "pairwise.complete")
    moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
  }

  if (cluster) {
    ii <- hclust(dist(moduleTraitCor))$order
    jj <- hclust(dist(t(moduleTraitCor)))$order
  } else {
    ii <- rev(order(rownames(moduleTraitCor)))
    jj <- order(colnames(moduleTraitCor))
  }
  moduleTraitCor <- moduleTraitCor[ii, jj]
  moduleTraitPvalue <- moduleTraitPvalue[ii, jj]

  if (justdata) {
    return(moduleTraitCor)
  }

  ## Will display correlations and their p-values
  if (pstar) {
    textPv <- cut(moduleTraitPvalue,
      breaks = c(-1, 0.001, 0.01, 0.05, 99),
      labels = c("★★★", "★★", "★", "")
    )
  } else {
    textPv <- paste0("(", signif(moduleTraitPvalue, 1), ")")
  }

  textMatrix <- paste0(signif(moduleTraitCor, 2), "\n", textPv)
  dim(textMatrix) <- dim(moduleTraitCor)
  if (setpar) par(mar = c(8, 8, 3, 3))
  if (is.null(main)) main <- "Module-trait relationships"
  ## Display the correlation values within a heatmap plot
  WGCNA::labeledHeatmap(
    Matrix = moduleTraitCor,
    xLabels = colnames(moduleTraitCor),
    yLabels = rownames(moduleTraitCor),
    ySymbols = rownames(moduleTraitCor),
    colorLabels = FALSE,
    colors = WGCNA::blueWhiteRed(50),
    textMatrix = textMatrix,
    setStdMargins = FALSE,
    cex.text = 0.7,
    zlim = c(-1, 1),
    main = main
  )
}

#' Plot membership correlation vs gene signficance (correlation with
#' trait) to discover biomarkers/driver genes.
#'
#' @export
wgcna.plotMMvsGS <- function(wgcna, module, trait, abs = TRUE, par = TRUE,
                             plotlib = "base") {
  ## module="ME3";trait="activated=act"
  moduleGenes <- wgcna$me.genes[[module]]
  nSamples <- nrow(wgcna$datExpr)

  ## Module membership correlation (with p-values)
  if ("stats" %in% names(wgcna)) {
    moduleMembership <- wgcna$stats$moduleMembership
    MMPvalue <- wgcna$stats$MMPvalue
  } else {
    moduleMembership <- as.data.frame(cor(wgcna$datExpr, wgcna$net$MEs, use = "p"))
    MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(moduleMembership), nSamples))
  }

  ## Gene-trait significance (trait correlation) (with p-values)
  if ("stats" %in% names(wgcna)) {
    traitSignificance <- wgcna$stats$traitSignificance
    GSPvalue <- wgcna$stats$GSPvalue
  } else {
    traitSignificance <- as.data.frame(cor(wgcna$datExpr, wgcna$datTraits, use = "p"))
    GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(traitSignificance), nSamples))
  }

  x <- (moduleMembership[moduleGenes, module])
  y <- (traitSignificance[moduleGenes, trait])
  if (abs == TRUE) {
    x <- abs(x)
    y <- abs(y)
  }
  ##
  px <- MMPvalue[moduleGenes, module]
  py <- GSPvalue[moduleGenes, trait]
  qx <- p.adjust(px, method = "fdr")
  qy <- p.adjust(py, method = "fdr")
  is.sig <- (qx < 0.05 & qy < 0.05)
  sigx <- (qx < 0.05)
  sigy <- (qy < 0.05)
  ii <- which(is.sig)
  qv <- quantile(x[ii], prob = 0.1)[1]
  qh <- quantile(y[ii], prob = 0.1)[1]

  pos <- cbind(x, y)
  rownames(pos) <- moduleGenes
  is.sig1 <- c("notsig", "onesig", "sig")[1 + 1 * sigx + 1 * sigy]
  hi1 <- NULL
  ## hi1 <- head(rownames(pos),10)
  col1 <- c("grey70", "grey20", "red2")

  if (par) par(mfrow = c(1, 1), mar = c(5, 5, 3, 2))
  if (plotlib == "ggplot") {
    pgx.scatterPlotXY.GGPLOT(
      pos,
      var = is.sig1, hilight = hi1, col = col1,
      xlab = paste("Module membership in", module, "module"),
      ylab = paste("Gene significance for trait", trait),
      title = paste("Module membership vs. gene significance\n"),
      cex.title = 0.9,
      girafe = FALSE
    )
  } else if (plotlib == "girafe") {
    pgx.scatterPlotXY.GGPLOT(
      pos,
      var = is.sig1, hilight = hi1, col = col1,
      xlab = paste("Module membership in", module, "module"),
      ylab = paste("Gene significance for trait", trait),
      title = paste("Module membership vs. gene significance\n"),
      cex.title = 0.7, cex.axis = 0.7,
      girafe = TRUE
    )
  } else {
    ii <- which(is.sig1 == "notsig")
    verboseScatterplot(
      x[-ii], y[-ii],
      xlab = paste("Module membership in", module, "module"),
      ylab = paste("Gene significance for trait", trait),
      main = paste("Module membership vs. gene significance\n"),
      cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = col1[1]
    )
    ii <- which(is.sig1 == "onesig")
    points(x[ii], y[ii], col = col1[2])
    ii <- which(is.sig1 == "sig")
    points(x[ii], y[ii], col = col1[3])
    abline(v = qv, h = qh, col = "darkred")
  }
}


#' Plot cluster dendrogram with eigengenes and traits.
#'
#' @export
wgcna.plotEigenGeneClusterDendrogram <- function(wgcna, add_traits = TRUE,
                                                 main = NULL) {
  # Matrix with eigengenes and traits
  MET <- wgcna$net$MEs
  if (add_traits) {
    MET <- cbind(MET, wgcna$datTraits)
  }
  MET <- WGCNA::orderMEs(MET)
  if (NCOL(MET) <= 2) MET <- cbind(MET, MET) ## error if ncol(MET)<=2 !!!!
  if (is.null(main)) main <- "Eigengene cluster dendrogram"
  WGCNA::plotEigengeneNetworks(
    MET, main,
    marDendro = c(0, 4, 2, 0),
    plotHeatmaps = FALSE
  )
}

#' Plot the adjacency correlation heatmap matrix of eigengenes with or
#' without traits. This can show how traits cluster together with the
#' eigengenes.
#'
#' @export
wgcna.plotEigenGeneAdjacencyHeatmap <- function(wgcna, add_traits = TRUE,
                                                marx = 1, main = NULL,
                                                justdata = FALSE) {
  # Matrix with eigengenes and traits
  MET <- wgcna$net$MEs
  if (add_traits) {
    MET <- cbind(MET, wgcna$datTraits)
  }
  colnames(MET) <- paste0(" ", colnames(MET))
  MET <- WGCNA::orderMEs(MET)
  if (NCOL(MET) <= 2) MET <- cbind(MET, MET) ## error if ncol(MET)<=2 !!!!

  if (justdata) {
    R <- (1 + cor(MET)) / 2
    return(R)
  }

  # Plot the correlation heatmap matrix (note: this plot will overwrite
  # the dendrogram plot)
  if (is.null(main)) main <- "Eigengene adjacency heatmap"
  WGCNA::plotEigengeneNetworks(
    MET, main,
    marHeatmap = c(8 * marx, 10 * marx, 2, 2),
    plotDendrograms = FALSE,
    colorLabels = TRUE,
    xLabelsAngle = 45
  )
}

#' @export
wgcna.plotEigenGeneGraph <- function(wgcna, add_traits = TRUE, main = NULL) {
  ## require(igraph)
  net <- wgcna$net
  MET <- net$MEs
  if (add_traits) {
    MET <- cbind(MET, wgcna$datTraits)
  }
  if (NCOL(MET) <= 2) MET <- cbind(MET, MET) ## error if ncol(MET)<=2 !!!!

  ## Recalculate MEs with color as labels
  clust <- hclust(dist(t(scale(MET))))
  clust
  phylo <- ape::as.phylo(clust)
  gr <- igraph::as.igraph(phylo, directed = FALSE)

  is.node <- grepl("Node", igraph::V(gr)$name)
  module.name <- igraph::V(gr)$name
  module.size <- table(wgcna$net$labels)
  module.size <- module.size / mean(module.size)
  module.size

  igraph::V(gr)$label <- igraph::V(gr)$name
  igraph::V(gr)$label[is.node] <- NA
  igraph::V(gr)$color <- wgcna$me.colors[module.name]
  igraph::V(gr)$size <- 20 * (module.size[module.name])**0.4
  igraph::V(gr)$size[is.na(igraph::V(gr)$size)] <- 0

  ## par(mfrow = c(1, 1), mar = c(1, 1, 1, 1) * 0)
  igraph::plot.igraph(
    gr,
    layout = igraph::layout.kamada.kawai,
    vertex.label.cex = 0.8,
    edge.width = 3
  )
  if (!is.null(main)) title(main, line = -1.5)
}


#' Plot Multi-dimensional scaling (MDS) of centered data matrix.
#'
#' @export
wgcna.plotMDS <- function(wgcna, main = NULL, scale = FALSE) {
  cc <- wgcna.labels2colors(wgcna$net$color)
  pc <- svd(t(scale(wgcna$datExpr, scale = scale)), nv = 2)$u[, 1:2]
  # pc <- svd(t(scale(wgcna$datExpr)),nv=1)$u[,1:2]
  colnames(pc) <- c("MDS-x", "MDS-y")
  if (is.null(main)) main <- "MDS of features"
  plot(pc, col = cc, main = main)
}

#'
#'
#' @export
wgcna.plotFeatureUMAP <- function(wgcna, nhub = 3, method = "clust",
                                  scale = FALSE, main = NULL,
                                  plotlib = "base") {
  if (method == "clust" && "clust" %in% names(wgcna)) {
    pos <- wgcna$clust[["umap2d"]]
  } else if (method == "umap") {
    cX <- t(scale(wgcna$datExpr, scale = scale)) ## WGCNA uses correlation
    pos <- uwot::umap(cX)
    colnames(pos) <- c("UMAP-x", "UMAP-y")
    rownames(pos) <- colnames(wgcna$datExpr)
    ##  } else if(method=="mds") {
  } else {
    pos <- svd(t(scale(wgcna$datExpr, scale = scale)), nv = 2)$u[, 1:2]
    colnames(pos) <- c("MDS-x", "MDS-y")
    rownames(pos) <- colnames(wgcna$datExpr)
  }

  if (is.null(main)) main <- "Feature UMAP colored by module"

  ## get top hub genes
  mm <- wgcna$stats$moduleMembership
  hubgenes <- apply(mm, 2, function(x) head(names(sort(-x)), 3), simplify = FALSE)
  hubgenes
  sel <- which(names(hubgenes) != "MEgrey")
  hubgenes <- unlist(hubgenes[sel])
  col1 <- wgcna$net$colors
  genes1 <- names(which(col1 != "grey"))
  pgx.scatterPlotXY(
    pos,
    var = col1, col = sort(unique(col1)),
    hilight = genes1, hilight2 = hubgenes,
    cex.lab = 1.2, label.clusters = FALSE,
    title = main,
    plotlib = plotlib
  )
}

#' Plot module significance.
#'
#' @export
wgcna.plotModuleSignificance <- function(wgcna, trait, main = NULL, abs = FALSE) {
  ## cc <- paste0("ME",wgcna$net$color)
  cc <- wgcna.labels2colors(wgcna$net$color)
  if ("stats" %in% names(wgcna)) {
    traitSignificance <- wgcna$stats$traitSignificance
  } else {
    traitSignificance <- as.data.frame(cor(wgcna$datExpr, wgcna$datTraits, use = "p"))
    names(traitSignificance) <- names(wgcna$datTraits)
    rownames(traitSignificance) <- colnames(wgcna$datExpr)
  }
  geneSig <- traitSignificance[, trait]
  if (is.null(main)) main <- paste("Module significance with", trait)
  if (abs) geneSig <- abs(geneSig)
  WGCNA::plotModuleSignificance(
    geneSig,
    colors = cc, main = main, boxplot = FALSE
  )
}

#'
#'
#' @export
wgcna.plotSampleDendroAndColors <- function(wgcna,
                                            what = c("me", "traits", "both")[3],
                                            clust.expr = TRUE,
                                            main = NULL, justdata = FALSE) {
  MET0 <- wgcna$net$MEs
  MET <- MET0[, 0]
  if (any(what %in% c("me", "both"))) {
    MET <- cbind(MET, MET0)
  }
  if (any(what %in% c("traits", "both"))) {
    MET <- cbind(MET, wgcna$datTraits)
  }
  if (NCOL(MET) <= 2) MET <- cbind(MET, MET) ## error if ncol(MET)<=2 !!!!

  ## Recalculate MEs with color as labels
  if (clust.expr) {
    sampleTree <- hclust(dist(wgcna$datExpr), method = "average")
  } else {
    sampleTree <- hclust(dist(scale(MET0)), method = "average")
  }
  ii <- sampleTree$order
  jj <- hclust(dist(t(scale(MET))))$order
  colors <- WGCNA::numbers2colors(MET[, jj])

  if (justdata) {
    return(MET)
  }

  if (is.null(main)) {
    if (what == "me") main <- "Sample dendrogram and module heatmap"
    if (what == "traits") main <- "Sample dendrogram and trait heatmap"
    if (what == "both") main <- "Sample dendrogram and module+traits heatmap"
  }

  ## Plot the dendrogram and the module colors underneath
  plotDendroAndColors(
    dendro = sampleTree,
    colors = colors,
    groupLabels = colnames(MET)[jj],
    dendroLabels = rownames(MET),
    hang = 0.03,
    addGuide = FALSE,
    guideHang = 0.05,
    marAll = c(0.2, 7, 1.5, 0.5),
    main = main
  )
}


#' @export
wgcna.plotLabeledCorrelationHeatmap <- function(R, nSamples, setpar = TRUE,
                                                cluster = TRUE,
                                                main = NULL, justdata = FALSE,
                                                pstar = TRUE) {
  ## Define numbers of genes and samples
  Pvalue <- corPvalueStudent(R, nSamples)
  if (cluster) {
    R0 <- R
    R0[is.na(R0)] <- 0
    ii <- hclust(dist(R0))$order
    jj <- hclust(dist(t(R0)))$order
  } else {
    ii <- rev(order(rownames(R)))
    jj <- order(colnames(R))
  }
  R <- R[ii, jj]
  Pvalue <- Pvalue[ii, jj]

  if (justdata) {
    return(R)
  }

  ## Will display correlations and their p-values
  if (pstar) {
    textPv <- cut(Pvalue,
      breaks = c(-1, 0.001, 0.01, 0.05, 99),
      labels = c("★★★", "★★", "★", "")
    )
  } else {
    textPv <- paste0("(", signif(Pvalue, 1), ")")
  }

  textMatrix <- paste0(signif(R, 2), "\n", textPv)
  textMatrix[which(is.na(R))] <- ""
  dim(textMatrix) <- dim(R)
  if (setpar) par(mar = c(8, 8, 3, 3))
  if (is.null(main)) main <- "Correlation heatmap"
  ## Display the correlation values within a heatmap plot
  WGCNA::labeledHeatmap(
    Matrix = R,
    xLabels = colnames(R),
    yLabels = rownames(R),
    xSymbols = colnames(R),
    ySymbols = rownames(R),
    colorLabels = FALSE,
    colors = WGCNA::blueWhiteRed(50),
    textMatrix = textMatrix,
    setStdMargins = FALSE,
    cex.text = 0.7,
    zlim = c(-1, 1),
    main = main
  )
}


#' Plot module hub genes as graph in circle layout. This shows the
#' main connectivity structure for the modules. It can be used for
#' (edge) preservation analysis for comparing different group of
#' samples..
#'
#' @export
wgcna.plotModuleHubGenes <- function(wgcna, modules = NULL,
                                     alpha = 0.5, setpar = TRUE) {
  if (is.null(modules)) {
    modules <- colnames(wgcna$stats$moduleMembership)
  }
  modules <- intersect(modules, colnames(wgcna$stats$moduleMembership))
  if (length(modules) == 0) {
    message("ERROR. no valid modules")
    return(NULL)
  }

  if (setpar) {
    nr <- floor(length(modules)**0.5)
    nc <- ceiling(length(modules) / nr)
    nr
    nc
    par(mfrow = c(nr, nc), mar = c(0, 1, 2.5, 1))
  }
  for (k in modules) {
    mm <- wgcna$stats$moduleMembership[, k]
    mm.score <- head(sort(mm, decreasing = TRUE), 30)
    topgenes <- names(mm.score)
    A <- cor(wgcna$datExpr[, topgenes])
    diag(A) <- NA
    A <- (A - min(A, na.rm = TRUE)) / (max(A, na.rm = TRUE) - min(A, na.rm = TRUE))
    gr <- igraph::graph_from_adjacency_matrix(
      A,
      mode = "undirected", weighted = TRUE, diag = FALSE
    )
    norm.mm.score <- (mm.score - min(mm.score)) / (max(mm.score) - min(mm.score))
    clr <- sub("ME", "", k)
    if (!is.na(as.integer(clr))) clr <- as.integer(clr)
    if (clr == "black") clr <- "grey40"
    plot(gr,
      layout = igraph::layout_in_circle,
      edge.width = 6 * igraph::E(gr)$weight**8,
      vertex.size = 5 + 15 * norm.mm.score,
      vertex.color = adjustcolor(clr, alpha.f = alpha),
      vertex.frame.color = clr
    )
    title(paste("Module", k), line = 0.33)
  }
}

#' @title Map module colors to rainbow palette
#'
#' @description
#' Maps the module colors from a WGCNA network to a rainbow palette.
#'
#' @param net A WGCNA network object.
#'
#' @details
#' This function takes a WGCNA network object and maps the module colors
#' to a rainbow palette based on the modules' hierarchical clustering order.
#'
#' It extracts the hierarchical clustering dendrogram order, gets the number
#' of unique module colors, ranks the module color means based on the
#' dendrogram order, and assigns rainbow colors accordingly.
#'
#' This allows easier visualization and interpretation of modules in the
#' standard rainbow palette order.
#'
#' @return
#' A named character vector mapping the original module colors to rainbow colors.
#'
#' @export
labels2rainbow <- function(net) {
  hc <- net$dendrograms[[1]]
  nc <- length(unique(net$colors))
  n <- length(net$colors)
  ii <- rep(NA, n)
  ii[net$goodGenes] <- hc$order
  col1 <- WGCNA::labels2colors(net$colors)
  col.rnk <- rank(tapply(1:n, col1[ii], mean))
  new.col <- grDevices::rainbow(nc)[col.rnk]
  names(new.col) <- names(col.rnk)
  new.col["grey"] <- "#AAAAAA"
  new.col <- new.col[col1]
  names(new.col) <- net$colors
  return(new.col)
}


#' Filter color vector by minimum KME and mergeCutHeight. Set color of
#' features with KME smaller than minKME to grey (or 0) group. Merge
#' similar modules with (module) correlation larger than
#' (1-mergeCutHeight) together.
#'
#' @export
wgcna.filterColors <- function(X, colors, minKME=0.3, mergeCutHeight=0.15,
                               minmodsize = 20, ntop=-1 ) {
  ##minKME=0.3;mergeCutHeight=0.15;minmodsize=20;ntop=-1

  sX <- X + 1e-8*matrix(rnorm(length(X)),nrow(X),ncol(X))
  sX <- t(scale(t(sX)))

  ## get singular vectors and correct sign
  vv <- tapply(1:nrow(sX), colors, function(i) svd(sX[i,],nv=1)$v[,1])
  mm <- tapply(1:nrow(sX), colors, function(i) colMeans(sX[i,]))
  vv.sign <- mapply( function(a,b) sign(cor(a,b)), mm, vv)
  vv <- mapply( function(a,b) a*b, vv, vv.sign, SIMPLIFY=FALSE)

  kme <- rep(NA,nrow(X))
  names(kme) <- rownames(X)
  names(colors) <- rownames(X)

  grey.val <- NULL
  is.color <- mean(colors %in% WGCNA::standardColors(435)) > 0.8  
  if(is.numeric(colors)) {
    colors <- as.integer(colors)
    grey.val <- 0
  } else {
    colors <- as.character(colors)
    grey.val <- "---"
    if(is.color) grey.val <- "grey"
  }
  names(colors) <- rownames(X)  
  new.colors <- colors

  if(minKME > 0) {
    i=1
    for(i in 1:length(vv)) {
      ii <- which(colors == names(vv)[i])
      r <- cor(t(X[ii,]), vv[[i]])[,1]
      max(r)
      jj <- ii[which(r < minKME)]
      if(length(jj)) {
        new.colors[jj] <- NA
      }
      kme[ii] <- r
    }  
    new.colors[is.na(new.colors)] <- grey.val
  }

  ## merge groups
  if(mergeCutHeight > 0) {
    mx <- rowmean(X, new.colors)
    rr <- cor(t(mx))
    diag(rr) <- 0
    merge.idx <- which(rr > (1 - mergeCutHeight), arr.ind=TRUE)
    if(nrow(merge.idx)>0) {
      i=1
      for(i in 1:nrow(merge.idx)) {
        aa <- rownames(rr)[merge.idx[i,]]
        jj <- which(new.colors  %in% aa)
        max.color <- names(which.max(table(new.colors[jj])))
        new.colors[jj] <- max.color
      }
    }
  }

  ## remove small groups
  modsize <- table(new.colors)
  modsize
  if( min(modsize) < minmodsize ) {
    small.mod <- names(which(modsize < minmodsize))
    sel <- which( new.colors %in% small.mod)
    new.colors[sel] <- NA
  }

  if(ntop>0) {
    keep <- tapply( names(kme), new.colors, function(i) head(names(sort(-kme[i])),ntop) )
    keep <- unlist(keep)
    not.keep <- setdiff(names(kme), keep)
    if(length(not.keep)) new.colors[not.keep] <- NA
  }

  new.colors[which(is.na(new.colors))] <- grey.val
  ##if(!is.numeric(colors)) new.colors <- factor(new.colors)

  return(new.colors)
}

