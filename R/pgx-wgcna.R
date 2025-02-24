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
    ngenes = 2000,
    gset.filter = "PATHWAY|HALLMARK|^GO|^C[1-9]"
    ) {

  ##minmodsize=10;power=NULL;cutheight=0.25;deepsplit=2;ngenes=2000;networktype="signed";tomtype="signed";numericlabels=FALSE;ngenes=2000;gset.filter=NULL

  samples <- pgx$samples
  ## no dot pheno
  samples <- samples[, grep("^[.]",colnames(samples),invert=TRUE),drop=FALSE]  
  
  wgcna <- wgcna.compute(
    X = pgx$X,
    samples = samples,
    minmodsize = minmodsize,   # default: min(20,...)
    power = power,             # default: 12 (for signed)
    cutheight = cutheight,     # default: 0.15
    deepsplit = deepsplit,     # default: 2
    minKME = minKME,           # default: 0.30
    networktype = networktype, # default: unsigned (but signed is better...)
    tomtype = tomtype,         # default: signed                  
    numericlabels = numericlabels,
    ngenes = ngenes)

  ##---------------------------------------------------
  ## compute clustering based on TOM matrix
  ##---------------------------------------------------
  dissTOM <- 1 - wgcna$TOM
  rownames(dissTOM) <- colnames(dissTOM) <- colnames(wgcna$datExpr)
  clust <- pgx.clusterBigMatrix(dissTOM, methods = c("umap", "tsne", "pca"), dims = c(2))
  if ("cluster.genes" %in% names(pgx)) {
    posx <- pgx$cluster.genes$pos[["umap2d"]]
    clust[["umap2d"]] <- posx[colnames(wgcna$datExpr), ]
  }

  ##----------------------------------------------------
  ## Do geneset analysis 
  ##----------------------------------------------------
  res.gse <- wgcna.compute_enrichment(
    wgcna, pgx, 
    method = c("fisher","gsetcor","mmcor"),
    ntop = 1000,
    filter = gset.filter) 
  
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
                          prefix = "ME"
                          ) {

  #minmodsize=20;power=12;cutheight=0.15;deepsplit=2;ngenes=2000;networktype="signed";tomtype="signed";numericlabels=FALSE;prefix="ME";minKME=0.3;reassignThreshold=1e-6

  nmissing <- sum(is.na(X))
  if (nmissing > 0) {
    message("Found ", nmissing, " missing values in X. Imputing prior to WGCNA.")
    ##X <- X[complete.cases(X), , drop = FALSE]
    X <- svdImpute2(X)
  }
  X <- X[!duplicated(rownames(X)), ]

  if( ngenes > 0 && nrow(X) > ngenes ) {
    sdx <- matrixStats::rowSds(X, na.rm = TRUE)
    X <- X[sdx > 0.1 * mean(sdx, na.rm = TRUE), ] ## filter low SD??
    X <- X[order(-matrixStats::rowSds(X, na.rm = TRUE)), ]
    X <- utils::head(X, ngenes)
  }

  message("[wgcna.compute] dim(X) = ", paste(dim(X),collapse=' x '))
  message("[wgcna.compute] dim(samples) = ", paste(dim(samples),collapse=' x '))
  
  datExpr <- t(X)

  if(is.null(power)) {
    ## Estimate best power
    message("[wgcna.compute] estimating optimal power...")        
    powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
    sft <- WGCNA::pickSoftThreshold(
      datExpr,
      powerVector = powers,
      networkType = networktype,
      verbose = 0
    )
    power <- sft$powerEstimate
    if(is.na(power)) power <- 20
  }
  
  ## adapt for small datasets (also done in WGCNA package)
  minmodsize = min( minmodsize, ncol(datExpr)/2 )
  
  message("[wgcna.compute] computing blockwise Modules...")        
  message("[wgcna.compute] minmodsize = ", minmodsize)
  message("[wgcna.compute] number of features = ", nrow(X))
  message("[wgcna.compute] minKME = ", minKME)      
  message("[wgcna.compute] power = ", power)
  message("[wgcna.compute] mergeCutHeight = ", cutheight)

  WGCNA::enableWGCNAThreads()
  cor <- WGCNA::cor  ## needed...
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
    verbose = 0
  )
  cor <- stats::cor

  ##  prefix="ME"
  table(net$colors)
  names(net$MEs) <- sub("^ME",prefix,names(net$MEs))
  net$labels <- paste0(prefix,net$colors)

  ## clean up traits matrix
  datTraits <- data.frame(samples, check.names=FALSE)
  ## no dates please...
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
  TOM <- NULL
  TOM <- WGCNA::TOMsimilarityFromExpr(
    datExpr,
    power = power,
    TOMType = tomtype,
    networkType = networktype,
    verbose = 0
  )
  rownames(TOM) <- colnames(TOM) <- colnames(datExpr)

  ## compute module eigenvectors (loading matrix)
  MVs <- list()
  for(clr in unique(net$colors)) {
    ii <- which(net$colors == clr)
    mX <- datExpr[,ii]
    mX <- scale(mX)  ## NOTE: seems WGCNA is using full scaling
    sv1 <- irlba::irlba(mX, nv=1, nu=1)
    sv <- rep(0, ncol(datExpr))
    sv[ii] <- sv1$v[,1] * sv1$d[1]
    MVs[[paste0(prefix,clr)]] <- sv
  }
  MVs <- as.matrix(do.call(cbind, MVs[names(net$MEs)]))
  rownames(MVs) <- colnames(datExpr)

  ## compute gene statistics
  stats <- wgcna.compute_geneStats(net, datExpr, datTraits, TOM=TOM) 
  
  ## construct results object
  results <- list(
      datExpr = datExpr,
      datTraits = datTraits,
      TOM = TOM,
      net = net,
      power = power,
      me.genes = me.genes,
      me.colors = me.colors,
      W = MVs,
      stats = stats
  )

  return(results)
}

wgcna.compute_geneStats <- function(net, datExpr, datTraits, TOM) {
  
  ## Define numbers of genes and samples
  nGenes = ncol(datExpr);
  nSamples = nrow(datExpr);  
  
  ## Recalculate MEs with color labels
  moduleTraitCor = cor(net$MEs, datTraits, use = "pairwise.complete");
  moduleTraitPvalue = WGCNA::corPvalueStudent(moduleTraitCor, nSamples);

  ## Module membership correlation (with p-values)
  moduleMembership = cor(datExpr, net$MEs, use = "p");
  MMPvalue = WGCNA::corPvalueStudent(as.matrix(moduleMembership), nSamples);
  
  ## Gene-trait significance (trait correlation) (with p-values)
  traitSignificance = cor(datExpr, datTraits, use = "p");
  GSPvalue = WGCNA::corPvalueStudent(as.matrix(traitSignificance), nSamples);

  ## Fold-change
  is.binary <- apply(datTraits,2,function(a) length(unique(a[!is.na(a)]))==2)
  Y <- datTraits[,which(is.binary)]
  lm <- lapply(Y, function(y) gx.limma( t(datExpr), y, lfc=0, fdr=1,
                                               sort.by='none', verbose=0))
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
    adj, mode="undirected", weighted=TRUE, diag=FALSE )
  geneCentrality <- rep(NA, nrow(adj))
  names(geneCentrality) <- rownames(adj)
  me.genes <- tapply( names(net$colors), net$colors, list)
  gg <- me.genes[[1]]
  for(gg in me.genes) {
    gr1 <- igraph::subgraph(gr, gg)
    ct <- igraph::page_rank(gr1, weights=NULL)$vector
    ct <- ct / mean(ct,na.rm=TRUE)
    geneCentrality[gg] <- ct
  }

  ## propotion variance explained (PVE) per module
  propVarExplained <- WGCNA::propVarExplained(
    datExpr, colors=net$colors, MEs=net$MEs)
  
  ##
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


##----------------------------------------------------
## Perform geneset analysis on modules
##----------------------------------------------------
wgcna.compute_enrichment <- function(wgcna, pgx,
                                     method = c("fisher","gsetcor","mmcor","wcor"),
                                     ntop = 1000,
                                     annot = NULL, GMT = NULL,
                                     filter = "^PATHWAY|^HALLMARK|^GO|^C[1-9]") {
  
  ## collapse features to symbol
  geneX <- t(wgcna$datExpr)
  symbol.col <- NULL
  ##annot=NULL;GMT=NULL
  if(!is.null(pgx)) {
    annot <- pgx$genes
    GMT <- pgx$GMT
    symbol.col <- intersect(c("symbol","gene_name"),colnames(annot))[1]
    gsetX <- pgx$gsetX
  }

  if(is.null(GMT)) {
    stop("FATAL. must supply GMT matrix")
    return(NULL)
  }
  if(is.null(gsetX)) {
    stop("FATAL. must supply GMT matrix")
    return(NULL)
  }
  if(!is.null(annot)) {
    geneX <- rename_by(geneX, annot, symbol.col)  
  }
  
  bg <- rownames(geneX)
  bg <- intersect(bg, rownames(GMT))
  G1 <- GMT[bg,, drop=FALSE]
  if(!is.null(filter)) {
    sel <- grep(filter,colnames(G1))
    G1 <- G1[,sel, drop=FALSE]
  }
  G1 <- G1[, which(Matrix::colSums(G1!=0) >= 4), drop=FALSE]
  gmt <- mat2gmt(G1)
  
  W <- as.matrix(wgcna$W)
  if(!is.null(annot)) {
    W <- rename_by(W, annot, symbol.col)
  }

  rho.list  <- list()
  pval.list  <- list()  

  ## Perform fisher-test on ME genes
  if("fisher" %in% method) {
    message("[wgcna.compute_enrichment] calculating Fisher tests...")
    i=1
    rho <- NULL
    pval <- NULL
    me.genes <- wgcna$me.genes
    for (i in 1:length(me.genes)) {
      gg <- me.genes[[i]]
      if(!is.null(annot)) {
        gg <- probe2symbol(me.genes[[i]], annot, query = symbol.col)
      }
      ## perform Fisher test
      rr <- try(gset.fisher(gg, gmt, background=bg, fdr=1, min.genes=4, verbose=0))
      if (!"try-error" %in% class(rr)) {
        rr <- rr[match(names(gmt),rownames(rr)),]
        if (is.null(rho)) {
          rho  <- cbind(rr$odd.ratio)
          pval <- cbind(rr$p.value)
          colnames(rho) <- colnames(pval) <- names(me.genes)[i]
        } else {          
          rho  <- cbind(rho, rr$odd.ratio)
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
    rho[is.infinite(rho)] <- 2*max(rho,na.rm=TRUE)
    rho[is.na(rho)] <- 0    
    pval[is.na(pval)] <- 1
    rho.list[['fisher']] <- rho
    pval.list[['fisher']] <- pval
  }

  ## Here we correlate geneset score (averageCLR) with the module
  ## eigengene (ME). This should select genesets correlated with the
  ## ME.
  if("gsetcor" %in% method) {
    message("[wgcna.compute_enrichment] calculating geneset correlation...")
    gsetX <- gsetX[colnames(G1),]
    ME <- as.matrix(wgcna$net$MEs)
    rc.rho <- cor(t(gsetX), ME)
    rc.pvalue <- cor.pvalue( rc.rho, n=nrow(ME))
    rho.list[['gsetcor']] <- rc.rho
    pval.list[['gsetcor']] <- rc.pvalue
  }
  
  ## Here we correlate genes with the module eigengene (ME) then do 
  ## a gset.rankcor() on the ME correlation.
  if("mmcor" %in% method) {
    message("[wgcna.compute_enrichment] calculating MM correlation...")    
    mm <- cor(t(geneX), wgcna$net$MEs)
    rc <- gset.rankcor(mm, G1, compute.p=TRUE)  ## NEEDS CHECK!!!
    rho.list[['mmcor']] <- rc$rho
    pval.list[['mmcor']] <- rc$p.value
  }

  ## Here we correlate genes with the module eigengene (ME) then do 
  ## a gset.rankcor() on the ME correlation.
  if("wcor" %in% method) {
    message("[wgcna.compute_enrichment] calculating W correlation...")    
    rc <- gset.rankcor(W, G1, compute.p=TRUE)  ## NEEDS CHECK!!!
    rho.list[['wcor']] <- rc$rho
    pval.list[['wcor']] <- rc$p.value
  }

  ## Compute meta rank and pval
  meta.p <- Reduce( pmax, pval.list )  ## NEED RETHINK!!!
  meta.q <- apply(meta.p, 2, p.adjust, method="fdr")
  rnk.list <- lapply( rho.list, function(x) apply(x,2,rank)/nrow(x))
  meta.rnk <- Reduce( '+', rnk.list ) / length(rnk.list) 

  ## create dataframe by module  
  message("[wgcna.compute_enrichment] creating dataframes...")
  gse.list <- list()
  i=1
  for(i in 1:ncol(meta.p)) {
    k <- colnames(meta.p)[i]
    pv <- sapply( pval.list, function(x) x[,i])
    colnames(pv) <- paste0("p.",colnames(pv))
    df <- data.frame(
      module = k,
      geneset = rownames(meta.p),
      score = meta.rnk[,i],
      p.value = meta.p[,i],
      q.value = meta.q[,i],      
      pv
    )
    df <- df[order(-df$score),]
    if(!is.null(ntop) && ntop>0) df <- head(df,n=ntop)
    gse.list[[k]] <- df
  }
  
  ## add genes
  gse.genes <- list()
  k = names(gse.list)[1]
  for(k in names(gse.list)) {
    gset <- rownames(gse.list[[k]])
    gg <- wgcna$me.genes[[k]]
    if(!is.null(annot)) {
      gg <- probe2symbol(gg, annot, query = symbol.col)
    }
    set.genes <- lapply( gmt[gset], function(s) intersect(s, gg))
    n0 <- sapply(gmt[gset], length)
    n1 <- sapply(set.genes, length)
    gse.genes[[k]] <- sort(table(unlist(set.genes)),decreasing=TRUE)
    set.genes <- sapply(set.genes, function(g) paste(sort(g), collapse="|"))
    gse.list[[k]]$genes <- set.genes
    gse.list[[k]]$overlap <- paste0(n1,"/",n0)
  }

  ## compute gene enrichment
  message("[wgcna.compute_enrichment] computing gene enrichment...")
  gene.list <- list()
  k = names(gse.list)[1]
  for(k in names(gse.list)) {
    me <- gse.list[[k]]
    rnk <- me$score
    names(rnk) <- rownames(me)
    me.genes <- strsplit(me$genes,split='\\|')
    me.names <- unlist(mapply(rep, rownames(me), sapply(me.genes,length)))
    me.gmt <- tapply(me.names, unlist(me.genes), list)
    res <- fgsea::fgsea(me.gmt, rnk)
    res <- res[order(-res$NES),]
    gene.list[[k]] <- head(res,100)
  }
  
  res <- list(
    gsets = gse.list,
    genes = gene.list
  )
  
  return(res)
}


#' @export
wgcna.runConsensusWGCNA <- function( exprList,
                                    datTraits,
                                    power = 12,
                                    minKME = 0.8,
                                    mergeCutHeight = 0.25,
                                    deepSplit = 2
                                    ) {

  ##exprList <- list(Set1 = V1, Set2 = V2)
  multiExpr = WGCNA::list2multiData(lapply(exprList,Matrix::t))

  # The variable exprSize contains useful information about the sizes of all
  # of the data sets now we run automatic module detection procedure
  net.list <- list()
  for(k in names(multiExpr)) {
    message("[wgcna.runConsensusWGCNA] computing WGCNA for ", k)
    net.list[[k]] = WGCNA::blockwiseModules(
      datExpr = multiExpr[[k]]$data,
      maxBlockSize = 5000,
      power = power,
      networkType = "signed",
      TOMType = "signed",
      minModuleSize = 20,
      deepSplit = 2,
      mergeCutHeight = 0.25, 
      numericLabels = FALSE,
      minKMEtoStay = minKME,
      saveTOMs = FALSE,
      verbose = 0
    )
  }

  # now we run automatic consensus module detection 
  message("[wgcna.runConsensusWGCNA] computing consensus modules...")
  cons = WGCNA::blockwiseConsensusModules(
    multiExpr[],
    maxBlockSize = 5000,
    power = power, 
    networkType = "signed",
    TOMType = "signed",
    minModuleSize = 20,
    deepSplit = deepSplit,
    mergeCutHeight = mergeCutHeight, 
    numericLabels = FALSE,
    minKMEtoStay = minKME,
    saveTOMs = FALSE,
    verbose = 1
  )

  ii <- which(cons$goodGenes)
  colors <- sapply(net.list, function(net) net$colors)
  colors <- cbind( Consensus=cons$colors, colors )[ii,]
  
  if(0) {
    dendro = cons$dendrograms[[1]];
    plotDendroAndColors(
      dendro,
      colors = colors,
      c("Consensus", names(net.list)),
      dendroLabels = FALSE, hang = 0.03,
      addGuide = TRUE, guideHang = 0.05,
      main = "Consensus gene dendrogram and module colors")
  }

  ## create module-trait matrices for each set
  Z.list <- list()
  k=1
  for(k in names(cons$multiME)) {
    M <- (cons$multiME[[k]][[1]])
    Z.list[[k]] <- cor(M, datTraits[rownames(M),], use="pairwise")
  }
  z.check <- sapply(Z.list, function(z) colSums(z!=0,na.rm=TRUE)>0)
  sel <- names(which(rowMeans(z.check)==1))
  sel
  Z.list <- lapply(Z.list, function(z) z[,sel,drop=FALSE])
  Z.list

  ## create consensus module-trait matrix
  mdim <- sapply(exprList,ncol)
  pv <- mapply(function(z,n) corPvalueStudent(z,n), Z.list, mdim, SIMPLIFY=FALSE)
  all.sig <- Reduce('*', lapply(pv, function(p) 1*(p<0.05)))
  all.sig
  all.pos <- Reduce('*', lapply(Z.list, function(z) sign(z)==1))
  all.neg <- Reduce('*', lapply(Z.list, function(z) sign(z)==-1))  
  disconcordant <- !(all.pos | all.neg)
  disconcordant
  zsign <- sign(Reduce('+', lapply(Z.list, sign)))
  consZ <- Reduce(pmin, lapply(Z.list,abs)) * zsign
  consZ[disconcordant] <- NA  
  ##consZ[!all.sig] <- NA
  ydim <- sapply(exprList,ncol)
  
  res <- list(
    cons = cons,
    consZ = consZ,
    dendro = cons$dendrograms[[1]],    
    colors = colors,
    zlist = Z.list,
    ydim = ydim
  )
  
  res
}


#'
#'
#' @export
wgcna.getGeneStats <- function(wgcna, module, trait, plot=TRUE,
                               showallmodules=TRUE, col=NULL,
                               main=NULL) {
  # module="MEblue";trait="activated=act"
  p1 <- c("moduleMembership","MMPvalue")
  p2 <- c("traitSignificance","GSPvalue","foldChange","foldChangePvalue")
  p3 <- c("geneCentrality")
  nrow <- ncol(wgcna$datExpr)
  df <- wgcna$stats[[p1[1]]][,0]  ## empty data.frame
  rownames(df) <- colnames(wgcna$datExpr)
  if(!is.null(module)) {
    A1 <- sapply( wgcna$stats[p1], function(x) x[,module])
    df <- cbind(df, A1)
  }
  if(!is.null(trait)) {
    A2 <- sapply( wgcna$stats[p2], function(x) x[,trait])
    df <- cbind(df, A2)
  }
  A3 <- wgcna$stats[[p3]]
  df <- cbind(df, centrality=A3)

  sel <- c("moduleMembership","traitSignificance","foldChange","centrality")
  sel <- intersect(sel, colnames(df))
  df1 <- pmax( as.matrix(df[,sel]), 1e-8)
  score <- exp(rowMeans(log(df1)))

  labels <- wgcna$net$labels
  df <- data.frame( module = labels, score=score, df )
  df <- df[order(-df$score),]
  if(!is.null(module) && !showallmodules) {
    sel <- which( df$module == module )
    df <- df[sel,,drop=FALSE]
  }
  
  if(plot) {
    cols <- c("moduleMembership","traitSignificance","foldChange","centrality")
    cols <- intersect(cols, colnames(df))
    df1 <- df[,cols]
    col1 <- wgcna.labels2colors(wgcna$net$colors[rownames(df1)])
    if(!is.null(col)) col1 <- col
    pairs(df1, col=col1)
    if(is.null(main)) {
      main <- paste("Gene significance for module",module,"and trait",trait)
    }
    title(main, line=3, cex.main=1.15)
  }

  df
}

##=========================================================================
## PLOTTING FUNCTIONS
##=========================================================================

#'
#'
#' @export
wgcna.plotTOM <- function(wgcna, justdata=FALSE) {

  datExpr <- wgcna$datExpr
  MEs <- wgcna$net$MEs
  moduleColors <- wgcna.labels2colors(wgcna$net$colors)  

  ## Topological overlap dissimilarity matrix
  dissTOM <- 1 - wgcna$TOM
  rownames(dissTOM) <- colnames(dissTOM) <- colnames(datExpr)

  ## clustering wgcnaults
  geneTree <- wgcna$net$dendrograms[[1]]
  
  ## Taking the dissimilarity to a power, say 10, makes the plot
  ## more informative by effectively changing the color palette;
  ## setting the diagonal to NA also improves the clarity of the
  ## plot
  plotDiss <- dissTOM^7
  diag(plotDiss) <- NA
  myheatcol <- gplots::colorpanel(250, "red", "orange", "lemonchiffon")
  myheatcol <- gplots::colorpanel(250, "lemonchiffon", "orange", "red")

  if(justdata) return(plotDiss)
  
  par(oma = c(2, 0, 0, 0))
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
    plotDiss,
    geneTree,
    moduleColors,
    #col = myheatcol,
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

}


#'
#'
#' @export
wgcna.plotDendroAndColors <- function(wgcna, main=NULL, unmerged=FALSE) {

  colors <- wgcna.labels2colors(wgcna$net$colors)
  groupLabels <- "Module colors"
  geneTree = wgcna$net$dendrograms[[1]]
  if(unmerged) {
    colors <- cbind(colors,
      wgcna.labels2colors(wgcna$net$unmergedColors))
      groupLabels <- c( "Merged colors", "Unmerged colors")
  }
  
  if(is.null(main)) main <- "Gene dendrogram and module colors"
  ## Plot the dendrogram and the module colors underneath
  plotDendroAndColors(
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
  if(all(is.numeric(colors))) {
    colors <- WGCNA::labels2colors(colors, ...)
    return(colors)
  }
  stdColors <- c("grey",WGCNA::standardColors())
  if(all(colors %in% stdColors)) {
    return(colors)
  }
  icolors <- as.integer(factor(as.character(colors)))
  colors <- WGCNA::standardColors()[icolors]
  return(colors)
}


#' @export
wgcna.plotModuleTraitHeatmap <- function(wgcna, setpar=TRUE, cluster=FALSE,
                                         main = NULL, justdata=FALSE,
                                         pstar = TRUE) {
  
  ## Define numbers of genes and samples
  nGenes = ncol(wgcna$datExpr);
  nSamples = nrow(wgcna$datExpr);
  if("stats" %in% names(wgcna)) {
    moduleTraitCor = wgcna$stats$moduleTraitCor
    moduleTraitPvalue <- wgcna$stats$moduleTraitPvalue
  } else {
    moduleTraitCor = cor(wgcna$net$MEs, wgcna$datTraits, use = "pairwise.complete");
    moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
  }
  
  if(cluster) {
    ii <- hclust(dist(moduleTraitCor))$order
    jj <- hclust(dist(t(moduleTraitCor)))$order
  } else {
    ii <- rev(order(rownames(moduleTraitCor)))
    jj <- order(colnames(moduleTraitCor))  
  }
  moduleTraitCor <- moduleTraitCor[ii,jj]
  moduleTraitPvalue <- moduleTraitPvalue[ii,jj]    

  if(justdata) return(moduleTraitCor)
  
  ## Will display correlations and their p-values
  if(pstar) {
    textPv = cut(moduleTraitPvalue, breaks=c(-1,0.001,0.01,0.05,99),
                 labels = c("★★★","★★","★",""))
  } else {
    textPv = paste0("(",signif(moduleTraitPvalue, 1), ")")
  }
  
  textMatrix = paste0(signif(moduleTraitCor, 2), "\n", textPv)
  dim(textMatrix) = dim(moduleTraitCor)
  if(setpar) par(mar = c(8, 8, 3, 3));
  if(is.null(main)) main <- "Module-trait relationships"
  ## Display the correlation values within a heatmap plot
  WGCNA::labeledHeatmap(Matrix = moduleTraitCor,
                        xLabels = colnames(moduleTraitCor),
                        yLabels = rownames(moduleTraitCor),
                        ySymbols = rownames(moduleTraitCor),
                        colorLabels = FALSE,
                        colors = WGCNA::blueWhiteRed(50),
                        textMatrix = textMatrix,
                        setStdMargins = FALSE,
                        cex.text = 0.7,
                        zlim = c(-1,1),
                        main = main)
  
}

#' Plot membership correlation vs gene signficance (correlation with
#' trait) to discover biomarkers/driver genes.
#' 
#' @export
wgcna.plotMMvsGS <- function(wgcna, module, trait, abs=TRUE, par=TRUE,
                             plotlib = "base") {
  ##module="ME3";trait="activated=act"
  moduleGenes = wgcna$me.genes[[module]]
  nSamples = nrow(wgcna$datExpr)
  
  ## Module membership correlation (with p-values)
  if("stats" %in% names(wgcna)) {
    moduleMembership <- wgcna$stats$moduleMembership
    MMPvalue <- wgcna$stats$MMPvalue
  } else {
    moduleMembership = as.data.frame(cor(wgcna$datExpr, wgcna$net$MEs, use = "p"));
    MMPvalue = as.data.frame(corPvalueStudent(as.matrix(moduleMembership), nSamples));
  }

  ## Gene-trait significance (trait correlation) (with p-values)
  if("stats" %in% names(wgcna)) {
    traitSignificance <- wgcna$stats$traitSignificance
    GSPvalue <- wgcna$stats$GSPvalue
  } else {
    traitSignificance = as.data.frame(cor(wgcna$datExpr, wgcna$datTraits, use = "p"));
    GSPvalue = as.data.frame(corPvalueStudent(as.matrix(traitSignificance), nSamples));
  }
    
  x <- (moduleMembership[moduleGenes, module])
  y <- (traitSignificance[moduleGenes, trait])
  if(abs==TRUE) {
    x <- abs(x)
    y <- abs(y)
  }
  ## 
  px <- MMPvalue[moduleGenes, module]
  py <- GSPvalue[moduleGenes, trait]
  qx <- p.adjust(px, method="fdr")
  qy <- p.adjust(py, method="fdr")
  is.sig <- ( qx < 0.05 & qy < 0.05 )
  sigx <- (qx < 0.05)
  sigy <- (qy < 0.05)  
  ii <- which(is.sig)
  qv <- quantile(x[ii],prob=0.1)[1]
  qh <- quantile(y[ii],prob=0.1)[1]

  pos <- cbind(x,y)
  rownames(pos) <- moduleGenes
  is.sig1 <- c("notsig","onesig","sig")[1 + 1*sigx + 1*sigy]
  hi1 <- NULL
  ##hi1 <- head(rownames(pos),10)
  col1 <- c("grey70","grey20","red2")

  if(par) par(mfrow = c(1,1), mar=c(5,5,3,2));
  if(plotlib == "ggplot") {
    pgx.scatterPlotXY.GGPLOT(
      pos, var=is.sig1, hilight=hi1, col=col1,
      xlab = paste("Module membership in", module, "module"),
      ylab = paste("Gene significance for trait",trait),
      title = paste("Module membership vs. gene significance\n"),
      cex.title = 0.9, 
      girafe=FALSE)
  } else if(plotlib == "girafe") {
    pgx.scatterPlotXY.GGPLOT(
      pos, var=is.sig1, hilight=hi1, col=col1,
      xlab = paste("Module membership in", module, "module"),
      ylab = paste("Gene significance for trait",trait),
      title = paste("Module membership vs. gene significance\n"),
      cex.title = 0.7, cex.axis = 0.7,
      girafe=TRUE)
  } else {
    ii <- which( is.sig1 == "notsig")
    verboseScatterplot(
      x[-ii], y[-ii],
      xlab = paste("Module membership in", module, "module"),
      ylab = paste("Gene significance for trait",trait),
      main = paste("Module membership vs. gene significance\n"),
      cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = col1[1])
    ii <- which( is.sig1 == "onesig")
    points( x[ii], y[ii], col=col1[2] )
    ii <- which( is.sig1 == "sig")
    points( x[ii], y[ii], col=col1[3] )    
    abline(v=qv, h=qh, col="darkred")  
  }  
}


#' Plot cluster dendrogram with eigengenes and traits.
#' 
#' @export
wgcna.plotEigenGeneClusterDendrogram <- function(wgcna, add_traits=TRUE,
                                                 main = NULL) {
  # Matrix with eigengenes and traits
  MET <- wgcna$net$MEs
  if(add_traits) {
    MET <- cbind(MET, wgcna$datTraits)
  }
  MET = WGCNA::orderMEs(MET)
  if (NCOL(MET) <= 2) MET <- cbind(MET, MET) ## error if ncol(MET)<=2 !!!!  
  if(is.null(main)) main <- "Eigengene cluster dendrogram"
  WGCNA::plotEigengeneNetworks(
    MET, main,
    marDendro = c(0,4,2,0),
    plotHeatmaps = FALSE)

}

#' Plot the adjacency correlation heatmap matrix of eigengenes with or
#' without traits. This can show how traits cluster together with the
#' eigengenes.
#' 
#' @export
wgcna.plotEigenGeneAdjacencyHeatmap <- function(wgcna, add_traits=TRUE,
                                                marx=1, main=NULL,
                                                justdata=FALSE) {
  # Matrix with eigengenes and traits
  MET <- wgcna$net$MEs
  if(add_traits) {
    MET <- cbind(MET, wgcna$datTraits)
  }
  colnames(MET) <- paste0(" ",colnames(MET))
  MET = WGCNA::orderMEs(MET)
  if (NCOL(MET) <= 2) MET <- cbind(MET, MET) ## error if ncol(MET)<=2 !!!!

  if(justdata) {
    R <- (1 + cor(MET))/2
    return(R)
  }
  
  # Plot the correlation heatmap matrix (note: this plot will overwrite
  # the dendrogram plot)
  if(is.null(main)) main <- "Eigengene adjacency heatmap"
  WGCNA::plotEigengeneNetworks(
    MET, main,
    marHeatmap = c(8*marx,10*marx,2,2),
    plotDendrograms = FALSE,
    colorLabels = TRUE,
    xLabelsAngle = 45)
  
}

#' @export
wgcna.plotEigenGeneGraph <- function(wgcna, add_traits=TRUE, main=NULL) {
  ##require(igraph)
  net <- wgcna$net
  MET <- net$MEs
  if(add_traits) {
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
  
  ##par(mfrow = c(1, 1), mar = c(1, 1, 1, 1) * 0)
  igraph::plot.igraph(
    gr,
    layout = igraph::layout.kamada.kawai,
    vertex.label.cex = 0.8,
    edge.width = 3
  )
  if(!is.null(main)) title(main, line=-1.5)  
}


#' Plot Multi-dimensional scaling (MDS) of centered data matrix.
#'
#' @export
wgcna.plotMDS <- function(wgcna, main=NULL, scale=FALSE) {
  cc <- wgcna.labels2colors(wgcna$net$color)
  pc <- svd(t(scale(wgcna$datExpr,scale=scale)),nv=2)$u[,1:2]
  #pc <- svd(t(scale(wgcna$datExpr)),nv=1)$u[,1:2]
  colnames(pc) <- c("MDS-x","MDS-y")
  if(is.null(main)) main <- "MDS of features"
  plot(pc, col=cc, main=main)
}

#'
#'
#' @export
wgcna.plotFeatureUMAP <- function(wgcna, nhub=3, method="clust",
                                  scale=FALSE, main=NULL,
                                  plotlib="base") {

  if(method=="clust" && "clust" %in% names(wgcna)) {
    pos <- wgcna$clust[['umap2d']]
  } else if(method=="umap") {
    cX <- t(scale(wgcna$datExpr,scale=scale))  ## WGCNA uses correlation
    pos <- uwot::umap(cX)
    colnames(pos) <- c("UMAP-x","UMAP-y")
    rownames(pos) <- colnames(wgcna$datExpr)    
    ##  } else if(method=="mds") {
  } else {    
    pos <- svd(t(scale(wgcna$datExpr,scale=scale)),nv=2)$u[,1:2]
    colnames(pos) <- c("MDS-x","MDS-y")
    rownames(pos) <- colnames(wgcna$datExpr)
  }

  if(is.null(main)) main="Feature UMAP colored by module"

  ## get top hub genes
  mm <- wgcna$stats$moduleMembership
  hubgenes <- apply(mm, 2, function(x) head(names(sort(-x)),3),simplify=FALSE)
  hubgenes
  sel <- which(names(hubgenes)!="MEgrey")
  hubgenes <- unlist(hubgenes[sel])
  col1 <- wgcna$net$colors
  genes1 <- names(which(col1!="grey"))
  pgx.scatterPlotXY(
    pos, var=col1, col=sort(unique(col1)),
    hilight=genes1, hilight2=hubgenes,
    cex.lab=1.2, label.clusters = FALSE,
    title = main, 
    plotlib=plotlib)
    
}

#' Plot module significance.
#'
#' @export
wgcna.plotModuleSignificance <- function(wgcna, trait, main=NULL, abs=FALSE) {
  ##cc <- paste0("ME",wgcna$net$color)
  cc <- wgcna.labels2colors(wgcna$net$color)
  if("stats" %in% names(wgcna)) {
    traitSignificance <- wgcna$stats$traitSignificance
  } else {
    traitSignificance = as.data.frame(cor(wgcna$datExpr, wgcna$datTraits, use = "p"))
    names(traitSignificance) = names(wgcna$datTraits)
    rownames(traitSignificance) <- colnames(wgcna$datExpr)
  }
  geneSig <- traitSignificance[,trait]
  if(is.null(main)) main <- paste("Module significance with",trait)
  if(abs) geneSig <- abs(geneSig) 
  WGCNA::plotModuleSignificance(
    geneSig, colors=cc, main = main, boxplot=FALSE)
}

#'
#'
#' @export
wgcna.plotSampleDendroAndColors <- function(wgcna,
                                            what=c("me","traits","both")[3],
                                            clust.expr=TRUE,
                                            main=NULL, justdata=FALSE) {

  MET0 <- wgcna$net$MEs
  MET <- MET0[,0]
  if(any(what %in% c("me","both"))) {
     MET <- cbind(MET, MET0)
  }
  if(any(what %in% c("traits","both"))) {
    MET <- cbind(MET, wgcna$datTraits)
  }
  if (NCOL(MET) <= 2) MET <- cbind(MET, MET) ## error if ncol(MET)<=2 !!!!  
  
  ## Recalculate MEs with color as labels
  if(clust.expr) {
    sampleTree <- hclust(dist(wgcna$datExpr), method="average")
  } else {
    sampleTree <- hclust(dist(scale(MET0)), method="average")
  }
  ii <- sampleTree$order
  jj <- hclust(dist(t(scale(MET))))$order  
  colors <- WGCNA::numbers2colors(MET[,jj])

  if(justdata) {
    return(MET)
  }
  
  if(is.null(main)) {
    if(what=="me") main <- "Sample dendrogram and module heatmap"
    if(what=="traits") main <- "Sample dendrogram and trait heatmap"
    if(what=="both") main <- "Sample dendrogram and module+traits heatmap"    
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
                                                main = NULL, justdata=FALSE,
                                                pstar = TRUE) {
  
  ## Define numbers of genes and samples
  Pvalue = corPvalueStudent(R, nSamples);
  if(cluster) {
    R0 <- R
    R0[is.na(R0)] <- 0
    ii <- hclust(dist(R0))$order
    jj <- hclust(dist(t(R0)))$order
  } else {
    ii <- rev(order(rownames(R)))
    jj <- order(colnames(R))  
  }
  R <- R[ii,jj]
  Pvalue <- Pvalue[ii,jj]    

  if(justdata) return(R)
  
  ## Will display correlations and their p-values
  if(pstar) {
    textPv = cut(Pvalue, breaks=c(-1,0.001,0.01,0.05,99),
                 labels = c("★★★","★★","★",""))
  } else {
    textPv = paste0("(",signif(Pvalue, 1), ")")
  }
  
  textMatrix = paste0(signif(R, 2), "\n", textPv)
  textMatrix[which(is.na(R))] <- ""
  dim(textMatrix) = dim(R)
  if(setpar) par(mar = c(8, 8, 3, 3));
  if(is.null(main)) main <- "Correlation heatmap"
  ## Display the correlation values within a heatmap plot
  WGCNA::labeledHeatmap(Matrix = R,
                        xLabels = colnames(R),
                        yLabels = rownames(R),
                        xSymbols = colnames(R),
                        ySymbols = rownames(R),                        
                        colorLabels = FALSE,
                        colors = WGCNA::blueWhiteRed(50),
                        textMatrix = textMatrix,
                        setStdMargins = FALSE,
                        cex.text = 0.7,
                        zlim = c(-1,1),
                        main = main)

}

