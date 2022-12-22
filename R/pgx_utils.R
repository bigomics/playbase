logCPM <- function(counts, total=1e6, prior=1) {
  ## Transform to logCPM (log count-per-million) if total counts is
  ## larger than 1e6, otherwise scale to previous avarage total count.
  ##
  ##
  if(is.null(total)) {
    ##total <- nrow(counts)
    ##total <- mean(colSums(counts1!=0)) ## avg. number of expr genes
    total0 <- mean(Matrix::colSums(counts,na.rm=TRUE)) ## previous sum
    total <- ifelse( total0 < 1e6, total0, 1e6 )
    message("[logCPM] setting column sums to = ",round(total,2))
  }
  if(any(class(counts)=="dgCMatrix")) {
    ## fast/sparse calculate CPM
    cpm <- counts
    cpm[is.na(cpm)] <- 0  ## OK??
    cpm@x <- total * cpm@x / rep.int(Matrix::colSums(cpm), diff(cpm@p))  ## fast divide by columns sum
    cpm@x <- log2(prior + cpm@x)
    return(cpm)
  } else {
    cpm <- t(t(counts) / Matrix::colSums(counts,na.rm=TRUE)) * total
    x <- log2(prior + cpm)
    return(x)
  }
}

pgx.clusterSamples2 <- function(pgx, methods=c("pca","tsne","umap"), dims=c(2,3),
                                reduce.sd=1000, reduce.pca=50, perplexity=30,
                                center.rows=TRUE, scale.rows=FALSE,
                                X=NULL, umap.pkg="uwot", replace.orig=TRUE )
{
  if(!is.null(X)) {
    message("using provided X matrix...")
  } else if(!is.null(pgx$X)) {
    message("using pgx$X matrix...")
    X <- pgx$X
  } else {
    message("using logCPM(pgx$counts)...")
    ## X <- log2(1 + pgx$counts)
    X <- logCPM(pgx$counts, total=NULL)
  }
  dim(X)
  ##sdx <- apply(X,1,sd)
  ##X <- X[head(order(-sdx),reduce.sd),]
  ##if(center.rows) X <- X  - rowMeans(X)
  ##if(scale.rows)  X <- X / (1e-6+apply(X,1,sd))
  ##if(rank.tf) X <- scale(apply(X,2,rank))  ## works nicely

  dim(X)
  clust.pos <- pgx.clusterBigMatrix(
    X, methods = methods,
    dims = dims,
    perplexity = perplexity,
    center.features = center.rows,
    scale.features = scale.rows,
    reduce.sd = reduce.sd,
    reduce.pca = reduce.pca,
    find.clusters = FALSE,
    umap.pkg = umap.pkg
  )
  names(clust.pos)
  ##clust.index <- paste0("c",clust.pos$membership)
  clust.index <- clust.pos$membership
  clust.pos$membership <- NULL
  table(clust.index)

  if(replace.orig) {
    message("[pgx.clusterSamples2] update tsne2d/tsne3d and 'cluster' pheno...")
    pgx$samples$cluster <- clust.index
    pgx$tsne2d <- clust.pos[["tsne2d"]]
    pgx$tsne3d <- clust.pos[["tsne3d"]]
  } else {
    message("[pgx.clusterSamples2] skipping tsne2d/tsne3d update...")
  }

  pgx$cluster <- NULL
  pgx$cluster$pos  <- clust.pos
  pgx$cluster$index <- clust.index

  message("[pgx.clusterSamples2] done!")
  pgx
}


pgx.inferGender <- function(X, gene_name=NULL) {
  ## List of cell cycle markers, from Tirosh et al, 2015
  ##
  ##cc.genes <- readLines(con = "../opt/seurat/regev_lab_cell_cycle_genes.txt")
  if(is.null(gene_name)) gene_name <- toupper(sub(".*:","",rownames(X)))
  if(0) {
    X <- log2(1+ngs$counts)
    gene_name = rownames(X)
    x.genes <- rownames(ngs$genes)[grep("X",ngs$genes$chr)]
    y.genes <- rownames(ngs$genes)[grep("Y",ngs$genes$chr)]
    x.genes
    y.genes
  } else {
    y.genes = intersect(c("DDX3Y","RPS4Y1","USP9Y","KDM5D"),gene_name)
    y.genes
    x.genes = intersect(c("XIST"),gene_name)
    x.genes
  }
  if( length(y.genes)==0 && length(x.genes)==0 ) {
    cat("warning:: could not determine sex. missing some X/Y marker genes\n")
    sex <- rep(NA, ncol(X))
    return(sex)
  }
  sex <- rep(NA, ncol(X))
  if( length(y.genes)>0 && length(x.genes)>0 ) {
    x.expr <- colMeans(X[match(x.genes, gene_name),,drop=FALSE])
    y.expr <- colMeans(X[match(y.genes, gene_name),,drop=FALSE])
    x.expr
    y.expr
    mean.expr <- colMeans(X)
    ##plot(x.expr, y.expr)
    sex <- rep(NA, ncol(X))
    sex <- ifelse( x.expr > mean.expr & y.expr < mean.expr, "F", sex)
    sex <- ifelse( y.expr > mean.expr & x.expr < mean.expr, "M", sex)
    sex
    return(sex)
  }
  return(sex)
}

pgx.inferCellCyclePhase <- function(counts)
{
  ## List of cell cycle markers, from Tirosh et al, 2015
  ##
  ##cc.genes <- readLines(con = "../opt/seurat/regev_lab_cell_cycle_genes.txt")
  cc.genes = strsplit("MCM5 PCNA TYMS FEN1 MCM2 MCM4 RRM1 UNG GINS2 MCM6 CDCA7 DTL PRIM1 UHRF1 MLF1IP HELLS RFC2 RPA2 NASP RAD51AP1 GMNN WDR76 SLBP CCNE2 UBR7 POLD3 MSH2 ATAD2 RAD51 RRM2 CDC45 CDC6 EXO1 TIPIN DSCC1 BLM CASP8AP2 USP1 CLSPN POLA1 CHAF1B BRIP1 E2F8 HMGB2 CDK1 NUSAP1 UBE2C BIRC5 TPX2 TOP2A NDC80 CKS2 NUF2 CKS1B MKI67 TMPO CENPF TACC3 FAM64A SMC4 CCNB2 CKAP2L CKAP2 AURKB BUB1 KIF11 ANP32E TUBB4B GTSE1 KIF20B HJURP CDCA3 HN1 CDC20 TTK CDC25C KIF2C RANGAP1 NCAPD2 DLGAP5 CDCA2 CDCA8 ECT2 KIF23 HMMR AURKA PSRC1 ANLN LBR CKAP5 CENPE CTCF NEK2 G2E3 GAS2L3 CBX5 CENPA",split=" ")[[1]]
  s_genes <- cc.genes[1:43]
  g2m_genes <- cc.genes[44:97]

  ## Create our Seurat object and complete the initalization steps
  rownames(counts) <- toupper(rownames(counts))  ## mouse...
  ##counts1 <- cbind(counts,counts,counts,counts,counts,counts)
  obj <- Seurat::CreateSeuratObject(counts)
  obj <- Seurat::NormalizeData(obj, verbose=0)
  suppressWarnings( obj <- Seurat::CellCycleScoring(obj, s.features=s_genes,
                                                    g2m.features=g2m_genes, set.ident=TRUE))
  ## view cell cycle scores and phase assignments
  ##head(x = obj@meta.data)
  ##table(obj@meta.data$Phase)
  s.score <- obj@meta.data$S.Score
  g2m.score <- obj@meta.data$G2M.Score
  phase <- obj@meta.data$Phase
  if(is.null(phase) || length(phase)==0) return(NULL)
  return(phase)
}

compute.cellcycle.gender <- function(ngs, rna.counts=ngs$counts)
{
  pp <- rownames(rna.counts)
  is.mouse = (mean(grepl("[a-z]",gsub(".*:|.*\\]","",pp))) > 0.8)
  is.mouse
  if(!is.mouse) {
    if(1) {
      message("estimating cell cycle (using Seurat)...")
      ngs$samples$cell.cycle <- NULL
      ngs$samples$.cell.cycle <- NULL
      ##counts <- ngs$counts
      counts <- rna.counts
      rownames(counts) <- toupper(ngs$genes[rownames(counts),"gene_name"])
      res <- try(pgx.inferCellCyclePhase(counts) )  ## can give bins error
      if(class(res)!="try-error") {
        ngs$samples$.cell_cycle <- res
        table(ngs$samples$.cell_cycle)
      }
    }
    if(!(".gender" %in% colnames(ngs$samples) )) {
      message("estimating gender...")
      ngs$samples$.gender <- NULL
      X <- log2(1+rna.counts)
      gene_name <- ngs$genes[rownames(X),"gene_name"]
      ngs$samples$.gender <- pgx.inferGender( X, gene_name )
      table(ngs$samples$.gender)
    } else {
      message("gender already estimated. skipping...")
    }
    Matrix::head(ngs$samples)
  }
  return(ngs)
}

ngs.getGeneAnnotation <- function(genes) {
  hs.genes <- unique(unlist(as.list(org.Hs.eg.db::org.Hs.egSYMBOL)))
  mm.genes <- unique(unlist(as.list(org.Mm.eg.db::org.Mm.egSYMBOL)))

  ## if multiple genes, take first
  genes0 <- genes ## save
  genes <- sapply(genes, function(s) strsplit(s, split = "[;,\\|]")[[1]][1])

  is.human <- mean(genes %in% hs.genes) > mean(genes %in% mm.genes)
  is.mouse <- mean(genes %in% hs.genes) < mean(genes %in% mm.genes)

  txlen <- SYMBOL <- CHRLOC <- MAP <- NULL
  if (is.human) {
    GENE.TITLE <- unlist(as.list(org.Hs.eg.db::org.Hs.egGENENAME))
    SYMBOL <- unlist(as.list(org.Hs.eg.db::org.Hs.egSYMBOL))
    names(GENE.TITLE) <- SYMBOL
    CHRLOC <- as.list(org.Hs.eg.db::org.Hs.egCHRLOC)
    CHRLOC <- CHRLOC[match(names(SYMBOL), names(CHRLOC))]
    MAP <- as.list(org.Hs.eg.db::org.Hs.egMAP)
    MAP <- MAP[match(names(SYMBOL), names(MAP))]
    MAP <- sapply(MAP, paste, collapse = "|")
    names(CHRLOC) <- SYMBOL
    names(MAP) <- SYMBOL

    ## BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
    tx <- GenomicFeatures::transcriptLengths(txdb)
    TXLEN <- tapply(tx$tx_len, tx$gene_id, mean, na.rm = TRUE)
    TXLEN <- TXLEN[match(names(SYMBOL), names(TXLEN))]
    names(TXLEN) <- SYMBOL

    ## get gene biotype
    daf <- GenomicFeatures::transcripts(
      EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86,
      columns = c("gene_name", "gene_biotype"),
      return.type = "DataFrame"
    )
    GENE.BIOTYPE <- daf$gene_biotype
    names(GENE.BIOTYPE) <- daf$gene_name
  }
  if (is.mouse) {
    GENE.TITLE <- unlist(as.list(org.Mm.eg.db::org.Mm.egGENENAME))
    SYMBOL <- unlist(as.list(org.Mm.eg.db::org.Mm.egSYMBOL))
    names(GENE.TITLE) <- SYMBOL
    CHRLOC <- as.list(org.Mm.eg.db::org.Mm.egCHRLOC)
    CHRLOC <- CHRLOC[match(names(SYMBOL), names(CHRLOC))]
    names(CHRLOC) <- SYMBOL
    MAP <- NULL ## no map for mouse???

    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
    tx <- GenomicFeatures::transcriptLengths(txdb)
    TXLEN <- tapply(tx$tx_len, tx$gene_id, mean, na.rm = TRUE)
    TXLEN <- TXLEN[match(names(SYMBOL), names(TXLEN))]
    names(TXLEN) <- SYMBOL

    daf <- GenomicFeatures::transcripts(
      EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79,
      columns = c("gene_name", "gene_biotype"),
      return.type = "DataFrame"
    )
    GENE.BIOTYPE <- daf$gene_biotype
    names(GENE.BIOTYPE) <- daf$gene_name
  }
  remove(txdb)
  remove(tx)

  gene_title <- gene_biotype <- chrom <- loc <- map <- txlen <- NULL
  gene_title <- GENE.TITLE[genes]
  chrloc0 <- CHRLOC[genes]
  loc <- sapply(chrloc0, "[", 1) ## take only first entry !!!
  loc[sapply(loc, is.null)] <- NA
  loc <- abs(as.integer(unlist(loc)))
  chrom <- sapply(chrloc0, function(s) names(s)[1])
  chrom[sapply(chrom, is.null)] <- NA
  chrom <- as.vector(unlist(chrom))
  txlen <- round(TXLEN[genes])
  gene_biotype <- GENE.BIOTYPE[genes]

  map <- paste0("chr", chrom)
  map <- sub("chrNA", NA, map)
  if (!is.null(MAP)) map <- MAP[genes]

  ## get protein info
  annot <- data.frame(
    gene_name = genes,
    gene_title = gene_title,
    gene_biotype = gene_biotype,
    chr = chrom,
    pos = as.integer(loc),
    tx_len = as.integer(txlen),
    map = map
  )

  rownames(annot) <- genes
  annot
}

alias2hugo <- function(s, org = NULL, na.orig = TRUE) {
  hs.symbol <- unlist(as.list(org.Hs.eg.db::org.Hs.egSYMBOL))
  mm.symbol <- unlist(as.list(org.Mm.eg.db::org.Mm.egSYMBOL))
  if (is.null(org)) {
    is.human <- mean(s %in% hs.symbol, na.rm = TRUE) > mean(s %in% mm.symbol,
                                                            na.rm = TRUE)
    org <- ifelse(is.human, "hs", "mm")
  }
  nna <- which(!is.na(s) & s != "" & s != " ")
  s1 <- trimws(s[nna])
  hugo <- NULL
  if (org == "hs") {
    eg <- sapply(AnnotationDbi::mget(s1, envir = org.Hs.eg.db::org.Hs.egALIAS2EG,
                                     ifnotfound = NA), "[", 1)
    eg[is.na(eg)] <- "unknown"
    hugo <- sapply(AnnotationDbi::mget(eg, envir = org.Hs.eg.db::org.Hs.egSYMBOL,
                                       ifnotfound = NA), "[", 1)
  } else if (org == "mm") {
    eg <- sapply(AnnotationDbi::mget(s1, envir = org.Mm.eg.db::org.Mm.egALIAS2EG,
                                     ifnotfound = NA), "[", 1)
    eg[is.na(eg)] <- "unknown"
    hugo <- sapply(AnnotationDbi::mget(eg, envir = org.Mm.eg.db::org.Mm.egSYMBOL,
                                       ifnotfound = NA), "[", 1)
  } else {
    stop("[alias2hugo] invalid organism")
  }
  jj <- which(is.na(hugo))
  if (na.orig && length(jj)) hugo[jj] <- s1[jj]
  hugo0 <- rep(NA, length(s))
  hugo0[nna] <- hugo
  return(hugo0)
}

contrastAsLabels <- function(contr.matrix, as.factor = FALSE) {
  contrastAsLabels.col <- function(contr, contr.name) {
    grp1 <- gsub(".*[:]|_vs_.*", "", contr.name)
    grp0 <- gsub(".*_vs_|@.*", "", contr.name)
    x <- rep(NA, length(contr))
    x[which(contr < 0)] <- grp0
    x[which(contr > 0)] <- grp1
    if (as.factor) x <- factor(x, levels = c(grp0, grp1))
    x
  }
  K <- data.frame(contr.matrix[, 0])
  i <- 1
  for (i in 1:ncol(contr.matrix)) {
    contr <- contr.matrix[, i]
    contr.name <- colnames(contr.matrix)[i]
    k1 <- contrastAsLabels.col(contr, contr.name)
    K <- cbind(K, k1)
  }
  ## colnames(K) <- sub("[:].*",colnames(contr.matrix))
  colnames(K) <- colnames(contr.matrix)
  rownames(K) <- rownames(contr.matrix)
  K
}

makeDirectContrasts <- function(Y, ref, na.rm = TRUE) {
  ## check enough levels
  nlevel <- apply(Y, 2, function(y) length(unique(y)))
  if (any(nlevel < 2)) {
    notlevely <- colnames(Y)[which(nlevel < 2)]
    message("warning:: not enough levels: ", notlevely)
  }
  ii <- which(nlevel > 1)
  ii
  if (length(ii) == 0) {
    message("warning:: no valid phenotypes")
    return(NULL)
  }


  Y <- Y[, ii, drop = FALSE]
  ref <- ref[ii]

  ## make contrast
  exp.matrix <- makeDirectContrasts000(Y = Y, ref = ref,
                                       na.rm = na.rm, warn = FALSE)
  exp.matrix <- sign(exp.matrix)
  no.vs <- grep("_vs_|_VS_", colnames(exp.matrix), invert = TRUE)
  no.vs
  if (length(no.vs) > 0) {
    colnames(exp.matrix)[no.vs] <- paste0(colnames(exp.matrix)[no.vs],
                                          ":Y_vs_N")
  }
  exp.matrix0 <- exp.matrix
  if (all(grepl("_vs_|_VS_", colnames(exp.matrix0)))) {
    exp.matrix0 <- contrastAsLabels(exp.matrix0)
  }
  group <- pgx.getConditions(exp.matrix0)
  table(group)
  if (length(levels(group)) > 0.5 * nrow(exp.matrix)) {
    warning("contrast matrix looks degenerate.
        consider removing a contrast.\n")
  }

  contr.matrix <- exp.matrix[which(!duplicated(group)), , drop = FALSE]
  rownames(contr.matrix) <- group[which(!duplicated(group))]

  list(contr.matrix = contr.matrix, group = group, exp.matrix = exp.matrix)
}

makeDirectContrasts000 <- function(Y, ref, na.rm = TRUE, warn = FALSE) {
  if (NCOL(Y) == 1) Y <- data.frame(Y = Y)

  ## check
  all <- c("all", "other", "others", "rest")
  full <- c("*", "full")
  has.ref <- rep(NA, ncol(Y))
  for (i in 1:ncol(Y)) has.ref[i] <- (ref[i] %in% Y[, i] || ref[i] %in% c(all, full))
  has.ref
  if (!all(has.ref)) {
    stop("ERROR:: reference ", which(!has.ref), " not in phenotype matrix\n")
    return(NULL)
  }

  contr.matrix <- c()
  if (length(ref) < ncol(Y)) ref <- Matrix::head(rep(ref, 99), ncol(Y))
  ref.pattern <- "wt|contr|ctr|untreat|normal|^neg|ref|^no$|^0$|^0h$|scrambl|none|dmso|vehicle"
  i <- 1
  for (i in 1:ncol(Y)) {
    m1 <- NULL
    ref1 <- ref[i]
    x <- as.character(Y[, i])
    x[is.na(x) | x == "NA"] <- "_"
    detect.ref <- any(grepl(ref.pattern, x, ignore.case = TRUE))
    detect.ref
    if (is.na(ref1) && detect.ref) {
      ref1 <- grep(ref.pattern, x, ignore.case = TRUE, value = TRUE)
      ref1 <- sort(ref1)[1]
    }
    cref <- as.character(ref1)
    m1 <- model.matrix(~ 0 + x)
    colnames(m1) <- sub("^x", "", colnames(m1))
    if (ref1 %in% full) {
      levels <- names(table(x))
      levels <- setdiff(levels, c(NA, "NA"))
      levels
      if (length(levels) > 1) {
        cc <- makeFullContrasts(levels)
        m1 <- m1[, rownames(cc)] %*% cc
      }
    } else if (!is.na(ref1) && !(ref1 %in% all)) {
      m1 <- m1 - m1[, cref] ## +1/-1 encoding
      m1 <- m1[, which(colnames(m1) != cref), drop = FALSE] ## remove refvsref...
      m1 <- m1[, !colnames(m1) %in% c("NA", "_"), drop = FALSE]
      colnames(m1) <- paste0(colnames(m1), "_vs_", ref1)
    } else if (!is.na(ref1) && (ref1 %in% all)) {
      m1 <- t(t(m1 == 1) / Matrix::colSums(m1 == 1) - t(m1 == 0) / Matrix::colSums(m1 == 0))
      m1 <- m1[, !colnames(m1) %in% c("NA", "_"), drop = FALSE]
      colnames(m1) <- paste0(colnames(m1), "_vs_others")
    } else {
      stop("Fatal error with ref in making direct contrasts")
    }
    if (!is.null(m1)) {
      mm <- gsub("[: ]", "_", colnames(Y)[i])
      colnames(m1) <- paste0(mm, ":", colnames(m1))
      contr.matrix <- cbind(contr.matrix, m1)
    }
  }

  ## take out any empty comparisons
  contr.matrix <- contr.matrix[, which(Matrix::colSums(contr.matrix != 0) > 0), drop = FALSE]

  ## normalize to zero mean and symmetric sum-to-one. Any NA to zero.
  for (i in 1:ncol(contr.matrix)) {
    m <- contr.matrix[, i]
    m[is.na(m)] <- 0
    contr.matrix[, i] <- 1 * (m > 0) / sum(m > 0) - 1 * (m < 0) / sum(m < 0)
  }
  rownames(contr.matrix) <- rownames(Y)
  sign(contr.matrix)
}

#' Convert probe names to symbol names
#'
#' More explanation is needed for this function
#'
#' @param probes unknown.
#' @param type unknown.
#' @param org unknown.
#' @param keep.na unknown.
#'
#' @return unknown.
#'
#' @examples
#' x <- 1
.probe_to_symbol <- function(probes, type = NULL, org = "human", keep.na = FALSE) {
  ## strip postfix for ensemble codes
  if (mean(grepl("^ENS", probes)) > 0.5) {
    probes <- gsub("[.].*", "", probes)
  }

  if (is.null(type)) {
    hs.list <- list(
      "human.ensembl" = unlist(as.list(org.Hs.eg.db::org.Hs.egENSEMBL)),
      "human.ensemblTRANS" = unlist(as.list(org.Hs.eg.db::org.Hs.egENSEMBLTRANS)),
      "human.refseq" = unlist(as.list(org.Hs.eg.db::org.Hs.egREFSEQ)),
      "human.accnum" = unlist(as.list(org.Hs.eg.db::org.Hs.egACCNUM)),
      "human.uniprot" = unlist(as.list(org.Hs.eg.db::org.Hs.egUNIPROT)),
      "human.symbol" = unlist(as.list(org.Hs.eg.db::org.Hs.egSYMBOL))
    )

    mm.list <- list(
      "mouse.ensembl" = unlist(as.list(org.Mm.eg.db::org.Mm.egENSEMBL)),
      "mouse.ensemblTRANS" = unlist(as.list(org.Mm.eg.db::org.Mm.egENSEMBLTRANS)),
      "mouse.refseq" = unlist(as.list(org.Mm.eg.db::org.Mm.egREFSEQ)),
      "mouse.accnum" = unlist(as.list(org.Mm.eg.db::org.Mm.egACCNUM)),
      "mouse.uniprot" = unlist(as.list(org.Mm.eg.db::org.Mm.egUNIPROT)),
      "mouse.symbol" = unlist(as.list(org.Mm.eg.db::org.Mm.egSYMBOL))
    )

    id.list <- c(hs.list, mm.list)
    mx <- sapply(id.list, function(id) mean(probes %in% id))
    max.mx <- max(mx, na.rm = TRUE)
    mx0 <- names(mx)[which.max(mx)]
    org <- sub("[.].*", "", mx0)
    type <- sub(".*[.]", "", mx0)
    if (max.mx < 0.5 && max.mx > 0) {
      warning("Low mapping ratio: r= ", max.mx)
    }
    if (max.mx == 0) {
      warning("Zero mapping ratio: r= ")
      type <- NULL
    }
  }
  if (is.null(type)) {
    warning("invalid type: ", type, "\n")
    return(NULL)
  }
  if (!type %in% c(
    "ensembl", "ensemblTRANS", "unigene", "refseq",
    "accnum", "uniprot", "symbol"
  )) {
    warning("invalid type: ", type, "\n")
    return(NULL)
  }

  if (type == "symbol") {
    if (any(grep(" /// ", probes))) {
      symbol0 <- strsplit(probes, split = " /// ")
    } else if (any(grep("[;,]", probes))) {
      symbol0 <- strsplit(probes, split = "[;,\\|]")
    } else {
      symbol0 <- probes
    }
  } else {
    if (org == "human") {
      symbol0 <- AnnotationDbi::mapIds(
        org.Hs.eg.db::org.Hs.eg.db,
        probes, "SYMBOL", toupper(type)
      )
    }
    if (org == "mouse") {
      symbol0 <- AnnotationDbi::mapIds(
        org.Mm.eg.db::org.Mm.eg.db,
        probes, "SYMBOL", toupper(type)
      )
    }
  }

  ## Unrecognize probes
  nna <- which(is.na(names(symbol0)))
  if (length(nna)) names(symbol0)[nna] <- probes[nna]

  ## What to do with unmapped/missing symbols????
  symbol <- sapply(symbol0, "[", 1) ## takes first symbol only!!!
  isnull <- which(sapply(symbol, is.null))
  symbol[isnull] <- NA
  if (keep.na) {
    sel.na <- which(is.na(symbol))
    symbol[sel.na] <- probes[sel.na]
  }
  symbol <- unlist(symbol)
  names(symbol) <- NULL

  return(symbol)
}
