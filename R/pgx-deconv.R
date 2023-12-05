##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


#' @title Infer cell types from gene expression
#'
#' @param counts Numeric gene expression matrix with genes in rows.
#' @param low.th Filtering threshold to remove lowly expressed genes.
#' @param add.unknown Add "unknown" cell type if no match.
#' @param min.prob Minimum probability threshold for calling a cell type.
#' @param min.count Minimum count threshold for including a gene.
#' @param scalex Scale the expression values before inference.
#' @param normalize.mat Normalize the expression matrix before inference.
#' @param method Inference method to use (NNLM, LDA, ...).
#' @param collapse Method for collapsing predictions from multiple markers.
#' @param markers Cell type marker genes, if provided.
#'
#' @return Matrix of predicted cell type probabilities for each sample.
#'
#' @description Infer cell types based on known marker genes and their expression levels.
#'
#' @details This function identifies cell types based on the expression of known marker genes.
#' It filters the expression matrix, matches to the marker genes, and runs inference using NNLM or LDA.
#' The output is a matrix of predicted probabilities for each sample and cell type.
#' Parameters allow customizing the marker genes, expression filters, and inference method.
#'
#' @export
pgx.inferCellType <- function(counts, low.th = 0.01, add.unknown = FALSE,
                              min.prob = 0.2, min.count = 3, scalex = FALSE,
                              normalize.mat = TRUE, method = "NNLM",
                              collapse = "max", markers = NULL) {
  ## infer cell type from markers
  if (is.null(markers)) {
    message("[pgx.inferCellType] using database: 'signature-immuneMeta.csv'")
    M <- playdata::SIGNATURE_IMMUNEMETA
    M <- as.matrix(M)
  } else {
    message("[pgx.inferCellType] using provided markers list")
    marker.genes <- gsub("[-+]*$", "", unlist(markers))
    marker.genes <- sort(unique(marker.genes))
    M <- matrix(0, nrow = length(marker.genes), ncol = length(markers))
    dimnames(M) <- list(marker.genes, names(markers))
    k <- 7
    k <- 1
    for (k in 1:ncol(M)) {
      mm <- markers[[k]]

      m.neg <- sub("[-]$", "", grep("[-]$", mm, value = TRUE))
      m.pos <- setdiff(gsub("[+-]$", "", mm), m.neg)
      if (length(m.pos)) M[match(m.pos, rownames(M)), k] <- +1
      if (length(m.neg)) M[match(m.neg, rownames(M)), k] <- -1
    }
    M <- pmax(M, 0) ## sorry.. no negative markers yet in this function
    M <- as.matrix(M)
  }

  ## Filter count matrix
  X <- counts
  X <- X[(rowMeans(X >= min.count) > low.th), ] ## OK???

  ## Match matrices
  rownames(X) <- toupper(rownames(X))
  rownames(M) <- toupper(rownames(M))
  gg <- intersect(rownames(M), rownames(X))
  X1 <- X[gg, ]
  M1 <- M[gg, ]
  M1 <- M1[, Matrix::colSums(M1 != 0) > 0, drop = FALSE]
  if (scalex) X1 <- X1 / (1 + rowMeans(X1))

  ## run deconvolution algorithm

  out <- pgx.deconvolution(X1,
    ref = M1, methods = method,
    add.unknown = add.unknown, normalize.mat = normalize.mat
  )

  P <- out$results[[method]] ## choose specified method
  rownames(P) <- colnames(X1)
  P[is.na(P)] <- 0

  ## Collapse to single cell.type
  if (collapse == "sum") {
    P <- tapply(1:ncol(P), colnames(P), function(i) rowSums(P[, i, drop = FALSE]))
  } else if (collapse == "mean") {
    P <- tapply(1:ncol(P), colnames(P), function(i) rowMeans(P[, i, drop = FALSE]))
  } else {
    P <- tapply(1:ncol(P), colnames(P), function(i) apply(P[, i, drop = FALSE], 1, max))
  }
  P <- do.call(cbind, P)

  ## Get maximum probability cell.type
  P <- P / (1e-6 + rowSums(P, na.rm = TRUE))
  pmax <- apply(P, 1, max, na.rm = TRUE)
  sel.dodgy <- (pmax < min.prob) ## not reliable
  celltype <- colnames(P)[max.col(P)]
  if (any(sel.dodgy)) celltype[which(sel.dodgy)] <- NA

  ## collapse small groups to 'other_cells'
  low.ct <- names(which(table(celltype) < low.th * length(celltype)))
  celltype[celltype %in% low.ct] <- "other_cells"
  names(celltype) <- colnames(counts)

  res <- list(celltype = celltype, probs = P)
  return(res)
}






#' @title Simplify cell types names
#'
#' @description Simplifies cell types names from LM22, DICE, ImmProt and
#' ImmunoStates to standardized classification names.
#'
#' @param ct Character vector of cell type names.
#' @param low.th Threshold to filter low abundance cell types.
#'
#' @details This function takes a character vector of cell type names
#' from deconvolution methods like LM22, DICE, ImmProt, ImmunoStates etc.
#' It simplifies the names using pattern matching to map to standardized
#' cell type categories like "T_cells", "B_cells", "Macrophages" etc.
#'
#' Cell types below the abundance threshold `low.th` are removed.
#'
#' @return A character vector with simplified cell type names.
#'
#' @examples
#' \dontrun{
#' ct <- c("B.cells.Naive", "T.cells.CD8", "Macrophages.M2")
#' ct_simple <- pgx.simplifyCellTypes(ct)
#' print(ct_simple)
#' }
#' @export
pgx.simplifyCellTypes <- function(ct, low.th = 0.01) {
  ## Simplifies cell types names from LM22, DICE, ImmProt and
  ## ImmunoStates to standardized classification names.
  ##

  ## LM22
  ct[grep("^B.cells", ct)] <- "B_cells"
  ct[grep("^T.cells", ct)] <- "T_cells"
  ct[grep("^Macrophages", ct)] <- "Macrophages"
  ct[grep("^NK", ct)] <- "NK_cells"
  ct[grep("^Dendritic", ct)] <- "Dendritic_cells"
  ct[grep("^Mast", ct)] <- "Mast_cells"
  ct[grep("^Plasma", ct)] <- "Plasma_cells"
  ct[grep("Eosinophils|Neutrophils", ct)] <- "Granulocytes"

  ## DICE
  ct[grep("CD4|CD8|TH[12]|TFH|TREG|THSTAR", ct)] <- "T_cells"
  ct[grep("^M2$", ct)] <- "Macrophages"
  ct[grep("^B_CELL", ct)] <- "B_cells"
  ct[grep("^MONOCYTES$", ct)] <- "Monocytes"

  ## ImmProt
  ct[grep("^Bmem|^Bnav|^Bplasma", ct)] <- "B_cells"
  ct[grep("^T4|^Th17|^T8|Th[12]|Tregs", ct)] <- "T_cells"
  ct[grep("^NK", ct)] <- "NK_cells"
  ct[grep("^mDC|^pDC", ct)] <- "Dendritic_cells"
  ct[grep("^M2$", ct)] <- "Macrophages"
  ct[grep("^MOim|^MOnc|^MOcl", ct)] <- "Monocytes"
  ct[grep("^BSph|^ESph|^NTph", ct)] <- "Granulocytes"

  ## ImmunoStates
  ct[grep("B_cell", ct)] <- "B_cells"
  ct[grep("T_cell", ct)] <- "T_cells"
  ct[grep("macrophage", ct)] <- "Macrophages"
  ct[grep("natural.killer", ct)] <- "NK_cells"
  ct[grep("dendritic", ct)] <- "Dendritic_cells"
  ct[grep("plasma", ct)] <- "Plasma_cells"
  ct[grep("MAST", ct)] <- "Mast_cells"
  ct[grep("monocyte", ct)] <- "Monocytes"
  ct[grep("eosinophil|neutrophil|basophil", ct)] <- "Granulocytes"
  ct[grep("hemato", ct)] <- "other_cells"

  ## otherCells (from EPIC)
  ct[grep("otherCells", ct)] <- "other_cells"

  ## collapse low frequency celltype to "other"
  low.ct <- names(which(table(ct) < low.th * length(ct)))
  ct[ct %in% low.ct] <- "other_cells"

  return(ct)
}


#' @title Purify bulk gene expression profiles
#'
#' @param expr Gene expression matrix with genes in rows and samples
#' in columns.
#' @param ref Index of reference samples for estimating contamination.
#' @param k Number of components for NMF.
#' @param method NMF method, either 1 for normal cell contamination or
#' 2 for tumor cell purification.
#'
#' @return A list containing:
#' \itemize{
#'   \item{\code{purified}}{Purified gene expression matrix.}
#'   \item{\code{contaminant}}{Estimated contaminating profile.}
#'   \item{\code{alpha}}{Estimated sample purity.}
#' }
#'
#' @description
#' Estimates and removes contaminating gene expression signals from bulk
#' RNA-seq data using non-negative matrix factorization (NMF).
#'
#' @details This function takes a gene expression matrix and reference samples.
#' It uses NMF to estimate the contaminating expression profile
#' and sample purity. Method 1 assumes the contaminant is normal cells, method 2
#' assumes it is tumor cells. The estimated contaminating profile is subtracted
#' to obtain a purified expression profile.
#'
#' The reference samples should represent the contaminating component (normal or tumor).
#' The number of components k determines the complexity of the estimated contaminating
#' profile. A purity vector and the contaminating profile matrix are returned along with
#' the purified data.
#'
#' @examples
#' \dontrun{
#' # Purify normal cell contamination
#' purified <- pgx.purify(expr, ref = 1:5, method = 1)
#'
#' # Purify tumor cell contamination
#' purified <- pgx.purify(expr, ref = 6:10, method = 2)
#' }
#' @export
pgx.purify <- function(X, ref, k = 3, method = 2) {
  ## ----------------------------------------------------------------------
  ## Using NNLM
  ## ----------------------------------------------------------------------

  if (method == 1) {
    ## contaminating profiles (e.g. stromal, normal cells)
    normalX <- cbind(const = 1, X[, ref, drop = FALSE]) ## contamination (e.g. stromal, normal cells)
    colnames(normalX) <- paste0("N", 1:ncol(normalX))

    ## compute proportion of tumour content using NNMF
    res.nmf <- NNLM::nnmf(X, k = k, init = list(W0 = normalX), check.k = FALSE)
    alpha <- with(res.nmf, Matrix::colSums(W[, 1:k] %*% H[1:k, ]) / Matrix::colSums(W %*% H))

    ## estimate "pure" matrix
    x.total <- res.nmf$W[, ] %*% res.nmf$H[, ]
    x.purified <- res.nmf$W[, 1:k] %*% res.nmf$H[1:k, ]

    x.contaminant <- pmax(X - x.purified, 0)
  } else if (method == 2) {
    jj <- setdiff(1:ncol(X), ref)
    tumorX <- cbind(const = 1, X[, jj, drop = FALSE])
    colnames(tumorX) <- paste0("T", 1:ncol(tumorX))

    ## compute proportion of contaminant content using NNMF
    res.nmf <- NNLM::nnmf(X, k = k, init = list(W0 = tumorX), check.k = FALSE)
    beta <- with(res.nmf, Matrix::colSums(W[, 1:k] %*% H[1:k, ]) / Matrix::colSums(W %*% H))
    alpha <- (1 - beta)

    ## estimate "pure" matrix
    x.contaminant <- res.nmf$W[, 1:k] %*% res.nmf$H[1:k, ]
    x.purified <- pmax(X - x.contaminant, 0)
  } else {
    dbg("fatal error:: unknown method\n")
    stop()
  }

  res <- list(
    purified = x.purified,
    contaminant = x.contaminant,
    alpha = alpha
  )
  return(res)
}

#' @title Infer Cell Cycle Phase
#'
#' @param counts Gene expression matrix with genes in rows and cells in columns.
#'
#' @return Seurat object containing cell cycle scores and phase assignments.
#'
#' @description Infers cell cycle phase based on expression of known marker genes.
#'
#' @details This function takes a gene expression matrix and calculates cell cycle
#' scores and assigns cell cycle phases using the Seurat package. It uses curated lists
#' of S phase and G2/M phase marker genes from Tirosh et al. 2015. The gene expression
#' matrix is normalized then scored based on the phase marker genes. The output Seurat
#' object contains the scores and assigned phases for each cell.
#'
#' @export
pgx.inferCellCyclePhase <- function(counts) {
  ## List of cell cycle markers, from Tirosh et al, 2015
  ##

  rownames(counts) <- toupper(rownames(counts)) ## mouse...
  cc.genes <- strsplit("MCM5 PCNA TYMS FEN1 MCM2 MCM4 RRM1 UNG GINS2 MCM6 CDCA7 DTL PRIM1 UHRF1 MLF1IP HELLS RFC2 RPA2 NASP RAD51AP1 GMNN WDR76 SLBP CCNE2 UBR7 POLD3 MSH2 ATAD2 RAD51 RRM2 CDC45 CDC6 EXO1 TIPIN DSCC1 BLM CASP8AP2 USP1 CLSPN POLA1 CHAF1B BRIP1 E2F8 HMGB2 CDK1 NUSAP1 UBE2C BIRC5 TPX2 TOP2A NDC80 CKS2 NUF2 CKS1B MKI67 TMPO CENPF TACC3 FAM64A SMC4 CCNB2 CKAP2L CKAP2 AURKB BUB1 KIF11 ANP32E TUBB4B GTSE1 KIF20B HJURP CDCA3 HN1 CDC20 TTK CDC25C KIF2C RANGAP1 NCAPD2 DLGAP5 CDCA2 CDCA8 ECT2 KIF23 HMMR AURKA PSRC1 ANLN LBR CKAP5 CENPE CTCF NEK2 G2E3 GAS2L3 CBX5 CENPA", split = " ")[[1]]
  cc.genes <- cc.genes[cc.genes %in% rownames(counts)] # Subset for small datasets
  s_genes <- cc.genes[1:43]
  g2m_genes <- cc.genes[44:97]
  ## Create our Seurat object and complete the initalization steps

  obj <- Seurat::CreateSeuratObject(counts)
  obj <- Seurat::NormalizeData(obj, verbose = 0)
  suppressWarnings(obj <- Seurat::CellCycleScoring(obj,
    s.features = s_genes,
    g2m.features = g2m_genes,
    set.ident = TRUE
  ))
  ## view cell cycle scores and phase assignments
  s.score <- obj@meta.data$S.Score
  g2m.score <- obj@meta.data$G2M.Score
  phase <- obj@meta.data$Phase
  if (is.null(phase) || length(phase) == 0) {
    return(NULL)
  }
  return(phase)
}


#' @title Score cell cycle phase from gene expression
#'
#' @param expr Gene expression matrix with genes in rows and cells in columns.
#'
#' @return DataFrame with cell cycle phase scores and assignments for each cell.
#'
#' @description
#' Scores cell cycle phase for single cells based on expression of known marker genes.
#'
#' @details This function takes a gene expression matrix and calculates cell cycle phase
#' scores and assignments for each cell using curated gene lists. It calculates S and G2M
#' scores based on phase-specific marker genes from Tirosh et al. 2015. The difference between
#' S and G2M scores is used to assign cell cycle phases (G1, S or G2M).
#'
#' The output data frame contains the scores and assigned phases for each cell. This provides
#' an overview of cell cycle state in the data based on expression of canonical markers.
#'
#' @export
pgx.scoreCellCycle <- function(counts) {
  ## List of cell cycle markers, from Tirosh et al, 2015
  ##

  cc.genes <- strsplit("MCM5 PCNA TYMS FEN1 MCM2 MCM4 RRM1 UNG GINS2 MCM6 CDCA7 DTL PRIM1 UHRF1 MLF1IP HELLS RFC2 RPA2 NASP RAD51AP1 GMNN WDR76 SLBP CCNE2 UBR7 POLD3 MSH2 ATAD2 RAD51 RRM2 CDC45 CDC6 EXO1 TIPIN DSCC1 BLM CASP8AP2 USP1 CLSPN POLA1 CHAF1B BRIP1 E2F8 HMGB2 CDK1 NUSAP1 UBE2C BIRC5 TPX2 TOP2A NDC80 CKS2 NUF2 CKS1B MKI67 TMPO CENPF TACC3 FAM64A SMC4 CCNB2 CKAP2L CKAP2 AURKB BUB1 KIF11 ANP32E TUBB4B GTSE1 KIF20B HJURP CDCA3 HN1 CDC20 TTK CDC25C KIF2C RANGAP1 NCAPD2 DLGAP5 CDCA2 CDCA8 ECT2 KIF23 HMMR AURKA PSRC1 ANLN LBR CKAP5 CENPE CTCF NEK2 G2E3 GAS2L3 CBX5 CENPA", split = " ")[[1]]
  s_genes <- cc.genes[1:43]
  g2m_genes <- cc.genes[44:97]

  ## Create our Seurat object and complete the initalization steps
  rownames(counts) <- toupper(rownames(counts)) ## mouse...

  obj <- Seurat::CreateSeuratObject(counts)
  obj <- Seurat::NormalizeData(obj, verbose = 0)
  suppressWarnings(obj <- Seurat::CellCycleScoring(obj, s_genes, g2m_genes, set.ident = TRUE))

  ## view cell cycle scores and phase assignments
  s.score <- obj@meta.data$S.Score
  g2m.score <- obj@meta.data$G2M.Score
  diff.score <- s.score - g2m.score
  phase <- obj@meta.data$Phase
  df <- data.frame(
    phase = phase, s_score = s.score, g2m_score = g2m.score,
    diff_score = diff.score
  )
  if (nrow(df) == 0) {
    return(NULL)
  }
  return(df)
}


#' @title Infer sample gender from gene expression
#'
#' @param expr Numeric gene expression matrix with genes in rows.
#' @param method Method for gender inference ("pca" or "ttest").
#' @param pc Number of principal components to use if method="pca".
#'
#' @return Vector of predicted gender labels ("F" or "M") for each sample.
#'
#' @description Infers sample gender based on gene expression of Y chromosome genes.
#'
#' @details This function takes a gene expression matrix and infers the gender of each
#' sample using either a t-test (default) or PCA approach. The t-test compares expression
#' of Y chromosome genes between groups of known male and female samples. The PCA approach
#' uses the first few principal components of Y gene expression as a gender signature.
#'
#' The output is a vector of predicted gender labels "F" or "M" for each sample.
#'
#' @export
pgx.inferGender <- function(X, gene_name = NULL) {
  ## List of cell cycle markers, from Tirosh et al, 2015
  ##

  if (is.null(gene_name)) gene_name <- toupper(sub(".*:", "", rownames(X)))
  y.genes <- intersect(c("DDX3Y", "RPS4Y1", "USP9Y", "KDM5D"), gene_name)
  x.genes <- intersect(c("XIST"), gene_name)
  if (length(y.genes) == 0 && length(x.genes) == 0) {
    dbg("warning:: could not determine sex. missing some X/Y marker genes\n")
    sex <- rep(NA, ncol(X))
    return(sex)
  }
  sex <- rep(NA, ncol(X))
  if (length(y.genes) > 0 && length(x.genes) > 0) {
    x.expr <- colMeans(X[match(x.genes, gene_name), , drop = FALSE])
    y.expr <- colMeans(X[match(y.genes, gene_name), , drop = FALSE])
    mean.expr <- colMeans(X)

    sex <- rep(NA, ncol(X))
    sex <- ifelse(x.expr > mean.expr & y.expr < mean.expr, "F", sex)
    sex <- ifelse(y.expr > mean.expr & x.expr < mean.expr, "M", sex)
    return(sex)
  }
  return(sex)
}


#' @title Run multiple deconvolution methods
#'
#' @param counts Numeric gene expression matrix with genes in rows and samples in columns.
#' @param refmat List of reference expression matrices or a single matrix.
#' @param methods Character vector of deconvolution methods to run.
#'
#' @return List containing results from each deconvolution method.
#'
#' @description Runs multiple deconvolution methods using different reference profiles.
#'
#' @details This function runs deconvolution across multiple methods like CIBERSORT, DeconRNASeq, EPIC, etc.
#' It takes a gene expression matrix and one or more reference matrices as input. The reference matrix/matrices
#' represent cell type-specific expression profiles.
#'
#' Deconvolution is performed for each reference separately using the specified methods. The results
#' from each method are returned as nested lists, with the outer list split by reference profile name.
#'
#' @export
pgx.multipleDeconvolution <- function(counts, refmat,
                                      methods = c(
                                        "I-NNLS", "CIBERSORT", "DCQ", "DeconRNAseq", "EPIC", "NNLM",
                                        "cor", "SingleR"
                                      )) {
  timings <- c()
  results <- list()

  if (!inherits(refmat, "list")) {
    refmat <- list("reference" = refmat)
  }

  refnames <- names(refmat)
  i <- 1
  for (i in 1:length(refmat)) {
    message("[pgx.multipleDeconvolution] computing for ", refnames[i])
    ref <- refmat[[i]]
    res <- pgx.deconvolution(counts, ref = ref, methods = methods)

    if (!is.null(res)) {
      m <- names(refmat)[i]
      results[[m]] <- res$results
      timings <- rbind(timings, res$timings)
    }
  }

  timings0 <- NULL
  if (!is.null(timings)) {
    timings0 <- apply(timings, 2, function(x) tapply(x, rownames(timings), sum))
  }
  res2 <- list(results = results, timings = timings0)
  return(res2)
}


#' @title Deconvolve cell type proportions from bulk RNA-seq data
#'
#' @description Deconvolves bulk RNA-seq data into cell type proportions
#' using reference gene expression signatures.
#'
#' @param pgx A PGX object containing bulk RNA-seq data.
#' @param ref A reference gene expression matrix with cell types in columns.
#' @param method Deconvolution method to use. Options are NNLM, CIBERSORT, etc.
#' @param add.unknown Logical for whether to add an 'unknown' cell type.
#' @param normalize.mat Logical for whether to normalize input matrices.
#'
#' @details This function takes a PGX object containing bulk RNA-seq data, plus a
#' reference gene expression matrix with cell types in columns. It runs
#' deconvolution to estimate the proportions of each cell type in the bulk
#' samples.
#'
#' The \code{method} parameter specifies the algorithm, such as NNLM or CIBERSORT.
#' An additional 'unknown' cell type can be added to soak up residuals.
#' Input matrices can be normalized before deconvolution with \code{normalize.mat}.
#'
#' @return The function returns a matrix of estimated cell type proportions,
#' with bulk RNA-seq samples in rows and cell types in columns.
#'
#' @export
pgx.deconvolution <- function(X, ref,
                              methods = c(
                                "I-NNLS", "CIBERSORT", "DCQ", "DeconRNAseq", "EPIC", "NNLM",
                                "cor", "SingleR"
                              ),
                              add.unknown = FALSE, normalize.mat = TRUE) {
  if (max(X) < 50 || min(X) < 0) {
    dbg("WARNING:: pgx.deconvolution: is X really counts? (not logarithmic)\n")
  }


  ## clean up matrix, remove duplicate names
  mat <- as.matrix(X)
  rownames(mat) <- gsub(".*:", "", rownames(mat)) ## strip prefix
  rownames(mat) <- toupper(rownames(mat)) ## handle mouse??
  mat <- mat[order(-rowMeans(mat)), , drop = FALSE]
  mat <- as.matrix(mat[!duplicated(rownames(mat)), , drop = FALSE])

  ref <- as.matrix(ref)
  rownames(ref) <- toupper(rownames(ref))
  ref <- ref[order(-rowMeans(ref)), , drop = FALSE]
  ref <- as.matrix(ref[!duplicated(rownames(ref)), , drop = FALSE])

  ## Add "unknown" class to reference matrix
  if (add.unknown) {
    gg <- intersect(rownames(ref), rownames(mat))
    x1 <- log(1 + ref[gg, , drop = FALSE])
    y1 <- log(1 + rowMeans(mat[gg, , drop = FALSE]))
    x1 <- cbind(offset = 1, x1)

    ## compute residual matrix by substracting all possible linear
    ## combinations of reference.

    x1 <- ref[gg, , drop = FALSE]
    y1 <- rowMeans(mat[gg, , drop = FALSE])
    cf <- NNLM::nnlm(x1, cbind(y1))$coefficients
    cf[is.na(cf)] <- 0
    resid <- pmax(y1 - x1 %*% cf, 0) ## residual vector
    resx <- rep(0, nrow(ref))
    names(resx) <- rownames(ref)
    resx[gg] <- resid
    ref <- cbind(ref, "other_cells" = resx)
  }

  ## normalize all matrices to CPM
  if (normalize.mat) {
    ref <- t(t(ref) / sqrt(1e-6 + Matrix::colSums(ref**2, na.rm = TRUE))) * 1e6
    mat <- t(t(mat) / sqrt(1e-6 + Matrix::colSums(mat**2, na.rm = TRUE))) * 1e6
  }

  ## add small noise, some methods need it...
  ref <- ref + 1e-2 * matrix(stats::rnorm(length(ref)), nrow(ref), ncol(ref))
  mat <- mat + 1e-2 * matrix(stats::rnorm(length(mat)), nrow(mat), ncol(mat))
  ref <- pmax(ref, 0)
  mat <- pmax(mat, 0)

  gg <- intersect(rownames(ref), rownames(mat))
  if (length(gg) < 10) {
    warning("WARNING:: pgx.deconvolution: no enough marker genes")
    return(NULL)
  }

  ## conform??

  timings <- list()
  results <- list()
  CIBERSORT.code <- "/opt/CIBERSORT/CIBERSORTmat.R"
  if ("CIBERSORT" %in% methods && file.exists(CIBERSORT.code)) {
    ## CIBERSORT
    source(CIBERSORT.code)
    dbg("starting deconvolution using CIBERSORT...\n")
    ciber.out <- NULL
    stime <- system.time(
      # The f CIBERSORT is likely from CIBERSORT package but better confirm
      try(ciber.out <- CIBERSORT(ref, mat, perm = 0, QN = FALSE))
    )
    if (!is.null(ciber.out)) {
      timings[["CIBERSORT"]] <- stime
      dbg("deconvolution using CIBERSORT took", stime[3], "s\n")
      ciber.out <- ciber.out[, !(colnames(ciber.out) %in% c("P-value", "Correlation", "RMSE"))]
      results[["CIBERSORT"]] <- ciber.out
    } else {
      dbg("WARNING:: CIBERSORT failed\n")
    }
  }

  if ("EPIC" %in% methods) {
    ## EPIC

    dbg("[pgx.deconvolution] calculating EPIC...")

    out <- NULL
    gg <- intersect(rownames(ref), rownames(mat))
    ref1 <- ref
    if ("other_cells" %in% colnames(ref1)) {
      ref1 <- ref1[, setdiff(colnames(ref1), "other_cells")]
    }
    ref.list <- list(refProfiles = ref1, sigGenes = gg)
    mat1 <- mat
    colnames(mat1) <- 1:ncol(mat) ## EPIC doesnt like duplicated column names...
    stime <- system.time(try(
      out <- EPIC::EPIC(bulk = mat1, reference = ref.list)
    ))
    remove(mat1)
    if (!is.null(out)) {
      dbg(paste0("deconvolution using EPIC took", stime[3], "s\n"))
      timings[["EPIC"]] <- stime
      out.mat <- out[["cellFractions"]]
      colnames(out.mat) <- sub("otherCells", "other_cells", colnames(out.mat))
      rownames(out.mat) <- colnames(mat)
      results[["EPIC"]] <- out.mat
    } else {
      dbg("WARNING:: EPIC fail (no pun intended...)\n")
    }
  }

  if (FALSE && "DeconRNAseq" %in% methods) {
    ## IK17.04.2023 ************ BROKEN *******************
    ## ---- needs psych & pcaMethods inside namespace----
    ## DeconRNAseq
    if ("package:Seurat" %in% search()) detach("package:Seurat", unload = TRUE)
    dbg("[pgx.deconvolution] calculating DeconRNAseq...")

    ## DeconRNASeq need psych and pcaMethods, so we temporarily
    ## load the library...

    drs <- NULL
    stime <- system.time(suppressMessages(suppressWarnings(
      drs <- try(DeconRNASeq::DeconRNASeq(
        data.frame(mat, check.names = FALSE),
        data.frame(ref, check.names = FALSE)
      )$out.all)
    )))
    ## ... and quickly remove these


    timings[["DeconRNAseq"]] <- stime
    if (!is.null(drs) && !inherits(drs, "try-error")) {
      dbg("deconvolution using DeconRNAseq took", stime[3], "s\n")
      rownames(drs) <- colnames(mat)
      results[["DeconRNAseq"]] <- drs
    } else {
      dbg("*** WARNING ***:: DeconRNAseq failed\n")
    }
  }

  ## --------- DCQ from ComICS ---------------------
  if ("DCQ" %in% methods) {
    ## DCQ seems to work in logX, so we use log-transform

    dbg("[pgx.deconvolution] calculating DCQ...")

    res.dcq <- NULL
    stime <- system.time(
      res.dcq <- try(
        ComICS::dcq(
          reference_data = log2(1 + as.matrix(ref)),
          mix_data = log2(1 + as.matrix(mat)), ## log data OK??

          marker_set = cbind(intersect(rownames(ref), rownames(mat))),
          alpha_used = 0.05, lambda_min = 0.2, number_of_repeats = 3,
          precent_of_data = 1.0
        )
      )
    )
    if (!is.null(res.dcq) && !inherits(res.dcq, "try-error")) {
      timings[["DCQ"]] <- stime
      dbg("deconvolution using DCQ took", stime[3], "s\n")
      results[["DCQ"]] <- res.dcq$average
    } else {
      dbg("WARNING:: DCQ failed\n")
    }
  }

  if ("I-NNLS" %in% methods) {
    ## !!!!!!!!!!!!!! WARNING:: needs more testing/validation !!!!!!!!!!!!!
    ## ----- Constrained iterative non-negative least squares (Abbas et al. 2009) ----
    dbg("[pgx.deconvolution] calculating I-NNLS...")

    GetFractions.Abbas <- function(XX, y, w = NA) {
      ## XX is immune expression data
      ## y is cancer expression data
      ss.remove <- c()
      ss.names <- colnames(XX)
      while (TRUE) {
        if (length(ss.remove) == 0) {
          tmp.XX <- XX
        } else {
          if (is.null(ncol(tmp.XX))) {
            return(rep(0, ncol(XX)))
          }
          tmp.XX <- tmp.XX[, -ss.remove, drop = FALSE]
        }
        if (length(ss.remove) > 0) {
          ss.names <- ss.names[-ss.remove]
          if (length(ss.names) == 0) {
            return(rep(0, ncol(XX)))
          }
        }
        if (is.na(w[1])) {
          tmp <- stats::lsfit(tmp.XX, y, intercept = FALSE)
        } else {
          tmp <- stats::lsfit(tmp.XX, y, w, intercept = FALSE)
        }
        if (ncol(tmp.XX) == 1) {
          tmp.beta <- tmp$coefficients[1]
        } else {
          tmp.beta <- tmp$coefficients[1:(ncol(tmp.XX) + 0)]
        }
        if (min(tmp.beta > 0)) break ## break if coefs are all positive
        ss.remove <- which.min(tmp.beta) ## removes most negative coeff
      }
      tmp.F <- rep(0, ncol(XX))
      names(tmp.F) <- colnames(XX)
      tmp.F[ss.names] <- tmp.beta
      return(tmp.F)
    }
    gg <- intersect(rownames(ref), rownames(mat))
    XX <- ref[gg, , drop = FALSE]
    YY <- mat[gg, , drop = FALSE]
    res.abbas <- NULL
    stime <- system.time(
      res.abbas <- try(apply(YY, 2, function(y) GetFractions.Abbas(XX, y, w = NA)))
    )
    timings[["I-NNLS"]] <- stime
    dbg("deconvolution using I-NNLS took", stime[3], "s\n")
    if (!is.null(res.abbas) && !inherits(res.abbas, "try-error")) {
      rownames(res.abbas) <- colnames(ref)
      results[["I-NNLS"]] <- t(res.abbas)
    }
  }

  ## Own NNLM (non-negative linear modeling)...
  if ("NNLM" %in% methods) {
    ## NNLM

    dbg("[pgx.deconvolution] calculating NNLM...")

    x1 <- log2(1 + ref[gg, , drop = FALSE])
    x2 <- log2(1 + mat[gg, , drop = FALSE])
    x1 <- cbind(offset = 1, x1)
    stime <- system.time(
      cf <- NNLM::nnlm(x1, x2)$coefficients[-1, , drop = FALSE]
    )
    timings[["NNLM"]] <- stime
    dbg("deconvolution using NNLM took", stime[3], "s\n")
    results[["NNLM"]] <- t(cf)

    ## very much the same as I-NNLS
    results[["NNLM.lin"]] <- t(NNLM::nnlm(ref[gg, , drop = FALSE], mat[gg, , drop = FALSE])$coefficients)

    r1 <- apply(ref[gg, , drop = FALSE], 2, rank, na.last = "keep")
    r2 <- apply(mat[gg, , drop = FALSE], 2, rank, na.last = "keep")
    r1 <- cbind(offset = 1, r1)
    cf <- NNLM::nnlm(r1, r2)$coefficients[-1, , drop = FALSE]
    results[["NNLM.rnk"]] <- t(cf)
  }

  ## Simple (rank) correlation
  if ("cor" %in% methods) {
    dbg("[pgx.deconvolution] calculating cor...")

    r1 <- apply(mat[gg, , drop = FALSE], 2, rank, na.last = "keep")
    r2 <- apply(ref[gg, , drop = FALSE], 2, rank, na.last = "keep")
    stime <- system.time(
      cf <- stats::cor(r1, r2, use = "pairwise")
    )
    timings[["cor"]] <- stime
    dbg("deconvolution using COR took", stime[3], "s\n")
    results[["cor"]] <- cf
  }

  if ("SingleR" %in% methods) {
    dbg("[pgx.deconvolution] calculating SingleR...")
    stime <- system.time(
      sr1 <- SingleR::SingleR(test = mat, ref = ref, labels = colnames(ref))
    )
    timings[["SingleR"]] <- stime
    dbg("deconvolution using SingleR took", stime[3], "s\n")
    results[["SingleR"]] <- sr1$scores
  }
  ## clean up
  results <- results[which(!sapply(results, is.null))]
  results <- lapply(results, function(x) {
    x[is.na(x)] <- 0
    x
  })

  dbg("[pgx.deconvolution] calculating meta values...")

  ## meta
  if (length(results) > 1) {
    jj <- colnames(ref)
    norm.results <- lapply(results, function(x) x[, jj, drop = FALSE] / (1e-8 + rowSums(x[, jj, drop = FALSE])))
    lognorm.results <- lapply(norm.results, function(x) log(0.001 + pmax(x, 0)))
    res.meta1 <- Reduce("+", norm.results) / length(norm.results)
    res.meta2 <- exp(Reduce("+", lognorm.results) / length(lognorm.results))
    results[["meta"]] <- res.meta1 / (1e-8 + rowSums(res.meta1, na.rm = TRUE))
    results[["meta.prod"]] <- res.meta2 / (1e-8 + rowSums(res.meta2, na.rm = TRUE))
  }

  timings0 <- do.call(rbind, timings)
  res2 <- list(results = results, timings = timings0)

  dbg("[pgx.deconvolution] done!")

  return(res2)
}
