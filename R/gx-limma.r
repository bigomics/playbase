##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


########################################################################
## Compute LIMMA from matrices or from GCT/CLS file
##
########################################################################

#' Differential expression analysis using limma
#'
#' @title Differential expression analysis using limma
#'
#' @description Performs differential expression analysis on a gene expression matrix using limma.
#' Handles single sample case by duplicating sample. Auto-detects reference level.
#'
#' @param X Gene expression matrix with genes in rows and samples in columns.
#' @param pheno Phenotype vector or factor for samples.
#' @param B Optional batch covariate matrix.
#' @param remove.na Logical for removing samples with missing phenotype.
#' @param fdr FDR threshold for identifying differentially expressed genes.
#' @param compute.means Logical for computing group means.
#' @param lfc Log fold change cutoff.
#' @param max.na Max proportion of missing values allowed for a gene.
#' @param ref Character vector of possible reference levels.
#' @param trend Logical for fitting a trend model.
#' @param verbose Verbosity level.
#'
#' @details This function performs differential expression analysis on the gene expression matrix \code{X} using limma.
#' It handles filtering, model design matrices, and output formatting.
#' The phenotype vector \code{pheno} is used to define the linear model.
#' Batch effects can be accounted for by providing a \code{B} matrix.
#' The function auto-detects the reference level or the user can provide possible values in \code{ref}.
#' It returns a data frame containing the limma analysis results.
#'
#' @return Data frame with limma results.
#'
#' @export
gx.limma <- function(X, pheno, B = NULL, remove.na = TRUE,
                     fdr = 0.05, compute.means = TRUE, lfc = 0.20,
                     max.na = 0.20, ref = c(
                       "ctrl", "ctr", "control", "dmso", "nt", "0", "0h", "0hr",
                       "non", "no", "not", "neg", "negative", "ref", "veh", "vehicle",
                       "wt", "wildtype", "untreated", "normal", "false", "healthy"
                     ),
                     trend = FALSE, method=1, verbose = 1) {
  if (sum(duplicated(rownames(X))) > 0) {
    cat("WARNING:: matrix has duplicated rownames\n")
  }

  if (!is.null(B) && NCOL(B) == 1) {
    B <- matrix(B, ncol = 1)
    rownames(B) <- rownames(pheno)
    colnames(B) <- "batch"
  }

  ## detect single sample case
  is.single <- (max(table(pheno)) == 1)
  if (is.single) {
    cat("WARNING:: no replicates, duplicating samples...\n")
    X <- cbind(X, X)
    pheno <- c(pheno, pheno)
    if (!is.null(B)) B <- rbind(B, B)
  }

  ## filter probes and samples??
  ii <- which(rowMeans(is.na(X)) <= max.na)
  jj <- 1:ncol(X)
  if (remove.na && any(is.na(pheno))) {
    jj <- which(!is.na(pheno))
    if (verbose > 0) message(sum(is.na(pheno) > 0), "with missing phenotype\n")
  }
  X0 <- X[ii, jj, drop = FALSE]
  pheno0 <- as.character(pheno[jj])
  X0 <- X0[!(rownames(X0) %in% c(NA, "", "NA")), ]
  B0 <- NULL
  if (!is.null(B)) B0 <- B[jj, , drop = FALSE]

  if (verbose > 0) {
    cat("analyzing", ncol(X0), "samples\n")
    cat("table.pheno: ", table(pheno), "samples\n")
    cat("testing", nrow(X0), "features\n")
    if (!is.null(B0)) cat("including", ncol(B0), "batch covariates\n")
  }

  ## auto-detect reference
  pheno.ref <- c()
  ref.detected <- FALSE
  ref <- toupper(ref)

  is.ref <- (toupper(pheno0) %in% toupper(ref))
  ref.detected <- (sum(is.ref) > 0 && sum(!is.ref) > 0)
  ref.detected

  if (ref.detected) {
    pheno.ref <- unique(pheno0[which(toupper(pheno0) %in% toupper(ref))])
    if (verbose > 0) cat("setting reference to y=", pheno.ref, "\n")
    levels <- c(pheno.ref, sort(setdiff(unique(pheno0), pheno.ref)))
  } else {
    if (verbose > 0) cat("WARNING: could not auto-detect reference\n")
    levels <- as.character(sort(unique(pheno0)))
    if (verbose > 0) cat("setting reference to first class", levels[1], "\n")
  }
  if (length(levels) != 2) {
    stop("gx.limma::fatal error:only two class comparisons. Please use gx.limmaF().")
    return
  }

  ## setup model and perform LIMMA. See LIMMA userguide p41 ("Two groups").
  ## https://bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pdf
  if(method==1) {
    ## first method without  contrast matrix
    design <- cbind(1, pheno0 == levels[2])
    colnames(design) <- c("intercept", "main_vs_ref")
    if (!is.null(B0)) {
      if (verbose > 0) cat("augmenting design matrix with:", paste(colnames(B0)), "\n")
      sel <- which(colMeans(B0 == 1) < 1) ## take out any constant term
      design <- cbind(design, B0[, sel, drop = FALSE])
    }
    fit <- limma::lmFit(X0, design)
    fit <- limma::eBayes(fit, trend = trend)
    top <- limma::topTable(fit, coef = "main_vs_ref", number = nrow(X0), sort.by='none')
  } else {
    ## second possible method with explicit contrast matrix
    design <- cbind(1*(pheno0 == levels[1]), 1*(pheno0 == levels[2]))
    design <- as.matrix(design)
    colnames(design) <- c("ref", "main")    
    if (!is.null(B0)) {
      if (verbose > 0) cat("augmenting design matrix with:", paste(colnames(B0)), "\n")
      sel <- which(colMeans(B0 == 1) < 1) ## take out any constant term
      design <- cbind(design, B0[, sel, drop = FALSE])
    }
    fit <- limma::lmFit(X0, design)
    contr.matrix <- limma::makeContrasts( main_vs_ref = main-ref, levels=design)
    fit2 <- limma::contrasts.fit(fit, contr.matrix)
    fit2 <- limma::eBayes(fit2, trend = trend)
    top <- limma::topTable(fit2, coef = "main_vs_ref", number = Inf, sort.by='none')
  }

  ## give rownames
  if ("ID" %in% colnames(top)) {
    rownames(top) <- top$ID
    top$ID <- NULL
  }
  top <- top[rownames(X0), ]

  ## only significant
  top <- top[which(top$adj.P.Val <= fdr & abs(top$logFC) >= lfc), ]
  if (verbose > 0) cat("found", nrow(top), "significant at fdr=", fdr, "and minimal FC=", lfc, "\n")

  if (compute.means && nrow(top) > 0) {
    avg <- t(apply(
      X0[rownames(top), ], 1,
      function(x) tapply(x, pheno0, mean, na.rm = TRUE)
    ))
    avg <- avg[, as.character(levels), drop = FALSE]
    colnames(avg) <- paste0("AveExpr.", colnames(avg))
    top <- cbind(top, avg)
  }
  top$B <- NULL

  if (is.single) {
    top$P.Value <- NA
    top$adj.P.Val <- NA
    top$t <- NA
  }

  ## reorder on fold change
  top <- top[order(abs(top$logFC), decreasing = TRUE), ]


  ## unlist???

  return(top)
}


#' Differential expression analysis with limma
#'
#' @param X Numeric gene expression matrix with genes in rows and samples in columns.
#' @param pheno Data frame with phenotype data for samples. Must have column named 'group'.
#' @param B Data frame with batch data for samples. Default is NULL.
#' @param fdr FDR threshold for significance. Default is 0.05.
#' @param compute.means Logical indicating whether to compute group means. Default is TRUE.
#' @param lfc Log fold change threshold. Default is 0.2.
#' @param max.na Maximum missing value fraction for gene filtering. Default is 0.2.
#' @param ref Character vector of reference group names to use as baseline. Default is common control names.
#' @param trend Logical indicating whether to fit a trend model. Default is FALSE.
#' @param verbose Verbosity level. Default is 1.
#'
#' @return List with differential expression results, including:
#' \itemize{
#'   \item tab - Data frame with stats for all genes
#'   \item top - Data frame with stats for top significant genes
#'   \item fstats - Data frame with F statistics for all genes
#'   \item means - Data frame with mean expression by group
#' }
#'
#' @details This function performs differential expression analysis on \code{X} using limma.
#' It handles filtering, model design matrices, and output formatting.
#'
#' @examples
#' \dontrun{
#' # TODO
#' }
#' @export
gx.limmaF <- function(X, pheno, B = NULL, fdr = 0.05, compute.means = TRUE, lfc = 0.20,
                      max.na = 0.20, ref = c(
                        "ctrl", "ctr", "control", "dmso", "nt", "0", "0h", "0hr",
                        "non", "no", "not", "neg", "negative", "ref", "veh", "vehicle",
                        "wt", "wildtype", "untreated", "normal", "false", "healthy"
                        ),
                      trend = FALSE, method=1, verbose = 1) {
  if (sum(duplicated(rownames(X))) > 0) {
    cat("matrix has duplicated rownames. please remove.\n")
  }
  ## detect single sample case
  is.single <- (max(table(pheno)) == 1)
  if (is.single) {
    cat("warning: no replicates, no stats...\n")
    X <- cbind(X, X)
    pheno <- c(pheno, pheno)
  }
  if (!is.null(B) && NCOL(B) == 1) {
    B <- matrix(B, ncol = 1)
    rownames(B) <- rownames(pheno)
    colnames(B) <- "batch"
  }

  ## filter probes and samples
  ii <- which(rowMeans(is.na(X)) <= max.na)
  jj <- which(!is.na(pheno))
  if (verbose > 0) cat(sum(is.na(pheno) > 0), "with missing phenotype\n")
  X0 <- X[ii, jj]
  pheno0 <- as.character(pheno[jj])
  X0 <- X0[!(rownames(X0) %in% c(NA, "", "NA")), ]
  B0 <- NULL
  if (!is.null(B)) B0 <- B[jj, , drop = FALSE]

  if (verbose > 0) {
    cat("analyzing", ncol(X0), "samples\n")
    cat("testing", nrow(X0), "features\n")
    cat("in", length(unique(pheno0)), "groups\n")
    if (!is.null(B0)) cat("including", ncol(B0), "batch covariates\n")
  }

  ## auto-detect reference
  pheno.ref <- c()
  ref.detected <- FALSE
  ref <- toupper(ref)

  is.ref <- (toupper(pheno0) %in% toupper(ref))
  ref.detected <- (sum(is.ref) > 0 && sum(!is.ref) > 0)
  ref.detected

  if (ref.detected) {
    pheno.ref <- unique(pheno0[which(toupper(pheno0) %in% toupper(ref))])
    if (verbose > 0) cat("setting reference to y=", pheno.ref, "\n")
    levels <- c(pheno.ref, sort(setdiff(unique(pheno0), pheno.ref)))
    pheno1 <- stats::relevel(factor(pheno0), ref = levels[1])
  } else {
    if (verbose > 0) cat("WARNING: could not auto-detect reference\n")
    levels <- as.character(sort(unique(pheno0)))
    if (verbose > 0) cat("setting reference to first class", levels[1], "\n")
    pheno1 <- stats::relevel(factor(pheno0), ref = levels[1])
  }
  if (0 && length(levels) != 2) {
    stop("gx.limma::fatal error:only two class comparisons")
    return
  }

  if(method==1) {
    ## this tests implicitly against the reference level (REF), It
    ## uses no contrast matrix.
    design <- stats::model.matrix(~pheno1)
    colnames(design)
    colnames(design)[2:ncol(design)] <- paste0(levels(pheno1)[-1], "_vs_", levels(pheno1)[1])
    if (!is.null(B0)) {
      if (verbose > 0) cat("augmenting design matrix with:", paste(colnames(B0)), "\n")
      sel <- which(colMeans(B0 == 1) < 1) ## take out any constant term
      design <- cbind(design, B0[, sel, drop = FALSE])
    }
    fit <- limma::lmFit(X0, design)
    fit <- limma::eBayes(fit, trend = trend)
    top <- limma::topTable(fit, number = Inf, sort.by='none')
  }
  if(method==2) {
    ## this tests implicitly against the reference level (REF), It
    ## uses explicit contrast matrix.
    design <- stats::model.matrix(~0+pheno1)
    colnames(design)
    colnames(design) <- sub("^pheno1","",colnames(design))
    if (!is.null(B0)) {
      if (verbose > 0) cat("augmenting design matrix with:", paste(colnames(B0)), "\n")
      sel <- which(colMeans(B0 == 1) < 1) ## take out any constant term
      design <- cbind(design, B0[, sel, drop = FALSE])
    }
    fit <- limma::lmFit(X0, design)
    ct <- makeDirectContrasts( data.frame(pheno1,check.names=FALSE),
      ref=levels(pheno1)[1])   
    contrast.matrix <- ct$contr.matrix    
    grp2pheno.mat <- table(ct$group, pheno1)
    grp2pheno.mat <- grp2pheno.mat[rownames(contrast.matrix),]
    grp2pheno <- colnames(grp2pheno.mat)[max.col(grp2pheno.mat)]
    rownames(contrast.matrix) <- grp2pheno
    contrast.matrix <- contrast.matrix[colnames(design),]
    fit2 <- limma::contrasts.fit(fit, contrast.matrix)
    fit2 <- limma::eBayes(fit2, trend = trend)
    top  <- limma::topTable(fit2, number = Inf, sort.by='none')
  }
  if(method==3) {
    ## this tests implicitly all comparisons, It uses explicit full
    ## contrast matrix.
    design <- stats::model.matrix(~0+pheno1)
    colnames(design)
    colnames(design) <- sub("^pheno1","",colnames(design))
    if (!is.null(B0)) {
      if (verbose > 0) cat("augmenting design matrix with:", paste(colnames(B0)), "\n")
      sel <- which(colMeans(B0 == 1) < 1) ## take out any constant term
      design <- cbind(design, B0[, sel, drop = FALSE])
    }
    fit <- limma::lmFit(X0, design)
    contrast.matrix <- makeFullContrasts(pheno1)
    contrast.matrix <- contrast.matrix[colnames(design),]
    fit2 <- limma::contrasts.fit(fit, contrast.matrix)
    fit2 <- limma::eBayes(fit2, trend = trend)
    top <- limma::topTable(fit2, number = Inf, sort.by='none')
  }
  
  ## clean-up
  cols <- c("logFC","AveExpr","F","P.Value","adj.P.Val")
  top <- top[, intersect(colnames(top), cols)]
  top$B <- NULL
  if ("ID" %in% colnames(top)) {
    rownames(top) <- top$ID
    top$ID <- NULL
  }
  top <- top[, setdiff(colnames(top), colnames(design)), drop = FALSE]
  top <- top[rownames(X0), ]

  ## compute averages
  avg <- do.call(cbind, tapply(1:ncol(X0), pheno1, function(i) {
    rowMeans(X0[, i, drop = FALSE])
  }))

  if (!"logFC" %in% colnames(top)) {
    maxFC <- apply(avg, 1, max, na.rm = TRUE) - apply(avg, 1, min, na.rm = TRUE)
    top$logFC <- NULL
    top <- cbind(logFC = maxFC, top)
    rownames(top) <- rownames(X0)
  }

  ## only significant
  top <- top[which(top$adj.P.Val <= fdr & abs(top$logFC) >= lfc), ]
  if (verbose > 0) cat("found", nrow(top), "significant at fdr=", fdr, "and minimal FC=", lfc, "\n")

  if (compute.means && nrow(top) > 0) {
    avg1 <- avg[rownames(top), ]
    colnames(avg1) <- paste0("AveExpr.", colnames(avg1))
    top <- cbind(top, avg1)
  }
  top$B <- NULL

  if (is.single) {
    top$P.Value <- NA
    top$adj.P.Val <- NA
    top$t <- NA
    top$F <- NA
  }

  ## reorder on fold change
  top <- top[order(abs(top$logFC), decreasing = TRUE), ]


  ## unlist

  return(top)
}






#' Differential expression analysis for paired data
#'
#' @title Differential expression analysis using limma for paired samples
#'
#' @description Performs differential expression analysis using limma on paired or repeated measurement data.
#' Designed for pre-post, treated-untreated, or matched designs with two measurements per sample.
#'
#' @param X Numeric matrix of gene expression values (genes in rows, samples in columns).
#' @param pheno Phenotype factor vector indicating the sample classes.
#' @param pair Factor vector indicating which samples are paired.
#' @param fdr FDR threshold for identifying differentially expressed genes.
#' @param lfc Log2 fold change cutoff for differential expression.
#' @param ref Character vector of reference levels to use as baseline in contrasts.
#' @param compute.means Logical indicating whether to append mean expression values.
#' @param trend Logical indicating whether to fit a mean-variance trend.
#'
#' @details This function handles the limma model matrix and contrasts for paired designs.
#' The \code{pair} vector indicates which samples are paired measurements.
#' A differential expression analysis is performed between the phenotype classes.
#' Results are filtered by FDR and fold change thresholds.
#'
#' @return Data frame containing limma analysis results.
#'
#' @export
gx.limma.paired <- function(X, pheno, pair, fdr = 0.05, lfc = 0.20,
                            ref = c(
                              "ctrl", "ctr", "control", "dmso", "nt", "0", "0h", "0hr",
                              "non", "no", "not", "neg", "negative", "ref", "veh", "vehicle",
                              "wt", "wildtype", "untreated", "normal", "false", "healthy"
                            ),
                            compute.means = FALSE, trend = FALSE) {
  ## LIMMA

  cat("Paired LIMMA\n")
  cat("analyzing", ncol(X), "samples\n")
  X <- X[!(rownames(X) %in% c(NA, "", "NA")), ]
  pheno <- as.character(pheno)
  pair <- as.character(pair)

  ## check
  a <- as.character(pair)
  b <- as.character(pheno)
  v1 <- sort(unique(a))
  v2 <- sort(unique(b))
  if (!is.null(ref) && sum(v1 %in% ref) > 0) {
    r1 <- sort(intersect(v1, ref))
    cat("setting reference to", r1, "\n")
    v1 <- c(r1, setdiff(v1, r1))
  }
  if (!is.null(ref) && sum(v2 %in% ref) > 0) {
    r2 <- sort(intersect(v2, ref))
    cat("setting reference to", r2, "\n")
    v2 <- c(r2, setdiff(v2, r2))
  }
  a <- factor(a, levels = v1)
  b <- factor(b, levels = v2)
  v1 <- levels(a)
  v2 <- levels(b)
  cat("pairs:", paste(v1), "\n")
  cat("phenotypes:", paste(v2), "\n")
  if (length(v2) > 2) {
    cat("gx.limma.paired:: fatal error: only two-groups implemented\n")
    return(NULL)
  }
  if (length(v1) < 2) {
    cat("gx.limma.paired:: fatal error: no pairs\n")
    return(NULL)
  }

  ## setup LIMMA (paired t-test)

  design <- stats::model.matrix(~ a + b)

  ## perform fitting
  fit0 <- limma::lmFit(X, design)
  fit2 <- limma::eBayes(fit0, trend = trend)

  ## extract toptable
  bcoef <- grep("^b", colnames(fit2$coefficients))
  top <- limma::topTable(fit2, coef = bcoef, number = nrow(X))
  if (colnames(top)[1] == "ID") { ## old style
    rownames(top) <- top[, "ID"]
    top$ID <- NULL
  }
  top <- top[rownames(X), ]

  ## only significant
  if (fdr < 1) {
    kk <- which(top$adj.P.Val <= fdr & abs(top$logFC) >= lfc)
    top <- top[kk, ]
  }
  cat("found", nrow(top), "significant at fdr=", fdr, "and minimum logFC=", lfc, "\n")

  ## compute means if requested
  X.m <- NULL
  if (compute.means && nrow(top) > 0) {
    ff <- paste(pair, pheno, sep = ".")
    X.m <- t(apply(
      X[rownames(top), ], 1,
      function(x) tapply(x, ff, mean, na.rm = TRUE)
    ))
    top$AveExpr <- NULL
    top$AveExpr <- X.m
  } else {
    ff <- pheno
    X.m <- t(apply(
      X[rownames(top), ], 1,
      function(x) tapply(x, ff, mean, na.rm = TRUE)
    ))
    top$AveExpr <- NULL
    top$AveExpr <- X.m
  }

  ## fold-change
  jj <- order(-abs(top$logFC))
  top <- top[jj, ]
  if (!is.null(X.m)) X.m <- X.m[rownames(top), ]

  ## results
  top$B <- NULL
  return(top)
}

## two-factorial design







#' Differential expression analysis using limma for RNA-seq data
#'
#' @title Differential expression analysis with limma for RNA-seq data
#'
#' @description Performs differential expression analysis on RNA-seq count data
#'  using limma-voom. Designed for gene-level count data from RNA-seq
#' experiments.
#'
#' @param countdata Numeric matrix of RNA-seq counts with genes in rows and
#' samples in columns.
#' @param y Factor or vector indicating the phenotype for each sample.
#' @param method Method for normalization and dispersion estimation. Either
#' "edgeR" or "DESeq2".
#'
#' @details This function takes RNA-seq count data and a phenotype factor as input.
#'  It uses limma-voom to fit a linear model and identify differentially expressed
#' genes. The count data is normalized and dispersion is estimated using either edgeR or DESeq2 methods.
#'
#' @return Data frame containing limma analysis results.
#'
#' @examples
#' \dontrun{
#' set.seed(151)
#' gx <- matrix(sample(1:100, 100 * 20, replace = TRUE), 100, 20)
#' rownames(gx) <- replicate(100, sample(paste0(LETTERS, 1:50), 1))
#' colnames(gx) <- sample(letters, 20)
#' y <- sample(letters[1:4], 20, replace = TRUE)
#' out <- playbase::seq_limma(gx, y)
#' #'
#' }
#' @export
seq_limma <- function(countdata, y, method = "edgeR") {
  ## https://bioinformatics-core-shared-training.github.io/RNAseq-R/rna-seq-de.nb.html
  if (min(countdata) < 0 || !all(countdata %% 1 == 0)) {
    cat("WARNING:: input X should be integer counts! Proceed on own risk\n")
  }

  ## Identify genes with at least 0.5 CPM in at least 2 samples
  myCPM <- edgeR::cpm(countdata)
  thresh <- myCPM > 0.5
  keep <- rowSums(thresh) >= 2
  ## Subset the rows of countdata to keep the more highly expressed genes
  counts.keep <- countdata[keep, ]

  ## Convert to an edgeR object
  group <- factor(y)
  dgeObj <- edgeR::DGEList(counts.keep, group = group)

  ## Perform TMM normalisation
  dgeObj <- edgeR::calcNormFactors(dgeObj, method = "TMM")

  ## Define design matrix
  design <- stats::model.matrix(~group)

  ## Estimating the dispersion

  dgeObj <- edgeR::estimateGLMCommonDisp(dgeObj, design)
  dgeObj <- edgeR::estimateGLMTrendedDisp(dgeObj)
  dgeObj <- edgeR::estimateTagwiseDisp(dgeObj)

  ## Fit the linear model
  fit <- edgeR::glmFit(dgeObj, design)
  names(fit)

  ## Conduct likelihood ratio tests
  res.test <- edgeR::glmLRT(fit, coef = 2)
  toptable <- edgeR::topTags(res.test, n = nrow(countdata))@.Data[[1]]
  Matrix::head(toptable)



  ## calculate means
  logcpm <- edgeR::cpm(dgeObj, log = TRUE)
  xmean <- c()
  for (y0 in unique(y)) {
    m1 <- rowMeans(logcpm[, which(y == y0)], na.rm = TRUE)
    xmean <- cbind(xmean, m1)
  }
  colnames(xmean) <- paste0("mean.", unique(y))
  Matrix::head(xmean)
  xmean <- cbind(mean = rowMeans(xmean), xmean)
}
