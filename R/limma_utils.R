# gx.limma
pgx_limma <- function(X, pheno, B = NULL, remove.na = TRUE,
                      fdr = 0.05, compute.means = TRUE, lfc = 0.20,
                      max.na = 0.20, ref = REF.CLASS, trend = FALSE, verbose = 1) {
  if (0) {
    fdr <- 0.05
    compute.means <- TRUE
    lfc <- 0.20
    ref <- REF.CLASS
    max.na <- 0.2
    trend <- FALSE
    verbose <- 1
    B <- NULL
  }
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
  ## is.ref <- grepl(paste(ref,collapse="|"),pheno0)
  is.ref <- (toupper(pheno0) %in% toupper(ref))
  ref.detected <- (sum(is.ref) > 0 && sum(!is.ref) > 0)
  ref.detected

  ## if(!is.null(ref) && sum( toupper(pheno0) %in% ref)>0 ) {
  if (ref.detected) {
    pheno.ref <- unique(pheno0[which(toupper(pheno0) %in% toupper(ref))])
    if (verbose > 0) cat("setting reference to y=", pheno.ref, "\n")
    bb <- c(pheno.ref, sort(setdiff(unique(pheno0), pheno.ref)))
  } else {
    if (verbose > 0) cat("WARNING: could not auto-detect reference\n")
    bb <- as.character(sort(unique(pheno0)))
    if (verbose > 0) cat("setting reference to first class", bb[1], "\n")
  }
  if (length(bb) != 2) {
    stop("gx.limma::fatal error:only two class comparisons. Please use gx.limmaF().")
    return
  }

  ## setup model and perform LIMMA
  design <- cbind(1, pheno0 == bb[2])
  colnames(design) <- c("WT", "2vs1")
  d1 <- colnames(design)[1]
  d2 <- colnames(design)[2]

  if (!is.null(B0)) {
    if (verbose > 0) cat("augmenting design matrix with:", paste(colnames(B0)), "\n")
    sel <- which(colMeans(B0 == 1) < 1) ## take out any constant term
    design <- cbind(design, B0[, sel, drop = FALSE])
  }

  fit <- limma::lmFit(X0, design)
  fit <- limma::eBayes(fit, trend = trend)
  top <- limma::topTable(fit, coef = d2, number = nrow(X0))
  if ("ID" %in% colnames(top)) {
    rownames(top) <- top$ID
    top$ID <- NULL
  }
  top <- top[rownames(X0), ]
  Matrix::head(top)

  ## only significant
  top <- top[which(top$adj.P.Val <= fdr & abs(top$logFC) >= lfc), ]
  if (verbose > 0) cat("found", nrow(top), "significant at fdr=", fdr, "and minimal FC=", lfc, "\n")

  if (compute.means && nrow(top) > 0) {
    avg <- t(apply(
      X0[rownames(top), ], 1,
      function(x) tapply(x, pheno0, mean, na.rm = TRUE)
    ))
    avg <- avg[, as.character(bb), drop = FALSE]
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
  ## colnames(top) <-   sub("logFC","logR",colnames(top))

  ## unlist???
  ## top = do.call(cbind, top)
  return(top)
}
