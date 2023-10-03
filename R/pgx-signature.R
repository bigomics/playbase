##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

## ================================================================================
## ========================= CONNECTIVITY FUNCTIONS ===============================
## ================================================================================


#' @title Compute connectivity scores
#'
#' @param pgx PGX object containing gene expression data
#' @param sigdb Signature database HDF5 file
#' @param ntop Number of top correlated genes to use for scoring. Default 1000.
#' @param contrasts Contrasts to compute signature for. Default is all contrasts.
#' @param remove.le Remove leading edge genes before scoring. Default FALSE.
#' @param inmemory Compute correlations in memory instead of on disk. Default FALSE.
#'
#' @return List of data frames with connectivity scores for each contrast.
#'
#' @description Compute connectivity scores between a signature and datasets.
#'
#' @details This function takes a signature from the PGX object (fold changes for a given contrast)
#' and computes connectivity scores against the signature database. The signature is correlated
#' against all datasets in the database to produce a ranked list of connectivity scores.
#'
#' The top \code{ntop} correlated genes are used for computing the connectivity score.
#' Scores are computed for each contrast unless a subset is specified with \code{contrasts}.
#'
#' Setting \code{remove.le = TRUE} removes the leading edge genes from the signature before
#' computing connectivity. This reduces bias from highly weighted genes.
#'
#' By default, correlations are computed on disk using chunks. Set \code{inmemory = TRUE}
#' to compute in memory, which may be faster for small datasets.
#'
#' @export
pgx.computeConnectivityScores <- function(pgx, sigdb, ntop = 1000, contrasts = NULL,
                                          remove.le = FALSE, inmemory = FALSE) {
  meta <- pgx.getMetaFoldChangeMatrix(pgx, what = "meta")

  is.h5ref <- grepl("h5$", sigdb)
  if (!is.h5ref) {
    cat("[pgx.computeConnectivityScores] ERROR: must be H5 formatted file\n")
    return(NULL)
  }
  if (inmemory) {
    cat("[pgx.computeConnectivityScores] *** using in-memory ***\n")
  }

  h5.file <- NULL
  if (file.exists(sigdb)) {
    h5.file <- sigdb
  }
  if (is.null(h5.file)) {
    cat("[pgx.computeConnectivityScores] ERROR: could not open H5 file\n")
    return(NULL)
  }

  if (is.null(contrasts)) {
    contrasts <- colnames(meta$fc)
  }
  contrasts <- intersect(contrasts, colnames(meta$fc))
  F1 <- meta$fc[, contrasts, drop = FALSE]

  if (inmemory) {
    scores <- pgx.correlateSignatureH5.inmemory(
      F1,
      h5.file = h5.file,
      nsig = 100,
      ntop = ntop,
      nperm = 9999
    )
  } else {
    scores <- list()
    ct <- colnames(F1)[1]
    for (ct in colnames(F1)) {
      fc <- F1[, ct]
      names(fc) <- rownames(meta$fc)
      names(fc) <- toupper(names(fc)) ## for MOUSE!!
      res <- pgx.correlateSignatureH5(
        fc,
        h5.file = h5.file,
        nsig = 100,
        ntop = ntop,
        nperm = 9999
      )
      scores[[ct]] <- res
    }
  }
  if (is.null(names(scores))) names(scores) <- contrasts

  ## remove leadingEdge (take too much memory!!!)
  if (remove.le) {
    for (j in 1:length(scores)) scores[[j]]$leadingEdge <- NULL
  }

  return(scores)
}


#' @title Correlate signature with datasets in HDF5 file
#'
#' @param F Numeric vector of signature fold changes
#' @param h5.file Path to HDF5 file containing dataset expression matrices
#' @param nsig Number of signature genes to use for correlation. Default 100.
#' @param ntop Number of top correlated genes to use for GSEA. Default 1000.
#' @param nperm Number of permutations for GSEA. Default 1000.
#'
#' @return Data frame with correlation results and GSEA enrichment scores.
#'
#' @description Correlates a signature (fold change vector) with datasets in an HDF5 file using in-memory computations instead of on-disk.
#' It is also possible to run the operation on-disk via \code{pgx.correlateSignatureH5}.
#'
#' @details This function computes the correlation and gene set enrichment between a
#' signature fold change vector \code{F} and dataset expression matrices stored in an
#' HDF5 file at \code{h5.file}.
#'
#' It first reads the top \code{nsig} signature genes into memory.
#' Then it computes pairwise Pearson correlation between the signature and each dataset.
#' It also runs pre-ranked GSEA using the \code{fgsea} package to compute enrichment scores.
#'
#' The main difference from \code{\link{pgx.correlateSignatureH5}} is that correlations
#' and GSEA are computed in memory instead of using on-disk chunks. This is faster for
#' small datasets but requires loading the full expression matrix into memory.
#'
#' @return Data frame with the correlation coefficient, p-value, and normalized enrichment score (NES) for each dataset.
#'
#' @export
pgx.correlateSignatureH5.inmemory <- function(F, h5.file, nsig = 100, ntop = 1000, nperm = 1000) {
  if (NCOL(F) == 1 && inherits(F, "numeric")) {
    rn <- names(F)
    F <- matrix(F, ncol = 1)
    rownames(F) <- rn
  }

  if (is.null(rownames(F))) stop("F must have rownames")
  ## mouse... mouse...
  rownames(F) <- toupper(rownames(F))

  ## or instead compute correlation on top100 fc genes (read from file)
  rn <- rhdf5::h5read(h5.file, "data/rownames")
  cn <- rhdf5::h5read(h5.file, "data/colnames")

  ## Entire matrix in memory????
  matG <- rhdf5::h5read(h5.file, "data/matrix") ### whole matrix!!!!
  matG[which(matG < -999999)] <- NA
  rownames(matG) <- rn
  colnames(matG) <- cn

  mem1 <- round(utils::object.size(matG) / 1e9, 2)
  cat("[pgx.correlateSignatureH5] utils::object.size(matG)=", mem1, "Gb\n") ## gigabytes....

  ## ---------------------------------------------------------------
  ## Compute simple correlation between query profile and signatures
  ## ---------------------------------------------------------------
  res <- list()
  i <- 1
  for (i in 1:ncol(F)) {
    gg <- intersect(rownames(F), rn)
    fc1 <- sort(F[gg, i])
    gg <- unique(names(c(Matrix::head(fc1, nsig), Matrix::tail(fc1, nsig))))

    remove(fc1)
    row.idx <- match(gg, rn)

    rG <- matG[row.idx, , drop = FALSE]
    rG <- apply(rG, 2, rank, na.last = "keep")
    dimnames(rG) <- list(rn[row.idx], cn)

    ## this FC signature
    fc <- F[, i]
    rfc <- rank(fc[gg], na.last = "keep") ## rank correlation??
    rG[is.na(rG)] <- 0
    rfc[is.na(rfc)] <- 0
    rho <- stats::cor(rG, rfc, use = "pairwise")[, 1]

    remove(rG, rfc)

    ## --------------------------------------------------
    ## test tops signatures using fGSEA
    ## --------------------------------------------------

    sel <- Matrix::head(names(sort(-abs(rho))), ntop)
    sel.idx <- match(sel, cn)
    sig100.up <- rhdf5::h5read(h5.file, "signature/sig100.up",
      index = list(NULL, sel.idx)
    )
    sig100.dn <- rhdf5::h5read(h5.file, "signature/sig100.dn",
      index = list(NULL, sel.idx)
    )


    ## combine up/down into one (unsigned GSEA test)
    gmt <- rbind(sig100.up, sig100.dn)
    gmt <- unlist(apply(gmt, 2, list), recursive = FALSE)
    names(gmt) <- cn[sel.idx]

    ## use entire fc vector
    suppressMessages(suppressWarnings(
      res1 <- fgsea::fgseaSimple(gmt, abs(fc), nperm = nperm) ## really unsigned???
    ))

    ## ---------------------------------------------------------------
    ## Combine correlation+GSEA by combined score (NES*rho)
    ## ---------------------------------------------------------------
    jj <- match(res1$pathway, names(rho))
    res1$rho <- rho[jj]
    res1$R2 <- rho[jj]**2
    res1$score <- (res1$R2 * res1$NES)

    fn <- colnames(F)[i]
    res[[fn]] <- res1[order(res1$score, decreasing = TRUE), , drop = FALSE]

    gc()
  }
  remove(matG)
  gc()

  return(res)
}


#' @describeIn  pgx.correlateSignatureH5.inmemory computes correlation and gene set enrichment between a
#' signature and datasets in an HDF5 file using on-disk chunked computations
#' @param h5.file: HDF5 file of reference expression signatures
#' @param ntop:    number of top signatures (in abs(rho)) to report
#' @param nperm:   number of permuations for fGSEA
#'
#' @export
pgx.correlateSignatureH5 <- function(fc, h5.file, nsig = 100, ntop = 1000, nperm = 10000) {
  if (is.null(names(fc))) stop("fc must have names")
  ## mouse... mouse...
  names(fc) <- toupper(names(fc))

  ## or instead compute correlation on top100 fc genes (read from file)
  rhdf5::h5closeAll()
  rn <- rhdf5::h5read(h5.file, "data/rownames")
  cn <- rhdf5::h5read(h5.file, "data/colnames")

  ## ---------------------------------------------------------------
  ## Compute simple correlation between query profile and signatures
  ## ---------------------------------------------------------------
  gg <- intersect(names(fc), rn)
  fc <- fc[gg]
  fc1 <- sort(fc[gg])
  gg <- unique(names(c(Matrix::head(fc1, nsig), Matrix::tail(fc1, nsig))))
  row.idx <- match(gg, rn)
  rhdf5::h5closeAll()
  G <- rhdf5::h5read(h5.file, "data/matrix", index = list(row.idx, 1:length(cn)))
  G[which(G < -999999)] <- NA
  dimnames(G) <- list(rn[row.idx], cn)

  ## rank correlation??
  rG  <- apply(G[gg, ], 2, rank, na.last = "keep")
  rfc <- rank(fc[gg], na.last = "keep")
  rG[is.na(rG)] <- 0
  rfc[is.na(rfc)] <- 0
  suppressWarnings(rho <- stats::cor(rG, rfc, use = "pairwise")[, 1])

  remove(G, rG, rfc)
  
  ## --------------------------------------------------
  ## test tops signatures using fGSEA
  ## --------------------------------------------------

  sel <- Matrix::head(names(sort(-abs(rho))), ntop)
  sel.idx <- match(sel, cn)

  ## if we have less than 100 genes, we should make smaller GMT sets!
  nsig <- min(100, round(length(fc) / 5))
  idx <- list(1:nsig, sel.idx)
  sig100.up <- rhdf5::h5read(h5.file, "signature/sig100.up", index = idx)
  sig100.dn <- rhdf5::h5read(h5.file, "signature/sig100.dn", index = idx)
  colnames(sig100.up) <- cn[sel.idx]
  colnames(sig100.dn) <- cn[sel.idx]
  
  ## combine up/down into one (unsigned GSEA test)
  gmt <- rbind(sig100.up, sig100.dn)
  gmt <- unlist(apply(gmt, 2, list), recursive = FALSE)
  names(gmt) <- cn[sel.idx]

  suppressMessages(suppressWarnings(
    res <- fgsea::fgseaSimple(gmt, abs(fc), nperm = nperm, scoreType = "pos")
  )) ## really unsigned???

  ## --------------------------------------------------
  ## Fisher test
  ## --------------------------------------------------
  fc.up <- fc[fc>0]
  fc.dn <- fc[fc<0]
  top.up <- head(names(sort(-fc.up)), 3*nsig)  ## RETHINK!
  top.dn <- head(names(sort(+fc.dn)), 3*nsig)  
  top.fc <- c(top.up, top.dn)
  bg <- intersect(names(fc), rn)  
  stats <- playbase::gset.fisher(top.fc, gmt, background=bg, fdr=1, min.genes=0, nmin=0)
  or.max <- max(stats$odd.ratio[!is.infinite(stats$odd.ratio)])
  stats$odd.ratio[is.infinite(stats$odd.ratio)] <- max(99, 2*or.max)

  ## --------------------------------------------------
  ## Sparse GSET matrix
  ## --------------------------------------------------
  bg <- intersect(names(fc), rn)
  gmt100.up <- unlist(apply(sig100.up, 2, list), recursive = FALSE)
  gmt100.dn <- unlist(apply(sig100.dn, 2, list), recursive = FALSE)    
  G1 <- playbase::gmt2mat(gmt100.up, bg=bg)
  G2 <- playbase::gmt2mat(gmt100.dn, bg=bg)
  G1 <- G1[match(bg,rownames(G1)),]
  G2 <- G2[match(bg,rownames(G2)),colnames(G1)] 
  G <- G1 - G2
  dim(G)
  remove(G1,G2)
##  tau <- qlcMatrix::corSparse( G[bg,], cbind(fc[bg]) )[,1]
  tau <- qlcMatrix::cosSparse( G[bg,], cbind(fc[bg]) )[,1]
  jj <- match(res$pathway, colnames(G))
  res$tau <- tau[jj]
  
  ## ---------------------------------------------------------------
  ## Combine correlation+GSEA by combined score (NES*rho)
  ## ---------------------------------------------------------------
  jj <- match(res$pathway, names(rho))
  res$rho <- rho[jj]
  res$R2 <- rho[jj]**2

  ii <- match(res$pathway, rownames(stats))
  res$odd.ratio <- stats$odd.ratio[ii]
  res$overlap   <- stats$overlap[ii]
  
  res$score <- res$R2 * res$NES * res$odd.ratio * abs(res$tau)
  res <- res[order(abs(res$score), decreasing = TRUE), ]

  return(res)
}


#' @title Create CREEDS signature database
#'
#' @param gmt.files Character vector of paths to CREEDS GMT files containing gene signatures.
#' @param h5.file Path to HDF5 file to write signature database to.
#' @param update.only Logical indicating whether to only update existing signatures. Default is FALSE.
#'
#' @return NULL. The CREEDS signature database is written to the specified HDF5 file.
#'
#' @description Reads gene signatures from CREEDS GMT files and writes to an HDF5 database.
#'
#' @details This function takes a set of CREEDS GMT files containing up/down regulated gene signatures.
#' It extracts the signatures and writes them to an HDF5 file in a standard format that can be used for
#' gene set enrichment analysis.
#'
#' Each GMT file contains multiple up and down regulated signatures for a perturbation experiment.
#' The signatures are named using the perturbation and direction, e.g. "DrugX-up".
#'
#' If update.only=FALSE, any existing signatures in the HDF5 file will be overwritten.
#' If TRUE, only signatures not already present will be added.
#'
#' @export
pgx.createCreedsSigDB <- function(gmt.files, h5.file, update.only = FALSE) {
  h5exists <- function(h5.file, obj) {
    xobjs <- apply(rhdf5::h5ls(h5.file)[, 1:2], 1, paste, collapse = "/")
    obj %in% gsub("^/|^//", "", xobjs)
  }

  if (update.only && h5exists(h5.file, "data/matrix")) {
    X <- rhdf5::h5read(h5.file, "data/matrix")
    rn <- rhdf5::h5read(h5.file, "data/rownames")
    cn <- rhdf5::h5read(h5.file, "data/colnames")
    rownames(X) <- rn
    colnames(X) <- cn
  } else {
    ## --------------------------------------------------
    ## make big FC signature matrix
    ## --------------------------------------------------
    F <- list()
    sig100.dn <- list()
    sig100.up <- list()
    cat("reading gene lists from", length(gmt.files), "gmt files ")
    i <- 1
    for (i in 1:length(gmt.files)) {
      if (!file.exists(gmt.files[i])) next()
      cat(".")
      try.error <- try(gmt <- read.gmt(gmt.files[i], add.source = TRUE))
      if (inherits(try.error, "try-error")) next()

      j1 <- grep("-up ", names(gmt))
      j2 <- grep("-dn ", names(gmt))
      f1 <- lapply(gmt[j1], function(gg) {
        x <- length(gg):1
        names(x) <- gg
        x
      })
      f2 <- lapply(gmt[j2], function(gg) {
        x <- -length(gg):-1
        names(x) <- gg
        x
      })

      s1 <- gmt[j1]
      s2 <- gmt[j2]

      ff <- lapply(1:length(f1), function(i) c(f1[[i]], f2[[i]]))
      sig.names <- sub("-up", "", names(f1))
      prefix <- gsub(".*/|single_|_perturbations|.gmt|_signatures", "", gmt.files[i])
      sig.names <- paste0("[CREEDS:", prefix, "] ", sig.names)

      names(s1) <- names(s2) <- names(ff) <- sig.names
      sig100.up <- c(sig100.up, s1)
      sig100.dn <- c(sig100.dn, s2)
      F <- c(F, ff)
    }
    cat("\n")

    genes <- as.vector(unlist(sapply(F[], names)))
    genes <- sort(unique(toupper(genes)))

    ## Filter out genes (not on known chromosomes...)
    gannot <- ngs.getGeneAnnotation(genes)
    sel <- which(!is.na(gannot$chr))
    genes <- sort(genes[sel])

    X <- lapply(F, function(x) x[match(genes, names(x))])
    X <- do.call(cbind, X)
    rownames(X) <- genes
    remove(F)

    h5.file
    pgx.saveMatrixH5(X, h5.file, chunk = c(nrow(X), 1))

    na100 <- rep(NA, 100)
    msig100.up <- sapply(sig100.up, function(g) Matrix::head(c(intersect(g, genes), na100), 100))
    msig100.dn <- sapply(sig100.dn, function(g) Matrix::head(c(intersect(g, genes), na100), 100))

    if (!h5exists(h5.file, "signature")) rhdf5::h5createGroup(h5.file, "signature")
    rhdf5::h5write(msig100.up, h5.file, "signature/sig100.up") ## can write list??
    rhdf5::h5write(msig100.dn, h5.file, "signature/sig100.dn") ## can write list???
    remove(sig100.up, sig100.dn, msig100.up, msig100.dn)

    ## check NA!!! sometimes it is set to large negative
    rhdf5::h5ls(h5.file)
    X <- rhdf5::h5read(h5.file, "data/matrix")
    X[which(X < -999999)] <- NA
    rhdf5::h5write(X, h5.file, "data/matrix") ## can write list??
    rhdf5::h5closeAll()
  }

  ## --------------------------------------------------
  ## Precalculate t-SNE/UMAP
  ## --------------------------------------------------

  if (!update.only || !h5exists(h5.file, "clustering")) {
    X[is.na(X)] <- 0
    pos <- pgx.clusterBigMatrix(
      abs(X), ## on absolute foldchange!!
      methods = c("pca", "tsne", "umap"),
      dims = c(2, 3),
      reduce.sd = 2000,
      reduce.pca = min(200, round(ncol(X) / 3))
    )

    if (!h5exists(h5.file, "clustering")) rhdf5::h5createGroup(h5.file, "clustering")
    rhdf5::h5ls(h5.file)
    rhdf5::h5write(pos[["pca2d"]], h5.file, "clustering/pca2d")
    rhdf5::h5write(pos[["pca3d"]], h5.file, "clustering/pca3d")
    rhdf5::h5write(pos[["tsne2d"]], h5.file, "clustering/tsne2d")
    rhdf5::h5write(pos[["tsne3d"]], h5.file, "clustering/tsne3d")
    rhdf5::h5write(pos[["umap2d"]], h5.file, "clustering/umap2d")
    rhdf5::h5write(pos[["umap3d"]], h5.file, "clustering/umap3d")
  }

  rhdf5::h5closeAll()
} ## end of pgx.createCreedsSigDB


#' Create a signature database from PGX files
#'
#' @title Create signature database
#' @description Reads fold change matrices from multiple PGX files and combines them into a single signature database stored in an H5 file.
#'
#' @param h5.file Path to the H5 file to store the signature database.
#' @param pgx.files Character vector of paths to PGX files to read fold change matrices from.
#' @param update.only Logical indicating if only new PGX files should be added. Default is FALSE.
#'
#' @details This function reads fold change matrices from multiple PGX experiment result files (.pgx).
#' It extracts the fold change values, normalizes the gene names, and combines the matrices into a single signature database.
#' This combined matrix is then stored in an H5 file for later use.
#'
#' @return The combined signature database matrix is saved to the specified H5 file.
#'
#' @export
pgx.createSignatureDatabaseH5 <- function(h5.file, pgx.files, update.only = FALSE) {
  h5exists <- function(h5.file, obj) {
    xobjs <- apply(rhdf5::h5ls(h5.file)[, 1:2], 1, paste, collapse = "/")
    obj %in% gsub("^/|^//", "", xobjs)
  }

  if (update.only && h5exists(h5.file, "data/matrix")) {
    X <- rhdf5::h5read(h5.file, "data/matrix")
    rn <- rhdf5::h5read(h5.file, "data/rownames")
    cn <- rhdf5::h5read(h5.file, "data/colnames")
    rownames(X) <- rn
    colnames(X) <- cn
  } else {
    ## --------------------------------------------------
    ## make big FC signature matrix
    ## --------------------------------------------------
    F <- list()
    cat("[pgx.createSignatureDatabaseH5] reading FC from", length(pgx.files), "pgx files ")
    i <- 1
    for (i in 1:length(pgx.files)) {
      if (!file.exists(pgx.files[i])) next()
      cat(".")
      pgx <- try(local(get(load(pgx.files[i], verbose = 0)))) ## override any name
      if (inherits(pgx, "try-error")) next()
      meta <- pgx.getMetaFoldChangeMatrix(pgx, what = "meta")
      rownames(meta$fc) <- toupper(rownames(meta$fc)) ## mouse-friendly
      pgxfile <- gsub(".*[/]|[.]pgx$", "", pgx.files[i])
      colnames(meta$fc) <- paste0("[", pgxfile, "] ", colnames(meta$fc))
      F[[pgxfile]] <- meta$fc
    }
    cat("\n")

    genes <- as.vector(unlist(sapply(F, rownames)))
    genes <- sort(unique(toupper(genes)))
    F <- lapply(F, function(x) x[match(genes, rownames(x)), , drop = FALSE])
    X <- do.call(cbind, F)
    rownames(X) <- genes

    ## Filter out genes (not on known chromosomes...)
    genes <- rownames(X)
    gannot <- ngs.getGeneAnnotation(genes)
    sel <- which(!is.na(gannot$chr))
    X <- X[sel, , drop = FALSE]
    remove(F)
  }

  pgx.createSignatureDatabaseH5.fromMatrix(
    X = X,
    h5.file = h5.file,
    update.only = update.only
  )
}


#' Create signature database from matrix
#'
#' @title Create signature database from matrix
#'
#' @param h5.file Path to HDF5 file to write signature database to
#' @param X Numeric matrix of gene expression data
#' @param update.only Logical indicating whether to update existing signatures only. Default is FALSE.
#'
#' @return NULL. The signature database is written to the HDF5 file.
#'
#' @description Reads a gene expression matrix and saves it to an HDF5 database file. Also extracts per-dataset top up/down regulated signatures.
#'
#' @details This function takes a gene expression matrix \code{X} and writes it to an HDF5 file at \code{h5.file} in a standard format.
#'
#' It also calculates the top 100 up and down regulated genes (signatures) for each dataset in \code{X}. These signatures are saved in the HDF5 file.
#'
#' If \code{update.only=TRUE}, existing signatures in the HDF5 file will not be overwritten. Only new signatures from \code{X} will be added.
#'
#' @export
pgx.createSignatureDatabaseH5.fromMatrix <- function(h5.file, X, update.only = FALSE) {
  if (file.exists(h5.file)) unlink(h5.file)
  dbg("[pgx.createSignatureDatabaseH5.fromMatrix] saving data matrix...")
  X <- as.matrix(X)
  pgx.saveMatrixH5(X, h5.file, chunk = c(nrow(X), 1))

  ## --------------------------------------------------
  ## Calculate top100 gene signatures
  ## --------------------------------------------------

  if (!update.only || !h5exists(h5.file, "signature")) {
    dbg("[pgx.createSignatureDatabaseH5.fromMatrix] Creating top-100 signatures...")
    rn <- rhdf5::h5read(h5.file, "data/rownames")
    cn <- rhdf5::h5read(h5.file, "data/colnames")
    rhdf5::h5ls(h5.file)
    X[is.na(X)] <- 0
    orderx <- apply(X, 2, function(x) {
      idx <- order(x)
      list(
        DN = utils::head(idx, 100),
        UP = rev(utils::tail(idx, 100))
      )
    })
    sig100.dn <- sapply(orderx, "[[", "DN")
    sig100.dn <- apply(sig100.dn, 2, function(i) rn[i])
    sig100.up <- sapply(orderx, "[[", "UP")
    sig100.up <- apply(sig100.up, 2, function(i) rn[i])

    dbg("[pgx.createSignatureDatabaseH5.fromMatrix] saving top-100 signatures...")
    if (!h5exists(h5.file, "signature")) rhdf5::h5createGroup(h5.file, "signature")
    rhdf5::h5write(sig100.dn, h5.file, "signature/sig100.dn") ## can write list???
    rhdf5::h5write(sig100.up, h5.file, "signature/sig100.up") ## can write list??

    remove(orderx)
    remove(sig100.dn)
    remove(sig100.up)
  }

  ## --------------------------------------------------
  ## Precalculate t-SNE/UMAP
  ## --------------------------------------------------
  if (!update.only || !h5exists(h5.file, "clustering")) {
    if (!h5exists(h5.file, "clustering")) rhdf5::h5createGroup(h5.file, "clustering")
    rhdf5::h5ls(h5.file)
    dbg("[pgx.createSignatureDatabaseH5.fromMatrix] compute clustering...")

    pos <- pgx.clusterBigMatrix(
      abs(X), ## on absolute foldchange!!
      methods = c("pca", "tsne", "umap"),
      dims = c(2, 3),
      reduce.sd = 2000,
      reduce.pca = min(200, round(ncol(X) / 3))
    )

    rhdf5::h5write(pos[["pca2d"]], h5.file, "clustering/pca2d") ## can write list??
    rhdf5::h5write(pos[["pca3d"]], h5.file, "clustering/pca3d") ## can write list??
    rhdf5::h5write(pos[["tsne2d"]], h5.file, "clustering/tsne2d") ## can write list??
    rhdf5::h5write(pos[["tsne3d"]], h5.file, "clustering/tsne3d") ## can write list??
    rhdf5::h5write(pos[["umap2d"]], h5.file, "clustering/umap2d") ## can write list??
    rhdf5::h5write(pos[["umap3d"]], h5.file, "clustering/umap3d") ## can write list??
  }

  ## --------------------------------------------------
  ## Add enrichment signatures
  ## --------------------------------------------------
  pgx.addEnrichmentSignaturesH5(h5.file, X = X, mc.cores = 0, methods = "rankcor")

  ## done!
  rhdf5::h5closeAll()
}


#' @title Add Enrichment Signatures to HDF5 File
#'
#' @description
#' Adds gene set enrichment results as signatures to an HDF5 file.
#'
#' @param h5.file The HDF5 file path to write enrichment results to.
#' @param X The gene expression matrix. If NULL, matrix is read from HDF5 file.
#' @param mc.cores Number of cores for parallel processing. Default 0.
#' @param methods Enrichment methods to run. Default c("gsea", "rankcor")[2].
#'
#' @details
#' This function runs gene set enrichment analysis on a gene expression matrix X,
#' using the specified methods. It writes the enrichment results as signatures to
#' the specified HDF5 file under the "enrichment" group.
#'
#' Parallel processing can be enabled by setting mc.cores to use multiple CPU cores.
#' By default, it runs only the "rankcor" method. Alternative methods like "gsea" can
#' be specified in the methods parameter.
#'
#' @return
#' No value is returned. Enrichment signatures are written to the HDF5 file.
#'
#' @export
pgx.addEnrichmentSignaturesH5 <- function(h5.file, X = NULL, mc.cores = 0,
                                          methods = c("gsea", "rankcor")[2]) {
  h5exists <- function(h5.file, obj) {
    xobjs <- apply(rhdf5::h5ls(h5.file)[, 1:2], 1, paste, collapse = "/")
    obj %in% gsub("^/|^//", "", xobjs)
  }

  if (is.null(X)) {
    X <- rhdf5::h5read(h5.file, "data/matrix")
    rn <- rhdf5::h5read(h5.file, "data/rownames")
    cn <- rhdf5::h5read(h5.file, "data/colnames")
    rownames(X) <- rn
    colnames(X) <- cn
    X[which(X < -999999)] <- NA
  }

  ## ---------------------- ONLY HALLMARK FOR NOW -----------------------

  G <- playdata::GSET_SPARSEG_XL
  sel <- grep("HALLMARK|C[1-9]|^GO", rownames(G))
  sel <- grep("HALLMARK", rownames(G))
  genes <- intersect(colnames(G), rownames(X))
  G <- G[sel, genes, drop = FALSE]
  X <- X[genes, ]
  gmt <- apply(G, 1, function(x) colnames(G)[which(x != 0)])

  if (!h5exists(h5.file, "enrichment")) {
    rhdf5::h5createGroup(h5.file, "enrichment")
  }
  if (h5exists(h5.file, "enrichment/genesets")) {
    rhdf5::h5delete(h5.file, "enrichment/genesets")
  }
  rhdf5::h5write(names(gmt), h5.file, "enrichment/genesets")

  if ("gsea" %in% methods) {
    cat("[pgx.addEnrichmentSignaturesH5] starting fGSEA for", length(gmt), "gene sets...")
    i <- 1
    F1 <- lapply(colnames(X), function(i) {
      xi <- X[, i]
      xi <- xi[!is.na(xi)]

      xi <- xi + 1e-3 * stats::rnorm(length(xi))
      suppressMessages(suppressWarnings(
        res1 <- fgsea::fgseaSimple(gmt, xi, nperm = 10000, nproc = mc.cores)
      ))
      r <- res1$NES
      names(r) <- res1$pathway
      r
    })
    F1 <- do.call(cbind, F1)
    F1 <- F1[match(names(gmt), rownames(F1)), , drop = FALSE]
    F1[is.na(F1)] <- 0
    colnames(F1) <- colnames(X)
    if (h5exists(h5.file, "enrichment/GSEA")) rhdf5::h5delete(h5.file, "enrichment/GSEA")
    rhdf5::h5write(F1, h5.file, "enrichment/GSEA")
  }

  if ("rankcor" %in% methods) {
    F3 <- gset.rankcor(X, Matrix::t(G))$rho
    F3 <- F3[match(names(gmt), rownames(F3)), , drop = FALSE]
    F3[is.na(F3)] <- 0
    if (h5exists(h5.file, "enrichment/rankcor")) rhdf5::h5delete(h5.file, "enrichment/rankcor")
    rhdf5::h5write(F3, h5.file, "enrichment/rankcor")
  }

  rhdf5::h5ls(h5.file)
  rhdf5::h5closeAll()

  cat("[pgx.addEnrichmentSignaturesH5] done!\n")
}


## -------------------------------------------------------------------
## Pre-calculate geneset expression with different methods
## -------------------------------------------------------------------


#' @title Compute gene set expression
#'
#' @description This function calculates gene set expression scores from an input gene expression matrix using different methods like GSVA, ssGSEA, etc.
#'
#' @param X Numeric gene expression matrix with genes in rows and samples in columns.
#' @param gmt List of gene sets, each containing gene names.
#' @param method Character vector of methods to use for computing gene set expression. Options are "gsva", "ssgsea", "spearman", "average". Default uses all.
#' @param min.size Minimum gene set size required. Sets with fewer genes are removed. Default 10.
#' @param center Logical indicating whether gene expression should be centered first. Default is TRUE.
#'
#' @details This function takes a gene expression matrix X and gene sets defined in gmt. It calculates gene set expression scores using the specified methods.
#'
#' By default it computes GSVA, ssGSEA, spearman correlation and mean expression scores. X is centered by gene means if center=TRUE. Small gene sets are filtered out.
#'
#' The output is a list containing the gene set expression matrices from each method.
#'
#' @return A list with gene set expression matrices, named by method.
#'
#' @export
pgx.computeGeneSetExpression <- function(X, gmt, method = NULL,
                                         min.size = 10, center = TRUE) {
  ALL.METHODS <- c("gsva", "spearman", "average")
  ALL.METHODS <- c("gsva", "ssgsea", "spearman", "average")
  if (is.null(method)) {
    method <- ALL.METHODS
  }
  ## this is important!!! centering on genes (GSVA does)
  if (center) {
    X <- X - rowMeans(X, na.rm = TRUE)
  }

  gmt.size <- sapply(gmt, function(x) sum(x %in% rownames(X)))
  gmt <- gmt[gmt.size >= min.size]

  S <- list()
  if ("gsva" %in% method) {
    S[["gsva"]] <- GSVA::gsva(X, gmt, method = "gsva")
  }
  if ("ssgsea" %in% method) {
    S[["ssgsea"]] <- GSVA::gsva(X, gmt, method = "ssgsea", min.sz = 1)
  }
  if (any(method %in% c("spearman", "average"))) {
    gg <- rownames(X)
    G <- gmt2mat(gmt, bg = gg)
    if ("spearman" %in% method) {
      rho <- t(G[gg, ]) %*% scale(apply(X[gg, ], 2, rank)) / sqrt(nrow(X) - 1)
      rho[is.na(rho)] <- 0
      S[["spearman"]] <- rho
    }
    if ("average" %in% method) {
      avg.X <- t(G[gg, ]) %*% X[gg, ] / Matrix::colSums(G[gg, ])
      avg.X[is.na(avg.X)] <- 0
      S[["average"]] <- avg.X
    }
  }

  ## compute meta score
  S1 <- lapply(S, function(x) apply(x, 2, rank)) ## rank by sample
  S[["meta"]] <- scale(Reduce("+", S1) / length(S1))
  gs <- Reduce(intersect, lapply(S, rownames))
  S <- lapply(S, function(x) x[gs, ])

  return(S)
}

## ================================================================================
## ========================= SIGDB H5 FUNCTIONS ===================================
## ================================================================================

## IK 2.6.23: going to rewrite following functions without global
## SIGDB.DIR. For the moment seems these functions are duplicated in


#' @title Get the SIGDB Directory
#'
#' @description This function returns the path to the SIGDB directory, which is used to store data related to the SIGDB database.
#'
#' @details The function first checks if the `SIGDB.DIR` variable exists in the global environment.
#' If it does not exist, the function checks if the `FILESX` variable exists in the global environment.
#' If `FILESX` exists, the function constructs the path to the SIGDB directory by concatenating `FILESX` and a subdirectory named "sigdb".
#' The resulting path is returned as a character vector of length 2, containing both `FILESX` and the full path to the SIGDB directory.
#' If `SIGDB.DIR` already exists in the global environment, its value is returned directly.
#'
#' @return A character vector containing the path to the SIGDB directory.
#'
getSIGDB.DIR <- function() {
  if (!exists("SIGDB.DIR") && exists("FILESX")) {
    SIGDB.DIR <- c(FILESX, file.path(FILESX, "sigdb"))
  } else {
    SIGDB.DIR
  }
}


#' @title Get Connectivity Contrasts from a SIGDB File
#'
#' @description This function reads the connectivity contrasts from a specified SIGDB file and returns them as a character vector.
#'
#' @param sigdb A character string specifying the name of the SIGDB file.
#' @param path An optional character string specifying the path to the directory containing the SIGDB file (default is NULL).
#'
#' @details The function first checks if a path to the SIGDB file is provided.
#' If a path is provided, it constructs the full path to the SIGDB file by concatenating the path and the file name.
#' The function then checks if the specified SIGDB file exists.
#' If the file exists, it uses the `h5read` function from the rhdf5 package to read the connectivity contrasts from the "data/colnames" dataset in the SIGDB file.
#'
#' @return A character vector containing the connectivity contrasts from the specified SIGDB file, or NULL if the file does not exist.
#'
#' @export
sigdb.getConnectivityContrasts <- function(sigdb, path = NULL) {
  if (!is.null(path)) {
    sigdb <- file.path(path, sigdb)
  }
  if (!file.exists(sigdb)) {
    return(NULL)
  }
  rhdf5::h5read(sigdb, "data/colnames")
}


#' @title Get elements from a sigDB
#'
#' @param sigdb Path to signature database HDF5 file
#' @param select Character vector of contrasts to return
#' @param genes Character vector of genes to return
#' @param path Optional path prefix for sigdb file
#'
#' @return Matrix of connectivity scores
#'
#' @description Get the connectivity matrix from a signature database HDF5 file.
#'
#' @details This function extracts the connectivity matrix from a signature database
#' stored in HDF5 format. The HDF5 file should contain the matrix under /data/matrix.
#'
#' You can filter the returned matrix by specifying a subset of contrasts to return
#' with the \code{select} parameter. Use the \code{genes} parameter to return only a
#' subset of genes.
#'
#' If the sigdb file is not in the working directory, set the \code{path} parameter.
#'
#' @export
sigdb.getConnectivityMatrix <- function(sigdb, select = NULL, genes = NULL, path = NULL) {
  if (sigdb == "" || is.null(sigdb)) {
    warning("[getConnectivityMatrix] ***WARNING*** sigdb=", sigdb)
    return(NULL)
  }
  if (!is.null(path)) {
    sigdb <- file.path(path, sigdb)
  }
  if (!file.exists(sigdb)) {
    warning("[getConnectivityMatrix] ***WARNING*** file ", sigdb, "not found")
    return(NULL)
  }


  if (grepl("csv$", sigdb)) {
    X <- utils::read.csv(sigdb, row.names = 1, check.names = FALSE)
    X <- as.matrix(X)
    X <- X[, colMeans(is.na(X)) < 0.99, drop = FALSE] ## omit empty columns
    if (!is.null(genes)) X <- X[intersect(genes, rownames(X)), , drop = FALSE]
    if (!is.null(select)) X <- X[, intersect(select, colnames(X))]
  } else if (grepl("h5$", sigdb)) {
    cn <- rhdf5::h5read(sigdb, "data/colnames")
    rn <- rhdf5::h5read(sigdb, "data/rownames")
    rowidx <- 1:length(rn)
    colidx <- 1:length(cn)
    if (!is.null(genes)) rowidx <- match(intersect(genes, rn), rn)
    if (!is.null(select)) colidx <- match(intersect(select, cn), cn)
    nr <- length(rowidx)
    nc <- length(colidx)
    dbg("[sigdb.getConnectivityMatrix] *WARNING* reading large H5 file:", nr, "x", nc, "")
    X <- rhdf5::h5read(sigdb, "data/matrix", index = list(rowidx, colidx))
    rownames(X) <- rn[rowidx]
    colnames(X) <- cn[colidx]
  } else {
    cat("[getConnectivityMatrix] WARNING: could not retrieve matrix\n")
  }
  return(X)
}


#' @describeIn sigdb.getConnectivityMatrix This function reads the enrichment data stored in a signature database in HDF5 format.
#' It extracts either the gene set enrichment analysis (GSEA) results or the rank
#' correlation analysis results.
#' @param which Method to use for enrichment, either "gsea" or "rankcor"
#'
#' @export
sigdb.getEnrichmentMatrix <- function(sigdb, select = NULL, path = NULL,
                                      which = c("gsea", "rankcor")[1]) {
  if (sigdb == "" || is.null(sigdb)) {
    dbg("[getEnrichmentMatrix] *WARNING* could not find sigdb=", sigdb)
    return(NULL)
  }
  if (!grepl("h5$", sigdb)) {
    stop("getEnrichmentMatrix:: only for H5 database files")
    return(NULL)
  }
  if (!is.null(path)) {
    sigdb <- file.path(path, sigdb)
  }
  if (!file.exists(sigdb)) {
    warning("[sigdb.getEnrichmentMatrix] *WARNING* file ", sigdb, "not found")
    return(NULL)
  }

  h5exists <- function(h5.file, obj) {
    xobjs <- apply(rhdf5::h5ls(h5.file)[, 1:2], 1, paste, collapse = "/")
    obj %in% gsub("^/|^//", "", xobjs)
  }

  cn <- rhdf5::h5read(sigdb, "data/colnames")
  has.gs <- h5exists(sigdb, "enrichment/genesets")
  if (!has.gs) {
    dbg("[getEnrichmentMatrix] WARNING: PGX object has no enrichment results")
    return(NULL)
  }

  rn <- rhdf5::h5read(sigdb, "enrichment/genesets")
  rowidx <- 1:length(rn)
  colidx <- 1:length(cn)
  if (!is.null(select)) colidx <- match(intersect(select, cn), cn)

  has.gsea <- h5exists(sigdb, "enrichment/GSEA")
  has.rankcor <- h5exists(sigdb, "enrichment/rankcor")
  Y <- NULL
  if (is.null(Y) && "gsea" %in% which && has.gsea) {
    Y <- rhdf5::h5read(sigdb, "enrichment/GSEA", index = list(rowidx, colidx))
  }
  if (is.null(Y) && "rankcor" %in% which && has.rankcor) {
    Y <- rhdf5::h5read(sigdb, "enrichment/rankcor", index = list(rowidx, colidx))
  }
  if (is.null(Y) || nrow(Y) == 0) {
    return(NULL)
  }
  rownames(Y) <- rn[rowidx]
  colnames(Y) <- cn[colidx]

  return(Y)
}

#' @describeIn sigdb.getConnectivityMatrix reads the up/down regulated gene signatures from a signature database HDF5 file.
#' @export
sigdb.getSignatureMatrix <- function(sigdb, path = NULL) {
  if (sigdb == "" || is.null(sigdb)) {
    dbg("[getEnrichmentMatrix] ***WARNING*** sigdb=", sigdb)
    return(NULL)
  }
  if (!grepl("h5$", sigdb)) {
    stop("getEnrichmentMatrix:: only for H5 database files")
  }
  if (!is.null(path)) {
    sigdb <- file.path(path, sigdb)
  }
  if (!file.exists(sigdb)) {
    warning("[sigdb.getEnrichmentMatrix] *WARNING* file ", sigdb, "not found")
    return(NULL)
  }

  rhdf5::h5ls(sigdb)
  cn <- rhdf5::h5read(sigdb, "data/colnames")
  dn <- rhdf5::h5read(sigdb, "signature/sig100.dn")
  up <- rhdf5::h5read(sigdb, "signature/sig100.up")
  colnames(dn) <- cn
  colnames(up) <- cn

  list(up = up, dn = dn)
}

#' @describeIn sigdb.getConnectivityMatrix removes a dataset and its associated data from a signature database HDF5 file.
#' @export
sigdb.removeDataset <- function(h5.file, pgxname) {
  ## delete columns from H5 file
  dd <- rhdf5::h5ls(h5.file)
  cn <- rhdf5::h5read(h5.file, "data/colnames")
  del <- grep(paste0("\\[", pgxname, "\\]"), cn)
  slots <- sub("^/*", "", paste0(dd$group, "/", dd$name))
  h5.removeRows <- function(h5.file, slot, rows) {
    if (!slot %in% slots) {
      return()
    }
    X <- rhdf5::h5read(h5.file, slot)
    if (length(dim(X)) == 1) X <- X[-rows]
    if (length(dim(X)) > 1) X <- X[-rows, , drop = FALSE]
    rhdf5::h5delete(h5.file, slot)
    rhdf5::h5write(X, h5.file, slot)
  }
  h5.removeCols <- function(h5.file, slot, cols, chunk = TRUE) {
    if (!slot %in% slots) {
      return()
    }
    X <- rhdf5::h5read(h5.file, slot)
    X <- X[, -cols, drop = FALSE]
    rhdf5::h5delete(h5.file, slot)
    if (!chunk) {
      rhdf5::h5write(X, h5.file, slot)
    } else {
      rhdf5::h5createDataset(
        h5.file, slot,
        c(nrow(X), ncol(X)),
        chunk = c(nrow(X), 1),
        level = 7
      )
      rhdf5::h5write(
        X,
        file = h5.file,
        name = "data/matrix",
        index = list(1:nrow(X), 1:ncol(X))
      )
    }
  }

  rhdf5::h5ls(h5.file)
  if (!is.null(del) && length(del) > 0) {
    h5.removeRows(h5.file, "clustering/pca2d", del)
    h5.removeRows(h5.file, "clustering/pca3d", del)
    h5.removeRows(h5.file, "clustering/tsne2d", del)
    h5.removeRows(h5.file, "clustering/tsne3d", del)
    h5.removeRows(h5.file, "clustering/umap2d", del)
    h5.removeRows(h5.file, "clustering/umap3d", del)
    h5.removeRows(h5.file, "data/colnames", del)
    h5.removeCols(h5.file, "data/matrix", del)
    h5.removeCols(h5.file, "enrichment/gsea", del)
    h5.removeCols(h5.file, "enrichment/rankcor", del)
    h5.removeCols(h5.file, "signature/sig100.dn", del)
    h5.removeCols(h5.file, "signature/sig100.up", del)
  }
}


## ================================================================================
## =============================== END OF FILE ====================================
## ================================================================================
