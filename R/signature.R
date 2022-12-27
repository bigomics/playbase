# pgx.computeConnectivityScores
compute_connectivity_scores <- function(ngs, sigdb, ntop = 1000, contrasts = NULL,
                                        remove.le = FALSE, inmemory = FALSE) {
  meta <- get_meta_fold_change_matrix(ngs, what = "meta")
  colnames(meta$fc)

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
    cat("[pgx.computeConnectivityScores] ERROR: could not H5 file\n")
    return(NULL)
  }

  if (is.null(contrasts)) {
    contrasts <- colnames(meta$fc)
  }
  contrasts <- intersect(contrasts, colnames(meta$fc))

  if (inmemory) {
    F <- meta$fc[, contrasts]
    scores <- correlate_signature_h5_inmemory(
      meta$fc,
      h5.file = h5.file,
      nsig = 100, ntop = ntop, nperm = 9999
    )
  } else {
    scores <- list()
    ct <- contrasts[1]
    for (ct in contrasts) {
      fc <- meta$fc[, ct]
      names(fc) <- rownames(meta$fc)
      names(fc) <- toupper(names(fc)) ## for MOUSE!!
      res <- correlate_signature_h5(
        fc,
        h5.file = h5.file,
        nsig = 100, ntop = ntop, nperm = 9999
      )
      dim(res)
      scores[[ct]] <- res
    }
  }
  if (is.null(names(scores))) names(scores) <- contrasts

  ## remove leadingEdge (take too much memory!!!)
  if (remove.le) {
    for (j in 1:length(scores)) scores[[j]]$leadingEdge <- NULL
  }

  names(scores)
  return(scores)
}

# pgx.correlateSignatureH5.inmemory
correlate_signature_h5_inmemory <- function(F, h5.file, nsig = 100, ntop = 1000, nperm = 1000) {
  ##
  ##
  ##
  ##

  if (NCOL(F) == 1 && class(F) == "numeric") {
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

  mem1 <- round(object.size(matG) / 1e9, 2)
  cat("[pgx.correlateSignatureH5] object.size(matG)=", mem1, "Gb\n") ## gigabytes....

  ## ---------------------------------------------------------------
  ## Compute simple correlation between query profile and signatures
  ## ---------------------------------------------------------------
  res <- list()
  i <- 1
  for (i in 1:ncol(F)) {
    gg <- intersect(rownames(F), rn)
    fc1 <- sort(F[gg, i])
    gg <- unique(names(c(Matrix::head(fc1, nsig), Matrix::tail(fc1, nsig))))
    ## gg <- intersect(gg,rn)
    remove(fc1)
    row.idx <- match(gg, rn)
    length(row.idx)

    ## G <- rhdf5::h5read(h5.file, "data/matrix", index=list(row.idx,1:length(cn)))  ### SLOW!!!
    rG <- matG[row.idx, , drop = FALSE]
    rG <- apply(rG, 2, rank, na.last = "keep")
    dim(rG)
    dimnames(rG) <- list(rn[row.idx], cn)

    ## this FC signature
    fc <- F[, i]
    ## rG  <- apply( G[gg,], 2, rank, na.last="keep" )
    rfc <- rank(fc[gg], na.last = "keep") ## rank correlation??
    ## rho <- stats::cor(rG, rfc, use="pairwise")[,1]
    rG[is.na(rG)] <- 0 ## NEED RETHINK: are missing values to be treated as zero???
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
    ## head(sig100.up,2)

    ## combine up/down into one (unsigned GSEA test)
    gmt <- rbind(sig100.up, sig100.dn)
    gmt <- unlist(apply(gmt, 2, list), recursive = FALSE)
    names(gmt) <- cn[sel.idx]
    length(gmt)

    ## use entire fc vector
    system.time(res1 <- fgsea::fgseaSimple(gmt, abs(fc), nperm = nperm)) ## really unsigned???
    dim(res1)

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

  if (0) {
    res$rho.p <- cor.pvalue(res$rho, n = length(gg))
    res$meta.p <- apply(res[, c("pval", "rho.p")], 1, function(p) metap::sumz(p)$p)
    res <- res[order(res$meta.p), ]
  }
  remove(matG)
  gc()

  return(res)
}

# pgx.correlateSignatureH5
correlate_signature_h5 <- function(fc, h5.file, nsig = 100, ntop = 1000, nperm = 10000,
                                   h5.data = "data/matrix", h5.rn = "data/rownames",
                                   h5.cn = "data/colnames") {
  ##
  ##
  ##
  ##


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
  fc1 <- sort(fc[gg])
  gg <- unique(names(c(Matrix::head(fc1, nsig), Matrix::tail(fc1, nsig))))
  ## gg <- intersect(gg,rn)
  length(gg)
  row.idx <- match(gg, rn)
  rhdf5::h5closeAll()
  G <- rhdf5::h5read(h5.file, "data/matrix", index = list(row.idx, 1:length(cn)))
  dim(G)
  ## head(G[,1])
  G[which(G < -999999)] <- NA
  ## G[is.na(G)] <- 0  ## NEED RETHINK: are missing values to be treated as zero???
  dim(G)
  dimnames(G) <- list(rn[row.idx], cn)

  ## rank correlation??
  rG <- apply(G[gg, ], 2, rank, na.last = "keep")
  rfc <- rank(fc[gg], na.last = "keep")
  ## rho <- stats::cor(rG, rfc, use="pairwise")[,1]
  rG[is.na(rG)] <- 0 ## NEED RETHINK: are missing values to be treated as zero???
  rfc[is.na(rfc)] <- 0
  suppressWarnings(rho <- stats::cor(rG, rfc, use = "pairwise")[, 1])

  remove(G, rG, rfc)

  ## --------------------------------------------------
  ## test tops signatures using fGSEA
  ## --------------------------------------------------

  sel <- Matrix::head(names(sort(-abs(rho))), ntop)
  sel.idx <- match(sel, cn)
  sig100.up <- rhdf5::h5read(h5.file, "signature/sig100.up",
    index = list(1:100, sel.idx)
  )
  sig100.dn <- rhdf5::h5read(h5.file, "signature/sig100.dn",
    index = list(1:100, sel.idx)
  )
  ## head(sig100.up,2)

  ## combine up/down into one (unsigned GSEA test)
  gmt <- rbind(sig100.up, sig100.dn)
  gmt <- unlist(apply(gmt, 2, list), recursive = FALSE)
  names(gmt) <- cn[sel.idx]
  length(gmt)

  ## system.time( res <- fgsea::fgsea(gmt, fc, nperm=10000))
  system.time(res <- fgsea::fgseaSimple(gmt, abs(fc), nperm = nperm)) ## really unsigned???
  dim(res)

  ## ---------------------------------------------------------------
  ## Combine correlation+GSEA by combined score (NES*rho)
  ## ---------------------------------------------------------------
  jj <- match(res$pathway, names(rho))
  res$rho <- rho[jj]
  res$R2 <- rho[jj]**2
  res$score <- res$R2 * res$NES
  res <- res[order(res$score, decreasing = TRUE), ]

  if (0) {
    res$rho.p <- cor.pvalue(res$rho, n = length(gg))
    res$meta.p <- apply(res[, c("pval", "rho.p")], 1, function(p) metap::sumz(p)$p)
    res <- res[order(res$meta.p), ]
  }

  Matrix::head(res)
  return(res)
}
