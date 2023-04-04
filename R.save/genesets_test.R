# compute.testGenesets
test_genesets <- function(pgx, max.features = 1000, lib.dir = "../lib",
                          test.methods = c("gsva", "camera", "fgsea"),
                          remove.outputs = TRUE) {
  ## Rewritten 24.12.2019. Now much faster but needs gset-sparseG-XL
  ## precomputed.
  ##
  ##
  if (!"X" %in% names(pgx)) {
    stop("FATAL : object must have normalized matrix X")
  }

  ## -----------------------------------------------------------
  ## Load huge geneset matrix
  ## -----------------------------------------------------------
  G <- playbase::GSET_SPARSEG_XL # now an internal dataset
  G <- Matrix::t(G)

  ## -----------------------------------------------------------
  ## Filter genes
  ## -----------------------------------------------------------

  ## filter genes only in dataset

  ##
  ## GENE.TITLE = unlist(as.list(org.Hs.egGENENAME))
  ## genes = Matrix::head(as.character(unlist(as.list(org.Hs.egSYMBOL))),1000)
  genes <- unique(as.character(pgx$genes$gene_name))
  genes <- toupper(genes) ## handle mouse genes...
  G <- G[rownames(G) %in% genes, ]

  ## -----------------------------------------------------------
  ## Filter gene sets
  ## -----------------------------------------------------------

  ## filter gene sets on size
  gmt.size <- Matrix::colSums(G != 0)
  size.ok <- (gmt.size >= 15 & gmt.size <= 1000)
  G <- G[, which(size.ok)]

  ## -----------------------------------------------------------
  ## create the full GENE matrix (always collapsed by gene)
  ## -----------------------------------------------------------

  single.omics <- !any(grepl("\\[", rownames(pgx$counts)))
  single.omics <- TRUE ## !!! for now...
  if (single.omics) {
    ## normalized matrix
    X <- pgx$X
  } else {
    data.type <- gsub("\\[|\\].*", "", rownames(pgx$counts))
    jj <- which(data.type %in% c("gx", "mrna"))
    if (length(jj) == 0) {
      stop("FATAL. could not find gx/mrna values.")
    }
    X <- pgx$X[jj, ]
  }

  ## if reduced samples
  ss <- rownames(pgx$model.parameters$exp.matrix)
  X <- X[, ss, drop = FALSE]

  ## -----------------------------------------------------------
  ## create the GENESETxGENE matrix
  ## -----------------------------------------------------------
  gg <- toupper(rownames(X)) ## accomodate for mouse...
  ii <- intersect(gg, rownames(G))
  G <- G[ii, , drop = FALSE]
  xx <- setdiff(gg, rownames(G))
  matX <- Matrix::Matrix(0, nrow = length(xx), ncol = ncol(G), sparse = TRUE)
  rownames(matX) <- xx
  colnames(matX) <- colnames(G)
  G <- rbind(G, matX)
  G <- G[match(gg, rownames(G)), , drop = FALSE]
  rownames(G) <- rownames(X) ## original name (e.g. mouse)

  ## -----------------------------------------------------------
  ## Prioritize gene sets by fast rank-correlation
  ## -----------------------------------------------------------
  if (is.null(max.features)) max.features <- 20000
  if (max.features < 0) max.features <- 20000

  if (max.features > 0) {
    ## Reduce gene sets by selecting top varying genesets. We use the
    ## very fast sparse rank-correlation for approximate single sample
    ## geneset activation.
    cX <- X - rowMeans(X, na.rm = TRUE) ## center!
    gsetX <- qlcMatrix::corSparse(G[, ], apply(cX[, ], 2, rank))
    grp <- pgx$model.parameters$group
    gsetX.bygroup <- NULL
    if (!is.null(grp)) {
      gsetX.bygroup <- Matrix::t(apply(gsetX, 1, function(x) tapply(x, grp, mean)))
      sdx <- apply(gsetX.bygroup, 1, sd)
    } else {
      sdx <- apply(gsetX, 1, sd)
    }
    names(sdx) <- colnames(G)
    jj <- Matrix::head(order(-sdx), max.features)
    must.include <- "hallmark|kegg|^go|^celltype"
    jj <- unique(c(jj, grep(must.include, colnames(G), ignore.case = TRUE)))
    jj <- jj[order(colnames(G)[jj])]
    G <- G[, jj, drop = FALSE]
    remove(gsetX)
    remove(gsetX.bygroup)
  }

  ## -----------------------------------------------------------
  ## get design and contrast matrix
  ## -----------------------------------------------------------
  design <- pgx$model.parameters$design
  contr.matrix <- pgx$model.parameters$contr.matrix

  ALL.GSET.METHODS <- c(
    "fisher", "ssgsea", "gsva", "spearman", "camera", "fry",
    "fgsea", "gsea.permPH", "gsea.permGS", "gseaPR"
  )

  ## -----------------------------------------------------------
  ## Run methods
  ## -----------------------------------------------------------

  ## convert to gene list
  gmt <- lapply(apply(G != 0, 2, which), names)
  Y <- pgx$samples
  gc()

  # gset.fitContrastsWithAllMethods(
  gset.meta <- fit_geneset_contrasts_with_all_methods(
    gmt = gmt, X = X, Y = Y, G = G,
    design = design, ## genes=GENES,
    contr.matrix = contr.matrix, methods = test.methods,
    mc.threads = 1, mc.cores = NULL, batch.correct = TRUE
  )

  rownames(gset.meta$timings) <- paste("[test.genesets]", rownames(gset.meta$timings))

  pgx$timings <- rbind(pgx$timings, gset.meta$timings)
  pgx$gset.meta <- gset.meta

  pgx$gsetX <- pgx$gset.meta$matrices[["meta"]] ## META or average FC??!
  pgx$GMT <- G[, rownames(pgx$gsetX)]

  ## -----------------------------------------------------------------------
  ## ------------------------ clean up -------------------------------------
  ## -----------------------------------------------------------------------

  ## remove large outputs... (uncomment if needed)
  if (remove.outputs) {
    pgx$gset.meta$outputs <- NULL
  }

  remove(X)
  remove(Y)
  remove(G)
  remove(gmt)

  return(pgx)
}
