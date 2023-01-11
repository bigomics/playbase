# compute.extra
compute_extra <- function(pgx,
                          extra = c("meta.go", "deconv", "infer",
                                    "drugs", "wordcloud"),
                          lib.dir = '../lib', sigdb = NULL) {
  compute_extra_helper(pgx, extra = extra, lib.dir = lib.dir, sigdb = sigdb)
}

# pgx.computeExtra
compute_extra_helper <- function(pgx, extra = c(
                                   "meta.go", "deconv", "infer",
                                   "drugs", "wordcloud"
                                 ),
                                 lib.dir = '../lib', sigdb = NULL) {
  timings <- c()

  libx.dir <- '../libx'

  EXTRA.MODULES = c("meta.go", "deconv", "infer", "drugs", ## "graph",
                    "connectivity", "wordcloud")

  extra <- intersect(extra, EXTRA.MODULES)
  if (length(extra) == 0) {
    return(pgx)
  }

  ## detect if it is single or multi-omics
  single.omics <- !any(grepl("\\[", rownames(pgx$counts)))
  if (single.omics) {
    rna.counts <- pgx$counts
  } else {
    data.type <- gsub("\\[|\\].*", "", rownames(pgx$counts))
    jj <- which(data.type %in% c("gx", "mrna"))
    length(jj)
    if (length(jj) == 0) {
      stop("FATAL. could not find gx/mrna values.")
    }
    rna.counts <- pgx$counts[jj, ]
    ## rownames(rna.counts) <- gsub(".*:|.*\\]","",rownames(rna.counts))
    is.logged <- (min(rna.counts, na.rm = TRUE) < 0 || max(rna.counts, na.rm = TRUE) < 50)
    if (is.logged) {
      rna.counts <- 2**rna.counts
    }
  }

  if ("meta.go" %in% extra) {
    tt <- system.time({
      pgx$meta.go <- compute_core_go_graph(pgx, fdr = 0.20)
    })
    timings <- rbind(timings, c("meta.go", tt))
  }

  if ("deconv" %in% extra) {
    tt <- system.time({
      pgx <- compute_deconvolution(
        pgx,
        lib.dir = lib.dir[1], rna.counts = rna.counts,
        full = FALSE
      )
    })
    timings <- rbind(timings, c("deconv", tt))
  }

  if ("infer" %in% extra) {
    tt <- system.time({
      pgx <- compute_cellcycle_gender(pgx, rna.counts = rna.counts)
    })
    timings <- rbind(timings, c("infer", tt))
  }

  if ("drugs" %in% extra) {
    pgx$drugs <- NULL ## reset??
    cmap.dir <- file.path(libx.dir, "cmap")
    if (!dir.exists(cmap.dir)) {
      cmap.dir <- file.path(lib.dir, "cmap") ## look for default lib
    }
    if (!dir.exists(cmap.dir)) {
      warning("Warning:: missing CMAP files. Skipping drug connectivity analysis!")
    }

    if (dir.exists(cmap.dir)) {
      tt <- system.time({
        pgx <- compute_drug_activity_enrichment(pgx, cmap.dir)
      })
      timings <- rbind(timings, c("drugs", tt))

      tt <- system.time({
        pgx <- compute_drug_sensitivity_enrichment(pgx, cmap.dir)
      })
      timings <- rbind(timings, c("drugs-sx", tt))
    }
  }

  if ("graph" %in% extra) {
    tt <- system.time({
      pgx <- compute_omics_graphs(pgx)
    })
    timings <- rbind(timings, c("graph", tt))
  }

  if ("wordcloud" %in% extra) {
    tt <- system.time({
      res <- calculate_word_cloud(pgx, progress = NULL, pg.unit = 1)
    })
    timings <- rbind(timings, c("wordcloud", tt))
    pgx$wordcloud <- res
    remove(res)
  }

  if ("connectivity" %in% extra) {
    if (is.null(sigdb)) {
      lib.dir2 <- unique(c(lib.dir, libx.dir)) ### NEED BETTER SOLUTION!!!
      sig.dir <- c(SIGDB.DIR, lib.dir2)
      sigdb <- dir(sig.dir, pattern = "^sigdb-.*h5$", full.names = TRUE)
      sigdb
    }

    db <- sigdb[1]
    for (db in sigdb) {
      if (file.exists(db)) {
        ntop <- 1000
        ## in memory for many comparisons
        meta <- get_meta_fold_change_matrix(pgx, what = "meta")
        inmemory <- ifelse(ncol(meta$fc) > 50, TRUE, FALSE)
        tt <- system.time({
          scores <- compute_connectivity_scores(
            pgx, db,
            ntop = ntop, contrasts = NULL,
            remove.le = TRUE, inmemory = inmemory
          )
        })
        timings <- rbind(timings, c("connectivity", tt))

        db0 <- sub(".*/", "", db)
        pgx$connectivity[[db0]] <- scores
        remove(scores)
      }
    }
    names(pgx$connectivity)
  }

  ## ------------------------------------------------------
  ## pretty collapse all timings
  ## ------------------------------------------------------
  ## timings0 <- do.call(rbind, timings)
  timings <- as.matrix(timings)
  rownames(timings) <- timings[, 1]
  timings0 <- apply(as.matrix(timings[, -1, drop = FALSE]), 2, as.numeric)
  if (nrow(timings) == 1) {
    timings0 <- matrix(timings0, nrow = 1)
    colnames(timings0) <- colnames(timings)[-1]
    rownames(timings0) <- rownames(timings)
  }
  rownames(timings0) <- rownames(timings)
  timings0 <- apply(timings0, 2, function(x) tapply(x, rownames(timings0), sum))
  if (is.null(nrow(timings0))) {
    cn <- names(timings0)
    rn <- unique(rownames(timings))
    timings0 <- matrix(timings0, nrow = 1)
    colnames(timings0) <- cn
    rownames(timings0) <- rn[1]
  }
  rownames(timings0) <- paste("[extra]", rownames(timings0))

  pgx$timings <- rbind(pgx$timings, timings0)
  return(pgx)
}

# compute.deconvolution
compute_deconvolution <- function(pgx, lib.dir, rna.counts = pgx$counts, full = FALSE) {
  ## list of reference matrices
  refmat <- list()
  readSIG <- function(f) read.csv(file.path(lib.dir, "sig", f), row.names = 1, check.names = FALSE)
  LM22 <- read.csv(file.path(lib.dir, "sig/LM22.txt"), sep = "\t", row.names = 1)
  refmat[["Immune cell (LM22)"]] <- LM22
  refmat[["Immune cell (ImmProt)"]] <- readSIG("immprot-signature1000.csv")
  refmat[["Immune cell (DICE)"]] <- readSIG("DICE-signature1000.csv")
  refmat[["Immune cell (ImmunoStates)"]] <- readSIG("ImmunoStates_matrix.csv")
  refmat[["Tissue (HPA)"]] <- readSIG("rna_tissue_matrix.csv")
  refmat[["Tissue (GTEx)"]] <- readSIG("GTEx_rna_tissue_tpm.csv")
  refmat[["Cell line (HPA)"]] <- readSIG("HPA_rna_celline.csv")
  refmat[["Cell line (CCLE)"]] <- readSIG("CCLE_rna_celline.csv")
  refmat[["Cancer type (CCLE)"]] <- readSIG("CCLE_rna_cancertype.csv")

  ## list of methods to compute
  methods <- c("DCQ", "DeconRNAseq", "I-NNLS", "NNLM", "cor", "CIBERSORT", "EPIC")

  if (full == FALSE) {
    ## Fast methods, subset of references
    sel <- c(
      "Immune cell (LM22)", "Immune cell (ImmunoStates)",
      "Immune cell (DICE)", "Immune cell (ImmProt)",
      "Tissue (GTEx)", "Cell line (HPA)", "Cancer type (CCLE)"
    )
    refmat <- refmat[intersect(sel, names(refmat))]
    methods <- c("DCQ", "DeconRNAseq", "I-NNLS", "NNLM", "cor")
  }

  ## counts <- pgx$counts
  counts <- rna.counts
  rownames(counts) <- toupper(pgx$genes[rownames(counts), "gene_name"])
  res <- pgx_multiple_deconvolution(counts, refmat = refmat, method = methods)

  pgx$deconv <- res$results
  if (!is.null(res$timings)) {
    rownames(res$timings) <- paste0("[deconvolution]", rownames(res$timings))
    res$timings
    pgx$timings <- rbind(pgx$timings, res$timings)
  }
  remove(refmat)
  remove(res)

  return(pgx)
}

# compute.cellcycle.gender
compute_cellcycle_gender <- function(pgx, rna.counts = pgx$counts) {
  pp <- rownames(rna.counts)
  is.mouse <- (mean(grepl("[a-z]", gsub(".*:|.*\\]", "", pp))) > 0.8)
  if (!is.mouse) {
    pgx$samples$cell.cycle <- NULL
    pgx$samples$.cell.cycle <- NULL
    counts <- rna.counts
    rownames(counts) <- toupper(pgx$genes[rownames(counts), "gene_name"])
    res <- try(infer_cell_cycle_phase(counts)) ## can give bins error
    if (class(res) != "try-error") {
      pgx$samples$.cell_cycle <- res
    }
    if (!(".gender" %in% colnames(pgx$samples))) {
      pgx$samples$.gender <- NULL
      X <- log2(1 + rna.counts)
      gene_name <- pgx$genes[rownames(X), "gene_name"]
      pgx$samples$.gender <- infer_gender(X, gene_name)
    }
    Matrix::head(pgx$samples)
  }
  return(pgx)
}

# compute.drugActivityEnrichment
compute_drug_activity_enrichment <- function(pgx, cmap.dir) {
  ## -------------- drug enrichment
  ref.db <- dir(cmap.dir, pattern = "^L1000-.*rds$")
  if (length(ref.db) == 0) {
    warning("Missing drug activity database")
    return(pgx)
  }
  names(ref.db) <- sub("-", "/", gsub("_.*", "", ref.db))
  f <- ref.db[1]

  for (i in 1:length(ref.db)) {
    f <- ref.db[i]
    message("[compute.drugActivityEnrichment] reading L1000 reference: ", f)
    X <- readRDS(file = file.path(cmap.dir, f))
    xdrugs <- gsub("[_@].*$", "", colnames(X))
    ndrugs <- length(table(xdrugs))
    message("number of profiles: ", ncol(X))
    message("number of drugs: ", ndrugs)
    is.drug <- grepl("activity|drug|ChemPert", f, ignore.case = TRUE)

    NPRUNE <- 250
    fname <- names(ref.db)[i]
    out1 <- compute_drug_enrichment(
      pgx, X, xdrugs,
      methods = c("GSEA", "cor"),
      nmin = 3, nprune = NPRUNE, contrast = NULL
    )

    if (is.null(out1)) {
      warning("compute_drug_enrichment failed!")
      next()
    }

    ## --------------- attach annotation
    annot0 <- NULL
    if (is.drug) {
      annot0 <- read.csv(file.path(cmap.dir, "L1000_repurposing_drugs.txt"),
        sep = "\t", comment.char = "#"
      )
      annot0$drug <- annot0$pert_iname
      rownames(annot0) <- annot0$pert_iname
    } else {
      ## gene perturbation OE/LIG/SH
      dd <- rownames(out1[["GSEA"]]$X)
      d1 <- dd
      d2 <- sub("-.*", "", dd)
      annot0 <- data.frame(drug = dd, moa = d1, target = d2)
      rownames(annot0) <- dd
    }

    annot0 <- annot0[match(rownames(out1[["GSEA"]]$X), rownames(annot0)), ]
    rownames(annot0) <- rownames(out1[["GSEA"]]$X)

    ## --------------- attach results to object
    db <- names(ref.db)[i]
    pgx$drugs[[db]] <- out1[["GSEA"]]
    pgx$drugs[[db]][["annot"]] <- annot0[, c("drug", "moa", "target")]
    pgx$drugs[[db]][["clust"]] <- out1[["clust"]]
    pgx$drugs[[db]][["stats"]] <- out1[["stats"]]
  }

  remove(X)
  remove(xdrugs)
  return(pgx)
}

# compute.drugSensitivityEnrichment
compute_drug_sensitivity_enrichment <- function(pgx, cmap.dir) {
  ref.db <- dir(cmap.dir, pattern = "sensitivity.*rds$")
  if (length(ref.db) == 0) {
    warning("Missing drug sensitivity database")
    return(pgx)
  }
  names(ref.db) <- sub("-", "/", gsub("_.*", "", ref.db))
  ref <- ref.db[1]
  for (i in 1:length(ref.db)) {
    ref <- ref.db[i]
    X <- readRDS(file = file.path(cmap.dir, ref))
    xdrugs <- gsub("[@_].*$", "", colnames(X))

    NPRUNE <- 250
    out1 <- compute_drug_enrichment(
      pgx, X, xdrugs,
      methods = c("GSEA", "cor"),
      nmin = 10, nprune = NPRUNE, contrast = NULL
    )

    if (!is.null(out1)) {
      ## attach annotation
      db <- sub("-.*", "", ref)
      annot0 <- read.csv(file.path(cmap.dir, paste0(db, "-drugs.csv")))
      rownames(annot0) <- annot0$drug
      annot0 <- annot0[match(rownames(out1[["GSEA"]]$X), rownames(annot0)), ]
      rownames(annot0) <- rownames(out1[["GSEA"]]$X)

      s1 <- names(ref.db)[i]
      pgx$drugs[[s1]] <- out1[["GSEA"]]
      pgx$drugs[[s1]][["annot"]] <- annot0[, c("moa", "target")]
      pgx$drugs[[s1]][["clust"]] <- out1[["clust"]]
      pgx$drugs[[s1]][["stats"]] <- out1[["stats"]]
    }
  } ## end of for rr

  remove(X)
  remove(xdrugs)
  return(pgx)
}

# compute.omicsGraphs
compute_omics_graphs <- function(pgx) {
  pgx$omicsnet <- create_omics_graph(pgx)
  pgx$pathscores <- compute_path_scores(pgx$omicsnet, strict.pos = FALSE)

  ## compute reduced graph
  pgx$omicsnet.reduced <- reduce_omics_graph(pgx)
  pgx$pathscores.reduced <- compute_path_scores(pgx$omicsnet.reduced, strict.pos = FALSE)
  return(pgx)
}
