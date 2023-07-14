##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

#' Compute Extra Analysis
#'
#' This function computes additional analysis based on the input data, such as GO core graph,
#' deconvolution, phenotype inference, drug activity enrichment, OmicsGraphs, WordCloud statistics,
#' connectivity scores, and wgcna.
#'
#' @param ngs An object containing the input data for analysis.
#' @param extra A character vector specifying which additional analyses to perform.
#' @param sigdb A character vector specifying the path to the sigdb-.h5 files for connectivity scores.
#' @param libx.dir The directory where the sigdb-.h5 files are located.
#' @return An updated object with additional analysis results.
#' @export
compute_extra <- function(ngs, extra = c(
                            "meta.go", "deconv", "infer", "drugs", ## "graph",
                            "connectivity", "wordcloud", "wgcna"
                          ), sigdb = NULL, libx.dir = NULL) {
  timings <- c()

  if (length(extra) == 0) {
    return(ngs)
  }

  ## detect if it is single or multi-omics
  single.omics <- !any(grepl("\\[", rownames(ngs$counts)))
  single.omics
  if (single.omics) {
    message(">>> computing extra for SINGLE-OMICS")
    rna.counts <- ngs$counts
  } else {
    message(">>> computing extra for MULTI-OMICS")
    data.type <- gsub("\\[|\\].*", "", rownames(ngs$counts))
    jj <- which(data.type %in% c("gx", "mrna"))
    length(jj)
    if (length(jj) == 0) {
      stop("FATAL. could not find gx/mrna values.")
    }
    rna.counts <- ngs$counts[jj, ]
    is.logged <- (min(rna.counts, na.rm = TRUE) < 0 ||
      max(rna.counts, na.rm = TRUE) < 50)
    if (is.logged) {
      message("expression data seems log. undoing logarithm")
      rna.counts <- 2**rna.counts
    }
  }

  if ("meta.go" %in% extra) {
    message(">>> Computing GO core graph...")
    tt <- system.time({
      ngs$meta.go <- pgx.computeCoreGOgraph(ngs, fdr = 0.20)
    })
    timings <- rbind(timings, c("meta.go", tt))
    message("<<< done!")
  }

  if ("deconv" %in% extra) {
    message(">>> computing deconvolution")
    tt <- system.time({
      ngs <- compute_deconvolution(
        ngs,
        rna.counts = rna.counts,
        full = FALSE
      )
    })
    timings <- rbind(timings, c("deconv", tt))
    message("<<< done!")
  }

  if ("infer" %in% extra) {
    message(">>> inferring extra phenotypes...")
    tt <- system.time({
      ngs <- compute_cellcycle_gender(ngs, rna.counts = rna.counts)
    })
    timings <- rbind(timings, c("infer", tt))
    message("<<< done!")
  }

  if ("drugs" %in% extra) {
    ngs$drugs <- NULL ## reset??

    message(">>> Computing drug activity enrichment...")
    tt <- system.time({
      ngs <- compute_drugActivityEnrichment(ngs, libx.dir = libx.dir)
    })
    timings <- rbind(timings, c("drugs", tt))

    if (!is.null(libx.dir)) {
      message(">>> Computing drug sensitivity enrichment...")
      tt <- system.time({
        ngs <- compute_drugSensitivityEnrichment(ngs, libx.dir)
      })
      timings <- rbind(timings, c("drugs-sx", tt))
    } else {
      message(">>> Skipping drug sensitivity enrichment (no libx.dir)...")
    }
    message("<<< done!")
  }

  if ("graph" %in% extra) {
    message(">>> computing OmicsGraphs...")
    tt <- system.time({
      ngs <- compute_omicsGraphs(ngs)
    })
    timings <- rbind(timings, c("graph", tt))
    message("<<< done!")
  }

  if ("wordcloud" %in% extra) {
    message(">>> computing WordCloud statistics...")
    tt <- system.time({
      res <- pgx.calculateWordCloud(ngs, progress = NULL, pg.unit = 1)
    })
    timings <- rbind(timings, c("wordcloud", tt))
    ngs$wordcloud <- res
    remove(res)
    message("<<< done!")
  }

  # I THINK THIS REQUIRES libx.dir TO BE SET TO FIND sigdb-.h5 FILES (-Nick)
  if ("connectivity" %in% extra) {
    # try to find sigdb in libx dir if not specified
    if (!is.null(libx.dir) || !is.null(sigdb)) {
      message(">>> Computing connectivity scores...")
      if (is.null(sigdb)) {
        sigdb <- dir(file.path(libx.dir, "sigdb"), pattern = "^sigdb-.*h5$", full.names = TRUE)
      }

      db <- sigdb[1]
      for (db in sigdb) {
        if (file.exists(db)) {
          ntop <- 10000
          ntop <- 1000
          message("computing connectivity scores for ", db)
          ## in memory for many comparisons
          meta <- pgx.getMetaFoldChangeMatrix(ngs, what = "meta")
          inmemory <- ifelse(ncol(meta$fc) > 50, TRUE, FALSE) ## NEED RETHINK!! reverse?
          inmemory
          tt <- system.time({
            scores <- pgx.computeConnectivityScores(
              ngs, db,
              ntop = ntop, contrasts = NULL,
              remove.le = TRUE, inmemory = inmemory
            )
          })
          timings <- rbind(timings, c("connectivity", tt))

          db0 <- sub(".*/", "", db)
          ngs$connectivity[[db0]] <- scores
          remove(scores)
        }
      }
    } else {
      message(">>> Skipping connectivity scores (no libx.dir)...")
    }
  }


  if ("wgcna" %in% extra) {
    message(">>> Computing wgcna...")
    ngs$wgcna <- playbase::pgx.wgcna(ngs)
  }

  ## ------------------------------------------------------
  ## pretty collapse all timings
  ## ------------------------------------------------------
  #
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

  ngs$timings <- rbind(ngs$timings, timings0)
  message("<<< done!")

  return(ngs)
}

## -------------- deconvolution analysis --------------------------------

#' Compute Deconvolution
#'
#' This function performs deconvolution analysis on the input RNA expression data using various reference matrices
#' and methods. It estimates the abundance of different cell types or tissue components present in the data.
#'
#' @param ngs An object containing the input data for analysis.
#' @param rna.counts A matrix or data frame of RNA expression counts. Defaults to the counts in the input object.
#' @param full A logical value indicating whether to use the full set of reference matrices and methods (TRUE),
#'   or a subset of faster methods and references (FALSE).
#' @return An updated object with deconvolution results.
#' @export
compute_deconvolution <- function(ngs, rna.counts = ngs$counts, full = FALSE) {
  ## list of reference matrices
  refmat <- list()
  refmat[["Immune cell (LM22)"]] <- playdata::LM22
  refmat[["Immune cell (ImmProt)"]] <- playdata::IMMPROT_SIGNATURE1000
  refmat[["Immune cell (DICE)"]] <- playdata::DICE_SIGNATURE1000
  refmat[["Immune cell (ImmunoStates)"]] <- playdata::IMMUNOSTATES_MATRIX
  refmat[["Tissue (HPA)"]] <- playdata::RNA_TISSUE_MATRIX
  refmat[["Tissue (GTEx)"]] <- playdata::GTEX_RNA_TISSUE_TPM
  refmat[["Cell line (HPA)"]] <- playdata::HPA_RNA_CELLINE
  refmat[["Cell line (CCLE)"]] <- playdata::CCLE_RNA_CELLINE
  refmat[["Cancer type (CCLE)"]] <- playdata::CCLE_RNA_CANCERTYPE

  ## list of methods to compute
  methods <- c("DCQ", "DeconRNAseq", "I-NNLS", "NNLM", "cor", "CIBERSORT", "EPIC")

  if (full == FALSE) {
    ## Faster methods, subset of references
    sel <- c(
      "Immune cell (LM22)", "Immune cell (ImmunoStates)",
      "Immune cell (DICE)", "Immune cell (ImmProt)",
      "Tissue (GTEx)", "Cell line (HPA)", "Cancer type (CCLE)"
    )
    refmat <- refmat[intersect(sel, names(refmat))]
    methods <- c("DCQ", "DeconRNAseq", "I-NNLS", "NNLM", "cor")
  }

  counts <- rna.counts
  rownames(counts) <- toupper(ngs$genes[rownames(counts), "gene_name"])
  res <- pgx.multipleDeconvolution(counts, refmat = refmat, method = methods)

  ngs$deconv <- res$results
  if (!is.null(res$timings)) {
    rownames(res$timings) <- paste0("[deconvolution]", rownames(res$timings))
    res$timings
    ngs$timings <- rbind(ngs$timings, res$timings)
  }
  remove(refmat)
  remove(res)

  return(ngs)
}

## -------------- infer sample characteristics --------------------------------

#' Compute Cell Cycle and Gender Inference
#'
#' This function performs cell cycle phase inference and gender estimation based on the input RNA expression data.
#' The cell cycle phase inference is performed using the Seurat package.
#'
#' @param ngs An object containing the input data for analysis.
#' @param rna.counts A matrix or data frame of RNA expression counts.
#'   Defaults to the counts in the input object.
#' @return An updated object with cell cycle and gender inference results.
#' @export
compute_cellcycle_gender <- function(ngs, rna.counts = ngs$counts) {
  pp <- rownames(rna.counts)
  is.mouse <- (mean(grepl("[a-z]", gsub(".*:|.*\\]", "", pp))) > 0.8)
  is.mouse
  if (!is.mouse) {
    message("estimating cell cycle (using Seurat)...")
    ngs$samples$cell.cycle <- NULL
    ngs$samples$.cell.cycle <- NULL

    counts <- rna.counts
    rownames(counts) <- toupper(ngs$genes[rownames(counts), "gene_name"])
    res <- try(pgx.inferCellCyclePhase(counts)) ## can give bins error
    if (class(res) != "try-error") {
      ngs$samples$.cell_cycle <- res
      table(ngs$samples$.cell_cycle)
    }
    if (!(".gender" %in% colnames(ngs$samples))) {
      message("estimating gender...")
      ngs$samples$.gender <- NULL
      X <- log2(1 + rna.counts)
      gene_name <- ngs$genes[rownames(X), "gene_name"]
      ngs$samples$.gender <- pgx.inferGender(X, gene_name)
      table(ngs$samples$.gender)
    } else {
      message("gender already estimated. skipping...")
    }
    Matrix::head(ngs$samples)
  }
  return(ngs)
}


#' Compute Drug Activity Enrichment
#'
#' This function performs drug activity enrichment analysis based on the input RNA expression data.
#' It uses drug activity databases to compute enrichment scores and attach the results to the input object.
#' The drug activity databases include L1000_ACTIVITYS_N20D1011 and L1000_GENE_PERTURBATION.
#'
#' @param ngs An object containing the input data for analysis.
#' @param libx.dir The directory path where the drug activity databases are located.
#'   This is required if calling the function compute_full_drugActivityEnrichment.
#' @return An updated object with drug activity enrichment results.
#' @export
compute_drugActivityEnrichment <- function(ngs, libx.dir = NULL) {
  ## -------------- drug enrichment
  # get drug activity databases

  if (is.null(libx.dir)) message("WARNING Need libx.dir for full drugActivityEnrichment")

  if (!is.null(libx.dir) && dir.exists(libx.dir)) {
    ## scan for extra connectivity reference files in libx
    cmap.dir <- file.path(libx.dir, "cmap")
    db.files <- dir(cmap.dir, pattern = "L1000-activity.*rds$|L1000-gene.*rds$")
    db.files
    ref.db <- lapply(db.files, function(f) readRDS(file.path(cmap.dir, f)))
    names(ref.db) <- sub("-", "/", gsub("_.*", "", db.files))
  } else {
    ## get the default 'light' version of the drug CMAP
    ref.db <- list(
      "L1000/activity" = playdata::L1000_ACTIVITYS_N20D1011
    )
  }

  for (i in seq_along(ref.db)) {
    f <- names(ref.db)[i]
    message("[compute_drugActivityEnrichment] computing activity CMAP for ", f)

    X <- ref.db[[i]]
    xdrugs <- gsub("[_@].*$", "", colnames(X))
    ndrugs <- length(table(xdrugs))
    is.drug <- grepl("activity|drug|ChemPert", f, ignore.case = TRUE)

    out1 <- playbase::pgx.computeDrugEnrichment(
      obj = ngs,
      X = X,
      xdrugs = xdrugs,
      methods = c("GSEA", "cor"),
      nmin = 10,
      nprune = 1000,
      contrast = NULL
    )

    if (is.null(out1)) {
      message("[compute_drugActivityEnrichment] WARNING:: pgx.computeDrugEnrichment failed!")
      next()
    }

    ## --------------- attach annotation
    annot0 <- NULL
    if (is.drug) {
      annot0 <- playdata::L1000_REPURPOSING_DRUGS
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
    Matrix::head(annot0)

    ## --------------- attach results to object
    db <- names(ref.db)[i]
    ngs$drugs[[db]] <- out1[["GSEA"]]
    ngs$drugs[[db]][["annot"]] <- annot0[, c("drug", "moa", "target")]
    ngs$drugs[[db]][["clust"]] <- out1[["clust"]]
    ngs$drugs[[db]][["stats"]] <- out1[["stats"]]
  }

  return(ngs)
}


#' Compute Drug Sensitivity Enrichment
#'
#' This function performs drug sensitivity enrichment analysis based on the input RNA expression data.
#' It uses drug sensitivity databases to compute enrichment scores and attach the results to the input object.
#' The drug sensitivity databases are located in the specified directory (libx.dir).
#'
#' @param ngs An object containing the input data for analysis.
#' @param libx.dir The directory path where the drug sensitivity databases are located.
#' @return An updated object with drug sensitivity enrichment results.
#' @export
compute_drugSensitivityEnrichment <- function(ngs, libx.dir = NULL) {
  if (is.null(libx.dir) || !dir.exists(libx.dir)) {
    return(ngs)
  }

  cmap.dir <- file.path(libx.dir, "cmap")
  ref.db <- dir(cmap.dir, pattern = "sensitivity.*rds$")
  if (length(ref.db) == 0) {
    message("[compute_drugSensitivityEnrichment] Warning:: missing drug sensitivity database")
    return(ngs)
  }
  names(ref.db) <- sub("-", "/", gsub("_.*", "", ref.db))
  ref.db
  ref <- ref.db[1]
  for (i in seq_along(ref.db)) {
    ref <- ref.db[i]
    message("[compute_drugSensitivityEnrichment] computing sensitivity CMAP for ", ref)
    X <- readRDS(file = file.path(cmap.dir, ref))
    xdrugs <- gsub("[@_].*$", "", colnames(X))
    length(table(xdrugs))
    dim(X)

    out1 <- playbase::pgx.computeDrugEnrichment(
      ngs, X, xdrugs,
      methods = c("GSEA", "cor"),
      nmin = 10,
      nprune = 1000,
      contrast = NULL
    )

    if (!is.null(out1)) {
      ## attach annotation
      db <- sub("-.*", "", ref)
      annot0 <- read.csv(file.path(cmap.dir, paste0(db, "-drugs.csv")))
      Matrix::head(annot0)
      rownames(annot0) <- annot0$drug
      annot0 <- annot0[match(rownames(out1[["GSEA"]]$X), rownames(annot0)), ]
      rownames(annot0) <- rownames(out1[["GSEA"]]$X)
      dim(annot0)

      s1 <- names(ref.db)[i]
      ngs$drugs[[s1]] <- out1[["GSEA"]]
      ngs$drugs[[s1]][["annot"]] <- annot0[, c("moa", "target")]
      ngs$drugs[[s1]][["clust"]] <- out1[["clust"]]
      ngs$drugs[[s1]][["stats"]] <- out1[["stats"]]
    }
  } ## end of for rr

  return(ngs)
}

## ------------------ Omics graphs --------------------------------

#' Compute Omics Graphs
#'
#' This function computes omics graphs and path scores based on the input data.
#' It creates the omics graph using the input object and calculates path scores using the omics graph.
#' It also computes a reduced graph and path scores based on the reduced graph.
#'
#' @param ngs An object containing the input data for analysis.
#' @return An updated object with omics graphs and path scores.
#' @export
compute_omicsGraphs <- function(ngs) {
  ngs$omicsnet <- pgx.createOmicsGraph(ngs)
  ngs$pathscores <- pgx.computePathscores(ngs$omicsnet, strict.pos = FALSE)

  ## compute reduced graph
  ngs$omicsnet.reduced <- pgx.reduceOmicsGraph(ngs)
  ngs$pathscores.reduced <- pgx.computePathscores(ngs$omicsnet.reduced, strict.pos = FALSE)
  return(ngs)
}
