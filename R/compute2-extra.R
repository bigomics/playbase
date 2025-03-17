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
#' @param pgx An object containing the input data for analysis.
#' @param extra A character vector specifying which additional analyses to perform.
#' @param sigdb A character vector specifying the path to the sigdb-.h5 files for connectivity scores.
#' @param libx.dir The directory where the sigdb-.h5 files are located.
#' @return An updated object with additional analysis results.
#' @export
compute_extra <- function(pgx, extra = c(
                            "meta.go", "infer", "deconv", "drugs", ## "graph",
                            "connectivity", "wordcloud", "wgcna", "mofa"
                          ), sigdb = NULL, pgx.dir = "./data", libx.dir = "./libx",
                          user_input_dir = getwd()) {
  timings <- c()

  if (length(extra) == 0) {
    return(pgx)
  }
  if (!is.null(pgx.dir) && !dir.exists(pgx.dir)) pgx.dir <- NULL
  if (!is.null(libx.dir) && !dir.exists(libx.dir)) libx.dir <- NULL

  rna.counts <- pgx$counts

  # If working on non-human species, use homologs
  if (!all(is.na(pgx$genes$human_ortholog))) {
    rownames(rna.counts) <- probe2symbol(rownames(rna.counts), pgx$genes, query = "human_ortholog")
  }

  if ("meta.go" %in% extra) {
    message(">>> Computing GO core graph...")
    tt <- system.time({
      pgx$meta.go <- tryCatch(
        pgx.computeCoreGOgraph(pgx, fdr = 0.20),
        error = function(e) {
          write(as.character(e), file = paste0(user_input_dir, "/ERROR_METAGO"))
          return(NULL)
        }
      )
    })
    timings <- rbind(timings, c("meta.go", tt))
    message("<<< done!")
  }

  if ("deconv" %in% extra) {
    message(">>> computing deconvolution")
    tt <- system.time({
      pgx <- tryCatch(
        {
          compute_deconvolution(
            pgx,
            rna.counts = rna.counts,
            full = FALSE
          )
        },
        error = function(e) {
          write(as.character(e), file = paste0(user_input_dir, "/ERROR_DECONVOLUTIION"))
          return(pgx)
        }
      )
    })
    timings <- rbind(timings, c("deconv", tt))
    message("<<< done!")
  }

  if ("infer" %in% extra) {
    message(">>> inferring extra phenotypes...")
    tt <- system.time({
      pgx <- tryCatch(
        compute_cellcycle_gender(pgx, rna.counts = rna.counts),
        error = function(e) {
          write(as.character(e), file = paste0(user_input_dir, "/ERROR_INFERENCE"))
          return(pgx)
        }
      )
    })
    timings <- rbind(timings, c("infer", tt))
    message("<<< done!")
  }

  if ("drugs" %in% extra) {
    pgx$drugs <- NULL ## reset??

    message(">>> Computing drug activity enrichment...")
    tt <- system.time({
      pgx <- tryCatch(
        {
          compute_drugActivityEnrichment(pgx, libx.dir = libx.dir)
        },
        error = function(e) {
          write(as.character(e), file = paste0(user_input_dir, "/ERROR_DRUG_ACTIVITY"))
          return(pgx)
        }
      )
    })
    timings <- rbind(timings, c("drugs", tt))

    if (!is.null(libx.dir)) {
      message(">>> Computing drug sensitivity enrichment...")
      tt <- system.time({
        pgx <- tryCatch(
          {
            compute_drugSensitivityEnrichment(pgx, libx.dir)
          },
          error = function(e) {
            write(as.character(e), file = paste0(user_input_dir, "/ERROR_DRUG_SENSITIVITY"))
            return(pgx)
          }
        )
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
      pgx <- tryCatch(
        {
          compute_omicsGraphs(pgx)
        },
        error = function(e) {
          write(as.character(e), file = paste0(user_input_dir, "/ERROR_GRAPH"))
          return(pgx)
        }
      )
    })
    timings <- rbind(timings, c("graph", tt))
    message("<<< done!")
  }

  if ("wordcloud" %in% extra) {
    message(">>> computing WordCloud statistics...")
    tt <- system.time({
      res <- tryCatch(pgx.calculateWordCloud(pgx, progress = NULL, pg.unit = 1),
        error = function(e) {
          write(as.character(e), file = paste0(user_input_dir, "/ERROR_WORDCLOUD"))
          return(NULL)
        }
      )
    })
    timings <- rbind(timings, c("wordcloud", tt))
    pgx$wordcloud <- res
    remove(res)
    message("<<< done!")
  }

  # This requires libx.dir to be set (or sigdb passed) to find sigdb-.h5 files
  if ("connectivity" %in% extra) {
    # try to find sigdb in libx dir if not specified
    if (!is.null(libx.dir) || !is.null(pgx.dir) || !is.null(sigdb)) {
      message(">>> Computing connectivity scores...")
      if (is.null(sigdb)) {
        sigdb <- NULL
        if (!is.null(pgx.dir)) {
          ## make sure h5 file is up-to-date
          pgxinfo.updateDatasetFolder(pgx.dir, force = FALSE, update.sigdb = TRUE)
          user.sigdb <- file.path(pgx.dir, "datasets-sigdb.h5")
          sigdb <- c(sigdb, user.sigdb)
        }
        if (!is.null(libx.dir)) {
          libx.sigdb <- dir(file.path(libx.dir, "sigdb"),
            pattern = "h5$",
            full.names = TRUE
          )
          ## do not follow symlinks because they can just be old names/aliases
          is.symlink <- (Sys.readlink(libx.sigdb) != "")
          libx.sigdb <- libx.sigdb[!is.symlink]
          libx.sigdb
          sigdb <- c(sigdb, libx.sigdb)
        }
      }

      db <- sigdb[1]
      for (db in sigdb) {
        if (file.exists(db)) {
          message("computing connectivity scores for ", db)
          tt <- system.time({
            scores <- tryCatch(
              {
                pgx.computeConnectivityScores(
                  pgx,
                  db,
                  ntop = 200,
                  contrasts = NULL,
                  remove.le = TRUE
                )
              },
              error = function(e) {
                write(as.character(e), file = paste0(user_input_dir, "/ERROR_CONNECTIVITY"))
                return(NULL)
              }
            )
            if (!is.null(scores)) {
              db0 <- sub(".*/", "", db)
              pgx$connectivity[[db0]] <- scores
            }
            remove(scores)
          })
          timings <- rbind(timings, c("connectivity", tt))
        }
      }
    } else {
      message(">>> Skipping connectivity scores (no libx.dir)...")
    }
  }

  if ("wgcna" %in% extra) {
    message(">>> Computing wgcna...")
    tt <- system.time({
      tryCatch(
        {
          pgx$wgcna <- pgx.wgcna(pgx)
        },
        error = function(e) {
          write(as.character(e), file = paste0(user_input_dir, "/ERROR_WGCNA"))
          return(NULL)
        }
      )
    })
    timings <- rbind(timings, c("wgcna", tt))
  }

  if ("mofa" %in% extra) {
    message(">>> Computing MOFA...")
    tt <- system.time({
      tryCatch(
        {
          pgx$mofa <- pgx.compute_mofa(
            pgx,
            kernel = "MOFA",
            ntop = 2000,
            gset.ntop = 2000,
            numfactors = 10,
            add_gsets = FALSE
          )
        },
        error = function(e) {
          write(as.character(e), file = paste0(user_input_dir, "/ERROR_MOFA"))
          return(NULL)
        }
      )
    })
    timings <- rbind(timings, c("mofa", tt))
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

  pgx$timings <- rbind(pgx$timings, timings0)
  message("<<< done!")

  return(pgx)
}

## -------------- deconvolution analysis --------------------------------

#' Compute Deconvolution
#'
#' This function performs deconvolution analysis on the input RNA expression data using various reference matrices
#' and methods. It estimates the abundance of different cell types or tissue components present in the data.
#'
#' @param pgx An object containing the input data for analysis.
#' @param rna.counts A matrix or data frame of RNA expression counts. Defaults to the counts in the input object.
#' @param full A logical value indicating whether to use the full set of reference matrices and methods (TRUE),
#'   or a subset of faster methods and references (FALSE).
#' @return An updated object with deconvolution results.
#'
#' @examples
#' \dontrun{
#' pgx <- playdata::GEIGER_PGX
#' deconv <- compute_deconvolution(pgx)
#' }
#' @export
compute_deconvolution <- function(pgx, rna.counts = pgx$counts, full = FALSE) {
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
  res <- pgx.multipleDeconvolution(counts, refmat = refmat, methods = methods)

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

## -------------- infer sample characteristics --------------------------------

#' Compute Cell Cycle and Gender Inference
#'
#' This function performs cell cycle phase inference and gender estimation based on the input RNA expression data.
#' The cell cycle phase inference is performed using the Seurat package.
#'
#' @param pgx An object containing the input data for analysis.
#' @param rna.counts A matrix or data frame of RNA expression counts.
#'   Defaults to the counts in the input object.
#' @return An updated object with cell cycle and gender inference results.
#'
#' @examples
#' \dontrun{
#' pgx <- playbase::GEIGER_PGX
#' deconv <- compute_cellcycle_gender(pgx)
#' }
#' @export
compute_cellcycle_gender <- function(pgx, rna.counts = pgx$counts) {
  if (!is.null(pgx$organism)) {
    is.human <- (tolower(pgx$organism) == "human")
  } else {
    is.human <- (tolower(pgx.getOrganism(pgx)) == "human")
  }
  if (is.human) {
    message("estimating cell cycle (using Seurat)...")
    pgx$samples$cell.cycle <- NULL
    pgx$samples$.cell.cycle <- NULL

    if (!(".cell_cycle" %in% colnames(pgx$samples))) {
      counts <- rna.counts
      message("length(unique(rownames(counts))): ", length(unique(rownames(counts))))
      message("dim(counts): ", dim(counts))
      ## In multi-species now use symbol, and deduplicate in case
      ## use retains feature as "gene_name/rowname"
      rownames(counts) <- toupper(pgx$genes[rownames(counts), "symbol"])
      counts <- counts[which(!is.na(rownames(counts))), ]
      if (any(duplicated(rownames(counts)))) {
        message("Deduplicate counts for cell cycle and gender inference")
        counts <- playbase::rowmean(counts, group = rownames(counts))
        counts[which(is.nan(counts))] <- NA
      }
      res <- try(pgx.inferCellCyclePhase(counts))
      if (!inherits(res, "try-error")) {
        pgx$samples$.cell_cycle <- res
      }
    } else {
      message("cell cycle already estimated. skipping...")
    }

    if (!(".gender" %in% colnames(pgx$samples))) {
      message("estimating gender...")
      pgx$samples$.gender <- NULL
      X <- log2(1 + rna.counts)
      gene_symbol <- pgx$genes[rownames(X), "symbol"] # Use gene-symbol also for gender
      pgx$samples$.gender <- pgx.inferGender(X, gene_symbol)
    } else {
      message("gender already estimated. skipping...")
    }
  }
  return(pgx)
}


#' Compute Drug Activity Enrichment
#'
#' This function performs drug activity enrichment analysis based on the input RNA expression data.
#' It uses drug activity databases to compute enrichment scores and attach the results to the input object.
#' The drug activity databases include L1000_ACTIVITYS_N20D1011 and L1000_GENE_PERTURBATION.
#'
#' @param pgx An object containing the input data for analysis.
#' @param libx.dir The directory path where the drug activity databases are located.
#'   This is required if calling the function compute_full_drugActivityEnrichment.
#' @return An updated object with drug activity enrichment results.
#' @export
compute_drugActivityEnrichment <- function(pgx, libx.dir = NULL) {
  ## -------------- drug enrichment
  # get drug activity databases

  if (is.null(libx.dir)) message("WARNING Need libx.dir for full drugActivityEnrichment")

  if (!is.null(libx.dir) && dir.exists(libx.dir)) {
    ## scan for extra connectivity reference files in libx
    cmap.dir <- file.path(libx.dir, "cmap")
    db.files <- dir(cmap.dir, pattern = "L1000-activity.*rds$|L1000-gene.*rds$")
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
    l1000genes <- rownames(X)
    if (!pgx$organism %in% c("Human", "human")) {
      rowid <- data.table::chmatch(rownames(X), pgx$genes$human_ortholog, nomatch = NA)
      rownames(X) <- pgx$genes$human_ortholog[rowid]
      X <- X[!is.na(rowid), , drop = FALSE]
    }

    xdrugs <- gsub("[_@].*$", "", colnames(X))
    ndrugs <- length(table(xdrugs))
    is.drug <- grepl("activity|drug|ChemPert", f, ignore.case = TRUE)

    out1 <- NULL
    out1 <- pgx.computeDrugEnrichment(
      pgx = pgx,
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

    ## --------------- attach results to object
    db <- names(ref.db)[i]
    pgx$drugs[[db]] <- out1[["GSEA"]]
    pgx$drugs[[db]][["annot"]] <- annot0[, c("drug", "moa", "target")]
    pgx$drugs[[db]][["clust"]] <- out1[["clust"]]
    pgx$drugs[[db]][["stats"]] <- out1[["stats"]]
  }

  return(pgx)
}


#' Compute Drug Sensitivity Enrichment
#'
#' This function performs drug sensitivity enrichment analysis based on the input RNA expression data.
#' It uses drug sensitivity databases to compute enrichment scores and attach the results to the input object.
#' The drug sensitivity databases are located in the specified directory (libx.dir).
#'
#' @param pgx An object containing the input data for analysis.
#' @param libx.dir The directory path where the drug sensitivity databases are located.
#' @return An updated object with drug sensitivity enrichment results.
#' @export
compute_drugSensitivityEnrichment <- function(pgx, libx.dir = NULL) {
  if (is.null(libx.dir) || !dir.exists(libx.dir)) {
    return(pgx)
  }

  cmap.dir <- file.path(libx.dir, "cmap")
  ref.db <- dir(cmap.dir, pattern = "sensitivity.*rds$")
  if (length(ref.db) == 0) {
    message("[compute_drugSensitivityEnrichment] Warning:: missing drug sensitivity database")
    return(pgx)
  }
  names(ref.db) <- sub("-", "/", gsub("_.*", "", ref.db))

  ref.db
  ref <- ref.db[1]
  for (i in seq_along(ref.db)) {
    # i=2
    ref <- ref.db[i]
    message("[compute_drugSensitivityEnrichment] computing sensitivity CMAP for ", ref)
    X <- readRDS(file = file.path(cmap.dir, ref))
    xdrugs <- gsub("[@_].*$", "", colnames(X))

    out1 <- pgx.computeDrugEnrichment(
      pgx = pgx,
      X = X,
      xdrugs = xdrugs,
      methods = c("GSEA", "cor"),
      nmin = 10,
      nprune = 1000,
      contrast = NULL
    )

    if (!is.null(out1)) {
      ## attach annotation
      db <- sub("-.*", "", ref)
      annot0 <- utils::read.csv(file.path(cmap.dir, paste0(db, "-drugs.csv")))
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

  return(pgx)
}

## ------------------ Omics graphs --------------------------------

#' Compute Omics Graphs
#'
#' This function computes omics graphs and path scores based on the input data.
#' It creates the omics graph using the input object and calculates path scores using the omics graph.
#' It also computes a reduced graph and path scores based on the reduced graph.
#'
#' @param pgx An object containing the input data for analysis.
#' @return An updated object with omics graphs and path scores.
#' @export
compute_omicsGraphs <- function(pgx) {
  pgx$omicsnet <- pgx.createOmicsGraph(pgx)
  pgx$pathscores <- pgx.computePathscores(pgx$omicsnet, strict.pos = FALSE)

  ## compute reduced graph
  pgx$omicsnet.reduced <- pgx.reduceOmicsGraph(pgx)
  pgx$pathscores.reduced <- pgx.computePathscores(pgx$omicsnet.reduced, strict.pos = FALSE)
  return(pgx)
}
