##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

#' Title
#'
#' @param ngs value
#' @param extra value
#' @param sigdb value
#'
#' @return
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
    ## rownames(rna.counts) <- gsub(".*:|.*\\]","",rownames(rna.counts))
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
      ngs <- compute_drugActivityEnrichment(ngs)
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
          inmemory <- ifelse(ncol(meta$fc) > 50, TRUE, FALSE)
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
    ngs$wgcna <- playbase::pgx.wgcna(pgx)
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

  ngs$timings <- rbind(ngs$timings, timings0)
  message("<<< done!")

  return(ngs)
}

## -------------- deconvolution analysis --------------------------------

#' Title
#'
#' @param ngs value
#' @param rna.counts value
#' @param full value
#'
#' @return
#' @export
#'
#' @examples
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
  ## methods = DECONV.METHODS
  methods <- c("DCQ", "DeconRNAseq", "I-NNLS", "NNLM", "cor", "CIBERSORT", "EPIC")
  ## methods <- c("NNLM","cor")

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

  ## counts <- ngs$counts
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

#' Title
#'
#' @param ngs value
#' @param rna.counts value
#'
#' @return
#' @export
#'
#' @examples
compute_cellcycle_gender <- function(ngs, rna.counts = ngs$counts) {
  pp <- rownames(rna.counts)
  is.mouse <- (mean(grepl("[a-z]", gsub(".*:|.*\\]", "", pp))) > 0.8)
  is.mouse
  if (!is.mouse) {
    if (1) {
      message("estimating cell cycle (using Seurat)...")
      ngs$samples$cell.cycle <- NULL
      ngs$samples$.cell.cycle <- NULL
      ## counts <- ngs$counts
      counts <- rna.counts
      rownames(counts) <- toupper(ngs$genes[rownames(counts), "gene_name"])
      res <- try(pgx.inferCellCyclePhase(counts)) ## can give bins error
      if (class(res) != "try-error") {
        ngs$samples$.cell_cycle <- res
        table(ngs$samples$.cell_cycle)
      }
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


#' Title
#'
#' @param ngs value
#'
#' @return
#' @export
#'
#' @examples
compute_drugActivityEnrichment <- function(ngs, libx.dir = NULL) {
  ## -------------- drug enrichment
  # get drug activity databases

  #ngs=pgx
  #libx.dir = params$libx.dir
  
  if (is.null(libx.dir)) message('WARNING Need libx.dir if you call compute_full_drugActivityEnrichment')
  
  cmap.dir <- file.path(libx.dir, "cmap")
  file.gene.db <- dir(cmap.dir, pattern = "n8m20g5812.*rds$")

  if (length(gene.db)>1) message('WARNING multiple gene.db files found. Using first one.')

  gene.db <- readRDS(file.path(cmap.dir,file.gene.db[1]))

  View(gene.db)

  

  if(file.exists(gene.db)){
    ref.db <- list(
    "L1000_ACTIVITYS_N20D1011" = playdata::L1000_ACTIVITYS_N20D1011,
    "L1000_GENE_PERTURBATION" = gene.db
  )

  }
  

  for (i in 1:length(ref.db)) {
    f <- names(ref.db)[i]
    message("[compute_drugActivityEnrichment] reading L1000 reference: ", f)
    X <- ref.db[[i]]
    xdrugs <- gsub("[_@].*$", "", colnames(X))
    ndrugs <- length(table(xdrugs))
    message("number of profiles: ", ncol(X))
    message("number of drugs: ", ndrugs)
    is.drug <- grepl("activity|drug|ChemPert", f, ignore.case = TRUE)

    NPRUNE <- 250
    fname <- names(ref.db)[i]
    out1 <- pgx.computeDrugEnrichment(
      ngs, X, xdrugs,
      methods = c("GSEA", "cor"),
      nmin = 3, nprune = NPRUNE, contrast = NULL
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

    if (1) {
      annot0 <- annot0[match(rownames(out1[["GSEA"]]$X), rownames(annot0)), ]
      rownames(annot0) <- rownames(out1[["GSEA"]]$X)
      Matrix::head(annot0)
    }

    ## --------------- attach results to object
    db <- names(ref.db)[i]
    ## ngs$drugs[["activity/L1000"]]  <- res.mono[["GSEA"]]
    ## ngs$drugs[["activity/L1000"]][["annot"]] <- annot0[,c("drug","moa","target")]
    ngs$drugs[[db]] <- out1[["GSEA"]]
    ngs$drugs[[db]][["annot"]] <- annot0[, c("drug", "moa", "target")]
    ngs$drugs[[db]][["clust"]] <- out1[["clust"]]
    ngs$drugs[[db]][["stats"]] <- out1[["stats"]]
  }

  remove(X)
  remove(xdrugs)
  return(ngs)
}

#' Title
#'
#' @param ngs value
#' @param cmap.dir value
#'
#' @return
#' @export
#'
#' @examples
compute_drugSensitivityEnrichment <- function(ngs, libx.dir) {

  if (is.null(libx.dir)) message('ERROR: Need libx.dir if you call compute_drugSensitivityEnrichment')
  cmap.dir <- file.path(libx.dir, "cmap")

  ref.db <- dir(cmap.dir, pattern = "sensitivity.*rds$")
  if (length(ref.db) == 0) {
    message("[compute_drugSensitivityEnrichment] Warning:: missing drug sensitivity database")
    return(ngs)
  }
  names(ref.db) <- sub("-", "/", gsub("_.*", "", ref.db))
  ref.db
  ref <- ref.db[1]
  for (i in 1:length(ref.db)) {
    ref <- ref.db[i]
    X <- readRDS(file = file.path(cmap.dir, ref))
    xdrugs <- gsub("[@_].*$", "", colnames(X))
    length(table(xdrugs))
    dim(X)

    NPRUNE <- -1
    NPRUNE <- 250
    out1 <- pgx.computeDrugEnrichment(
      ngs, X, xdrugs,
      methods = c("GSEA", "cor"),
      nmin = 10, nprune = NPRUNE, contrast = NULL
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

  names(ngs$drugs)

  remove(X)
  remove(xdrugs)
  return(ngs)
}

## ------------------ Omics graphs --------------------------------

#' Title
#'
#' @param ngs value
#'
#' @return
#' @export
#'
#' @examples
compute_omicsGraphs <- function(ngs) {
  ## gr1$layout <- gr1$layout[igraph::V(gr1)$name,]  ## uncomment to keep entire layout
  ngs$omicsnet <- pgx.createOmicsGraph(ngs)
  ngs$pathscores <- pgx.computePathscores(ngs$omicsnet, strict.pos = FALSE)

  ## compute reduced graph
  ngs$omicsnet.reduced <- pgx.reduceOmicsGraph(ngs)
  ngs$pathscores.reduced <- pgx.computePathscores(ngs$omicsnet.reduced, strict.pos = FALSE)
  ## save(ngs, file=rda.file)
  return(ngs)
}
