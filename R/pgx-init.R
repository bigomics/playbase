##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

#' @title Initialize and validate a PGX object
#'
#' @description
#' Validates and initializes a PGX object by performing checks,
#' conversions, and adding default parameters.
#'
#' @param pgx A PGX object to initialize.
#'
#' @details
#' This function performs several validation checks on the PGX object:
#'
#' - Checks that counts and samples data frames are valid
#' - Checks that model parameters are defined
#' - Converts counts to matrix
#' - Converts model matrix to matrix
#' - Defines group labels if not provided
#' - Converts contrasts to labeled matrix form
#'
#' It will throw errors if required components like counts, samples, or groups are missing.
#' Additional default parameters like model formula and contrasts are added if not already defined.
#'
#' @return
#' The initialized PGX object.
#'
#' @examples
#' \dontrun{
#' # TODO
#' }
#' @export
pgx.initialize <- function(pgx) {
  ## ---------------------------------------------------------------------
  ## This function must be called after creation of a PGX object
  ## and include some cleaning up and updating some internal
  ## structures to keep compatibility with new/old versions.
  ## ---------------------------------------------------------------------
  message("[pgx.initialize] initializing pgx object")

  ## ----------------- check object
  obj.needed <- c(
    "samples", "genes", "GMT", "gx.meta",
    "model.parameters", "tsne2d", "X"
  )
  all(obj.needed %in% names(pgx))
  if (!all(obj.needed %in% names(pgx))) {
    obj.missing <- setdiff(obj.needed, names(pgx))
    msg <- paste("invalid pgx object. missing parts in object: ", obj.missing)
    message("[pgx-init.R] *** WARNING ***", msg)
    return(NULL)
  }

  if (is.null(pgx$version)) {
    # this is needed in case the species is human, and we dont have the
    # homolog column or if we have an old pgx which will ensure consistency
    # between old and new pgx
    pgx$genes$gene_name <- as.character(pgx$genes$gene_name)
    pgx$genes$gene_title <- as.character(pgx$genes$gene_title)
    pgx$genes$human_ortholog <- toupper(as.character(pgx$genes$gene_name))
    pgx$genes$feature <- as.character(rownames(pgx$genes))
    pgx$genes$symbol <- pgx$genes$gene_name
    col_order <- c(
      "feature", "symbol", "human_ortholog",
      "gene_title", "gene_name", colnames(pgx$genes)
    )
    col_order <- col_order[!duplicated(col_order)]
    pgx$genes <- pgx$genes[, col_order, drop = FALSE]
    if (!"chr" %in% colnames(pgx$genes)) {
      pgx$genes$chr <- 0
    }
  }

  ## for COMPATIBILITY: if no counts, estimate from X
  if (is.null(pgx$counts)) {
    message("WARNING:: no counts table. estimating from X\n")
    pgx$counts <- pmax(2**pgx$X - 1, 0)
    k <- grep("lib.size|libsize", colnames(pgx$samples))[1]
    if (length(k) > 0) {
      libsize <- pgx$samples[colnames(pgx$counts), k]
      libsize <- libsize / colSums(pgx$counts, na.rm = TRUE)
      pgx$counts <- Matrix::t(Matrix::t(pgx$counts) * libsize)
    }
  }
  pgx$counts <- as.matrix(pgx$counts)
  if (!is.null(pgx$X)) pgx$X <- as.matrix(pgx$X)

  ## ----------------------------------------------------------------
  ## model parameters
  ## ----------------------------------------------------------------
  has.design <- !is.null(pgx$model.parameters$design)
  has.expmatrix <- !is.null(pgx$model.parameters$exp.matrix)
  if (!"group" %in% names(pgx$model.parameters) && has.design) {
    ii <- max.col(pgx$model.parameters$design)
    group <- colnames(pgx$model.parameters$design)[ii]
    names(group) <- rownames(pgx$model.parameters$design)
    pgx$model.parameters$group <- group
  }
  if (!"group" %in% names(pgx$model.parameters) && has.expmatrix) {
    group <- pgx.getConditions(pgx$model.parameters$exp.matrix, nmax = 0)
    names(group) <- rownames(pgx$model.parameters$exp.matrix)
    pgx$model.parameters$group <- group
  }

  if (is.null(pgx$model.parameters$group)) {
    stop("[pgx.initialize] FATAL: group is null!!!")
  }

  ## ----------------------------------------------------------------
  ## Convert to labeled contrast matrix (new style)
  ## ----------------------------------------------------------------

  ## check if 'contrasts' is present in pgx object. Older pgx might not have it.
  if (!"contrasts" %in% names(pgx) && "model.parameters" %in% names(pgx)) {
    pgx$contrasts <- contrastAsLabels(pgx$model.parameters$contr.matrix)
  }

  ## check if 'contrasts' is sample-wise and label-format (new format)
  is.numlev <- all(unique(pgx$contrasts) %in% c(NA, "", -1, 0, 1))
  is.samplewise <- all(rownames(pgx$contrasts) == rownames(pgx$samples))
  if ("contrasts" %in% names(pgx) && (!is.samplewise || is.numlev)) {
    #    design <- pgx$model.parameters$design
    #    expmat <- pgx$model.parameters$exp.matrix
    #    contr.mat <- pgx$model.parameters$contr.matrix
    new.contr <- pgx$contrasts
    is.numlev <- all(unique(new.contr) %in% c(NA, "", -1, 0, 1))
    is.numlev <- is.numlev && (-1 %in% new.contr) ## must have -1 !!
    if (is.numlev) {
      new.contr <- contrastAsLabels(new.contr)
    }
    is.groupwise <- all(rownames(new.contr) %in% pgx$samples$group)
    is.groupwise
    if (is.groupwise) {
      grp <- as.character(pgx$samples$group)
      ii <- match(grp, rownames(new.contr))
      new.contr <- new.contr[ii, , drop = FALSE]
      rownames(new.contr) <- rownames(pgx$samples)
    }
    pgx$contrasts <- new.contr
  }
  # If contrasts is present and numeric, got to run contrastAsLabels
  if (is.numeric(pgx$contrasts)) {
    pgx$contrasts <- contrastAsLabels(pgx$contrasts)
  }


  ## ----------------------------------------------------------------
  ## Tidy up phenotype matrix (important!!!): get numbers/integers
  ## into numeric, categorical into factors....
  ## ----------------------------------------------------------------
  pgx$samples <- utils::type.convert(pgx$samples, as.is = TRUE) ## autoconvert to datatypes
  pgx$samples <- pgx$samples[, which(colMeans(is.na(pgx$samples)) < 1), drop = FALSE]

  is.num <- sapply(pgx$samples, class) %in% c("numeric", "integer")
  numlev <- apply(pgx$samples, 2, function(x) length(unique(x[!is.na(x)])))
  is.numfac <- (is.num & numlev <= 3)

  if (any(is.numfac)) {
    for (i in which(is.numfac)) pgx$samples[, i] <- as.character(pgx$samples[, i])
  }

  ## clean up: pgx$Y is a cleaned up pgx$samples
  kk <- grep("batch|lib.size|norm.factor|repl|donor|clone|sample|barcode",
    colnames(pgx$samples),
    invert = TRUE, value = TRUE
  )
  kk <- grep("lib.size|norm.factor|donor|clone|barcode",
    colnames(pgx$samples),
    invert = TRUE, value = TRUE
  )

  ## *****************************************************************
  ## ******************NEED RETHINK***********************************
  ## *****************************************************************

  ## IK: ONLY categorical variables for the moment!!! pgx$Y is somehow
  ## DEPRECATED. Better use pgx$samples directly. Idea was to
  ## 'cleanup' the samples matrix from numerical and id columns.
  pgx$Y <- pgx$samples[colnames(pgx$X), kk, drop = FALSE]
  pgx$Y <- utils::type.convert(pgx$Y, as.is = TRUE) ## autoconvert to datatypes
  ny1 <- nrow(pgx$Y) - 1
  k1 <- pgx.getCategoricalPhenotypes(pgx$Y, min.ncat = 2, max.ncat = ny1) ## exclude
  k2 <- grep("OS.survival|cluster|condition|group", colnames(pgx$Y), value = TRUE) ## must include
  kk <- sort(unique(c(k1, k2)))
  pgx$Y <- pgx$Y[, kk, drop = FALSE]

  ## -----------------------------------------------------------------------------
  ## intersect and filter gene families (convert species to human gene sets)
  ## -----------------------------------------------------------------------------
  # Here we use the homologs when available, instead of gene_name
  genes <- ifelse(!is.na(pgx$genes$human_ortholog),
    pgx$genes$human_ortholog,
    pgx$genes$gene_name
  )

  if (is.null(pgx$organism)) {
    pgx$organism <- playbase::pgx.getOrganism(pgx)
  }

  # Check if human ortholog is empty, if it is
  # 1) run getHumanOrtholog (maybe it failed on pgx.compute bc server was unreachable)
  # 2) if still empty, grag the symbols toUpper
  if (all(is.na(pgx$genes$human_ortholog))) {
    genes_ho <- getHumanOrtholog(pgx$organism, pgx$genes$symbol)$human
    if (all(is.na(genes_ho))) {
      pgx$genes$human_ortholog <- toupper(pgx$genes$symbol)
    } else {
      pgx$genes$human_ortholog <- genes_ho
    }
  }

  if (pgx$organism %in% c("Human", "human") | !is.null(pgx$version)) {
    pgx$families <- lapply(playdata::FAMILIES, function(x) {
      intersect(x, genes)
    })
  } else {
    pgx$families <- lapply(playdata::FAMILIES, function(x, genes, annot_table) {
      x <- intersect(x, genes)
      x <- annot_table$symbol[match(x, annot_table$human_ortholog)]
      return(x)
    }, genes = genes, annot_table = pgx$genes)
  }
  famsize <- sapply(pgx$families, length)
  pgx$families <- pgx$families[which(famsize >= 10)]

  all.genes <- sort(unique(pgx$genes$symbol))
  pgx$families[["<all>"]] <- all.genes

  ## -----------------------------------------------------------------------------
  ## Recompute geneset meta.fx as average fold-change of genes
  ## -----------------------------------------------------------------------------

  if (pgx$organism != "No organism" || nrow(pgx$GMT) > 0) {
    message("[pgx.initialize] Recomputing geneset fold-changes")
    nc <- length(pgx$gset.meta$meta)
    for (i in 1:nc) {
      gs <- pgx$gset.meta$meta[[i]]
      fc <- pgx$gx.meta$meta[[i]]$meta.fx
      names(fc) <- rownames(pgx$gx.meta$meta[[i]])
      # If use does not collapse by gene
      if (!all(names(fc) %in% pgx$genes$symbol)) {
        names(fc) <- probe2symbol(names(fc), pgx$genes, "symbol", fill_na = TRUE)
        fc <- fc[names(fc) != ""]
        fc <- fc[!is.na(names(fc))]
      }
      gg <- intersect(names(fc), rownames(pgx$GMT))
      G1 <- Matrix::t(pgx$GMT[gg, rownames(gs)])
      mx <- (G1 %*% fc[gg])[, 1]
      pgx$gset.meta$meta[[i]]$meta.fx <- mx
    }
  } else {
    message("[pgx.initialize] No genematrix found")
  }

  ## -----------------------------------------------------------------------------
  ## Recode survival
  ## -----------------------------------------------------------------------------
  pheno <- colnames(pgx$Y)
  ## DLBCL coding
  if (("OS.years" %in% pheno && "OS.status" %in% pheno)) {
    message("found OS survival data")
    event <- (pgx$Y$OS.status %in% c("DECEASED", "DEAD", "1", "yes", "YES", "dead"))
    pgx$Y$OS.survival <- ifelse(event, pgx$Y$OS.years, -pgx$Y$OS.years)
  }

  ## cBioportal coding
  if (("OS_MONTHS" %in% pheno && "OS_STATUS" %in% pheno)) {
    message("[pgx.initialize] found OS survival data\n")
    event <- (pgx$Y$OS_STATUS %in% c("DECEASED", "DEAD", "1", "yes", "YES", "dead"))
    pgx$Y$OS.survival <- ifelse(event, pgx$Y$OS_MONTHS, -pgx$Y$OS_MONTHS)
  }

  ## -----------------------------------------------------------------------------
  ## Check if clustering is done
  ## -----------------------------------------------------------------------------
  message("[pgx.initialize] Check if clustering is done...")
  if (!"cluster.genes" %in% names(pgx)) {
    message("[pgx.initialize] clustering genes...")
    pgx <- pgx.clusterGenes(pgx, methods = "umap", dims = c(2), level = "gene")
    pgx$cluster.genes$pos <- lapply(pgx$cluster.genes$pos, pos.compact)
  }
  if (!"cluster.gsets" %in% names(pgx) && length(pgx$gsetX) > 0) {
    message("[pgx.initialize] clustering genesets...")
    pgx <- pgx.clusterGenes(pgx, methods = "umap", dims = c(2), level = "geneset")
    pgx$cluster.gsets$pos <- lapply(pgx$cluster.gsets$pos, pos.compact)
  }
  if (!"cluster" %in% names(pgx)) {
    message("[pgx.initialize] clustering samples...")
    pgx <- pgx.clusterSamples(
      pgx,
      dims = c(2, 3),
      perplexity = NULL,
      methods = c("pca", "tsne", "umap")
    )
  }

  ## -----------------------------------------------------------------------------
  ## Remove redundant???
  ## -----------------------------------------------------------------------------
  message("[pgx.initialize] Remove redundant phenotypes...")
  if (".gender" %in% colnames(pgx$Y) &&
    any(c("gender", "sex") %in% tolower(colnames(pgx$Y)))) {
    pgx$Y$.gender <- NULL
  }

  ## -----------------------------------------------------------------------------
  ## Keep compatible with OLD formats
  ## -----------------------------------------------------------------------------
  message("[pgx.initialize] Keep compatible OLD formats...")
  if (any(c("mono", "combo") %in% names(pgx$drugs))) {
    dd <- pgx$drugs[["mono"]]
    aa1 <- pgx$drugs[["annot"]]
    if (is.null(aa1)) {
      aa1 <- playdata::L1000_REPURPOSING_DRUGS
      aa1$drug <- aa1$pert_iname
      rownames(aa1) <- aa1$pert_iname
    }
    dd[["annot"]] <- aa1
    pgx$drugs[["activity/L1000"]] <- dd
    if ("combo" %in% names(pgx$drugs)) {
      dd2 <- pgx$drugs[["combo"]]
      combo <- rownames(dd2$X)
      aa2 <- pgx.createComboDrugAnnot(combo, aa1)
      dd2[["annot"]] <- aa2
      pgx$drugs[["activity-combo/L1000"]] <- dd2
    }
    pgx$drugs$mono <- NULL
    pgx$drugs$annot <- NULL
    pgx$drugs$combo <- NULL
  }

  ## -----------------------------------------------------------------------------
  ## remove large deprecated outputs from objects
  ## -----------------------------------------------------------------------------
  message("[pgx.initialize] Removing deprecated objects...")
  pgx$gx.meta$outputs <- NULL
  pgx$gset.meta$outputs <- NULL
  pgx$gmt.all <- NULL

  message("[pgx.initialize] done!")
  return(pgx)
}
