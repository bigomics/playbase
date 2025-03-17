##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


#' Compute Test Genesets
#'
#' Computes test genesets for the given multi-omics data.
#'
#' @param pgx A data object representing multi-omics data.
#' @param custom.geneset Custom geneset object.
#' @param test.methods Character vector with the methods to use for testing.
#'  Currently supported: "gsva", "camera", "fgsea".
#' @param remove.outputs Logical indicating whether to remove large outputs.
#'
#' @return The updated \code{pgx} object with computed test genesets.
#'
#' @export
compute_testGenesets <- function(pgx,
                                 custom.geneset = list(gmt = NULL, info = NULL),
                                 test.methods = c("gsva", "camera", "fgsea"),
                                 remove.outputs = TRUE) {
  if (!"X" %in% names(pgx)) {
    stop("[compute_testGenesets] FATAL : object must have normalized matrix X")
  }

  if (is.null(pgx$genes$human_ortholog)) {
    # this is needed in case the species is human, and we dont have the homolog column or if we have an old pgx
    # which will ensure consistency between old and new pgx
    pgx$genes$human_ortholog <- NA
  }

  ## -----------------------------------------------------------
  ## get design and contrast matrix, and get gene list
  ## -----------------------------------------------------------

  # Get X matrix
  X <- pgx$X
  design <- pgx$model.parameters$design
  contr.matrix <- pgx$model.parameters$contr.matrix

  ALL.GSET.METHODS <- c(
    "fisher", "ssgsea", "gsva", "spearman", "camera", "fry",
    "fgsea", "gsea.permPH", "gsea.permGS", "gseaPR"
  )

  ## convert to gene list
  G <- pgx$GMT
  gmt <- playbase::mat2gmt(G)

  ## if reduced samples
  ss <- rownames(pgx$model.parameters$exp.matrix)
  if (!is.null(ss)) {
    X <- X[, ss, drop = FALSE]
  }

  ## -----------------------------------------------------------
  ## Collapse X by gene
  ## -----------------------------------------------------------
  dbg("creating temporary GENE matrix by SYMBOL...")

  X <- pgx$X
  if (!all(rownames(X) %in% pgx$genes$symbol)) {
    X <- rename_by(X, pgx$genes, "symbol", unique = TRUE)
  }

  ## -----------------------------------------------------------
  ## Run methods
  ## -----------------------------------------------------------
  message(">>> Testing gene sets with methods:", test.methods, "\n")

  if (length(intersect(rownames(X), rownames(G))) == 0) {
    message("[compute_testGenesets] ERROR : matrix X and G do not overlap")
    return(pgx)
  }

  Y <- pgx$samples
  gc()

  gset.meta <- gset.fitContrastsWithAllMethods(
    gmt = gmt,
    X = X,
    Y = Y,
    G = G,
    design = design,
    contr.matrix = contr.matrix,
    methods = test.methods,
    mc.threads = 1,
    mc.cores = NULL,
    batch.correct = TRUE
  )

  rownames(gset.meta$timings) <- paste("[test.genesets]", rownames(gset.meta$timings))
  pgx$timings <- rbind(pgx$timings, gset.meta$timings)
  pgx$gset.meta <- gset.meta

  pgx$gsetX <- pgx$gset.meta$matrices[["meta"]] ## META or average FC?
  pgx$GMT <- G[, rownames(pgx$gsetX)]


  ## -------------------------------------------------------
  ## calculate gset info and store as pgx$gset.meta
  ## -------------------------------------------------------

  gset.size <- Matrix::colSums(pgx$GMT != 0)
  gset.size.raw <- playdata::GSET_SIZE

  # combine standard genesets with custom genesets size vector
  if (!is.null(custom.geneset$gmt)) {
    gset.size.raw <- c(gset.size.raw, custom.geneset$info$GSET_SIZE)
  }

  gset.idx <- match(names(gset.size), names(gset.size.raw))
  gset.size.raw <- gset.size.raw[gset.idx]
  names(gset.size.raw) <- names(gset.size)
  gset.fraction <- gset.size / gset.size.raw

  pgx$gset.meta$info <- data.frame(
    gset.size.raw = gset.size.raw,
    gset.size = gset.size,
    gset.fraction = gset.fraction
  )

  ## -----------------------------------------------------------------------
  ## ------------------------ clean up -------------------------------------
  ## -----------------------------------------------------------------------

  ## remove large outputs... (uncomment if needed)
  if (remove.outputs) {
    pgx$gset.meta$outputs <- NULL
  }

  remove(X, Y, G, gmt)
  gc()
  return(pgx)
}



#' Clean GMT
#'
#' Cleans the GMT (Gene Matrix Transposed) object by modifying its names and removing duplicates.
#'
#' @param gmt.all A list representing the GMT object.
#' @param gmt.db A character vector representing the GMT database.
#'
#' @return The cleaned GMT object.
#'
#' @export
clean_gmt <- function(gmt.all, gmt.db) {
  gmt.db <- toupper(gmt.db)
  for (i in 1:length(gmt.all)) {
    names(gmt.all[[i]]) <- sub("\\(GO:", "(GO_", names(gmt.all[[i]]))
    names(gmt.all[[i]]) <- gsub("%", "_", names(gmt.all[[i]])) # substitute % sign in wikipathways
    names(gmt.all[[i]]) <- sub(":", "", names(gmt.all[[i]]))
    names(gmt.all[[i]]) <- paste0(toupper(gmt.db[i]), ":", names(gmt.all[[i]]))
  }
  j0 <- grep("_up", names(gmt.all))
  j1 <- grep("_down", names(gmt.all))
  for (i in j0) {
    names(gmt.all[[i]]) <- paste0(names(gmt.all[[i]]), " (up)")
  }
  for (i in j1) {
    names(gmt.all[[i]]) <- paste0(names(gmt.all[[i]]), " (down)")
  }
  names(gmt.all) <- NULL
  gmt.all <- unlist(gmt.all, recursive = FALSE, use.names = TRUE)

  ## get rid of trailing numeric values
  gmt.all <- lapply(gmt.all, function(x) gsub("[,].*", "", x))



  ## order by length and take out duplicated sets (only by name)
  gmt.all <- gmt.all[order(-sapply(gmt.all, length))]
  gmt.all <- gmt.all[!duplicated(names(gmt.all))]

  gmt.all <- gmt.all[order(names(gmt.all))]

  return(gmt.all)
}

#' Create Sparse Geneset Matrix
#'
#' Creates a sparse matrix representation of genesets from a given GMT file, filtered by size and gene frequency.
#'
#' @param gmt.all The list of genesets in GMT format.
#' @param min.geneset.size The minimum size of a geneset to be included in the matrix. Default is 15.
#' @param max.geneset.size The maximum size of a geneset to be included in the matrix. Default is 500.
#' @param min_gene_frequency The minimum frequency of a gene to be included in the matrix. Default is 10.
#'
#' @return A sparse matrix representing the filtered genesets.
#'
#' @export
createSparseGenesetMatrix <- function(
    gmt.all,
    min.geneset.size = 15,
    max.geneset.size = 500,
    min_gene_frequency = 10,
    all_genes = NULL,
    annot = NULL,
    filter_genes = TRUE) {
  # WARNING #
  # This function is used in playbase and playdata to generate curated GMT. Do not change it without testing it in both packages to ensure reproducibility.

  if (is.null(all_genes)) {
    all_genes <- unique(unlist(gmt.all))
  }
  all_genes <- sort(all_genes)

  ## ------------- filter by size
  gmt.size <- sapply(gmt.all, length)
  gmt.all <- gmt.all[which(gmt.size >= min.geneset.size & gmt.size <= max.geneset.size)]

  ## ------------- filter genes by minimum frequency and chrom
  genes.table <- table(unlist(gmt.all))
  genes <- names(which(genes.table >= min_gene_frequency))
  genes <- intersect(genes, all_genes)

  if (filter_genes == TRUE) {
    genes <- genes[grep("^LOC|RIK$", genes, invert = TRUE)]
  }

  if (!is.null(annot) && filter_genes == TRUE) {
    annot <- annot[genes, ]
    annot <- annot[annot$chr %in% c(1:22, "X", "Y"), ]
    genes <- genes[!is.na(annot$chr)]
  }

  ## Filter genesets with permitted genes (official and min.sharing)
  gmt.all <- lapply(gmt.all, function(s) intersect(s, genes))
  gmt.all <- gmt.all[which(gmt.size >= min.geneset.size & gmt.size <= max.geneset.size)] # legacy
  ## build huge sparsematrix gene x genesets
  genes <- sort(genes)
  idx.j <- lapply(gmt.all, function(s) match(s, genes))
  idx.i <- lapply(1:length(gmt.all), function(i) rep(i, length(idx.j[[i]])))

  ii <- unlist(idx.i)
  jj <- unlist(idx.j)

  G <- Matrix::sparseMatrix(
    i = ii, j = jj, x = rep(1, length(ii)),
    dims = c(length(gmt.all), length(genes))
  )
  colnames(G) <- genes
  rownames(G) <- names(gmt.all)

  # remove NA rows
  G <- G[!is.na(rownames(G)), , drop = FALSE]

  return(G)
}

