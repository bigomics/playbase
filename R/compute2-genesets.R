##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

#' Normalize Matrix by Row
#'
#' Normalizes a matrix by dividing each row by the sum of its elements.
#'
#' @param G The matrix to be normalized.
#'
#' @return The normalized matrix.
#'
#' @export
normalize_matrix_by_row <- function(G) {
  # efficient normalization using linear algebra
  row_sums <- Matrix::rowSums(G)
  D <- Matrix::Diagonal(x = 1 / row_sums)
  G_scaled <- D %*% G
  rownames(G_scaled) <- rownames(G)
  return(G_scaled)
}

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
  G <- pgx$G
  gmt <- playbase::mat2gmt(G)

  ## if reduced samples
  ss <- rownames(pgx$model.parameters$exp.matrix)
  if (!is.null(ss)) {
    X <- X[, ss, drop = FALSE]
  }

  ## -----------------------------------------------------------
  ## Collapse X by gene
  ## -----------------------------------------------------------

  X <- pgx$X
  single.omics <- TRUE ## !!! for now...
  if (single.omics) {
    ## normalized matrix
    X <- X
  } else {
    data.type <- gsub("\\[|\\].*", "", rownames(pgx$counts))
    jj <- which(data.type %in% c("gx", "mrna"))
    if (length(jj) == 0) {
      stop("[pgx.add_GMT] FATAL. could not find gx/mrna values.")
    }
    X <- X[jj, ]
  }

  if (!all(rownames(X) %in% pgx$genes$symbol)) {
    X <- rename_by(X, pgx$genes, "symbol")
    X <- X[!rownames(X) == "", , drop = FALSE]
    if (any(duplicated(rownames(X)))) {
      X <- log2(rowsum(2**X, rownames(X)))
    }
  }

  ## -----------------------------------------------------------
  ## Run methods
  ## -----------------------------------------------------------
  cat(">>> Testing gene sets with methods:", test.methods, "\n")

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

  all_genes <- sort(all_genes)

  ## ------------- filter by size
  gmt.size <- sapply(gmt.all, length)

  gmt.all <- gmt.all[which(gmt.size >= min.geneset.size & gmt.size <= max.geneset.size)]

  ## ------------- filter genes by minimum frequency and chrom
  genes.table <- table(unlist(gmt.all))
  genes <- names(which(genes.table >= min_gene_frequency))

  annot <- annot

  if (filter_genes == TRUE) {
    genes <- genes[grep("^LOC|RIK$", genes, invert = TRUE)]
    genes <- intersect(genes, all_genes)
  }

  if (filter_genes == TRUE) {
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


#' Merge Sparse Matrix
#'
#' Merges two sparse matrices by combining their columns and ensuring the same gene order.
#'
#' @param m1 The first sparse matrix.
#' @param m2 The second sparse matrix.
#'
#' @return A merged sparse matrix with combined columns and matching gene order.
#'
#' @export
merge_sparse_matrix <- function(m1, m2) {
  num_cols1 <- ncol(m1)
  num_cols2 <- ncol(m2)

  gene_vector <- unique(c(colnames(m1), colnames(m2)))

  if (num_cols1 < length(gene_vector)) {
    num_missing_cols <- length(gene_vector) - num_cols1
    zero_cols <- Matrix::Matrix(0, nrow = nrow(m1), ncol = num_missing_cols, sparse = TRUE)
    genes_missing_in_m1 <- setdiff(gene_vector, colnames(m1))
    colnames(zero_cols) <- genes_missing_in_m1

    m1 <- cbind(m1, zero_cols)

    m1 <- m1[, gene_vector]
  }

  if (num_cols2 < length(gene_vector)) {
    num_missing_cols <- length(gene_vector) - num_cols2
    zero_cols <- Matrix::Matrix(0, nrow = nrow(m2), ncol = num_missing_cols, sparse = TRUE)
    genes_missing_in_m2 <- setdiff(gene_vector, colnames(m2))
    colnames(zero_cols) <- genes_missing_in_m2
    m2 <- cbind(m2, zero_cols)
    m2 <- m2[, gene_vector]
  }

  combined_gmt <- Matrix::rbind2(m1, m2)

  # if rows have duplicated names, then sum them and keep only one row
  if (any(duplicated(rownames(combined_gmt)))) {
    dup_rows <- unique(rownames(combined_gmt)[duplicated(rownames(combined_gmt))])
    summed_rows <- lapply(dup_rows, function(x) Matrix::colSums(combined_gmt[rownames(combined_gmt) == x, ]))

    # convert summed rows to a matrix
    summed_rows <- do.call(rbind, summed_rows)

    # add names to summed rows
    rownames(summed_rows) <- dup_rows

    # remove duplicated rows
    combined_gmt <- combined_gmt[!rownames(combined_gmt) %in% dup_rows, ]

    # add summed rows

    combined_gmt <- Matrix::rbind2(combined_gmt, summed_rows)
  }
  return(combined_gmt)
}
