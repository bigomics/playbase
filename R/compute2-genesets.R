##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


#' Compute Test Genesets
#'
#' Computes test genesets for the given multi-omics data.
#'
#' @param pgx A data object representing multi-omics data.
#' @param max.features Maximum number of features to consider.
#' @param custom.geneset Custom geneset object.
#' @param test.methods Methods to use for testing.
#' @param remove.outputs Logical indicating whether to remove large outputs.
#'
#' @return The updated \code{pgx} object with computed test genesets.
#'
#' @export
compute_testGenesets <- function(pgx,
                                 max.features = 1000,
                                 custom.geneset = list(gmt = NULL, info = NULL),
                                 test.methods = c("gsva", "camera", "fgsea"),
                                 remove.outputs = TRUE) {
  if (!"X" %in% names(pgx)) {
    stop("[compute_testGenesets] FATAL : object must have normalized matrix X")
  }

  if(0) {
    max.features = 1000;custom.geneset = list(gmt = NULL, info = NULL);
    test.methods = c("gsva", "camera", "fgsea");remove.outputs = TRUE
  }
  
  if (is.null(pgx$genes$human_ortholog)) {
    # this is needed in case the species is human, and we dont have the homolog column or if we have an old pgx
    # which will ensure consistency between old and new pgx
    pgx$genes$human_ortholog <- NA
  }

  ## -----------------------------------------------------------
  ## Load huge geneset matrix
  ## -----------------------------------------------------------
  G <- Matrix::t(playdata::GSETxGENE)

  ## -----------------------------------------------------------
  ## Filter genes
  ## -----------------------------------------------------------

  ## filter genes by gene or homologous, if it exists
  # replace "" to NA in pgx$genes$human_ortholog
  if (pgx$organism != "Human") {
    human_genes <- ifelse(!is.na(pgx$genes$human_ortholog),
      pgx$genes$human_ortholog,
      pgx$genes$symbol
    )
  } else {
    human_genes <- pgx$genes$symbol
  }
  genes <- pgx$genes$symbol

  # Change HUMAN gene names to species symbols if NOT human and human_ortholog column is NOT all NA
  G <- G[rownames(G) %in% human_genes,, drop = FALSE]

  if (pgx$organism != "Human" && !all(is.na(pgx$genes$human_ortholog))) {
    rownames(G) <- pgx$genes$symbol[match(rownames(G), pgx$genes$human_ortholog)]
  }

  # Normalize G after removal of genes
  G <- playbase::normalize_cols(G)

  ## -----------------------------------------------------------
  ## Filter gene sets
  ## -----------------------------------------------------------

  ## filter gene sets on size
  cat("Filtering gene sets on size...\n")
  gmt.size <- Matrix::colSums(G != 0)
  size.ok <- which(gmt.size >= 15 & gmt.size <= 400)

  # If dataset is too small that size.ok == 0, then select top 100
  ## if (length(size.ok) == 0) {
  ##   size.ok <- head(sample(1:ncol(G)),100)
  ## }
  G <- G[, size.ok, drop = FALSE]

  ## -----------------------------------------------------------
  ## Add random genesets
  ## -----------------------------------------------------------

  add.gmt <- NULL

  rr <- sample(15:400,100)
  gg <- rownames(pgx$X)
  random.gmt <- lapply( rr, function(n) head(sample(gg),min(n,length(gg)/2)))
  names(random.gmt) <- paste0("TEST:random_geneset.",1:length(random.gmt))
  add.gmt <- random.gmt

  ## -----------------------------------------------------------
  ## Add custom genesets
  ## -----------------------------------------------------------
  
  if (!is.null(custom.geneset$gmt)) {
    add.gmt <- c( add.gmt, custom.geneset$gmt )
  }
  
  if (!is.null(add.gmt)) {
    # convert gmt standard to SPARSE matrix
    custom_gmt <- playbase::createSparseGenesetMatrix(
##    gmt.all = custom.geneset$gmt,
      gmt.all = add.gmt,
      min.geneset.size = 3,
      max.geneset.size = 9999,
      min_gene_frequency = 1,
      all_genes = pgx$all_genes,
      annot = pgx$genes,
      filter_genes = FALSE
    )
    custom_gmt <- Matrix::t(custom_gmt)    
    custom_gmt <- custom_gmt[rownames(custom_gmt) %in% genes,, drop = FALSE]
    custom_gmt <- playbase::normalize_cols(custom_gmt)

    G <- Matrix::t(playbase::merge_sparse_matrix(m1 = Matrix::t(G), m2 = Matrix::t(custom_gmt)))
    remove(custom_gmt)
  }

  ## -----------------------------------------------------------
  ## create the full GENE matrix (always collapsed by gene)
  ## -----------------------------------------------------------

  X <- pgx$X
  if (!all(rownames(X) %in% pgx$genes$symbol)) {
    X <- rename_by(X, pgx$genes, "symbol")
    X <- X[!rownames(X) == "", , drop = FALSE]
    if (any(duplicated(rownames(X)))) {
      X <- log2(rowsum(2**X, rownames(X)))
    }
  }

  ## if reduced samples
  ss <- rownames(pgx$model.parameters$exp.matrix)
  if (!is.null(ss)) {
    X <- X[, ss, drop = FALSE]
  }

  ## -----------------------------------------------------------
  ## create the GENESETxGENE matrix
  ## -----------------------------------------------------------
  cat("Matching gene set matrix...\n")
  gg <- rownames(X)
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
    cat("Reducing gene set matrix...\n")
    ## Reduce gene sets by selecting top varying genesets. We use the
    ## very fast sparse rank-correlation for approximate single sample
    ## geneset activation.
    cX <- X - rowMeans(X, na.rm = TRUE) ## center!
    cX <- apply(cX, 2, rank)
    gsetX <- qlcMatrix::corSparse(G, cX) ## slow!
    grp <- pgx$model.parameters$group
    gsetX.bygroup <- NULL
    ## If groups/conditions are present we calculate the SD by group
    if (!is.null(grp)) {
      gsetX.bygroup <- tapply(1:ncol(gsetX), grp, function(i) rowMeans(gsetX[, i, drop = FALSE]))
      gsetX.bygroup <- do.call(cbind, gsetX.bygroup)
      sdx <- apply(gsetX.bygroup, 1, stats::sd)
    } else {
      sdx <- matrixStats::rowSds(gsetX)
    }
    rm(gsetX)
    rm(gsetX.bygroup)
    names(sdx) <- colnames(G)
    jj <- Matrix::head(order(-sdx), max.features)
    must.include <- "hallmark|kegg|^go|^celltype|^pathway|^custom"
    jj <- unique(c(jj, grep(must.include, colnames(G), ignore.case = TRUE)))
    jj <- jj[order(colnames(G)[jj])] ## sort alphabetically
    G <- G[, jj, drop = FALSE]
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
  cat(">>> Testing gene sets with methods:", test.methods, "\n")

  ## convert to gene list
  gmt <- playbase::mat2gmt(G)
  Y <- pgx$samples
  gc()

  gset.meta <- playbase::gset.fitContrastsWithAllMethods(
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

  remove(X)
  remove(Y)
  remove(G)
  remove(gmt)
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
