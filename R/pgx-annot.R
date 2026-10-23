##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

#' Get probetype annotation for organism and datatype. For multi-omics
#' probe names must be prefixed with data type.
#'
#' Forwards to `playbase.annotation::getProbeAnnotation()`.
#' @export
getProbeAnnotation <- function(organism,
                               probes,
                               datatype,
                               meth_type = NULL,
                               probetype = "",
                               ortholog_species = "Human",
                               annot_table = NULL) {
  playbase.annotation::getProbeAnnotation(
    organism = organism,
    probes = probes,
    datatype = datatype,
    meth_type = meth_type,
    probetype = probetype,
    ortholog_species = ortholog_species,
    annot_table = annot_table
  )
}


#' Cleanup annotation
#'
cleanupAnnotation <- function(genes) {
  if (is.null(genes)) {
    return(NULL)
  }

  ## add missing columns if needed, then reorder
  columns <- c(
    "feature", "symbol", "ortholog", "gene_title", ## "gene_biotype",
    ## "map", "pos", "tx_len",
    "chr", "source", "gene_name"
  )
  missing.cols <- setdiff(columns, colnames(genes))
  missing.cols
  for (a in missing.cols) genes[[a]] <- NA
  #  genes <- genes[, columns]
  #  colnames(genes) <- columns

  # gene_name should ALWAYS be assigned to feature for compatibility
  # with gene_name legacy implementation
  genes$gene_name <- genes$feature

  # add space after ; to conform with playbase <= 1.3.2
  genes$gene_title <- gsub(";[ ]*", "; ", genes$gene_title)

  # trim whitespace
  char.cols <- which(sapply(genes, class) == "character")
  for (k in char.cols) {
    genes[[k]] <- trimws(genes[[k]])
  }

  # rename protein-coding to protein_coding to confirm with playbase <= v1.3.2
  ## genes$gene_biotype <- sub("protein-coding", "protein_coding", genes$gene_biotype)

  # replace NA in gene_ortholog by "" to conform with old
  # pgx objects. For collapsing to symbol this is important.
  genes$ortholog[is.na(genes$ortholog)] <- ""
  if (!is.null(genes$human_ortholog)) genes$human_ortholog[is.na(genes$human_ortholog)] <- ""

  # replace NA or empty symbol by "{feature}" so there is always a readable name
  ii <- which(genes$symbol %in% c(NA, "", "-"))
  genes$symbol[ii] <- paste0("{", genes$feature[ii], "}")
  genes$gene_title[ii] <- "Uknown feature"

  # if organism is human, ortholog should be NA (matching old
  # playbase annot). NEED RETHINK (this is not very consistent).
  if (is.null(genes$ortholog)) genes$ortholog <- NA

  ## reorder
  ordered.cols <- c(columns, setdiff(colnames(genes), columns))
  genes <- genes[, ordered.cols]

  ## Attempt: remove "pos", "tx_len"
  ##  keep <- colnames(genes)[!colnames(genes) %in% c("pos", "tx_len")]
  ##  genes <- genes[, keep]

  genes
}


## ================================================================================
## ================== GET ORTHOLOG FUNCTIONS ======================================
## ================================================================================

#' @title Get human ortholog from given symbols of organism by
#'   aggregation of multiple methods. This package needs internet
#'   connection.
#'
#' @description Forwards to `playbase.annotation::getHumanOrtholog()`.
#'
#' @export
getHumanOrtholog <- function(organism, symbols,
                             ortho.methods = c("homologene","gprofiler","babelgene",
                               "gprofiler2","uppercase"),
                             verbose = 1) {
  playbase.annotation::getHumanOrtholog(
    organism = organism,
    symbols = symbols,
    ortho.methods = ortho.methods,
    verbose = verbose
  )
}
