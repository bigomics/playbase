##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


#' @title Get gene annotation data
#'
#' @description Retrieves additional gene annotation data for a Salmon quantification object.
#'
#' @param keys Character vector of row names in the Salmon quant matrix.
#' @param keytype The type of identifiers in `keys` (e.g. 'ENSEMBL').
#' @param gencode Data frame with GENCODE gene annotation data.
#'
#' @details This function takes a Salmon quantification object and retrieves additional
#' gene annotation data for the quantified genes. It extracts the row names from the
#' Salmon counts matrix as the `keys`. It then uses these keys to retrieve corresponding
#' annotation from the provided `gencode` data frame, matching on the gene identifiers.
#'
#' Additional annotation columns added are:
#' - gene_id: GENCODE gene identifier
#' - gene_type: GENCODE gene type
#' - gene_name: Gene name
#' - chr: Chromosome
#' - start: Gene start position
#' - stop: Gene end position
#'
#' The additional annotation is added to the `genes` slot of the Salmon object.
#'
#' @return The Salmon quantification object with additional gene annotation data added.
#'
#' @export
ngs.getGeneAnnot <- function(keys, keytype, gencode) {
  ## add more gene annotation
  ## build gene annotation (expects pure ENSEMBL.ID in rows)
  idx <- match(keys, sub("[.].*", "", gencode$gene_id)) ## assumes single match
  kk <- c("gene_id", "gene_type", "gene_name", "chr", "start", "stop")
  gencode <- gencode[idx, kk]
  dim(gencode)


  ## add annotation using org.Hs.eg.db (NEED RETHINK ON MULTIPLE MATCHES)
  biomaRt::columns(org.Hs.eg.db::org.Hs.eg.db)

  sel.keys <- c("ENTREZID", "GENENAME")
  org.annot <- plotly::select(org.Hs.eg.db::org.Hs.eg.db,
    keys = keys, keytype = keytype, columns = sel.keys
  )
  dim(org.annot)
  idx <- match(keys, org.annot$ENSEMBL) ## assumes single match
  org.annot <- org.annot[idx, ]


  genes <- data.frame(gencode, org.annot)
  rownames(genes) <- keys
  Matrix::head(genes)
  return(genes)
}
