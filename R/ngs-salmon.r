##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##




#' @title Import Salmon Abundance Estimates
#'
#' @description This function imports Salmon transcript abundance estimates from Salmon SF files
#' and converts them to a DGEList object for use with edgeR.
#'
#' @param sf.files A named character vector of Salmon SF files.
#' @param count.type The type of counts to extract, either "TPM", "NumReads" or "lengthScaledTPM".
#' @param organism The organism name, either "Hsapiens" or "Mmusculus".
#' @param txOut Whether to output transcript-level estimates in addition to gene-level.
#'
#' @details This function takes a named character vector of Salmon SF files, extracts the abundance
#' estimates, converts them to a DGEList object with integer counts scaled to library size.
#' It requires the EnsDb and org packages for the corresponding organism.
#'
#' @return A DGEList object containing the abundance estimates, ready for use with edgeR.
#'
#' @export
ngs.tximportSalmon <- function(sf.files, count.type = "lengthScaledTPM", organism = "Hsapiens",
                               txOut = FALSE) {
  if (is.null(names(sf.files))) stop("sf.files must be named!")
  if (!all(file.exists(sf.files))) stop("missing SF files")

  if (organism == "Hsapiens") {
    edb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
    org <- org.Hs.eg.db::org.Hs.eg.db
  }
  if (organism == "mouse") {
    edb <- EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79
    org <- org.Mm.eg.db::org.Mm.eg.db
  }

  ## ------------------------------------------------------------
  ## get transcript annotation
  ## ------------------------------------------------------------
  ## then the transcript to gene file from Ensembl
  listColumns(edb)
  daf <- GenomicFeatures::transcripts(edb,
    columns = c(
      "tx_id", "gene_id", "entrezid",
      "gene_name", "gene_biotype", "name"
    ),
    return.type = "DataFrame"
  )
  dim(daf)
  Matrix::head(daf)
  annot_genes <- AnnotationDbi::select(org,
    keytype = "ENSEMBL", keys = daf$gene_id,
    columns = c("SYMBOL", "REFSEQ", "ENTREZID", "GENENAME")
  )
  if (1) {
    ## Add REFSEQ????
    cat("quering biomaRt...\n")

    ensembl <- biomaRt::useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
    bm.res <- biomaRt::getBM(
      attributes = c("refseq_mrna", "ensembl_gene_id", "ensembl_transcript_id", "external_gene_name"),
      filters = "ensembl_transcript_id",
      values = daf$tx_id,
      mart = ensembl
    )
    Matrix::head(bm.res)
    dim(bm.res)
    daf$refseq2 <- NULL
    daf$refseq <- bm.res[match(daf$tx_id, bm.res[, "ensembl_transcript_id"]), "refseq_mrna"]
    Matrix::head(daf)
  }

  dim(annot_genes)
  Matrix::head(annot_genes)
  daf$gene_title <- annot_genes$GENENAME[match(daf$gene_id, annot_genes$ENSEMBL)]
  Matrix::head(daf, 200)

  ## ------------------------------------------------------------
  ## Import Salmon files using tximport
  ## ------------------------------------------------------------
  tx2gene <- daf[, c("tx_id", "gene_id")] ## map to Ensemble gene ID

  ## now import all files and collapse to gene. The 'lengthScaleTPM'
  ## is essential for LIMMA/VOOM and Deseq2 handles this fine (see
  ## code for DESeqDataSetFromTximport)


  if (txOut == FALSE) {
    ## collapse gene level
    txi <- tximport::tximport(sf.files,
      type = "salmon",
      countsFromAbundance = count.type, txOut = FALSE,
      tx2gene = tx2gene, ignoreTxVersion = TRUE
    )
    daf0 <- daf
    daf <- daf[match(rownames(txi$counts), daf$gene_id), ]

    ## collapse also transcript-level annotation
    tx.ids <- tapply(tx2gene$tx_id, tx2gene$gene_id, paste, collapse = ",")
    daf$tx_id <- tx.ids[match(daf$gene_id, names(tx.ids))] ## replace with TX list
    if (!is.null(daf$refseq)) {
      rfq.ids <- tapply(daf0$refseq, daf0$gene_id, function(x) paste(setdiff(x, ""), collapse = ","))
      daf$refseq <- rfq.ids[match(daf$gene_id, names(rfq.ids))] ## replace with TX list
    }
    remove(daf0)
  } else {
    ## transcript level
    txi <- tximport::tximport(sf.files,
      type = "salmon",
      countsFromAbundance = count.type, txOut = TRUE,
      tx2gene = NULL, ignoreTxVersion = TRUE
    )
    tx.id <- sub("[.][0-9]*$", "", rownames(txi$counts))
    daf <- daf[match(tx.id, daf$tx_id), ]
  }

  txi$genes <- daf[, c("tx_id", "gene_id", "refseq", "gene_name", "gene_biotype", "gene_title")]
  rownames(txi$genes) <- rownames(txi$counts)
  return(txi)
}

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