##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##



#' Title
#'
#' @param sampleTable value
#' @param gencode value
#' @param txi value
#'
#' @return
#' @export
#'
#' @examples
ngs.rawSalmon <- function(sampleTable, gencode, txi = NULL) {
  if (is.null(txi)) {
    if (!("sf.file" %in% colnames(sampleTable))) stop("need sf.files in table")
    sf.files <- sampleTable$sf.files
    names(sf.files) <- rownames(sampleTable)
    txi <- ngs.tximportSalmon(sf.files)
  } else {
    if (!all(rownames(sampleTable) == colnames(txi$counts))) {
      stop("provided txi does not match sampleTable")
    }
  }

  keys <- rownames(txi$counts)
  geneTable <- ngs.getGeneAnnot(keys, keytype = "ENSEMBL", gencode)

  ## we like this..
  probegene <- paste0(rownames(geneTable), ":", geneTable$gene_name)
  Matrix::head(probegene)
  rownames(txi$counts) <- probegene
  rownames(txi$abundance) <- probegene
  rownames(geneTable) <- probegene
  all(rownames(txi$counts) == rownames(geneTable))

  raw <- edgeR::DGEList(round(txi$counts), group = NULL) ## we like integer counts...
  raw$group <- NULL
  raw$samples <- sampleTable
  raw$genes <- geneTable

  return(raw)
}


#' Title
#'
#' @param sf.files value
#' @param count.type value
#' @param organism value
#' @param txOut value
#'
#' @return
#' @export
#'
#' @examples
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

#' Title
#'
#' @param keys value
#' @param keytype value
#' @param gencode value
#'
#' @return
#' @export
#'
#' @examples
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
