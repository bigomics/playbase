##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

#' @title Get gene annotation data
#'
#' @description This function retrieves gene annotation data including gene
#' symbols, titles, chromosome locations, lengths, and biotypes.
#'
#' @param is.human Logical indicating if organism is human or mouse. Default
#' is TRUE.
#' @param SYMBOL Vector of gene symbols.
#' @param GENE.TITLE Named vector mapping gene symbols to gene titles.
#' @param CHRLOC Named vector mapping gene symbols to chromosome locations.
#' @param TXLEN Named vector mapping gene symbols to transcript lengths.
#' @param MAP Data frame with chromosome map (for human only).
#' @param GENE.BIOTYPE Named vector mapping gene symbols to gene biotypes
#' (for human only).
#'
#' @details This function first checks if the organism is human or mouse based on
#' the is.human parameter. For human, it retrieves gene symbols, titles, chromosome
#' locations, transcript lengths and biotypes using various Bioconductor annotation packages.
#' For mouse, it retrieves symbols, titles, chromosome locations, and transcript lengths.
#' The data is stored in a series of named vectors mapped to the gene symbols.
#'
#' @return A data.frame containing gene annotation data mapped to gene symbols. The columns
#' are: "gene_name", "gene_title", "gene_biotype", "chr", "pos", "tx_len", "map".
#'
#' @examples
#' \dontrun{
#' d <- get_mini_example_data()
#' genes <- sample(rownames(d$counts), 10)
#' annotation <- ngs.getGeneAnnotation(genes)
#' print(annotation)
#' }
#' @export
ngs.getGeneAnnotation <- function(genes) {
  hs.genes <- unique(unlist(as.list(org.Hs.eg.db::org.Hs.egSYMBOL)))
  mm.genes <- unique(unlist(as.list(org.Mm.eg.db::org.Mm.egSYMBOL)))

  ## if multiple genes, take first
  genes0 <- genes ## save
  genes <- sapply(genes, function(s) strsplit(s, split = "[;,\\|]")[[1]][1])

  is.human <- mean(genes %in% hs.genes) > mean(genes %in% mm.genes)
  is.mouse <- mean(genes %in% hs.genes) < mean(genes %in% mm.genes)
  is.human
  is.mouse

  txlen <- SYMBOL <- CHRLOC <- MAP <- NULL
  if (is.human) {
    message("detected organism: human")
    GENE.TITLE <- unlist(as.list(org.Hs.eg.db::org.Hs.egGENENAME))
    SYMBOL <- unlist(as.list(org.Hs.eg.db::org.Hs.egSYMBOL))
    names(GENE.TITLE) <- SYMBOL
    CHRLOC <- as.list(org.Hs.eg.db::org.Hs.egCHRLOC)
    CHRLOC <- CHRLOC[match(names(SYMBOL), names(CHRLOC))]
    MAP <- as.list(org.Hs.eg.db::org.Hs.egMAP)
    MAP <- MAP[match(names(SYMBOL), names(MAP))]
    MAP <- sapply(MAP, paste, collapse = "|")
    names(CHRLOC) <- SYMBOL
    names(MAP) <- SYMBOL

    #
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
    tx <- GenomicFeatures::transcriptLengths(txdb)
    TXLEN <- tapply(tx$tx_len, tx$gene_id, mean, na.rm = TRUE)
    TXLEN <- TXLEN[match(names(SYMBOL), names(TXLEN))]
    names(TXLEN) <- SYMBOL

    ## get gene biotype
    daf <- GenomicFeatures::transcripts(
      EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86,
      columns = c("gene_name", "gene_biotype"),
      return.type = "DataFrame"
    )
    GENE.BIOTYPE <- daf$gene_biotype
    names(GENE.BIOTYPE) <- daf$gene_name
  }
  if (is.mouse) {
    message("detected organism: mouse")
    GENE.TITLE <- unlist(as.list(org.Mm.eg.db::org.Mm.egGENENAME))
    SYMBOL <- unlist(as.list(org.Mm.eg.db::org.Mm.egSYMBOL))
    names(GENE.TITLE) <- SYMBOL
    CHRLOC <- as.list(org.Mm.eg.db::org.Mm.egCHRLOC)
    CHRLOC <- CHRLOC[match(names(SYMBOL), names(CHRLOC))]
    names(CHRLOC) <- SYMBOL
    MAP <- NULL ## no map for mouse???

    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
    tx <- GenomicFeatures::transcriptLengths(txdb)
    TXLEN <- tapply(tx$tx_len, tx$gene_id, mean, na.rm = TRUE)
    TXLEN <- TXLEN[match(names(SYMBOL), names(TXLEN))]
    names(TXLEN) <- SYMBOL

    daf <- GenomicFeatures::transcripts(
      EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79,
      columns = c("gene_name", "gene_biotype"),
      return.type = "DataFrame"
    )
    GENE.BIOTYPE <- daf$gene_biotype
    names(GENE.BIOTYPE) <- daf$gene_name
  }
  remove(txdb)
  remove(tx)

  gene_title <- gene_biotype <- chrom <- loc <- map <- txlen <- NULL
  gene_title <- GENE.TITLE[genes]
  chrloc0 <- CHRLOC[genes]
  loc <- sapply(chrloc0, "[", 1) ## take only first entry !!!
  loc[sapply(loc, is.null)] <- NA
  loc <- abs(as.integer(unlist(loc)))
  chrom <- sapply(chrloc0, function(s) names(s)[1])
  chrom[sapply(chrom, is.null)] <- NA
  chrom <- as.vector(unlist(chrom))
  txlen <- round(TXLEN[genes])
  gene_biotype <- GENE.BIOTYPE[genes]

  map <- paste0("chr", chrom)
  map <- sub("chrNA", NA, map)
  if (!is.null(MAP)) map <- MAP[genes]

  ## get protein info
  ## fill me
  annot <- data.frame(
    gene_name = genes,
    gene_title = gene_title,
    gene_biotype = gene_biotype,
    chr = chrom,
    pos = as.integer(loc),
    tx_len = as.integer(txlen),
    map = map
  )

  rownames(annot) <- genes
  annot
}
