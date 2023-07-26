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

#' Detect organism from NGS object
#'
#' @title Detect organism from NGS object
#'
#' @description Detects the organism for an NGS object by
#' looking at the case of gene identifiers.
#'
#' @param ngs An NGS object.
#'
#' @details This function detects the organism for an NGS object
#' by calculating the ratio of lowercase characters in the
#' rownames of the counts matrix. It assumes human data will have a
#' high proportion of uppercase ENSEMBL identifiers, while mouse
#' data will have more lowercase.
#'
#' It calculates the ratio of lowercase characters in the first 100
#' characters of the rownames of ngs$counts.
#' If the ratio is greater than 0.5 it assigns 'mouse', otherwise 'human'.
#'
#' @return A character string indicating the detected organism ('human' or 'mouse').
#'
#' @examples
#' \dontrun{
#' # TODO
#' }
#' @export
ngs.detectOrganism <- function(ngs) {
  lowcase.ratio <- mean(grepl("[a-z]", substring(rownames(ngs$counts), 2, 100)))
  c("human", "mouse")[1 + 1 * (lowcase.ratio > 0.5)]
}


#' @title Match gene names to feature identifiers
#'
#' @description Matches gene names in a character vector to the
#' corresponding feature identifiers in an NGS object.
#'
#' @param ngs An NGS object containing feature annotation data.
#' @param genes Character vector of gene names to match.
#'
#' @details This function takes a character vector of gene names and
#' finds the matching feature identifiers in the feature annotation data
#' contained in the \code{ngs$genes} data frame of an NGS object. It
#' matches the \code{gene_name} column in \code{ngs$genes}.
#'
#' Matching is case-insensitive. The function returns the feature
#' identifiers from the rownames of \code{ngs$genes} that match the
#' provided gene names.
#'
#' @return Character vector of matching feature identifiers
#' (typically Ensembl IDs).
#'
#' @examples
#' \dontrun{
#' # TODO
#' }
#'
#' @export
ngs.matchFeatures <- function(ngs, genes) {
  jj <- match(toupper(genes), toupper(ngs$genes$gene_name))
  rownames(ngs$genes)[jj]
}

#' Collapse duplicate genes in NGS dataset
#'
#' @param ngs An NGS dataset object
#'
#' @return The NGS dataset object with collapsed genes
#'
#' @details
#' This function collapses duplicate gene names in an NGS dataset by summing their counts.
#'
#' It first identifies genes that occur only once or more than once in the dataset.
#' For uniquely occurring genes, it keeps them as is.
#' For duplicate genes, it sums their counts by gene name.
#'
#' The gene and count matrices in the NGS dataset object are then updated accordingly.
#'
#' @export
ngs.collapseByGene <- function(ngs) {
  gene <- as.character(ngs$genes$gene_name)
  p1 <- names(which(table(gene) == 1))
  p2 <- names(which(table(gene) > 1))
  if (length(p2) == 0) {
    gene <- as.character(ngs$genes$gene_name)
    rownames(ngs$genes) <- gene
    rownames(ngs$counts) <- gene
    return(ngs)
  }
  j1 <- which(gene %in% p1)
  j2 <- which(gene %in% p2)
  x1 <- ngs$counts[j1, , drop = FALSE]
  rownames(x1) <- ngs$genes$gene_name[j1]
  x2 <- apply(ngs$counts[j2, , drop = FALSE], 2, function(x) tapply(x, gene[j2], sum))
  if (length(p2) == 1) {
    x2 <- matrix(x2, nrow = 1)
    rownames(x2) <- p2
    colnames(x2) <- colnames(x1)
  }
  x1 <- rbind(x1, x2)
  x1 <- x1[!(rownames(x1) %in% c(NA, "", "NA")), , drop = FALSE]
  x1 <- x1[order(rownames(x1)), , drop = FALSE]
  ngs$genes <- ngs$genes[match(rownames(x1), ngs$genes$gene_name), ]
  ngs$counts <- x1
  rownames(ngs$genes) <- rownames(ngs$counts) <- rownames(x1)
  return(ngs)
}
