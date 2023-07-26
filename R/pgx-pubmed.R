##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

NCORE <- function() {
  parallel::detectCores(all.tests = FALSE, logical = TRUE) / 2
}


#' @title Get gene context from PubMed
#'
#' @description 
#' Retrieves PubMed article information related to a gene and keyword context.
#'
#' @param gene Gene name or symbol. 
#' @param keyword Optional keyword or phrases to search for context.
#'
#' @details
#' This function searches PubMed for articles related to the input \code{gene}, 
#' optionally filtered by \code{keyword} context phrases.
#' 
#' It first searches PubMed for the gene name and synonyms. It then filters these 
#' articles based on matches to the \code{keyword} terms.  The PubMed ids (PMIDs) 
#' of matching articles are returned.
#'
#' Statistics are calculated to assess enrichment of the \code{keyword} within articles 
#' matching the \code{gene}. These include a Fisher's exact test p-value, and word 
#' co-occurrence statistics.
#'
#' @return 
#' A list containing:
#' \itemize{
#' \item rifs: Character vector of PMIDs for articles matching the gene and keyword filter.
#' \item table: Contingency table used for calculating p-value.  
#' \item p.value: Fisher's exact test p-value.
#' \item context: Word context scores.
#' }
#'
#' @export
pmid.getGeneContext <- function(gene, keyword) {
  gene1 <- c(gene, sub("([0-9])", "-\\1", gene))
  gene1 <- paste0("^", gene1, "$|^", gene1, "[-]")
  gene1 <- paste(gene1, collapse = "|")

  if (gene %in% biomaRt::keys(org.Hs.eg.db::org.Hs.egALIAS2EG)) {
    gname <- get(get(gene, org.Hs.eg.db::org.Hs.egALIAS2EG), org.Hs.eg.db::org.Hs.egGENENAME)
    gname <- gsub("[, -]", ".", gname)
    gene1 <- paste0(gene1, "|", gname)
  } else if (gene %in% biomaRt::keys(org.Mm.eg.db::org.Mm.egALIAS2EG)) {
    gname <- get(get(gene, org.Mm.eg.db::org.Mm.egALIAS2EG), org.Mm.eg.db::org.Mm.egGENENAME)
    gname <- gsub("[, -]", ".", gname)
    gene1 <- paste0(gene1, "|", gname)
  }

  rif.words <- colnames(GENERIF.MATRIX)
  ii <- grepl(gene1, rif.words, ignore.case = TRUE)
  match0 <- rowSums(GENERIF.MATRIX[, ii, drop = FALSE]) > 0
  match1 <- rep(1, length(match0))
  if (length(keyword) > 0) {
    i <- 1
    for (i in 1:length(keyword)) {
      jj <- grepl(keyword[i], rif.words, ignore.case = TRUE)
      match1 <- match1 & (rowSums(GENERIF.MATRIX[, jj, drop = FALSE]) > 0)
    }
  }
  sel <- ((match0 * match1) > 0)
  rif.hits <- rownames(GENERIF.MATRIX)[sel]
  rif.hits <- rif.hits[!duplicated(rif.hits)]

  ## calculate P-value for this keyword
  A <- table(gene = match0, keyword = match1)
  pv <- NA
  if (nrow(A) == 2 && ncol(A) == 2) {
    pv <- fisher.test(A, alternative = "greater")$p.value
  }

  context1 <- NULL
  match2 <- ((match0 * match1) > 0)
  m0 <- Matrix::colSums(match2 * (GENERIF.MATRIX != 0))
  m1 <- Matrix::colSums(GENERIF.MATRIX != 0)
  pp <- corpora::fisher.pval(m0, sum(match2) + 1, m1, nrow(GENERIF.MATRIX) + 1, alternative = "greater")
  pp <- sort(pp)
  qq <- p.adjust(pp)
  qq <- sort(qq)
  context1 <- Matrix::head(qq[qq < 1], 100)

  out <- list(rifs = rif.hits, table = A, p.value = pv, context = context1)
  return(out)
}


#' @title Get PubMed context for a gene
#'
#' @description 
#' Retrieves PubMed abstract text containing a gene name and context words.
#'
#' @param gene Gene name or symbol to search for.
#' @param context Context words or phrases to search for along with gene.
#'
#' @details
#' This function searches PubMed for articles containing the specified gene name 
#' and context words. It uses the EUtils API to search titles and abstracts in 
#' PubMed for the gene name and context terms.
#'
#' It returns a list containing the PubMed IDs, article titles, and abstract 
#' excerpts (RIF strings) containing both the gene and context.
#'
#' @return
#' A list with the PubMed IDs, article titles, and abstract excerpts containing
#' the gene name and context words.
#'
#' @examples
#' \dontrun{
#' res <- pmid.getPubMedContext("TP53", "DNA damage")
#' }
#' @export
pmid.getPubMedContext <- function(gene, context) {
  res <- EUtilsSummary(
    paste0(gene, "[sym] AND ", context),
    type = "esearch", db = "pubmed", datetype = "pdat",
    mindate = 2000, maxdate = 2099, retmax = 1000
  )
  QueryCount(res)
  if (QueryCount(res) == 0) {
    return(NULL)
  }
  gene1 <- c(gene, sub("([0-9])", "-\\1", gene), sub("([0-9])", "[ ]\\1", gene))
  gene1 <- paste(gene1, collapse = "|")
  gname <- get(get(toupper(gene), org.Hs.eg.db::org.Hs.egALIAS2EG), org.Hs.eg.db::org.Hs.egGENENAME)
  gene2 <- paste0(gene1, "|", gsub("[ -]", ".", gname))
  extractRIF <- function(a) {
    s <- strsplit(a, split = "[.;:]")[[1]]
    hit <- grepl(gene2, s, ignore.case = TRUE) & grepl(context, s, ignore.case = TRUE)
    if (!any(hit)) {
      return(NULL)
    }
    s <- paste(s[hit], collapse = ".")
    sub("^[ ]*", "", s)
  }
  fetch <- EUtilsGet(res)
  tt <- ArticleTitle(fetch)
  aa <- AbstractText(fetch)
  tt2 <- paste(tt, ".", aa)
  pp <- PMID(fetch)
  rif <- lapply(tt2, extractRIF)
  rif[sapply(rif, is.null)] <- NA
  rif <- unlist(rif)
  rif[!is.na(rif)] <- paste0(rif, " (PMID:", pp, ")")[!is.na(rif)]
  out <- list(pmid = pp, title = tt, rif = rif)
  return(out)
}


#' @title Build PMID annotation matrix
#'
#' @description 
#' Builds an annotation matrix mapping PubMed IDs to gene symbols 
#' based on mappings in org.Hs.eg.db.
#'
#' @param None
#' 
#' @details
#' This function retrieves PubMed ID to Entrez Gene ID mappings from 
#' org.Hs.eg.db and converts them to a sparse matrix mapping PMIDs to 
#' gene symbols. It filters to PMIDs associated with <=10 genes.
#'
#' It collapses duplicate PMID mappings and builds a sparse matrix 
#' with rows as PMIDs, columns as gene symbols, and entries indicating 
#' an association between a PMID and symbol.
#'
#' @return 
#' A sparse matrix with rows as PMIDs, columns as gene symbols, and 
#' entries indicating a mapping between a PMID and symbol.
#'
#' @export
pmid.buildMatrix <- function() {
  pmid <- as.list(org.Hs.eg.db::org.Hs.egPMID2EG)
  symbol <- as.list(org.Hs.eg.db::org.Hs.egSYMBOL)
  eg <- names(symbol)
  symbol <- sapply(symbol, "[", 1)
  names(symbol) <- as.character(eg)
  ngene <- sapply(pmid, length)
  pmid <- pmid[which(ngene <= 10)]

  ## collapse duplicates
  pmid.gg <- sapply(pmid, paste, collapse = ",")
  idx <- tapply(names(pmid), pmid.gg, function(x) x)
  pmid <- pmid[sapply(idx, "[", 1)]
  idx <- lapply(idx, function(x) paste0("PMID:", x))
  names(pmid) <- sapply(idx, paste, collapse = ",")

  ## build PMID2SYMBOL matrix

  idx0 <- parallel::mclapply(1:length(pmid), function(i) cbind(i, which(eg %in% pmid[[i]])),
    mc.cores = NCORE()
  )
  idx <- do.call(rbind, idx0)
  P <- Matrix::sparseMatrix(
    i = idx[, 1], j = idx[, 2], x = rep(1, nrow(idx)),
    dims = c(length(idx0), length(eg))
  )
  rownames(P) <- names(pmid)
  colnames(P) <- symbol
  P <- P[, which(Matrix::colSums(P) > 0)]
  return(P)
}


#' @title Build gene association graph from PubMed identifiers
#'
#' @description Builds an undirected gene association graph from PubMed 
#' identifiers (PMIDs) mapped to gene symbols.
#' 
#' @param P Sparse matrix mapping PMIDs to gene symbols 
#'
#' @details This function takes as input a sparse matrix \code{P} that maps 
#' PubMed identifiers (PMIDs) to gene symbols. 
#' It filters \code{P} to only include PMIDs mapped to 2-10 genes.
#' 
#' It then calculates a co-occurrence matrix \code{M} between genes by multiplying 
#' \code{P} and its transpose. \code{M[i,j]} indicates the number of PMIDs shared 
#' between gene i and gene j.
#'
#' An undirected graph is constructed from \code{M} using the igraph package, 
#' with genes as nodes and edge weights proportional to their PMID co-occurrences.
#'
#' @return An igraph undirected graph object representing the gene association network.
#'
#' @export
pmid.buildGraph <- function(P) {
  P <- P[which(Matrix::rowSums(P) <= 10), ]
  P <- P[which(Matrix::rowSums(P) >= 2), ]
  P <- P[, which(Matrix::colSums(P) > 0)]
  P[1:10, 1:10]

  ## create graph from overlap
  M <- P[, ] %*% t(P[, ])
  diag(M) <- 0
  object.size(M)
  gr <- igraph::graph_from_adjacency_matrix(M,
    mode = "undirected",
    diag = FALSE, weighted = TRUE
  )
  igraph::V(gr)$name <- rownames(M)
  gr <- igraph::subgraph.edges(gr, which(igraph::E(gr)$weight > 0))

  P <- P[igraph::V(gr)$name, ]
  pmids <- parallel::mclapply(igraph::V(gr)$name, function(x) gsub("PMID:", "", strsplit(x, split = ",")[[1]]))
  nref <- sapply(pmids, length)
  vgenes <- parallel::mclapply(1:nrow(P), function(i) names(which(P[i, ] != 0)), mc.cores = NCORE())
  vgenes2 <- unlist(sapply(vgenes, paste, collapse = ","))
  igraph::V(gr)$size <- nref
  igraph::V(gr)$genes <- vgenes
  igraph::V(gr)$pmid <- pmids
  return(gr)
}

#' @title Annotate edges with shared genes in a Pubmed network graph 
#'
#' @description Annotates edges in a Pubmed citation network graph with the genes shared between connected articles.
#' 
#' @param gr An igraph network graph object generated from Pubmed data
#'
#' @details This function takes a network graph generated from Pubmed citation data, with articles as nodes.
#' It extracts the genes associated with each article node using the V(gr)$genes attribute.
#' For each edge, it finds the intersecting genes between the connected nodes. 
#' The number of shared genes is assigned as the edge weight.
#' The shared gene symbols are assigned to a new edge attribute E(gr)$genes.
#'
#' @return The input graph object gr with edge weights and gene annotations added.
#'
#' @export
pmid.annotateEdges <- function(gr) {
  ee <- igraph::get.edges(gr, igraph::E(gr))
  g1 <- igraph::V(gr)[ee[, 1]]$genes
  g2 <- igraph::V(gr)[ee[, 2]]$genes
  shared.genes <- mapply(intersect, g1, g2)
  nshared <- sapply(shared.genes, length)
  igraph::E(gr)$genes <- shared.genes
  igraph::E(gr)$weight <- nshared
  return(gr)
}


#' @title Extract gene-specific subgraph from PubMed network  
#'
#' @param gr An igraph network object generated from PubMed data 
#' @param gene Character vector of gene symbols to extract network for
#' @param nmin Minimum number of genes in connected components to keep
#'
#' @return An igraph object containing the gene-specific subgraph
#'
#' @description Extracts a subgraph related to the specified genes from a PubMed network graph.
#'
#' @details This function takes a PubMed network graph \code{gr} and a vector of \code{gene} symbols. 
#' It extracts the nodes containing those genes and the edges between them.  
#' Connected components with fewer than \code{nmin} genes are removed.
#'
#' Edge weights and gene annotations are added to the subgraph.
#' 
#' @export
pmid.extractGene <- function(gr, gene, nmin = 3) {
  jj <- c()
  for (g in gene) {
    j1 <- which(sapply(igraph::V(gr)$genes, function(s) (gene %in% s)))
    jj <- c(jj, j1)
  }
  jj <- unique(jj)
  ngene1 <- sapply(igraph::V(gr)$genes[jj], length)
  gr1 <- igraph::induced_subgraph(gr, jj)
  if (verbose > 0) cat("annotating edges...\n")
  gr1 <- pmid.annotateEdges(gr1)
  nshared <- unlist(parallel::mclapply(igraph::E(gr1)$genes, length, mc.cores = NCORE()))
  ee <- which(!(nshared == 1 & sapply(igraph::E(gr1)$genes, "[", 1) %in% gene))
  gr1 <- igraph::subgraph.edges(gr1, ee, delete.vertices = TRUE)
  cmp <- igraph::components(gr1)
  jj <- which(cmp$membership %in% which(cmp$csize >= nmin))
  gr1 <- igraph::induced_subgraph(gr1, jj)
  return(gr1)
}


#' Link to PubMed article
#'
#' @title Create HTML link to PubMed article
#'
#' @param s PubMed ID (PMID) to link to. 
#' 
#' @return HTML link to the PubMed article page.
#'
#' @description Generates an HTML hyperlink for a given PubMed ID (PMID).
#'
#' @details This function takes a PubMed ID as input and returns 
#' an HTML hyperlink that will open the article page on PubMed website.
#' 
#' The PMID is embedded in the link URL along with link text containing 
#' the PMID. This allows easy creation of hyperlinks to PubMed articles.
#'
#' @examples
#' \dontrun{
#' link <- pubmedlink("12345678")
#' cat(link)
#' }
#' @export
pubmedlink <- function(s) {
  paste0(
    "<a href='https://www.ncbi.nlm.nih.gov/pubmed/", s,
    "' target='_blank'>PMID:", s, "</a>"
  )
}
