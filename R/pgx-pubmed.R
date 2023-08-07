##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

NCORE <- function() {
  parallel::detectCores(all.tests = FALSE, logical = TRUE) / 2
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
