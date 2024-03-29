% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pgx-pubmed.R
\name{pmid.annotateEdges}
\alias{pmid.annotateEdges}
\title{Annotate edges with shared genes in a Pubmed network graph}
\usage{
pmid.annotateEdges(gr)
}
\arguments{
\item{gr}{An igraph network graph object generated from Pubmed data}
}
\value{
The input graph object gr with edge weights and gene annotations added.
}
\description{
Annotates edges in a Pubmed citation network graph with the genes shared between connected articles.
}
\details{
This function takes a network graph generated from Pubmed citation data, with articles as nodes.
It extracts the genes associated with each article node using the V(gr)$genes attribute.
For each edge, it finds the intersecting genes between the connected nodes.
The number of shared genes is assigned as the edge weight.
The shared gene symbols are assigned to a new edge attribute E(gr)$genes.
}
