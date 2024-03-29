% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute2-genesets.R
\name{createSparseGenesetMatrix}
\alias{createSparseGenesetMatrix}
\title{Create Sparse Geneset Matrix}
\usage{
createSparseGenesetMatrix(
  gmt.all,
  min.geneset.size = 15,
  max.geneset.size = 500,
  min_gene_frequency = 10,
  all_genes = NULL,
  annot = NULL,
  filter_genes = TRUE
)
}
\arguments{
\item{gmt.all}{The list of genesets in GMT format.}

\item{min.geneset.size}{The minimum size of a geneset to be included in the matrix. Default is 15.}

\item{max.geneset.size}{The maximum size of a geneset to be included in the matrix. Default is 500.}

\item{min_gene_frequency}{The minimum frequency of a gene to be included in the matrix. Default is 10.}
}
\value{
A sparse matrix representing the filtered genesets.
}
\description{
Creates a sparse matrix representation of genesets from a given GMT file, filtered by size and gene frequency.
}
