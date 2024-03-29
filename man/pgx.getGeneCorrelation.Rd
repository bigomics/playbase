% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pgx-correlation.R
\name{pgx.getGeneCorrelation}
\alias{pgx.getGeneCorrelation}
\title{Get correlation of a gene with other genes}
\usage{
pgx.getGeneCorrelation(gene, xref)
}
\arguments{
\item{gene}{Character string specifying the target gene name.}

\item{xref}{List of gene expression matrices to compute correlation against.}
}
\value{
Matrix of correlation values between the target gene and other genes.
}
\description{
Computes correlation between a given gene and all other genes using
external gene expression data.
}
\details{
This function takes a gene name and a list of gene expression matrices.
It calculates the correlation between the target gene and all other genes in
each expression data matrix.

The function returns a matrix with genes in rows and datasets in columns.
Each column contains the correlation values between the target gene and
other genes in that dataset.
}
