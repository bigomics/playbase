% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pgx-wgcna.R
\name{pgx.wgcna}
\alias{pgx.wgcna}
\title{WGCNA network construction and module detection}
\usage{
pgx.wgcna(
  pgx,
  minmodsize = 30,
  power = 6,
  cutheight = 0.25,
  deepsplit = 2,
  ngenes = 1000
)
}
\arguments{
\item{pgx}{PGX object containing gene expression data}

\item{minmodsize}{Minimum module size cutoff}

\item{power}{Soft thresholding power for network construction}

\item{cutheight}{Cut height for module dendrogram}

\item{deepsplit}{Number of splits for module dendrogram}

\item{ngenes}{Number of genes to use (most variable)}
}
\value{
List containing WGCNA network and module results
}
\description{
Constructs a weighted gene coexpression network and detects
modules using WGCNA on a PGX object.
}
\details{
This function takes a PGX object containing gene expression data.
It constructs a weighted gene coexpression network using the WGCNA package.
Soft thresholding power is set by \code{power}.

Modules are detected by cutting the module dendrogram at \code{cutheight} and
with \code{deepsplit} splits. Only the \code{ngenes} most variable genes are used.

The output is a list containing the WGCNA network object and module results,
including module assignments, colors, and summary statistics.
}
