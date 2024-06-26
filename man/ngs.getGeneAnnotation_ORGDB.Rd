% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pgx-ensembl.R
\name{ngs.getGeneAnnotation_ORGDB}
\alias{ngs.getGeneAnnotation_ORGDB}
\title{Get gene annotation data}
\usage{
ngs.getGeneAnnotation_ORGDB(organism, probes, probe_type)
}
\arguments{
\item{organism}{Organism name (e.g. "Human", "Mouse", "Rat").}

\item{probes}{Character vector of gene/transcript identifiers.}

\item{probe_type}{Type of identifiers provided in probes (e.g. "ENSEMBL").}
}
\value{
Data frame containing annotation data for the provided identifiers.
}
\description{
Retrieves gene annotation data from an organism-specific
annotation package.
}
\details{
Queries an annotation package for the specified organism
using the provided identifiers and probe type. Returns a cleaned and
reformatted data frame with feature identifiers as row names.
}
\examples{
\dontrun{
probes <- c("ENSMUSG00000051951", "ENSMUSG00000033845")
annot <- ngs.getGeneAnnotation_ORGDB(probes, "ENSEMBL", organism = "Mouse")
}
}
