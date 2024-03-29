% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pgx-drugs.R
\name{pgx.computeDrugEnrichment}
\alias{pgx.computeDrugEnrichment}
\title{Compute drug enrichment from expression data}
\usage{
pgx.computeDrugEnrichment(
  obj,
  X,
  xdrugs,
  drug_info = NULL,
  methods = c("GSEA", "cor"),
  nmin = 15,
  nprune = 250,
  contrast = NULL
)
}
\arguments{
\item{X}{Numeric gene expression matrix}

\item{methods}{Methods for enrichment (rank, gsea)}

\item{nmin}{Minimum number of genes per drug set}

\item{drugs}{Character vector of drug names}
}
\value{
List with enrichment results for each drug
}
\description{
Computes enrichment of drug gene sets in expression data using
rank correlation or GSEA.
}
\details{
This function takes a gene expression matrix and a set of drugs. It extracts
gene sets for each drug from MSigDB. Enrichment is computed by:

1) Rank correlation of drug set ranks with experiment ranks.
2) GSEA using the drug sets as gene sets.

Drugs with fewer than nmin genes are filtered out. The output is a list
containing the enrichment results for each drug.
}
