% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gx-limma.r
\name{gx.limmaF}
\alias{gx.limmaF}
\title{Differential expression analysis with limma}
\usage{
gx.limmaF(
  X,
  pheno,
  B = NULL,
  fdr = 0.05,
  compute.means = TRUE,
  lfc = 0.2,
  max.na = 0.2,
  sort.by = "FC",
  ref = c("ctrl", "ctr", "control", "dmso", "nt", "0", "0h", "0hr", "non", "no", "not",
    "neg", "negative", "ref", "veh", "vehicle", "wt", "wildtype", "untreated", "normal",
    "false", "healthy"),
  trend = FALSE,
  method = 1,
  verbose = 1
)
}
\arguments{
\item{X}{Numeric gene expression matrix with genes in rows and samples in columns.}

\item{pheno}{Data frame with phenotype data for samples. Must have column named 'group'.}

\item{B}{Data frame with batch data for samples. Default is NULL.}

\item{fdr}{FDR threshold for significance. Default is 0.05.}

\item{compute.means}{Logical indicating whether to compute group means. Default is TRUE.}

\item{lfc}{Log fold change threshold. Default is 0.2.}

\item{max.na}{Maximum missing value fraction for gene filtering. Default is 0.2.}

\item{ref}{Character vector of reference group names to use as baseline. Default is common control names.}

\item{trend}{Logical indicating whether to fit a trend model. Default is FALSE.}

\item{verbose}{Verbosity level. Default is 1.}
}
\value{
List with differential expression results, including:
\itemize{
  \item tab - Data frame with stats for all genes
  \item top - Data frame with stats for top significant genes
  \item fstats - Data frame with F statistics for all genes
  \item means - Data frame with mean expression by group
}
}
\description{
Differential expression analysis with limma
}
\details{
This function performs differential expression analysis on \code{X} using limma.
It handles filtering, model design matrices, and output formatting.
}
\examples{
\dontrun{
# TODO
}
}
