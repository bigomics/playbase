% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pgx-singlecell.R
\name{pgx.poolCells}
\alias{pgx.poolCells}
\title{Pool single cell counts into pseudo-bulk groups}
\usage{
pgx.poolCells(
  counts,
  ncells,
  groups = NULL,
  stats = "sum",
  clust.method = "umap",
  prior = 1,
  X = NULL,
  meta = NULL,
  verbose = TRUE
)
}
\arguments{
\item{counts}{Single cell count matrix with cells as columns.}

\item{ncells}{Number of cells to sample per group. If NULL, takes all cells in each group.}

\item{groups}{Group assignments for each cell. If NULL, a single group is used.}

\item{stats}{Aggregation method for pooling counts within each group (sum, mean, median, etc)}

\item{clust.method}{Clustering method used to assign groups if groups=NULL (umap, pca, etc)}

\item{prior}{Per-group prior used for pooling. Default is 1 (no prior).}

\item{X}{Optional log-expression matrix for clustering if groups=NULL.}

\item{meta}{Optional cell metadata for clustering if groups=NULL.}

\item{verbose}{Print progress messages?}
}
\value{
Matrix of pooled counts for each group
}
\description{
This function pools single cell RNA-seq count matrices into
pseudo-bulk groups. Cells are assigned to groups based on metadata or
clustering. Counts within each group are aggregated using summation or other statistics.
}
\details{
This function takes a single cell count matrix and pools the counts into pseudo-bulk groups.
If groups are provided, cells are assigned to these groups. Otherwise clustering is performed using umap/pca
on the log-expression matrix X to infer groups. Within each group, cell counts are aggregated into a pseudo-bulk profile
using summation or other statistics. A per-group prior can be used to normalize pooling.
The output is a pooled count matrix with groups as columns.
}
