% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pgx-cluster.R
\name{pgx.clusterMatrix.DEPRECATED}
\alias{pgx.clusterMatrix.DEPRECATED}
\title{Cluster rows of a matrix}
\usage{
pgx.clusterMatrix.DEPRECATED(
  X,
  perplexity = 30,
  dims = c(2, 3),
  ntop = 1000,
  npca = 50,
  prefix = "c",
  row.center = TRUE,
  row.scale = FALSE,
  find.clusters = TRUE,
  kclust = 1,
  clust.detect = c("louvain", "hclust"),
  method = c("tsne", "umap", "pca")
)
}
\arguments{
\item{X}{Numeric matrix with rows as features to cluster}

\item{perplexity}{Perplexity parameter for tSNE and UMAP}

\item{dims}{Output dimensions for tSNE/UMAP (2 or 3)}

\item{ntop}{Number of top variable rows to use}

\item{npca}{Number of PCA components for initial reduction}

\item{prefix}{Prefix for cluster variable names}

\item{row.center}{Center rows before clustering?}

\item{row.scale}{Scale rows before clustering?}

\item{find.clusters}{Detect optimal number of clusters?}

\item{kclust}{Number of clusters for k-means}

\item{clust.detect}{Clustering method for optimal clusters (louvain, hclust)}

\item{method}{Embedding method (tsne, umap, pca)}
}
\value{
Updated matrix with cluster assignments
}
\description{
Clusters the rows of a matrix using dimensionality reduction and clustering.
Designed for large matrices where tSNE/UMAP cannot be run directly.
}
\details{
This function takes a matrix and clusters the rows using a combination of
PCA, tSNE/UMAP, and k-means clustering.

PCA is first used to reduce the dimensions if the matrix is large. The top
variable rows are selected before reduction. tSNE or UMAP is then applied
to further reduce the dimensions to 2 or 3. k-means clustering finally
groups the rows into clusters.

Multiple parameters control perplexity, centering/scaling, number of clusters, etc.
By default it tries to detect the optimal number of clusters.
}
