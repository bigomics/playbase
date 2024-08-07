% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gx-heatmap.r
\name{gx.heatmap}
\alias{gx.heatmap}
\title{Create a heatmap of a gene expression matrix}
\usage{
gx.heatmap(
  gx,
  values = NULL,
  clust.method = "ward.D2",
  dist.method = "pearson",
  col.dist.method = "euclidean",
  plot.method = "heatmap.2",
  col = (grDevices::colorRampPalette(c("royalblue3", "grey90", "indianred3")))(64),
  softmax = FALSE,
  scale = "row",
  verbose = 1,
  symm = FALSE,
  col.annot = NULL,
  row.annot = NULL,
  annot.ht = 1,
  annot.cex = 1,
  nmax = 1000,
  cmax = NULL,
  show_colnames = TRUE,
  indent.names = FALSE,
  ...
)
}
\arguments{
\item{gx}{Gene expression matrix with genes in rows and samples in columns}

\item{values}{value}

\item{clust.method}{Hierarchical clustering method for rows and columns}

\item{dist.method}{Distance metric for clustering genes}

\item{col.dist.method}{Distance metric for clustering samples}

\item{plot.method}{Heatmap plotting method (heatmap.2 or ComplexHeatmap)}

\item{col}{Colors to use for heatmap}

\item{softmax}{value}

\item{scale}{Scaling method for the data matrix}

\item{verbose}{value}

\item{symm}{value}

\item{col.annot}{Annotations for samples}

\item{row.annot}{Annotations for genes}

\item{annot.ht}{value}

\item{nmax}{Maximum number of genes to include}

\item{cmax}{value}

\item{show_colnames}{value}

\item{indent.names}{value}

\item{...}{Additional arguments passed to heatmap plotting function}
}
\value{
A heatmap grob object
}
\description{
Create a heatmap of a gene expression matrix
}
\details{
This function subsets the input matrix to the top variable genes,
clusters genes and samples using the specified methods, and visualizes the
result as a heatmap using heatmap.2 or ComplexHeatmap.

Gene and sample clustering is optional. The data can be scaled by rows or
columns before clustering. A color palette is used to map values to colors.
}
\examples{
\dontrun{
set.seed(123)
gx <- matrix(rnorm(250)^2, nrow = 25, ncol = 10)
rownames(gx) <- sample(LETTERS, 25)
colnames(gx) <- sample(letters, 10)
p <- gx.heatmap(gx)
}
}
