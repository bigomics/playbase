% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pgx-compute.R
\name{pgx.createPGX}
\alias{pgx.createPGX}
\title{Create a PGX object
This function creates a pgx object, which is the core object in the
OmicsPlayground.}
\usage{
pgx.createPGX(
  counts,
  samples,
  contrasts,
  organism,
  custom.geneset = NULL,
  annot_table = NULL,
  max.genesets = 5000,
  name = "Data set",
  datatype = "RNA-seq",
  probe_type = NULL,
  creator = "unknown",
  description = "No description provided.",
  X = NULL,
  impX = NULL,
  norm_method = "CPM",
  is.logx = NULL,
  batch.correct = TRUE,
  auto.scale = TRUE,
  filter.genes = TRUE,
  prune.samples = FALSE,
  only.known = TRUE,
  only.hugo = TRUE,
  convert.hugo = FALSE,
  only.proteincoding = TRUE,
  remove.xxl = TRUE,
  remove.outliers = TRUE,
  add.gmt = TRUE,
  settings = list()
)
}
\arguments{
\item{counts}{Matrix of count data with genes as rows and samples as columns.}

\item{samples}{Data frame containing sample information.}

\item{contrasts}{Data frame defining sample contrasts.}

\item{organism}{Default "Human", it indicates the species used
for the gene annotation table and the probe to symbol conversion.}

\item{custom.geneset}{Custom gene sets to test, as a named list with gmt and info elements.}

\item{max.genesets}{Maximum number of gene sets to test. Default is 5000.}

\item{X}{(Optional) Matrix of normalized expression data. If NULL, will be calculated from counts.}

\item{is.logx}{Logical indicating if count matrix is already log-transformed. If NULL, guessed automatically.}

\item{batch.correct}{Logical indicating whether to perform batch correction. Default is TRUE.}

\item{auto.scale}{Logical indicating whether to automatically scale/center genes. Default is TRUE.}

\item{filter.genes}{Logical indicating whether to filter lowly expressed genes. Default is TRUE.}

\item{prune.samples}{Logical indicating whether to remove samples without contrasts. Default is FALSE.}

\item{only.known}{Logical indicating whether to keep only known genes. Default is TRUE.}

\item{only.hugo}{Logical indicating whether to convert symbols to HUGO names. Default is TRUE.}

\item{convert.hugo}{Logical indicating whether to convert symbols to HUGO names. Default is TRUE.}

\item{only.proteincoding}{Logical indicating whether to keep only protein-coding genes. Default is TRUE.}
}
\value{
List. PGX object containing input data and parameters.
}
\description{
Create a PGX object
This function creates a pgx object, which is the core object in the
OmicsPlayground.
}
\details{
pgx.createPGX creates a pgx object with the following slots:

- `name`: Name of the dataset
- `organism`: Organism for the dataset
- `version`: Dataset version
- `date`: Date the dataset was created
- `creator`: Creator of the dataset
- `datatype`: Type of data (e.g. RNA-seq, microarray)
- `description`: Description of the dataset
- `samples`: Sample metadata
- `counts`: Raw count matrix
- `contrasts`: Contrast matrix
- `X`: Normalized expression matrix
- `total_counts`: Total counts per sample
- `counts_multiplier`: Counts multiplier for each sample
- `genes`: Gene annotation data.frame (initially NULL)
- `all_genes`: Full list of genes
- `probe_type`: Probe type according to biomaRt classification(e.g. ensemble_id)
- `GMT`: Gene set matrix
}
