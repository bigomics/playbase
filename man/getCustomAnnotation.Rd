% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pgx-annot.R
\name{getCustomAnnotation}
\alias{getCustomAnnotation}
\title{Custom Gene Annotation}
\usage{
getCustomAnnotation(probes, custom_annot)
}
\arguments{
\item{custom_annot}{data.frame with custom annotation data. If provided,
it has to contain at least the columns "feature", "symbol", "gene_name". Also,
the features has to match the rownames of the counts provided.}

\item{counts}{A counts matrix}
}
\value{
The pgx object with custom gene annotation added/appended. The gene annotation
table has the same format as the one returned by pgx.gene_table(). However, the
columns human_ortholog, gene_title, gene_biotype, chr, pos, tx_len, map, source are filled
with default values.
}
\description{
Adds custom gene annotation table to a pgx object
}
\details{
This function allows adding a gene annotation data.frame to a pgx object when
the user has not provided an organism or it's not known.  The custom_annot data.frame
should contain gene IDs that match the pgx object genes, plus any additional columns
of annotation data.

The id_type parameter specifies the type of ID used in custom_annot to match genes.
Possible options are "symbol", "ensembl_gene_id", etc. By default it will try to match
on the "symbol" field.

Any columns in custom_annot that match existing pgx gene annotation columns will
overwrite the original data. New columns will be appended.
}
\examples{
\dontrun{
custom_annot <- data.frame(
  feature = c("A1", "A2", "A3"),
  symbol = c("TP53", "MYC", "EGFR"),
  gene_name = c("A1", "A2", "A3")
)

pgx <- getCustomAnnotation(counts, custom_annot)
}
}
