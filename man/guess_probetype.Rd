% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pgx-ensembl.R
\name{guess_probetype}
\alias{guess_probetype}
\title{Guess probe type}
\usage{
guess_probetype(probes, organism, for.biomart = FALSE)
}
\arguments{
\item{probes}{Character vector of probes to detect type for.}

\item{mart}{BioMart object specifying the database to use for mapping.}

\item{verbose}{Logical indicating whether to print progress messages.}
}
\value{
The probe type with the best mapping performance. One of:
\itemize{
 \item ensembl_transcript_id
 \item ensembl_transcript_id_version
 \item ensembl_gene_id
 \item ensembl_gene_id_version
 \item uniprot_gn_id
 \item refseq_peptide
 \item refseq_mrna
 \item hgnc_symbol
}
}
\description{
This function tries to automatically detect the probe type of a set
of input probes using regular exprssion or by testing mapping
success against different identifier types in a BioMart database.
}
\details{
This function subsamples a subset of the input probes and tries mapping
them to each identifier type in the specified BioMart database. It returns
the type with the maximum number of mapped probes.
}
\examples{
\dontrun{
library(playbase)
library(biomaRt)
# list databases
listEnsembl(version = 110)
# use genes database
ensembl <- useEnsembl(biomart = "genes", version = 110)
# here we see a list of 214 species available in ensembl
datasets <- listDatasets(ensembl)
# here we can select between 214 species
dataset_hsa <- searchDatasets(mart = ensembl, pattern = "hsapiens")
ensembl_species <- useDataset(dataset = dataset_hsa$dataset, mart = ensembl)
probes <- c(
  "ENSG00000230915.1", "ENSG00000275728.1", "ENSG00000277599.1",
  "ENSG00000186163.9", "ENSG00000164823.11", "ENSG00000274234.1",
  "ENSG00000282461.1", "ENSG00000283056.1", "ENSG00000239021.1",
  "ENSG00000214268.2", "ENSG00000206687.1", "ENSG00000171148.14",
  "ENSG00000250027.1", "ENSG00000244217.1", "ENSG00000103502.14",
  "ENSG00000213178.3", "ENSG00000235059.5", "ENSG00000204555.3",
  "ENSG00000221044.2", "ENSG00000267162.1"
)
type <- detect_probe(probes, mart = ensembl_species)
}

}
