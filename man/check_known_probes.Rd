% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pgx-ensembl.R
\name{check_known_probes}
\alias{check_known_probes}
\title{Known probes}
\usage{
check_known_probes(probes, probe_types_to_check = NULL)
}
\arguments{
\item{probes}{A character vector of probe identifiers.}

\item{probe_types_to_check}{A character vector of probe types to check.}
}
\value{
A character vector of probe types to check mapping for.
}
\description{
This function checks if a vector of probe identifiers matches known patterns
for common identifier types like Ensembl IDs. It returns a vector of probe
types that should be checked for mapping based on the matches.
}
\details{
This function checks `probes` for the following known patterns:

- Ensembl IDs for higher animals (ENS*)
- Drosophila melanogaster Ensembl IDs (FBgn*)
- Caenorhabditis elegans Ensembl IDs (WBGene*)

If a high proportion (>50%) of `probes` match one of these patterns,
the corresponding Ensembl ID types are returned in `probe_types_to_check`.
}
\examples{
probes <- c("ENSG000001", "ENST0001", "FBgn001", "WBGene0001")
check_known_probes(probes, NULL)

}
