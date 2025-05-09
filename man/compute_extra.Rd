% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute2-extra.R
\name{compute_extra}
\alias{compute_extra}
\title{Compute Extra Analysis}
\usage{
compute_extra(
  pgx,
  extra = c("meta.go", "infer", "deconv", "drugs", "connectivity", "wordcloud", "wgcna",
    "mofa"),
  sigdb = NULL,
  pgx.dir = "./data",
  libx.dir = "./libx",
  user_input_dir = getwd()
)
}
\arguments{
\item{pgx}{An object containing the input data for analysis.}

\item{extra}{A character vector specifying which additional analyses to perform.}

\item{sigdb}{A character vector specifying the path to the sigdb-.h5 files for connectivity scores.}

\item{libx.dir}{The directory where the sigdb-.h5 files are located.}
}
\value{
An updated object with additional analysis results.
}
\description{
This function computes additional analysis based on the input data, such as GO core graph,
deconvolution, phenotype inference, drug activity enrichment, OmicsGraphs, WordCloud statistics,
connectivity scores, and wgcna.
}
