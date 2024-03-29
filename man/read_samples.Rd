% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pgx-read.R
\name{read_samples}
\alias{read_samples}
\title{Read samples data from file}
\usage{
read_samples(file)
}
\arguments{
\item{file}{Path to input sample data file. Should be a matrix with samples as rows and metadata as columns.}
}
\value{
dataframe the file with the data
}
\description{
Read samples data from file
}
\details{
This function reads the sample matrix with the meta-data information
of the counts and converts it to a dataframe.
}
\examples{
samples <- read_samples(playbase::example_file("samples.csv"))
}
