% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pgx-read.R
\name{read_contrasts}
\alias{read_contrasts}
\title{Read contrasts data from file}
\usage{
read_contrasts(file)
}
\arguments{
\item{file}{string. path to file}
}
\value{
matrix. the file with the data
}
\description{
Read contrasts data from file
}
\examples{
contrasts <- read_contrasts(playbase::example_file("contrasts.csv"))
}
