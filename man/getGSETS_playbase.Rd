% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pgx-functions.R
\name{getGSETS_playbase}
\alias{getGSETS_playbase}
\title{Get gene sets from playbase data}
\usage{
getGSETS_playbase(gsets = NULL, pattern = NULL)
}
\arguments{
\item{pattern}{Pattern to match gene set names. Default NULL matches all sets.}

\item{sets}{Character vector of gene set names to retrieve. Default NULL retrieves all sets.}
}
\value{
A named list containing the gene sets matching the criteria.
}
\description{
Retrieves gene sets from the playbase data package matching a pattern.
Allows optionally including additional custom gene sets.
}
\details{
This function extracts gene sets from the playbase data package.
It returns gene sets matching the provided pattern.
If a sets parameter is provided, only those specific sets are returned.
By default it returns all gene sets in playbase.
}
