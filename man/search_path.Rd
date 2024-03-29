% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pgx-functions.R
\name{search_path}
\alias{search_path}
\title{Search for file in system path}
\usage{
search_path(paths, file)
}
\arguments{
\item{file}{Character string of filename to search for}

\item{path}{Character vector of directories to search}
}
\value{
Full path of file if found, otherwise NULL
}
\description{
Search for file in system path
}
\details{
Searches the directories in \code{path} to find the file specified by \code{file}.
Returns the full path if the file is found in one of the directories.
Returns NULL if the file is not found.
}
