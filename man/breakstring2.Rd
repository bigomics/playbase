% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pgx-functions.R
\name{breakstring}
\alias{breakstring}
\alias{breakstring2}
\title{Break up a long string into multiple lines}
\usage{
breakstring(s, n, nmax = 999, force = FALSE, brk = "\\n")

breakstring2(s, n, brk = "\\n", nmax = 999)
}
\arguments{
\item{s}{The input string}

\item{n}{The maximum line length}

\item{nmax}{The maximum length to process. Longer strings are truncated. Default is 999.}

\item{force}{logical not used}

\item{brk}{The string to insert for line breaks. Default is "\\n".}
}
\value{
A string with line breaks inserted.
}
\description{
Break up a long string into multiple lines
}
\details{
This function takes a long string \code{s} and inserts line breaks \code{brk} every \code{n} characters.
It first truncates the input string to a maximum length \code{nmax} if longer.
The string is then split into words, and words are added to each line until it reaches length \code{n}.
A line break \code{brk} is inserted, and the process continues until all words have been processed.
}
\section{Functions}{
\itemize{
\item \code{breakstring()}: Breaks a character string into substrings of a
specified length, separated by a specified line break character.

}}
