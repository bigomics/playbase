% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pgx-functions.R
\name{shortstring}
\alias{shortstring}
\alias{shortstring0}
\title{Truncate a string to a maximum length}
\usage{
shortstring(s, n, dots = 1)

shortstring0(s, n, dots = 1)
}
\arguments{
\item{s}{The string to truncate.}

\item{n}{The maximum length of the truncated string.}

\item{dots}{The number of characters to devote to the ellipses.}
}
\value{
The truncated string.
}
\description{
Truncates a string to a maximum length by replacing the middle part with ellipses.
}
\details{
This function truncates long strings to a specified maximum length
by replacing the middle part with ellipses.

It first converts the input string to UTF-8 encoding and removes any HTML entities.
If the string is already shorter than the maximum length, it is returned unchanged.
Otherwise, it truncates the beginning and end of the string to fit within the max length when concatenated with the ellipses.

The number of characters to use for the ellipses is set by the \code{dots}
parameter. By default, it uses one dot per 10 characters of max length.
}
\section{Functions}{
\itemize{
\item \code{shortstring()}: The shortstring function truncates a string to a specified length.

}}
\examples{
\dontrun{
shortstring0("This is a very long string", 20)
# Returns "This is a ...string"

shortstring0("Short string", 20)
# Returns "Short string"
}
}
