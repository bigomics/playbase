% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pgx-init.R
\name{pgx.initialize}
\alias{pgx.initialize}
\title{Initialize and validate a PGX object}
\usage{
pgx.initialize(pgx)
}
\arguments{
\item{pgx}{A PGX object to initialize.}
}
\value{
The initialized PGX object.
}
\description{
Validates and initializes a PGX object by performing checks,
conversions, and adding default parameters.
}
\details{
This function performs several validation checks on the PGX object:

- Checks that counts and samples data frames are valid
- Checks that model parameters are defined
- Converts counts to matrix
- Converts model matrix to matrix
- Defines group labels if not provided
- Converts contrasts to labeled matrix form

It will throw errors if required components like counts, samples, or groups are missing.
Additional default parameters like model formula and contrasts are added if not already defined.
}
\examples{
\dontrun{
# TODO
}
}
