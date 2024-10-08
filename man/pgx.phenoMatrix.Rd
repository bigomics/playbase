% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pgx-functions.R, R/pgx-pheno.R
\name{pgx.phenoMatrix}
\alias{pgx.phenoMatrix}
\title{Create a Phenotype Matrix}
\usage{
pgx.phenoMatrix(pgx, phenotype)

pgx.phenoMatrix(pgx, phenotype)
}
\arguments{
\item{pgx}{A list representing a PGX object containing the data to be analyzed.}

\item{phenotype}{A character string specifying the name of the phenotype to use for creating the matrix.}
}
\value{
A numeric matrix representing the phenotype matrix created from the specified phenotype.

A numeric matrix representing the phenotype matrix created from the specified phenotype.
}
\description{
This function creates a phenotype matrix from a PGX object.

This function creates a phenotype matrix from a PGX object.
}
\details{
This function takes a PGX object `pgx` and a character string `phenotype` as input and
creates a phenotype matrix using the specified phenotype. The phenotype matrix is created by extracting
the specified phenotype from the `samples` data frame in the `pgx` object and converting it into a model matrix.
The resulting phenotype matrix is returned as a numeric matrix.

This function takes a PGX object `pgx` and a character string `phenotype` as input and
creates a phenotype matrix using the specified phenotype. The phenotype matrix is created by extracting
the specified phenotype from the `samples` data frame in the `pgx` object and converting it into a model matrix.
The resulting phenotype matrix is returned as a numeric matrix.
}
