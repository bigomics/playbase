% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pgx-deconv.R
\name{pgx.inferGender}
\alias{pgx.inferGender}
\title{Infer sample gender from gene expression}
\usage{
pgx.inferGender(X, gene_name = NULL)
}
\arguments{
\item{expr}{Numeric gene expression matrix with genes in rows.}

\item{method}{Method for gender inference ("pca" or "ttest").}

\item{pc}{Number of principal components to use if method="pca".}
}
\value{
Vector of predicted gender labels ("F" or "M") for each sample.
}
\description{
Infers sample gender based on gene expression of Y chromosome genes.
}
\details{
This function takes a gene expression matrix and infers the gender of each
sample using either a t-test (default) or PCA approach. The t-test compares expression
of Y chromosome genes between groups of known male and female samples. The PCA approach
uses the first few principal components of Y gene expression as a gender signature.

The output is a vector of predicted gender labels "F" or "M" for each sample.
}
