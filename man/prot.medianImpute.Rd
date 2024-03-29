% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pgx-proteomics.R
\name{prot.medianImpute}
\alias{prot.medianImpute}
\title{Median Imputation of Proteomics Data}
\usage{
prot.medianImpute(X, groups)
}
\arguments{
\item{X}{A numeric matrix of proteomics data, where rows represent proteins and columns represent samples.}

\item{groups}{A factor or character vector indicating the group membership of each sample (column) in X.}
}
\value{
A numeric matrix of the same dimensions as X, with missing and zero values replaced by the median of their respective group.
}
\description{
This function performs median imputation on proteomics data, replacing missing or zero values with the median of the non-missing values in the same group.
}
\details{
The function first replaces any missing values in X with zeros.
Then, for each protein (row), the function calculates the median of the non-missing values in each group, using a log-transformation to handle skewed data.
The function then replaces any zero values in X with the calculated median for their respective group.
}
