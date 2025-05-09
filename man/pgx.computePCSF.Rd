% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pgx-pcsf.R
\name{pgx.computePCSF}
\alias{pgx.computePCSF}
\title{Compute PCSF solution from pgx object}
\usage{
pgx.computePCSF(
  pgx,
  contrast,
  level = "gene",
  ntop = 250,
  ncomp = 3,
  beta = 1,
  rm.negedge = TRUE,
  use.corweight = TRUE,
  dir = "both",
  ppi = c("STRING", "GRAPHITE"),
  gset.rho = 0.8,
  gset.filter = NULL,
  as.name = NULL
)
}
\description{
Compute PCSF solution from pgx object for given contrast.
}
