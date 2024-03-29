% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ui-colors.R
\name{omics_pal_d}
\alias{omics_pal_d}
\title{Return function to interpolate a discrete omics color palette}
\usage{
omics_pal_d(palette = "default", reverse = FALSE)
}
\arguments{
\item{palette}{Character name of palette in omics_palettes}

\item{reverse}{Boolean indicating whether the palette should be reversed}
}
\description{
Return function to interpolate a discrete omics color palette
}
\examples{
omics_pal_d()
omics_pal_d("light")
omics_pal_d("mUTED", reverse = TRUE)

}
