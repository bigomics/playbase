% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ui-colors.R
\name{scale_fill_omics_c}
\alias{scale_fill_omics_c}
\title{Fill scale constructor for continuous omics color palettes}
\usage{
scale_fill_omics_c(palette = "brand_blue", reverse = FALSE, ...)
}
\arguments{
\item{palette}{Character name of palette in omics_palettes}

\item{reverse}{Boolean indicating whether the palette should be reversed}

\item{...}{Additional arguments passed to discrete_scale() or
scale_fill_gradientn(), used respectively when discrete is TRUE
or FALSE}
}
\description{
Fill scale constructor for continuous omics color palettes
}
\examples{
library(ggplot2)
ggplot(mpg, aes(displ, cty, fill = hwy)) +
  geom_point(shape = 21, size = 4, stroke = 1) +
  scale_fill_omics_c()
ggplot(iris, aes(Sepal.Width, Sepal.Length, fill = Petal.Width)) +
  geom_point(shape = 21, size = 4, stroke = 1) +
  scale_fill_omics_c("blue_red_grey", reverse = TRUE)

}
