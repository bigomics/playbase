% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ui-colors.R
\name{scale_fill_omics_d}
\alias{scale_fill_omics_d}
\title{Fill scale constructor for discrete omics color palettes}
\usage{
scale_fill_omics_d(palette = "default", reverse = FALSE, ...)
}
\arguments{
\item{palette}{Character name of palette in omics_palettes}

\item{reverse}{Boolean indicating whether the palette should be reversed}

\item{...}{Additional arguments passed to discrete_scale() or
scale_fill_gradientn(), used respectively when discrete is TRUE
or FALSE}
}
\description{
Fill scale constructor for discrete omics color palettes
}
\examples{
library(ggplot2)
ggplot(mpg, aes(class, fill = class)) +
  geom_bar() +
  scale_fill_omics_d()
ggplot(mpg, aes(class, fill = class)) +
  geom_bar() +
  scale_fill_omics_d("mUTED", reverse = TRUE)
ggplot(iris, aes(Species, Sepal.Width, fill = Species)) +
  geom_jitter(shape = 21, size = 4, stroke = 1) +
  scale_fill_omics_d("super_light")

}
