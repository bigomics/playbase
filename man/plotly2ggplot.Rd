% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pgx-plotting.R
\name{plotly2ggplot}
\alias{plotly2ggplot}
\title{Convert a plotly plot to a ggplot object}
\usage{
plotly2ggplot(
  plot,
  width = NULL,
  height = NULL,
  scale = 1,
  hjust = 0,
  vjust = 0
)
}
\arguments{
\item{plot}{The plotly plot object to convert}

\item{width}{The width of the resulting ggplot, in pixels. Default NULL.}

\item{height}{The height of the resulting ggplot, in pixels. Default NULL.}

\item{scale}{The scaling factor to apply. Default 1.}

\item{hjust}{The horizontal adjustment to apply. Default 0.}

\item{vjust}{The vertical adjustment to apply. Default 0.}
}
\value{
A ggplot object representing the converted plotly plot.
}
\description{
Convert a plotly plot to a ggplot object
}
\details{
This function takes a plotly plot object and converts it to a ggplot object.
The plotly plot is rendered to a temporary PNG file using the orca() function. This PNG
is then read in and added as an annotation to a blank ggplot. Optional parameters allow
adjusting the width, height, scaling, and justification of the resulting ggplot.
}
