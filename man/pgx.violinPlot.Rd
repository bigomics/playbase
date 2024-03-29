% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pgx-plotting.R
\name{pgx.violinPlot}
\alias{pgx.violinPlot}
\title{Violin plot with scatter dots}
\usage{
pgx.violinPlot(
  x,
  y,
  group = NULL,
  xlab = "",
  ylab = "",
  srt = 0,
  cex.lab = 1,
  cex.main = 1.1,
  jitter = 0.015,
  vcol = "grey85",
  main = NULL,
  plotlib = "base",
  maxbee = NULL,
  ...
)
}
\arguments{
\item{x}{Factor for x-axis groups}

\item{y}{Numeric vector of values to plot}

\item{group}{Optional grouping variable}

\item{xlab}{Label for x-axis}

\item{ylab}{Label for y-axis}

\item{srt}{Rotation angle for x-axis labels}

\item{cex.lab}{Expansion factor for axis labels}

\item{cex.main}{Expansion factor for main title}

\item{jitter}{Jittering factor for scatter dots}

\item{vcol}{Violin fill color}

\item{plotlib}{Plotting library to use (base, ggplot, plotly)}

\item{maxbee}{Maximum number of dots to show in violins}

\item{...}{Other arguments passed to plotting functions}
}
\value{
A violin plot
}
\description{
Generate a violin plot visualization with embedded scatter dots
}
\details{
This function takes a factor \code{x} and numeric vector \code{y} and generates a violin plot.
An optional \code{group} variable can be used to produce grouped/facetted violins.
Scatter dots are added to each violin with jittering controlled by \code{jitter}.
The \code{maxbee} parameter limits the number of dots.
The plot can be rendered using base, ggplot2 or plotly graphics depending on \code{plotlib}.
}
