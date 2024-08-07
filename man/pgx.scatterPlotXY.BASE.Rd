% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pgx-plotting.R
\name{pgx.scatterPlotXY.BASE}
\alias{pgx.scatterPlotXY.BASE}
\title{Scatter plot with base graphics}
\usage{
pgx.scatterPlotXY.BASE(
  pos,
  var = NULL,
  type = NULL,
  col = NULL,
  title = "",
  zlim = NULL,
  zlog = FALSE,
  zsym = FALSE,
  softmax = FALSE,
  pch = 20,
  cex = NULL,
  cex.lab = 1,
  cex.title = 1.2,
  cex.legend = 1,
  zoom = 1,
  legend = TRUE,
  bty = "o",
  legend.ysp = 0.85,
  legend.pos = "bottomright",
  lab.pos = NULL,
  repel = TRUE,
  xlab = NULL,
  ylab = NULL,
  xlim = NULL,
  ylim = NULL,
  dlim = 0.02,
  hilight2 = hilight,
  hilight.cex = NULL,
  lab.xpd = TRUE,
  hilight = NULL,
  hilight.col = NULL,
  hilight.lwd = 0.8,
  label.clusters = FALSE,
  cex.clust = 1.5,
  tstep = 0.1,
  rstep = 0.1,
  na.color = "#AAAAAA55",
  tooltip = NULL,
  theme = NULL,
  set.par = TRUE,
  axt = "s",
  xaxs = TRUE,
  yaxs = TRUE,
  labels = NULL,
  label.type = NULL,
  opacity = 1
)
}
\arguments{
\item{var}{Numeric vector of values to plot on x/y axes.}

\item{type}{Variable type, either "continuous" or "factor".}

\item{col}{Colors for points.}

\item{title}{Plot title.}

\item{zlim}{Value limits for coloring points by zvar.}

\item{zlog}{Log scale zvar colors?}

\item{zsym}{Symmetrize zvar colors?}

\item{softmax}{Softmax transform zvar colors?}

\item{pch}{Point shapes.}

\item{cex}{Point expansion factor.}

\item{cex.lab}{Label size.}

\item{cex.title}{Title size.}

\item{cex.legend}{Legend size.}

\item{zoom}{Zoom factor.}

\item{legend}{Show legend?}

\item{bty}{Legend box type.}

\item{legend.ysp}{Legend y-spacing.}

\item{legend.pos}{Legend position.}

\item{lab.pos}{Label positions, matrix with x, y.}

\item{repel}{Repel overlapping labels?}

\item{xlab, ylab}{Axis labels.}

\item{xlim, ylim}{Axis limits.}

\item{dlim}{Axis limit expansion fraction.}

\item{hilight2}{Secondary highlight indices.}

\item{hilight.cex}{Size for highlight points.}

\item{lab.xpd}{Allow labels outside plot area?}

\item{hilight}{Indices to highlight.}

\item{hilight.col}{Color for highlighted points.}

\item{hilight.lwd}{Line width for highlighted points.}

\item{label.clusters}{Label clusters?}

\item{cex.clust}{Cluster label size.}

\item{tstep, rstep}{Label repel increments.}

\item{tooltip}{Tooltip text.}

\item{theme}{Theme parameters.}

\item{set.par}{Reset par?}

\item{axt, xaxs, yaxs}{Axis type, style, limits.}

\item{labels}{Point labels.}

\item{label.type}{Label type.}

\item{opacity}{Opacity for all elements.}
}
\value{
None. Plot is produced as a side-effect.
}
\description{
Scatter plot with base graphics
}
