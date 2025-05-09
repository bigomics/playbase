% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pgx-plotting.R
\name{pgx.scatterPlotXY.GGPLOT}
\alias{pgx.scatterPlotXY.GGPLOT}
\title{Scatter plot with ggplot2}
\usage{
pgx.scatterPlotXY.GGPLOT(
  pos,
  var = NULL,
  type = NULL,
  col = NULL,
  cex = NULL,
  cex.lab = 0.8,
  cex.title = 1.2,
  cex.clust = 1.5,
  cex.legend = 1,
  cex.axis = 1,
  gridcolor = NULL,
  bgcolor = NULL,
  zoom = 1,
  legend = TRUE,
  bty = "n",
  hilight = NULL,
  zlim = NULL,
  cmin = NULL,
  cmax = NULL,
  zlog = FALSE,
  softmax = FALSE,
  zsym = FALSE,
  xlab = NULL,
  ylab = NULL,
  xlim = NULL,
  ylim = NULL,
  hilight2 = hilight,
  hilight.col = NULL,
  hilight.lwd = 0.3,
  hilight2.lwd = 1,
  hilight.cex = NULL,
  na.color = "#AAAAAA55",
  opacity = 1,
  label.clusters = FALSE,
  labels = NULL,
  legend.ysp = 0.85,
  legend.pos = "bottomright",
  tooltip = NULL,
  theme = NULL,
  set.par = TRUE,
  label.type = c("text", "box"),
  base_size = 11,
  title = NULL,
  barscale = 0.8,
  axis = TRUE,
  box = TRUE,
  guide = "legend",
  girafe = FALSE
)
}
\arguments{
\item{pos}{Data frame containing x and y coordinates.}

\item{var}{Grouping variable.}

\item{type}{Variable type, either "continuous" or "factor".}

\item{col}{Colors.}

\item{cex}{Point expansion factor.}

\item{cex.lab}{Label size factor.}

\item{cex.title}{Title size factor.}

\item{cex.clust}{Cluster label size factor.}

\item{cex.legend}{Legend text size factor.}

\item{cex.axis}{Axis text size factor.}

\item{gridcolor}{Grid line color.}

\item{bgcolor}{Background color.}

\item{zoom}{Zoom factor.}

\item{legend}{Show legend? Logical.}

\item{bty}{Legend box type.}

\item{hilight}{Indices of points to highlight.}

\item{zlim}{Color scale limits.}

\item{cmin, cmax}{Min and max color values.}

\item{zlog}{Log scale z axis? Logical.}

\item{softmax}{Apply softmax to z axis? Logical.}

\item{zsym}{Symmetrize z axis? Logical.}

\item{xlab, ylab}{Axis labels.}

\item{xlim, ylim}{Axis limits.}

\item{hilight2}{Secondary highlight indices.}

\item{hilight.col}{Color for highlighted points.}

\item{hilight.lwd}{Line width for highlighted points.}

\item{hilight.cex}{Size for highlighted points.}

\item{opacity}{Opacity factor, between 0 and 1.}

\item{label.clusters}{Logical, label clusters?}

\item{labels}{Point labels.}

\item{legend.ysp, legend.pos}{Legend position controls.}

\item{tooltip}{Tooltip text.}

\item{theme}{Theme parameters.}

\item{set.par}{Reset par? Logical.}

\item{label.type}{Label type, "text" or "box".}

\item{base_size}{Base point size.}

\item{title}{Plot title.}

\item{barscale}{Bar width scaling for box labels.}

\item{axis, box, guide}{Axis, box, and legend controls.}
}
\value{
A ggplot2 scatterplot object.
}
\description{
Scatter plot with ggplot2
}
\examples{
\dontrun{
x <- stats::rnorm(100)
y <- stats::rnorm(100)
df <- data.frame(x = x, y = y)
p <- pgx.scatterPlotXY.GGPLOT(df)
}
}
