% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pgx-functions.R
\name{param.class}
\alias{param.class}
\title{Get Classes of Parameters}
\usage{
param.class(A)
}
\arguments{
\item{A}{A data frame representing the input data.}
}
\value{
A character vector representing the classes of the parameters in
the input data frame.
}
\description{
This function returns the classes of the parameters in a data frame.
}
\details{
The function takes a data frame `A` as input and applies the `class`
function to each column of the data frame using the `sapply` function.
The resulting vector of classes is returned.
}
\examples{
\dontrun{
# example code
A <- data.frame(x = c(1, 2, 3), y = c("a", "b", "a"))
param.class(A)
}
}
