#' Example pgx object
#'
#' This is an example pgx object already created. It comes from
#' Geiger2016-arginine.
#'
#' @format ## `GEIGER_PGX`
#' A pgx object
#' \describe{
#'   \item{name}{PGX name}
#'   \item{description}{PGX description}
#'   \item{date}{PGX data created}
#'   ...
#' }
#' @source Geiger2016
"GEIGER_PGX"


#' Get path to omp example dataset(s)
#'
#' `playbase` comes bundled with a number of sample files in its `inst/extdata`
#' directory. This function make them easy to access. This function was
#' taken from tidyverse/readr.
#'
#' @param file string. Name of file. If `NULL`, the example files will
#'   be listed.
#' @export
#' @examples
#' example_file()
#' example_file("counts.csv")
example_file <- function(file = NULL) {
  if (is.null(file)) {
    dir(system.file("extdata", package = "playbase"))
  } else {
    system.file("extdata", file, package = "playbase", mustWork = TRUE)
  }
}
