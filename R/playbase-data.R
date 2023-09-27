#' Get path to omp example dataset(s)
#'
#' `playbase` comes bundled with a number of sample files in its `inst/extdata`
#' directory. This function make them easy to access. This function was
#' taken from tidyverse/readr.
#'
#' @param file string. Name of file. If `NULL`, the example files will
#'   be listed.

#' @examples
#' example_file()
#' example_file("counts.csv")
#' @export
example_file <- function(file = NULL) {
  if (is.null(file)) {
    dir(system.file("extdata", package = "playbase"))
  } else {
    system.file("extdata", file, package = "playbase", mustWork = TRUE)
  }
}

#' Example sample data from unknown GEO dataset
#'
#' @format ## `SAMPLES`
#' data.frame
#' @source unknown
"SAMPLES"

#' Example contrasts data from unknown GEO dataset
#'
#' @format ## `CONTRASTS`
#' data.frame with rows and columns as contrasts, values are -1, 0 or 1.
#' @source unknown
"CONTRASTS"

#' Example counts data from unknown GEO dataset
#'
#' @format ## `COUNTS`
#' data.frame with genes as rows and samples as columns
#' @source unknown
"COUNTS"

#' Checks performed by the pgx.CHECK function
#'
#' @format ## `data.frame`
#' rows are checks, columns are description of the check performed.
"PGX_CHECKS"


#' Example GMT (geneset matrix transpose) of genes targeted by microRNA
#'
#' @format ## `GMT`
#' list of genes (targeted by a microRNA), and the list name is the microRNA
#' @source C3: regulatory target gene sets, MIR: microRNA targets from https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
"EXAMPLE_GMT"

#' Table for species ID conversin conversion from scientific name to Ensembl db ID
#'
#' @format ## `SPECIES_TABLE`
#' data.frame with columns "dataset"	"description"	"version"	"species_name"
#' @source BiomaRt
"EXAMPLE_GMT"
