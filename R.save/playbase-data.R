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


#' Some kind of sparse matrix dataset?
#'
#' No idea what this data is but it is used in `test_genesets`.
#'
#' @format ## `GSET_SPARSEG_XL`
#' A dgCMatrix / Matrix object
#' @source unknown
"GSET_SPARSEG_XL"


#' A dataset used by `compute_deconvolution`
#'
#' No idea what this data is but it is used in `compute_deconvolution`.
#'
#' @format ## `CCLE_RNA_CANCERTYPE`
#' @source unknown
"CCLE_RNA_CANCERTYPE"

#' A dataset used by `compute_deconvolution`
#'
#' No idea what this data is but it is used in `compute_deconvolution`.
#'
#' @format ## `CCLE_RNA_CELLINE`
#' @source unknown
"CCLE_RNA_CELLINE"

#' A dataset used by `compute_deconvolution`
#'
#' No idea what this data is but it is used in `compute_deconvolution`.
#'
#' @format ## `DICE_SIGNATURE1000`
#' @source unknown
"DICE_SIGNATURE1000"

#' A dataset used by `compute_deconvolution`
#'
#' No idea what this data is but it is used in `compute_deconvolution`.
#'
#' @format ## `GTEX_RNA_TISSUE_TPM`
#' @source unknown
"GTEX_RNA_TISSUE_TPM"

#' A dataset used by `compute_deconvolution`
#'
#' No idea what this data is but it is used in `compute_deconvolution`.
#'
#' @format ## `HPA_RNA_CELLINE`
#' @source unknown
"HPA_RNA_CELLINE"

#' A dataset used by `compute_deconvolution`
#'
#' No idea what this data is but it is used in `compute_deconvolution`.
#'
#' @format ## `IMMPROT_SIGNATURE1000`
#' @source unknown
"IMMPROT_SIGNATURE1000"

#' A dataset used by `compute_deconvolution`
#'
#' No idea what this data is but it is used in `compute_deconvolution`.
#'
#' @format ## `IMMUNOSTATES_MATRIX`
#' @source unknown
"IMMUNOSTATES_MATRIX"

#' A dataset used by `compute_deconvolution`
#'
#' No idea what this data is but it is used in `compute_deconvolution`.
#'
#' @format ## `LM22`
#' @source unknown
"LM22"

#' A dataset used by `compute_deconvolution`
#'
#' No idea what this data is but it is used in `compute_deconvolution`.
#'
#' @format ## `RNA_TISSUE_MATRIX`
#' @source unknown
"RNA_TISSUE_MATRIX"

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
