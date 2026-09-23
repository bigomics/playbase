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
"SPECIES_TABLE"


#' GENES_TABLE dataset
#'
#' A dataset containing the gene annotation table for the example dataset
#'
#' @format A data frame with 6826 rows and 10 columns:
#' \describe{
#'   \item{symbol}{Description of column1}
#'   \item{gene_title}{Description of column2}
#'   \item{chr}{Description of column2}
#'   \item{pos}{Description of column2}
#'   \item{tx_len}{Description of column2}
#'   \item{map}{Description of column2}
#'   \item{gene_name}{Description of column2}
#'   \item{feature}{Description of column2}
#' }
#' @source Proteome profiles of activated vs resting human naive T cells at different times (Geiger et al., Cell 2016).
"GENES_TABLE"


#' GENES_TABLE dataset
#'
#' A PGX object after running pgx.createPGX() on the mini example counts,
#' samples and contrast data
#'
#' @format A PGX object
#' @source Proteome profiles of activated vs resting human naive T cells at different times (Geiger et al., Cell 2016).
"PGX_CREATE"
