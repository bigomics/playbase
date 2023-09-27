##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

#' Detect probe type
#'
#' This function tries to automatically detect the probe type of a set of
#' input probes by testing mapping success against different identifier
#' types in a BioMart database.
#'
#' @param probes Character vector of probes to detect type for.
#' @param mart BioMart object specifying the database to use for mapping.
#' @param verbose Logical indicating whether to print progress messages.
#'
#' @return The probe type with the best mapping performance. One of:
#' \itemize{
#'  \item ensembl_transcript_id
#'  \item ensembl_transcript_id_version
#'  \item ensembl_gene_id
#'  \item ensembl_gene_id_version
#'  \item uniprot_gn_id
#'  \item refseq_peptide
#'  \item refseq_mrna
#'  \item hgnc_symbol
#' }
#'
#' @details
#' This function subsamples a subset of the input probes and tries mapping
#' them to each identifier type in the specified BioMart database. It returns
#' the type with the maximum number of mapped probes.
#'
#' @examples
#' \dontrun{
#' library(playbase)
#' library(biomaRt)
#' # list databases
#' listEnsembl(version = 110)
#'  # use genes database
#' ensembl <- useEnsembl(biomart="genes", version = 110)
#' # here we see a list of 214 species available in ensembl
#' datasets <- listDatasets(ensembl)
#' # here we can select between 214 species
#' dataset_hsa <- searchDatasets(mart = ensembl, pattern = "hsapiens")
#' ensembl_species <- useDataset(dataset = dataset_hsa$dataset, mart = ensembl)
#' probes <- c("ENSG00000230915.1", "ENSG00000275728.1", "ENSG00000277599.1",
#'  "ENSG00000186163.9", "ENSG00000164823.11", "ENSG00000274234.1",
#'  "ENSG00000282461.1", "ENSG00000283056.1", "ENSG00000239021.1",
#'  "ENSG00000214268.2", "ENSG00000206687.1", "ENSG00000171148.14",
#'  "ENSG00000250027.1", "ENSG00000244217.1", "ENSG00000103502.14",
#'  "ENSG00000213178.3", "ENSG00000235059.5", "ENSG00000204555.3",
#'  "ENSG00000221044.2", "ENSG00000267162.1")
#' type <- detect_probe(probes, mart = ensembl_species)
#' }
#'
#' @export
detect_probe <- function(probes, mart = NULL, verbose = TRUE){

  # Check mart
  if (is.null(mart)) {
    stop("Mart not found. Please specify a BioMart database to use.")
  }
  # Prepare inputs
  if (verbose) {
    message("[createPGX] Guessing probe type...")
  }
  clean_probes <- probes[!is.na(probes)]
  n <- length(clean_probes)
  # if number of probes above 10, keep only 10 random probes

  if (n > 100L) n2 <- 100L else n2 <- n
  subsample <- sample(1:n, n2)


  # Vector with input types to check
  probe_types_to_check <- c("ensembl_gene_id",
                            "ensembl_transcript_id",
                            "ensembl_transcript_id_version",
                            "hgnc_symbol",
                            "ensembl_gene_id_version",
                            "uniprot_gn_id",
                            "refseq_peptide",
                            "refseq_mrna"
                            )
  subset_probes <- clean_probes[subsample]

  # often we see multiples probes at once
  # initially we can try to detect as many probes as possible from each probe type
  # the type that guess better, is the one that will be used
  # this approach still has issues, as sometimes we have mixed probe types in on study.


  probe_check <- sapply(probe_types_to_check, FUN = function(x) {
    # TODO write a while loop that stops when we first get 20 IDs mapped, o
    # continue until the end if not
    tryCatch({
      tmp <- biomaRt::getBM(attributes = x,
                            filters = x,
                            values = subset_probes,
                            mart = mart)
      Sys.sleep(2)  # Sleep time to prevent bounce from ensembl for consecutive calls
      out <- nrow(tmp)
      return(out)
    }, error = function(e) {
      return(0)
    }
    )
  })

  # Check matches and return if winner
  if (all(probe_check == 0)) {
    #TODO the probe2symbol and detect_probe code should be used in data-preview,
    # and we should warn the user in case no matches are found
    stop("Probe type not found, please, check your probes")
  } else {
    probe_type <- names(which.max(probe_check))
  }
  return(probe_type)
}


#' Get gene annotation data
#'
#' Retrieves gene annotation information from BioMart for a set of input
#' gene/transcript identifiers.
#'
#' @param probes Character vector of gene/transcript identifiers to retrieve annotation for.
#' @param probe_type Character specifying the type of input identifiers. If NULL,
#' it will be automatically detected. Options are "ensembl_gene_id", "ensembl_transcript_id", etc.
#' @param mart BioMart object specifying the database to query.
#' @param verbose Logical indicating whether to print status messages.
#'
#' @return Data frame with gene annotation data for the input identifiers. Columns are:
#' \itemize{
#'   \item \code{gene_name}: Gene name
#'   \item \code{hgnc_symbol}: Gene symbol
#'   \item \code{gene_title}: Gene description
#'   \item \code{gene_biotype}: Gene biotype
#'   \item \code{chr}: Chromosome
#'   \item \code{pos}: Transcript start position
#'   \item \code{tx_len}: Transcript length
#'   \item \code{map}: Chromosome band
#' }
#'
#' @details This function queries BioMart to retrieve key gene annotation data for
#' a set of input gene/transcript identifiers. It can detect the identifier
#' type automatically if not provided.
#'
#'
#' @examples
#' \dontrun{
#' probes <- c("ENSG00000142192", "ENST00000288602")
#' mart <- biomaRt::useMart("ensembl")
#' result <- ngs.getGeneAnnotation(probes, mart)
#' head(result)
#' }
#' @export
ngs.getGeneAnnotation <- function(probes,
                                  probe_type = NULL,
                                  mart,
                                  verbose = TRUE) {

  # Prepare inputs
  if (verbose) {
    message("[createPGX] Filling genes information...")
  }

  clean_probes <- probes[!duplicated(probes)]

  if (is.null(probe_type)) {
    probe_type <- detect_probe(probes, mart)
  }

  # Select attributes
  attr_call <- c(
    probe_type,
    "external_gene_name",  # gene_name
    "hgnc_symbol",         # Hugo name
    "description",         # gene_title
    "gene_biotype",        # gene_biotype
    "chromosome_name",     # chr
    "transcript_start",    # pos
    "transcript_length",   # tx_len
    "band"                 # map
  )
  attr_call <- attr_call[!duplicated(attr_call)]

  # Get the gene annotations
  annot <- biomaRt::getBM(
    attributes = attr_call,
    filters = probe_type,
    values = probes,
    mart = mart
  )
  annot <- data.table::data.table(annot)
  data.table::setkeyv(annot, probe_type)
  return(annot)
}


#' Map probe identifiers to gene symbols
#'
#' This function converts a vector of probe identifiers to
#' standard HGNC gene symbols using an annotation lookup table.
#'
#' @param probes Character vector of probe IDs to convert.
#' @param annot_table Data frame with columns "probe_type" and "hgnc_symbol".
#'   The probe_type matches the type of IDs in probes.
#'
#' @return Character vector of mapped HGNC gene symbols.
#'
#' @details The annot_table should contain a column with the probe IDs
#'   (matching type of probes input) and a column with the corresponding HGNC
#'   gene symbols. This function matches the input probes to the table
#'   to retrieve the gene symbols. Unmatched probes are returned as is.
#'
#' \dontrun{
#' probes <- c("ENSG00000142192", "ENST00000288602")
#' annot_table <- data.frame(
#'   ensembl_gene_id = c("ENSG00000142192", "ENSG00000099977"),
#'   hgnc_symbol = c("EGFR", "CDKN2A")
#' )
#' symbols <- probe2symbol(probes, annot_table)
#' }
#'
#' @export
probe2symbol <- function(probes, annot_table, query = "external_gene_name") {

  # Prepare inputs
  probe_type <- colnames(annot_table[, 1])
  probe_cols <- c(probe_type, query)

  # Match annot_table
  annot <- annot_table[, .SD,
                       .SDcols = probe_cols][probes,
                                             on = probe_type,
                                             mult = "first"]

  # Deal with NA
  annot[is.na(annot[[query]])|annot[[query]] ==  "",
        (query) := .SD,
        .SDcols = probe_type]

  # Return queryed col
  return(annot[[query]])
}