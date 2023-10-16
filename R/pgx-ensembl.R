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
                            "ensembl_peptide_id",
                            "external_gene_name",
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
#'    \item \code{feature}: The probe identifier.
#'   \item \code{gene_name}: HUman readable gene name.
#'   \item \code{human_homolog}: Gene symbol for human. Only present if working with non-human dataset.
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
                                  mart = NULL,
                                  verbose = TRUE) {
  # Check mart
  if (is.null(mart)) {
    stop("Mart not found. Please specify a BioMart database to use.")
  }

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
    values = clean_probes,
    mart = mart
  )
  annot <- data.table::data.table(annot)

  # Get homologs if working with non-human dataset
  # This should come as separate call because attr belong to diff. page
  if (!mart@dataset == "hsapiens_gene_ensembl") {
    annot_homologs <- biomaRt::getBM(
      attributes = c(probe_type, "hsapiens_homolog_associated_gene_name"),
      filters = probe_type,
      values = clean_probes,
      mart = mart
    )
    annot_homologs <- data.table::data.table(annot_homologs)
    data.table::setnames(annot_homologs, 
                        old = "hsapiens_homolog_associated_gene_name", 
                        new = "human_ortholog")
    annot <- annot[annot_homologs, on = probe_type]
  }

  # Join with clean_probes vector
  out <- annot[clean_probes, on = probe_type, mult = "first"]

  # Renaming for backwards compatibility
  if (probe_type != "external_gene_name") {
    new_names <- c("feature", 
                  "symbol",
                  "gene_title",
                  "gene_biotype",
                  "chr", 
                  "pos",
                  "tx_len", 
                  "map")
  } else {
    new_names <- c("feature", 
                  "gene_title",
                  "gene_biotype",
                  "chr", 
                  "pos",
                  "tx_len", 
                  "map")
    out[, symbol := external_gene_name]
  }
  data.table::setnames(out, old = attr_call, new = new_names)
  
  # Reorder columns and rows
  if ("human_ortholog" %chin% colnames(out)) {
    col_order <- c("feature", 
                   "symbol", 
                   "human_ortholog",
                   "gene_title",
                   "gene_biotype")
  } else {
    col_order <- c("feature", 
                   "symbol", 
                   "gene_title",
                   "gene_biotype")
  }
  data.table::setcolorder(out, col_order)
  data.table::setkeyv(out, "feature")
  
  # Keep it for back compatibility
  out[, gene_name := feature]

  # Return as data.frame and in the same order as input probes
  out <- as.data.frame(out)
  rownames(out) <- out$feature
  out <- out[probes, , drop = FALSE]
  return(out)
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
#' @import data.table
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
#' @import data.table
#' @export
probe2symbol <- function(probes, annot_table, query = "symbol", fill_na = FALSE) {

  # Prepare inputs
  query_col <- annot_table[probes, query]

  # Deal with NA
  if (fill_na) {
    query_col <- ifelse(is.na(query_col)|query_col ==  "",
          query_col,
          probes)
  }

  # Return queryed col
  return(query_col)
}


#' Retrieve gene annotation table
#'
#' @description Retrieves a gene annotation table for the given organism 
#' from Ensembl using biomaRt. Adds the table to the PGX object.
#' 
#' @param pgx PGX object with a counts table.
#' @param organism Char. Organism name. For more info see \code{\link{playbase::SPECIES_TABLE}}. 
#'
#' @return Updated PGX object with gene annotation table
#' 
#'
#' @details Queries the Ensembl database to get a gene annotation table 
#' containing external gene IDs mapped to Ensembl IDs. Handles retries in case
#' of temporary Ensembl API errors.
#'
#'
#' @examples
#' \dontrun{
#' pgx <- list()
#' pgx$counts <- matrix(rnorm(4), nrow = 2)
#' rownames(pgx$counts) <- c("ENSG00000142192", "ENSG00000288602")
#' pgx <- pgx.gene_table(pgx, "Human")
#' }
#' @export
pgx.gene_table <- function(pgx, organism) {

  # Safety checks
  stopifnot(is.list(pgx))
  stopifnot(is.character(organism))
  
  # Init vals
  genes <- NULL
  counter <- 0
  counts <- pgx$counts
  probes <- rownames(counts)
  probe_type <- NA_character_
  species_info <- playbase::SPECIES_TABLE[species_name == organism]

  # Some species appear in more than one mart, select ensembl only to avoid confusion
  if (nrow(species_info) > 1) {
    species_info <- species_info[mart == "ensembl"]
    species_info <- species_info[1, ]
  }

  # Use while loop for retries
  while (is.null(genes) && counter <= 5) {
    
    # Set waiter so that we can make multiple calls with waiting time
    Sys.sleep(counter * 60)
    
    # lock ensembl to version 110 (latest) and genes dataset
    if (species_info$mart == "ensembl") {
      ensembl <- biomaRt::useEnsembl(biomart = "genes", host = species_info$host, version = species_info$version)
      # lock ensembl to species
      ensembl <- biomaRt::useDataset(dataset = species_info$dataset, mart = ensembl)
      
    } else {
      ensembl <- biomaRt::useEnsemblGenomes(
        biomart = species_info$mart,
        dataset = species_info$dataset)

      ensembl <- biomaRt::useDataset(dataset = species_info$dataset, mart = ensembl)
    }
    
    # Get probe type
    if (is.na(probe_type)) {
      probe_type <- detect_probe(probes, ensembl)
    }

    # Get gene table
    genes <- ngs.getGeneAnnotation(
      probes = probes,
      probe_type = probe_type,
      mart = ensembl)

    all_genes <- biomaRt::getBM(attributes = "external_gene_name", mart = ensembl)
    all_genes <- all_genes[, 1]

  }

  # Return data
  pgx$genes <- genes
  pgx$all_genes <- all_genes
  pgx$probe_type <- probe_type

  return(pgx)
}
