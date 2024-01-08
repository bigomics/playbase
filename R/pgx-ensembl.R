##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

#' Retrieve gene annotation table
#'
#' @description Retrieves a gene annotation table for the given organism
#' from Ensembl using biomaRt. Adds the table to the PGX object.
#'
#' @param pgx PGX object with a counts table.
#' @param organism Char. Organism name. For more info see \code{\link{playbase::SPECIES_TABLE}}.
#' @param annot_table Custom annotation table. See \code{\link{playbase::pgx.custom_annotation}}.
#' @param use_biomart Logical. If TRUE, use biomaRt to retrieve gene annotation.
#' @return Updated PGX object with gene annotation table
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
pgx.addGeneAnnotation <- function(pgx, organism = NULL, annot_table = NULL, use_biomart = NULL) {
  # Safety checks
  stopifnot(is.list(pgx))
  probes <- rownames(pgx$counts)

  if (is.null(organism) && !is.null(pgx$organism)) {
    organism <- pgx$organism
  }
  if (is.null(organism) && !is.null(pgx$organism)) {
    organism <- guess_organism(probes)
  }
  if (is.null(organism)) {
    stop("could not determine organism. please specify")
  }
  # Get gene table
  genes <- ngs.getGeneAnnotation(
    pgx = pgx, 
    probes = probes,
    annot_table = annot_table,
    organism = organism,
    use_biomart = use_biomart ## auto-select
  )

  # Return data
  pgx$genes <- genes

  return(pgx)
}


#' Get gene annotation data
#'
#' Retrieves gene annotation information from BioMart for a set of input
#' gene/transcript identifiers.
#'
#' @param probes Character vector of gene/transcript identifiers to retrieve annotation for.
#' @param organism Organism name, e.g. "hsapiens_gene_ensembl".
#' @param probe_type Character specifying the type of input identifiers. If NULL,
#' it will be automatically detected. Options are "ensembl_gene_id", "ensembl_transcript_id", etc.
#' @param mart BioMart object specifying the database to query.
#' @param verbose Logical indicating whether to print status messages.
#'
#' @return Data frame with gene annotation data for the input identifiers. Columns are:
#' \itemize{
#'    \item \code{feature}: The probe identifier.
#'   \item \code{sybmol}: Human readable gene name.
#'   \item \code{human_homolog}: Gene symbol for human. Only present if working with non-human dataset.
#'   \item \code{gene_title}: Gene description
#'   \item \code{gene_biotype}: Gene biotype
#'   \item \code{chr}: Chromosome
#'   \item \code{pos}: Transcript start position
#'   \item \code{tx_len}: Transcript length
#'   \item \code{map}: Chromosome band
#'   \item \code{gene_name}: equivalent to the rownames. Kept for back compatibility
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
ngs.getGeneAnnotation <- function(pgx, probes, organism = NULL, annot_table = NULL, use_biomart = NULL) {
  if (is.null(organism)) {
    organism <- guess_organism(probes)
  }
  if (is.null(organism)) {
    warning("[getGeneAnnotation] ERROR could not detect organism")
    return(NULL)
  }

  message("[ngs.getGeneAnnotation] organism = ", organism)

  is.primary_organism <- (tolower(organism) %in% c("human", "mouse", "rat"))
  if (is.null(use_biomart) && is.primary_organism) {
    use_biomart <- FALSE
  }
  if (is.null(use_biomart) && !is.primary_organism) {
    use_biomart <- TRUE
  }

  genes <- NULL

  ## first try ORGDB for human, mouse, rat
  if (is.null(genes) && is.primary_organism && !use_biomart) {
    message("[getGeneAnnotation] >>> annotating genes using ORGDB libraries")
    probe_type <- guess_probetype(probes, for.biomart = FALSE)
    if (is.null(probe_type)) stop("probe_type is NULL")
    message("[ngs.getGeneAnnotation] probe_type = ", probe_type)
    genes <- ngs.getGeneAnnotation_ORGDB(
      probes = probes,
      organism = organism,
      probe_type = probe_type
    )
  }

  ## try biomaRt for the rest
  if (is.null(genes) && organism != "No organism") {
    message("[ngs.getGeneAnnotation] >>> annotating genes using biomaRt")
    mart <- use_mart(organism)
    if (is.null(mart)) {
      message("[ngs.getGeneAnnotation] FAIL : could not connect to mart")
    }
    probe_type <- detect_probe.DEPRECATED(probes, mart)
    message("[ngs.getGeneAnnotation] probe_type = ", probe_type)
    if (is.null(probe_type)) stop("probe_type is NULL")
    genes <- ngs.getGeneAnnotation_BIOMART(
      probes = probes,
      organism = organism,
      probe_type = as.character(probe_type),
      mart = mart
    )
  } else if (organism == "No organism") {
    genes <- pgx.custom_annotation(pgx, custom_annot = annot_table)
  }

  if (is.null(genes)) {
    warning("[getGeneAnnotation] ERROR : could not create gene annotation")
    return(NULL)
  }
  return(genes)
}

ngs.getGeneAnnotation_BIOMART <- function(
    probes,
    organism,
    mart = NULL,
    probe_type = NULL,
    verbose = TRUE) {
  # Check mart
  if (is.null(mart)) {
    mart <- use_mart(organism)
    if (is.null(mart)) stop("ERROR : could not establish connection to MART.")
  }

  # Prepare inputs
  if (verbose) {
    message("[ngs.getGeneAnnotation_BIOMART] Retrieving gene annotations...")
  }

  clean_probes <- probes[!duplicated(probes)]

  if (is.null(probe_type)) {
    probe_type <- detect_probe.DEPRECATED(probes, mart)
  }

  # Select attributes
  attr_call <- c(
    probe_type,
    "external_gene_name", # gene_name
    "description", # gene_title
    "gene_biotype", # gene_biotype
    "chromosome_name", # chr
    "transcript_start", # pos
    "transcript_length", # tx_len
    "band" # map
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
      attributes = c("external_gene_name", "hsapiens_homolog_associated_gene_name"),
      filters = "external_gene_name",
      values = annot[!is.na(external_gene_name), external_gene_name],
      mart = mart
    )
    annot_homologs <- data.table::data.table(annot_homologs)
    data.table::setnames(annot_homologs,
      old = "hsapiens_homolog_associated_gene_name",
      new = "human_ortholog"
    )
    annot <- annot[annot_homologs, on = "external_gene_name", mult = "first"]
  }

  # Join with clean_probes vector
  out <- annot[clean_probes, on = probe_type, mult = "first"]

  # Renaming for backwards compatibility
  if (probe_type != "external_gene_name") {
    new_names <- c(
      "feature",
      "symbol",
      "gene_title",
      "gene_biotype",
      "chr",
      "pos",
      "tx_len",
      "map"
    )
  } else {
    new_names <- c(
      "feature",
      "gene_title",
      "gene_biotype",
      "chr",
      "pos",
      "tx_len",
      "map"
    )
    out[, symbol := external_gene_name]
  }
  data.table::setnames(out, old = attr_call, new = new_names)

  # Reorder columns and rows (remove human_ortholog if not present)
  col_order <- c("feature", "symbol", "human_ortholog", "gene_title", "gene_biotype")
  col_order <- col_order[col_order %in% colnames(out)]
  data.table::setcolorder(out, col_order)
  data.table::setkeyv(out, "feature")

  # Take out the source info from gene_title
  out[, c("gene_title", "source") :=
    data.table::tstrsplit(gene_title, "\\[", keep = 1:2)]
  out[, source := gsub("\\]", "", source)]

  if (organism == "Saccharomyces cerevisiae") {
    out[, gene_title := data.table::tstrsplit(gene_title, ";", keep = 1)]
  }
  out[, gene_title := trimws(gene_title, which = "right")]

  # Switch NA to empty string
  character_cols <- names(out)[sapply(out, is.character)]
  out[, (character_cols) :=
    lapply(.SD, function(x) data.table::fifelse(is.na(x), "", x)),
  .SDcols = character_cols
  ]

  # Keep it for back compatibility
  out[, gene_name := feature]

  # Return as data.frame and in the same order as input probes
  out <- as.data.frame(out)
  rownames(out) <- out$feature
  out <- out[probes, , drop = FALSE]
  return(out)
}


#' @title Get gene annotation data
#' @description Retrieves gene annotation data from an organism-specific
#' annotation package.
#'
#' @param probes Character vector of gene/transcript identifiers.
#' @param probe_type Type of identifiers provided in probes (e.g. "ENSEMBL").
#' @param organism Organism name (e.g. "Human", "Mouse", "Rat").
#'
#' @return Data frame containing annotation data for the provided identifiers.
#'
#' @details Queries an annotation package for the specified organism
#' using the provided identifiers and probe type. Returns a cleaned and
#' reformatted data frame with feature identifiers as row names.
#'
#' @examples
#' \dontrun{
#' probes <- c("ENSMUSG00000051951", "ENSMUSG00000033845")
#' annot <- ngs.getGeneAnnotation_ORGDB(probes, "ENSEMBL", organism = "Mouse")
#' }
#' @export
ngs.getGeneAnnotation_ORGDB <- function(probes, organism, probe_type) {
  organism <- tolower(organism)
  if (is.null(probe_type)) {
    stop("must provide probe_type")
  }

  # Get org database and columns request
  if (organism == "human") {
    org_db <- org.Hs.eg.db::org.Hs.eg.db
    cols_req <- c("SYMBOL", "GENENAME", "CHR", "CHRLOC", "MAP", "GENETYPE")
  } else if (organism == "mouse") {
    org_db <- org.Mm.eg.db::org.Mm.eg.db
    cols_req <- c("SYMBOL", "GENENAME", "CHR", "CHRLOC", "GENETYPE")
  } else if (organism == "rat") {
    org_db <- org.Rn.eg.db::org.Rn.eg.db
    cols_req <- c("SYMBOL", "GENENAME", "CHR", "CHRLOC", "GENETYPE")
  }

  if (!probe_type %in% AnnotationDbi::keytypes(org_db)) {
    warning("[ngs.getGeneAnnotation_ORGDB] ERROR : probe_type not in keytypes: ", probe_type)
    warning(
      "[ngs.getGeneAnnotation_ORGDB] keytypes available: ",
      paste(AnnotationDbi::keytypes(org_db), collapse = " ")
    )
    return(NULL)
  }

  ## discard version numbers if ENSEMBL/ENSEMBLTRANS
  clean.probes <- probes
  if (grepl("^ENSEMBL", probe_type)) {
    clean.probes <- sub("[.][0-9]+", "", probes)
  }

  # Call for annotation table
  suppressWarnings(
    d <- AnnotationDbi::select(
      org_db,
      keys = clean.probes,
      keytype = probe_type,
      columns = cols_req
    )
  )

  d <- data.table::as.data.table(d)

  # Add human ortholog and map for non-human organisms
  # Here, we use the old strategy of capitalise the symbol
  if (organism %in% c("mouse", "rat")) {
    d$human_ortholog <- toupper(d$SYMBOL)
    d$MAP <- d$CHR
  }

  ## if ENSEMBL get original probe names with version
  if (probe_type == "ENSEMBL") {
    d$ENSEMBL <- probes[match(d$ENSEMBL, clean.probes)]
  }
  if (probe_type == "ENSEMBLTRANS") {
    d$ENSEMBLTRANS <- probes[match(d$ENSEMBLTRANS, clean.probes)]
  }

  # Rename cols, add extra cols, reorder cols and rows
  if (probe_type == "SYMBOL") {
    d[, feature := SYMBOL]
    old_names <- c("feature", "SYMBOL", "GENENAME", "CHR", "CHRLOC", "MAP", "GENETYPE")
  } else {
    old_names <- c(probe_type, "SYMBOL", "GENENAME", "CHR", "CHRLOC", "MAP", "GENETYPE")
  }
  data.table::setnames(d,
    old = old_names,
    new = c("feature", "symbol", "gene_title", "chr", "pos", "map", "gene_biotype")
  )

  d$tx_len <- 1000
  d$source <- "local db"
  d$CHRLOCCHR <- NULL
  d$gene_name <- d$feature

  col_order <- c("feature", "symbol", "human_ortholog", "gene_title", "gene_biotype")
  col_order <- col_order[col_order %in% colnames(d)]
  data.table::setcolorder(d, col_order)
  data.table::setkeyv(d, "feature")

  # Return data.frame with rownames
  d <- as.data.frame(d)

  # Current DF do not admit multiple rownames, since CHRLOC can have more than one
  # match per gene/probe, we remove duplicates
  d <- d[!duplicated(d$feature), ]
  rownames(d) <- d$feature
  ii <- match(probes, rownames(d))
  d <- d[ii, , drop = FALSE]
  return(d)
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
#' @examples
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
    query_col <- data.table::fifelse(query_col == "" | is.na(query_col),
      yes = probes,
      no = query_col
    )
  }

  # Return queryed col
  return(query_col)
}

use_mart <- function(organism) {
  organism <- capitalize(organism) ## in utils.R
  message("[use_mart] connecting to bioMART server for organism ", organism)
  species_info <- playbase::SPECIES_TABLE[species_name == organism]
  # Some species appear in more than one mart, select ensembl only to avoid confusion
  if (nrow(species_info) > 1) {
    species_info <- species_info[mart == "ensembl"]
    species_info <- species_info[1, ]
  }
  ensembl <- biomaRt::useEnsembl(biomart = "ensembl")
  mart <- biomaRt::useDataset(dataset = species_info$dataset, mart = ensembl)
  return(mart)
}

#' Guess probe type
#'
#' This function tries to automatically detect the probe type of a set
#' of input probes using regular exprssion or by testing mapping
#' success against different identifier types in a BioMart database.
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
#' # use genes database
#' ensembl <- useEnsembl(biomart = "genes", version = 110)
#' # here we see a list of 214 species available in ensembl
#' datasets <- listDatasets(ensembl)
#' # here we can select between 214 species
#' dataset_hsa <- searchDatasets(mart = ensembl, pattern = "hsapiens")
#' ensembl_species <- useDataset(dataset = dataset_hsa$dataset, mart = ensembl)
#' probes <- c(
#'   "ENSG00000230915.1", "ENSG00000275728.1", "ENSG00000277599.1",
#'   "ENSG00000186163.9", "ENSG00000164823.11", "ENSG00000274234.1",
#'   "ENSG00000282461.1", "ENSG00000283056.1", "ENSG00000239021.1",
#'   "ENSG00000214268.2", "ENSG00000206687.1", "ENSG00000171148.14",
#'   "ENSG00000250027.1", "ENSG00000244217.1", "ENSG00000103502.14",
#'   "ENSG00000213178.3", "ENSG00000235059.5", "ENSG00000204555.3",
#'   "ENSG00000221044.2", "ENSG00000267162.1"
#' )
#' type <- detect_probe(probes, mart = ensembl_species)
#' }
#'
#' @export
guess_probetype <- function(probes, organism = NULL, for.biomart = FALSE) {
  ## 1. determine probe type using regular expression
  probe_type <- xbioc::idtype(probes)
  if (probe_type == "") probe_type <- NULL

  ## 2. match with human/mouse/rat genesets
  if (is.null(probe_type)) {
    ## matches SYMBOL for human, mouse and rat
    symbol <- as.list(org.Hs.eg.db::org.Hs.egSYMBOL)
    avg.match <- mean(toupper(probes) %in% symbol)
    if (avg.match > 0.5) probe_type <- "SYMBOL"
  }

  ## 3. check if they are proteins
  if (is.null(probe_type)) {
    type.regex <- list(
      "UNIPROT" = "^[OPQ][0-9]",
      "REFSEQ"  = "^N[MP]_[0-9]+$"
    )
    avg.match <- sapply(type.regex, function(s) mean(grepl(s, probes)))
    if (any(avg.match > 0.5)) {
      probe_type <- names(which.max(avg.match))
    }
  }

  KEYTYPES <- c("ENSEMBL", "ENSEMBLPROT", "ENSEMBLTRANS", "ENTREZID", "REFSEQ", "SYMBOL", "UNIPROT")
  if (!probe_type %in% KEYTYPES) {
    warning("[guess_probetype] ERROR : unsupported probe_type: ", probe_type)
    warning("[guess_probetype] keytypes available: ", KEYTYPES)
    return(NULL)
  }

  ## for biomart we have different nomenclature
  if (for.biomart) {
    keytype2biomart <- c(
      "ENSEMBL" = "ensembl_gene_id",
      "ENSEMBLTRANS" = "ensembl_transcript_id",
      "ENSEMBLPROT" = "ensembl_peptide_id",
      "SYMBOL" = "external_gene_name",
      "UNIPROT" = "uniprot_gn_id",
      ##    "REFSEQPROT"   = "refseq_peptide",
      "REFSEQ" = "refseq_mrna"
    )
    probe_type <- keytype2biomart[probe_type]

    ## add version if versioned
    if (probe_type %in% c("ensembl_gene_id", "ensembl_transcript_id")) {
      has.version <- (mean(grepl("[.][0-9]+$", probes)) > 0.8)
      if (has.version) probe_type <- paste0(probe_type, "_version")
    }
  }

  if (is.null(probe_type) || probe_type == "") {
    warning("[guess_probe] Could not guess probe_type. Please provide.")
  } else {
    message("[guess_probe] auto-detected probe_type = ", probe_type)
  }
  probe_type
}

#' @export
guess_organism <- function(probes) {
  org.regex <- list(
    "Human" = "^ENSG|^ENST|^[A-D,F-Z]{3,}",
    "Mouse" = "^ENSMUS|^[A-Z][a-z]{2,}",
    "Rat" = "^ENSRNO|^[A-Z][a-z]{2,}"
  )
  org.match <- function(probes, org) {
    mean(grepl(org.regex[[org]], probes))
  }
  avg.match <- sapply(names(org.regex), function(g) org.match(probes, g))
  avg.match

  if (any(avg.match > 0.33)) {
    organism <- names(which.max(avg.match))
  } else {
    organism <- NULL
  }
  if (is.null(organism)) {
    warning("[guess_organism] Could not auto-detect organism.")
  } else {
    message("[guess_organism] auto-detected organism = ", organism)
  }
  organism
}

#' @export
id2symbol <- function(probes, organism = "human") {
  if (is.null(organism)) {
    organism <- guess_organism(probes)
    if (is.null(organism)) {
      stop("could not determine organism. please specify.")
    }
  }
  ## this auto-selects using ORG.DB or BIOMARRT
  genes <- ngs.getGeneAnnotation(probes, organism = organism, use_biomart = NULL)
  genes <- genes[match(probes, rownames(genes)), ] ## just to be sure
  ## just return the symbol
  genes$symbol
}


#' @title Custom Gene Annotation
#'
#' @description Adds custom gene annotation table to a pgx object
#'
#' @param pgx pgx object
#' @param custom_annot data.frame with custom annotation data. If provided, 
#' it has to contain at least the columns "feature", "symbol", "gene_name". Also,
#' the features has to match the rownames of the counts provided.
#' 
#'
#' @details This function allows adding a gene annotation data.frame to a pgx object when
#' the user has not provided an organism or it's not known.  The custom_annot data.frame 
#' should contain gene IDs that match the pgx object genes, plus any additional columns 
#' of annotation data.
#' 
#' The id_type parameter specifies the type of ID used in custom_annot to match genes.
#' Possible options are "symbol", "ensembl_gene_id", etc. By default it will try to match
#' on the "symbol" field.
#'
#' Any columns in custom_annot that match existing pgx gene annotation columns will 
#' overwrite the original data. New columns will be appended.
#'
#' @return The pgx object with custom gene annotation added/appended. The gene annotation
#' table has the same format as the one returned by pgx.gene_table(). However, the 
#' columns human_ortholog, gene_title, gene_biotype, chr, pos, tx_len, map, source are filled
#' with default values.
#' 
#' @examples
#' \dontrun{
#'  custom_annot <- data.frame(
#'    feature = c("A1", "A2", "A3"), 
#'    symbol = c("TP53", "MYC", "EGFR"),
#'    gene_name = c("A1", "A2", "A3")
#'  )
#'  
#'  pgx <- pgx.custom_annotation(pgx, custom_annot)
#' }
#' @export
pgx.custom_annotation <- function(pgx, custom_annot = NULL) {

  message("[pgx.custom_annotation] Adding custom annotation table...")
  # If the user has provided a custom gene table, check it and use it
  if (!is.null(custom_annot)) {
    required_cols <- c(
      "feature",
      "symbol",
      "gene_name"
    )

    if (!all(required_cols %in% colnames(custom_annot))) {
      missing_cols <- required_cols[!required_cols %in% colnames(custom_annot)]
      stop("Custom gene table must contain the following columns: ", 
          paste0(required_cols, collapse = ", "), "\ncols missing: ", paste0(missing_cols, collapse = ", "))
    }

    # add extra cols if not present
    message("[pgx.custom_annotation] Filling annotation table...")
    extra_cols <- c("human_ortholog", "gene_title", "gene_biotype",  
                    "chr", "pos", "tx_len", "map", "source"
                    )
    for (col_i in extra_cols) {
      if (!col_i %in% colnames(custom_annot)) {
        custom_annot[[col_i]] <- switch(col_i,
          "human_ortholog" = "",
          "gene_title" = "unknown",
          "gene_biotype" = "unknown", 
          "chr" = "unknown",
          "pos" = 0,
          "tx_len" = 0,
          "map" = "1",
          "source" = "custom"
          )
      }
    }

    # Conform annotation table to pgx$counts
    annot_genes <- sum(rownames(pgx$counts) %in% custom_annot$feature) 
    annot_fraction <- annot_genes/ nrow(pgx$counts)
    
    if (annot_fraction > .5) {
      # filter annotated table by pgx$counts rownames using match
      custom_annot <- custom_annot[match(rownames(pgx$counts), custom_annot$feature), ]
    } else {
      stop("[pgx.custom_annotation] Not enought annoated genes. Be sure 
        custom_annot$feature matches counts rownames")
    }

  } else {
    # Create custom gene table from counts rownames
    message("[pgx.custom_annotation] Creating annotation table from counts rownames...")
    custom_annot <- data.frame(
      feature = rownames(pgx$counts),
      symbol = rownames(pgx$counts),
      gene_name = rownames(pgx$counts),
      human_ortholog = "",
      gene_title = "unknown",
      gene_biotype = "unknown", 
      chr = "unknown",
      pos = 0,
      tx_len = 0,
      map = "1",
      source = "custom"
    )
  }

  rownames(custom_annot) <- rownames(pgx$counts)

  return(custom_annot)
}


#' @title Detect probe type from probe set
#' @description Detects the most likely probe type (ENSEMBL, SYMBOL, etc)
#' for a set of gene/transcript identifiers.
#'
#' @param probes Character vector of gene/transcript identifiers.
#' @param organism Organism name (e.g. "Human", "Mouse", "Rat").
#'
#' @return The probe type with the most matches to the input identifiers.
#'
#' @details Checks the input identifiers against known gene/transcript
#' identifiers from an organism-specific annotation package to determine the
#' most likely probe type. Useful for identifying the type of identifiers
#' in a gene expression matrix when not explicitly provided.
#'
#' @examples
#' \dontrun{
#' probes <- probes <- c("NM_001081979", "NM_001081980", "NM_001081981", "NM_001081982", "NM_001081983")
#' organism <- "Mouse"
#' type <- detect_probetype_ORGDB(probes, organism)
#' }
#' @export

detect_probetype_ORGDB <- function(probes, organism) {
  warning("DEPRECATED. Please use guess_probetype")

  # Get org database
  if (tolower(organism) == "human") {
    org_db <- org.Hs.eg.db::org.Hs.eg.db
  } else if (tolower(organism) == "mouse") {
    org_db <- org.Mm.eg.db::org.Mm.eg.db
  } else if (tolower(organism) == "rat") {
    org_db <- org.Rn.eg.db::org.Rn.eg.db
  }

  # Probe types
  keytypes <- c(
    "ENSEMBL", "ENSEMBLTRANS", "SYMBOL", "REFSEQ", "UNIPROT", "ACCNUM"
  )

  ##  key_matches <- vector("character", length(keytypes))
  key_matches <- rep(0L, length(keytypes))
  names(key_matches) <- keytypes

  ## discard version numbers if ENSEMBL
  if (mean(grepl("^ENS", probes)) > 0.8) {
    probes <- sub("[.][0-9]+", "", probes)
  }

  # Subset probes if too many
  if (length(probes) > 100) {
    probes <- probes[as.integer(seq(1, length(probes), 100))]
  }

  ## remove versioning postfix from ensembl
  if( mean(grepl("^ENST",probes)) > 0.5 ) {
    probes <- sub("[.][0-9]+$","",probes)  
  }
    
  # Iterate over probe types
  for (key in keytypes) {
    n <- 0
    probe_matches <- data.frame(NULL)
    try(
      probe_matches <- AnnotationDbi::select(org_db,
        keys = probes,
        keytype = key,
        columns = key
      ),
      silent = TRUE
    )
    n <- nrow(probe_matches)
    key_matches[key] <- n
  }

  ## Return top match
##  key_matches    
  if (all(key_matches == 0)) {
    stop("Probe type not found, please, check your probes")
  } else {
    top_match <- names(which.max(key_matches))
  }

  return(top_match)
}


## ================================================================================
## ========================= DEPRECATED ===========================================
## ================================================================================


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
#' # use genes database
#' ensembl <- useEnsembl(biomart = "genes", version = 110)
#' # here we see a list of 214 species available in ensembl
#' datasets <- listDatasets(ensembl)
#' # here we can select between 214 species
#' dataset_hsa <- searchDatasets(mart = ensembl, pattern = "hsapiens")
#' ensembl_species <- useDataset(dataset = dataset_hsa$dataset, mart = ensembl)
#' probes <- c(
#'   "ENSG00000230915.1", "ENSG00000275728.1", "ENSG00000277599.1",
#'   "ENSG00000186163.9", "ENSG00000164823.11", "ENSG00000274234.1",
#'   "ENSG00000282461.1", "ENSG00000283056.1", "ENSG00000239021.1",
#'   "ENSG00000214268.2", "ENSG00000206687.1", "ENSG00000171148.14",
#'   "ENSG00000250027.1", "ENSG00000244217.1", "ENSG00000103502.14",
#'   "ENSG00000213178.3", "ENSG00000235059.5", "ENSG00000204555.3",
#'   "ENSG00000221044.2", "ENSG00000267162.1"
#' )
#' type <- detect_probe(probes, mart = ensembl_species)
#' }
#'
#' @export
detect_probe.DEPRECATED <- function(probes, mart = NULL, verbose = TRUE) {
  warning("DEPRECATED. Please use guess_probetype")

  # Check mart
  if (is.null(mart)) {
    stop("[detect_probe] Mart not found. Specify a BioMart database to use.")
  }
  # Prepare inputs
  if (verbose) message("[detect_probe] guessing probe type...")

  clean_probes <- probes[!is.na(probes)]

  n <- length(clean_probes)
  # if number of probes above 10, keep only 10 random probes
  if (n > 100L) n2 <- 100L else n2 <- n
  subsample <- sample(1:n, n2)
  subset_probes <- clean_probes[subsample]

  # Vector with input types to check
  probe_types_to_check <- c(
    "ensembl_gene_id",
    "ensembl_transcript_id",
    "ensembl_transcript_id_version",
    "ensembl_peptide_id",
    "external_gene_name",
    "ensembl_gene_id_version",
    "uniprot_gn_id",
    "refseq_peptide",
    "refseq_mrna"
  )


  # often we see multiples probes at once
  # initially we can try to detect as many probes as possible from each probe type
  # the type that guess better, is the one that will be used
  # this approach still has issues, as sometimes we have mixed probe types in on study.

  ## !!!!!NEED REFACTORING!!!!! : this is so bad... (IK)
  probe_check <- rep(0L, length(probe_types_to_check))
  names(probe_check) <- probe_types_to_check
  for(i in seq_along(probe_types_to_check)) {
    probe_check[i] <- tryCatch(
      {
        tmp <- biomaRt::getBM(
          attributes = probe_types_to_check[i],
          filters = probe_types_to_check[i],
          values = subset_probes,
          mart = mart
        )
        Sys.sleep(5) # Sleep time to prevent bounce from ensembl for consecutive calls
        out <- nrow(tmp)
        out
      },
      error = function(e) {
        return(0)
      }
    )
    # If more than 50 probes are found, stop
    if(probe_check[i] > 50) {
      break
    }
  }

  # Check matches and return if winner
  if (all(probe_check == 0)) {
    # TODO the probe2symbol and detect_probe code should be used in data-preview,
    # and we should warn the user in case no matches are found
    stop("[detect_probe]Probe type not found, please, check your probes")
  } else {
    probe_type <- names(which.max(probe_check))
  }
  if (verbose) message("[detect_probe] detected probe type = ", probe_type)

  return(probe_type)
}


id2symbol.DEPRECATED <- function(probes, type = NULL, org = "human", keep.na = FALSE) {
  require(org.Hs.eg.db)
  require(org.Mm.eg.db)
  require(org.Rn.eg.db)

  ## strip postfix for ensemble codes
  if (mean(grepl("^ENS", probes)) > 0.5) {
    probes <- gsub("[.].*", "", probes)
  }

  if (is.null(type) || is.null(org)) {
    hs.list <- list(
      "human.ensembl" = unlist(as.list(org.Hs.egENSEMBL)),
      "human.ensemblTRANS" = unlist(as.list(org.Hs.egENSEMBLTRANS)),
      # "human.unigene" = unlist(as.list(org.Hs.egUNIGENE)),
      "human.refseq" = unlist(as.list(org.Hs.egREFSEQ)),
      "human.accnum" = unlist(as.list(org.Hs.egACCNUM)),
      "human.uniprot" = unlist(as.list(org.Hs.egUNIPROT)),
      "human.symbol" = unlist(as.list(org.Hs.egSYMBOL))
    )

    mm.list <- list(
      "mouse.ensembl" = unlist(as.list(org.Mm.egENSEMBL)),
      "mouse.ensemblTRANS" = unlist(as.list(org.Mm.egENSEMBLTRANS)),
      # "mouse.unigene" = unlist(as.list(org.Mm.egUNIGENE)),
      "mouse.refseq" = unlist(as.list(org.Mm.egREFSEQ)),
      "mouse.accnum" = unlist(as.list(org.Mm.egACCNUM)),
      "mouse.uniprot" = unlist(as.list(org.Mm.egUNIPROT)),
      "mouse.symbol" = unlist(as.list(org.Mm.egSYMBOL))
    )

    rn.list <- list(
      "rat.ensembl" = unlist(as.list(org.Rn.egENSEMBL)),
      "rat.ensemblTRANS" = unlist(as.list(org.Rn.egENSEMBLTRANS)),
      # "rat.unigene" = unlist(as.list(org.Rn.egUNIGENE)),
      "rat.refseq" = unlist(as.list(org.Rn.egREFSEQ)),
      "rat.accnum" = unlist(as.list(org.Rn.egACCNUM)),
      "rat.uniprot" = unlist(as.list(org.Rn.egUNIPROT)),
      "rat.symbol" = unlist(as.list(org.Rn.egSYMBOL))
    )

    id.list <- c(hs.list, mm.list, rn.list)
    mx <- sapply(id.list, function(id) mean(probes %in% id))
    mx
    org <- type <- NULL
    max.mx <- max(mx, na.rm = TRUE)
    mx0 <- names(mx)[which.max(mx)]
    org <- sub("[.].*", "", mx0)
    type <- sub(".*[.]", "", mx0)
    message("[id2symbol] mapped ", format(100 * max.mx, digits = 2), "% of probes")
    if (max.mx < 0.5 && max.mx > 0) {
      message("[id2symbol] WARNING! low mapping ratio: r= ", max.mx)
    }
    if (max.mx == 0) {
      message("[id2symbol] WARNING! zero mapping ratio: r= ")
      type <- NULL
    }
    org
    type
  }

  ## checks
  if (is.null(type)) {
    cat("id2symbol: missing type: type = NULL\n")
    return(NULL)
  }
  if (!type %in% c("ensembl", "ensemblTRANS", "unigene", "refseq", "accnum", "uniprot", "symbol")) {
    cat("id2symbol: invalid type: ", type, "\n")
    return(NULL)
  }

  if (is.null(org)) {
    cat("id2symbol: missing organism: org = NULL\n")
    return(NULL)
  }
  if (!tolower(org) %in% c("human", "mouser", "rat")) {
    cat("id2symbol: not supported organism: org = ", org, "\n")
    return(NULL)
  }


  cat("[id2symbol] organism = ", org, "\n")
  cat("[id2symbol] probe.type = ", type, "\n")
  type

  if (type == "symbol") {
    cat("id2symbol: probe is already symbol\n")
    if (any(grep(" /// ", probes))) {
      symbol0 <- strsplit(probes, split = " /// ")
    } else if (any(grep("[;,]", probes))) {
      symbol0 <- strsplit(probes, split = "[;,\\|]")
    } else {
      symbol0 <- probes
    }
  } else {
    org
    if (org == "human") {
      symbol0 <- AnnotationDbi::mapIds(org.Hs.eg.db, probes, "SYMBOL", toupper(type))
    }
    if (org == "mouse") {
      symbol0 <- AnnotationDbi::mapIds(org.Mm.eg.db, probes, "SYMBOL", toupper(type))
    }
    if (org == "rat") {
      symbol0 <- AnnotationDbi::mapIds(org.Rn.eg.db, probes, "SYMBOL", toupper(type))
    }
  }

  ## Unrecognize probes
  nna <- which(is.na(names(symbol0)))
  length(nna)
  if (length(nna)) names(symbol0)[nna] <- probes[nna]

  ## What to do with unmapped/missing symbols????
  symbol <- sapply(symbol0, "[", 1) ## takes first symbol only!!!
  isnull <- which(sapply(symbol, is.null))
  symbol[isnull] <- NA
  if (keep.na) {
    sel.na <- which(is.na(symbol))
    symbol[sel.na] <- probes[sel.na]
  }
  symbol <- unlist(symbol)
  names(symbol) <- NULL
  Matrix::head(symbol)

  symbol
}


pgx.gene_table.OLD <- function(pgx, organism) {
  # Safety checks
  stopifnot(is.list(pgx))
  stopifnot(is.character(organism))

  # Init vals
  genes <- NULL
  ensembl_dataset <- NULL
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
  counter <- 0
  while (is.null(genes) && counter <= 5) {
    message(paste0("[pgx.gene_table] attempt to annotate genes, counter = ", counter + 1))
    # Set waiter so that we can make multiple calls with waiting time
    Sys.sleep(counter * 10)

    # Get ensembl
    if (is.null(ensembl_dataset)) {
      try({
        ensembl <- biomaRt::useEnsembl(biomart = "ensembl")
        ensembl_dataset <- biomaRt::useDataset(dataset = species_info$dataset, mart = ensembl)
      })
    }

    # Get probe type
    if (!is.null(ensembl_dataset) & is.na(probe_type)) {
      probe_type <- detect_probe.DEPRECATED(probes, ensembl_dataset)
    }

    # Get gene table
    if (!is.na(probe_type) & !is.null(ensembl_dataset)) {
      genes <- ngs.getGeneAnnotation_BIOMART(
        probes = probes,
        organism = organism,
        probe_type = probe_type,
        mart = ensembl_dataset
      )

      all_genes <- biomaRt::getBM(attributes = "external_gene_name", mart = ensembl_dataset)
      all_genes <- all_genes[, 1]
    }
    counter <- counter + 1
  }

  if (counter > 5 && is.null(genes)) {
    message("[pgx.gene_table] WARNING. could not reach ensembl server to get gene annotation")
    return(pgx)
  }

  # Return data
  pgx$genes <- genes
  pgx$all_genes <- all_genes ## Why need this?? unneeded baggage.. (IK)
  pgx$probe_type <- probe_type ## Why need this?? unneeded baggage.. (IK)

  return(pgx)
}
