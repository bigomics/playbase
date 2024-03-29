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
#' pgx <- pgx.pgx.addGeneAnnotation(pgx, "Human")
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
ngs.getGeneAnnotation <- function(probes, pgx = NULL, organism = NULL, annot_table = NULL, use_biomart = NULL) {
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
    probe_type <- detect_probetype.BIOMART(probes, mart)
    message("[ngs.getGeneAnnotation] probe_type = ", probe_type)
    if (is.null(probe_type)) stop("probe_type is NULL")
    genes <- ngs.getGeneAnnotation_BIOMART(
      probes = probes,
      organism = organism,
      probe_type = as.character(probe_type),
      mart = mart
    )
  } else if (organism == "No organism" && !is.null(pgx)) {
    genes <- pgx.custom_annotation(counts = pgx$counts, custom_annot = annot_table)
  }

  if (is.null(genes)) {
    warning("[getGeneAnnotation] ERROR : could not create gene annotation")
    return(NULL)
  }
  return(genes)
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
guess_probetype <- function(probes, organism = "", for.biomart = FALSE) {
  best.match <- function(type.regex, avg.min = 0.33) {
    avg.match <- sapply(type.regex, function(s) mean(grepl(s, probes)))
    probe_type <- ""
    if (any(avg.match > 0.5)) {
      probe_type <- names(which.max(avg.match))
    }
    probe_type
  }
  probe_type <- ""

  ## 0. Matching: worm, fly, yeast
  type.regex <- list(
    "ENSEMBL" = "^WBGene|^FBgn0|^Y[A-Z][LR]" ## worm, fly, yeast
  )
  probe_type <- best.match(type.regex, 0.33)


  ## 1. determine probe type using regular expression
  if (probe_type == "") {
    probe_type <- xbioc::idtype(probes)
  }

  ## 2. match with human/mouse/rat genesets
  if (probe_type == "") {
    ## matches SYMBOL for human, mouse and rat
    symbol <- as.list(org.Hs.eg.db::org.Hs.egSYMBOL)
    avg.match <- mean(toupper(probes) %in% symbol)
    avg.match
    if (avg.match > 0.3) probe_type <- "SYMBOL"
  }

  ## 3. check if they are proteins
  if (probe_type == "") {
    type.regex <- list(
      "UNIPROT" = "^[OPQ][0-9]",
      "REFSEQ"  = "^N[MP]_[0-9]+$"
    )
    probe_type <- best.match(type.regex, 0.33)
  }

  KEYTYPES <- c("ENSEMBL", "ENSEMBLPROT", "ENSEMBLTRANS", "ENTREZID", "REFSEQ", "SYMBOL", "UNIPROT")
  if (!probe_type %in% KEYTYPES) {
    warning("[guess_probetype] ERROR : unsupported probe_type: ", probe_type)
    warning("[guess_probetype] keytypes available: ", paste(KEYTYPES, collapse = " "))
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
      ## "REFSEQPROT"   = "refseq_peptide",
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
    "Human" = "^ENSG|^ENST|^[A-Z]+[0-9]*$",
    "Mouse" = "^ENSMUS|^[A-Z][a-z]{2,}",
    "Rat" = "^ENSRNO|^[A-Z][a-z]{2,}",
    "Caenorhabditis elegans" = "^WBGene",
    "Drosophila melanogaster" = "^FBgn0",
    "Saccharomyces cerevisiae" = "^Y[A-P][RL]"
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
#' @param counts A counts matrix
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
#' custom_annot <- data.frame(
#'   feature = c("A1", "A2", "A3"),
#'   symbol = c("TP53", "MYC", "EGFR"),
#'   gene_name = c("A1", "A2", "A3")
#' )
#'
#' pgx <- pgx.custom_annotation(counts, custom_annot)
#' }
#' @export
pgx.custom_annotation <- function(counts, custom_annot = NULL) {
  message("[pgx.custom_annotation] Adding custom annotation table...")
  # If the user has provided a custom gene table, check it and use it

  annot_genes <- sum(rownames(counts) %in% custom_annot$feature)

  annot_map <- list(
    "human_ortholog" = "",
    "gene_title" = "unknown",
    "gene_biotype" = "unknown",
    "chr" = "unknown",
    "pos" = 0,
    "tx_len" = 0,
    "map" = "1",
    "source" = "custom"
  )

  required_cols <- c(
    "feature",
    "symbol",
    "gene_name"
  )

  # this will be used at the end to order df columns
  table_col_order <- c(required_cols, names(annot_map))

  # legacy code but maybe this could be removed in the future...
  required_in_annot <- all(required_cols %in% colnames(custom_annot))

  if (!is.null(custom_annot) && annot_genes > 1 && required_in_annot) {
    # remove all NA columns, otherwise the for loop below will not work
    custom_annot <- custom_annot[, !apply(custom_annot, 2, function(x) all(is.na(x)))]

    # identify missing columns and fill them with annot_map
    missing_cols <- setdiff(names(annot_map), names(custom_annot))

    custom_annot[missing_cols] <- annot_map[missing_cols]

    # filter annotated table by counts using match
    custom_annot <- custom_annot[match(rownames(counts), custom_annot$feature), ]

    # if row was missing from annotation table (NA from match call above), input NA based on rownames(counts)

    rownames(custom_annot) <- rownames(counts)

    custom_annot$feature <- ifelse(is.na(custom_annot$feature), rownames(custom_annot), custom_annot$feature)
    custom_annot$symbol <- ifelse(is.na(custom_annot$symbol), rownames(custom_annot), custom_annot$symbol)
    custom_annot$gene_name <- ifelse(is.na(custom_annot$gene_name), rownames(custom_annot), custom_annot$gene_name)

    # Fill NA values with corresponding values from annot_map

    res <- lapply(names(annot_map), function(x) {
      ifelse(is.na(custom_annot[[x]]), annot_map[[x]], custom_annot[[x]])
    })

    names(res) <- names(annot_map)

    res <- as.data.frame(res)

    custom_annot[, names(annot_map)] <- res[, names(annot_map)]
  } else {
    # Create custom gene table from counts rownames
    message("[pgx.custom_annotation] Creating annotation table from counts rownames...")
    custom_annot <- data.frame(
      feature = rownames(counts),
      symbol = rownames(counts),
      gene_name = rownames(counts),
      human_ortholog = "",
      gene_title = "unknown",
      gene_biotype = "unknown",
      chr = "unknown",
      pos = 0,
      tx_len = 0,
      map = "1",
      source = "custom"
    )
    rownames(custom_annot) <- rownames(counts)
  }

  custom_annot <- custom_annot[, table_col_order]

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
#' probes <- c("NM_001081979", "NM_001081980", "NM_001081981")
#' organism <- "Mouse"
#' type <- detect_probetype.ORGDB(probes, organism)
#' }
#' @export

detect_probetype.ORGDB <- function(probes, organism) {
  warning("DEPRECATED. Please use guess_probetype")

  ## Get org database
  org_db <- NULL
  if (tolower(organism) == "human") {
    org_db <- org.Hs.eg.db::org.Hs.eg.db
  } else if (tolower(organism) == "mouse") {
    org_db <- org.Mm.eg.db::org.Mm.eg.db
  } else if (tolower(organism) == "rat") {
    org_db <- org.Rn.eg.db::org.Rn.eg.db
  }

  if (is.null(org_db)) {
    return(NULL)
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
  if (mean(grepl("^ENST", probes)) > 0.5) {
    probes <- sub("[.][0-9]+$", "", probes)
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
detect_probetype.BIOMART <- function(probes, mart = NULL, verbose = TRUE) {
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
    "ensembl_gene_id_version",
    "ensembl_transcript_id",
    "ensembl_transcript_id_version",
    "ensembl_peptide_id",
    "external_gene_name",
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
  for (i in seq_along(probe_types_to_check)) {
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
    if (probe_check[i] > 50) {
      break
    }
  }

  # Check matches and return if winner
  if (all(probe_check == 0)) {
    # TODO the probe2symbol and detect_probe code should be used in data-preview,
    # and we should warn the user in case no matches are found
    stop("[detect_probe] Probe type not found, please, check your probes")
  } else {
    probe_type <- names(which.max(probe_check))
  }
  if (verbose) message("[detect_probe] detected probe type = ", probe_type)

  return(probe_type)
}


#' Known probes
#'
#' This function checks if a vector of probe identifiers matches known patterns
#' for common identifier types like Ensembl IDs. It returns a vector of probe
#' types that should be checked for mapping based on the matches.
#'
#' @param probes A character vector of probe identifiers.
#' @param probe_types_to_check A character vector of probe types to check.
#'
#' @return A character vector of probe types to check mapping for.
#'
#' @details
#' This function checks `probes` for the following known patterns:
#'
#' - Ensembl IDs for higher animals (ENS*)
#' - Drosophila melanogaster Ensembl IDs (FBgn*)
#' - Caenorhabditis elegans Ensembl IDs (WBGene*)
#'
#' If a high proportion (>50%) of `probes` match one of these patterns,
#' the corresponding Ensembl ID types are returned in `probe_types_to_check`.
#'
#' @examples
#' probes <- c("ENSG000001", "ENST0001", "FBgn001", "WBGene0001")
#' check_known_probes(probes, NULL)
#'
#' @export
check_known_probes <- function(probes, probe_types_to_check = NULL) {
  # Check higher animals ENSEMBL notation
  if (sum(grepl("^ENS*", probes)) > length(probes) * 0.5) {
    probe_types_to_check <- c(
      "ensembl_gene_id",
      "ensembl_transcript_id",
      "ensembl_transcript_id_version",
      "ensembl_peptide_id"
    )
  }

  # Check D. melanogaster notation
  if (sum(grepl("^FBgn*", probes)) > length(probes) * 0.5) {
    probe_types_to_check <- c(
      "ensembl_gene_id",
      "ensembl_transcript_id",
      "ensembl_transcript_id_version",
      "ensembl_peptide_id"
    )
  }

  # Check C. elegans notation
  if (sum(grepl("^WBGene*", probes)) > length(probes) * 0.5) {
    probe_types_to_check <- c(
      "ensembl_gene_id",
      "ensembl_transcript_id",
      "ensembl_transcript_id_version",
      "ensembl_peptide_id"
    )
  }

  return(probe_types_to_check)
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
ngs.getGeneAnnotation_BIOMART <- function(
    probes,
    organism,
    probe_type = NULL,
    mart = NULL,
    verbose = TRUE) {
  # Check mart
  if (is.null(mart)) {
    stop("[ngs.getGeneAnnotation_BIOMART] Mart not found. Please specify a BioMart database.")
  }

  # Prepare inputs
  if (verbose) {
    message("[ngs.getGeneAnnotation_BIOMART] Filling genes information...")
  }

  clean_probes <- probes[!duplicated(probes)]

  if (is.null(probe_type)) {
    probe_type <- detect_probe(probes, mart)
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

  # Reorder columns and rows
  if ("human_ortholog" %chin% colnames(out)) {
    col_order <- c(
      "feature",
      "symbol",
      "human_ortholog",
      "gene_title",
      "gene_biotype"
    )
  } else {
    col_order <- c(
      "feature",
      "symbol",
      "gene_title",
      "gene_biotype"
    )
  }
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
      probe_type <- detect_probetype.BIOMART(probes, ensembl_dataset)
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
