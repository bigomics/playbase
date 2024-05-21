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
#' pgx <- pgx.addGeneAnnotation(pgx, "Human")
#' }
#' @export
pgx.addGeneAnnotation <- function(pgx, organism = NULL, annot_table = NULL, use_biomart = NULL) {
  # Safety checks
  stopifnot(is.list(pgx))
  probes <- rownames(pgx$counts)

  if (is.null(organism) && !is.null(pgx$organism)) {
    organism <- pgx$organism
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
ngs.getGeneAnnotation <- function(probes, organism, pgx = NULL,
                                  annot_table = NULL, use_biomart = NULL) {
  if (is.null(organism)) {
    warning("[getGeneAnnotation] Please specify organism")
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
      organism = organism,
      probes = probes,
      probe_type = probe_type
    )
  }

  ## try biomaRt for the rest
  if (is.null(genes) && organism != "No organism") {
    genes <- ngs.getGeneAnnotation_ANNOTHUB(
      organism = organism,
      probes = probes,
      probe_type = NULL,
      verbose = FALSE
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
ngs.getGeneAnnotation_ORGDB <- function(organism, probes, probe_type) {
  organism <- tolower(organism)
  if (is.null(probe_type)) {
    probe_type <- guess_probetype(probes, organism, for.biomart = FALSE)
  }
  message("probe_type = ", probe_type)

  if (organism == "homo sapiens") organism <- "human"
  if (organism == "mus musculus") organism <- "mouse"
  if (organism == "rattus norvegicus") organism <- "rat"

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
  } else {
    stop("ERROR: organism", organism, "not supported by ORGDB")
    return(NULL)
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
  d$source <- org_db$packageName
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

  if (is.null(d$human_ortholog)) {
    d$human_ortholog <- NA
  }

  return(d)
}

#' Get gene annotation data using AnnotationHub/OrgDb
#'
#' Retrieves gene annotation information from AnnotationHub for a set of input
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
#'   \item \code{feature}: The probe identifier.
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
ngs.getGeneAnnotation_ANNOTHUB <- function(
    organism,
    probes,
    probe_type = NULL,
    verbose = TRUE) {
  # Prepare inputs
  if (verbose) {
    message("[ngs.getGeneAnnotation_ANNOTHUB] Retrieving gene annotation...")
  }
  ##  require(AnnotationHub)
  ##  require(GO.db)

  if (tolower(organism) == "human") organism <- "Homo sapiens"
  if (tolower(organism) == "mouse") organism <- "Mus musculus"
  if (tolower(organism) == "rat") organism <- "Rattus norvegicus"
  organism

  ## Load the annotation resource.
  suppressMessages({
    ah <- AnnotationHub::AnnotationHub()
    cat("querying AnnotationHub for", organism, "\n")
    ahDb <- AnnotationHub::query(ah, pattern = c(organism, "OrgDb"))

    ## select on exact organism name
    ahDb <- ahDb[which(tolower(ahDb$species) == tolower(organism))]
    k <- length(ahDb)
    cat("selecting database for", ahDb$species[k], "\n")
    orgdb <- tryCatch(
      {
        ahDb[[k]]
      },
      error = function(e) {
        message("An error occurred: ", e, ". Retrying with force=TRUE.")
        ahDb[[k, force = TRUE]]
      }
    )
  })

  if (is.null(probes)) {
    probes <- keys(orgdb)
  }
  probes0 <- probes
  names(probes) <- probes0

  ## clean up probes
  probes <- probes[!is.na(probes) & probes != ""]
  is.ensembl <- mean(grepl("^ENS", probes)) > 0.5
  if (is.ensembl) probes <- sub("[.][0-9]+$", "", probes)
  if (sum(duplicated(probes)) > 0) {
    message("WARNING: duplicated probes")
    probes <- probes[!duplicated(probes)]
  }

  if (is.null(probe_type)) {
    probe_type <- detect_probetype.ANNOTHUB(organism, probes, ah = ah)
  }
  message("detected probe_type = ", probe_type)

  ## --------------------------------------------
  ## retrieve table
  ## --------------------------------------------
  keytypes(orgdb)
  cols <- c(
    "SYMBOL", "GENENAME", "GENETYPE", "ENTREZID",
    ## "ALIAS", "ACCNUM","REFSEQ",  ## balloon warning!!!
    "ENSEMBL", "ENSEMBLPROT", "UNIPROT", "GENETYPE", "MAP", "MGI"
  )
  mean_transcript <- mean(grepl("ENS[A-Z]*T[0]", probes))
  mean_transcript
  is_transcript <- (mean_transcript > 0.33)
  if (is_transcript) {
    ## add transcript/peptide level
    cols <- c(cols, "ENSEMBLTRANS")
  }
  cols <- intersect(cols, keytypes(orgdb))

  cat("get gene annotation columns:", cols, "\n")
  message("retrieving annotation for ", length(probes), " features...")
  annot <- AnnotationDbi::select(orgdb, keys = probes, columns = cols, keytype = probe_type)

  annot$SYMBOL[is.na(annot$SYMBOL)] <- ""

  ## match annotation table to probes
  cat("got", length(unique(annot$SYMBOL)), "unique SYMBOLs...\n")
  annot <- annot[match(probes, annot[, probe_type]), ]
  annot$PROBE <- names(probes) ## original probes

  ## --------------------------------------------
  ## get human ortholog using 'orthogene'
  ## --------------------------------------------
  ortho.map <- orthogene::map_species(method = "gprofiler")
  head(ortho.map)
  cat("\ngetting human orthologs...\n")
  if (organism == "Homo sapiens") {
    annot$ORTHOGENE <- annot$SYMBOL
  } else if (organism %in% ortho.map$scientific_name) {
    ortho.out <- orthogene::convert_orthologs(
      gene_df = unique(annot$SYMBOL),
      input_species = organism,
      output_species = "human",
      non121_strategy = "drop_both_species",
      method = "gprofiler"
    )
    ii <- match(annot$SYMBOL, ortho.out$input_gene)
    annot$ORTHOGENE <- rownames(ortho.out)[ii]
  } else {
    message("WARNING: ", organism, " not found in orthogene database. please check name.")
    annot$ORTHOGENE <- NA
  }

  ## Return as standardized data.frame and in the same order as input
  ## probes.
  annot$SOURCE <- ahDb$dataprovider
  annot.cols <- c(
    "PROBE", "SYMBOL", "ORTHOGENE", "GENENAME", "GENETYPE",
    "MAP", "CHR", "POS", "TXLEN", "SOURCE", "SYMBOL"
  )
  missing.cols <- setdiff(annot.cols, colnames(annot))
  missing.cols
  out <- annot
  for (a in missing.cols) out[[a]] <- NA
  out <- out[, annot.cols]
  new.names <- c(
    "feature", "symbol", "human_ortholog", "gene_title", "gene_biotype",
    "map", "chr", "pos", "tx_len", "source", "gene_name"
  )
  colnames(out) <- new.names
  out <- as.data.frame(out)
  if (!all(probes0 %in% out$feature)) {
    message("WARNING: not all probes could be annotated")
  }
  out <- out[match(probes0, out$feature), , drop = FALSE]

  # add space after ; to conform with playbase <= 1.3.2
  out$gene_title <- gsub(";", "; ", out$gene_title)

  # rename protein-coding to protein_coding to confirm with playbase <= v1.3.2
  out$gene_biotype <- sub("protein-coding", "protein_coding", out$gene_biotype)

  # replace NA in gene_ortholog by "" to conform with old pgx objects
  out$human_ortholog[is.na(out$human_ortholog)] <- ""

  # if organism is human, human_ortholog should be NA (matching old playbase annot)
  if (is.null(out$human_ortholog)) {
    out$human_ortholog <- NA
  }


  rownames(out) <- out$feature
  return(out)
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

## ================================================================================
## ========================= FUNCTIONS ============================================
## ================================================================================


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
  ## organism <- capitalize(organism) ## in utils.R
  organism <- sub("[H|h]omo sapiens", "Human", organism)
  organism <- sub("[M|m]us musculus", "Mouse", organism)
  organism <- sub("[R|r]attus norvegicus", "Rat", organism)
  message("[use_mart] connecting to bioMART server for organism ", organism)
  species_info <- playbase::SPECIES_TABLE[tolower(species_name) == tolower(organism)]
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
guess_probetype <- function(probes, organism, for.biomart = FALSE) {
  detect_probetype.MATCH(probes = probes, organism = organism, for.biomart = for.biomart)
}

detect_probetype.MATCH <- function(probes, organism = "", for.biomart = FALSE) {
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

  ## 1. match with human/mouse/rat genesets
  if (probe_type == "") {
    ## matches SYMBOL for human, mouse and rat
    symbol <- as.list(org.Hs.eg.db::org.Hs.egSYMBOL)
    avg.match <- mean(toupper(probes) %in% symbol)
    avg.match
    if (avg.match > 0.3) probe_type <- "SYMBOL"
  }

  ## 2. determine probe type using regular expression
  if (probe_type == "") {
    ## probe_type <- xbioc::idtype(probes)
    idtype.table <- table(sapply(head(sample(probes), 1000), xbioc::idtype))
    probe_type <- names(which.max(idtype.table))
  }

  ## 3. check if they are proteins
  if (probe_type == "") {
    type.regex <- list(
      "UNIPROT" = "^[OPQ][0-9]",
      "REFSEQ"  = "^N[MP]_[0-9]+$"
    )
    probe_type <- best.match(type.regex, 0.33)
  }

  ## 4. if probetype == genebank, replace by uniprot
  if (probe_type == "GENBANK") {
    probe_type <- "UNIPROT"
  }

  KEYTYPES <- c("ENSEMBL", "ENSEMBLPROT", "ENSEMBLTRANS", "ENTREZID", "REFSEQ", "SYMBOL", "UNIPROT")
  if (!probe_type %in% KEYTYPES) {
    warning("[detect_probetype] ERROR : unsupported probe_type: ", probe_type)
    warning("[detect_probetype] keytypes available: ", paste(KEYTYPES, collapse = " "))
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
    warning("[detect_probe] Could not guess probe_type. Please provide.")
  } else {
    message("[detect_probe] auto-detected probe_type = ", probe_type)
  }
  probe_type
}

#' @title Detect probe type from probe set
#' @export
detect_probetype.ANNOTHUB <- function(organism, probes, ah = NULL, nprobe = 100) {
  ##  require(AnnotationHub)
  ##  require(GO.db)

  if (tolower(organism) == "human") organism <- "Homo sapiens"
  if (tolower(organism) == "mouse") organism <- "Mus musculus"
  if (tolower(organism) == "rat") organism <- "Rattus norvegicus"
  organism

  ## Load the annotation resource.
  if (is.null(ah)) {
    ah <- AnnotationHub::AnnotationHub()
  }
  suppressMessages({
    cat("querying AnnotationHub for", organism, "\n")

    ahDb <- AnnotationHub::query(ah, pattern = c(organism, "OrgDb"))

    ## select on exact organism name
    ahDb <- ahDb[which(tolower(ahDb$species) == tolower(organism))]
    k <- length(ahDb)
    cat("selecting database for", ahDb$species[k], "\n")
    orgdb <- ahDb[[k]] ## last one, newest version
  })

  ## clean up probes
  probes <- probes[!is.na(probes) & probes != ""]
  if (sum(duplicated(probes)) > 0) {
    message("WARNING: duplicated probes")
    probes <- unique(probes)
  }

  ## get probe types for organism
  keytypes <- c(
    "ENSEMBL", "ENSEMBLTRANS", "SYMBOL", "REFSEQ", "UNIPROT",
    "ACCNUM", "ENTREZID"
  )
  keytypes <- c(
    "SYMBOL", "ENTREZID", "ACCNUM", "REFSEQ",
    "ENSEMBL", "ENSEMBLTRANS", "MGI",
    "ENSEMBLPROT", "UNIPROT"
  )
  keytypes <- intersect(keytypes, keytypes(orgdb))

  key_matches <- rep(0L, length(keytypes))
  names(key_matches) <- keytypes

  ## clean up probes
  probes <- probes[!is.na(probes) & probes != ""]
  if (sum(duplicated(probes)) > 0) {
    message("WARNING: duplicated probes")
    probes <- unique(probes)
  }

  ## discard version numbers if ENSEMBL
  if (mean(grepl("^ENS", probes)) > 0.5) {
    probes <- sub("[.][0-9]+$", "", probes)
  }

  ## Subset probes if too many
  if (length(probes) > nprobe) {
    if (nprobe > length(probes)) nprobe <- length(probes)

    # get random probes for query
    probes <- sample(probes, nprobe)
  }

  # Iterate over probe types
  for (key in keytypes) {
    # key = keytypes[4]
    n <- 0
    probe_matches <- data.frame(NULL)
    try(
      probe_matches <- AnnotationDbi::select(
        orgdb,
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
  top_match <- NULL
  if (all(key_matches == 0)) {
    message("WARNING:: Probe type not found, please check your probes")
    return(NULL)
  } else {
    top_match <- names(which.max(key_matches))
    message("Guessed probe type = ", top_match)
  }

  return(top_match)
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
  if (grepl("homo sapiens", tolower(organism))) organism <- "Human"
  if (grepl("mus musculus", tolower(organism))) organism <- "Mouse"
  if (grepl("rattus norvegicus", tolower(organism))) organism <- "Rat"

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
    message("WARNING: unsupported organism")
    return(NULL)
  }

  # Probe types
  keytypes <- c(
    "ENSEMBL", "ENSEMBLTRANS", "SYMBOL", "REFSEQ", "UNIPROT", "ACCNUM", "ENTREZID"
  )

  ##  key_matches <- vector("character", length(keytypes))
  key_matches <- rep(0L, length(keytypes))
  names(key_matches) <- keytypes

  ## discard version numbers if ENSEMBL
  if (mean(grepl("^ENS", probes)) > 0.5) {
    probes <- sub("[.][0-9]+", "", probes)
  }

  # Subset probes if too many
  if (length(probes) > 100) {
    probes <- probes[as.integer(seq(1, length(probes), 100))]
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
  key_matches
  top_match <- NULL
  if (all(key_matches == 0)) {
    message("WARNING:: Probe type not found, please check your probes")
  } else {
    top_match <- names(which.max(key_matches))
    message("Guessed probe type = ", top_match)
  }

  return(top_match)
}


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
  ##  warning("DEPRECATED. Please use guess_probetype")

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
    "entrezgene_id",
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


#' @export
id2symbol <- function(probes, organism) {
  if (is.null(organism)) {
    stop("could not determine organism. please specify.")
  }
  ## this auto-selects using ORG.DB or BIOMART
  genes <- ngs.getGeneAnnotation(probes, organism = organism, use_biomart = NULL)
  genes <- genes[match(probes, rownames(genes)), ] ## just to be sure
  ## just return the symbol
  genes$symbol
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
check_known_probes.NOTUSED <- function(probes, probe_types_to_check = NULL) {
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
