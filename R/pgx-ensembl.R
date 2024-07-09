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
pgx.addGeneAnnotation <- function(pgx, organism = NULL, annot_table = NULL) {
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
    organism = organism
  )

  # Return data
  pgx$genes <- genes

  return(pgx)
}


#' Get gene annotation data using AnnotationHub
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
ngs.getGeneAnnotation <- function(
    organism,
    probes,
    probe_type = NULL,
    pgx = NULL,
    annot_table = NULL,
    verbose = TRUE) {

  if (is.null(organism)) {
    warning("[getGeneAnnotation] Please specify organism")
    return(NULL)
  }

  if (verbose) {
    message("[ngs.getGeneAnnotation] Retrieving gene annotation...")
  }

  if (tolower(organism) == "human") organism <- "Homo sapiens"
  if (tolower(organism) == "mouse") organism <- "Mus musculus"
  if (tolower(organism) == "rat") organism <- "Rattus norvegicus"

  genes <- NULL

  ## get correct OrgDb database for this organism
  orgdb <- getOrgDb(organism, ah = NULL)

  if (is.null(probes)) {
    probes <- keys(orgdb)
  }
  probes0 <- probes
  names(probes) <- probes0

  ## clean up probes
  probes <- probes[!is.na(probes) & probes != ""]
  probes <- sapply(strsplit(probes, split = ";"), head, 1) ## take first
  is.ensembl <- mean(grepl("^ENS", probes)) > 0.5
  if (is.ensembl) probes <- sub("[.][0-9]+$", "", probes)

  if (is.null(probe_type)) {
    probe_type <- playbase::detect_probetype(organism, probes, ah = ah)
  }
  message("detected probe_type = ", probe_type)
  if (is.null(probe_type)) {
    message("ERROR: could not determine probe_type ")
    return(NULL)
  }

  ## --------------------------------------------
  ## retrieve table
  ## --------------------------------------------
  cols <- c("SYMBOL", "GENENAME", "GENETYPE", "MAP")
  cols <- intersect(cols, keytypes(orgdb))

  cat("get gene annotation columns:", cols, "\n")
  message("retrieving annotation for ", length(probes), " ", probe_type, " features...")

  suppressMessages( suppressWarnings( 
    annot <- AnnotationDbi::select(orgdb, keys = probes, columns = cols, keytype = probe_type)
  ))

  # some organisms do not provide symbol but rather gene name (e.g. yeast)
  if (!"SYMBOL" %in% colnames(annot)) {
    annot$SYMBOL <- annot$GENENAME
  }

  annot$SYMBOL[is.na(annot$SYMBOL)] <- ""

  ## match annotation table to probes
  cat("got", length(unique(annot$SYMBOL)), "unique SYMBOLs...\n")
  annot <- annot[match(probes, annot[, probe_type]), ]
  annot$PROBE <- names(probes) ## original probe names

  ## --------------------------------------------
  ## get human ortholog using 'orthogene'
  ## --------------------------------------------
  ortho.map <- orthogene::map_species(method = "gprofiler")
  head(ortho.map)
  cat("\ngetting human orthologs...\n")
  has.ortho <- organism %in% ortho.map$scientific_name
  has.symbol <- "SYMBOL" %in% colnames(annot)
  if (organism == "Homo sapiens") {
    annot$ORTHOGENE <- annot$SYMBOL
  } else if (has.ortho && has.symbol) {

    ortho.out <- orthogene::convert_orthologs(
      gene_df = unique(annot$SYMBOL),
      input_species = organism,
      output_species = "human",
      non121_strategy = "drop_both_species",
      method = "gprofiler"
    )

    if (dim(ortho.out)[1] == 0) {
      ortho.out <- orthogene::convert_orthologs(
        gene_df = unique(annot$SYMBOL),
        input_species = organism,
        output_species = "human",
        non121_strategy = "drop_both_species",
        method = "homologene"
      )
    }

    if (dim(ortho.out)[1] == 0) {
      ortho.out <- orthogene::convert_orthologs(
        gene_df = unique(annot$SYMBOL),
        input_species = organism,
        output_species = "human",
        non121_strategy = "drop_both_species",
        method = "babelgene"
      )
    }

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
    "MAP", "CHR", "POS", "TXLEN", "SOURCE"
  )
  missing.cols <- setdiff(annot.cols, colnames(annot))
  missing.cols
  genes <- annot
  for (a in missing.cols) genes[[a]] <- NA
  genes <- genes[, annot.cols]
  new.names <- c(
    "feature", "symbol", "human_ortholog", "gene_title", "gene_biotype",
    "map", "chr", "pos", "tx_len", "source"
  )
  colnames(genes) <- new.names
  genes <- as.data.frame(genes)

  # gene_name should ALWAYS be assigned to feature for compatibility
  # with gene_name legacy implementation
  genes$gene_name <- genes$feature

  if (!all(probes0 %in% genes$feature)) {
    message("WARNING: not all probes could be annotated")
  }
  genes <- genes[match(probes0, genes$feature), , drop = FALSE]

  # add space after ; to conform with playbase <= 1.3.2
  genes$gene_title <- gsub(";", "; ", genes$gene_title)

  # rename protein-coding to protein_coding to confirm with playbase <= v1.3.2
  genes$gene_biotype <- sub("protein-coding", "protein_coding", genes$gene_biotype)

  # replace NA in gene_ortholog by "" to conform with old pgx objects
  genes$human_ortholog[is.na(genes$human_ortholog)] <- ""

  # if organism is human, human_ortholog should be NA (matching old playbase annot)
  if (is.null(genes$human_ortholog)) {
    genes$human_ortholog <- NA
  }

  rownames(genes) <- probes0

  # annotation table is mandatory for No organism (until server side can handle missing genesets)
  if (organism == "No organism" && !is.null(pgx)) {
    genes <- pgx.custom_annotation(counts = pgx$counts, custom_annot = annot_table)
  }

  if (is.null(genes)) {
    warning("[getGeneAnnotation] ERROR : could not create gene annotation")
    return(NULL)
  }
  return(genes)
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


## not exported
getOrgDb <- function(organism, ah = NULL) {

  if (tolower(organism) == "human") organism <- "Homo sapiens"
  if (tolower(organism) == "mouse") organism <- "Mus musculus"
  if (tolower(organism) == "rat") organism <- "Rattus norvegicus"

  if(!is.null(ah)) {
    all_species <- getAllSpecies(ah)
  } else {
    ## If organism is in localHub we select localHub=TRUE because
    ## this is faster. Otherwise switch to online Hub
    suppressMessages(
      ah <- AnnotationHub::AnnotationHub(localHub=TRUE)
    )
    local_species <- getAllSpecies(ah)  ## orgDb species only
    if(tolower(organism) %in% tolower(local_species)) {
      message("[selectAnnotationHub] organism '",organism,"' in local Hub")
      all_species <- local_species
    } else {
      message("[selectAnnotationHub] querying online Hub...")
      ah <- AnnotationHub::AnnotationHub(localHub=FALSE)  
      all_species <- getAllSpecies(ah)
    }
  }

  if(!tolower(organism) %in% tolower(all_species)) {
    message("WARNING: organism '",organism,"' not in AnnotationHub")
    return(NULL)
  }

  ## correct capitalization
  species <- all_species[which(tolower(all_species) == tolower(organism))]
  
  suppressMessages({
    message("querying AnnotationHub for '", organism, "'\n")
    ahDb <- AnnotationHub::query(ah, pattern = c(organism, "OrgDb"))

    ## select on exact organism name
    ahDb <- ahDb[which(tolower(ahDb$species) == tolower(organism))]
    k <- length(ahDb)  ## latest of multiple
    message("selecting database for '", ahDb$species[k], "'\n")

    message("retrieving annotation...\n")
    orgdb <- tryCatch({
      ahDb[[k]]
    },
    error = function(e) {
      message("An error occurred: ", e, ". Retrying with force=TRUE.")
      ahDb[[k, force = TRUE]]
    })
  })
  
  return(orgdb)  
}


#' @title Detect probe type from probe set
#' @export
detect_probetype <- function(organism, probes, ah = NULL, nprobe = 100) {

  if (tolower(organism) == "human") organism <- "Homo sapiens"
  if (tolower(organism) == "mouse") organism <- "Mus musculus"
  if (tolower(organism) == "rat") organism <- "Rattus norvegicus"
  organism

  ## get correct OrgDb database for organism
  orgdb <- getOrgDb(organism, ah = ah)
  if(is.null(orgdb)) {
    message("WARNING: unsupported organism '",organism,"'\n")
    return(NULL)
  }
  
  ## get probe types for organism
  keytypes <- c(
    "SYMBOL", "ENSEMBL", "UNIPROT",
    "GENENAME", "MGI",    
    "ENSEMBLTRANS", "ENSEMBLPROT",
    "ACCNUM", "REFSEQ", "ENTREZID"
  )
  keytypes <- intersect(keytypes, keytypes(orgdb))
  key_matches <- rep(0L, length(keytypes))
  names(key_matches) <- keytypes

  ## clean up probes
  probes <- probes[!is.na(probes) & probes != ""]
  probes <- sapply(strsplit(probes, split = ";"), head, 1) ## take first
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
  key = keytypes[1]
  for (key in keytypes) {
    probe_matches <- data.frame(NULL)

    # add symbol and genename on top of key as they will be used to
    # count the real number of probe matches
    key2 <- c(key, c("SYMBOL", "GENENAME"))
    key2 <- intersect(key2, keytypes)
    suppressMessages( suppressWarnings( try(
      probe_matches <- AnnotationDbi::select(
        orgdb,
        keys = probes,
        keytype = key,
        columns = key2
      ),
      silent = TRUE
    )))

    # set empty character to NA, as we only count not-NA to define probe type
    probe_matches[probe_matches == ""] <- NA
    # check which probe types (genename, symbol) return the most matches
    n1 <- n2 <- 0
    if ("SYMBOL" %in% colnames(probe_matches)) n1 <- sum(!is.na(probe_matches[, "SYMBOL"]))
    if ("GENENAME" %in% colnames(probe_matches)) n2 <- sum(!is.na(probe_matches[, "GENENAME"]))
    matchratio <- max(n1, n2) / length(probes)
    key_matches[key] <- matchratio

    ## stop search prematurely if matchratio > 95%
    if(matchratio > 0.95) break()
  }

  ## Return top match
  ##  key_matches
  top_match <- NULL
  if (all(key_matches == 0)) {
    message("WARNING:: Probe type not found, use one of the following probe types: ", paste(keytypes, collapse = " "))
    return(NULL)
  } else {
    top_match <- names(which.max(key_matches))
    message("Guessed probe type = ", top_match)
  }

  return(top_match)
}

#' @title Get all species in AnnotationHub/OrgDB
#'
#' @export
getAllSpecies <- function(ah = NULL) {

  if (is.null(ah)) {
    ah <- AnnotationHub::AnnotationHub() ## make global??
  }
  db <- AnnotationHub::query(ah, "OrgDb")
  sort(unique(mcols(db)$species))
}


#' @title Get species table in AnnotationHub/OrgDB
#'
#' @export
getSpeciesTable <- function(ah = NULL) {

  if (is.null(ah)) {
    ah <- AnnotationHub::AnnotationHub() ## make global??
  }
  ah.tables <- AnnotationHub::query(ah, "OrgDb")
  
  variables <- c(
    "ah_id", "species", "description", "rdatadateadded", "rdataclass",
    "title", "taxonomyid", "coordinate_1_based", "preparerclass", "sourceurl",
    "dataprovider", "genome", "maintainer", "tags", "sourcetype"
  )

  # Iterate through each variable and store it as a table
  tables <- lapply(variables, function(var) {
    table <- eval(parse(text = paste0("ah.tables$", var)))
  })
  tables <- do.call(cbind, tables)
    
  colnames(tables) <- variables
  names(tables) <- variables
  return(tables)
}

#' @title Check if probes are valid for organism
#'
#' @return TRUE    if probes match any probetype of organism
#' @return FALSE   if probes do not any probetype of organism
#'
#' @export
checkProbes <- function(organism, probes, ah = NULL) {
  probe_type <- detect_probetype(organism, probes, ah = ah)
  if (is.null(probe_type)) {
    message("[checkProbes] WARNING: could not validate probes")
    return(FALSE)
  }
  message("[checkProbes] detected probe_type = ", probe_type)
  return(TRUE)
}
