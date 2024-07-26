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

  if (organism == "No organism") {
    # annotation table is mandatory for 'No organism' (until server side
    # can handle missing genesets)
    genes <- getCustomAnnotation(
      probes = probes,
      custom_annot = annot_table
    )
  }

  if (organism != "No organism") {
    # Get gene table
    genes <- ngs.getGeneAnnotation(
      organism = organism,
      probes = probes
    )
  }

  # Add to pgx object
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
    annot_table = NULL,
    use.ah = NULL,
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
  orgdb <- getOrgDb(organism, use.ah = use.ah)

  if (is.null(probes)) {
    probes <- AnnotationDbi::keys(orgdb)
  }
  probes0 <- probes
  names(probes) <- probes0

  ## clean up probes
  probes <- probes[!is.na(probes) & probes != ""]
  probes <- sapply(strsplit(probes, split = ";"), head, 1) ## take first
  is.ensembl <- mean(grepl("^ENS", probes)) > 0.5
  if (is.ensembl) probes <- sub("[.][0-9]+$", "", probes)

  is.uniprot <- mean(grepl("^[QP][0-9]*", probes)) > 0.8
  if (is.uniprot) probes <- sub("-[0-9]+", "", probes)

  if (is.null(probe_type)) {
    probe_type <- detect_probetype(organism, probes, orgdb = NULL)
  }
  if (is.null(probe_type)) {
    message("ERROR: could not determine probe_type ")
    return(NULL)
  }

  ## --------------------------------------------
  ## retrieve table
  ## --------------------------------------------
  cols <- c("SYMBOL", "GENENAME", "GENETYPE", "MAP")
  cols <- intersect(cols, AnnotationDbi::keytypes(orgdb))

  cat("get gene annotation columns:", cols, "\n")
  message("retrieving annotation for ", length(probes), " ", probe_type, " features...")

  suppressMessages(suppressWarnings(
    annot <- AnnotationDbi::select(orgdb,
      keys = probes, columns = cols,
      keytype = probe_type
    )
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

  # Dog id should be canis lupus familiaris and not cannis familiaris, as in ah
  if (organism == "Canis familiaris") {
    org_orthogene <- "Canis lupus familiaris"
  } else {
    org_orthogene <- organism
  }

  ## get human ortholog using 'orthogene'
  ortho.map <- orthogene::map_species(method = "gprofiler")
  head(ortho.map)
  cat("\ngetting human orthologs...\n")
  has.ortho <- org_orthogene %in% ortho.map$scientific_name
  has.symbol <- "SYMBOL" %in% colnames(annot)
  if (organism == "Homo sapiens") {
    annot$ORTHOGENE <- annot$SYMBOL
  } else if (has.ortho && has.symbol) {
    ortho.out <- orthogene::convert_orthologs(
      gene_df = unique(annot$SYMBOL),
      input_species = org_orthogene,
      output_species = "human",
      non121_strategy = "drop_both_species",
      method = "gprofiler"
    )

    if (dim(ortho.out)[1] == 0) {
      ortho.out <- orthogene::convert_orthologs(
        gene_df = unique(annot$SYMBOL),
        input_species = org_orthogene,
        output_species = "human",
        non121_strategy = "drop_both_species",
        method = "homologene"
      )
    }

    if (dim(ortho.out)[1] == 0) {
      ortho.out <- orthogene::convert_orthologs(
        gene_df = unique(annot$SYMBOL),
        input_species = org_orthogene,
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
  annot$SOURCE <- "OrgDb"
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

  if (!all(probes0 %in% genes$feature)) {
    message("WARNING: not all probes could be annotated")
  }
  genes <- genes[match(probes0, genes$feature), , drop = FALSE]

  if (is.null(genes)) {
    warning("[getGeneAnnotation] ERROR : could not create gene annotation")
    return(NULL)
  }

  ## in case there were duplicated probe names we must make them unique
  rownames(genes) <- make_unique(probes0) ## in pgx-functions.R

  ## clean up
  genes <- cleanupAnnotation(genes)

  return(genes)
}


#' Cleanup annotation
#'
cleanupAnnotation <- function(genes) {
  ## add missing columns if needed, then reorder
  columns <- c(
    "feature", "symbol", "human_ortholog", "gene_title", "gene_biotype",
    "map", "chr", "pos", "tx_len", "source", "gene_name"
  )
  missing.cols <- setdiff(columns, colnames(genes))
  missing.cols
  for (a in missing.cols) genes[[a]] <- NA
  genes <- genes[, columns]
  colnames(genes) <- columns

  # gene_name should ALWAYS be assigned to feature for compatibility
  # with gene_name legacy implementation
  genes$gene_name <- genes$feature

  # add space after ; to conform with playbase <= 1.3.2
  genes$gene_title <- gsub(";[ ]*", "; ", genes$gene_title)

  # rename protein-coding to protein_coding to confirm with playbase <= v1.3.2
  genes$gene_biotype <- sub("protein-coding", "protein_coding", genes$gene_biotype)

  # replace NA in gene_ortholog by "" to conform with old pgx objects
  genes$human_ortholog[is.na(genes$human_ortholog)] <- ""

  # if organism is human, human_ortholog should be NA (matching old playbase annot)
  if (is.null(genes$human_ortholog)) {
    genes$human_ortholog <- NA
  }

  genes
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
#' pgx <- getCustomAnnotation(counts, custom_annot)
#' }
#' @export
getCustomAnnotation <- function(probes, custom_annot = NULL) {
  message("[getCustomAnnotation] Adding custom annotation table...")
  # If the user has provided a custom gene table, check it and use it

  num_annot <- sum(probes %in% custom_annot$feature)

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

  if (!is.null(custom_annot) && num_annot > 1 && required_in_annot) {
    # remove all NA columns, otherwise the for loop below will not work
    custom_annot <- custom_annot[, !apply(custom_annot, 2, function(x) all(is.na(x)))]

    # identify missing columns and fill them with annot_map
    missing_cols <- setdiff(names(annot_map), names(custom_annot))

    custom_annot[missing_cols] <- annot_map[missing_cols]

    # filter annotated table by probes using match
    custom_annot <- custom_annot[match(probes, custom_annot$feature), ]

    # if row was missing from annotation table (NA from match call above), input NA based on probes

    rownames(custom_annot) <- probes

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
    # Create custom gene table from probe names
    message("[getCustomAnnotation] Creating annotation table from probe names...")
    custom_annot <- data.frame(
      feature = probes,
      symbol = probes,
      gene_name = probes,
      human_ortholog = "",
      gene_title = "unknown",
      gene_biotype = "unknown",
      chr = "unknown",
      pos = 0,
      tx_len = 0,
      map = "1",
      source = "custom"
    )
    rownames(custom_annot) <- probes
  }

  custom_annot <- custom_annot[, table_col_order]
  custom_annot <- cleanupAnnotation(custom_annot)
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
.getOrgDb <- function(organism, use.ah = NULL) {
  if (tolower(organism) == "human") organism <- "Homo sapiens"
  if (tolower(organism) == "mouse") organism <- "Mus musculus"
  if (tolower(organism) == "rat") organism <- "Rattus norvegicus"
  organism
  
  if(is.null(use.ah) || !use.ah) {
    if(organism == "Homo sapiens" && require("org.Hs.eg.db", quietly=TRUE)) {
      message("[getOrgDb] returning org.Hs.eg.db for '", organism, "'\n")
      return(org.Hs.eg.db::org.Hs.eg.db)
    }
    if(organism == "Mus musculus" && require("org.Mm.eg.db", quietly=TRUE)) {
      message("[getOrgDb] returning org.Mm.eg.db for '", organism, "'\n")    
      return(org.Mm.eg.db::org.Mm.eg.db)
    }
    if(organism == "Rattus norvegicus" && require("org.Rn.eg.db", quietly=TRUE)) {
      message("[getOrgDb] returning org.Rn.eg.db for '", organism, "'\n")
      return(org.Rn.eg.db::org.Rn.eg.db)
    }
  }
  
  ah <- AnnotationHub::AnnotationHub()
  all_species <- getAllSpecies(ah)
  if (!tolower(organism) %in% tolower(all_species)) {
    message("WARNING: organism '", organism, "' not in AnnotationHub")
    return(NULL)
  }

  ## correct capitalization
  species <- all_species[which(tolower(all_species) == tolower(organism))]

  message("querying AnnotationHub for '", organism, "'\n")
  suppressMessages({
    ahDb <- AnnotationHub::query(ah, pattern = c(organism, "OrgDb"))
  })

  ## select on exact organism name
  ahDb <- ahDb[which(tolower(ahDb$species) == tolower(organism))]
  k <- length(ahDb) ## latest of multiple
  message("selecting database for '", ahDb$species[k], "'\n")
  
  message("retrieving annotation...\n")
  orgdb <- tryCatch({
    ahDb[[k]]
  },
  error = function(e) {
    message("An error occurred: ", e, ". Retrying with force=TRUE.")
    ahDb[[k, force = TRUE]]
  })

  return(orgdb)
}

#' @export
getOrgDb <- function(organism, use.ah = NULL) {
  if (tolower(organism) == "human") organism <- "Homo sapiens"
  if (tolower(organism) == "mouse") organism <- "Mus musculus"
  if (tolower(organism) == "rat") organism <- "Rattus norvegicus"
  organism
  orgdb <- .getOrgDb(organism, use.ah = use.ah) 
  if(is.null(orgdb)) return(NULL)
  
  suppressMessages({
    check.org <- grep("ORGANISM",capture.output(orgdb),value=TRUE)
  })
  check.org <- sub(".*ORGANISM: ","",check.org)
  check.org
  if( is.null(check.org) || check.org != organism ) {
    message("[getOrgDb] ***WARNING***: AnnotationHub is corrupt! removing cache")
    ah <- AnnotationHub::AnnotationHub(localHub=TRUE)
    AnnotationHub::removeCache(ah, ask=FALSE)
    orgdb <- .getOrgDb(organism, use.ah = use.ah) 
  }
  orgdb
}



#' @title Detect probe type from probe set
#' @export
detect_probetype <- function(organism, probes, orgdb = NULL,
                             nprobe = 100, use.ah = NULL) {
  if (tolower(organism) == "human") organism <- "Homo sapiens"
  if (tolower(organism) == "mouse") organism <- "Mus musculus"
  if (tolower(organism) == "rat") organism <- "Rattus norvegicus"
  
  ## get correct OrgDb database for organism
  if(is.null(orgdb)) {
    orgdb <- getOrgDb(organism, use.ah = use.ah)
  }
  if (is.null(orgdb)) {
    message("[detect_probetype] ERROR: unsupported organism '", organism, "'\n")
    return(NULL)
  }

  ## get probe types for organism
  keytypes <- c(
    "SYMBOL", "ENSEMBL", "UNIPROT", "ENTREZID",
    "GENENAME", "MGI",
    "ENSEMBLTRANS", "ENSEMBLPROT",
    "ACCNUM", "REFSEQ"
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

  ## discard isoform if UNIPROT
  is.uniprot <- mean(grepl("^[QP][0-9]*", probes)) > 0.8
  if (is.uniprot) probes <- sub("-[0-9]+", "", probes)

  ## Subset probes if too many
  if (length(probes) > nprobe) {
    if (nprobe > length(probes)) nprobe <- length(probes)
    # get random probes for query
    probes <- sample(probes, nprobe)
  }

  # Iterate over probe types
  key <- keytypes[1]
  for (key in keytypes) {
    probe_matches <- data.frame(NULL)

    # add symbol and genename on top of key as they will be used to
    # count the real number of probe matches
    key2 <- c(key, c("SYMBOL", "GENENAME"))
    key2 <- intersect(key2, keytypes)
    suppressMessages(suppressWarnings(try(
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
    if (matchratio > 0.95) break()
  }
  
  ## Return top match
  ##  key_matches
  top_match <- NULL
  if (all(key_matches == 0)) {
    message("WARNING:: Probe type not found. Valid probe types: ", paste(keytypes, collapse = " "))
    return(NULL)
  } else {
    top_match <- names(which.max(key_matches))
    message("Guessed probe type = ", top_match)
  }

  return(top_match)
}

#' @title Get human ortholog from given symbols of organism by using
#'   orthogene package. This package needs internet connection.
#'
#' @export
getHumanOrtholog <- function(organism, symbols) {
  ## test if orthogene server is reachable
  res <- try(orthogene::map_genes("CDK1"))
  if ("try-error" %in% class(res)) {
    message("[getHumanOrtholog] failed to contact server")
    df <- data.frame(symbols, "human" = NA)
    colnames(df)[1] <- organism
    rownames(df) <- NULL
    return(NULL)
  }

  # Dog id should be 'Canis lupus familiaris' (in orthogene) and not
  # Canis familiaris (in AH)
  if (organism == "Canis familiaris") {
    organism <- "Canis lupus familiaris"
  }

  clean2 <- function(s) {
    paste(head(strsplit(s, split = "[ _-]")[[1]], 2), collapse = " ")
  }

  organism <- gsub("[_]", " ", organism) ## orgDB names have sometimes underscores
  ortho.methods <- c("gprofiler", "homologene", "babelgene")
  has.ortho <- sapply(ortho.methods, function(m) !is.null(orthogene::map_species(organism, method = m, verbose = FALSE)))
  has.ortho

  ## build clean species list (only 'Genus species')
  ortho.species <- lapply(ortho.methods, function(m) {
    orthogene::map_species(method = m, verbose = FALSE)$scientific_name
  })
  ortho.clean <- lapply(ortho.species, function(s) {
    g <- sapply(s, clean2)
    names(g) <- s
    g
  })
  organism2 <- clean2(organism)
  has.clean <- sapply(ortho.clean, function(s) organism2 %in% s)
  has.clean

  ## build Genus list
  genus <- strsplit(organism, split = " ")[[1]][1]
  ortho.genus <- lapply(ortho.species, function(s) {
    g <- gsub(" .*|\\[|\\]", "", s)
    names(g) <- s
    g
  })
  has.genus <- sapply(ortho.genus, function(s) genus %in% s)
  has.genus

  ## select orthogene method (DB) and best matching species. If a
  ## species does not match, then we take the first matching species
  ## with the same Genus that is in the DB.
  if (any(has.ortho)) {
    sel <- head(which(has.ortho), 1)
    method <- ortho.methods[sel]
    ortho_organism <- orthogene::map_species(organism, method = method, verbose = FALSE)
    message("auto matching: ", organism)
  } else if (any(has.clean)) {
    ## if no exact species match, we try a clean2
    sel <- head(which(has.clean), 1)
    method <- ortho.methods[sel]
    ortho_organism <- head(names(which(ortho.clean[[sel]] == organism2)), 1)
    message("cleaned matching: ", organism2)
  } else if (any(has.genus)) {
    ## if no species match, we take the first Genus hit
    sel <- head(which(has.genus), 1)
    method <- ortho.methods[sel]
    ortho_organism <- head(names(which(ortho.genus[[sel]] == genus)), 1)
    message("genus matching: ", genus)
  } else {
    ## no match
    method <- NULL
    ortho_organism <- NULL
  }

  message("method = ", method)
  message("organism = ", organism)
  message("ortho_organism = ", ortho_organism)

  human.genes <- playdata::GENE_SYMBOL
  looks.human <- mean(toupper(symbols) %in% human.genes)
  looks.human
  message("looks.human = ", looks.human)

  if (organism == "Homo sapiens") {
    orthogenes <- symbols
  } else if (!is.null(method)) {
    ortho.out <- orthogene::convert_orthologs(
      gene_df = unique(symbols),
      input_species = ortho_organism,
      output_species = "human",
      non121_strategy = "drop_both_species",
      method = method,
      verbose = FALSE
    )
    ii <- match(symbols, ortho.out$input_gene)
    orthogenes <- rownames(ortho.out)[ii]
  } else if (looks.human > 0.5) {
    message("WARNING: symbols look 'human', using capitalized symbols")
    orthogenes <- toupper(symbols)
    orthogenes[which(!orthogenes %in% human.genes)] <- NA
  } else {
    message("WARNING: ", organism, " not found in orthogene database.")
    orthogenes <- rep(NA, length(symbols))
  }

  df <- data.frame(symbols, "human" = orthogenes)
  colnames(df)[1] <- organism
  return(df)
}


#' @title Show some probe types for selected organism
#'
#' @export
showProbeTypes <- function(organism, keytypes = NULL, use.ah = NULL, n = 10) {
  if (tolower(organism) == "human") organism <- "Homo sapiens"
  if (tolower(organism) == "mouse") organism <- "Mus musculus"
  if (tolower(organism) == "rat") organism <- "Rattus norvegicus"
  organism

  message(paste("retrieving probe types for", organism, "..."))

  ## get correct OrgDb database for organism
  orgdb <- getOrgDb(organism, use.ah = use.ah)
  if (is.null(orgdb)) {
    message("ERROR: unsupported organism '", organism, "'\n")
    return(NULL)
  }

  ## get probe types for organism
  if (is.null(keytypes)) {
    keytypes <- c(
      "SYMBOL", "ENSEMBL", "UNIPROT", "ENTREZID",
      "GENENAME", "MGI",
      "ENSEMBLTRANS", "ENSEMBLPROT",
      "ACCNUM", "REFSEQ"
    )
  }
  keytypes0 <- keytypes
  keytypes <- intersect(keytypes, keytypes(orgdb))

  if (length(keytypes) == 0) {
    message("ERROR: no valid keytypes in: ", keytypes0)
    return(NULL)
  }

  ## example probes
  keytype0 <- "ENTREZID"
  probes <- try(head(keys(orgdb, keytype = keytype0), n))
  if("try-error" %in% class(probes)) {
    keytype0 <- "SYMBOL"
    probes <- try(head(keys(orgdb, keytype = keytype0), n))
  }
  
  ## Iterate over probe types
  key_matches <- list()
  key <- keytypes[1]
  for (key in keytypes) {
    ## add symbol and genename on top of key as they will be used to
    ## count the real number of probe matches
    probe_matches <- try(
      suppressMessages(suppressWarnings(
        AnnotationDbi::select(
          orgdb,
          keys = probes,
          keytype = keytype0,
          columns = key
        ))),
      silent = TRUE
    )
    if(!"try-error" %in% class(probe_matches)) {
      ## set empty character to NA, as we only count not-NA to define probe type
      types <- probe_matches[, key]
      types <- setdiff(types, c("", NA))
      key_matches[[key]] <- head(types, n)
    }
  }

  return(key_matches)
}


#' @title Get all species in AnnotationHub/OrgDB
#'
#' @export
getAllSpecies <- function(ah = NULL) {
  if (is.null(ah)) {
    ah <- AnnotationHub::AnnotationHub() ## make global??
  }
  db <- AnnotationHub::query(ah, "OrgDb")
  sort(unique(AnnotationHub::mcols(db)$species))
}


#' @title Get species table in AnnotationHub/OrgDB
#'
#' @export
getSpeciesTable <- function(ah = NULL) {
  if (is.null(ah)) {
    ah <- AnnotationHub::AnnotationHub(localHub=FALSE) ## make global??
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

#' Get GO gene sets for organism directly from
#' AnnotationHub/OrgDB. Restrict to genes as background.
#'
#' export
getOrganismGO <- function(organism, use.ah = NULL, orgdb = NULL) {
  if (tolower(organism) == "human") organism <- "Homo sapiens"
  if (tolower(organism) == "mouse") organism <- "Mus musculus"
  if (tolower(organism) == "rat") organism <- "Rattus norvegicus"

  ## Load the annotation resource.
  if(is.null(orgdb)) {
    orgdb <- getOrgDb(organism, use.ah = use.ah)
  }
  
  go.gmt <- list()
  AnnotationDbi::keytypes(orgdb)
  if (!"GOALL" %in% AnnotationDbi::keytypes(orgdb)) {
    cat("WARNING:: missing GO annotation in database!\n")
  } else {
    ## create GO annotets
    cat("\nCreating GO annotation from AnnotationHub...\n")
    ont_classes <- c("BP", "CC", "MF")
    k <- "BP"
    for (k in ont_classes) {
      suppressMessages(suppressWarnings(
        go_id <- AnnotationDbi::mapIds(orgdb,
          keys = k, keytype = "ONTOLOGY",
          column = "GO", multiVals = "list"
        )[[1]]
      ))
      go_id <- unique(go_id)
      suppressMessages(suppressWarnings(      
        sets <- AnnotationDbi::mapIds(orgdb,
          keys = go_id, keytype = "GOALL",
          column = "SYMBOL", multiVals = "list"
        )
      ))
      
      ## get GO title
      sets <- sets[which(names(sets) %in% keys(GO.db::GOTERM))]      
      go <- sapply(GO.db::GOTERM[names(sets)], Term)
      new_names <- paste0("GO_", k, ":", go, " (", sub("GO:", "GO_", names(sets)), ")")
      names(sets) <- new_names

      ## add to list
      go.gmt <- c(go.gmt, sets)
    }
  }
  go.gmt
}


## ==================== using orthogene =====================
## WIP. This seems much faster than AnnotHub. There are about 700
## species supported. Online connection to server is needed but we are
## already using remote AnnotHub and orthogene for ortholog
## matching. The advantage is that probe type detection is not needed
## because orthogene seems to detect is automatically.

#' export
getGeneAnnotation.ORTHOGENE <- function(
    organism,
    probes,
    verbose = TRUE) {
  species <- orthogene::map_species(organism, method = "gprofiler", verbose = FALSE)
  species

  if (is.null(species)) {
    message("ERROR: unknown organism ", organism)
    return(NULL)
  }

  gene.out <- orthogene::map_genes(
    genes = probes,
    species = species,
    verbose = FALSE
  )
  head(gene.out)
  gene.out <- gene.out[match(probes, gene.out$input), ]

  ortho.out <- orthogene::convert_orthologs(
    gene_df = probes,
    input_species = species,
    output_species = "human",
    non121_strategy = "drop_both_species",
    method = "gprofiler",
    verbose = FALSE
  )
  head(ortho.out)
  ortholog <- rownames(ortho.out)[match(probes, ortho.out$input_gene)]
  ortholog[grep("^NA[.][1-9]", ortholog)] <- NA

  df <- data.frame(
    feature = probes,
    symbol = gene.out$name,
    human_ortholog = ortholog,
    gene_title = sub(" \\[.*", "", gene.out$description),
    gene_biotype = NA,
    map = NA,
    chr = NA,
    pos = NA,
    tx_len = NA,
    source = gene.out$namespace,
    gene_name = probes
  )
  head(df)

  return(df)
}

#' export
getAllSpecies.ORTHOGENE <- function() {
  M <- orthogene::map_species(method = "gprofiler")
  M$scientific_name
}
