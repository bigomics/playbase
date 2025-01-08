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
    genes <- getGeneAnnotation(
      organism = organism,
      probes = probes
    )
  }

  ## if annot_table is provided we override our annotation and append
  ## extra columns
  if (!is.null(annot_table)) {
    genes <- genes[, setdiff(colnames(genes), colnames(annot_table))]
    annot_table <- annot_table[match(rownames(genes), rownames(annot_table)), ]
    genes <- cbind(genes, annot_table)
  }

  ## cleanup entries and reorder columns
  genes <- cleanupAnnotation(genes)

  # Add to pgx object
  pgx$genes <- genes

  return(pgx)
}

# old function call
ngs.getGeneAnnotation <- function(...) {
  getGeneAnnotation(...)
}

#' @export
getGeneAnnotation <- function(
    organism,
    probes,
    use.ah = NULL,
    verbose = TRUE) {
  annot <- NULL

  if (tolower(organism) == "human") organism <- "Homo sapiens"
  if (tolower(organism) == "mouse") organism <- "Mus musculus"
  if (tolower(organism) == "rat") organism <- "Rattus norvegicus"
  if (tolower(organism) == "dog") organism <- "Canis familiaris"

  ## first annotate with ANNOTHUB
  info("[getGeneAnnotation] annotating with ANNOTHUB")
  annot <- getGeneAnnotation.ANNOTHUB(
    organism = organism,
    probes = probes,
    use.ah = use.ah,
    verbose = verbose
  )

  ## fallback with ORTHOGENE
  missing <- (is.na(annot$symbol) | annot$symbol == "")
  if (is.null(annot) || any(missing)) {
    info(
      "[getGeneAnnotation] annotating", sum(missing),
      "missing features with ORTHOGENE"
    )
    if (any(missing)) {
      missing.probes <- probes[which(missing)]
    } else {
      missing.probes <- probes
    }
    missing.annot <- try(getGeneAnnotation.ORTHOGENE(
      organism = organism,
      probes = missing.probes,
      verbose = verbose
    ))
    if (!"try-error" %in% class(missing.annot) && nrow(missing.annot)) {
      ## replace missing entries
      missing.annot <- missing.annot[, colnames(annot)]
      jj <- match(missing.probes, probes)
      annot[jj, ] <- missing.annot
    }
  }

  ## clean up
  annot <- cleanupAnnotation(annot)

  return(annot)
}

getOrthoSpecies <- function(organism) {
  S <- playbase::SPECIES_TABLE
  df <- data.frame(rownames(S), S[, c("species", "species_name", "ortho_species")])
  match <- colSums(apply(df, 2, tolower) == tolower(organism), na.rm = TRUE)
  if (all(match == 0)) {
    return(NULL)
  }
  k <- which.max(match)
  sel <- match(tolower(organism), tolower(df[, k]))
  if (length(sel) == 0) {
    return(NULL)
  }
  df[sel, "ortho_species"]
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
#' @details This function queries AnnotHub to retrieve key gene annotation data for
#' a set of input gene/transcript identifiers. It can detect the identifier
#' type automatically if not provided.
#'
#'
#' @examples
#' \dontrun{
#' probes <- c("ENSG00000142192", "ENST00000288602")
#' result <- getGeneAnnotation(organism, probes)
#' head(result)
#' }
#' @export
getGeneAnnotation.ANNOTHUB <- function(
    organism,
    probes,
    use.ah = NULL,
    probe_type = NULL,
    verbose = TRUE) {
  if (is.null(organism)) {
    warning("[getGeneAnnotation.ANNOTHUB] Please specify organism")
    return(NULL)
  }

  if (verbose) {
    message("[getGeneAnnotationy] Retrieving gene annotation...")
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

  if (is.null(probe_type)) {
    probe_type <- detect_probetype(organism, probes, orgdb = NULL)
    if (is.null(probe_type)) {
      message("ERROR: could not determine probe_type. Please specify. ")
      annot <- data.frame(feature = probes, symbol = "")
      annot <- cleanupAnnotation(annot)
      return(annot)
    }
  }

  ## clean up probe names from suffixes
  ## probes <- clean_probe_names(probes)
  probes <- match_probe_names(probes, orgdb, probe_type)

  ## --------------------------------------------
  ## retrieve table
  ## --------------------------------------------
  cols <- c("SYMBOL", "GENENAME", "GENETYPE", "MAP")
  cols <- intersect(cols, AnnotationDbi::keytypes(orgdb))

  if (organism %in% c("Mus musculus", "Rattus norvegicus")) {
    cols <- unique(c(cols, "ENTREZID"))
  }

  cat("get gene annotation columns:", cols, "\n")
  message("retrieving annotation for ", length(probes), " ", probe_type, " features...")

  suppressMessages(suppressWarnings(
    annot <- AnnotationDbi::select(
      orgdb,
      keys = probes,
      columns = cols,
      keytype = probe_type
    )
  ))

  ## Attempt to retrieve chr map via org.Mm.egCHRLOC / org.Rn.egCHRLOC.
  if (organism %in% c("Mus musculus", "Rattus norvegicus")) {
    if (organism == "Mus musculus") {
      library(org.Mm.eg.db)
      chrloc <- org.Mm.egCHRLOC
    }
    if (organism == "Rattus norvegicus") {
      library(org.Rn.eg.db)
      chrloc <- org.Rn.egCHRLOC
    }
    mapped_genes <- as.list(chrloc[mappedkeys(chrloc)])
    cm <- intersect(as.character(annot$ENTREZID), names(mapped_genes))
    mapped_genes <- mapped_genes[cm]
    locs <- unlist(lapply(mapped_genes, function(x) names(x[1])))
    jj <- match(names(locs), annot$ENTREZID)
    annot$MAP <- NA
    annot$MAP[jj] <- unname(locs)
    cls <- setdiff(colnames(annot), "ENTREZID")
    annot <- annot[, cls, drop = FALSE]
  }

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
  ## second pass for missing symbols
  ## --------------------------------------------
  is.missing <- (is.na(annot$SYMBOL) | annot$SYMBOL == "")
  missing.probes <- probes[which(is.missing)] ## probes match annot!
  missing.probes <- missing.probes[!is.na(missing.probes)]
  if (length(missing.probes)) {
    dbg(
      "[getGeneAnnotation.ANNOTHUB] retrying missing",
      length(missing.probes), "symbols..."
    )
    suppressWarnings(suppressMessages(
      missing.probe_type <- detect_probetype(organism, missing.probes, orgdb = orgdb)
    ))
    missing.probe_type
    if (!is.null(missing.probe_type)) {
      missing.probes1 <- match_probe_names(missing.probes, orgdb, missing.probe_type)
      suppressMessages(suppressWarnings(
        missing.annot <- AnnotationDbi::select(orgdb,
          keys = missing.probes1,
          columns = cols,
          keytype = missing.probe_type
        )
      ))
      head(missing.annot)
      missing.key <- missing.annot[, missing.probe_type]
      missing.annot$PROBE <- names(missing.probes[match(missing.key, missing.probes1)])
      jj <- match(missing.annot$PROBE, probes)
      colnames(missing.annot) <- colnames(annot)
      annot[jj, ] <- missing.annot
    }
  }

  ## get human ortholog using 'orthogene'
  cat("\ngetting human orthologs...\n")
  ortho_organism <- getOrthoSpecies(organism)
  annot$ORTHOGENE <- getHumanOrtholog(ortho_organism, annot$SYMBOL)$human

  ## Return as standardized data.frame and in the same order as input
  ## probes.
  pkgname <- orgdb$packageName
  if (length(pkgname) == 0) pkgname <- "OrgDb"
  annot$SOURCE <- pkgname[1]

  annot.cols <- c(
    "PROBE", "SYMBOL", "ORTHOGENE", "GENENAME",
    ##  "GENETYPE", "MAP", "CHR", "POS", "TXLEN", "SOURCE"
    "MAP", "SOURCE"
  )
  missing.cols <- setdiff(annot.cols, colnames(annot))
  missing.cols
  genes <- annot
  for (a in missing.cols) genes[[a]] <- NA
  genes <- genes[, annot.cols]
  new.names <- c(
    "feature", "symbol", "human_ortholog", "gene_title",
    ## "gene_biotype", "map", "chr", "pos", "tx_len", "source"
    "chr", "source"
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

  ## in case there were duplicated probe names we _must_ make them
  ## unique??? IK: really?? or should we remove duplicates?
  rownames(genes) <- make_unique(probes0) ## in pgx-functions.R

  return(genes)
}


#' Cleanup probe names from postfixes or version numbers
#'
#' @export
clean_probe_names <- function(probes, sep = "._") {
  probes0 <- probes
  probes[is.na(probes)] <- ""
  sel <- grep("[;]", probes)
  if (length(sel)) {
    probes[sel] <- sapply(strsplit(probes[sel], split = ";"), head, 1) ## take first
  }

  ## strip away anything after a 'dot' or 'underscore'
  probes <- sub(paste0("[", sep, "].*"), "", probes)

  ## is.ensembl <- mean(grepl("^ENS", probes)) > 0.5
  ## if (is.ensembl) {
  ##   probes <- sub("[.][0-9]+$", "", probes) ## strip version number
  ## }

  ## If UNIPROT we also strip isoform extension (orgDb does not like it)
  is.uniprot <- mean(grepl("^[QP][0-9]*", probes)) > 0.8
  if (is.uniprot) {
    ## probes <- sub("[.][0-9]+$", "", probes) ## strip phosphosite
    probes <- sub("-[0-9]+", "", probes) ## strip isoform
  }
  ##  names(probes) <- probes0
  probes
}

#' Match dirty probe names to clean key names
#'
#' @export
match_probe_names <- function(probes, org, probe_type = NULL) {
  if (is.character(org)) org <- getOrgDb(org)
  if (is.null(probe_type)) {
    probe_type <- detect_probetype(organism = "custom", probes, orgdb = org)
  }
  all.keys <- keys(org, probe_type)
  tsub <- function(s) gsub("[-:;.]|\\[|\\]", ".", s)
  ii <- match(toupper(tsub(probes)), toupper(tsub(all.keys)))
  table(is.na(ii))
  new.probes <- all.keys[ii]
  if (sum(is.na(new.probes))) {
    jj <- which(is.na(new.probes))
    new.probes[jj] <- probes[jj]
    jj.probes1 <- clean_probe_names(probes[jj], sep = "._")
    jj.probes2 <- clean_probe_names(probes[jj], sep = "._-")
    ii1 <- match(toupper(tsub(jj.probes1)), toupper(tsub(all.keys)))
    ii2 <- match(toupper(tsub(jj.probes2)), toupper(tsub(all.keys)))
    ii <- ifelse(!is.na(ii1), ii1, ii2)
    if (any(!is.na(ii))) {
      k <- which(!is.na(ii))
      new.probes[jj[k]] <- all.keys[ii[k]]
    }
  }
  names(new.probes) <- probes
  new.probes
}


#' Cleanup annotation
#'
cleanupAnnotation <- function(genes) {
  ## add missing columns if needed, then reorder
  columns <- c(
    "feature", "symbol", "human_ortholog", "gene_title", ## "gene_biotype",
    ## "map", "pos", "tx_len",
    "chr", "source", "gene_name"
  )
  missing.cols <- setdiff(columns, colnames(genes))
  missing.cols
  for (a in missing.cols) genes[[a]] <- NA
  #  genes <- genes[, columns]
  #  colnames(genes) <- columns

  # gene_name should ALWAYS be assigned to feature for compatibility
  # with gene_name legacy implementation
  genes$gene_name <- genes$feature

  # add space after ; to conform with playbase <= 1.3.2
  genes$gene_title <- gsub(";[ ]*", "; ", genes$gene_title)

  # rename protein-coding to protein_coding to confirm with playbase <= v1.3.2
  ## genes$gene_biotype <- sub("protein-coding", "protein_coding", genes$gene_biotype)

  # replace NA in gene_ortholog by "" to conform with old pgx objects
  genes$human_ortholog[is.na(genes$human_ortholog)] <- ""

  # if organism is human, human_ortholog should be NA (matching old playbase annot)
  if (is.null(genes$human_ortholog)) {
    genes$human_ortholog <- NA
  }

  ## reorder
  ordered.cols <- c(columns, setdiff(colnames(genes), columns))
  genes <- genes[, ordered.cols]

  ## Attempt: remove "pos", "tx_len"
  ##  keep <- colnames(genes)[!colnames(genes) %in% c("pos", "tx_len")]
  ##  genes <- genes[, keep]

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
    #    "gene_biotype" = "unknown",
    "chr" = "unknown",
    #    "pos" = 0,
    #    "tx_len" = 0,
    #    "map" = "1",
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
      #      gene_biotype = "unknown",
      chr = "unknown",
      #      pos = 0,
      #      tx_len = 0,
      #      map = "1",
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

  if (is.null(use.ah) || !use.ah) {
    if (organism == "Homo sapiens" && require("org.Hs.eg.db", quietly = TRUE)) {
      return(org.Hs.eg.db::org.Hs.eg.db)
    }
    if (organism == "Mus musculus" && require("org.Mm.eg.db", quietly = TRUE)) {
      return(org.Mm.eg.db::org.Mm.eg.db)
    }
    if (organism == "Rattus norvegicus" && require("org.Rn.eg.db", quietly = TRUE)) {
      return(org.Rn.eg.db::org.Rn.eg.db)
    }
    if (organism == "Plasmodium falciparum" && require("org.Pf.plasmo.db", quietly = TRUE)) {
      return(org.Pf.plasmo.db::org.Pf.plasmo.db)
    }
  }

  ah <- AnnotationHub::AnnotationHub()
  all_species <- allSpecies()
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
  orgdb <- tryCatch(
    {
      ahDb[[k]]
    },
    error = function(e) {
      message("An error occurred: ", e, ". Retrying with force=TRUE.")
      ahDb[[k, force = TRUE]]
    }
  )

  return(orgdb)
}

#'
#'
#' @export
getOrgDb <- function(organism, use.ah = NULL) {
  if (tolower(organism) == "human") organism <- "Homo sapiens"
  if (tolower(organism) == "mouse") organism <- "Mus musculus"
  if (tolower(organism) == "rat") organism <- "Rattus norvegicus"
  organism
  orgdb <- .getOrgDb(organism, use.ah = use.ah)
  if (is.null(orgdb)) {
    return(NULL)
  }

  ## Extra check for validity of database
  suppressMessages({
    check.org <- grep("ORGANISM", capture.output(orgdb), value = TRUE)
  })
  check.org <- sub(".*ORGANISM: ", "", check.org)
  check.org
  if (is.null(check.org) || check.org != organism) {
    message("[getOrgDb] ***WARNING***: AnnotationHub is corrupt! removing cache")
    ah <- AnnotationHub::AnnotationHub(localHub = TRUE)
    AnnotationHub::removeCache(ah, ask = FALSE)
    orgdb <- .getOrgDb(organism, use.ah = use.ah)
  }
  orgdb
}


#' @title Detect probe type from probe set
#' @export
detect_probetype <- function(organism, probes, orgdb = NULL,
                             nprobe = 1000, use.ah = NULL, datatype = NULL,
                             verbose = TRUE) {
  if (tolower(organism) == "human") organism <- "Homo sapiens"
  if (tolower(organism) == "mouse") organism <- "Mus musculus"
  if (tolower(organism) == "rat") organism <- "Rattus norvegicus"

  if (!is.null(datatype) && datatype == "metabolomics") {
    probe_type <- mx.detect_probetype(probes)
    return(probe_type)
  }

  if (!is.null(datatype) && datatype == "multi-omics") {
    mx.probes <- sub("^mx:","",grep("^mx:", probes, value=TRUE))
    px.probes <- sub("^px:","",grep("^px:", probes, value=TRUE))
    gx.probes <- sub("^gx:","",grep("^gx:", probes, value=TRUE))    
    gx.probe_types=px.probe_types=mx.probe_types=NULL
    if(length(gx.probes)) gx.probe_types <- detect_probetype(organism, gx.probes)
    if(length(px.probes)) px.probe_types <- detect_probetype(organism, px.probes)    
    if(length(mx.probes)) mx.probe_types <- mx.detect_probetype(mx.probes)    
    probe_type <- c(gx=gx.probe_types, px=px.probe_types, mx=mx.probe_types)
    return(probe_type)
  }

  ## get correct OrgDb database for organism
  if (is.null(orgdb)) {
    orgdb <- getOrgDb(organism, use.ah = use.ah)
  }
  if (is.null(orgdb)) {
    if (verbose) message("[detect_probetype] ERROR: unsupported organism '", organism, "'\n")
    return(NULL)
  }

  ## get probe types for organism
  keytypes <- c(
    "SYMBOL", "ENSEMBL", "ACCNUM", "UNIPROT", "GENENAME",
    "MGI", "TAIR", ## organism specific
    "ENSEMBLTRANS", "ENSEMBLPROT",
    "REFSEQ", "ENTREZID"
  )
  keytypes <- intersect(keytypes, keytypes(orgdb))
  key_matches <- rep(0L, length(keytypes))
  names(key_matches) <- keytypes

  ## clean up probes
  probes <- probes[!is.na(probes) & probes != ""]
  probes <- sapply(strsplit(probes, split = ";"), head, 1) ## take first
  probes <- unique(probes)

  ## Subset probes if too many
  if (length(probes) > nprobe) {
    if (nprobe > length(probes)) nprobe <- length(probes)
    # get random probes for query
    probes <- sample(probes, nprobe)
  }

  probes0 <- probes
  ## try different cleaning methods
  probes <- clean_probe_names(probes)
  probes1 <- clean_probe_names(probes, sep = "._-")
  probesx <- unique(c(probes0, probes, probes1))

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
        keys = probesx,
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
    matchratio <- max(n1, n2) / (1e-4 + nrow(probe_matches))
    key_matches[key] <- matchratio

    ## stop search prematurely if matchratio > 99%
    if (matchratio > 0.99) break()
  }

  ## Return top match
  ##  key_matches
  top_match <- NULL
  if (all(key_matches == 0)) {
    if (verbose) {
      message("head.probes = ", paste(head(probes), collapse = " "))
      message("WARNING: Probe type not found. Valid probe types: ", paste(keytypes, collapse = " "))
    }
    return(NULL)
  } else {
    top_match <- names(which.max(key_matches))
  }

  return(top_match)
}

#' @export
collapse_by_humansymbol <- function(obj, annot) {
  annot <- cbind(annot, rownames = rownames(annot))
  target <- c("human_ortholog", "symbol", "gene_name", "rownames")
  target <- intersect(target, colnames(annot))
  complete_targets <- lapply(target, function(x) {
    sum(is.na(annot[, x]) | annot[, x] %in% c("")) < 1
  }) |> unlist()
  target <- target[complete_targets]
  if (length(target) == 0) {
    message("[map_humansymbol] WARNING: could not find symbol mapping column.")
    return(obj)
  } else {
    ## call rename_by with target column
    map.obj <- rename_by(obj, annot_table = annot, new_id = target[1])
  }
  if (!is.null(dim(map.obj))) rownames(map.obj) <- toupper(rownames(map.obj))
  if (is.null(dim(map.obj))) names(map.obj) <- toupper(names(map.obj))
  map.obj
}


#' @title Get human ortholog from given symbols of organism by using
#'   orthogene package. This package needs internet connection.
#'
#' @export
getHumanOrtholog <- function(organism, symbols) {
  ## test if orthogene server is reachable
  res <- try(orthogene::map_genes("CDK1", verbose = FALSE))
  if ("try-error" %in% class(res)) {
    message("[getHumanOrtholog] failed to contact server")
    df <- data.frame(symbols, "human" = NA)
    colnames(df)[1] <- organism
    rownames(df) <- NULL
    return(NULL)
  }

  ## map to correct orthogene species name, if not
  ## done. SPECIES_TABLE$species are annothub names,
  ## SPECIES_TABLE$ortho_species are matched orthogene/gprofiler
  ## names.
  ortho_organism <- getOrthoSpecies(organism)

  orthogenes <- NULL
  ortho.out <- try(orthogene::convert_orthologs(
    gene_df = c("---", unique(symbols[!is.na(symbols)])),
    input_species = ortho_organism,
    output_species = "human",
    non121_strategy = "drop_both_species",
    method = "gprofiler",
    verbose = FALSE
  ))

  if (!"try-error" %in% class(ortho.out)) {
    ii <- match(symbols, ortho.out$input_gene)
    orthogenes <- rownames(ortho.out)[ii]
  }

  if (is.null(orthogenes)) {
    message("WARNING: could not find orthogene for ", organism)
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
  if (!is.null(keytypes) && keytypes[1] == "*") {
    keytypes <- keytypes(orgdb)
  }
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
  keytypes

  if (length(keytypes) == 0) {
    message("ERROR: no valid keytypes in: ", keytypes0)
    return(NULL)
  }

  ## example probes
  keytype0 <- "ENTREZID"
  suppressMessages(suppressWarnings(
    probes <- try(head(keys(orgdb, keytype = keytype0), n))
  ))
  if ("try-error" %in% class(probes)) {
    keytype0 <- setdiff(keytypes, "ENTREZID")[1]
    probes <- try(head(keys(orgdb, keytype = keytype0), n))
  }
  keytype0

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
        )
      )),
      silent = TRUE
    )
    if (!"try-error" %in% class(probe_matches)) {
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
allSpecies <- function(col = c("species", "ortho_species", "species_name")[1]) {
  M <- data.frame(playbase::SPECIES_TABLE)
  species <- as.character(M[, col])
  names(species) <- M[, "taxonomyid"]
  species
}

#' Return all species that are supported by the ANNOTHUB annotation
#' engine.
#'
#' @return character vector of species names
#'
#' @export
allSpecies.ANNOTHUB <- function() {
  M <- getSpeciesTable(ah = NULL)
  M <- data.frame(M)
  M <- M[M$rdataclass == "OrgDb", ]
  species <- as.character(M[, "species"])
  names(species) <- M[, "taxonomyid"]
  species
}

#' Return all species that are supported by the ORTHOGENE annotation
#' engine.
#'
#' @return character vector of species names
#'
#' @export
allSpecies.ORTHOGENE <- function() {
  M <- orthogene::map_species(method = "gprofiler", verbose = FALSE)
  species <- M[, "scientific_name"]
  names(species) <- M[, "taxonomy_id"]
  species
}

#' @title Get species table in AnnotationHub/OrgDB
#'
#' @export
getSpeciesTable <- function(ah = NULL) {
  if (is.null(ah)) {
    ah <- AnnotationHub::AnnotationHub(localHub = FALSE) ## make global??
  }
  ah.tables <- AnnotationHub::query(ah, "OrgDb")

  variables <- c(
    "ah_id", "species", "description", "rdatadateadded", "rdataclass",
    "title", "taxonomyid", "coordinate_1_based", "preparerclass", "sourceurl",
    "dataprovider", "genome", "maintainer", "tags", "sourcetype"
  )
  variables <- c(
    "ah_id", "species", "description", "rdatadateadded", "rdataclass",
    "title", "taxonomyid", ## "coordinate_1_based", "preparerclass", "sourceurl",
    ## "dataprovider", "genome", "maintainer", "tags",
    "sourcetype"
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
#' @export
getOrganismGO <- function(organism, use.ah = NULL, orgdb = NULL) {
  if (tolower(organism) == "human") organism <- "Homo sapiens"
  if (tolower(organism) == "mouse") organism <- "Mus musculus"
  if (tolower(organism) == "rat") organism <- "Rattus norvegicus"

  ## Load the annotation resource.
  if (is.null(orgdb)) {
    orgdb <- getOrgDb(organism, use.ah = use.ah)
  }

  go.gmt <- list()
  AnnotationDbi::keytypes(orgdb)
  ont_classes <- c("BP", "CC", "MF")
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
      sets <- lapply(sets, function(s) unique(s))
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

#' @export
getGeneAnnotation.ORTHOGENE <- function(
    organism,
    probes,
    verbose = TRUE) {
  ## correct organism names different from OrgDb
  if (organism == "Canis familiaris") {
    organism <- "Canis lupus familiaris"
  }

  ## map given name to official species name
  ## species <- try(orthogene::map_species(organism, method = "gprofiler", verbose = FALSE))
  ## species
  ## if ("try-error" %in% class(species)) {
  ##   message("[getGeneAnnotation.ORTHOGENE] *WARNING* could not connect to server")
  ##   return(NULL)
  ## }
  species <- getOrthoSpecies(organism)
  if (is.null(species)) {
    message("ERROR: unknown organism ", organism)
    return(NULL)
  }

  probes1 <- clean_probe_names(probes)

  gene.out <- orthogene::map_genes(
    genes = probes1,
    species = species,
    verbose = FALSE
  )
  head(gene.out)
  gene.out <- gene.out[match(probes1, gene.out$input), ]

  ortholog <- getHumanOrtholog(organism, gene.out$name)$human
  genebuild <- gprofiler2::get_version_info()$genebuild

  df <- data.frame(
    feature = probes,
    symbol = gene.out$name,
    human_ortholog = ortholog,
    gene_title = sub(" \\[.*", "", gene.out$description),
    #    gene_biotype = NA,
    #    map = NA,
    chr = NA,
    #    pos = NA,
    #    tx_len = NA,
    # source = gene.out$namespace,
    source = "gprofiler2",
    gene_name = probes
  )

  return(df)
}

#' Check if probes of organism are automatically recognized by
#' ORTHOGENE annotation engine.
#'
#' @return TRUE  if probes are recognized
#' @return FALSE if probes are not recognized
#'
#' @export
check_probetype.ORTHOGENE <- function(organism, probes) {
  map <- try(orthogene::map_genes(probes, species = organism, drop_na = FALSE, verbose = FALSE))
  if ("try-error" %in% class(map) || is.null(map)) {
    message("[check_probetype.ORTHOGENE] *WARNING* could not connect to server")
    return(NULL)
  }
  mean.mapped <- mean(!is.na(map$target))
  ## get correct OrgDb database for organism
  if (mean.mapped < 0.20) {
    return(FALSE)
  }
  return(TRUE)
}


#' Check if probes can be detected by Orthogene or AnnotHub/OrgDb
#' annotation engines.
#'
#' export
check_probetype <- function(organism, probes) {
  chk1 <- check_probetype.ORTHOGENE(organism, probes)
  if (!is.null(chk1) && chk1 == TRUE) {
    return(TRUE)
  }
  ## using AnnotHub/OrgDb
  chk2 <- detect_probetype(organism, probes)
  if (!is.null(chk2)) {
    return(TRUE)
  }
  return(FALSE)
}

#' Create new feature name by concatenating some columns of input
#' annotation table. Make all feature names unique.
#'
#' @param annot  some annotation dataframe
#' @param target vector of character. e.g. c("feature","_","symbol")
#'
#' export
combine_feature_names <- function(annot, target) {
  new.feature <- annot[, 0]
  for (i in 1:length(target)) {
    if (target[i] == 0 || target[i] == "rownames") {
      new.feature <- paste0(new.feature, rownames(annot))
    } else if (target[i] %in% colnames(annot)) {
      new.feature <- paste0(new.feature, annot[, target[i]])
    } else {
      ## some character
      new.feature <- paste0(new.feature, target[i])
    }
  }
  if (sum(duplicated(new.feature)) > 0) {
    message("[merge feature names] duplicated = ", sum(duplicated(new.feature)))
    new.feature <- make_unique(new.feature)
  }
  new.feature
}

#' @export
pgx.getGeneInfo <- function(pgx, gene) {
  feature <- pgx$genes[match(gene, pgx$genes$symbol), "feature"]
  ortholog <- pgx$genes[match(gene, pgx$genes$symbol), "human_ortholog"]
  datatype <- ifelse(is.null(pgx$datatype), "rna-seq", pgx$datatype)
  getOrgGeneInfo(pgx$organism, gene, feature, ortholog, datatype, as.link = TRUE)
}

#' @export
getOrgGeneInfo <- function(organism, gene, feature, ortholog, datatype, as.link = TRUE) {
  if (is.null(gene) || length(gene) == 0) {
    return(NULL)
  }
  if (is.na(gene) || gene == "") {
    return(NULL)
  }

  orgdb <- getOrgDb(organism, use.ah = NULL)
  cols <- c("SYMBOL", "UNIPROT", "GENENAME", "MAP", "OMIM", "PATH", "GO")
  cols <- intersect(cols, keytypes(orgdb))

  if (!"SYMBOL" %in% cols) {
    keytype <- detect_probetype(organism, gene)
  } else {
    keytype <- "SYMBOL"
  }

  ## return if gene is not known
  if (!gene %in% keys(orgdb, keytype)) {
    info <- list()
    info[["feature"]] <- feature
    info[["symbol"]] <- gene
    info[["organism"]] <- organism
    info[["summary"]] <- "(no info available)"
    return(info)
  }

  ## get info from different environments
  info <- lapply(cols, function(k) {
    tryCatch(
      {
        AnnotationDbi::select(
          orgdb,
          keys = gene,
          keytype = keytype,
          columns = k
        )[[k]]
      },
      error = function(w) {
        NULL
      }
    )
  })

  #  if(is.null(unlist(info))){
  #    return(NULL)
  #  }
  names(info) <- cols

  info[["ORGANISM"]] <- organism

  ## take out duplicates
  info <- lapply(info, unique)
  symbol <- info[[keytype]]
  uniprot <- info[["UNIPROT"]]
  if (length(uniprot) == 0) {
    this.uniprot <- NULL
  } else {
    this.uniprot <- uniprot[which(sapply(uniprot, function(p) grepl(p, feature)))]
    if (length(this.uniprot) == 0) this.uniprot <- uniprot[1]
  }

  if (as.link && length(symbol)) {
    gene.link <- "<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=GENE' target='_blank'>GENE</a>"
    gene.link <- sapply(symbol, function(s) gsub("GENE", s, gene.link))
    info[["SYMBOL"]] <- paste(gene.link, collapse = ", ")
  }

  if (as.link && length(uniprot)) {
    prot.link <- "<a href='https://www.uniprot.org/uniprotkb/UNIPROT' target='_blank'>UNIPROT</a>"
    prot.link <- sapply(uniprot, function(s) gsub("UNIPROT", s, prot.link))
    info[["UNIPROT"]] <- paste(prot.link, collapse = ", ")
  }

  ## create link to external databases: OMIM, GeneCards, Uniprot
  if (as.link) {
    genecards.link <- "<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=GENE' target='_blank'>GeneCards</a>"
    uniprot.link <- "<a href='https://www.uniprot.org/uniprotkb/UNIPROT' target='_blank'>UniProtKB</a>"
    genecards.link <- NULL
    uniprot.link <- NULL
    if (length(symbol)) genecards.link <- sub("GENE", symbol[1], genecards.link)
    if (length(this.uniprot)) uniprot.link <- sub("UNIPROT", this.uniprot, uniprot.link)
    info[["databases"]] <- paste(c(genecards.link, uniprot.link), collapse = ", ")
  }

  if (length(this.uniprot) && grepl("proteomics", datatype, ignore.case = TRUE)) {
    ## create link to PhosphoSitePlus
    phosphositeplus.link <- "<a href='https://www.phosphosite.org/simpleSearchSubmitAction.action?searchStr=GENE' target='_blank'>PhosphoSitePlus</a>"
    phosphositeplus.link <- "<a href='https://www.phosphosite.org/uniprotAccAction?id=UNIPROT' target='_blank'>PhosphoSitePlus</a>"
    ## phosphositeplus.link <- sub("GENE", symbol[1], phosphositeplus.link)
    phosphositeplus.link <- sub("UNIPROT", this.uniprot, phosphositeplus.link)
    info[["databases"]] <- paste(c(info[["databases"]], phosphositeplus.link), collapse = ", ")

    ## ## create links to PhosphoELM for proten and gene: db of S/T/Y phosphorylation sites
    ## phosphoELM.link <- "<a href='http://phospho.elm.eu.org/byAccession/UNIPROT' target='_blank'>PhosphoELM</a>"
    ## phosphoELM.link <- sub("UNIPROT", uniprot, phosphoELM.link)
    ## info[["databases"]] <- paste(c(info[["databases"]], phosphoELM.link), collapse = ", ")
  }

  ## create link to OMIM
  if (as.link && length(info[["OMIM"]])) {
    omim.link <- "<a href='https://www.omim.org/entry/OMIM' target='_blank'>OMIM</a>"
    info[["OMIM"]] <- sapply(info[["OMIM"]], function(x) gsub("OMIM", x, omim.link))
  }

  ## create link to KEGG
  if (as.link && length(info[["PATH"]])) {
    kegg.link <- "<a href='https://www.genome.jp/kegg-bin/show_pathway?map=hsaKEGGID&show_description=show' target='_blank'>KEGGNAME (KEGGID)</a>"
    for (i in 1:length(info[["PATH"]])) {
      kegg.id <- info[["PATH"]][[i]]
      kegg.id <- setdiff(kegg.id, NA)
      if (length(kegg.id) > 0) {
        kegg.name <- AnnotationDbi::mget(kegg.id, envir = KEGG.db::KEGGPATHID2NAME, ifnotfound = NA)[[1]]
        if (!is.na(kegg.name) && as.link) {
          info[["PATH"]][[i]] <- gsub("KEGGNAME", kegg.name, gsub("KEGGID", kegg.id, kegg.link))
        } else {
          info[["PATH"]][[i]] <- kegg.name
        }
      }
    }
  }

  ## create link to GO
  if (length(info[["GO"]]) && !is.na(info[["GO"]][1])) {
    ## sometimes GO.db is broken...
    suppressWarnings(try.out <- try(AnnotationDbi::Term(AnnotationDbi::mget("GO:0000001",
      envir = GO.db::GOTERM,
      ifnotfound = NA
    )[[1]])))
    go.ok <- (class(try.out) != "try-error")
    if (go.ok) {
      amigo.link <- "<a href='http://amigo.geneontology.org/amigo/term/GOID' target='_blank'>GOTERM (GOID)</a>"
      i <- 1
      for (i in 1:length(info[["GO"]])) {
        go_id <- info[["GO"]][i]

        term_id <- AnnotationDbi::mget(go_id, envir = GO.db::GOTERM, ifnotfound = NA)[[1]]
        if (class(term_id) == "GOTerms") {
          go_term <- AnnotationDbi::Term(term_id)
          if (as.link) {
            info[["GO"]][i] <- gsub("GOTERM", go_term, gsub("GOID", go_id, amigo.link))
          } else {
            info[["GO"]][i] <- go_term
          }
        } else {
          info[["GO"]][i] <- go_id
        }
      }
    } else {
      info[["GO"]] <- NULL
    }
  }

  ## pull summary
  info[["SUMMARY"]] <- "(no info available)"
  ## ortholog <- getHumanOrtholog(organism, symbol)$human

  if (!is.null(ortholog) && ortholog %in% names(playdata::GENE_SUMMARY)) {
    info[["SUMMARY"]] <- playdata::GENE_SUMMARY[ortholog]
    info[["SUMMARY"]] <- gsub("Publication Note.*|##.*", "", info[["SUMMARY"]])
  }

  ## rename
  tags <- c(
    "ORGANISM", "SYMBOL", "UNIPROT", "GENENAME", "MAP", "OMIM", "PATH",
    "GO", "SUMMARY", "databases"
  )
  info <- info[tags]
  names(info) <- c(
    "organism", "gene_symbol", "uniprot", "name", "map_location",
    "OMIM", "pathway", "GO", "summary", "databases"
  )

  return(info)
}

#' @export
getMetaboliteInfo <- function(organism = "Human", chebi) {
  
  if (is.null(chebi) || length(chebi) == 0) {
    return(NULL)
  }
  if (is.na(chebi) || chebi == "") {
    return(NULL)
  }

  metabolite_metadata <- playdata::METABOLITE_METADATA
  annot <- playdata::METABOLITE_ANNOTATION

  inf <- list()
  inf[["name"]] <- metabolite_metadata[metabolite_metadata$ID == chebi, "name"]
  inf[["summary"]] <- metabolite_metadata[metabolite_metadata$ID == chebi, "definition"]
  inf[["organism"]] <- organism
  
  # remove summary if it is null
  if (is.null(inf[["summary"]])) inf[["summary"]] <- "Summary not available for this metabolite."
  if (inf[["summary"]] == "null") inf[["summary"]] <- "Summary not available for this metabolite."

  # get annotation for a given chebi id
  annotation <- annot[annot$ID == chebi, ]
  
  # remove NA columns from annotation
  annotation <- annotation[, colSums(is.na(annotation)) < nrow(annotation)]
  cols <- colnames(annotation)[-1] ## exclude chebi IDS as we already have it

  ## get info from different environments
  res <- lapply(cols, function(k) {
    link <- NULL
    matched_id <- annotation[annotation$ID == chebi, k]
    if (k == "HMDB") link <- glue::glue("<a href='https://hmdb.ca/metabolites/{matched_id}' target='_blank'>{matched_id}</a>")
    if (k == "KEGG") link <- glue::glue("<a href='https://www.kegg.jp/dbget-bin/www_bget?{matched_id}' target='_blank'>{matched_id}</a>")
    if (k == "PubChem") link <- glue::glue("<a href='https://pubchem.ncbi.nlm.nih.gov/compound/{matched_id}' target='_blank'>{matched_id}</a>")
    if (k == "ChEBI") link <- glue::glue("<a href='https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:{matched_id}' target='_blank'>{matched_id}</a>")
    if (k == "PATHBANK") link <- glue::glue("<a href='https://moldb.wishartlab.com/molecules/{matched_id}/curation.html' target='_blank'>{matched_id}</a>")    
#    if (k == "METLIN") link <- matched_id # METLIN is offline at the time of this request, needs to be updated
#    if (k == "SMILES") link <- matched_id
    return(link)
  })

  # merge all info
  names(res) <- cols
  inf <- c(inf, res)

  ## create link to external databases: these libraries are not always
  ## available for a given chebi id
  hmdb.link <- NULL
  kegg.link <- NULL
  pubchem.link <- NULL
  if (!is.null(inf[["HMDB"]])) hmdb.link <- glue::glue("<a href='https://hmdb.ca/metabolites/{annotation[,'HMDB']}' target='_blank'>HMDB</a>")
  if (!is.null(inf[["KEGG"]])) kegg.link <- glue::glue("<a href='https://www.kegg.jp/dbget-bin/www_bget?{annotation[,'KEGG']}' target='_blank'>KEGG</a>")
  if (!is.null(inf[["PubChem"]])) pubchem.link <- glue::glue("<a href='https://pubchem.ncbi.nlm.nih.gov/compound/{annotation[,'PubChem']}' target='_blank'>PubChem</a>")

  # these libraries are always available
  chebi.link <- glue::glue("<a href='https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:{chebi}' target='_blank'>ChEBI</a>")
  reactome.link <- glue::glue("<a href='https://reactome.org/content/query?q=chebi%3A{chebi}' target='_blank'>Reactome</a>")

  inf[["databases"]] <- paste(c(hmdb.link, kegg.link, reactome.link, pubchem.link, chebi.link), collapse = ", ")

  return(inf)
}

#' Automatically detects species by trying to detect probetype from
#' list of test_species. Warning. bit slow.
#'
#' @export
detect_species_probetype <- function(
    probes,
    test_species = c("human", "mouse", "rat"),
    datatype = NULL) {

  probes <- unique(clean_probe_names(probes))
  ptype <- list()
  s="human"
  for (s in test_species) {
    ptype[[s]] <- detect_probetype(
      organism = s,
      probes = probes,
      use.ah = FALSE,
      datatype = datatype,
      verbose = FALSE
    )
  }
  ptype <- unlist(ptype)
  species <- sub("[.].*","",names(ptype))  ## remove datatype postfix
  out <- tapply( ptype, species, c)
  ## out <- list(
  ##   species = species,
  ##   probetype = as.character(ptype)
  ## )
  return(out)
}

#' Annotate phosphosite with residue symbol. Feature names must be of
#' form 'uniprot_position'.
#'
#' @export
annotate_phospho_residue <- function(features, detect.only = FALSE) {
  valid_name <- mean(grepl("[_.][1-9].*", features), na.rm = TRUE) > 0.9
  uniprot <- sub("[_.][1-9].*", "", features)
  positions <- strsplit(sub(".*[_.]", "", features), split = "[;/,]")

  P <- playdata::PHOSPHOSITE
  prot.match <- mean(uniprot %in% P$UniProt, na.rm = TRUE)
  pos.match <- mean(positions %in% P$Position, na.rm = TRUE)
  is_phospho <- (valid_name && prot.match > 0.50 && pos.match > 0.50)

  if (detect.only) {
    return(is_phospho)
  }

  if (is_phospho) {
    ## determine separators
    sep1.match <- sapply(c("_", "."), function(s) {
      sum(grepl(s, features, fixed = TRUE), na.rm = TRUE)
    })
    sep1 <- names(which.max(sep1.match))
    sel <- grep("[;/,]", features)
    sep2.match <- sapply(c(";", "/", ","), function(s) {
      sum(grepl(s, features[sel], fixed = TRUE), na.rm = TRUE)
    })
    sep2 <- names(which.max(sep2.match))
    P <- P[which(P$UniProt %in% uniprot), ]
    dim(P)
    P.id <- paste0(P$UniProt, "_", P$Position)
    F.id <- lapply(1:length(uniprot), function(i) {
      paste0(uniprot[i], "_", positions[[i]])
    })

    ## this takes a while...
    p.idx <- lapply(uniprot, function(p) which(P$UniProt == p))
    type <- sapply(1:length(positions), function(i) {
      jj <- match(positions[[i]], P$Position[p.idx[[i]]])
      tt <- P$Residue[p.idx[[i]][jj]]
      tt[is.na(tt)] <- "" ## not found
      tt
    })

    new.features <- sapply(1:length(features), function(i) {
      tt <- type[[i]]
      pp <- paste(paste0(tt, positions[[i]]), collapse = sep2)
      paste0(uniprot[i], sep1, pp)
    })
    features <- new.features
  }
  features
}


#' Convert probetype unsing annothub
#'
#' @export
convert_probetype <- function(organism, probes, target_id, from_id = NULL,
                              datatype = NULL, orgdb = NULL, verbose = TRUE) {
  if (tolower(organism) == "human") organism <- "Homo sapiens"
  if (tolower(organism) == "mouse") organism <- "Mus musculus"
  if (tolower(organism) == "rat") organism <- "Rattus norvegicus"

  if (!is.null(datatype) && datatype == "metabolomics") {
    ## new.probes <- mx.detect_probetype(probes)
    ## return(new.probes)
    message("WARNING: metabolomics not yet implemented")
    return(NULL)
  }

  ## get correct OrgDb database for organism
  if (is.null(orgdb)) {
    orgdb <- getOrgDb(organism)
  }
  if (is.null(orgdb)) {
    if (verbose) message("[convert_probetype] ERROR: unsupported organism '", organism, "'\n")
    return(NULL)
  }

  if (!target_id %in% keytypes(orgdb)) {
    message("[convert_probetype] invalid target probetype")
    return(NULL)
  }
  if (is.null(from_id)) {
    from_id <- detect_probetype(organism, probes, orgdb = orgdb, datatype = NULL)
  }
  from_id
  message("[convert_probetype] converting from ", from_id, " to ", target_id)

  suppressMessages(suppressWarnings(try(
    res <- AnnotationDbi::select(
      orgdb,
      keys = probes,
      keytype = from_id,
      columns = target_id
    ),
    silent = TRUE
  )))
  new.probes <- res[match(probes, res[, from_id]), target_id]
  return(new.probes)
}


#' Convert multi-omics probetype. Probe names *must  be prefixed with
#' data type unless classical transcriptomics/proteomics.
#'
#' @export
getProbeAnnotation <- function(organism, probes) {
  
  if(sum(grepl("[:]",probes))) {
    dtype <- sub(":.*","",probes)
  } else {
    ## no colon in names. try to guess by matching.
    ptype <- detect_probetype(organism, probes)
    mtype <- mx.detect_probetype(probes)    
    gx.types <- c(
      "SYMBOL", "ENSEMBL", "ACCNUM", "GENENAME",
      "MGI", "TAIR",  "ENSEMBLTRANS", "REFSEQ", "ENTREZID"
    )
    px.types <- c("UNIPROT", "ENSEMBLPROT")
    dx <- ifelse(ptype %in% px.types, "px", "gx")
    dtype <- rep(dx,length(probes))
  }
  table(dtype)
  dtype <- tolower(dtype)
  dtype <- sub(paste0("ensembl|symbol|hugo|gene|hgnc",".*"),"gx",dtype)
  dtype <- sub(paste0("uniprot|protein",".*"),"px",dtype)
  dtype <- sub(paste0("chebi|hmdb|kegg|pubchem",".*"),"mx",dtype)  
  table(dtype)

  ## populate with defaults
  symbol <- toupper(sub(".*:","",probes))
  annot <- list()
  if(any(dtype %in% c('gx','px'))) {
    ii <- which(dtype %in% c('gx','px'))
    pp <- sub(".*:","",probes[ii])
    aa <- getGeneAnnotation(organism, pp)
    head(aa)
    aa$data_type <- sub(":.*","",probes[ii])
    rownames(aa) <- probes[ii]
    aa$feature <- probes[ii]
    annot <- c(annot, list(aa))
  }
  if("mx" %in% dtype) {
    ii <- which(dtype == 'mx')
    pp <- sub(".*:","",probes[ii])
    aa <- getMetaboliteAnnotation(pp)
    head(aa)
    aa$data_type <- 'mx'
    rownames(aa) <- probes[ii]
    aa$feature <- probes[ii]
    annot <- c(annot, list(aa))
  }

  cols <- Reduce( intersect, lapply(annot,colnames))
  annot <- lapply(annot, function(a) a[,cols])
  annot <- do.call( rbind, annot )
  annot <- annot[match(probes, annot$feature),]
  rownames(annot) <- probes
  head(annot)

  ## fill NA
  annot$feature <- ifelse(is.na(annot$feature), probes, annot$feature)
  annot$symbol <- ifelse(is.na(annot$symbol), symbol, annot$symbol)
  annot$gene_name <- ifelse(is.na(annot$gene_name), probes, annot$gene_name)
  annot$data_type <- ifelse(is.na(annot$data_type), dtype, annot$data_type)    

  head(annot)
  tail(annot)  
  annot
}
