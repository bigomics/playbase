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
detect_probe <- function(probes, mart = NULL, verbose = TRUE) {
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

  probe_check <- sapply(probe_types_to_check, FUN = function(x) {
    tryCatch(
      {
        tmp <- biomaRt::getBM(
          attributes = x,
          filters = x,
          values = subset_probes,
          mart = mart
        )
        Sys.sleep(5) # Sleep time to prevent bounce from ensembl for consecutive calls
        out <- nrow(tmp)
        return(out)
      },
      error = function(e) {
        return(0)
      }
    )
  })

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
ngs.getGeneAnnotation <- function(
    probes,
    organism,
    mart = NULL,
    probe_type = NULL,
    verbose = TRUE) {
  # Check mart
  if (is.null(mart)) {
    stop("Mart not found. Please specify a BioMart database to use.")
  }

  # Prepare inputs
  if (verbose) {
    message("[ngs.getGeneAnnotation] Filling genes information...")
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
pgx.gene_table <- function(pgx, organism) {
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
      probe_type <- detect_probe(probes, ensembl_dataset)
    }

    # Get gene table
    if (!is.na(probe_type) & !is.null(ensembl_dataset)) {
      genes <- ngs.getGeneAnnotation(
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
  pgx$all_genes <- all_genes
  pgx$probe_type <- probe_type

  return(pgx)
}


#' @title Detect probe type from probe set
#' @description Detects the most likely probe type (ENSEMBL, SYMBOL, etc)
#' for a set of gene/transcript identifiers. This function has been deprecated
#' and is mantained for fallback purposes.
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
#' type <- detect_probe_DEPRECATED(probes, organism)
#' }
#' @export
detect_probe_DEPRECATED <- function(probes, organism) {
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

  # Return top match
  if (all(key_matches == 0)) {
    stop("Probe type not found, please, check your probes")
  } else {
    top_match <- names(which.max(key_matches))
  }

  return(top_match)
}


#' @title Get gene annotation data
#' @description Retrieves gene annotation data from an organism-specific
#' annotation package. This function has been deprecated
#' and is mantained for fallback purposes.
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
#' annot <- ngs.getGeneAnnotation_DEPRECATED(probes, "ENSEMBL", organism = "Mouse")
#' }
#' @export
ngs.getGeneAnnotation_DEPRECATED <- function(probes, probe_type, organism) {
  # Get org database and columns request
  if (organism == "Human") {
    org_db <- org.Hs.eg.db::org.Hs.eg.db
    cols_req <- c("SYMBOL", "GENENAME", "CHR", "CHRLOC", "MAP", "GENETYPE")
  } else if (organism == "Mouse") {
    org_db <- org.Mm.eg.db::org.Mm.eg.db
    cols_req <- c("SYMBOL", "GENENAME", "CHR", "CHRLOC", "GENETYPE")
  } else if (organism == "Rat") {
    org_db <- org.Rn.eg.db::org.Rn.eg.db
    cols_req <- c("SYMBOL", "GENENAME", "CHR", "CHRLOC", "GENETYPE")
  }

  ## discard version numbers if ENSEMBL
  orig.probes <- probes
  if (mean(grepl("^ENS", probes)) > 0.8) {
    probes <- sub("[.][0-9]+", "", probes)
  }

  # Call for annotation table
  suppressWarnings(d <- AnnotationDbi::select(org_db,
    keys = probes,
    keytype = probe_type,
    columns = cols_req
  ))

  d <- data.table::as.data.table(d)

  # Add human ortholog and map for non-human organisms
  # Here, we use the old strategy of capitalise the symbol
  if (organism %in% c("Mouse", "Rat")) {
    d$human_ortholog <- toupper(d$SYMBOL)
    d$MAP <- d$CHR
  }

  ## if ENSEMBL get original probe names with version
  if (probe_type == "ENSEMBL") {
    d$ENSEMBL <- orig.probes[match(d$ENSEMBL, probes)]
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
  return(d[probes, , drop = FALSE])
}


##type=NULL;org="human";keep.na=FALSE
feature2symbol <- function(probes, type=NULL, org="human", keep.na=FALSE) {
    require(org.Hs.eg.db)
    require(org.Mm.eg.db)
    require(org.Rn.eg.db)

    ## strip postfix for ensemble codes
    if(mean(grepl("^ENS",probes))>0.5) {
        probes <- gsub("[.].*","",probes)
    }

    if(is.null(type)) {

        hs.list <- list(
            "human.ensembl" = unlist(as.list(org.Hs.egENSEMBL)),
            "human.ensemblTRANS" = unlist(as.list(org.Hs.egENSEMBLTRANS)),
            #"human.unigene" = unlist(as.list(org.Hs.egUNIGENE)),
            "human.refseq"  = unlist(as.list(org.Hs.egREFSEQ)),
            "human.accnum"  = unlist(as.list(org.Hs.egACCNUM)),
            "human.uniprot" = unlist(as.list(org.Hs.egUNIPROT)),
            "human.symbol"  = unlist(as.list(org.Hs.egSYMBOL))
            )

        mm.list <- list(
            "mouse.ensembl" = unlist(as.list(org.Mm.egENSEMBL)),
            "mouse.ensemblTRANS" = unlist(as.list(org.Mm.egENSEMBLTRANS)),
            #"mouse.unigene" = unlist(as.list(org.Mm.egUNIGENE)),
            "mouse.refseq"  = unlist(as.list(org.Mm.egREFSEQ)),
            "mouse.accnum"  = unlist(as.list(org.Mm.egACCNUM)),
            "mouse.uniprot" = unlist(as.list(org.Mm.egUNIPROT)),
            "mouse.symbol"  = unlist(as.list(org.Mm.egSYMBOL))
        )

        rn.list <- list(
          "rat.ensembl" = unlist(as.list(org.Rn.egENSEMBL)),
          "rat.ensemblTRANS" = unlist(as.list(org.Rn.egENSEMBLTRANS)),
          #"rat.unigene" = unlist(as.list(org.Rn.egUNIGENE)),
          "rat.refseq"  = unlist(as.list(org.Rn.egREFSEQ)),
          "rat.accnum"  = unlist(as.list(org.Rn.egACCNUM)),
          "rat.uniprot" = unlist(as.list(org.Rn.egUNIPROT)),
          "rat.symbol"  = unlist(as.list(org.Rn.egSYMBOL))
        )

        id.list <- c(hs.list, mm.list, rn.list)
        mx <- sapply(id.list, function(id) mean(probes %in% id))
        mx
        org=type=NULL
        max.mx <- max(mx,na.rm=TRUE)
        mx0 <- names(mx)[which.max(mx)]
        org  <- sub("[.].*","",mx0)
        type <- sub(".*[.]","",mx0)
        message("[probe2symbol] mapped ",format(100*max.mx,digits=2),"% of probes")
        if(max.mx < 0.5 && max.mx>0) {
            message("[probe2symbol] WARNING! low mapping ratio: r= ",max.mx)
        }
        if(max.mx==0) {
            message("[probe2symbol] WARNING! zero mapping ratio: r= ")
            type = NULL
        }
        org
        type
    }
    if(is.null(type)) {
        cat("probe2symbol: invalid type: ",type,"\n")
        return(NULL)
    }
    if(!type %in% c("ensembl","ensemblTRANS","unigene","refseq","accnum","uniprot","symbol")) {
        cat("probe2symbol: invalid type: ",type,"\n")
        return(NULL)
    }

    cat("[probe2symbol] organism = ",org,"\n")
    cat("[probe2symbol] probe.type = ",type,"\n")
    type

    if(type=="symbol") {
        cat("probe2symbol: probe is already symbol\n")
        if(any(grep(" /// ",probes))) {
            symbol0 <- strsplit(probes, split=" /// ")
        } else if(any(grep("[;,]",probes))) {
            symbol0 <- strsplit(probes, split="[;,\\|]")
        } else {
            symbol0 <- probes
        }
        ## all.symbols <- NULL
        ## if(org=="human") all.symbols <- unlist(as.list(org.Hs.egSYMBOL))
        ## if(org=="mouse") all.symbols <- unlist(as.list(org.Mm.egSYMBOL))
        ## symbol0 <- lapply(symbol0, function(s) intersect(s,all.symbols))

    } else {
        org
        if(org=="human") {
            symbol0 <- AnnotationDbi::mapIds(org.Hs.eg.db, probes, 'SYMBOL', toupper(type))
        }
        if(org=="mouse") {
            symbol0 <- AnnotationDbi::mapIds(org.Mm.eg.db, probes, 'SYMBOL', toupper(type))
        }
        if(org=="rat") {
          symbol0 <- AnnotationDbi::mapIds(org.Rn.eg.db, probes, 'SYMBOL', toupper(type))
        }
    }

    ## Unrecognize probes
    nna <- which(is.na(names(symbol0)))
    length(nna)
    if(length(nna)) names(symbol0)[nna] <- probes[nna]

    ## What to do with unmapped/missing symbols????
    symbol <- sapply(symbol0,"[",1)  ## takes first symbol only!!!
    isnull <- which(sapply(symbol,is.null))
    symbol[isnull] <- NA
    if(keep.na) {
        sel.na <- which(is.na(symbol))
        symbol[sel.na] <- probes[sel.na]
    }
    symbol <- unlist(symbol)
    names(symbol) <- NULL
    Matrix::head(symbol)

    symbol
}
