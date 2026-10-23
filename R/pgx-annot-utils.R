##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

MAIN_ORGANISMS = c("Human","Mouse","Rat","Homo sapiens",
  "Rattus norvegicus","Mus musculus","hsapiens","mmusculus","rnorvegicus")


#' Normalize organism name to standard format
#'
#' @description Converts common organism name variants to their standard
#' scientific names. This ensures consistent organism naming across all
#' annotation functions.
#'
#' @param organism Character string with organism name (e.g., "human", "mouse", "dog")
#' @return Normalized organism name in standard format (e.g., "Homo sapiens")
#'
#' @examples
#' normalizeOrganism("human")
#' # Returns: "Homo sapiens"
#'
#' normalizeOrganism("dog")
#' # Returns: "Canis familiaris"
#'
#' @export
normalizeOrganism <- function(organism) {
  if (is.null(organism) || is.na(organism)) {
    return(organism)
  }
  org_lower <- tolower(organism)
  if (org_lower == "human") {
    return("Homo sapiens")
  }
  if (org_lower == "mouse") {
    return("Mus musculus")
  }
  if (org_lower == "rat") {
    return("Rattus norvegicus")
  }
  if (grepl("canis.*familiaris|^dog$", org_lower)) {
    return("Canis familiaris")
  }
  organism
}

#' @export
gene2uniprot <- function(genes, organism) {
  gp.organism <- .map_gprofiler_id(organism)
  out <- try(gprofiler2::gconvert(genes, organism = gp.organism, target = "UNIPROT_GN_ACC"))
  if (is.null(out) || "try-error" %in% class(out)) {
    return(NULL)
  }
  res <- tapply(out$target, out$input, function(s) {
    paste(setdiff(unique(s), c(NA, "")), collapse = ";")
  })
  ii <- match(genes, names(res))
  ## seems input is uppercase!
  ii <- ifelse(is.na(ii), match(toupper(genes), toupper(names(res))), ii)
  res[ii]
}

#' @export
uniprot2gene <- function(uniprots, organism) {
  gp.organism <- .map_gprofiler_id(organism)  
  out <- try(gprofiler2::gconvert(uniprots, organism = gp.organism, target = "ENSG"))
  if (is.null(out) || "try-error" %in% class(out)) {
    return(NULL)
  }
  res <- tapply(out$name, out$input, function(s) paste(setdiff(unique(s), c(NA, "")), collapse = ";"))
  res[uniprots]
}


#' non-greedy removal of numerical postfix. Postfix is defined as
#' numerical substring after . (dot)
#'
strip_postfix <- function(s) {
  stripFUN <- function(s) {
    #sub(paste0("[._].*$|[-][0-9.]+$"), "", s)
    sub(paste0("[.][0-9]+$"), "", s) 
  }
  ss <- strsplit(s, split = ";")
  ss <- lapply(ss, function(s) stripFUN(s))
  sapply(ss, paste, collapse = ";")
}


#' Cleanup probe names from postfixes or version numbers. Retains
#' prefix needed for multi-omics.
#'
.clean_probe_names <- function(probes, sep = ".-") {
  probes0 <- probes
  probes <- trimws(probes)  
  probes[is.na(probes)] <- ""
  ## strip multiple probes
  probes <- sub("[;].*", "", probes)
  ## strip away anything postfix after a 'dot' or 'underscore'
  probes <- sub(paste0("[", sep, "].*"), "", probes)
  names(probes) <- probes0
  return(probes)
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
probe2symbol <- function(probes, annot_table, query = "symbol",
                         key = NULL, fill_na = FALSE, add_datatype = FALSE) {

  # NULL annot_table: no mapping possible — return probes as-is (same behaviour as fill_na=TRUE)
  if (is.null(annot_table)) return(probes)

  # Prepare inputs. add extra matching columns.
  annot_table <- cbind(rownames = rownames(annot_table), annot_table)
  id.cols <- intersect(c("feature", "gene_name", "symbol"), colnames(annot_table))
  if (length(id.cols) > 0) {
    stripped_annot <- apply(annot_table[, id.cols, drop = FALSE], 2, function(a) sub("^[A-Za-z]+:", "", a))
    ## colnames(stripped_annot) <- paste0(colnames(stripped_annot),"_stripped")
    annot_table <- cbind(annot_table, stripped_annot)
  }

  probes1 <- setdiff(probes, c(NA, ""))
  if (is.null(key) || !key %in% colnames(annot_table)) {
    key <- which.max(apply(annot_table, 2, function(a) sum(probes1 %in% a)))
  }
  if (is.null(key)) {
    message("[probe2symbol] FATAL. could not get key column.")
    return(NULL)
  }

  query <- head(intersect(query, colnames(annot_table)), 1)
  if (length(query) == 0) {
    message("ERROR. no symbol column.")
    return(NULL)
  }

  # fall back on old gene_name
  if (query == "symbol" && !"symbol" %in% colnames(annot_table) &&
    "gene_name" %in% colnames(annot_table)) {
    query <- "gene_name"
  }

  # match query
  ii <- match(probes, annot_table[, key])
  query_col <- annot_table[ii, query]

  # Deal with NA
  if (fill_na) {
    query_col <- data.table::fifelse(query_col == "" | is.na(query_col),
      yes = probes,
      no = query_col
    )
  }

  # Prepend datatype if requested and available
  if (add_datatype && "data_type" %in% colnames(annot_table)) {
    datatype_col <- annot_table[ii, "data_type"]
    has_datatype <- !is.na(datatype_col) & datatype_col != ""
    # Check if query_col already has the datatype prefix
    already_has_prefix <- startsWith(query_col, paste0(datatype_col, ":"))
    should_add <- has_datatype & !already_has_prefix
    query_col <- ifelse(should_add, paste0(datatype_col, ":", query_col), query_col)
  }

  # Return queryed col
  return(query_col)
}

#'
#'
#' @export
getOrgDb <- function(organism, use.ah = NULL) {
  organism <- normalizeOrganism(organism)
  orgdb <- .getOrgDb(organism, use.ah = use.ah)
  if (is.null(orgdb)) {
    message("[getOrgDb] ERROR: could not get orgdb")
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


## not exported
.getOrgDb <- function(organism, use.ah = NULL) {
  organism <- normalizeOrganism(organism)

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
    if (organism == "Plasmodium falciparum" && requireNamespace("org.Pf.plasmo.db", quietly = TRUE)) {
      return(org.Pf.plasmo.db::org.Pf.plasmo.db)
    }
  }

  ah <- AnnotationHub::AnnotationHub()
  #  all_species <- allSpecies()
  #  if (!tolower(organism) %in% tolower(all_species)) {
  #    message("WARNING: organism '", organism, "' not in AnnotationHub")
  #    return(NULL)
  #  }

  message("querying AnnotationHub for '", organism, "'\n")
  suppressMessages({
    ahDb <- try(AnnotationHub::query(ah, pattern = c(organism, "OrgDb")))
  })

  if (length(ahDb) == 0 || inherits(ahDb, "try-error")) {
    message("WARNING: organism '", organism, "' not in AnnotationHub.")
    return(NULL)
  }

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


#' Rename features names of object to available human symbol by
#' ortholog or other 'human-like' uppercased annotation
#' columns. WARNING: does not necessarily keep original length.
#'
#' @export
collapse_by_humansymbol <- function(obj, annot) {
  annot <- cbind(annot, rownames = rownames(annot))
  target <- c("human_ortholog", "ortholog", "symbol", "gene_name", "rownames")
  target <- intersect(target, colnames(annot))
  if (length(target) == 0) {
    message("[collapse_by_humansymbol] WARNING: could not find symbol mapping column.")
    return(obj)
  } else {
    ## call rename_by with target column
    k <- target[1]
    sel.na <- which(annot[, k] %in% c(NA, "", "-", "---", "NA"))
    annot[sel.na, k] <- "---"
    annot[, k] <- toupper(annot[, k]) ## all uppercase??
    map.obj <- rename_by(obj, annot_table = annot, new_id = target[1])
  }
  if (!is.null(dim(map.obj))) rownames(map.obj) <- toupper(rownames(map.obj))
  if (is.null(dim(map.obj))) names(map.obj) <- toupper(names(map.obj))
  map.obj
}

#'
#' 
.map_gprofiler_id <- function(species) {

  if (is.null(species) || is.na(species) || !nzchar(species)) return(NULL)

  orgs <- jsonlite::fromJSON("https://biit.cs.ut.ee/gprofiler/api/util/organisms_list")

  ## exact match
  exact.species <- paste0("^",species,"$")
  i <- which(
    grepl( exact.species, orgs$id, ignore.case = TRUE) |
      grepl(exact.species, orgs$scientific_name, ignore.case = TRUE) |
      grepl( exact.species, orgs$display_name, ignore.case = TRUE)
  )

  ## internal match
  if(length(i)==0) {
    i <- which(grepl( species, orgs$scientific_name, ignore.case = TRUE) |
                 grepl( species, orgs$display_name, ignore.case = TRUE))
  }

  if(length(i)==0) return(NULL)

  ## Several strains/assemblies of the same species can match (e.g.
  ## "Caenorhabditis elegans" matches both the core "celegans" organism
  ## and WormBase ParaSite imports like "caelegprjna13758"). g:Profiler's
  ## core reference organisms use short, digit-free ids; strain/assembly
  ## imports carry a numeric accession in their id. Prefer the former
  ## when several exact matches tie, since that's the one gorth()/GO
  ## lookups actually support.
  if (length(i) > 1) {
    core <- i[!grepl("[0-9]", orgs[i,"id"])]
    if (length(core) > 0) i <- core
  }

  id <- orgs[i[1],"id"]
  return(id)
}

#' Create new feature name by concatenating some columns of input
#' annotation table. Make all feature names unique.
#'
#' @param annot  some annotation dataframe
#' @param target vector of character. e.g. c("feature","_","symbol")
#'
#' export
combine_feature_names <- function(annot, target) {
  annot$rownames <- rownames(annot)
  new.feature <- strsplit(annot[, target[1]], split = ";")
  for (i in 2:length(target)) {
    if (target[i] %in% colnames(annot)) {
      new.feature <- mapply(paste0, new.feature, annot[, target[i]])
    } else {
      ## some character
      new.feature <- mapply(paste0, new.feature, target[i])
    }
  }
  new.feature <- sapply(new.feature, paste0, collapse = ";")
  if (sum(duplicated(new.feature)) > 0) {
    new.feature <- make_unique(new.feature)
  }
  new.feature
}

#' @title Detect probe type from probe set
#' @description Forwards to `playbase.annotation::detect_probetype()`.
#' @export
detect_probetype <- function(organism, probes, datatype = NULL, 
                             nprobe = 1000, verbose = FALSE) {
  playbase.annotation::detect_probetype(
    organism = organism,
    probes = probes,
    datatype = datatype,
    nprobe = nprobe,
    verbose = verbose
  )
}

#' @title Detect probe type from probe set
#' @export
detect_probetype.ANNOTHUB <- function(organism, probes, orgdb = NULL,
                                      nprobe = 1000, use.ah = NULL, datatype = NULL,
                                      verbose = TRUE) {

  ## get correct OrgDb database for organism
  if (is.null(orgdb)) {
    orgdb <- getOrgDb(organism, use.ah = use.ah)
  }
  if (is.null(orgdb)) {
    if (verbose) message("[detect_probetype] ERROR: unsupported organism '", organism, "'\n")
    return(NULL)
  }

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

  ## try different cleaning methods. NEED RETHINK!!!! refseq has
  ## underscore!
  probes0 <- probes
  probes1 <- .clean_probe_names(probes)
  probesx <- unique(c(probes0, probes1))

  ## get probe types for organism
  keytypes <- c(
    "SYMBOL", "ENSEMBL", "ACCNUM", "UNIPROT", "GENENAME",
    "ALIAS", "MGI", "TAIR", ## organism specific
    "ENSEMBLTRANS", "ENSEMBLPROT",
    "REFSEQ", "ENTREZID"
  )
  keytypes <- intersect(keytypes, AnnotationDbi::keytypes(orgdb))
  key_matches <- rep(0L, length(keytypes))
  names(key_matches) <- keytypes

  ## Get all organism symbols
  org_annot <- AnnotationDbi::select(
    orgdb,
    keys = AnnotationDbi::keys(orgdb, "ENTREZID"),
    keytype = "ENTREZID",
    columns = intersect(c("SYMBOL", "GENENAME"), keytypes)
  )
  org_symbols <- NULL
  org_genenames <- NULL
  if ("SYMBOL" %in% colnames(org_annot)) org_symbols <- setdiff(org_annot[, "SYMBOL"], c("", NA))
  if ("GENENAME" %in% colnames(org_annot)) org_genenames <- setdiff(org_annot[, "GENENAME"], c("", NA))

  # Iterate over probe types
  key <- keytypes[1]
  for (key in keytypes) {
    probe_matches <- data.frame(NULL)
    # add symbol and genename on top of key as they will be used to
    # count the real number of probe matches
    key2 <- intersect(c(key, "SYMBOL", "GENENAME"), keytypes)
    suppressMessages(suppressWarnings(try(
      probe_matches <- AnnotationDbi::select(
        orgdb,
        keys = probesx,
        keytype = key,
        columns = key2
      ),
      silent = TRUE
    )))

    if (nrow(probe_matches) && ncol(probe_matches)) {
      ## extra check: if key is SYMBOL or GENENAME first column can be
      ## wrongly set as the key.
      if ("SYMBOL" %in% colnames(probe_matches) && !is.null(org_symbols)) {
        not.symbol <- !(probe_matches[, "SYMBOL"] %in% org_symbols)
        probe_matches[, "SYMBOL"][not.symbol] <- NA
      }
      if ("GENENAME" %in% colnames(probe_matches) && !is.null(org_genenames)) {
        not.gene <- !(probe_matches[, "GENENAME"] %in% org_genenames)
        probe_matches[, "GENENAME"][not.gene] <- NA
      }

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
  }
  key_matches <- round(key_matches, 4)
  key_matches

  ## Return top match key_matches
  top_match <- NULL
  if (all(key_matches == 0)) {
    message("WARNING: Probe type not found. Valid probe types: ",
      paste(keytypes, collapse = " "))
    return(NULL)
  }
  if (max(key_matches, na.rm = TRUE) < 0.01) {
    message("WARNING: Insufficient matching ratio. Max match = ",
      max(key_matches, na.rm = TRUE))
    return(NULL)
  }
  if (max(key_matches, na.rm = TRUE) < 0.50) {
    message("WARNING: Low matching ratio. Max match = ", max(key_matches, na.rm = TRUE))
  }
  top_match <- names(which.max(key_matches))
  return(top_match)
}

#' Check if probes can be detected by Orthogene or AnnotHub/OrgDb
#' annotation engines.
#'
#' export
check_probetype <- function(organism, probes, verbose=1) {
  chk1 <- .check_probetype.GPROFILER(organism, probes, min.map=0.20)  
  if (!is.null(chk1) && chk1 == TRUE) {
    if(verbose) message("organism/features supported by Gprofiler")
    return(TRUE)
  }
  ## using AnnotHub/OrgDb
  chk2 <- detect_probetype(organism, probes)
  if (!is.null(chk2)) {
    if(verbose) message("organism/features supported by AnnotHub")    
    return(TRUE)
  }
  if(verbose) message("Warning: organism/features not recognized")    
  return(FALSE)
}

.check_probetype.GPROFILER <- function(organism, probes, min.map=0.20) {
  gp.organism <- .map_gprofiler_id(organism)
  map <- try(gprofiler2::gconvert(probes, organism = gp.organism, target = "ENSG"))
  if ("try-error" %in% class(map) || is.null(map)) {
    message("[check_probetype.GPROFILER] *WARNING* organism not  recogized, or server not reachable")
    return(NULL)
  }
  mean.mapped <- mean(!is.na(map$target))
  ## get correct OrgDb database for organism
  if (mean.mapped < min.map) {
    message("[check_probetype.GPROFILER] *WARNING* too low mapping coverage")
    return(FALSE)
  }
  return(TRUE)
}


#' Automatically detects species by trying to detect probetype from
#' list of test_species. Warning. bit slow.
#'
#' Forwards to `playbase.annotation::check_species_probetype()`.
#'
#' @export
check_species_probetype <- function(
  probes,
  test_species = c("Human", "Mouse", "Rat"),
  datatype = NULL, annot.cols = NULL
) {
  playbase.annotation::check_species_probetype(
    probes = probes,
    test_species = test_species,
    datatype = datatype,
    annot.cols = annot.cols
  )
}

#' Annotate phosphosite with residue symbol. Feature names must be of
#' form 'uniprot_position'. NOTE!!! Annotation is currently done here
#' in feature name but it would be 'better' to add phosphosite
#' modification type in the pgx$genes general annotation table.
#'
#' @export
annotate_phospho_residue <- function(features, detect.only = FALSE) {
  valid_name <- mean(grepl("[_][A-Z]?[0-9]+", features), na.rm = TRUE) > 0.9
  valid_name
  uniprot <- sub("[_].*", "", features)
  positions <- gsub(".*[_][A-Za-z]?|[.].*", "", features)
  positions <- strsplit(positions, split = "[;/,]")

  P <- playdata::PHOSPHOSITE
  prot.match <- mean(uniprot %in% P$UniProt, na.rm = TRUE)
  pos.match <- mean(positions %in% P$Position, na.rm = TRUE)
  is_phospho <- (valid_name && prot.match > 0.50 && pos.match > 0.50)
  is_phospho

  if (detect.only) {
    return(is_phospho)
  }

  if (is_phospho) {
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

    ## determine separators for paste: sep1 for main position
    ## separator. sep2 for entries with multiple positions.
    sep1.match <- sapply(c("_", "."), function(s) {
      sum(grepl(s, features, fixed = TRUE), na.rm = TRUE)
    })
    sep1 <- names(which.max(sep1.match))
    sel <- grep("[;/,]", features)
    sep2.match <- sapply(c(";", "/", ","), function(s) {
      sum(grepl(s, features[sel], fixed = TRUE), na.rm = TRUE)
    })
    sep2 <- names(which.max(sep2.match))

    ## insert modification type in front of position
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
                              orgdb = NULL, verbose = TRUE) {
  organism <- normalizeOrganism(organism)

  ## get correct OrgDb database for organism
  if (is.null(orgdb)) {
    orgdb <- getOrgDb(organism)
  }
  if (is.null(orgdb)) {
    if (verbose) message("[convert_probetype] ERROR: unsupported organism '", organism, "'\n")
    return(NULL)
  }

  if (!target_id %in% AnnotationDbi::keytypes(orgdb)) {
    message("[convert_probetype] invalid target probetype")
    return(NULL)
  }
  if (is.null(from_id)) {
    from_id <- detect_probetype.ANNOTHUB(organism, probes, orgdb = orgdb, datatype = NULL)
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


