##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##


#' Merges any missing annotation in df with non-missing annotation of
#' annot_table.
#'
merge_annot_table <- function(df, df2, priority = 1) {
  #  df2 <- df2[match(rownames(df), rownames(df2)), ]
  #  rownames(df2) <- rownames(df)
  if (nrow(df) != nrow(df2)) stop("df and df2 not same size")
  if (inherits(df2, "matrix")) {
    df2 <- as.data.frame(df2)
  }
  ## add columns by filling missing values in df or df2
  cols <- union(colnames(df), colnames(df2))
  for (k in setdiff(cols, colnames(df))) df[[k]] <- NA
  for (k in setdiff(cols, colnames(df2))) df2[[k]] <- NA

  ## merge common columns by filling missing values in df2
  for (k in cols) {
    a <- df[, k]
    b <- df2[, k]
    na.chars <- c(NA, "NA", "", "-", "---", "unknown")
    if (priority == 1) {
      replace.a <- (a %in% na.chars) & !(b %in% na.chars)
    } else {
      replace.a <- !(b %in% na.chars) ## always replace
    }
    if (any(replace.a)) df[, k] <- ifelse(replace.a, b, a)
  }

  return(df)
}


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


#' non-greedy removal of prefixes. Prefix is defined as any
#' character substring (no spaces, no special chars) before the
#' matching colon character :.
#'
strip_prefix <- function(s) {
  stripFUN <- function(s) {
    sub("^[a-zA-Z]+:", "", s)
  }
  ss <- strsplit(s, split = ";")
  ss <- lapply(ss, function(s) stripFUN(s))
  sapply(ss, paste, collapse = ";")
}


#' Cleanup symbols names from postfixes and prefixes. Take only first
#' symbol. This is mostly used for symbol lookup tables that need one
#' clean symbol.
#'
.clean_symbols <- function(symbols) {
  strip_prefix(strip_postfix(sub(";.*", "", trimws(symbols))))
}

#' Cleanup probe names from postfixes or version numbers. Retains
#' prefix needed for multi-omics.
#'
.clean_probe_names <- function(probes, sep = ".-") {
  probes0 <- trimws(probes)
  probes[is.na(probes)] <- ""
  ## strip multiple probes
  probes <- sub("[;].*", "", probes)
  ## strip away anything postfix after a 'dot' or 'underscore'
  probes <- sub(paste0("[", sep, "].*"), "", probes)
  return(probes)
}

#' Match dirty probe names to clean key names
#'
#' @export
match_probe_names <- function(probes, orgdb, probe_type = NULL) {
  if (is.character(orgdb)) orgdb <- getOrgDb(orgdb)
  if (is.null(orgdb)) {
    message("[match_probe_names] ERROR could not get orgdb!")
    return(NULL)
  }
  if (is.null(probe_type)) {
    probe_type <- detect_probetype(organism = "custom", probes, orgdb = orgdb)
  }
  ## bail out if not annothub keytypes
  if (!probe_type %in% AnnotationDbi::keytypes(orgdb)) {
    return(probes)
  }

  probe.names <- names(probes)
  all.keys <- AnnotationDbi::keys(orgdb, probe_type)
  tsub <- function(s) gsub("[-:;.]|\\[|\\]", ".", s)
  ii <- match(toupper(tsub(probes)), toupper(tsub(all.keys)))
  table(is.na(ii))
  new.probes <- all.keys[ii]
  if (sum(is.na(new.probes))) {
    jj <- which(is.na(new.probes))
    new.probes[jj] <- probes[jj]
    jj.probes <- .clean_probe_names(probes[jj], sep = ".-")
    ii <- match(toupper(tsub(jj.probes)), toupper(tsub(all.keys)))
    if (any(!is.na(ii))) {
      k <- which(!is.na(ii))
      new.probes[jj[k]] <- all.keys[ii[k]]
    }
  }
  names(new.probes) <- probe.names
  new.probes
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
    if (organism == "Plasmodium falciparum" && require("org.Pf.plasmo.db", quietly = TRUE)) {
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


#' Rename features names of object to available human symbol by
#' human_ortholog or other 'human-like' uppercased annotation
#' columns. WARNING: does not necessarily keep original length.
#'
#' @export
collapse_by_humansymbol <- function(obj, annot) {
  annot <- cbind(annot, rownames = rownames(annot))
  target <- c("human_ortholog", "symbol", "gene_name", "rownames")
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

