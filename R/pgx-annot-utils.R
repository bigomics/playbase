##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

MAIN_ORGANISMS = c("Human","Mouse","Rat","Homo sapiens",
  "Rattus norvegicus","Mus musculus","hsapiens","mmusculus","rnorvegicus")


#' Show all aliases
#' 
getSpeciesAliases <- function(species) {
  S <- playbase::SPECIES_TABLE
  matches <- apply(S, 2, function(g) which(tolower(g) %in% tolower(species)))
  matches <- unique(unlist(matches))
  kk <- c("species_name","display_name","ah_species","gprofiler_species",
    "gprofiler_id")
  kk <- intersect(kk, colnames(S))
  unique(as.character(unlist(S[matches, kk])))
}

#' Merges any missing annotation in df with non-missing annotation of
#' df2.
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
    probe_type <- detect_probetype.ANNOTHUB(organism = "custom", probes, orgdb = orgdb)
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


#' Rename features names of object to available human symbol by
#' ortholog or other 'human-like' uppercased annotation
#' columns. WARNING: does not necessarily keep original length.
#'
#' @export
collapse_by_humansymbol <- function(obj, annot) {
  annot <- cbind(annot, rownames = rownames(annot))
  target <- c("ortholog", "symbol", "gene_name", "rownames")
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

#' Two-pass implemenation of AnnotationDbi::select that checks for
#' original probe names and then matches clean names.
#' 
AnnotationDbi_select_2pass <- function(orgdb, keys, columns, keytype,
                                       verbose=1) {

  ## First pass: 
  annot <- try(AnnotationDbi::select(
    orgdb,
    keys = keys,
    columns = columns,
    keytype = keytype
  ), silent = TRUE)

  dim(annot)
  ## determine symbol column
  symbol.col <- intersect(colnames(annot),c("SYMBOL","GENENAME","ALIAS"))[1]

  if(inherits(annot, "try-error")) {
    annot <- NULL
    which.missing <- 1:length(keys)
  } else {
    which.missing <- which(is.na(annot[[symbol.col]]))
  }
  
  ## Second pass
  if(length(which.missing)) {
    clean.keys <- .clean_probe_names(keys[which.missing])
    names(clean.keys) <- keys[which.missing]
    annot2 <- try(AnnotationDbi::select(
      orgdb,
      keys = clean.keys,
      columns = columns,
      keytype = keytype
    ), silent = TRUE)
    if(!inherits(annot2, "try-error")) {
      jj <- which(!is.na(annot2[[symbol.col]]))
      if(length(jj)) {
        pp <- match( annot2[,keytype], clean.keys )
        annot2[,keytype] <- names(clean.keys)[pp]
        annot <- rbind(annot, annot2[jj,])
      }
    }
  }
  dim(annot)
  
  ## collapse to original keys (ordered). There may be duplicates
  ## from 2-pass matching. Prefer non-NA entries
  symbol <- annot[[symbol.col]]
  annot <- annot[order(symbol, na.last=TRUE),]
  annot <- annot[match(keys, annot[,keytype]),,drop=FALSE]
  
  return(annot)
}

#' Maps a gmt list to symbol using annotation table. We go via
#' sparsematrix because it is much faster than list filtering.
#' 
gmt.map2symbol <- function(gmt, annot, target="symbol") {
  G1 <- gmt2mat(gmt)
  as.symbol <- map2symbol(annot=annot, genes=rownames(G1), target=target, na.rm=FALSE)
  jj <- which(!is.na(as.symbol))
  G1 <- G1[jj,]
  rownames(G1) <- as.symbol[jj]
  mat2gmt(G1)
}


#' Merge duplicated GO sets in a gmt collection by merging the
#' terms. This is often needed after merging to GO collections if
#' retrieved by different methods.
#' 
go.merge_duplicates <- function(gmt) {
  gmt.id <- gsub(".*\\(GO_|\\)$","",names(gmt))
  gmt.names <- names(gmt)
  names(gmt.names) <- gmt.id
  ndup <- sum(duplicated(gmt.id))
  message(paste("merging",ndup,"duplicated GO terms"))
  if(ndup == 0) return(gmt)
  id.dup <- gmt.id[which(duplicated(gmt.id))]
  id.one <- setdiff(gmt.id, id.dup)
  ## colllapse duplicates by set union  
  gmt1 <- gmt[which(gmt.id %in% id.one)]
  jj <- which(gmt.id %in% id.dup)
  gmt2 <- tapply(gmt[jj], gmt.id[jj], function(g) unique(unlist(g)))
  names(gmt2) <- gmt.names[names(gmt2)]
  gmt <- c(gmt1, gmt2)
  return(gmt)
}

#' @title Detect probe type from probe set
#' @export
detect_probetype <- function(organism, probes, datatype = NULL, 
                             nprobe = 1000, verbose = FALSE) {

  organism <- normalizeOrganism(organism)

  if (is.null(datatype) && all(grepl("[:]", probes))) {
    dbg("[detect_probetype] datatype is multi-omics?")
    datatype <- "multi-omics"
  }

  if (!is.null(datatype) && datatype %in% c("metabolomics","lipidomics")) {
    probe_type <- mx.detect_probetype(probes)
    return(probe_type)
  }

  if (!is.null(datatype) && datatype == "multi-omics") {
    mx.probes <- sub("^mx:", "", grep("^mx:", probes, value = TRUE))
    px.probes <- sub("^px:", "", grep("^px:", probes, value = TRUE))
    gx.probes <- sub("^gx:", "", grep("^gx:", probes, value = TRUE))
    gx.probe_types <- px.probe_types <- mx.probe_types <- NA
    if (length(gx.probes)) gx.probe_types <- detect_probetype(organism, gx.probes)
    if (length(px.probes)) px.probe_types <- detect_probetype(organism, px.probes)
    if (length(mx.probes)) mx.probe_types <- mx.detect_probetype(mx.probes)
    probe_type <- c(gx = gx.probe_types, px = px.probe_types, mx = mx.probe_types)
    dtypes <- sort(unique(sub(":.*", "", probes)))
    probe_type <- probe_type[dtypes]
    return(probe_type)
  }

  ptype <- detect_probetype.ANNOTHUB(
    organism = organism, probes = probes, orgdb = NULL,
    nprobe = nprobe, use.ah = NULL, datatype = datatype,
    verbose = verbose) 

  if(is.null(ptype) || is.na(ptype)) {
    ptype <- detect_probetype.GPROFILER(
      organism = organism, probes = probes, nprobe = nprobe,
      datatype = datatype, verbose = verbose) 
  }

  if(is.null(ptype) || is.na(ptype)) {
    return(NA)  ## expects NA for fail
  }
  return(ptype)
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

#' Detect/validate features with gprofiler
#'
detect_probetype.GPROFILER <- function(organism, probes, nprobe = 1000,
                                       datatype = NULL, verbose = TRUE) {
  gp.organism <- .map_gprofiler_id(organism)    

  if (length(probes) > nprobe) {
    probes <- sample(probes, nprobe)
  }
  probesx <- gsub(";.*|_.*|.*:","",probes)
  probes <- c(probes, probesx)
  gp.out <- tryCatch(
  {
    gprofiler2::gconvert(probes, organism = gp.organism, target = "UNIPROT_GN_ACC")
  },
  error = function(e) {
    return(NULL)
  }
  )
  return("GPROFILER2")
}


#' @title Get all species in AnnotationHub/OrgDB
#'
#' @export
allSpecies <- function(col = "species_name") {
  M <- data.frame(playbase::SPECIES_TABLE)
  col <- intersect(col, colnames(M))[1]
  if(length(col)==0) return(NULL)
  species <- as.character(M[, col])
  names(species) <- M[, "taxonomyid"]
  species
}

.getSpeciesTable.GPROFILER <- function() {
  jsonlite::fromJSON("https://biit.cs.ut.ee/gprofiler/api/util/organisms_list")
}

#' @title Get species table in AnnotationHub/OrgDB
#'
#' @export
.getSpeciesTable.ANNOTHUB <- function(ah = NULL) {
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
#' @export
check_species_probetype <- function(
  probes,
  test_species = c("Human", "Mouse", "Rat"),
  datatype = NULL, annot.cols = NULL
) {
  ## No check if custom
  custom_datatype <- !is.null(datatype) && tolower(datatype) %in% c("custom", "unknown", "")
  custom_organism <- any(tolower(test_species) %in% c("custom", "unknown", "no organism"))

  if (custom_datatype || custom_organism) {
    out <- rep("custom", length(test_species))
    names(out) <- test_species
    return(as.list(out))
  }

  probes <- unique(.clean_probe_names(probes))
  ## report possible probetype per organism
  ptype <- vector("list", length(test_species))
  names(ptype) <- test_species
  if (datatype == "metabolomics") {
    mx.type <- NA
    if (!is.null(annot.cols)) {
      mx.ids <- toupper(colnames(playdata::METABOLITE_ID)[-1])
      mx.ids <- c(mx.ids, paste0(mx.ids, "_ID"))
      has.id <- any(toupper(annot.cols) %in% mx.ids)
      if (has.id) {
        ids <- intersect(toupper(annot.cols), mx.ids)
        mx.type <- ids[1]
      }
    }
    if (all(is.na(mx.type))) {
      db <- mx.check_mapping(probes, check.first = TRUE)
      table(db)
      if (!all(is.na(db))) {
        mx.type <- names(which.max(table(db[!is.na(db)])))
      }
    }
    for (s in test_species) ptype[[s]] <- mx.type
  } else {
    s <- "Human"
    for (s in test_species) {
      ptype[[s]] <- detect_probetype(
        organism = s,
        probes = probes,
        datatype = datatype,
        verbose = FALSE
      )
    }
  }

  ## remove NA
  ptype <- ptype[!sapply(ptype, function(p) all(is.na(p)))]
  return(ptype)
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
                              datatype = NULL, orgdb = NULL, verbose = TRUE) {
  organism <- normalizeOrganism(organism)

  if (!is.null(datatype) && datatype == "metabolomics") {
    new.probes <- mx.convert_probe(probes, target_id = target_id)
    return(new.probes)
  }

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


#' Sort/dedupe a set of candidate ortholog symbols, dropping empty and
#' N/A placeholders. Returns a sorted character vector (possibly empty).
#'
#' @noRd
.clean_ortholog_set <- function(x) {
  x <- unlist(strsplit(as.character(x), ";"))
  sort(unique(setdiff(x, c("", NA, "NA", "N/A"))))
}
