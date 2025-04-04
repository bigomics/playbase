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
pgx.addGeneAnnotation <- function(pgx, annot_table = NULL) {
  # Safety checks
  stopifnot(is.list(pgx))

  dbg("[pgx.addGeneAnnotation] *** DEPRECATED ***")
  
  probes <- rownames(pgx$counts)
  datatype <- pgx$datatype
  organism <- pgx$organism
  probe_type <- pgx$probe_type
  
  genes <- getProbeAnnotation(
    organism,
    probes,
    datatype,
    probetype = probe_type,    
    annot_table = annot_table
  )

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

merge_annot_table <- function(df, annot_table) {
  kk <- unique(c(colnames(df),colnames(annot_table)))
  annot_table <- annot_table[match(rownames(df), rownames(annot_table)), ]
  rownames(annot_table) <- rownames(df)
  
  ## merge common columns by filling missing values in annot_table
  common.cols <- intersect(colnames(df),colnames(annot_table))
  for(k in common.cols) {
    a <- annot_table[,k]
    b <- df[,k]
    annot_table[,k] <- ifelse( a %in% c(NA,"NA","","-","---"), b, a)
  }
  
  ## merge data.frame, skip common columns
  df <- df[, setdiff(colnames(df), colnames(annot_table)), drop=FALSE]
  cbind(df, annot_table)[,kk]
}


#' Convert multi-omics probetype. Probe names *must  be prefixed with
#' data type unless classical transcriptomics/proteomics.
#'
#' @export
getProbeAnnotation <- function(organism,
                               probes,
                               datatype,
                               probetype = "",
                               annot_table = NULL
                               ) {
  
  if(is.null(datatype)) datatype <- "unknown"  
  if(is.null(probetype)) probetype <- "unknown"

  unknown.organism <- (tolower(organism) %in% c("no organism","custom","unkown"))
  unknown.datatype <- (datatype %in% c("custom","unkown"))
  unknown.probetype <- (probetype %in% c("custom","unkown"))  
  annot.unknown <- unknown.organism || unknown.datatype || unknown.probetype
  annot.unknown

  ## clean probe names
  probes <- trimws(probes)
  probes[probes=="" | is.na(probes)] <- 'NA'
  probes0 <- make_unique(probes)  ## make unique but do not clean
#  probes <- make_unique(clean_probe_names(probes0))  ## NEED RETHINK!! really???
  if(!is.null(annot_table)) {
    rownames(annot_table) <- make_unique(rownames(annot_table))
  }
  
  genes <- NULL
  if (annot.unknown) {
    # annotation table is mandatory for 'No organism' (until server side
    # can handle missing genesets)
    info("[getProbeAnnotation] annotating with custom annotation")
    genes <- getCustomAnnotation2( probes0, annot_table )
  } else if (datatype == "metabolomics") {
    dbg("[getProbeAnnotation] annotating for metabolomics")
    mx.check <- mx.check_mapping(
      probes, all.db=c("playdata","annothub","refmet"), check.first = TRUE)
    mx.check <- mean(!is.na(mx.check)) > 0.01
    mx.check
    if(mx.check) {
      ## Directly annotate if probes are recognized
      ##genes <- getMetaboliteAnnotation(probes, add_id=TRUE, annot_table=annot_table)
      genes <- getMetaboliteAnnotation(probes, add_id=TRUE, annot_table=NULL)      
    } else {
      ## Fallback on custom
      dbg("[getProbeAnnotation] WARNING: not able to map metabolomics probes")
    }
  } else if (datatype == "multi-omics") {
    dbg("[getProbeAnnotation] annotating for multi-omics")
    genes <- getMultiOmicsProbeAnnotation(organism, probes)
  } else {
    dbg("[getProbeAnnotation] annotating for transcriptomics")
    genes <- getGeneAnnotation( organism = organism, probes = probes )
  }
  
  ## final fallback is genes==NULL
  if(is.null(genes)) {
    dbg("[getProbeAnnotation] WARNING: fallback to UNKNOWN probes")
    genes <- getCustomAnnotation( probes0, custom_annot = NULL )
  }

  ## if annot_table is provided we override our annotation and append
  ## extra columns
  if (!is.null(genes) && !is.null(annot_table)) {
    dbg("[getProbeAnnotation] merging custom annotation table")
    genes <- merge_annot_table(genes, annot_table) 
  }

  ## restore original probe names
  rownames(genes) <- genes$feature <- probes0

  ## cleanup entries and reorder columns
  genes <- cleanupAnnotation(genes)

  return(genes)

}

#' Get gene annotation data using annothub or orthogene.
#'
#' @export
getGeneAnnotation <- function(
    organism,
    probes,
    use.ah = NULL,
    verbose = TRUE,
    methods = c("annothub", "gprofiler")
    ) {
  annot <- NULL

  if (tolower(organism) == "human") organism <- "Homo sapiens"
  if (tolower(organism) == "mouse") organism <- "Mus musculus"
  if (tolower(organism) == "rat") organism <- "Rattus norvegicus"
  if (tolower(organism) == "dog") organism <- "Canis familiaris"

  probes <- trimws(probes)
  probes[probes=="" | is.na(probes)] <- 'NA'
  
  if(mean(grepl("[:]",probes)) > 0.98) {
    message("[getGeneAnnotation] WARNING. stripping multi-omics prefix")
    probes <- sub("^[a-zA-Z]+:", "", probes)
  }

  # init empty (all missings)
  annot <- data.frame(feature=probes, stringsAsFactors=FALSE)
  missing <- rep(TRUE, length(probes))

  for (method in methods) {
    if (any(missing)) {
      # annotation for current method
      curr_annot <- try(switch(
        method,
        "annothub" = getGeneAnnotation.ANNOTHUB(
          organism = organism,
          probes = probes,
          use.ah = use.ah,
          verbose = verbose
        ),
        "gprofiler" = getGeneAnnotation.ORTHOGENE(
          organism = organism,
          probes = probes,
          verbose = verbose
        ),
        stop("Unknown method: ", method)
      ))
      if (inherits(curr_annot, "try-error")) next

      if (!is.null(curr_annot) && nrow(curr_annot) > 0) {
        # not all methods have the same columns
        new_cols <- setdiff(colnames(curr_annot), colnames(annot))
        if (length(new_cols) > 0) {
          for (col in new_cols) {
            annot[, col] <- rep(NA, nrow(annot))
          }
        }
        annot[missing, colnames(curr_annot)] <- curr_annot
        missing <- is.na(annot$symbol) | annot$symbol == ""
      }
    }
  }

  if (all(missing)) { # unsuccessful annotation
    annot <- NULL
  }
  
  ## clean up
  if(!is.null(annot)) {
    annot <- cleanupAnnotation(annot)
  }

  return(annot)
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
    second.pass = TRUE,
    verbose = TRUE) {
  if (is.null(organism)) {
    warning("[getGeneAnnotation.ANNOTHUB] Please specify organism")
    return(NULL)
  }

  if (verbose) {
    message("[getGeneAnnotation.ANNOTHUB] Retrieving gene annotation...")
  }

  if (tolower(organism) == "human") organism <- "Homo sapiens"
  if (tolower(organism) == "mouse") organism <- "Mus musculus"
  if (tolower(organism) == "rat") organism <- "Rattus norvegicus"

  genes <- NULL

  ## get correct OrgDb database for this organism
  orgdb <- getOrgDb(organism, use.ah = use.ah)
  if(is.null(orgdb)) {
    message("[getGeneAnnotation.ANNOTHUB] ERROR: orgdb == NULL: ", is.null(orgdb) )
    return(NULL)
  }
  
  if (is.null(probes)) {
    probes <- AnnotationDbi::keys(orgdb)
  }
  ## Backup probes as probes0, give names of probes the original name.
  probes[is.na(probes) | probes==""] <- "NA"
  probes0 <- make.unique(probes)
  names(probes) <- probes0

  ## clean up probe names from suffixes
  probes <- clean_probe_names(probes)

  if (is.null(probe_type)) {
    probe_type <- detect_probetype(organism, probes, orgdb = orgdb)
    probe_type
    if (is.null(probe_type) || is.na(probe_type) ) {
      message("ERROR: could not determine probe_type.")
      message("WARNING!!! returning empty annotation.")
      annot <- data.frame(feature = probes, symbol = "")
      annot <- cleanupAnnotation(annot)
      return(annot)
    }
  }

  if (probe_type == "GPROFILER") {
    dbg("[getGeneAnnotation.ANNOTHUB] probe_type = GPROFILER; skipping annothub")
    return(NULL)
  }
  
  ## Match to clean probe names (???)
  probes <- match_probe_names(probes, orgdb, probe_type)

  ## --------------------------------------------
  ## retrieve table
  ## --------------------------------------------
  cols <- c("SYMBOL", "GENENAME", "GENETYPE", "MAP", "ALIAS", "UNIPROT")
  #cols <- c("SYMBOL", "GENENAME", "GENETYPE", "MAP")  
  cols <- intersect(cols, AnnotationDbi::keytypes(orgdb))

  if (organism %in% c("Mus musculus", "Rattus norvegicus")) {
    cols <- unique(c(cols, "ENTREZID"))
  }

  dbg("[getGeneAnnotation.ANNOTHUB] annotation columns:", cols)
  dbg("[getGeneAnnotation.ANNOTHUB] probe_type =  ", probe_type)
  dbg("[getGeneAnnotation.ANNOTHUB] retrieving annotation for ",length(probes),"features.")  
  suppressMessages(suppressWarnings(
    annot <- AnnotationDbi::select(
      orgdb,
      keys = probes,
      columns = cols,
      keytype = probe_type
    )
  ))
  symbols <- AnnotationDbi::keys(orgdb, keytype="SYMBOL")
  
  # some organisms do not provide symbol but rather gene name (e.g. yeast)
  if (!"SYMBOL" %in% colnames(annot)) {
    annot$SYMBOL <- annot$GENENAME
    annot$GENENAME <- annot$ALIAS
  }
  if ("SYMBOL" %in% colnames(annot)) {
    not.symbols <- !(annot$SYMBOL %in% symbols)
    if(length(not.symbols)) annot$SYMBOL[not.symbols] <- NA
  }
  annot$ALIAS <- NULL
  annot$SYMBOL[is.na(annot$SYMBOL)] <- ""
  
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
  
  ## match annotation table to probes
  info("got", length(unique(annot$SYMBOL)), "unique SYMBOLs")
  key <- annot[, probe_type]
  dfA <- apply(annot,2,function(a) tapply(a,key,function(b)
    paste(setdiff(unique(b),c(NA,"")),collapse=';')))
  annot <- data.frame(dfA, check.names=FALSE)
  annot <- annot[match(probes,rownames(annot)),]
  rownames(annot) <- names(probes)
  annot$PROBE <- names(probes) ## original probe names
  
  ## ---------------------------------------------------------------------------------
  ## Second pass for missing symbols. Still trying annothub but
  ## missing symbols may map to different keytype. NEED RETHINK:
  ## DO WE REALLY NEED THIS???
  ## ------------------------------------------------------------------------------------
  is.missing <- (is.na(annot$SYMBOL) | annot$SYMBOL == "")
  missing.probes <- probes[which(is.missing)] ## probes match annot!
  missing.probes <- missing.probes[!is.na(missing.probes)]  
  length(missing.probes)
  if (second.pass && length(missing.probes)) {
    dbg("[getGeneAnnotation.ANNOTHUB] second pass: retrying missing",
        length(missing.probes), "symbols...")
    suppressWarnings(suppressMessages(
      missing.probe_type <- detect_probetype(organism, missing.probes, orgdb = orgdb)
    ))
    dbg("[getGeneAnnotation.ANNOTHUB] missing.probe_type=", missing.probe_type)

    ## only do second try if missing.probetype is different
    if (!is.null(missing.probe_type) && !is.na(missing.probe_type)
        && missing.probe_type != probe_type
        ) {      
      missing.probes1 <- match_probe_names(missing.probes, orgdb, missing.probe_type)
      suppressMessages(suppressWarnings(
        missing.annot <- AnnotationDbi::select(orgdb,
          keys = missing.probes1,
          columns = cols,
          keytype = missing.probe_type
        )
      ))
      missing.key <- missing.annot[, missing.probe_type]
      dfA <- apply(missing.annot,2,function(a) tapply(a,missing.key,function(b)
        paste(setdiff(unique(b),c(NA,"")),collapse=';')))
      missing.annot <- data.frame(dfA, check.names=FALSE)
      missing.annot <- missing.annot[match(missing.probes1,rownames(missing.annot)),]
      rownames(missing.annot) <- names(missing.probes)
      missing.annot$PROBE <- names(missing.probes)
      
      # some organisms do not provide SYMBOL but rather GENENAME (e.g. yeast)
      if (!"SYMBOL" %in% colnames(missing.annot)) {
        missing.annot$SYMBOL <- missing.annot$GENENAME
        missing.annot$GENENAME <- missing.annot$ALIAS
      }
      if ("SYMBOL" %in% colnames(missing.annot)) {
        not.symbols <- !(missing.annot$SYMBOL %in% symbols)
        if(length(not.symbols)) missing.annot$SYMBOL[not.symbols] <- NA
      }
      
      for(k in setdiff(colnames(annot),colnames(missing.annot))) {
        missing.annot[[k]] <- NA
      }
      kk <- match(colnames(annot),colnames(missing.annot))
      missing.annot <- missing.annot[,kk,drop=FALSE]
      jj <- match(missing.annot$PROBE, probes)
      ii <- which(!is.na(jj))
      annot[jj[ii], ] <- missing.annot[ii,]
    }
  }

  ## get human ortholog using 'orthogene'
  message("[getGeneAnnotation.ANNOTHUB] getting human orthologs...")
  ortho_organism <- getOrthoSpecies(organism, use="map")
  annot$ORTHOGENE <- getHumanOrtholog(ortho_organism, annot$SYMBOL)$human

  ## Return as standardized data.frame and in the same order as input
  ## probes.
  pkgname <- orgdb$packageName
  if (length(pkgname) == 0) pkgname <- "OrgDb"
  annot$SOURCE <- pkgname[1]
  
  annot.cols <- c(
    "PROBE", "SYMBOL", "ORTHOGENE", "UNIPROT", "GENENAME", 
    ## "GENETYPE", "MAP", "CHR", "POS", "TXLEN", "SOURCE"
    "MAP", "SOURCE"
  )
  missing.cols <- setdiff(annot.cols, colnames(annot))
  missing.cols

  ## create genes data.frame
  genes <- annot 
  for (a in missing.cols) genes[[a]] <- NA
  genes <- genes[, annot.cols]
  new.names <- c(
    "feature", "symbol", "human_ortholog", "uniprot", "gene_title", 
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

getOrthoSpecies <- function(organism, use=c("table","map")[1]) {
  if(use == "map") {
    species <- try(orthogene::map_species(organism, method = "gprofiler", verbose = FALSE))
  }
  if(use == "table") {  
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
    species <- df[sel, "ortho_species"]
  }
  species
}

#' @export
gene2uniprot <- function(genes, organism) {
  gp.organism <- orthogene::map_species(
    species = organism, method = "gprofiler", output_format = "id", verbose = FALSE)
  out <- try(gprofiler2::gconvert(genes, organism = gp.organism, target="UNIPROT_GN_ACC"))
  if(is.null(out) || "try-error" %in% class(out)) return(NULL)
  res <- tapply(out$target, out$input, function(s)
    paste(setdiff(unique(s),c(NA,"")),collapse=';'))
  ii <- match(genes, names(res))
  ## seems input is uppercase!
  ii <- ifelse(is.na(ii), match(toupper(genes),toupper(names(res))), ii)
  res[ii]
}

#' @export
uniprot2gene <- function(uniprots, organism) {
  gp.organism <- orthogene::map_species(
    species = organism, method = "gprofiler", output_format = "id", verbose = FALSE)
  out <- try(gprofiler2::gconvert(uniprots, organism = gp.organism, target="ENSG"))
  if(is.null(out) || "try-error" %in% class(out)) return(NULL)
  res <- tapply(out$name, out$input, function(s) paste(setdiff(unique(s),c(NA,"")),collapse=';'))
  res[uniprots]
}

#' Clean up inline duplicated features: eg: feature1;feature1;....
#'
#' @export
clean_dups_inline_probenames <- function(probes) {
  probes[is.na(probes)] <- ""
  probes <- sapply(strsplit(probes,split="[;,._]"),function(s) paste(unique(s),collapse=";"))
  return (probes)
}

#' non-greedy removal of numerical postfix. Postfix is defined as (1)
#' last numerical substring after - (minus), or (2) any substring
#' after special separators [_.]. 
#'
strip_postfix <- function(s) {  
  stripFUN <- function(s) {  
    sub(paste0("[._].*$|[-][0-9.]+$"),"",s)
  }
  ss <- strsplit(s,split=';')
  ss <- lapply(ss, function(s) stripFUN(s))
  sapply(ss, paste, collapse=';')
}


#' non-greedy removal of prefixes. Prefix is defined as any
#' alphanumerical substring (no spaces, no special chars) before the
#' matching colon character :.
#'
strip_prefix <- function(s) {  
  stripFUN <- function(s) {  
    sub("^[0-9a-zA-Z]+:","",s)
  }
  ss <- strsplit(s,split=';')
  ss <- lapply(ss, function(s) stripFUN(s))
  sapply(ss, paste, collapse=';')
}

#' Cleanup symbols names from postfixes and prefixes. Take only first
#' symbol. This is mostly used for symbol lookup tables that need one
#' clean symbol.
#'
#' @export
clean_symbols <- function(symbols) {
  strip_prefix(strip_postfix(sub(";.*","",trimws(symbols))))
}

#' Cleanup probe names from postfixes or version numbers. Retains
#' prefix needed for multi-omics.
#'
#' @export
clean_probe_names <- function(probes, sep = "_.-") {
  probes0 <- trimws(probes)
  probes[is.na(probes)] <- ""  
  ## strip multiple probes  
  probes <- sub("[;].*","",probes) 
  ## strip away anything postfix after a 'dot' or 'underscore'
  probes <- sub(paste0("[", sep, "].*"), "", probes)
  ##probes <- strip_postfix(probes)
  return(probes)
}

#' Match dirty probe names to clean key names
#'
#' @export
match_probe_names <- function(probes, orgdb, probe_type = NULL) {
  if (is.character(orgdb)) orgdb <- getOrgDb(orgdb)
  if (is.null(orgdb)) {
    message("[match_probe_names] ERROR could not get orgdb!" )
    return(NULL)
  }
  if (is.null(probe_type)) {
    probe_type <- detect_probetype(organism = "custom", probes, orgdb = orgdb)
  }
  probe.names <- names(probes)
  all.keys <- keys(orgdb, probe_type)
  tsub <- function(s) gsub("[-:;.]|\\[|\\]", ".", s)
  ii <- match(toupper(tsub(probes)), toupper(tsub(all.keys)))
  table(is.na(ii))
  new.probes <- all.keys[ii]
  if (sum(is.na(new.probes))) {
    jj <- which(is.na(new.probes))
    new.probes[jj] <- probes[jj]
    jj.probes <- clean_probe_names(probes[jj], sep = ".-")
    ii <- match(toupper(tsub(jj.probes)), toupper(tsub(all.keys)))
    if (any(!is.na(ii))) {
      k <- which(!is.na(ii))
      new.probes[jj[k]] <- all.keys[ii[k]]
    }
  }
  names(new.probes) <- probe.names
  new.probes
}


#' Cleanup annotation
#'
cleanupAnnotation <- function(genes) {

  if(is.null(genes)) {
    return(NULL)
  }
  
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
getCustomAnnotation <- function(probes, custom_annot) {
  message("[getCustomAnnotation] Adding custom annotation table...")
  # If the user has provided a custom gene table, check it and use it
  custom_annot <- data.frame(custom_annot, check.names = FALSE)

  if (!"feature" %in% colnames(custom_annot) && !is.null(rownames(custom_annot))) {
    custom_annot$feature <- rownames(custom_annot)
  }
  
  annot_map <- list(
    "human_ortholog" = "",
    "gene_title" = "unknown",
    "chr" = "unknown",
    "source" = "custom"
  )

  required_cols <- c("feature", "symbol", "gene_name")
  if (!"symbol" %in% colnames(custom_annot)) {
    custom_annot$symbol <- custom_annot$feature
  }
  if (!"gene_name" %in% colnames(custom_annot)) {
    custom_annot$gene_name <- custom_annot$feature
  }

  # this will be used at the end to order df columns
  table_col_order <- c(required_cols, names(annot_map))

  # legacy code but maybe this could be removed in the future...
  # required_in_annot <- all(required_cols %in% colnames(custom_annot))

  ##  if (!is.null(custom_annot) && num_annot > 1 && required_in_annot) {
  num_annot <- sum(probes %in% custom_annot$feature)
  if (!is.null(custom_annot) && num_annot > 1) {  
    message("[getCustomAnnotation] Cleaning custom annotation table...")
    # remove all NA columns, otherwise the for loop below will not work
    custom_annot <- custom_annot[, !apply(custom_annot, 2, function(x) all(is.na(x)))]

    # identify missing columns and fill them with annot_map
    missing_cols <- setdiff(names(annot_map), names(custom_annot))
    custom_annot[missing_cols] <- annot_map[missing_cols]

    # filter annotated table by probes using match
    custom_annot <- custom_annot[match(probes, custom_annot$feature), ]
    
    # if row was missing from annotation table (NA from match call
    # above), input NA based on probes
    rownames(custom_annot) <- probes
    nr <- nrow(custom_annot)
    if(is.null(custom_annot$feature)) custom_annot$feature <- rep(NA,nr)
    if(is.null(custom_annot$symbol)) custom_annot$symbol <- rep(NA,nr)
    if(is.null(custom_annot$gene_name)) custom_annot$gene_name <- rep(NA,nr)    
    cname <- rownames(custom_annot)
    custom_annot$feature <- ifelse(is.na(custom_annot$feature), cname, custom_annot$feature)
    custom_annot$symbol <- ifelse(is.na(custom_annot$symbol), cname, custom_annot$symbol)
    custom_annot$gene_name <- ifelse(is.na(custom_annot$gene_name), cname, custom_annot$gene_name)

    # Fill NA values with corresponding values from annot_map
    res <- lapply(names(annot_map), function(x) {
      ifelse(is.na(custom_annot[[x]]), annot_map[[x]], custom_annot[[x]])
    })    
    names(res) <- names(annot_map)
    res <- as.data.frame(res)
    res$source <- ifelse(res$source == "custom", "custom", paste0("custom+", res$source))
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
      chr = "unknown",
      source = "custom"
    )
    rownames(custom_annot) <- probes
  }
  
  custom_annot <- custom_annot[, table_col_order, drop=FALSE]
  custom_annot <- cleanupAnnotation(custom_annot)

  return(custom_annot)
}

#' Like getCustomAnnotation() but allows custom column names mapping
#' to feature, symbol and title. Also simplified
#' implementation. Should eventually replace getCustomAnnotation().
#' 
#' @export
getCustomAnnotation2 <- function(probes, custom_annot, feature.col='feature',
                                 symbol.col='symbol', gene_title.col='gene_title',
                                 ortholog.col='human_ortholog',
                                 extra.columns = TRUE) {

#  feature.col='feature';symbol.col='symbol';gene_title.col='gene_title';ortholog.col='human_ortholog';extra.columns = TRUE
  
  message("[getCustomAnnotation2] Adding custom annotation table...")
  # Create custom gene table from probe names
  message("[getCustomAnnotation2] Creating annotation table from probe names...")
  annot <- data.frame(
    feature = probes,
    symbol = probes,
    gene_name = probes,
    human_ortholog = NA,
    gene_title = "unknown",
    ## chr = NA,
    source = "custom"
  )
  rownames(annot) <- probes
  required.columns <- colnames(annot)
  
  # If the user has provided a custom gene table, check it and use it
  if(!is.null(custom_annot)) {

    custom_annot <- data.frame(custom_annot, check.names=FALSE)

    if(!feature.col %in% colnames(custom_annot)) {
      if(!is.null(rownames(custom_annot))) {
        custom_annot$rownames <- rownames(custom_annot)
      }
      fsum <- apply(custom_annot,2,function(a) mean(probes %in% a,na.rm=TRUE))
      feature.col <- NULL
      if(max(fsum)>0.9) feature.col <- names(which.max(fsum))
      if(length(feature.col)==0) {
        custom_annot$feature <- probes
        feature.col <- "feature"
      }
    }
    if(!symbol.col %in% colnames(custom_annot)) {
      symbol.col <- head(grep("symbol|name|gene|protein|alias",
        setdiff(colnames(custom_annot),feature.col),
        ignore.case=TRUE, value=TRUE),1)
      if(length(symbol.col)==0) symbol.col <- NA
    }
    if(!gene_title.col %in% colnames(custom_annot)) {
      gene_title.col <- head(grep("title|description|name",
        setdiff(colnames(custom_annot),c(feature.col,symbol.col)),
        ignore.case=TRUE, value=TRUE),1)
      if(length(gene_title.col)==0) gene_title.col <- NA
    }
    if(!ortholog.col %in% colnames(custom_annot)) {
      ortholog.col <- head(grep("ortholog|human|hgnc",
        setdiff(colnames(custom_annot),c(feature.col,symbol.col)),
        ignore.case=TRUE, value=TRUE),1)
      if(length(ortholog.col)==0) ortholog.col <- NA
    }
    dbg("[getCustomAnnotation2] feature.col = ", feature.col)
    dbg("[getCustomAnnotation2] symbol.col = ", symbol.col)
    dbg("[getCustomAnnotation2] title.col = ", gene_title.col)
    dbg("[getCustomAnnotation2] ortholog.col = ", ortholog.col)    

    features <- custom_annot[,feature.col]
    custom_annot <- custom_annot[match(probes, features),]

    # Rename columns
    newcols <- c("feature" = feature.col, "symbol" = symbol.col,
                 "gene_title" = gene_title.col,
                 "human_ortholog" = ortholog.col)
    newcols <- newcols[which(newcols != names(newcols))]
    newcols <- newcols[which(newcols %in% colnames(custom_annot))]
    if(length(newcols)) {
      custom_annot <- dplyr::rename(custom_annot, all_of(newcols))
    }
    
    if("feature" %in% colnames(custom_annot)) {
      annot$feature <- custom_annot[,"feature"]
    }
    if("symbol" %in% colnames(custom_annot)) {
      annot$symbol  <- custom_annot[,"symbol"]
    }
    if("gene_title" %in% colnames(custom_annot)) {    
      annot$gene_title <- custom_annot[,"gene_title"]
    }

    ##  if (!is.null(custom_annot) && num_annot > 1 && required_in_annot) {
    # remove all NA columns, otherwise the for loop below will not work
    custom_annot <- custom_annot[, colMeans(is.na(custom_annot))!=1, drop=FALSE]

    # identify missing columns and fill them with annot_map
    missing_cols <- setdiff( colnames(custom_annot), colnames(annot))
    missing_cols <- setdiff(missing_cols, c(NA))
    annot <- cbind(annot, custom_annot[,missing_cols,drop=FALSE])
  }

  if(!extra.columns) {
    sel <- (colnames(annot) %in% required.columns)
    annot <- annot[,sel]
  }
  
  message("[getCustomAnnotation2] Cleaning custom annotation table...")
  annot <- cleanupAnnotation(annot)
  return(annot)
}



## ================================================================================
## ================== GET ORTHOLOG FUNCTIONS ======================================
## ================================================================================

#' @title Get human ortholog from given symbols of organism by using
#'   orthogene package. This package needs internet connection.
#'
#' @export
getHumanOrtholog <- function(organism, symbols,
                             ortho.methods=c("homologene", "gprofiler", "babelgene"),
                             verbose=1) {

  orthogenes <- rep(NA, length(symbols))
  orthosource <- rep(NA, length(symbols))    

  ## clean symbols
  orig.symbols <- symbols
  symbols <- clean_symbols(symbols)
  
  ## Try mapping with orthogene's databases 
  ##ortho.methods <- c("gprofiler", "homologene", "babelgene") ## mapping methods
  ortho_organism <- getOrthoSpecies(organism)
  ortho.found <- FALSE
  i=1
  while(i <= length(ortho.methods) && !ortho.found) {
    ortho.out <- try(orthogene::convert_orthologs(
      gene_df = c("---", unique(symbols[!is.na(symbols)])),
      input_species = ortho_organism,
      output_species = "human",
      method = ortho.methods[i],
      non121_strategy = "drop_both_species",
      verbose = FALSE
    ), silent=TRUE)
    class(ortho.out)
    results.ok <- (!"try-error" %in% class(ortho.out) &&
                     inherits(ortho.out, "data.frame") &&
                     nrow(ortho.out) > 0)
    results.ok
    if (results.ok) {
      ii <- which(is.na(orthogenes))
      ##ii <- head(which(is.na(orthogenes)),10)      
      jj <- match(symbols[ii], ortho.out$input_gene)
      kk <- ii[which(!is.na(jj))]
      jj <- jj[which(!is.na(jj))]
      if(verbose>0) message("[getHumanOrtholog] mapping ",length(kk)," symbols with ", ortho.methods[i])
      orthogenes[kk] <- rownames(ortho.out)[jj]
      orthosource[kk] <- ortho.methods[i]
    } else {
      if(verbose>0) message("[getHumanOrtholog] failed lookup: ", ortho.methods[i])
    }
    ortho.found <- all(!is.na(orthogenes))
    i <- i + 1
  }

  table(is.na(orthogenes))
  mean.mapped <- round( 100*mean(!is.na(orthogenes)), digits=4)
  orthogene.failed <- (mean.mapped < 10.0)
  orthogene.failed
  if (orthogene.failed) {
    if(verbose>0) message("[getHumanOrtholog] ratio mapped using orthogene = ", mean.mapped,"%")
    if(verbose>0) message("[getHumanOrtholog] Trying biomart...")
    ## test if biomart is reachable
    ii <- which(is.na(symbols))
    res.biomart <- try(getHumanOrtholog.biomart(organism, symbols[ii]), silent=TRUE)
    class(res.biomart)
    if (!"try-error" %in% class(res.biomart)) {
      jj <- which(is.na(res.biomart))
      ii <- ii[jj]
      orthogenes[ii] <- res.biomart[jj]
      orthosource[ii] <- 'biomart'
    }
  } else {
    if(verbose>0) message("[getHumanOrtholog] skipping biomart...")
  }

  ## Map any missing symbols that look like human genes
  human.genes <- playdata::GENE_SYMBOL
  ii <- which(is.na(orthogenes) & toupper(symbols) %in% human.genes)
  orthogenes[ii] <- toupper(symbols[ii])

  mean.mapped <- round( 100*mean(!is.na(orthogenes)), digits=4)  
  if(verbose>0) message("[getHumanOrtholog] total ratio mapped  = ", mean.mapped,"%")
  
  ## return dataframe. First column organism symbols, second column
  ## human ortholog. NA if missing.
  df <- data.frame(input=orig.symbols, symbols, "human"=orthogenes, source=orthosource)
  colnames(df)[2] <- organism
  return(df)
}

#' @title Get human ortholog from given symbols of organism by using
#'   biomart package. This package needs internet connection.
#' This is an alternative to orthogene::map_genes failure.
#' Unfortunately, biomart is *much* less reliable than orthogene (often down).
#'
#' @export
getHumanOrtholog.biomart <- function(organism, symbols, verbose=1) {

  if (tolower(organism) == "human") organism <- "Homo sapiens"
  if (tolower(organism) == "mouse") organism <- "Mus musculus"
  if (tolower(organism) == "rat") organism <- "Rattus norvegicus"
  if (tolower(organism) %in% c("dog","canis familiaris")) organism <- "Canis LFamiliaris"  

  if(verbose>0) message("[getHumanOrtholog.biomart] Mapping ", organism, " genes with biomart.")
  require(biomaRt)
  s1 <- tolower(substring(organism,1,1))
  s2 <- tolower(strsplit(organism, " ")[[1]][2])
  organism0 <- paste0(s1, s2)
  organism0
  if(verbose>0) message("[getHumanOrtholog.biomart] Searching biomart for '", organism0,"'")
  dd <- listDatasets(useEnsembl(biomart = "genes"))
  hh <- grep(organism0, dd$dataset)
  hh
  if(length(hh)==0) {
    message("ERROR: ",paste0(organism, " not found in biomart databases. Exiting."))
    return(NULL)
  }

  dataset <- dd$dataset[hh]
  if(verbose>0) message("[getHumanOrtholog.biomart] found matching dataset '", dataset,"'")
  organism_mart <- NULL
  human_mart <- NULL    
  organism_mart <- biomaRt::useEnsembl(biomart = "genes", dataset = dataset)
  human_mart <- biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  mart1.ok <- (!is.null(organism_mart) && !"try-error" %in% class(organism_mart) &&
                 class(organism_mart) == "Mart")
  mart2.ok <- (!is.null(human_mart) && !"try-error" %in% class(human_mart) &&
                 class(human_mart) == "Mart")
  if(verbose>0) dbg("[getHumanOrtholog.biomart] organism mart OK -> ", mart1.ok)
  if(verbose>0) dbg("[getHumanOrtholog.biomart] human mart OK -> ", mart2.ok)    
  if(!(mart1.ok && mart2.ok)) {
    message("[getHumanOrtholog.biomart] ERROR: Could not create Marts...")
    return(NULL)
  }

  if (organism == "Mus musculus") {
    attrs <- c("ensembl_gene_id", "mgi_symbol")
    flt <- "mgi_symbol"
  } else {
    attrs <- c("ensembl_gene_id", "external_gene_name")
    flt <- "external_gene_name"
  }
  if(verbose>0) message("[getHumanOrtholog.biomart] Testing if biomart is reachable...")
  nz.symbols <- symbols[!is.na(symbols) & symbols!='']
  res.biomart <- try(biomaRt::getLDS(
    attributes = attrs,
    filters = flt,
    # values = "Xist",
    values = head(nz.symbols),
    mart = organism_mart,
    attributesL = c("ensembl_gene_id", "hgnc_symbol"),
    martL = human_mart
  ), silent=TRUE)
  class(res.biomart)
  if ("try-error" %in% class(res.biomart)) {
    message("[getHumanOrtholog] ERROR: biomart::getLDS failed to contact server")
    return(NULL)
  }

  ## Full query
  orthogenes <- try(biomaRt::getLDS(
    attributes = attrs,
    filters = flt,
    values = symbols,
    mart = organism_mart,
    attributesL = c("ensembl_gene_id", "hgnc_symbol"),
    martL = human_mart
  ), silent=TRUE)
  class(orthogenes)
  
  if ("try-error" %in% class(res.biomart)) {
    message("[getHumanOrtholog] ERROR: full biomart::getLDS query failed")
    return(NULL)
  }

  ## succesful
  return(orthogenes)
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
probe2symbol <- function(probes, annot_table, query = c("symbol", "gene_name"),
                         key = NULL, fill_na = FALSE) {
  # Prepare inputs
  annot_table <- cbind(rownames = rownames(annot_table), annot_table)
  probes1 <- setdiff(probes, c(NA, ""))
  if (is.null(key) || !key %in% colnames(annot_table)) {
    key <- which.max(apply(annot_table, 2, function(a) sum(probes1 %in% a)))
  }
  if (is.null(key)) {
    stop("[probe2symbol] FATAL. could not get key column.")
  }

  # match query
  ii <- match(probes, annot_table[, key])
  query <- intersect(query, colnames(annot_table))
  if (length(query) == 0) {
    stop("ERROR. no symbol column.")
  }
  query_col <- annot_table[ii, query[1]]

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


#' @title Detect probe type from probe set
#' @export
detect_probetype <- function(organism, probes, orgdb = NULL,
                             nprobe = 1000, use.ah = NULL, datatype = NULL,
                             verbose = TRUE) {
  if (tolower(organism) == "human") organism <- "Homo sapiens"
  if (tolower(organism) == "mouse") organism <- "Mus musculus"
  if (tolower(organism) == "rat") organism <- "Rattus norvegicus"

  if (is.null(datatype) && all(grepl("[:]", probes))) {
    dbg("[detect_probetype] datatype is multi-omics")
    datatype <- "multi-omics"
  }

  if (!is.null(datatype) && datatype == "metabolomics") {
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
    "ALIAS", "MGI", "TAIR", ## organism specific
    "ENSEMBLTRANS", "ENSEMBLPROT",
    "REFSEQ", "ENTREZID"
  )
  keytypes <- intersect(keytypes, AnnotationDbi::keytypes(orgdb))
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

  ## try different cleaning methods. NEED RETHINK!!!! refseq has
  ## underscore!
  probes0 <- probes
  probes1 <- clean_probe_names(probes)
  probesx <- unique(c(probes0, probes1))

  ## Get all organism symbols
  org_annot <- AnnotationDbi::select(
    orgdb,
    keys = keys(orgdb, "ENTREZID"),
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

  ## Return top match
  ##  key_matches
  top_match <- NULL
  if (all(key_matches == 0)) {
    if (verbose) {
      message("head.probes = ", paste(head(probes), collapse = " "))
      message("WARNING: Probe type not found. Valid probe types: ", paste(keytypes, collapse = " "))
    }
    # fallback before giving up; try gprofiler
    gp.organism <- orthogene::map_species(species = organism, method = "gprofiler", 
      output_format = "id", verbose = FALSE)
    gp.out <- gprofiler2::gconvert(probesx, organism = gp.organism, target="UNIPROT_GN_ACC")
    if (!is.null(gp.out)) {
      key_matches['GPROFILER'] <- length(unique(gp.out$input)) / length(probesx)
    }
  }


  if (max(key_matches, na.rm = TRUE) < 0.01) {
    message("WARNING: Insufficient matching ratio. Max match = ", max(key_matches, na.rm = TRUE))
    return(NA)
  }
  if (max(key_matches, na.rm = TRUE) < 0.50) {
    message("WARNING: Low matching ratio. Max match = ", max(key_matches, na.rm = TRUE))
  }
  top_match <- names(which.max(key_matches))
  return(top_match)
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
    sel.na <- which(annot[,k] %in% c(NA,"","-","---","NA"))
    annot[sel.na,k] <- '---'
    annot[,k] <- toupper(annot[,k])  ## all uppercase??
    map.obj <- rename_by(obj, annot_table = annot, new_id = target[1])
  }
  if (!is.null(dim(map.obj))) rownames(map.obj) <- toupper(rownames(map.obj))
  if (is.null(dim(map.obj))) names(map.obj) <- toupper(names(map.obj))
  map.obj
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
    message("[showProbeTypes] ERROR: unsupported organism '", organism, "'\n")
    return(NULL)
  }

  ## get probe types for organism
  if (!is.null(keytypes) && keytypes[1] == "*") {
    keytypes <- AnnotationDbi::keytypes(orgdb)
  }
  if (is.null(keytypes)) {
    keytypes <- c(
      "SYMBOL", "ENSEMBL", "UNIPROT", "ENTREZID",
      "GENENAME", "MGI", "TAIR",
      "ENSEMBLTRANS", "ENSEMBLPROT",
      "ACCNUM", "REFSEQ"
    )
  }
  keytypes0 <- keytypes
  keytypes <- intersect(keytypes, AnnotationDbi::keytypes(orgdb))
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
  if (is.null(orgdb)) {
    message("[getOrganismGO] WARNING unsupported organism? ", organism)
    return(NULL)
  }
  
  go.gmt <- list()
  AnnotationDbi::keytypes(orgdb)
  ont_classes <- c("BP", "CC", "MF")
  if (!"GOALL" %in% AnnotationDbi::keytypes(orgdb)) {
    message("WARNING:: missing GO annotation in database!\n")
  } else {
    ## create GO annotets
    message("Creating GO annotation from AnnotationHub...")
    ont_classes <- c("BP", "CC", "MF")
    k <- "BP"
    for (k in ont_classes) {
      gene.column <- intersect( c("SYMBOL","GENENAME","MGI","ALIAS"), columns(orgdb))
      gene.column <- head(gene.column,1)
      if(length(gene.column)>0) {
        suppressMessages(suppressWarnings(
          go_id <- AnnotationDbi::mapIds(
            orgdb, keys = k, keytype = "ONTOLOGY",
            column = "GO", multiVals = "list"
          )[[1]]
        ))
        go_id <- unique(go_id)
        suppressMessages(suppressWarnings(
          sets <- AnnotationDbi::mapIds(
            orgdb, keys = go_id, keytype = "GOALL",
            column = gene.column, multiVals = "list"
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

  ## map given organism to ORTHOGENE species name
  species <- getOrthoSpecies(organism, use="map")
  if ("try-error" %in% class(species)) {
    message("[getGeneAnnotation.ORTHOGENE] *WARNING* could not connect to server")
    return(NULL)
  }
  if (is.null(species)) {
    message("ERROR: unknown organism ", organism)
    return(NULL)
  }
  message("[getGeneAnnotation.ORTHOGENE] mapping to species: ", species)
  probes[is.na(probes) | probes==""] <- "NA"
  probes <- make.unique(probes)
  probes1 <- clean_probe_names(probes)  ## CHECK dots!!

  gene.out <- try(orthogene::map_genes(
    genes = probes1,
    species = species,
    verbose = FALSE
  ), silent=TRUE)

  ortholog <- toupper(probes1)
  if (!inherits(gene.out, "try-error")) {
    gene.out <- gene.out[match(probes1, gene.out$input), ]
    ortholog <- getHumanOrtholog(organism, gene.out$name)$human
  }

  ## get uniprot id
  gp.organism <- orthogene::map_species(species = species, method = "gprofiler", 
    output_format = "id", verbose = FALSE)
  gp.out <- gprofiler2::gconvert(probes1, organism = gp.organism, target="UNIPROT_GN_ACC")
  uniprot <- tapply(gp.out$target, gp.out$input, function(x) paste(x,collapse=";"))
  uniprot <- as.character(uniprot[match(probes1,names(uniprot))])
  
  df <- data.frame(
    feature = probes,
    symbol = gene.out$name,
    human_ortholog = ortholog,
    uniprot = as.character(uniprot),
    gene_title = sub(" \\[.*", "", gene.out$description),
    chr = NA,
    source = "gprofiler2",
    gene_name = probes
  )
  rownames(df) <- probes
  
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
  annot$rownames <- rownames(annot)
  new.feature <- strsplit(annot[,target[1]],split=';')
  for (i in 2:length(target)) {
    if (target[i] %in% colnames(annot)) {
      new.feature <- mapply(paste0, new.feature, annot[, target[i]])
    } else {
      ## some character
      new.feature <- mapply(paste0, new.feature, target[i])      
    }
  }
  new.feature <- sapply(new.feature, paste0, collapse=";")
  if (sum(duplicated(new.feature)) > 0) {
    new.feature <- make_unique(new.feature)
  }
  new.feature
}


#' @export
pgx.getFeatureInfo <- function(pgx, feature) {

  if (is.null(feature) || length(feature) == 0) {
    return(NULL)
  }
  if (is.na(feature) || feature == "") {
    return(NULL)
  }
  
  annot <- as.list( pgx$genes[feature,] )
  names(annot) <- sub("uniprot","protein",names(annot))
  annot$gene_name <- NULL
  annot$chr <- NULL
  annot$pos <- NULL
  annot$tx_len <- NULL

  datatype <- pgx$datatype
  if(pgx$datatype == "multi-omics") {
    ##dtype <- c("mx"="metabolomics","px"="proteomics","gx"="transcriptomics")
    datatype <- ifelse(grepl("^mx:|metabolomics",feature), "metabolomics", datatype)
  }
  
  if(!"protein" %in% names(annot) && datatype!="metabolomics") {
    annot[["protein"]] <- gene2uniprot(annot$symbol, pgx$organism)
  }
  
  if(annot$human_ortholog %in% names(playdata::GENE_SUMMARY)) {
    annot.summary <- playdata::GENE_SUMMARY[annot$human_ortholog]
    annot.summary <- gsub("Publication Note.*|##.*", "", annot.summary)
    annot[["summary"]] <- annot.summary
  }
  
##  annot[["GO"]] <- pgx.get_goterms(pgx, annot$symbol) 
  if(datatype == "metabolomics") {
    annot <- getMetaboliteInfo(
      organism = pgx$organism,
      id = annot$symbol,
      info = annot)
  } else {
    annot <- info.add_hyperlinks(
      annot, feature, datatype, 
      nm.symbol='symbol',
      nm.prot='protein',
      nm.ortholog='human_ortholog',
      as.link=TRUE, add.summary=FALSE
    )
  }
  
  if (!is.null(annot)) {
    annot <- annot[!duplicated(names(annot))]
    ## reorder
    nn1 <- intersect(
      c(
        "gene_symbol", "organism", "name", "map_location",
        "uniprot", "databases", "summary", names(annot)
      ),
      names(annot)
    )
    nn2 <- setdiff(names(annot), nn1)
    annot <- annot[c(nn1,nn2)]
    names(annot) <- sub("gene_symbol", "symbol", names(annot))
    names(annot) <- sub("uniprot", "protein", names(annot))
    names(annot) <- sub("map_location", "genome location", names(annot))
  } else {
    annot <- list()
    annot$summary <- "(no info available)"
  }
  
  return(annot)

}

#' Retrieve gene information from different databases (AnnotHub,
#' orgdb) and wrap hyperlinks.
#' 
#' @export
getOrgGeneInfo <- function(organism, gene, feature, ortholog, datatype,
                           as.link = TRUE) {
  if (is.null(gene) || length(gene) == 0) {
    return(NULL)
  }
  if (is.na(gene) || gene == "") {
    return(NULL)
  }

  orgdb <- getOrgDb(organism, use.ah = NULL)
  if (is.null(orgdb)) {
    message("[getOrgGeneInfo] WARNING could not get orgdb!" )
    return(NULL)
  }
  cols <- c("SYMBOL", "UNIPROT", "GENENAME", "MAP", "OMIM", "PATH", "GO")
  cols <- intersect(cols, AnnotationDbi::keytypes(orgdb))

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
  
  annot <- AnnotationDbi::select(
    orgdb,
    keys = gene,
    keytype = keytype,
    columns = cols
  )
  
  names(info) <- cols

  ## take out duplicates
  info <- lapply(info, unique)

  info[["ORGANISM"]] <- organism
  info <- info.add_hyperlinks(info, feature, datatype, as.link = TRUE)   

  ## rename
  tags <- c(
    "ORGANISM", "SYMBOL", "UNIPROT", "GENENAME", "MAP", "OMIM", "PATH",
    "GO", "SUMMARY", "DATABASES"
  )
  tags <- intersect(tags, names(info))
  info <- info[tags]
  new.names <- c(
    "ORGANISM" = "organism", "SYMBOL" = "gene_symbol", "UNIPROT" = "uniprot",
    "GENENAME" = "name", "MAP" = "map_location",
    "OMIM" = "OMIM", "PATH" = "pathway", "GO" = "GO",
    "SUMMARY" = "summary", "DATABASES" = "databases"
  )
  names(info) <- new.names[tags]
  return(info)
}

pgx.get_goterms <- function(pgx, symbol) {
  if(!symbol %in% rownames(pgx$GMT)) return(NULL)
  jj <- grep("^GO",colnames(pgx$GMT),value=TRUE)
  terms <- names(which( pgx$GMT[symbol,jj] !=0 ))
  go_ids <- paste0("GO:",gsub(".*\\(GO_|\\)","",terms))
  go_ids
}

info.add_hyperlinks <- function(info, feature, datatype,
                                nm.symbol='SYMBOL', nm.prot='UNIPROT',
                                nm.ortholog='ORTHOLOG',
                                as.link=TRUE, add.summary=TRUE) {

  ##nm.symbol='SYMBOL';nm.prot='UNIPROT';nm.ortholog='ORTHOLOG';as.link=TRUE;add.summary=TRUE
  ##nm.symbol='symbol';nm.prot='protein';nm.ortholog='human_ortholog';as.link=TRUE;add.summary=TRUE
  
  symbol  <- info[[nm.symbol]]
  ortholog  <- info[[nm.ortholog]]
  uniprot <- info[[nm.prot]]
  
  if (length(uniprot) == 0) {
    this.uniprot <- NULL
  } else if(length(uniprot) == 1 && uniprot[1]=="") {
    this.uniprot <- NULL
  } else {
    uniprot <- sort(unique(unlist(strsplit(uniprot, split=';'))))
    this.uniprot <- uniprot[which(sapply(uniprot, function(p) grepl(p, feature)))]
    if (length(this.uniprot) == 0) this.uniprot <- rev(uniprot)[1]
  }
  
  if (as.link && length(symbol)) {
    gene.link <- "<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=GENE' target='_blank'>GENE</a>"
    gene.link <- sapply(symbol, function(s) gsub("GENE", s, gene.link))
    info[[ nm.symbol ]] <- paste(gene.link, collapse = ", ")
  }
  
  if (as.link && length(ortholog)) {
    gene.link <- "<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=GENE' target='_blank'>GENE</a>"
    ortho.link <- sapply(ortholog, function(s) gsub("GENE", s, gene.link))
    info[[ nm.ortholog ]] <- paste(ortho.link, collapse = ", ")
  }
  
  if (as.link && length(uniprot)) {
    prot.link <- "<a href='https://www.uniprot.org/uniprotkb/UNIPROT' target='_blank'>UNIPROT</a>"
    prot.link <- sapply(uniprot, function(s) gsub("UNIPROT", s, prot.link))
    info[[ nm.prot ]] <- paste(prot.link, collapse = ", ")
  }

  ## create link to external databases: OMIM, GeneCards, Uniprot
  if (as.link) {
    genecards.link <- "<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=GENE' target='_blank'>GeneCards</a>"
    uniprot.link <- "<a href='https://www.uniprot.org/uniprotkb/UNIPROT' target='_blank'>UniProtKB</a>"
    if (length(symbol)) genecards.link <- sub("GENE", symbol[1], genecards.link)
    if (length(this.uniprot)) uniprot.link <- sub("UNIPROT", this.uniprot, uniprot.link)
    db.links <- paste(c(genecards.link, uniprot.link), collapse = ", ")
    if(db.links=="") db.links <- NULL
    info[["DATABASES"]] <- db.links
  }
  
  if (length(this.uniprot) && grepl("proteomics", datatype, ignore.case = TRUE)) {
    ## create link to PhosphoSitePlus
    phosphositeplus.link <- "<a href='https://www.phosphosite.org/simpleSearchSubmitAction.action?searchStr=GENE' target='_blank'>PhosphoSitePlus</a>"
    phosphositeplus.link <- "<a href='https://www.phosphosite.org/uniprotAccAction?id=UNIPROT' target='_blank'>PhosphoSitePlus</a>"
    ## phosphositeplus.link <- sub("GENE", symbol[1], phosphositeplus.link)
    phosphositeplus.link <- sub("UNIPROT", this.uniprot, phosphositeplus.link)
    info[["DATABASES"]] <- paste(c(info[["DATABASES"]], phosphositeplus.link), collapse = ", ")

    ## ## create links to PhosphoELM for proten and gene: db of S/T/Y phosphorylation sites
    ## phosphoELM.link <- "<a href='http://phospho.elm.eu.org/byAccession/UNIPROT' target='_blank'>PhosphoELM</a>"
    ## phosphoELM.link <- sub("UNIPROT", uniprot, phosphoELM.link)
    ## info[["DATABASES"]] <- paste(c(info[["DATABASES"]], phosphoELM.link), collapse = ", ")
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
      envir = GO.db::GOTERM, ifnotfound = NA )[[1]])))
    go.ok <- (!"try-error" %in% class(try.out))
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
  
  if(add.summary) {
    ## pull summary
    info[["SUMMARY"]] <- "(no info available)"
    ortholog <- info[[nm.ortholog]]
    if (!is.null(ortholog) && ortholog %in% names(playdata::GENE_SUMMARY)) {
      info[["SUMMARY"]] <- playdata::GENE_SUMMARY[ortholog]
      info[["SUMMARY"]] <- gsub("Publication Note.*|##.*", "", info[["SUMMARY"]])
    }
  }
  
  return(info)
}


#' Automatically detects species by trying to detect probetype from
#' list of test_species. Warning. bit slow.
#'
#' @export
check_species_probetype <- function(
    probes,
    test_species = c("Human", "Mouse", "Rat"),
    datatype = NULL, annot.cols = NULL) {

  ## No check if custom
  custom_datatype <- !is.null(datatype) && tolower(datatype) %in% c("custom","unknown","")
  custom_organism <- any(tolower(test_species) %in% c("custom","unknown","no organism"))

  if(custom_datatype || custom_organism) {
    out <- rep("custom", length(test_species))
    names(out) <- test_species
    return(as.list(out))
  }

  probes <- unique(clean_probe_names(probes))
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
        use.ah = FALSE,
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
  valid_name <- mean(grepl("[_][1-9]+", features), na.rm = TRUE) > 0.9
  valid_name
  uniprot <- sub("[_].*", "", features)
  positions <- gsub(".*[_]|[.].*", "", features)
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
  if (tolower(organism) == "human") organism <- "Homo sapiens"
  if (tolower(organism) == "mouse") organism <- "Mus musculus"
  if (tolower(organism) == "rat") organism <- "Rattus norvegicus"

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

#' Annotate multi-omics probetype. Probe names *must  be prefixed with
#' data type unless classical transcriptomics/proteomics.
#'
getMultiOmicsProbeAnnotation <- function(organism, probes) {
  
  if(all(grepl("^[A-Za-z]+:",probes))) {
    dtype <- sub(":.*","",probes)
  } else {
    ## no colon in names it is single type probes. try to guess by
    ## matching.
    ptype <- detect_probetype(organism, probes)
    mtype <- mx.detect_probetype(probes)
    dbg("[getMultiOmicsProbeAnnotation] ptype =", ptype)
    dbg("[getMultiOmicsProbeAnnotation] mtype =", mtype)
    gx.types <- c(
      "SYMBOL", "ENSEMBL", "ACCNUM", "GENENAME",
      "MGI", "TAIR", "ENSEMBLTRANS", "REFSEQ", "ENTREZID"
    )
    px.types <- c("UNIPROT", "ENSEMBLPROT")
    if (!is.na(ptype)) {
      dx <- ifelse(ptype %in% px.types, "px", "gx")
    } else if (!is.na(mtype)) {
      dx <- mtype
    } else {
      dx <- "custom"
    }
    info("[getMultiOmicsProbeAnnotation] detected as:", dx)
    dtype <- rep(dx, length(probes))
  }
  table(dtype)
  dtype <- tolower(dtype)
  dtype <- ifelse(grepl("ensembl|symbol|hugo|gene|hgnc", dtype), "gx", dtype)
  dtype <- ifelse(grepl("uniprot|protein", dtype), "px", dtype)
  dtype <- ifelse(grepl("chebi|hmdb|kegg|pubchem", dtype), "mx", dtype)
  table(dtype)
  dbg("[getMultiOmicsProbeAnnotation] dtypes =", unique(dtype))

  ## populate with defaults
  symbol <- toupper(sub("^[a-zA-Z]+:", "", probes))

  annot <- list()
  if (any(dtype %in% c("gx", "px"))) {
    ii <- which(dtype %in% c("gx", "px"))
    pp <- sub("^[a-zA-Z]+:", "", probes[ii])
    aa <- getGeneAnnotation(organism, pp)
    head(aa)
    aa$data_type <- sub(":.*", "", probes[ii])
    rownames(aa) <- probes[ii]
    aa$feature <- probes[ii]
    annot <- c(annot, list(aa))
  }
  if ("mx" %in% dtype) {
    ii <- which(dtype == "mx")
    hh <- grep("mx:NA$", probes[ii])
    if (any(hh)) ii <- ii[-hh]
    pp <- sub("^[a-zA-Z]+:", "", probes[ii])
    aa <- getMetaboliteAnnotation(pp)
    head(aa)
    aa$data_type <- "mx"
    rownames(aa) <- probes[ii]
    aa$feature <- probes[ii]
    annot <- c(annot, list(aa))
  }
  if ("custom" %in% dtype) {
    ii <- which(dtype == "custom")
    pp <- sub("^[a-zA-Z]+:", "", probes[ii])
    aa <- getCustomAnnotation(pp, custom_annot = NULL)
    head(aa)
    aa$data_type <- "custom"
    rownames(aa) <- probes[ii]
    aa$feature <- probes[ii]
    annot <- c(annot, list(aa))
  }

  cols <- Reduce(intersect, lapply(annot, colnames))
  annot <- lapply(annot, function(a) a[, cols])
  annot <- do.call(rbind, annot)
  annot <- annot[match(probes, annot$feature), ]
  rownames(annot) <- probes
  head(annot)

  ## fill NA
  annot$human_ortholog[which(annot$human_ortholog == "")] <- NA
  annot$feature <- ifelse(is.na(annot$feature), probes, annot$feature)
  annot$symbol <- ifelse(is.na(annot$symbol), symbol, annot$symbol)
  annot$human_ortholog <- ifelse(is.na(annot$human_ortholog), symbol, annot$human_ortholog)
  annot$gene_name <- ifelse(is.na(annot$gene_name), probes, annot$gene_name)
  annot$data_type <- ifelse(is.na(annot$data_type), dtype, annot$data_type)

  return(annot)
}

#' Annotate multi-species probetype. Probe names *must  be prefixed with
#' data type unless classical transcriptomics/proteomics.
#'
getMultiSpeciesProbeAnnotation <- function(probes, organisms, probetype,
                                           datatype = "rna-seq") {  
  
  ## FILL ME!!!
  annot <- data.frame()
  return(annot)
}
