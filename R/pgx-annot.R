##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
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


#' Get probetype annotation for organism and datatype. For multi-omics
#' probe names must be prefixed with data type.
#' @export
getProbeAnnotation <- function(organism,
                               probes,
                               datatype,
                               meth_type = NULL, 
                               probetype = "",
                               bridge = NULL,
                               bridge_species = NULL,
                               annot_table = NULL) {

  if (is.null(datatype)) datatype <- "unknown"
  if (is.null(probetype)) probetype <- "unknown"

  unknown.organism <- (tolower(organism) %in% c("no organism", "custom", "unkown"))
  unknown.datatype <- (datatype %in% c("custom", "unkown"))
  unknown.probetype <- (probetype %in% c("custom", "unkown"))
  annot.unknown <- unknown.organism || unknown.datatype || unknown.probetype
  organism <- normalizeOrganism(organism)

  if (datatype == "methylomics") {
    c1 <- is.null(meth_type)
    c2 <- !meth_type %in% c("450K array", "EPIC array")
    if (c1 | c2) meth_type = "450K array"
    genes <- annotate_methylomics(organism, probes, meth_type = meth_type)
    return(genes)
  }

  ## clean probe names
  probes <- trimws(probes)
  probes[probes == "" | is.na(probes)] <- "NA"
  probes0 <- make_unique(probes) ## make unique but do not clean
  if (!is.null(annot_table)) {
    rownames(annot_table) <- make_unique(rownames(annot_table))
  }    

  ## only first feature (warning: can make probe not unique). This
  ## assumes semicolumn is NEVER used for single feature names. 
  probes <- sub("[;].*","",probes) ## only first probe
  probes <- make_unique(probes)
  
  genes <- NULL
  if (annot.unknown) {
    # annotation table is mandatory for 'No organism' (until server side
    # can handle missing genesets)
    info("[getProbeAnnotation] annotating with custom annotation")
    genes <- getCustomAnnotation2(probes0, annot_table)
  } else if (datatype == "metabolomics") {
    mx.check <- mx.check_mapping(
      probes,
      all.db = c("playdata", "annothub", "refmet"), check.first = TRUE
    )
    mx.check <- mean(!is.na(mx.check)) > 0.01
    if (mx.check) {
      genes <- getMetaboliteAnnotation(
        probes,
        extra_annot = TRUE,
        annot_table = annot_table
      )
    } else {
      ## Fallback on custom
      dbg("[getProbeAnnotation] WARNING: not able to map metabolomics probes")
    }
  } else if (datatype == "lipidomics") {
    genes <- getLipidAnnotation(
      probes,
      extra_annot = TRUE,
      annot_table = annot_table
    )
  } else if (datatype == "multi-omics") {
    genes <- getMultiOmicsProbeAnnotation(organism, probes)
  } else {
    if (datatype == "proteomics") {
      is.phospho <- annotate_phospho_residue(probes, detect.only = TRUE)
      genes <- getGeneAnnotation(organism = organism, probes = probes,
        bridge = bridge, bridge_species = bridge_species,
        is.phospho = is.phospho)
    } else {
      genes <- getGeneAnnotation(
        organism = organism, probes = probes,
        bridge = bridge, bridge_species = bridge_species
      )
    }
  }

  ## final fallback is genes==NULL
  if (is.null(genes)) {
    dbg("[getProbeAnnotation] WARNING: fallback to UNKNOWN probes")
    genes <- getCustomAnnotation(probes0, custom_annot = NULL)
  }

  ## if annot_table is provided we (priority) override our annotation
  ## and append any extra columns.
  if (!is.null(genes) && !is.null(annot_table)) {
    dbg("[getProbeAnnotation] merging custom annotation table")
    colnames(annot_table) <- sub("^ortholog$", "human_ortholog",
      colnames(annot_table),
      ignore.case = TRUE
    )
    colnames(annot_table) <- sub("^Symbol$|^gene$|^gene_name$", "symbol",
      colnames(annot_table),
      ignore.case = TRUE
    )
    genes <- merge_annot_table(genes, annot_table, priority = 2)
  }

  genes <- genes[match(probes, genes$feature), ]
  rownames(genes) <- probes0
  genes <- cleanupAnnotation(genes)

  return(genes)

}


#' Get feature annotation (for RNAseq and proteomics) data using annothub
#' or gprofiler.
#' 
#' @export
getGeneAnnotation <- function(
  organism,
  probes,
  is.phospho = FALSE,
  bridge = NULL,
  bridge_species = NULL,
  use.ah = NULL,
  verbose = TRUE,
  methods = c("annothub", "gprofiler")
) {
  organism <- normalizeOrganism(organism)

  ## clean up probes name. be careful not to 'overclean'.
  probes0 <- make_unique(trimws(probes))
  probes <- trimws(probes)
  probes[probes == "" | is.na(probes)] <- "NA"

  ## only first feature! (warning: can make probes non-unique)
  probes <- sub("[;].*","",probes) ## only first probe

  if (mean(grepl("[:]", probes)) > 0.98) {
    message("[getGeneAnnotation] WARNING. stripping multi-omics prefix")
    probes <- sub("^[a-zA-Z0-9]+:", "", probes)
  }

  if (is.phospho) {
    old_probes <- probes
    probes <- sub("[_].*", "", probes)
  }

  # init empty (all missings)
  annot <- data.frame(feature = probes, stringsAsFactors = FALSE)
  rownames(annot) <- make_unique(probes0)
  missing <- rep(TRUE, length(probes))
  
  for (method in methods) {
    if (any(missing)) {
      # annotation for current method
      missing_probes <- probes[which(missing)]
      missing_annot <- try(switch(method,
        "annothub" = getGeneAnnotation.ANNOTHUB(
          organism = organism,
          probes = missing_probes,
          use.ah = use.ah,
          bridge = bridge,
          bridge_species = bridge_species,
          verbose = verbose
        ),
        "gprofiler" = getGeneAnnotation.GPROFILER(
          organism = organism,
          probes = missing_probes,
          bridge = bridge,
          bridge_species = bridge_species,
          verbose = verbose
        ),
        stop("Unknown method: ", method)
      ))

      annot_ok <- !inherits(missing_annot, "try-error") &&
        !is.null(missing_annot) && nrow(missing_annot) > 0

      if (annot_ok) {
        # not all methods have the same columns
        new_cols <- setdiff(colnames(missing_annot), colnames(annot))
        if (length(new_cols) > 0) {
          for (col in new_cols) annot[[col]] <- NA
        }
        mm <- merge_annot_table(annot[missing, ], missing_annot)
        annot[missing, ] <- mm[, colnames(annot)]
        missing <- is.na(annot$symbol) | annot$symbol == ""
      }
    }
  }

  if (all(missing)) { # unsuccessful annotation
    message("[getGeneAnnotation] WARNING. all missing??? missing.ratio=", mean(missing))
    annot <- NULL
  }
  
  if (verbose > 0) {
    mean.mapped <- round(100*mean(!is.na(annot$symbol)),2)
    mean.ortho <- round(100*mean(!is.na(annot$human_ortholog)),2)
    message("[getGeneAnnotation] mapping ratio SYMBOLS  = ", mean.mapped, "%")
    message("[getGeneAnnotation] mapping ratio ORTHOLOGS  = ", mean.ortho, "%")
  }
  
  ## clean up
  if (!is.null(annot)) {
    annot <- cleanupAnnotation(annot)
  }

  ## restore original phospho probe names
  if (is.phospho && !is.null(annot)) {
    annot$feature <- old_probes
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
  bridge = NULL,  ## NULL -> auto
  bridge_species = NULL,
  verbose = TRUE
) {
  if (is.null(organism)) {
    warning("[getGeneAnnotation.ANNOTHUB] Please specify organism")
    return(NULL)
  }

  if (verbose) {
    message("[getGeneAnnotation.ANNOTHUB] Retrieving gene annotation...")
  }

  organism <- normalizeOrganism(organism)
  genes <- NULL

  ## get correct OrgDb database for this organism
  orgdb <- getOrgDb(organism, use.ah = use.ah)
  if (is.null(orgdb)) {
    message("[getGeneAnnotation.ANNOTHUB] ERROR: orgdb == NULL: ", is.null(orgdb))
    return(NULL)
  } else {
    message(paste0("[getGeneAnnotation.ANNOTHUB] OrgDb for '",organism,"' retrieved..."))
  }

  if (is.null(probes)) probes <- AnnotationDbi::keys(orgdb)

  ## Backup probes as probes0, give names of probes the original name.
  probes[is.na(probes) | probes == ""] <- "NA"
  probes0 <- make_unique(probes)
  names(probes) <- probes0

  ## clean up probe names from suffixes
  ## if(clean_probes) {
  ##   probes <- .clean_probe_names(probes)
  ## }
  
  if (is.null(probe_type)) {
    probe_type <- detect_probetype(organism, probes, orgdb = orgdb)
    if (is.null(probe_type) || is.na(probe_type)) {
      message("ERROR: could not determine probe_type.")
      message("WARNING. returning empty annotation.")
      annot <- data.frame(feature = probes, symbol = "")
      annot <- cleanupAnnotation(annot)
      annot$symbol <- NA
      annot$gene_title <- NA      
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
  cols <- c("SYMBOL", "GENENAME", "GENETYPE", "ALIAS", "MAP", "UNIPROT")
  # cols <- c("SYMBOL", "GENENAME", "GENETYPE", "MAP")
  cols <- intersect(cols, AnnotationDbi::keytypes(orgdb))

  if (organism %in% c("Mus musculus", "Rattus norvegicus")) {
    cols <- unique(c(cols, "ENTREZID"))
  }

  suppressMessages(suppressWarnings(
    annot <- AnnotationDbi_select_2pass(
      orgdb,
      keys = probes,
      columns = cols,
      keytype = probe_type
    )
  ))

  # some organisms do not provide symbol but rather gene name (e.g. yeast)
  if ("SYMBOL" %in% cols) {
    symbols <- AnnotationDbi::keys(orgdb, keytype = "SYMBOL")
  } else if ("GENENAME" %in% cols) {
    symbols <- AnnotationDbi::keys(orgdb, keytype = "GENENAME")
  }

  if (!"SYMBOL" %in% colnames(annot)) {
    annot$SYMBOL <- annot$GENENAME
    annot$GENENAME <- annot$ALIAS
  } else {
    not.symbols <- !(annot$SYMBOL %in% symbols)
    if (length(not.symbols)) annot$SYMBOL[not.symbols] <- NA
  }

  ##  annot$ALIAS <- NULL
  annot$SYMBOL[is.na(annot$SYMBOL)] <- ""

  ## Attempt to retrieve chr map via org.Mm.egCHRLOC / org.Rn.egCHRLOC.
  if (organism %in% c("Mus musculus", "Rattus norvegicus")) {
    if (organism == "Mus musculus") {
      require(org.Mm.eg.db)
      chrloc <- org.Mm.eg.db::org.Mm.egCHRLOC
    }
    if (organism == "Rattus norvegicus") {
      require(org.Rn.eg.db)
      chrloc <- org.Rn.eg.db::org.Rn.egCHRLOC
    }
    mapped_genes <- as.list(chrloc[AnnotationDbi::mappedkeys(chrloc)])
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
  dfA <- apply(annot, 2, function(a) {
    tapply(a, key, function(b) {
      paste(setdiff(unique(b), c(NA, "")), collapse = ";")
    })
  })
  annot <- data.frame(dfA, check.names = FALSE)
  annot <- annot[match(probes, rownames(annot)), ]
  rownames(annot) <- make_unique(names(probes))
  annot$PROBE <- names(probes) ## original probe names

  ## -----------------------------------------------------------------------------
  ## Second pass for missing symbols. Still trying annothub but
  ## missing symbols may map to different keytype. NEED RETHINK:
  ## DO WE REALLY NEED THIS???
  ## -----------------------------------------------------------------------------
  is.missing <- (is.na(annot$SYMBOL) | annot$SYMBOL == "")
  missing.probes <- probes[which(is.missing)] ## probes match annot!
  missing.probes <- missing.probes[!is.na(missing.probes)]
  length(missing.probes)
  if (second.pass && length(missing.probes)) {
    missing.probe_type <- try(suppressWarnings(suppressMessages(
      detect_probetype(organism, missing.probes, orgdb = orgdb)
    )), silent = TRUE)
    if (inherits(missing.probe_type, "try-error")) {
      missing.probe_type <- NULL
    }

    ## only do second try if missing.probetype is different
    if (!is.null(missing.probe_type) &&
      !is.na(missing.probe_type) &&
      missing.probe_type != probe_type &&
      missing.probe_type %in% AnnotationDbi::keytypes(orgdb)
    ) {
      missing.probes1 <- match_probe_names(missing.probes, orgdb, missing.probe_type)
      suppressMessages(suppressWarnings(
        missing.annot <- AnnotationDbi_select_2pass(
          orgdb,
          keys = missing.probes1,
          columns = cols,
          keytype = missing.probe_type
        )
      ))
      
      missing.key <- missing.annot[, missing.probe_type]
      dfA <- apply(missing.annot, 2, function(a) {
        tapply(a, missing.key, function(b) {
          paste(setdiff(unique(b), c(NA, "")), collapse = ";")
        })
      })
      missing.annot <- data.frame(dfA, check.names = FALSE)
      missing.annot <- missing.annot[match(missing.probes1, rownames(missing.annot)), , drop = FALSE]
      rownames(missing.annot) <- names(missing.probes)
      missing.annot$PROBE <- names(missing.probes)

      # some organisms do not provide SYMBOL but rather GENENAME (e.g. yeast)
      if (!"SYMBOL" %in% colnames(missing.annot)) {
        missing.annot$SYMBOL <- missing.annot$GENENAME
        missing.annot$GENENAME <- missing.annot$ALIAS
      }
      if ("SYMBOL" %in% colnames(missing.annot)) {
        not.symbols <- !(missing.annot$SYMBOL %in% symbols)
        if (length(not.symbols)) missing.annot$SYMBOL[not.symbols] <- NA
      }

      for (k in setdiff(colnames(annot), colnames(missing.annot))) {
        missing.annot[[k]] <- NA
      }
      kk <- match(colnames(annot), colnames(missing.annot))
      missing.annot <- missing.annot[, kk, drop = FALSE]
      jj <- match(missing.annot$PROBE, probes)
      ii <- which(!is.na(jj))
      annot[jj[ii], ] <- missing.annot[ii, ]
    }
  }

  ## get human ortholog using 'orthogene'
  message("[getGeneAnnotation.ANNOTHUB] getting human orthologs...")
  ortho <- getHumanOrtholog(
    organism, annot$SYMBOL,
    bridge = bridge,
    bridge_species = bridge_species,
    verbose = 0)

  annot$ORTHOGENE <- ortho$human    ## single-valued
  annot$ORTHOGENES <- ortho$humans  ## all candidates, ";"-joined

  ## Return as standardized data.frame and in the same order as input
  ## probes.
  pkgname <- orgdb$packageName
  if (length(pkgname) == 0) pkgname <- "OrgDb"
  annot$SOURCE <- pkgname[1]

  annot.cols <- c(
    "PROBE", "SYMBOL", "ORTHOGENE", "ORTHOGENES", "UNIPROT", "GENENAME",
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
    "feature", "symbol", "human_ortholog", "human_orthologs", "uniprot", "gene_title",
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


#' Annotate using gprofiler. This seems much faster than
#' AnnotHub. There are about 700 species supported. Online connection
#' to server is needed but we are already using remote AnnotHub and
#' gprofiler for ortholog matching. The advantage is that probe type
#' detection is not needed because orthogene seems to detect is
#' automatically.
#' 
#' @export
getGeneAnnotation.GPROFILER <- function(
  organism,
  probes,
  bridge = NULL,
  bridge_species = NULL,  
  verbose = TRUE
) {
  ## correct organism names different from OrgDb
  organism <- sub("Canis familiaris", "Canis lupus familiaris", organism, fixed = TRUE)

  ## map given organism to ORTHOGENE species name
  #species <- try(.getGprofilerSpecies(organism))
  species <- try(.map_gprofiler_id(organism))
  if ("try-error" %in% class(species)) {
    message("[getGeneAnnotation.GPROFILER] *WARNING* could not connect to server")
    return(NULL)
  }
  if (is.null(species)) {
    message("ERROR: unknown organism ", organism)
    return(NULL)
  }
  message("Mapping to species: ", species)
  probes[is.na(probes) | probes == ""] <- "NA"
  
  out <- try(orthogene::map_genes(
    genes = probes,
    species = species,
    run_map_species = FALSE,  ## disable map and check
    verbose = FALSE
  ), silent = TRUE)

  message("Round 1: mapped ratio: ", round(100*mean(!is.na(out$name)),2),"%")
  
  ## Retry missing probes with clean probes
  ii <- which(is.na(out$name))
  length(ii)
  if(length(ii)) {
    clean.probes <- .clean_probe_names(probes[ii], sep='.-') 
    names(clean.probes) <- probes[ii]
    out2 <- try(orthogene::map_genes(
      genes = clean.probes,
      species = species,
      run_map_species = FALSE,  ## disable map and check
      verbose = FALSE
    ), silent = TRUE)

    out2$input <- names(clean.probes)[match(out2$input,clean.probes)]
    jj <- which(!is.na(out2$name))
    if(length(jj)) {
      ii <- match(out2$input[jj], out$input)
      out[ii,] <- out2[jj,,drop=FALSE]
    }
    message("Round 2: mapped ratio: ", round(100*mean(!is.na(out$name)),2),"%")
  }
  
  df <- data.frame(
    feature = probes,
    symbol = NA,
    human_ortholog = "",
    human_orthologs = "",
    uniprot = "",
    gene_title = "",
    chr = NA,
    source = NA,
    gene_name = probes
  )
  rownames(df) <- make_unique(probes)

  if (!inherits(out, "try-error")) {
    ## map to original probe names
    out <- out[match(probes, out$input), ]

    ## add extra uniprot id
    gp.out <- try(gprofiler2::gconvert(probes, organism = species, target = "UNIPROT_GN_ACC"))
    if (!is.null(gp.out) && !inherits(gp.out, "try-error")) {
      uniprot <- tapply(gp.out$target, gp.out$input, function(x) paste(x, collapse = ";"))
      uniprot <- as.character(uniprot[match(probes, names(uniprot))])
    } else {
      uniprot <- rep(NA, length(probes))
    }

    df$symbol <- out$name
    df$gene_title <- sub(" \\[.*", "", out$description)
    ortho <- getHumanOrtholog(
      organism, out$name,
      bridge = bridge,
      bridge_species = bridge_species,
      verbose = 0
    )
    df$human_ortholog <- ortho$human    ## single-valued
    df$human_orthologs <- ortho$humans  ## all candidates, ";"-joined
    df$uniprot <- uniprot
    df$source <- "gprofiler2"
  }

  return(df)
}

.getGprofilerSpecies <- function(organism, as = c("name", "id")[1]) {
  id <- .map_gprofiler_id(organism)
  if(as == 'id') return(id)
  S <- playbase::SPECIES_TABLE
  sel <- which(S$gprofiler_id == id)
  if(length(sel)==0) return(NULL)
  species <- S[sel, "gprofiler_species"]
  species
}

.getGprofilerID <- function(organism) {
  .getGprofilerSpecies(organism, "id")
}


#' Cleanup annotation
#'
cleanupAnnotation <- function(genes) {
  if (is.null(genes)) {
    return(NULL)
  }

  ## add missing columns if needed, then reorder
  columns <- c(
    "feature", "symbol", "human_ortholog", "human_orthologs", "gene_title", ## "gene_biotype",
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

  # trim whitespace
  char.cols <- which(sapply(genes, class) == "character")
  for (k in char.cols) {
    genes[[k]] <- trimws(genes[[k]])
  }

  # rename protein-coding to protein_coding to confirm with playbase <= v1.3.2
  ## genes$gene_biotype <- sub("protein-coding", "protein_coding", genes$gene_biotype)

  # replace NA in gene_ortholog by "" to conform with old
  # pgx objects. For collapsing to symbol this is important.
  genes$human_ortholog[is.na(genes$human_ortholog)] <- ""

  # replace NA or empty symbol by "{feature}" so there is always a readable name
  ii <- which(genes$symbol %in% c(NA, "", "-"))
  genes$symbol[ii] <- paste0("{", genes$feature[ii], "}")
  genes$gene_title[ii] <- "Uknown feature"

  # if organism is human, human_ortholog should be NA (matching old
  # playbase annot). NEED RETHINK (this is not very consistent).
  if (is.null(genes$human_ortholog)) genes$human_ortholog <- NA

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
    if (is.null(custom_annot$feature)) custom_annot$feature <- rep(NA, nr)
    if (is.null(custom_annot$symbol)) custom_annot$symbol <- rep(NA, nr)
    if (is.null(custom_annot$gene_name)) custom_annot$gene_name <- rep(NA, nr)
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

  custom_annot <- custom_annot[, table_col_order, drop = FALSE]
  custom_annot <- cleanupAnnotation(custom_annot)

  return(custom_annot)
}

#' Like getCustomAnnotation() but allows custom column names mapping
#' to feature, symbol and title. Also simplified
#' implementation. Should eventually replace getCustomAnnotation().
#'
#' @export
getCustomAnnotation2 <- function(probes, custom_annot, feature.col = "feature",
                                 symbol.col = "symbol", gene_title.col = "gene_title",
                                 ortholog.col = "human_ortholog",
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
  rownames(annot) <- make_unique(probes)
  required.columns <- colnames(annot)

  # If the user has provided a custom gene table, check it and use it
  if (!is.null(custom_annot)) {
    custom_annot <- data.frame(custom_annot, check.names = FALSE)

    if (!feature.col %in% colnames(custom_annot)) {
      if (!is.null(rownames(custom_annot))) {
        custom_annot$rownames <- rownames(custom_annot)
      }
      fsum <- apply(custom_annot, 2, function(a) mean(probes %in% a, na.rm = TRUE))
      feature.col <- NULL
      if (max(fsum) > 0.9) feature.col <- names(which.max(fsum))
      if (length(feature.col) == 0) {
        custom_annot$feature <- probes
        feature.col <- "feature"
      }
    }
    if (!symbol.col %in% colnames(custom_annot)) {
      symbol.col <- head(grep("symbol|name|gene|protein|alias",
        setdiff(colnames(custom_annot), feature.col),
        ignore.case = TRUE, value = TRUE
      ), 1)
      if (length(symbol.col) == 0) symbol.col <- NA
    }
    if (!gene_title.col %in% colnames(custom_annot)) {
      gene_title.col <- head(grep("title|description|name",
        setdiff(colnames(custom_annot), c(feature.col, symbol.col)),
        ignore.case = TRUE, value = TRUE
      ), 1)
      if (length(gene_title.col) == 0) gene_title.col <- NA
    }
    if (!ortholog.col %in% colnames(custom_annot)) {
      ortholog.col <- head(grep("ortholog|human|hgnc",
        setdiff(colnames(custom_annot), c(feature.col, symbol.col)),
        ignore.case = TRUE, value = TRUE
      ), 1)
      if (length(ortholog.col) == 0) ortholog.col <- NA
    }

    features <- custom_annot[, feature.col]
    custom_annot <- custom_annot[match(probes, features), ]

    # Rename columns
    newcols <- c(
      "feature" = feature.col, "symbol" = symbol.col,
      "gene_title" = gene_title.col,
      "human_ortholog" = ortholog.col
    )
    newcols <- newcols[which(newcols != names(newcols))]
    newcols <- newcols[which(newcols %in% colnames(custom_annot))]
    if (length(newcols)) {
      custom_annot <- dplyr::rename(custom_annot, all_of(newcols))
    }

    if ("feature" %in% colnames(custom_annot)) {
      annot$feature <- custom_annot[, "feature"]
    }
    if ("symbol" %in% colnames(custom_annot)) {
      annot$symbol <- custom_annot[, "symbol"]
    }
    if ("gene_title" %in% colnames(custom_annot)) {
      annot$gene_title <- custom_annot[, "gene_title"]
    }

    ##  if (!is.null(custom_annot) && num_annot > 1 && required_in_annot) {
    # remove all NA columns, otherwise the for loop below will not work
    custom_annot <- custom_annot[, colMeans(is.na(custom_annot)) != 1, drop = FALSE]

    # identify missing columns and fill them with annot_map
    missing_cols <- setdiff(colnames(custom_annot), colnames(annot))
    missing_cols <- setdiff(missing_cols, c(NA))
    annot <- cbind(annot, custom_annot[, missing_cols, drop = FALSE])
  }

  if (!extra.columns) {
    sel <- (colnames(annot) %in% required.columns)
    annot <- annot[, sel]
  }

  message("[getCustomAnnotation2] Cleaning custom annotation table...")
  annot <- cleanupAnnotation(annot)
  return(annot)
}


## ================================================================================
## ================== GET ORTHOLOG FUNCTIONS ======================================
## ================================================================================

##BRIDGE_SPECIES = c("dmelanogaster","drerio","celegans","scerevisiae","athaliana","ecoli")
BRIDGE_SPECIES = c("dmelanogaster","drerio","celegans","scerevisiae","athaliana")

## Auto-trigger the bridge search below this direct orthology ratio. The
## poorly annotated species the bridge is meant for sit well above the old
## 0.10 gate (trout ~0.16, hamster ~0.33) and so never bridged at all.
BRIDGE_THRESHOLD = 0.5

## Number of features probed to decide whether to bridge (see getHumanOrtholog)
BRIDGE_PROBE_SIZE = 500

## Below this fraction of resolvable identifiers we consider the input
## namespace unsupported for the organism and skip the bridge fan-out:
## bridging cannot work on IDs that cannot be resolved in the first place.
ID_RECOGNITION_MIN = 0.02

#' Fraction of input identifiers that g:Profiler can resolve for this
#' organism at all, independent of orthology. This separates the two
#' reasons a lookup comes back near-empty: an unsupported or obsolete
#' identifier namespace (nothing resolves) versus genuine absence of
#' orthologs (IDs resolve, but have no human counterpart).
#'
#' @noRd
.id_recognition_ratio <- function(genes, species, nprobe = 200) {
  genes <- unique(setdiff(genes, c("", NA, "NA", "N/A", "---")))
  if (!length(genes)) return(NA)
  if (length(genes) > nprobe) {
    genes <- genes[unique(round(seq(1, length(genes), length.out = nprobe)))]
  }
  cv <- try(gprofiler2::gconvert(query = genes, organism = species,
    target = "ENSG", numeric_ns = "ENTREZGENE_ACC", filter_na = FALSE),
    silent = TRUE)
  if (inherits(cv, "try-error") || is.null(cv) || !nrow(cv)) return(NA)
  cv$target[cv$target %in% c("N/A", "NA", "")] <- NA
  mean(genes %in% cv$input[!is.na(cv$target)])
}

#' Collapse a one-to-many (input, ortholog) table to one row per input,
#' with a single deterministic ortholog (first alphabetically) and the
#' full ";"-joined candidate set.
#'
#' @noRd
.collapse_orthologs <- function(df) {
  x <- tapply(df$ortholog, df$input, .clean_ortholog_set, simplify = FALSE)
  data.frame(
    input = names(x),
    ortholog = sapply(x, function(s) if (length(s)) s[1] else NA),
    orthologs = sapply(x, function(s) if (length(s)) paste(s, collapse = ";") else NA),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

#' @title Get human ortholog from given symbols of organism by
#'   aggregation of multiple methods. This package needs internet
#'   connection.
#'
#' @export
getHumanOrtholog <- function(organism, symbols,
                             ortho.methods = c("homologene","gprofiler","babelgene",
                               "gprofiler2","uppercase"),
                             bridge = NULL,
                             bridge_species = BRIDGE_SPECIES,
                             verbose = 1) {
  
  ## try also clean symbols
  symbols[is.na(symbols)] <- "NA"
  clean.symbols <- .clean_symbols(symbols)   ## NEED RETHINK!
  names(clean.symbols) <- symbols
  names(symbols) <- symbols
  if(is.null(bridge_species)) bridge_species <- BRIDGE_SPECIES
  
  ## Try mapping with orthogene's databases
  species_id <- .getGprofilerSpecies(organism, "id")
  species_id  

  i=1
  batch_size=500
  batch_size=1500

  ##genes <- c("---", unique(symbols[!is.na(symbols)]))
  genes <- c("---", unique(c(symbols,clean.symbols)))
  genes[is.na(genes)] <- "NA"

  ## Decide ONCE whether to bridge, before the batching loop below. The
  ## threshold used to be evaluated inside each convert_orthologs() call,
  ## so the same dataset could bridge on one batch and not on the next
  ## depending purely on feature order. The probe is an evenly spaced
  ## (deterministic) sample, so the decision is reproducible across runs.
  ## NOTE: only needed when batching. A single batch is one convert_orthologs()
  ## call that decides on the same threshold itself, so probing there would
  ## just do the whole lookup twice.
  if (is.null(bridge) && length(genes) > batch_size) {
    probe <- genes
    if (length(probe) > BRIDGE_PROBE_SIZE) {
      probe <- probe[unique(round(seq(1, length(probe), length.out = BRIDGE_PROBE_SIZE)))]
    }
    probe.out <- try(convert_orthologs(
      genes = probe,
      species = species_id,
      method = ortho.methods,
      bridge = FALSE,
      verbose = 0
    ), silent = TRUE)
    if (inherits(probe.out, "data.frame") && nrow(probe.out) > 0) {
      probe.ratio <- mean(!probe.out$ortholog %in% c(NA, "", "NA", "N/A"))
      bridge <- (probe.ratio < BRIDGE_THRESHOLD)
      if (verbose > 0) {
        message("[getHumanOrtholog] probe ratio mapped = ", round(probe.ratio, 4),
          " -> bridge = ", bridge)
      }
    } else {
      bridge <- FALSE
    }
  }

  ## Batch processing for large gene lists
  if (length(genes) > batch_size) {      
    
    n_genes <- length(genes)
    n_batches <- ceiling(n_genes / batch_size)
    if (verbose > 0) message("[getHumanOrtholog] processing ", n_genes,
      " genes in ", n_batches, " batches")
    
    batch_results <- list()
    b=1
    for (b in 1:n_batches) {
      start_idx <- (b - 1) * batch_size + 1
      end_idx <- min(b * batch_size, n_genes)
      genes_batch <- genes[start_idx:end_idx]
      
      batch_out <- convert_orthologs(
        genes = genes_batch,
        species = species_id,
        method = ortho.methods,
        bridge = bridge,
        bridge_species = bridge_species,
        verbose = verbose
      )
      
      ## Only keep successful results
      if (!"try-error" %in% class(batch_out) &&
            inherits(batch_out, "data.frame") &&
            nrow(batch_out) > 0) {
        batch_results[[length(batch_results) + 1]] <- batch_out
      }
    }
    
    ## Combine all successful batch results
    if (length(batch_results) > 0) {
      ortho.out <- do.call(rbind, batch_results)
    } else {
      ortho.out <- try(stop("All batches failed"), silent = TRUE)
    }
  } else {
    ## Standard processing for smaller gene lists
    ortho.out <- convert_orthologs(
      genes = genes,
      species = species_id,
      method = ortho.methods,
      bridge = bridge,
      bridge_species = bridge_species,
      verbose = verbose
    )
  }
  
  class(ortho.out)
  results.ok <- (!"try-error" %in% class(ortho.out) &&
                   inherits(ortho.out, "data.frame") &&
                   nrow(ortho.out) > 0)
  results.ok
  if (!results.ok) {
    if (verbose > 0) message("[getHumanOrtholog] failed lookup")
    ortho.out <- NULL
  }
  
  ## return dataframe. First column organism symbols, second column
  ## human ortholog. NA if missing.
  table(ortho.out$method)
  gg <- c(symbols, clean.symbols)
  ortho.out <- ortho.out[match(gg, ortho.out$input), ]
  ortho.out$symbol <- c(names(symbols), names(clean.symbols))

  df <- data.frame(
    symbol = ortho.out$symbol,
    input = ortho.out$input,
    human = ortho.out$ortholog,
    humans = ortho.out$orthologs,
    source = ortho.out$method)

  ## collapse duplicates
  df.concat <- function(x) tapply(x, df$symbol,
    function(s) paste(unique(setdiff(s,c("",NA,"NA","N/A"))), collapse=";"))
  ## NOTE: keep simplify=FALSE, otherwise a single unique symbol collapses
  ## the result to a plain vector and df2$human becomes NULL
  df2 <- do.call(cbind, apply(df, 2, df.concat, simplify = FALSE))
  df2[which(df2=="")] <- NA
  df2 <- data.frame(df2)
  df2 <- df2[match(symbols, df2$symbol),,drop=FALSE]

  ## A symbol and its cleaned form can map to different orthologs, so the
  ## collapse above may have produced a ";"-joined 'human' as well. Keep
  ## 'human' strictly single-valued (first alphabetically) and carry the
  ## full candidate set in 'humans'.
  hh <- lapply(paste(df2$human, df2$humans, sep=";"), .clean_ortholog_set)
  df2$human <- sapply(hh, function(s) if(length(s)) s[1] else NA)
  df2$humans <- sapply(hh, function(s) if(length(s)) paste(s, collapse=";") else NA)

  ## compute mapping ratio
  mean.mapped <- round(100 * mean(!is.na(df2$human)), digits = 4)
  if (verbose > 0) message("[getHumanOrtholog] total ratio mapped  = ", mean.mapped, "%")
  if (mean.mapped==0) message("[getHumanOrtholog] WARNING: no orthologs found!")

  ## A poor result has two very different causes, and reporting only the
  ## mapping ratio hides which one it is: identifiers that cannot be
  ## resolved for this organism at all (unsupported namespace, or an
  ## obsolete annotation release) versus resolvable identifiers that
  ## simply have no human counterpart. Only checked when the result is
  ## poor, so well-mapped organisms pay nothing.
  if (mean.mapped < 100 * BRIDGE_THRESHOLD) {
    recog <- .id_recognition_ratio(genes, species_id)
    if (!is.na(recog)) {
      message("[getHumanOrtholog] recognised identifiers = ",
        round(100 * recog, 2), "%")
      if (recog < ID_RECOGNITION_MIN) {
        message("[getHumanOrtholog] WARNING: these identifiers are not ",
          "recognised for ", organism, ". Check the feature ID type.")
      } else if (recog < BRIDGE_THRESHOLD) {
        message("[getHumanOrtholog] NOTE: most features could not be resolved ",
          "for ", organism, ". They may come from an obsolete annotation ",
          "release; re-annotating against a current release should improve this.")
      }
    }
  }

  df2$input <- NULL
  return(df2)
}


#'
#'
#' 
convert_orthologs <- function(genes, species, methods = c("homologene",
  "gprofiler","babelgene", "gprofiler2", "uppercase"), bridge = NULL,
  bridge_species = BRIDGE_SPECIES, target_species = "hsapiens", verbose = 1)
{

  res <- data.frame( input = genes, ortholog = NA, orthologs = NA, method=NA)
  ##species = .getGprofilerSpecies(species, "id")

  ## try all methods
  for(m in methods) {
    ii <- which( is.na(res$ortholog))
    out <- .query_orthologs(
      genes[ii], species, method = m,
      target_species = target_species,
      verbose = verbose
    )
    if(!is.null(out) && any(!is.na(out$ortholog)) ) {
      out <- out[!is.na(out$ortholog),,drop=FALSE]
      ## gorth() is one-to-many, so an input can occur on several rows.
      ## Collapse to one row per input: a single deterministic ortholog
      ## plus the full candidate set. Assigning the long table directly
      ## would silently keep whichever candidate came last.
      out <- .collapse_orthologs(out)
      out$method <- m
      if(nrow(out)>0) {
        jj <- match( out$input, res$input)
        res[jj,] <- out[,colnames(res)]
      }
    }
  }

  orth.ratio <- mean(!res$ortholog %in% c(NA,"","NA","N/A"))
  if(verbose) {
    message("[convert_orthologs] pass 1: orth.ratio = ", orth.ratio)
  }
  if(is.null(bridge)) {
    ## Auto enable bridge lookup if poorly mapped. NOTE: when called from
    ## getHumanOrtholog() the decision is taken once, above the batching
    ## loop, and passed in explicitly: deciding here per batch makes the
    ## same dataset bridge inconsistently depending on feature order.
    bridge <- (orth.ratio < BRIDGE_THRESHOLD)
  }
  
  ## Bridging requires the input identifiers to be resolvable for this
  ## organism at all. If (almost) none are recognised, every bridge species
  ## fails in turn -- five round-trips guaranteed to return nothing. Skip
  ## the fan-out and report what is actually wrong instead.
  if(bridge && !is.null(bridge_species)) {
    recog <- .id_recognition_ratio(genes, species)
    if(!is.na(recog) && recog < ID_RECOGNITION_MIN) {
      message("WARNING: none of the features are recognised identifiers for ",
        "this organism (", species, "). Skipping bridge search: the input ",
        "looks like an unsupported namespace or an obsolete annotation release.")
      bridge <- FALSE
    }
  }

  ## try with bridge species
  if(bridge && !is.null(bridge_species)) {
    for(b in bridge_species) {
      ii <- which(is.na(res$ortholog))
      if(length(ii)) {
        B <- .query_orthologs_bridged(
          genes[ii],
          species,
          bridge_species = b,
          target_species = target_species,
          verbose = verbose
        )
        B$method <- paste0("bridged:",b)
        B <- B[,c("input","ortholog","orthologs","method")]
        jj <- which(!is.na(B$ortholog))
        if(length(jj)) {
          res[ii[jj],] <- B[jj,]
        }
      }
    }

    if(verbose) {
      orth.ratio2 <- mean(!res$ortholog %in% c(NA,"","NA","N/A"))
      message("[convert_orthologs] pass 2: orth.ratio = ", orth.ratio2)
    }
    
  }
  
  res
}

#' Single method lookup of genes/features to human orthologs. For
#' multi-methods search use .convert_orthologs()
#'
#' 
.query_orthologs <- function(genes, species, method = c("homologene",
  "gprofiler","babelgene", "gprofiler2"), target_species = "hsapiens",
  verbose = 1)
{
  out <- NULL
  if(method == "gprofiler2") {
    out <- try(gprofiler2::gorth(query = genes, source_organism = species,
      target_organism = target_species, mthreshold = Inf, filter_na = FALSE,
      numeric_ns = "ENTREZGENE_ACC"), silent = TRUE)
    if(!inherits(out,"try-error")) {
      out$ortholog_name <-sub("N/A",NA,out$ortholog_name)
      out$ortholog_name <-sub("^NA$",NA,out$ortholog_name)
      out <- out[,c("input","ortholog_name")]
      colnames(out) <- c("input","ortholog")
    }
  }

  if(method %in% c("homologene","gprofiler","babelgene")) {
    ## homologene, gprofiler, babelgene
    out <- try(orthogene::convert_orthologs(
      gene_df = genes,
      input_species = species,
      output_species = target_species,
      method = method,
      non121_strategy = "drop_both_species",
      verbose = FALSE
    ), silent = TRUE)
    if(!inherits(out,"try-error")) {
      out <- data.frame( input = out[,"input_gene"],
        ortholog = rownames(out))
    }
  }

  if(method %in% c("uppercase")) {  
    ## Map any missing symbols that look like human genes
    human.genes <- playdata::GENE_SYMBOL
    ii <- which(!is.na(genes) & toupper(genes) %in% human.genes)
    if(length(ii)) {
      out <- data.frame( input = genes[ii], ortholog = toupper(genes[ii]))
    }
  }

  if(inherits(out,"try-error")) out <- NULL
  out
}

#' Bridged lookup of genes/features to human orthologs. Uses an
#' intermediate bridge organism for mapping to human ortholog. This
#' can be useful for poorly annotated species.
#' 
#' 
.query_orthologs_bridged <- function(genes, species, bridge_species,
                                     target_species = "hsapiens",
                                     verbose = 1)
{
  genes1 <- genes
  if(any(grepl("[;]",genes))) {
    if(verbose) message("WARNING: multiple gene features detected. taking only first gene.")
    genes1 <- sub(";.*","",genes) ## strip multiple genes
  }

  if(length(bridge_species)>1) {
    if(verbose) message("WARNING: multiple bridge species detected. taking only first species.")    
    bridge_species = bridge_species[1]
  }
  if(verbose) {
    message("input species = ", species)
    message("bridge species = ", bridge_species)
    message("target species = ", target_species)  
  }
  
  out1 <- try(gprofiler2::gorth(
    query = genes1,
    source_organism = species,
    target_organism = bridge_species,
    mthreshold = Inf, filter_na = FALSE, numeric_ns = "ENTREZGENE_ACC"),
    silent = TRUE)
  if(!inherits(out1,"try-error")) {
    out1$ortholog_name <- sub("^N/A$",NA, out1$ortholog_name)
  } else {
    out1 <- NULL
  }

  bridge_genes <- out1$ortholog_name
  out2 <- NULL
  if(any(!is.na(bridge_genes))) {
    jj <- which(!is.na(bridge_genes))
    out2 <- try(gprofiler2::gorth(
      query = bridge_genes[jj],
      source_organism = bridge_species,
      target_organism = target_species,
      mthreshold = Inf, filter_na = FALSE, numeric_ns = "ENTREZGENE_ACC"),
      silent = TRUE)
    if(!inherits(out2,"try-error")) {
      out2$ortholog_name <- sub("^N/A$",NA, out2$ortholog_name)
    } else {
      out2 <- NULL
    }
  }

  if(is.null(out1) || is.null(out2)) {
    message("WARNING: bridge conversion with ",bridge_species," failed")
    out <- data.frame(input = genes, bridge = NA, ortholog = NA, orthologs = NA )
    return(out)
  }
  
  ## mappings
  gene2bridge <- tapply(out1$ortholog_name, out1$input, list)
  gene2bridge <- gene2bridge[match(genes1, names(gene2bridge))]
  bridge2human <- tapply(out2$ortholog_name, out2$input, list)  
  ##bridge2human <- bridge2human[match(genes1, names(bridge2human))]

  ## mapped features
  bridge <- lapply(gene2bridge, function(s) sort(unique(unlist(s))))
  bridge <- sapply(bridge, paste, collapse=";")
  bridge[bridge %in% c("","N/A")] <- NA
  ## the bridge is one-to-many on both hops: report a single deterministic
  ## ortholog (first alphabetically) plus the full candidate set. Keeping
  ## only the ";"-joined value breaks every exact-match consumer downstream.
  ortho <- lapply(gene2bridge, function(s) {
    .clean_ortholog_set(unlist(bridge2human[s]))
  })
  ortholog  <- sapply(ortho, function(s) if(length(s)) s[1] else NA)
  orthologs <- sapply(ortho, function(s) if(length(s)) paste(s, collapse=";") else NA)

  table(!is.na(bridge))
  table(!is.na(ortholog))

  out <- data.frame(input = genes, bridge = bridge,
    ortholog = ortholog, orthologs = orthologs )
  out
}

#' Sort/dedupe a set of candidate ortholog symbols, dropping empty and
#' N/A placeholders. Returns a sorted character vector (possibly empty).
#'
#' @noRd
.clean_ortholog_set <- function(x) {
  x <- unlist(strsplit(as.character(x), ";"))
  sort(unique(setdiff(x, c("", NA, "NA", "N/A"))))
}
  

## ================================================================================
## ========================= FUNCTIONS ============================================
## ================================================================================


#' @title Detect probe type from probe set
#' @export
detect_probetype <- function(organism, probes, orgdb = NULL,
                             nprobe = 1000, use.ah = NULL, datatype = NULL,
                             verbose = TRUE) {
  organism <- normalizeOrganism(organism)

  if (is.null(datatype) && all(grepl("[:]", probes))) {
    dbg("[detect_probetype] datatype is multi-omics?")
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
    if (verbose) {
      message("head.probes = ", paste(head(probes), collapse = " "))
      message("WARNING: Probe type not found. Valid probe types: ", paste(keytypes, collapse = " "))
    }
    # fallback before giving up; try gprofiler to convert to UNIPROT
    gp.organism <- .map_gprofiler_id(organism)    
    gp.out <- tryCatch(
    {
      gprofiler2::gconvert(probesx, organism = gp.organism, target = "UNIPROT_GN_ACC")
    },
    error = function(e) {
      return(NULL)
    }
    )
    if (!is.null(gp.out)) { key_matches["GPROFILER"] <- length(unique(gp.out$input)) / length(probesx) }
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
  is.prefixed <- mean(grepl("^[A-Za-z]+:", probes)) > 0.8
  if (is.prefixed) {
    dtype <- sub(":.*", "", probes)
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
  dtype <- sub("^tx$", "gx", dtype)
  dtype <- ifelse(grepl("ensembl|symbol|hugo|gene|hgnc", dtype), "gx", dtype)
  dtype <- ifelse(grepl("uniprot|protein", dtype), "px", dtype)
  dtype <- ifelse(grepl("chebi|hmdb|kegg|pubchem|lipid|refmet", dtype), "mx", dtype)
  table(dtype)
  dtype[!dtype %in% c("gx", "px", "mx", "lx")] <- "custom"
  dbg("[getMultiOmicsProbeAnnotation] detected datatypes = ", unique(dtype))

  ## populate with defaults
  symbol <- sub("^[a-zA-Z]+:", "", probes)
  annot <- list()
  for (dt in unique(dtype)) {
    ii <- which(dtype == dt)
    pp <- sub("^[a-zA-Z]+:", "", probes[ii])
    aa <- NULL
    if (dt %in% c("gx", "px")) {
      aa <- getGeneAnnotation(organism, pp)
    }
    if (dt %in% c("mx")) {
      aa <- getMetaboliteAnnotation(
        pp,
        db = c("lipids", "refmet", "playdata", "annothub"),
        extra_annot = TRUE, annot_table = NULL,
        prefix.symbol = FALSE
      )
    }
    if (dt %in% c("lx")) {
      aa <- getLipidAnnotation(pp, annot_table = NULL)
    }
    if (dt %in% c("custom")) {
      aa <- getCustomAnnotation(pp, custom_annot = NULL)
    }
    aa$data_type <- dt
    aa$feature <- probes[ii]
    annot[[dt]] <- aa
  }

  ## Merge all annotation tables
  names(annot)
  ## cols <- Reduce(intersect, lapply(annot, colnames))
  cols <- Reduce(union, lapply(annot, colnames))
  k <- 1
  for (k in 1:length(annot)) {
    missing.cols <- setdiff(cols, colnames(annot[[k]]))
    for (m in missing.cols) annot[[k]][[m]] <- "-"
  }
  annot <- lapply(annot, function(a) a[, cols])
  annot <- do.call(rbind, annot)
  annot <- annot[match(probes, annot$feature), ]
  rownames(annot) <- make_unique(probes)
  head(annot)

  ## fill NA
  annot$symbol[annot$symbol %in% c("-", "")] <- NA
  symbolx <- paste0("{", symbol, "}")
  annot$human_ortholog[which(annot$human_ortholog == "")] <- NA
  annot$feature <- ifelse(is.na(annot$feature), probes, annot$feature)
  annot$symbol <- ifelse(is.na(annot$symbol), symbolx, annot$symbol)
  annot$human_ortholog <- ifelse(is.na(annot$human_ortholog), symbol, annot$human_ortholog)
  annot$gene_name <- ifelse(is.na(annot$gene_name), probes, annot$gene_name)
  annot$data_type <- ifelse(is.na(annot$data_type), dtype, annot$data_type)

  return(annot)
}


#' Return n example features (symbols) for given organism
#' 
#' @export
getExampleFeatures <- function(organism, n=20, db=c("gprofiler","orgdb")) {
  if(organism %in% MAIN_ORGANISMS && length(db)==2) {
    db <- c("orgdb","gprofiler")
  }
  f <- NULL
  for(d in db) {
    if(d == "orgdb") {
      f <- try(getExampleFeatures.ORGDB(organism, n=n, protein.coding=TRUE,
        type="SYMBOL"), silent=TRUE)
    }
    if(d == "gprofiler") {
      f <- try(getExampleFeatures.GPROFILER(organism, n=n), silent=TRUE)
    }
    if(inherits(f,"try-error")) f <- NULL
    if(!is.null(f)) break
  }
  return(f)
}

getExampleFeatures.ORGDB <- function(organism, n, protein.coding=TRUE, type="SYMBOL") {
  organism <- normalizeOrganism(organism)
  orgdb <- getOrgDb(organism, use.ah = NULL)
  if (is.null(orgdb)) {
    message("[getGeneAnnotation.ANNOTHUB] ERROR: orgdb == NULL: ", is.null(orgdb))
    return(NULL)
  } else {
    message(paste0("[getGeneAnnotation.ANNOTHUB] OrgDb for '",organism,"' retrieved..."))
  }

  cols <- c(type, "SYMBOL", "ALIAS", "GENETYPE")
  cols <- intersect(cols, AnnotationDbi::keytypes(orgdb))
  if("SYMBOL" %in% cols) {
    cols <- setdiff(cols, c("ALIAS"))    
  }
  
  ez <- AnnotationDbi::keys(orgdb, keytype="ENTREZID")
  ez <- head( sample(ez), 10*n)
  
  suppressMessages(suppressWarnings(
    annot <- AnnotationDbi::select(
      orgdb,
      keys = ez,
      columns = unique(cols),
      keytype = "ENTREZID"
    )
  ))

  if("GENETYPE" %in% cols && protein.coding) {
    annot <- annot[grep("protein", annot$GENETYPE),]
  }
  symbol.name <- cols[1]
  symbols <- unique(annot[,symbol.name])
  head(symbols, n)
}

getExampleFeatures.GPROFILER <- function(organism, n) {
  species_id <- .map_gprofiler_id(organism)  
  if(is.null(species_id)) {
    message("[getExampleFeatures.GPROFILER] unknown species")
    return(NULL)
  }
  query = c("GO:0008150")  ## biological process
  out <- try(gprofiler2::gconvert(query, organism=species_id,
    mthreshold=Inf, target="ENSG"))
  sample(out$name, n)
}

