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
                               ortholog_species = "Human",
                               annot_table = NULL) {
  if (is.null(datatype)) datatype <- "unknown"
  if (is.null(probetype)) probetype <- "unknown"
  if (is.null(ortholog_species)) ortholog_species <- "Human"  

  unknown.organism <- (tolower(organism) %in% c("no organism", "custom", "unkown"))
  unknown.datatype <- (datatype %in% c("custom", "unkown"))
  unknown.probetype <- (probetype %in% c("custom", "unkown"))
  annot.unknown <- unknown.organism || unknown.datatype || unknown.probetype
  organism <- normalizeOrganism(organism)
  
  if (datatype == "methylomics") {
    c1 <- is.null(meth_type)
    c2 <- !meth_type %in% c("450K array", "EPIC array")
    if (c1 | c2) meth_type <- "450K array"
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
    genes <- getMultiOmicsProbeAnnotation(
      organism = organism,
      probes = probes,
      ortholog_species = ortholog_species      
    )
  } else {
    if (datatype == "proteomics") {
      is.phospho <- annotate_phospho_residue(probes, detect.only = TRUE)
    } else {
      is.phospho <- FALSE
    }

    genes <- getGeneAnnotation(
      organism = organism,
      probes = probes,
      ortholog_species = ortholog_species,
      is.phospho = is.phospho
    )
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
    ## colnames(annot_table) <- sub("^ortholog$", "ortholog",
    ##   colnames(annot_table),
    ##   ignore.case = TRUE
    ## )
    colnames(annot_table) <- sub("^Symbol$|^gene$|^gene_name$", "symbol",
      colnames(annot_table),
      ignore.case = TRUE
    )
    kk <- match(rownames(genes), rownames(annot_table))
    annot_table <- annot_table[kk,]
    rownames(annot_table) <- rownames(genes)
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
  use.ah = NULL,
  methods = c("annothub", "gprofiler"),
  ortholog_species = "Human",
  verbose = TRUE
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
          verbose = verbose
        ),
        "gprofiler" = getGeneAnnotation.GPROFILER(
          organism = organism,
          probes = missing_probes,
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

  ## Ortholog lookup is shared by all backends and done once, here, on
  ## the merged symbol column. Doing it inside each backend meant one
  ## remote lookup per backend, each seeing only the probes that backend
  ## happened to annotate, instead of a single lookup over the final
  ## symbols.
  if (!is.null(annot)) {
    if (verbose > 0) message(paste("[getGeneAnnotation] getting",
      ortholog_species, "orthologs..."))
    ortho <- getOrtholog(
      symbols = annot$symbol,
      organism = organism,
      target_species = ortholog_species,
      verbose = 0
    )
    annot$ortholog <- ortho$ortholog    ## single-valued
    ##annot$orthologs <- ortho$orthologs  ## all candidates, ";"-joined
    annot$ortholog_description <- ortho$description

    ## A number of downstream resources (default genesets, TileDB,
    ## CMAP, GENE_SUMMARY, Reactome/WikiPathways, GTEx tissue,
    ## cross-dataset compare, ...) are keyed on human gene symbols
    ## regardless of the user-selected ortholog_species. Keep a
    ## guaranteed-human mapping alongside the species-specific one.
    ## Skip the extra remote lookup when the two targets are the same.
    if (tolower(ortholog_species) %in% c("human","hsapiens")) {
      annot$human_ortholog <- annot$ortholog
    } else {
      if (verbose > 0) message("[getGeneAnnotation] getting Human orthologs...")
      ortho.human <- getOrtholog(
        symbols = annot$symbol,
        organism = organism,
        target_species = "Human",
        verbose = 0
      )
      annot$human_ortholog <- ortho.human$ortholog
    }
  }

  if (verbose > 0) {
    mean.mapped <- round(100*mean(!is.na(annot$symbol)),2)
    mean.ortho <- round(100*mean(!is.na(annot$ortholog)),2)
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
getGeneAnnotation.ANNOTHUB <- function(
  organism,
  probes,
  use.ah = NULL,
  probe_type = NULL,
  second.pass = TRUE,
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
    probe_type <- detect_probetype(organism, probes)
    if (is.null(probe_type) || is.na(probe_type)) {
      message("ERROR: could not determine probe_type.")
      message("WARNING. returning empty annotation.")
      annot <- data.frame(feature = probes, symbol = "")
      annot <- cleanupAnnotation(annot)
      annot$symbol <- NA
      annot$gene_title <- NA      
      annot$ortholog <- NULL
      annot$orthologs <- NULL     
      annot$ortholog <- NULL
      annot$orthologs <- NULL     
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
      detect_probetype.ANNOTHUB(organism, missing.probes, orgdb = orgdb)
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

  ## NOTE: no ortholog lookup here, and no ortholog/orthologs
  ## columns in the output. getGeneAnnotation() does the lookup once, on
  ## the merged symbol column, and adds the columns there.

  ## Return as standardized data.frame and in the same order as input
  ## probes.
  pkgname <- orgdb$packageName
  if (length(pkgname) == 0) pkgname <- "OrgDb"
  annot$SOURCE <- pkgname[1]

  annot.cols <- c(
    "PROBE", "SYMBOL", "UNIPROT", "GENENAME",
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
    "feature", "symbol", "uniprot", "gene_title",
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
getGeneAnnotation.GPROFILER <- function(
  organism,
  probes,
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
  message("Mapping to gprofiler species: ", species)
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

    ## NOTE: no ortholog lookup here, and no ortholog/orthologs
    ## columns: getGeneAnnotation() does the lookup once, on the merged
    ## symbol column, and adds the columns there.
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
    "feature", "symbol", "ortholog", "gene_title", ## "gene_biotype",
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
  genes$ortholog[is.na(genes$ortholog)] <- ""
  if (!is.null(genes$human_ortholog)) genes$human_ortholog[is.na(genes$human_ortholog)] <- ""

  # replace NA or empty symbol by "{feature}" so there is always a readable name
  ii <- which(genes$symbol %in% c(NA, "", "-"))
  genes$symbol[ii] <- paste0("{", genes$feature[ii], "}")
  genes$gene_title[ii] <- "Uknown feature"

  # if organism is human, ortholog should be NA (matching old
  # playbase annot). NEED RETHINK (this is not very consistent).
  if (is.null(genes$ortholog)) genes$ortholog <- NA

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
#' columns ortholog, gene_title, gene_biotype, chr, pos, tx_len, map, source are filled
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
    "ortholog" = "",
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
      ortholog = "",
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
                                 ortholog.col = "ortholog",
                                 extra.columns = TRUE) {
  #  feature.col='feature';symbol.col='symbol';gene_title.col='gene_title';ortholog.col='ortholog';extra.columns = TRUE

  message("[getCustomAnnotation2] Adding custom annotation table...")
  # Create custom gene table from probe names
  message("[getCustomAnnotation2] Creating annotation table from probe names...")
  annot <- data.frame(
    feature = probes,
    symbol = probes,
    gene_name = probes,
    ortholog = NA,
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
      "ortholog" = ortholog.col
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

## Below this direct orthology ratio we consider the result poor and
## report *why* it is poor (see getHumanOrtholog).
LOW_ORTHOLOGY_THRESHOLD = 0.5

## Below this fraction of resolvable identifiers we consider the input
## namespace unsupported for the organism.
ID_RECOGNITION_MIN = 0.02

## Max identifiers per ortholog query. The remote services cap the size
## of a single request, so longer queries are chunked (see .query_orthologs).
ORTHOLOG_BATCH_SIZE = 1500

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
  ortholog = sapply(x, function(s) if (length(s)) s[1] else NA)
  orthologs = sapply(x, function(s) if (length(s)) paste(s, collapse = ";") else NA)
  df <- df[match(names(x),df$input), ]
  df$ortholog  <- unname(ortholog)
  df$orthologs <- unname(orthologs)  
  return(df)
}

#' @title Get human ortholog from given symbols of organism by
#'   aggregation of multiple methods. This package needs internet
#'   connection.
#'
#' @export
getHumanOrtholog <- function(organism, symbols,
                             ortho.methods = c("homologene","gprofiler","babelgene",
                               "gprofiler2","uppercase"),
                             verbose = 1) {
  df <- getOrtholog(
    symbols = symbols,
    organism = organism,
    target_species = "Human", 
    ortho.methods = ortho.methods,
    verbose = verbose)
  df
}


#' @title Get ortholog from given symbols of organism by
#'   aggregation of multiple methods. This package needs internet
#'   connection.
#'
#' @export
getOrtholog <- function(symbols, organism, target_species, 
                        ortho.methods = c("homologene","gprofiler","babelgene",
                          "gprofiler2","uppercase"),
                        verbose = 1) {

  ## try also clean symbols
  symbols[is.na(symbols)] <- "NA"
  clean.symbols <- .clean_symbols(symbols)   ## NEED RETHINK!
  names(clean.symbols) <- symbols
  names(symbols) <- symbols
  
  ## Try mapping with orthogene's databases
  species_id <- .getGprofilerSpecies(organism, "id")
  target_id <- .getGprofilerSpecies(target_species, "id")

  ## Degrade gracefully when either species can't be resolved (unknown
  ## name, typo, organism not covered by g:Profiler) instead of letting
  ## the NULL propagate into .convert_orthologs()/.map_gprofiler_id() and
  ## throw, which would abort the whole PGX computation.
  if (is.null(species_id) || is.null(target_id)) {
    unresolved <- c(
      if (is.null(species_id)) organism,
      if (is.null(target_id)) target_species
    )
    if (verbose > 0) {
      message("[getOrtholog] could not resolve species: ",
        paste(unresolved, collapse = ", "), " -- returning empty ortholog mapping")
    }
    return(data.frame(
      symbol = symbols,
      ortholog = NA_character_,
      orthologs = NA_character_,
      description = NA_character_,
      source = NA_character_,
      row.names = NULL
    ))
  }

  ##genes <- c("---", unique(symbols[!is.na(symbols)]))
  genes <- c("---", unique(c(symbols,clean.symbols)))
  genes[is.na(genes)] <- "NA"

  ## NOTE: no batching here. Long queries are chunked inside
  ## .query_orthologs(), which is where the request-size limit actually
  ## lives, so the method cascade below sees the whole gene set at once.
  ortho.out <- .convert_orthologs(
    genes = genes,
    species = species_id,
    target_species = target_id,
    methods = ortho.methods,
    verbose = verbose
  )
  
  class(ortho.out)
  results.ok <- (!"try-error" %in% class(ortho.out) &&
                   inherits(ortho.out, "data.frame") &&
                   nrow(ortho.out) > 0)
  results.ok
  if (!results.ok) {
    if (verbose > 0) message("[getOrtholog] failed lookup")
    ortho.out <- NULL
  }
  
  ## return dataframe. First column organism symbols, second column
  ## ortholog. NA if missing.
  table(ortho.out$method)
  gg <- c(symbols, clean.symbols)
  ortho.out <- ortho.out[match(gg, ortho.out$input), ]
  ortho.out$symbol <- c(names(symbols), names(clean.symbols))

  df <- data.frame(
    symbol = ortho.out$symbol,
    input = ortho.out$input,
    ortholog = ortho.out$ortholog,
    orthologs = ortho.out$orthologs,
    description = ortho.out$description,
    source = ortho.out$method
  )

  ## collapse duplicates
  df.concat <- function(x) tapply(x, df$symbol,
    function(s) paste(unique(setdiff(s,c("",NA,"NA","N/A"))), collapse=";"))
  ## NOTE: keep simplify=FALSE, otherwise a single unique symbol collapses
  ## the result to a plain vector and df2$ortholog becomes NULL
  df2 <- do.call(cbind, apply(df, 2, df.concat, simplify = FALSE))
  df2[which(df2=="")] <- NA
  df2 <- data.frame(df2)
  df2 <- df2[match(symbols, df2$symbol),,drop=FALSE]

  ## A symbol and its cleaned form can map to different orthologs, so the
  ## collapse above may have produced a ";"-joined 'ortholog' as well. Keep
  ## 'ortholog' strictly single-valued (first alphabetically) and carry the
  ## full candidate set in 'orthologs'.
  hh <- lapply(paste(df2$ortholog, df2$orthologs, sep=";"), .clean_ortholog_set)
  df2$ortholog <- sapply(hh, function(s) if(length(s)) s[1] else NA)
  df2$orthologs <- sapply(hh, function(s) if(length(s)) paste(s, collapse=";") else NA)

  ## compute mapping ratio
  mean.mapped <- round(100 * mean(!is.na(df2$ortholog)), digits = 4)
  if (verbose > 0) message("[getOrtholog] total ratio mapped  = ", mean.mapped, "%")
  if (mean.mapped==0) message("[getOrtholog] WARNING: no orthologs found!")

  ## A poor result has two very different causes, and reporting only the
  ## mapping ratio hides which one it is: identifiers that cannot be
  ## resolved for this organism at all (unsupported namespace, or an
  ## obsolete annotation release) versus resolvable identifiers that
  ## simply have no ortholog counterpart. Only checked when the result is
  ## poor, so well-mapped organisms pay nothing.
  if (mean.mapped < 100 * LOW_ORTHOLOGY_THRESHOLD) {
    recog <- .id_recognition_ratio(genes, species_id)
    if (!is.na(recog)) {
      message("[getOrtholog] recognised identifiers = ",
        round(100 * recog, 2), "%")
      if (recog < ID_RECOGNITION_MIN) {
        message("[getOrtholog] WARNING: these identifiers are not ",
          "recognised for ", organism, ". Check the feature ID type.")
      } else if (recog < LOW_ORTHOLOGY_THRESHOLD) {
        message("[getOrtholog] NOTE: most features could not be resolved ",
          "for ", organism, ". They may come from an obsolete annotation ",
          "release; re-annotating against a current release should improve this.")
      }
    }
  }

  df2$input <- NULL
  return(df2)
}


#' Internal helper function.
#'
#' 
.convert_orthologs <- function(genes, species, methods = c("homologene",
  "gprofiler","babelgene", "gprofiler2", "uppercase"),
  target_species = "Human", verbose = 1)
{

  res <- data.frame( input = genes, ortholog = NA, orthologs = NA,
    method = NA, description = NA)
  target_species = .getGprofilerSpecies(target_species, "id")

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
        kk <- match( colnames(res), colnames(out))
        res[jj,] <- out[,kk]
      }
    }
  }

  orth.ratio <- mean(!res$ortholog %in% c(NA,"","NA","N/A"))
  if(verbose) {
    message("[convert_orthologs] orth.ratio = ", orth.ratio)
  }
  res
}

#' Single method lookup of genes/features to human orthologs. For
#' multi-methods search use .convert_orthologs()
#'
#' 
.query_orthologs <- function(genes, species, method = c("homologene",
  "gprofiler","babelgene", "gprofiler2"), target_species = "Human",
  batch_size = ORTHOLOG_BATCH_SIZE, verbose = 1)
{
  ## The remote services cap the number of identifiers per request, so a
  ## long query is chunked here, at the single point that talks to them,
  ## rather than by every caller. Chunking per method (instead of above
  ## the method cascade in .convert_orthologs) means each method still
  ## sees the full set of genes the previous methods left unresolved, so
  ## a gene is queried once per method, not once per method per chunk.
  if(length(genes) > batch_size) {
    idx <- split(seq_along(genes), ceiling(seq_along(genes) / batch_size))
    if(verbose > 0) {
      message("[.query_orthologs] ", method, ": querying ", length(genes),
        " genes in ", length(idx), " batches")
    }
    target_species = .getGprofilerSpecies(target_species, "id")
    out <- lapply(idx, function(ii) {
      .query_orthologs(genes[ii], species, method = method,
        target_species = target_species, batch_size = batch_size,
        verbose = 0)
    })
    ## a failed batch yields NULL: keep whatever the others returned
    out <- out[!sapply(out, is.null)]
    if(!length(out)) return(NULL)
    out <- do.call(rbind, out)
    rownames(out) <- NULL
    return(out)
  }

  out <- NULL
  if(method == "gprofiler2") {
    out <- try(gprofiler2::gorth(query = genes, source_organism = species,
      target_organism = target_species, mthreshold = Inf, filter_na = FALSE,
      numeric_ns = "ENTREZGENE_ACC"), silent = TRUE)
    if(!inherits(out,"try-error") && nrow(out)) {
      out$ortholog_name <- sub("N/A",NA,out$ortholog_name)
      out$ortholog_name <- sub("^NA$",NA,out$ortholog_name)
      out$description <- sub("\\[Source.*","",out$description)
      out$description <- sub("N/A",NA,out$description)            
      out <- out[,c("input","ortholog_name","description")]
      colnames(out) <- c("input","ortholog","description")
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
    if(!inherits(out,"try-error") && nrow(out)) {
      out <- data.frame(
        input = out[,"input_gene"],
        ortholog = rownames(out),
        description = NA
      )
    }
  }

  if(method %in% c("uppercase")) {  
    ## Map any missing symbols that look like human genes
    human.genes <- playdata::GENE_SYMBOL
    ii <- which(!is.na(genes) & toupper(genes) %in% human.genes)
    if(length(ii)) {
      out <- data.frame(
        input = genes[ii],
        ortholog = toupper(genes[ii]),
        description = NA        
      )
    }
  }

  if(inherits(out,"try-error")) out <- NULL
  out
}
  

## ================================================================================
## ========================= FUNCTIONS ============================================
## ================================================================================


#' Annotate multi-omics probetype. Probe names *must  be prefixed with
#' data type unless classical transcriptomics/proteomics.
#'
getMultiOmicsProbeAnnotation <- function(organism, probes, ortholog_species) {
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
      aa <- getGeneAnnotation(
        organism = organism,
        probes = pp,
        ortholog_species = ortholog_species
      )
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
  annot$ortholog[which(annot$ortholog == "")] <- NA
  annot$feature <- ifelse(is.na(annot$feature), probes, annot$feature)
  annot$symbol <- ifelse(is.na(annot$symbol), symbolx, annot$symbol)
  annot$ortholog <- ifelse(is.na(annot$ortholog), symbol, annot$ortholog)
  if (!is.null(annot$human_ortholog)) {
    annot$human_ortholog[which(annot$human_ortholog == "")] <- NA
    annot$human_ortholog <- ifelse(is.na(annot$human_ortholog), symbol, annot$human_ortholog)
  }
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

