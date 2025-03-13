##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##



#' Check if RefMet server is alive
#' 
mx.ping_refmet <- function() {
  out <- try(RefMet::refmet_metadata("Tyrosine"))
  if("try-error" %in% class(out)) message("[mx.ping_refmet] WARNING! RefMet server is down")
  !("try-error" %in% class(out))
}

#' Check metabolite mapping in databases. For each probe return
#' matched database or NA if not found.
#'
#' @export
mx.check_mapping <- function(probes,
                             all.db = c("playdata", "annothub", "refmet"),
                             check.first = TRUE) {
  ## this goes through all databases and checks if there is a match
  src <- rep("", length(probes))
  refmet.alive <- mx.ping_refmet()
  if(!refmet.alive && "refmet" %in% all.db) {
    message("WARNING: RefMet server is not alive")
    all.db <- setdiff(all.db,"refmet")
  }
  
  for(db in all.db) {
    message("trying db: ", db)
    if (check.first) {
      jj <- which(src == "")
    } else {
      jj <- 1:length(probes)
    }
    if (length(jj)) {
      annot <- getMetaboliteAnnotation(probes[jj], db = db)
      ii <- which(!is.na(annot$source))
      jj <- jj[ii]
      src[jj] <- paste0(src[jj], "+", annot$source[ii])
    }
  }
  src
  src <- sub("^[+]","",src) ## remove plus sign at beginning
  src[which(src=='')] <- NA
  src
}

#' Detect metabolomics ID type.
#'
#' @export

mx.detect_probetype <- function(probes, min.match=0.2) {
  aa <- playdata::METABOLITE_ID
  probes <- setdiff(probes, c("", "-", NA))
  probes <- gsub("^[a-zA-Z]+:|[_.-].*", "", probes) ## strip any pre and postfix
  probes <- gsub("chebi:|chebi", "", probes, ignore.case = TRUE)
  match <- apply(aa, 2, function(s) mean(probes %in% s))

  if (max(match, na.rm = TRUE) > min.match) {
    ptype <- names(which.max(match))
    return(ptype)
  }

  ## otherwise check RefMet
  refmet.alive <- mx.ping_refmet()
  if(refmet.alive) {
    map <- mx.check_mapping(probes, all.db="refmet",check.first = TRUE)
    table(map)
    if(mean(!is.na(map)) > min.match) {
      ptype <- names(which.max(table(map[!is.na(map)])))
      return(ptype)
    }
  }
  message("[mx.detect_probetype] WARNING. could not detect probetype")
  return(NA)
}

#' Convert IDs to CHEBI using base R functions (DEPRECATED)
#'
#' @export
convert_probe_to_chebi <- function(probes, probe_type = NULL) {
  message("WARNING: deprecated. please use mx.convert_probe()")
  annot <- playdata::METABOLITE_ID
  valid_probe_types <- colnames(annot)
  # for id that are "", set it to na
  probes[probes == ""] <- NA
  probes <- gsub(".*:|[_.-].*", "", probes) ## strip any pre/postfix

  # check that probetype is valid
  if (is.null(probe_type)) {
    probes0 <- setdiff(probes, NA)
    nmatch <- apply(annot, 2, function(a) mean(probes0 %in% a))
    if (max(nmatch, na.rm = TRUE) > 0) {
      probe_type <- names(which.max(nmatch))
    }
    probe_type
  }
  probe_type <- match.arg(probe_type, valid_probe_types)

  if (probe_type == "ChEBI") {
    # keep only numbers in ids, as chebi ids are numeric
    chebi_ids <- gsub("[^0-9]", "", probes)
    # check that ChEBI ids are in the dictionary (for better results always use ChEBI own ids (ID col) not MetaboAnalyst derived)
    chebi_ids <- ifelse(chebi_ids %in% annot$ID, chebi_ids, NA)
  } else if (probe_type %in% valid_probe_types) {
    matches <- match(probes, annot[, probe_type])
    chebi_ids <- annot[matches, "ChEBI"]
  } else {
    chebi_ids <- rep(NA, length(probes))
  }

  return(chebi_ids)
}

#' Convert IDs to CHEBI using internal playdata annotation table. Can
#' handle different target columns.
#'
#' @export
mx.convert_probe <- function(probes, probe_type = NULL, target_id = "ID") {
  # for id that are "", set it to na
  probes[probes == ""] <- NA
  ## probes <- sub("[ _.-].*","",probes)  ## strip any postfix??? Not good for lipid names
  if (!is.null(probe_type) && probe_type == "ChEBI") {
    # keep only numbers in ids, as chebi ids are numeric
    probes <- gsub("[^0-9]", "", probes)
  }
  if (is.null(probe_type)) {
    probe_type <- mx.detect_probetype(probes)
  }
  if(is.null(probe_type) || is.na(probe_type)) {
    message("WARNING: could not determine probe_type. Check your names.")
    return(rep(NA,length(probes)))
  }
  # check that probetype is valid
  annot <- playdata::METABOLITE_ID
  if (probe_type == "ChEBI") {
    # keep only numbers in ids, as chebi ids are numeric
    probes <- gsub("[^0-9]", "", probes)
  }
  if (probe_type == "RefMet") {
    ## if probe_type is 'long name' we convert first to ChEBI using
    ## the RefMet server, if the server is alive.
    refmet.alive <- mx.ping_refmet()
    if(!refmet.alive) {
      message("WARNING: RefMet server is down. Cannot do conversion.")
      return(rep(NA,length(probes)))      
    }
    res <- RefMet::refmet_map_df(probes)  ## request on API server
    probes <- res$ChEBI_ID
    probe_type <- "ChEBI"
  }
  if(!probe_type %in% colnames(annot)) {
    message("FATAL ERROR: probetype not in annot. probetype=",probe_type)
    return(NULL)
  }
  ii <- match(probes, annot[, probe_type])
  ids <- annot[ii, target_id]
  ## jj <- which(is.na(ii))
  return(ids)
}

#'
#'
#' @export
getMetaboliteAnnotation <- function(probes, add_id=FALSE, 
                                    db = c("refmet","playdata","annothub"),
                                    annot_table = NULL ) {
  ##add_id=TRUE;db=c("refmet","playdata","annothub") 
  orig.probes <- probes

  ## strip multi-omics prefix
  has.prefix <- all(grepl("^[A-Za-z]+:", tolower(probes)))
  has.prefix
  if (has.prefix) {
    probes <- sub("^[A-Za-z]+:", "", probes)
  }
  names(orig.probes) <- probes

  ## check for cross annotation table
  has.id <- FALSE
  if (!is.null(annot_table)) {
    MX <- playdata::METABOLITE_ID
    mx.ids <- toupper(setdiff(colnames(MX), "ID"))
    mx.ids <- c(mx.ids, paste0(mx.ids, "_ID"))
    annot.cols <- toupper(colnames(annot_table))
    has.id <- !is.null(annot_table) && any(annot.cols %in% mx.ids)
    if (has.id) {
      id.match <- apply(annot_table, 2, function(a) {
        max(apply(MX, 2, function(m) sum(a %in% m)))
      })
      id.col <- which.max(id.match)
      dbg("[getProbeAnnotation] cross lookup using id.col =", colnames(annot_table)[id.col])
      probes <- annot_table[, id.col]
    }
  }
  has.id

  sum(duplicated(probes))
  if (sum(duplicated(probes))) {
    message("WARNING duplicated probes. result will not match length")
    probes <- probes[!duplicated(probes)]
  }
  probes[probes %in% c("", "-", "NA")] <- NA
  if (any(is.na(probes))) {
    message("WARNING NA probes. result will not match length")
    probes <- probes[!is.na(probes)]
  }

  colnames(playdata::METABOLITE_METADATA)
  COLS <- c(
    "ID", "feature", "name",
    "super_class", "main_class", "sub_class", "formula", "exactmass",
    "definition", "source"
  )

  metadata <- data.frame(matrix(NA, nrow = length(probes), ncol = length(COLS)))
  colnames(metadata) <- COLS
  ##  rownames(metadata) <- probes
  metadata$feature <- probes

  d <- db[1]
  for (d in db) {
    no.name <- any(is.na(metadata$name))
    no.name

    ## RefMet also handles metabolite/lipid long names, so this is
    ## convenient
    if( d == "refmet" && curl::has_internet() && no.name) {
      refmet.alive <- mx.ping_refmet()
      if(!refmet.alive && "refmet" %in% db) {
        message("WARNING: RefMet server is not alive. Skipping RefMet lookup")
      } else {
        message("[getMetaboliteAnnotation] annotating with RefMet server...")
        ii <- which(!is.na(probes) & is.na(metadata$name) )
        if(length(ii)) {
          probes1 <- probes[ii]
          res <- RefMet::refmet_map_df(probes1)  ## request on API server
          res$definition <- '-'   ## fill it??
          res$ID <- res$ChEBI_ID
          res$source <- ifelse(res$RefMet_ID!='-', "RefMet", NA)
          cols <- c("ID","Input.name","Standardized.name","Super.class","Main.class",
            "Sub.class","Formula","Exact.mass",
            "definition","source")
          res <- res[,cols]
          colnames(res) <- COLS
          ## only fill missing entries
          jj <- which( res != '-' & !is.na(res) & is.na(metadata[ii,]), arr.ind=TRUE)
          if(length(jj)) metadata[ii,][jj] <- res[jj]
        }
      }
    }

    ## this uses internal datatable. maybe we can get rid of it in
    ## future and use only online annotation.
    if (d == "playdata" && no.name) {
      # get annotation for probes
      message("[getMetaboliteAnnotation] annotating with METABOLITE_METADATA...")
      id <- mx.convert_probe(probes, target = "ID")
      no.name <- is.na(metadata$name) | is.na(metadata$definition) | is.na(metadata$ID)
      ii <- which(!id %in% c("", NA, "-") & no.name)
      if (length(ii)) {
        mm <- playdata::METABOLITE_METADATA
        mm <- mm[match(id[ii], mm$ID), ]
        mm$feature <- probes[ii]
        ## rownames(mm) <- probes[ii]
        mm$ID <- id[ii]
        mm$source <- ifelse(id[ii] %in% mm$ID, "ChEBI", NA)
        mm <- mm[, COLS]
        ## only fill missing entries
        jj <- which(!is.na(mm) & is.na(metadata[ii, , drop = FALSE]), arr.ind = TRUE)
        if (length(jj)) metadata[ii, ][jj] <- mm[jj]
      }
    }

    ## this uses metaboliteIDmapping through AnnotHub. The first time
    ## retrieval of the database can take some time for caching. Needs
    ## internet for first time download.
    if (d == "annothub" && curl::has_internet() && no.name) {
      message("[getMetaboliteAnnotation] annotating with AnnotHub ...")
      adb <- AnnotationHub::AnnotationHub()
      ah <- adb[["AH91792"]]
      ii <- which(is.na(metadata$name))
      nmatch <- sapply(ah, function(a) sum(probes %in% a))
      if (max(nmatch) == 0) {
        message("[getMetaboliteAnnotation] WARNING: could not match ID column")
      } else {
        match.col <- which.max(nmatch)
        jj <- match(probes[ii], ah[[match.col]])
        res <- as.data.frame(ah[jj, ])
        res$feature <- probes[ii]
        res$ID <- res$ChEBI
        res$Super.class <- "-"
        res$Main.class <- "-"
        res$Sub.class <- "-"
        res$Formula <- "-"
        res$Exact.mass <- "-"
        # res$INCHI_KEY <- "-"
        res$definition <- "-" ## fill it?
        res$lipidmaps <- "-"
        res$refmet <- "-"
        res$source <- ifelse(!is.na(jj), "metIDmap", NA)
        cols <- c(
          "ID", "feature", "Name", "Super.class", "Main.class",
          "Sub.class", "Formula", "Exact.mass", "definition", "source"
        )
        res <- res[, cols]
        colnames(res) <- COLS
        ## only fill missing entries
        jj <- which(res != "-" & !is.na(res) & is.na(metadata[ii, ]), arr.ind = TRUE)
        if (length(jj)) metadata[ii, ][jj] <- res[jj]
      }
    }
  } ## for db

  rownames(metadata) <- NULL
  df <- data.frame(
    feature = probes,
    symbol = metadata$ID,
    human_ortholog = metadata$ID,
    gene_title = metadata$name,
    source = metadata$source,
    gene_name = metadata$ID
  )
  rownames(df) <- as.character(probes)

  ## add ID table 
  if(add_id) {
    id_table <- playdata::METABOLITE_ID
    ii <- match( df$symbol, id_table$ID )
    kk <- setdiff(colnames(id_table), c("ID","NAME"))
    id_table <- id_table[ii,kk]
    colnames(id_table) <- paste0(sub("_ID","",colnames(id_table)),"_ID")
    if("ChEBI_ID" %in% colnames(id_table)) {
      ## fill missing ChEBI entries with internal CHEBI id (in symbol column)
      id_table$ChEBI_ID <- ifelse(is.na(id_table$ChEBI_ID), df$symbol, id_table$ChEBI_ID)
    }
    df <- cbind( df, id_table )   
  }

  if (has.id) {
    ii <- match(orig.probes, names(probes))
    df <- df[ii, ]
    df$feature <- orig.probes
    df$gene_name <- orig.probes
    rownames(df) <- orig.probes
  }

  return(df)
}


#' Given a ChEBI id this provides information about the metabolite
#' using internal annotation databases METABOLITE_METADATA and
#' METABOLITE_ID. Also provides links to external metabolite
#' databases.
#'
#' @export
getMetaboliteInfo <- function(organism = "Human", id) {
  if (is.null(id) || length(id) == 0) {
    return(NULL)
  }
  if (is.na(id) || id == "") {
    return(NULL)
  }

  metabolite_metadata <- playdata::METABOLITE_METADATA
  annot <- playdata::METABOLITE_ID

  inf <- list()
  inf[["name"]] <- metabolite_metadata[metabolite_metadata$ID == id, "name"]
  inf[["summary"]] <- metabolite_metadata[metabolite_metadata$ID == id, "definition"]
  inf[["organism"]] <- organism

  # remove summary if it is null
  if (is.null(inf[["summary"]])) inf[["summary"]] <- "Summary not available for this metabolite."
  if (inf[["summary"]] == "null") inf[["summary"]] <- "Summary not available for this metabolite."

  # get annotation for a given metabolite id
  annotation <- annot[annot$ID == id, ]

  # remove NA columns from annotation
  annotation <- annotation[, colSums(is.na(annotation)) < nrow(annotation)]
  cols <- colnames(annotation)[-1] ## exclude chebi IDS as we already have it

  ## get info from different environments
  res <- lapply(cols, function(k) {
    link <- NULL
    matched_id <- annotation[annotation$ID == id, k]
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
  chebi <- annotation[annotation$ID == id, "ChEBI"]
  chebi.link <- glue::glue("<a href='https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:{chebi}' target='_blank'>ChEBI</a>")
  reactome.link <- glue::glue("<a href='https://reactome.org/content/query?q=chebi%3A{chebi}' target='_blank'>Reactome</a>")

  inf[["databases"]] <- paste(c(hmdb.link, kegg.link, reactome.link, pubchem.link, chebi.link), collapse = ", ")

  return(inf)
}


#' @export
extend_metabolite_sets <- function(M) {
  ppi <- playdata::GRAPHITE_PPI
  ppi[, 1] <- ifelse(grepl("CHEBI", ppi[, 1]), ppi[, 1], paste0("SYMBOL:", ppi[, 1]))
  ppi[, 2] <- ifelse(grepl("CHEBI", ppi[, 2]), ppi[, 2], paste0("SYMBOL:", ppi[, 2]))
  sel <- which(grepl("CHEBI", ppi[, 1]) & grepl("CHEBI", ppi[, 2]) & ppi[, 3] <= 0.33)
  gr <- igraph::graph_from_edgelist(as.matrix(ppi[sel, 1:2]))
  ppi_mat <- as.matrix(gr)
  table(colnames(M) %in% rownames(ppi_mat))
  table(rownames(ppi_mat) %in% colnames(M))

  ## build neigborhood matrix for metabolites
  idx <- Matrix::which(M != 0, arr.ind = TRUE)
  x <- M[idx]
  features <- unique(c(colnames(M), rownames(ppi_mat)))
  dims <- c(nrow(M), length(features))
  dimnames <- list(rownames(M), features)
  M2 <- Matrix::sparseMatrix(idx[, 1], idx[, 2],
    x = x,
    dims = dims, dimnames = dimnames
  )

  ppi0 <- rbind("na" = 0, cbind("na" = 0, ppi_mat))
  ii <- match(colnames(M2), rownames(ppi0))
  ii[is.na(ii)] <- 1
  B <- ppi0[ii, ii]
  diag(B) <- 1
  colnames(B) <- rownames(B) <- colnames(M2)
  ## propagate neighbor
  extM <- M2 %*% B
  extM <- 1 * (extM != 0)
  extM

  ## row merge with original
  rownames(extM) <- paste(rownames(extM), "(extended)")
  rn <- unique(c(rownames(M), rownames(extM)))
  extM <- Matrix::t(merge_sparse_matrix(Matrix::t(M), Matrix::t(extM)))
  extM <- extM[match(rn, rownames(extM)), ]

  return(extM)
}
