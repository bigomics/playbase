##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


#' Detect metabolomics ID type.
#' 
#' @export
mx.detect_probetype <- function(probes, min.match=0.05) {
  aa <- playdata::METABOLITE_ID
  probes <- setdiff(probes,c("","-",NA))
  match <- apply(aa, 2, function(s) mean(probes %in% s))
  if (max(match) < min.match) {
    ## if no match we try stripping of non-numerical prefixes
    ##probes1 <- gsub("[^0-9]", "", probes)
    probes1 <- gsub("^[a-zA-Z:]+", "", probes)
    match2 <- apply(aa, 2, function(s) mean(probes1 %in% s))
    match <- match + match2
  }
  if (max(match, na.rm = TRUE) < min.match) {
    message("[mx.detect_probetype] WARNING. could not detect probetype")
    return(NA)
  }
  names(which.max(match))
}

#' Convert IDs to CHEBI using base R functions (DEPRECATED)
#' 
#' @export
convert_probe_to_chebi <- function(probes, probe_type=NULL) {
  message("WARNING: deprecated. please use mx.convert_probe()")
  annot <- playdata::METABOLITE_ID
  valid_probe_types <- colnames(annot)
  # for id that are "", set it to na
  probes[probes == ""] <- NA

  # check that probetype is valid
  if(is.null(probe_type)) {
    probes0 <- setdiff(probes,NA)
    nmatch <- apply(annot, 2, function(a) mean(probes0 %in% a))
    if(max(nmatch,na.rm=TRUE) > 0) {
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
    matches <- match(probes, annot[,probe_type])
    chebi_ids <- annot[matches, "ChEBI"]
  } else {
    chebi_ids <- rep(NA,length(probes))
  }

  return(chebi_ids)
}

#' Convert IDs to CHEBI using internal playdata annotation table. Can
#' handle different target columns.
#' 
#' @export
mx.convert_probe <- function(probes, probe_type = NULL, target_id = "ID") {
  if (is.null(probe_type)) {
    probe_type <- mx.detect_probetype(probes)
  }
  if(is.null(probe_type) || is.na(probe_type)) {
    message("WARNING: could not determine probe_type")
    return(rep(NA,length(probes)))
  }
  # check that probetype is valid
  annot <- playdata::METABOLITE_ID
  valid_probe_types <- colnames(annot)
  probe_type <- match.arg(probe_type, valid_probe_types)

  # for id that are "", set it to na
  probes[probes == ""] <- NA
  if (probe_type == "ChEBI") {
    # keep only numbers in ids, as chebi ids are numeric
    probes <- gsub("[^0-9]", "", probes)
  }
  matches <- match(probes, annot[, probe_type])
  ids <- annot[matches, target_id]
  return(ids)
}

#' 
#'
#' @export
getMetaboliteAnnotation <- function(probes, add_id=FALSE, probe_type = NULL, 
                                    db = c("refmet","playdata","annothub") ) {
  ##add_id=TRUE;db=c("refmet","playdata","annothub") 
  
  ## strip multi-omics prefix
  has.prefix <- all(grepl("^[A-Za-z]+:",tolower(probes)))
  has.prefix
  if(has.prefix) {
    probes <- sub("^[A-Za-z]+:","",probes)
  }

  if(sum(duplicated(probes))) {
    message("WARNING duplicated probes. result will not match length")
  }
  if(any(is.na(probes))) {
    message("WARNING NA probes. result will not match length")
  }
  probes <- probes[!duplicated(probes) & !is.na(probes)]  
  probes[ probes %in% c('','-','NA')] <- NA

  colnames(playdata::METABOLITE_METADATA)  
  COLS <- c("ID", "feature","name",
    "super_class","main_class", "sub_class","formula","exactmass",
    ## "chebi_id","pubchem_id","hmdb_id","kegg_id","lipidmaps_id","refmet_id",
    # "inchi_key",
    "definition","source")

  metadata <- data.frame(matrix(NA, nrow=length(probes), ncol=length(COLS)))
  colnames(metadata) <- COLS
##  rownames(metadata) <- probes
  metadata$feature <- probes

  d=db[1]
  for(d in db) {

    no.name <- any(is.na(metadata$name))
    no.name

    ## First pass using RefMet. RefMet also handles metabolite/lipid
    ## long names, so this is convenient
    if( d == "refmet" && curl::has_internet() && no.name) {
      message("[getMetaboliteAnnotation] annotating with RefMet server...")
      ii <- which(!is.na(probes) & is.na(metadata$name) )
      if(length(ii)) {
        res <- RefMet::refmet_map_df(probes[ii])  ## request on API server
        res$definition <- '-'
        res$ID <- res$ChEBI_ID
        res$source <- "RefMet:API"
        cols <- c("ID","Input.name","Standardized.name","Super.class","Main.class",
          "Sub.class","Formula","Exact.mass",
          ## "ChEBI_ID", "PubChem_CID", "HMDB_ID", "KEGG_ID", "LM_ID", "RefMet_ID",
          #"INCHI_KEY",
          "definition","source")
        res <- res[,cols]
        colnames(res) <- COLS
        ## only fill missing entries
        jj <- which( res != '-' & !is.na(res) & is.na(metadata[ii,]), arr.ind=TRUE)
        if(length(jj)) metadata[ii,][jj] <- res[jj]
      }
    }

    ## this uses internal datatable. maybe we can get rid of it in
    ## future and use only online annotation.
    if( d == "playdata" && no.name ) {
      # get annotation for probes
      message("[getMetaboliteAnnotation] annotating with METABOLITE_METADATA...")
      id <- mx.convert_probe(probes, target="ID")
      no.name <- is.na(metadata$name) | is.na(metadata$definition) | is.na(metadata$ID)
      ii <- which(!id %in% c('',NA,'-') & no.name)      
      if(length(ii)) {
        mm <- playdata::METABOLITE_METADATA
        mm <- mm[ match(id[ii], mm$ID), ]      
        mm$feature <- probes[ii]
        ##rownames(mm) <- probes[ii]
        mm$ID <- id[ii]
        mm$source <- "ChEBI+RefMet"      
        mm <- mm[,COLS]
        ## only fill missing entries
        jj <- which( !is.na(mm) & is.na(metadata[ii,,drop=FALSE]), arr.ind=TRUE)
        if(length(jj)) metadata[ii,][jj] <- mm[jj]      
      }
    }
    
    ## this uses metaboliteIDmapping through AnnotHub. The first time
    ## retrieval of the database can take some time for caching. Needs
    ## internet for first time download.
    if(d == "annothub" && curl::has_internet() && no.name) {
      message("[getMetaboliteAnnotation] annotating with AnnotHub ...")
      adb <- AnnotationHub::AnnotationHub()
      ah <- adb[["AH91792"]]
      ii <- which(is.na(metadata$name))
      nmatch <- sapply(ah, function(a) sum(probes %in% a))
      if( max(nmatch)== 0) {
        message("[getMetaboliteAnnotation] WARNING: could not match ID column")
      } else {
        match.col <- which.max(nmatch)
        jj <- match(probes[ii], ah[[match.col]])
        res <- as.data.frame(ah[jj,]) 
        res$feature <- probes[ii]
        res$ID <- res$ChEBI
        res$Super.class <- "-"
        res$Main.class <- "-"
        res$Sub.class <- "-"
        res$Formula <- "-"
        res$Exact.mass <- "-"
        #res$INCHI_KEY <- "-"
        res$definition <- "-"            
        res$lipidmaps <- "-"
        res$refmet <- "-"
        res$source <- "R:metaboliteIDmapping"
        cols <- c("ID","feature","Name","Super.class","Main.class",
          "Sub.class","Formula","Exact.mass", 
          ## "ChEBI","CID","HMDB","KEGG","lipidmaps","refmet",
          ##"INCHI_KEY",
          "definition", "source")
        res <- res[,cols]
        colnames(res) <- COLS
        ## only fill missing entries
        jj <- which( res != '-' & !is.na(res) & is.na(metadata[ii,]), arr.ind=TRUE)
        if(length(jj)) metadata[ii,][jj] <- res[jj]      
      }
    }

  } ## for db

  rownames(metadata) <- NULL
  id <- metadata$ID
  hmdb <- mx.convert_probe(id, probe_type="ChEBI", target_id="HMDB")
  
  df <- data.frame(
    feature = probes,
    symbol = id,
    human_ortholog = id,  ## or use HMDB???
    gene_title = metadata$name,
    source = metadata$source,
    gene_name = id
  )

  kk <- setdiff(colnames(metadata),c(colnames(df),"name","ID"))
  df <- cbind( df, metadata[,kk])   

  ## add ID table 
  if(add_id) {
    annot_table <- playdata::METABOLITE_ID
    ii <- match( id, annot_table$ID )
    kk <- setdiff(colnames(annot_table), c("ID","NAME"))
    id_table <- annot_table[ii,kk]
    colnames(id_table) <- paste0(sub("_ID","",colnames(id_table)),"_ID")
    df <- cbind( df, id_table )   
  }

  probes[is.na(probes)] <- "NA"
  rownames(df) <- make_unique(as.character(probes))
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
