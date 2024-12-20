##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


#' detect metabolomics ID type
#' @export
mx.detect_probetype <- function(probes, min.match=0.05) {
  aa <- playdata::METABOLITE_ANNOTATION
  probes <- setdiff(probes,c("","-",NA))
  match <- apply(aa, 2, function(s) mean(probes %in% s))
  if (max(match) < min.match) {
    ## if no match we try stripping of non-numerical prefixes
    ##probes1 <- gsub("[^0-9]", "", probes)
    probes1 <- gsub("^[a-zA-Z :_]+", "", probes)
    match2 <- apply(aa, 2, function(s) mean(probes1 %in% s))
    match <- match + match2
  }
  if (max(match, na.rm = TRUE) < min.match) {
    message("[mx.detect_probetype] WARNING. could not detect probetype")
    return(NULL)
  }
  names(which.max(match))
}

#' convert IDs to CHEBI using base R functions (DEPRECATED)
#' 
#' @export
convert_probe_to_chebi <- function(probes, probe_type=NULL) {

  annot <- playdata::METABOLITE_ANNOTATION
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

#' convert IDs to CHEBI using internal playdata annotation table
#' 
#' @export
mx.convert_probe <- function(probes, probe_type = NULL, target_id = "ChEBI") {
  if (is.null(probe_type)) {
    probe_type <- mx.detect_probetype(probes)
  }
  if(is.null(probe_type)) {
    message("WARNING: could not determine probe_type")
    return(NULL)
  }
  # check that probetype is valid
  annot <- playdata::METABOLITE_ANNOTATION
  valid_probe_types <- colnames(annot)
  probe_type <- match.arg(probe_type, valid_probe_types)

  # for id that are "", set it to na
  probes[probes == ""] <- NA
  if (probe_type == "ChEBI") {
    # keep only numbers in ids, as chebi ids are numeric
    probes <- gsub("[^0-9]", "", probes)
  }
  matches <- match(probes, annot[, probe_type])
  chebi_ids <- annot[matches, target_id]
  return(chebi_ids)
}

#' @export
getMetaboliteAnnotation <- function(probes, probe_type = NULL, id=FALSE,
                                    db = c("refmet","playdata") ) {
  ##probes = c("Ceramide","d18:0/26:2","citrate","Trilauroyl-glycerol","PE(aa-40:4)","HMDB00201","C00099","Octenoyl-L-carnitine")

  has.prefix <- all(grepl("[a-z]+:",tolower(probes)))
  has.prefix
  if(has.prefix) {
    probes <- sub(".*:","",probes)
  }

  if(sum(duplicated(probes))) {
    message("WARNING duplicated probes. result will not match length")
  }
  if(any(is.na(probes))) {
    message("WARNING NA probes. result will not match length")
  }
  probes <- probes[!duplicated(probes) & !is.na(probes)]
  
  colnames(playdata::METABOLITE_METADATA)
  
  COLS <- c("ID", "feature","name","super_class","main_class",
            "sub_class","formula","exactmass",
#            "chebi_id", "pubchem_cid",
#            "hmdb_id", "lipidmaps_id","kegg_id", "refmet_id",
            "inchi_key", "definition","source")

  metadata <- data.frame(matrix(NA, nrow=length(probes), ncol=length(COLS)))
  colnames(metadata) <- COLS
##  rownames(metadata) <- probes
  metadata$feature <- probes

  ## First pass using RefMet. RefMet also handles metabolite/lipid
  ## long names, so this is convenient.
  if( curl::has_internet() && "refmet" %in% db) {
    message("annotating with RefMet API server...")
    ##res <- mx.getRefMet(probes)
    ii <- which(!probes %in% c('',NA,'-'))
    res <- RefMet::refmet_map_df(probes[ii])  ## request on API server
    res$definition <- '-'
    res$ID <- res$ChEBI_ID
    res$source <- "RefMet:API"
    cols <- c("ID","Input.name","Standardized.name","Super.class","Main.class",
              "Sub.class","Formula","Exact.mass",
              ##"ChEBI_ID", "PubChem_CID", "HMDB_ID", "LM_ID","KEGG_ID", "RefMet_ID",
              "INCHI_KEY", "definition","source")
    res <- res[,cols]
    colnames(res) <- COLS
    ## only fill missing entries
    jj <- which( res != '-' & !is.na(res) & is.na(metadata[ii,]), arr.ind=TRUE)
    metadata[ii,][jj] <- res[jj]
  }

  if("playdata" %in% db ) {
    # get annotation for probes
    ##probes_chebi <- convert_probe_to_chebi(probes, probe_type)
    message("annotating with METABOLITE_METADATA...")
    chebi_id <- mx.convert_probe(probes)
    if(!all(is.na(chebi_id))) {
      mm <- playdata::METABOLITE_METADATA
      ii <- which(!chebi_id %in% c('',NA,'-'))
      match_ids <- match(chebi_id[ii], mm$ID)
      mm <- mm[match_ids, ]      
      colnames(mm) <- sub("DEFINITION","definition",colnames(mm))
      colnames(mm) <- sub("NAME","name",colnames(mm))      
      mm$feature <- probes[ii]
      ##rownames(mm) <- probes[ii]
      mm$ID <- chebi_id[ii]
      mm$source <- "ChEBI+RefMet"      
      mm <- mm[,COLS]
      ## only fill missing entries
      jj <- which( !is.na(mm) & is.na(metadata[ii,,drop=FALSE]), arr.ind=TRUE)
      if(length(jj)) metadata[ii,][jj] <- mm[jj]      
    }
  }

  rownames(metadata) <- NULL
  chebi_id <- metadata$ID
  
  df <- data.frame(
    feature = probes,
    symbol = chebi_id,
    human_ortholog = chebi_id,
    gene_title = metadata$name,
    source = metadata$source,
    gene_name = chebi_id
  )

  kk <- setdiff(colnames(metadata),c(colnames(df),"name","ID"))
  df <- cbind( df, metadata[,kk])   

  if(id) {
    annot_table <- playdata::METABOLITE_ANNOTATION
    ii <- match( chebi_id, annot_table$ID )
    kk <- setdiff(colnames(annot_table), c("ID","NAME"))
    df <- cbind( df, annot_table[ii,kk] )   
  }

  rownames(df) <- probes 
  return(df)
}

