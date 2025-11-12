##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


#' Check if RefMet server is alive
#'
mx.ping_refmet <- function() {
  out <- try(RefMet::refmet_metadata("Tyrosine"))
  if ("try-error" %in% class(out)) message("[mx.ping_refmet] WARNING! RefMet server is down")
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
  if (!refmet.alive && "refmet" %in% all.db) {
    message("WARNING: RefMet server is not alive")
    all.db <- setdiff(all.db, "refmet")
  }

  for (db in all.db) {
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
  src <- sub("^[+]", "", src) ## remove plus sign at beginning
  src[which(src == "")] <- NA
  src
}

#' Detect metabolomics ID type.
#'
#' @export
mx.detect_probetype <- function(probes, min.match = 0.2) {
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
  if (refmet.alive) {
    map <- mx.check_mapping(probes, all.db = "refmet", check.first = TRUE)
    table(map)
    if (mean(!is.na(map)) > min.match) {
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
  if (is.null(probe_type) || is.na(probe_type)) {
    message("WARNING: could not determine probe_type. Check your names.")
    return(rep(NA, length(probes)))
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
    if (!refmet.alive) {
      message("WARNING: RefMet server is down. Cannot do conversion.")
      return(rep(NA, length(probes)))
    }
    res <- RefMet::refmet_map_df(probes) ## request on API server
    probes <- res$ChEBI_ID
    probe_type <- "ChEBI"
  }
  if (is.null(probe_type) || !probe_type %in% colnames(annot)) {
    message("FATAL ERROR: probetype not in annot. probetype=", probe_type)
    return(NULL)
  }

  both.chebi <- (probe_type == "ChEBI" && target_id == "ID")
  both.hmdb <- (probe_type == "HMDB" && target_id == "HMDB")
  if (both.chebi || both.hmdb) {
    ids <- probes
  } else {
    ii <- match(probes, annot[, probe_type])
    ids <- annot[ii, target_id]
  }
  # Make sure NA are maintained (if there are NAs on annot, they get matched to random IDs sometimes) XEM
  na.probes <- is.na(probes)
  ids[na.probes] <- NA
  ids[ids == ""] <- NA
  return(ids)
}


#'
#'
#' @export
getLipidAnnotation <- function(probes,
                               extra_annot = TRUE,
                               annot_table = NULL
                               ) {
  annot <- getMetaboliteAnnotation(
    probes,
    db = c("lipids","refmet"),
    extra_annot = extra_annot,
    annot_table = annot_table,
    prefix.symbol = FALSE
  )
  return(annot)
}


#'
#'
#' @export
getMetaboliteAnnotation <- function(probes,
                                    db = c("lipids","refmet","playdata","annothub"),
                                    extra_annot = TRUE,
                                    annot_table = NULL,
                                    prefix.symbol = FALSE
                                    ) {

  ## save original probes
  orig.probes <- probes
  
  ## check for cross annotation table
  has.id <- FALSE
  if (!is.null(annot_table)) {
    ## this checks for the column in the provided annot_table that
    ## best matches any column in METABOLITE_ID and takes that column
    ## as new probe names.
    annot_table <- annot_table[match(orig.probes,rownames(annot_table)),]
    #rownames(annot_table) <- orig.probes    
    MX <- playdata::METABOLITE_ID
    id.match <- apply(annot_table, 2, function(a) {
      max(apply(MX, 2, function(m) mean(a %in% m, na.rm = TRUE)))
    })
    has.id <- max(id.match, na.rm = TRUE) > 0.8
    if (has.id) {
      ## If a good cross-lookup ID is found, take that column as new
      ## probe names.
      id.col <- which.max(id.match)
      probes <- annot_table[, id.col]
    }
  }

  ## strip multi-omics or datatype prefix
  probes <- trimws(probes)
  has.prefix <- all(grepl("^[A-Za-z]+:", probes))
  has.prefix
  if (has.prefix) {
    probes <- sub("^[A-Za-z]+:", "", probes)
  }
  probes <- trimws(probes)
  probes[probes %in% c("","-","NA",NA)] <- '-'
  ##names(probes) <- orig.probes
  names(orig.probes) <- probes

  ## check for duplicated probes
  sum(duplicated(probes))
  if (sum(duplicated(probes))) {
    message("WARNING duplicated probes. result may not match length")
    probes <- probes[!duplicated(probes)]
  }

  #colnames(playdata::METABOLITE_METADATA)
  COLS <- c(
    "feature", "symbol", "name",
    "super_class", "main_class", "sub_class", "formula", "exactmass",
    "definition", "source"
  )

  metadata <- data.frame(matrix(NA, nrow = length(probes), ncol = length(COLS)))
  colnames(metadata) <- COLS
  metadata$feature <- probes

  d <- db[1]
  for (d in db) {

    message("querying database ",d,"...")
    
    no.name <- any(is.na(metadata$name))
    no.name
    
    ## Annotate lipids using Refmet and rgoslin 
    if (d == "lipids" && curl::has_internet() && no.name) {
      ii <- which(!is.na(probes) & is.na(metadata$name) )
      aa <- mx.annotateLipids(probes[ii], two.pass=TRUE) 
      aa$definition <- '-'
      sel.col <- c("Input.name","Input.name","Standardized.name", 
        "Super.class", "Main.class", "Sub.class",
        "Formula", "Exact.mass", 'definition','Source')
      res <- aa[,sel.col]
      colnames(res) <- COLS
      
      ## only fill missing entries
      jj <- which(res != "-" & !is.na(res) & is.na(metadata[ii, ]), arr.ind = TRUE)
      if (length(jj)) {
        message("annotating ",nrow(jj)," entries with RefMet/rgoslin")
        metadata[ii, ][jj] <- res[jj]
      }
    }

    ## ----------------- RefMET -------------------------
    ## RefMet also handles metabolite/lipid long names, so this is

    ## convenient to call first.
    if( d == "refmet" && curl::has_internet() && no.name) {
      refmet.alive <- mx.ping_refmet()
      if (!refmet.alive && "refmet" %in% db) {
        message("WARNING: RefMet server is not alive. Skipping RefMet lookup")
      } else {
        message("[getMetaboliteAnnotation] annotating with RefMet server...")
        ii <- which(!is.na(probes) & is.na(metadata$name) )
        if(length(ii)) {
          id <- probes[ii]
          ## refmet wants prefix for chebi id's
          if(all(grepl("^[0-9]+",id))) id <- paste0("CHEBI:",id)
          res <- RefMet::refmet_map_df(id)  ## request on API server
          res$definition <- '-'   ## fill it??
          res$symbol <- res$ChEBI_ID
          res$source <- ifelse(res$RefMet_ID!='-', "RefMet", NA)
          cols <- c("Input.name","symbol","Standardized.name","Super.class","Main.class",
            "Sub.class","Formula","Exact.mass","definition","source")
          res <- res[,cols]
          colnames(res) <- COLS
          ## only fill missing entries
          jj <- which( res != '-' & !is.na(res) & is.na(metadata[ii,]), arr.ind=TRUE)
          if(length(jj)) {
            message("annotating ",nrow(jj)," entries with RefMet")
            metadata[ii,][jj] <- res[jj]
          }
        }
      }
    }

    ## This uses internal datatable but requires mapping to ChEBI
    ## exists. Maybe we can get rid of it in future and use only
    ## online annotation.
    if (d == "playdata" && no.name) {
      # get annotation for probes
      message("[getMetaboliteAnnotation] annotating with METABOLITE_METADATA...")
      id <- mx.convert_probe(probes, target = "ID")
      no.name <- is.na(metadata$name) | is.na(metadata$definition) | is.na(metadata$ID)
      ii <- which(!id %in% c("", NA, "-") & no.name)
      if (length(ii)) {
        M <- playdata::METABOLITE_METADATA
        M <- M[match(id[ii], M$ID), ]
        M$feature <- probes[ii]
        M$symbol <- id[ii]
        M$source <- ifelse(id[ii] %in% M$ID, "ChEBI", NA)
        M <- M[, COLS]
        ## only fill missing entries
        jj <- which(!is.na(M) & is.na(metadata[ii, , drop = FALSE]), arr.ind = TRUE)
        if (length(jj)) {
          message("annotating ",nrow(jj)," entries with METABOLITE_METADATA")
          metadata[ii, ][jj] <- M[jj]
        }
      }
    }

    ## this uses 'metaboliteIDmapping' through AnnotHub. The first
    ## time retrieval of the database can take some time for
    ## caching. Needs internet for first time download.
    if (d == "annothub" && curl::has_internet() && no.name) {
      message("[getMetaboliteAnnotation] annotating with AnnotHub ...")
      ah <- AnnotationHub::AnnotationHub()
      #AnnotationHub::query( ah, "metaboliteIDmapping")
      # AH79817 represents the original ID mapping containing 9 different ID formats
      # AH83115 mapping table which also includes common names for each compound
      # AH91792 current version of the mapping table that also accounts for tautomers
      #map <- ah[["AH91792"]]  ## 1012267 rows
      map <- ah[["AH83115"]]  ## 1116049 rows
      map <- map[setdiff(colnames(map),c("SID","CID"))] ## these are numeric
      nx <- apply(map, 1, function(m) sum(probes %in% m))
      sum(nx>0)
      if (all(nx == 0)) {
        message("[getMetaboliteAnnotation] WARNING: could not match any ID")
      } else {
        map1 <- as.data.frame(map[which(nx > 0),])
        ii <- which(is.na(metadata$name))
        res <- match.dataframe(probes[ii], map1)
        if(!is.null(res) && nrow(res)>0) {
          res$feature <- probes[ii]
          res$symbol <- probes[ii]
          res$Super.class <- "-"
          res$Main.class <- "-"
          res$Sub.class <- "-"
          res$Formula <- "-"
          res$Exact.mass <- "-"
          res$definition <- "-" ## fill it?
          res$lipidmaps <- "-"
          res$refmet <- "-"
          res$source <- ifelse( is.na(res$Name), '-', 'AH83115')
          cols <- c(
            "feature", "symbol", "Name", "Super.class", "Main.class",
            "Sub.class", "Formula", "Exact.mass", "definition", "source"
          )
          res <- res[, cols]
          colnames(res) <- COLS
          ## only fill missing entries
          jj <- which(res != "-" & !is.na(res) & is.na(metadata[ii, ]), arr.ind = TRUE)
          if (length(jj)) {
            message("annotating ",nrow(jj)," entries with AnnotHub")
            metadata[ii, ][jj] <- res[jj]
          }
        }
      }
    }
  } ## end of for db
  
  ## Harmonize lipid class names. Different databases have different
  ## namings.
  if("sub_class" %in% colnames(metadata)) {
    metadata$sub_class <- mx.harmonizeSubclassNames(metadata$sub_class)
  }
  if("main_class" %in% colnames(metadata)) {
    metadata$main_class <- mx.harmonizeMainclassNames(metadata$main_class)
  }
  if("super_class" %in% colnames(metadata)) {  
    metadata$super_class <- mx.harmonizeSuperclassNames(metadata$super_class)
  }
  
  ## This sets the default data.frame structure for metabolites. 
  rownames(metadata) <- NULL
  df <- data.frame(
    feature = probes,
    symbol = rep(NA,length(probes)),
    human_ortholog = rep(NA,length(probes)),
    gene_title = metadata$name,
    source = metadata$source,
    gene_name = metadata$feature,
    data_type = "metabolomics"
  )
  rownames(df) <- as.character(probes)
  
  if(extra_annot) {
    extra_cols <- setdiff(colnames(metadata),colnames(df))
    extra_cols <- setdiff(extra_cols, c("symbol","name","source","feature"))
    extra_cols
    if(length(extra_cols)) df <- cbind(df, metadata[,extra_cols])
  }

  ## add ID reference table. Need rethink METABOLITE_ID is not always
  ## complete.
  message("Creating ID conversion (mapping) table...")
  idtable <- mx.get_metabolite_mapping(
    probes,
    method=c("refmet","playdata","annhub")
  ) 
  df <- cbind(df, idtable)   

  ## For metabolomics we use any 'best' ID as symbol. Not only CHEBI
  ## but if CHEBI not exists any other ID. If no ID exists we fill
  ## with '{feature}'
  df$symbol <- idtable$mapping_ID
  df$human_ortholog <- df$HMDB_ID ## human
  
  ## Fill empty symbols with something. Many lipids are not mapped to
  ## our standard ChEBI id.
  fill_no_symbol = TRUE
  if(fill_no_symbol) {
    ii <- which(is.na(df$symbol))
    if(length(ii)) {
      df$symbol[ii] <- paste0("FEATURE:{",df$feature[ii],"}")
      df$name[ii] <- df$feature[ii]
      df$source[ii] <-"-"
    }
  }

  ## strip or leave the prefix of symbols??
  if(prefix.symbol == FALSE) {
    df$symbol <- sub("^[A-Za-z]+:", "", df$symbol)
  }

  ## Replace/add title and definition if we have more info in our
  ## METABOLITE_METADATA.
  message("Adding extra names and definitions...")
  M <- playdata::METABOLITE_METADATA
  chebi_id <- trimws(df$CHEBI_ID)
  M <- M[match(chebi_id,M$ID),]
  ii <- which(!is.na(M$name))
  same.name <- which(tolower(df$gene_title[ii]) == tolower(M$name[ii]))
  if(length(same.name)) {
    ## if names are same, discard first
    jj <- ii[same.name]
    df$gene_title[jj] <- '-'
  }
  df$gene_title[ii] <- ifelse(
    df$gene_title[ii] %in% c(NA,'NA','','-'),
    M$name[ii],
    paste(df$gene_title[ii], M$name[ii], sep=" | ")
  )
  ## we fill empty titles with the symbol
  jj <- which(df$gene_title %in% c(NA,'NA','','-'))  
  if(length(jj)) {
    ##df$gene_title[jj] <- '-'
    df$gene_title[jj] <- gsub("^\\{|\\}$","",df$symbol[jj]) ## remove brackets
  }
  
  ii <- which(!is.na(M$definition))
  df$definition[ii] <- ifelse(
    df$definition[ii] %in% c(NA,'NA','','-'),
    M$definition[ii],
    paste(df$definition[ii], M$definition[ii], sep=". ")
  )
  jj <- which(df$definition %in% c(NA,'NA','','-'))
  if(length(jj)) df$definition[jj] <- '-'
  
  ## Let's annotate the title with extra info
  if(extra_annot) {
    details <- paste0("(class=",df$sub_class,
      ";formula=",trimws(df$formula),";mass=",df$exactmass, ")")
    df$definition <- paste(df$definition, details)
  }
  
  ## remove unneeded columns
  df$input_ID <- NULL
  df$mapping_ID <- NULL    
  
  ## match original probe names
  df <- df[ match(names(orig.probes), df$feature), ]
  df$feature <- orig.probes
  df$gene_name <- orig.probes
  rownames(df) <- make_unique(orig.probes)

  if (extra_annot && !is.null(annot_table)) {
    extra_cols <- setdiff(colnames(annot_table),colnames(df))
    extra_cols <- setdiff( extra_cols, c("ID","symbol","name","source","feature"))
    df <- cbind(df, annot_table[,extra_cols])
  }

  ## conform empty cells
  for(k in colnames(df)) {
    jj <- which(df[[k]] %in% c(NA,'NA','','-'))
    if(length(jj)) df[[k]][jj] <- '-'
  }

  ## cleanup
  char.cols <- which(sapply(df,class) == "character")
  for(k in char.cols) df[[k]] <- trimws(df[[k]])
  
  rownames(df) <- NULL
  return(df)
}

mx.formula2chebi <- function(f, prefix=TRUE) {
  meta <- playdata::METABOLITE_METADATA
  ii <- match(f, meta[,"formula"])
  chebi <- meta[ii, "ID"]
  if(prefix) chebi <- paste("CHEBI:",chebi) ## chebi
  chebi
}

#' This uses different databases to create a ID mapping table for
#' given 'id' to other ID types (ChEBI, HMDB, LIPIDMAPS, KEGG,
#' REFMET).
#'
#' 
mx.get_metabolite_mapping <- function(met, method=c("refmet","playdata","annhub"),
                                      verbose=1) {

  COLS <- c("CHEBI","HMDB","LIPIDMAPS","KEGG","REFMET","Formula")
  map <- data.frame(matrix(NA, length(met), length(COLS)))
  colnames(map) <- COLS

  ## strip prefix. 
  orig.id <- met
  met <- sub("^[A-Za-z]+:","",met)
  prefix <- mx.prefix_id(met, return.prefix=TRUE)
  table(prefix)

  # retain prefix for chebi for refmet
  id2 <- ifelse(prefix=="CHEBI", paste0("CHEBI:",met), met)
  id2[id2 %in% c(NA,"")] <- "-"
    
  method
  for(m in method) {
    
    if(m == "refmet") {
      df <- RefMet::refmet_map_df(id2)
      sel <- c("ChEBI_ID","HMDB_ID","LM_ID","KEGG_ID","RefMet_ID","Formula")
      df <- df[,sel]
      colnames(df) <- COLS      
      ii <- which(is.na(map) & df != '-', arr.ind=TRUE)
      if(length(ii)) {
        if(verbose>0) message("mapping ", nrow(ii), " entries using RefMet")
        map[ii] <- df[ii]
      }
    }
    
    if(m == "playdata") {
      df <- playdata::METABOLITE_ID
      df$Formula <- '-'
      sel <- c("ChEBI","HMDB","LIPIDMAPS", "KEGG","REFMET","Formula")
      df <- df[,sel]
      colnames(df) <- COLS
      nmatch <- sum(met %in% unlist(df))
      if(nmatch > 0) {
        df <- match.dataframe(met, df)
        df[is.na(df)] <- '-'
        ii <- which(is.na(map) & df != '-', arr.ind=TRUE)
        if(length(ii)) {
          if(verbose>0) message("mapping ", nrow(ii), " entries using METABOLITE_ID")
          map[ii] <- df[ii]
        }
      }
    }

    if(m == "annhub") {
      suppressMessages(ah <- AnnotationHub::AnnotationHub())
      #AnnotationHub::query( ah, "metaboliteIDmapping")
      #AH79817 represents the original ID mapping containing 9 different ID formats
      #AH83115 mapping table which also includes common names for each compound
      #AH91792 current version of the mapping table that also accounts for tautomers
      #map <- ah[["AH91792"]]  ## 1012267 rows
      suppressMessages(df <- ah[["AH83115"]])  ## 1116049 rows
      df$LIPIDMAPS <- '-'
      df$Formula <- '-'      
      df$REFMET <- '-'
      sel <- c("ChEBI","HMDB","LIPIDMAPS", "KEGG","REFMET","Formula")
      df <- as.data.frame(df[,sel])
      colnames(df) <- COLS
      nmatch <- sum(met %in% unlist(df))
      if(nmatch > 0) {
        df1 <- match.dataframe(met, df, parallel=TRUE)
        df1[is.na(df1)] <- '-'
        ii <- which(is.na(map) & df1 != '-', arr.ind=TRUE)
        if(length(ii)) {
          if(verbose>0) message("mapping ", nrow(ii), " entries using AnnotationHub")
          map[ii] <- df1[ii]
        }
      }
    }
  }

  ## map contains complete mapping
  map[is.na(map)] <- '-'
  
  ## check and transfer if known ID
  ii <- which(prefix=="HMDB" & map$HMDB=='-')
  if(length(ii)) map$HMDB[ii] <- sub("HMDB:","",met[ii])
  ii <- which(prefix=="CHEBI" & map$CHEBI=='-')
  if(length(ii)) map$CHEBI[ii] <- sub("CHEBI:","",met[ii])
  ii <- which(prefix=="LIPIDMAPS" & map$LIPIDMAPS=='-')
  if(length(ii)) map$LIPIDMAPS[ii] <- sub("LIPIDMAPS:","",met[ii])
  ii <- which(prefix=="KEGG" & map$KEGG=='-')
  if(length(ii)) map$KEGG[ii] <- sub("KEGG:","",met[ii])
  ii <- which(prefix=="REFMET" & map$REFMET=='-')
  if(length(ii)) map$REFMET[ii] <- sub("REFMET:","",met[ii])

  ## determine optimal mapping ID
  mapx <- map[,order(colSums(map=='-'))]
  mapx$Formula <- NULL
  idx.argmax <- max.col(mapx != '-', ties.method='first')
  idx.type <- colnames(mapx)[idx.argmax]
  idx  <- mapx[cbind(1:nrow(mapx),idx.argmax)]
  idx.na <- rowSums(mapx != '-')==0
  idx <- paste0(idx.type,":",idx)
  idx[idx.na] <- NA
  #idx[idx.na] <- paste0("{",map$formula[idx.na],"}")
  #idx[idx.na] <- paste0("FORMULA:",map$formula[idx.na])  
  idx[idx.na] <- paste0("{",met[idx.na],"}")

  ## cleanup
  map$Formula <- NULL
  map1 <- cbind(input = met, mapping = idx, map)
  colnames(map1) <- paste0(colnames(map1),"_ID")
  
  return(map1)
}


#' Given a ChEBI id this provides information about the metabolite
#' using internal annotation databases METABOLITE_METADATA and
#' METABOLITE_ID. Also provides links to external metabolite
#' databases.
#'
#' @export
getMetaboliteInfo <- function(organism = "Human", id, info = NULL) {
  if (is.null(id) || length(id) == 0) {
    return(info)
  }
  if (is.na(id) || id == "") {
    return(info)
  }
  if (!id %in% playdata::METABOLITE_ID$ID) {
    message("[getMetaboliteInfo] unknown metabolite ID = ", id)
    return(info)
  }

  ## add metadata from tables. Append to info if provided.
  idx <- playdata::METABOLITE_ID[playdata::METABOLITE_ID$ID == id, ]
  metadata <- playdata::METABOLITE_METADATA[playdata::METABOLITE_METADATA$ID == id, ]
  if (is.null(info)) info <- list()
  info <- c(info, metadata, idx)
  info <- info[!duplicated(names(info))]
  info <- info[!sapply(info, function(s) all(is.na(s)))]

  ## substitute some names
  names(info) <- sub("definition", "summary", names(info))

  # remove summary if it is null
  if (is.null(info[["summary"]])) info[["summary"]] <- "Summary not available for this metabolite."
  if (info[["summary"]] == "null") info[["summary"]] <- "Summary not available for this metabolite."

  ## create links to external databases: these libraries are not
  ## always available for a given chebi id
  hmdb.link <- NULL
  kegg.link <- NULL
  pubchem.link <- NULL
  pathbank.link <- NULL
  refmet.link <- NULL
  lipidmaps.link <- NULL
  if (!is.null(info$HMDB)) hmdb.link <- glue::glue("<a href='https://hmdb.ca/metabolites/{info$HMDB}' target='_blank'>HMDB</a>")
  if (!is.null(info$KEGG)) kegg.link <- glue::glue("<a href='https://www.kegg.jp/dbget-bin/www_bget?{info$KEGG}' target='_blank'>KEGG</a>")
  if (!is.null(info$PubChem)) pubchem.link <- glue::glue("<a href='https://pubchem.ncbi.nlm.nih.gov/compound/{info$PubChem}' target='_blank'>PubChem</a>")
  if (!is.null(info$PATHBANK)) pathbank.link <- glue::glue("<a href='https://moldb.wishartlab.com/molecules/{info$PATHBANK}/curation.html' target='_blank'>Pathbank</a>")
  if (!is.null(info$REFMET)) refmet.link <- glue::glue("<a href='https://www.metabolomicsworkbench.org/databases/refmet/refmet_details.php?REFMET_ID={info$REFMET}' target='_blank'>RefMet</a>")
  if (!is.null(info$LIPIDMAPS)) lipidmaps.link <- glue::glue("<a href='https://dev.lipidmaps.org/databases/lmissd/{info$LIPIDMAPS}' target='_blank'>LipidMaps</a>")

  # these libraries are always available
  chebi <- info$ID
  chebi.link <- glue::glue("<a href='https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:{chebi}' target='_blank'>ChEBI</a>")
  reactome.link <- glue::glue("<a href='https://reactome.org/content/query?q=chebi%3A{chebi}' target='_blank'>Reactome</a>")

  ## collapse all links into 'databases' and remove individual entries
  info[["databases"]] <- paste(c(
    hmdb.link, chebi.link, kegg.link, reactome.link, pubchem.link,
    pathbank.link, refmet.link, lipidmaps.link
  ), collapse = ", ")
  id.cols <- colnames(playdata::METABOLITE_ID)
  id.cols <- c(id.cols, grep("_ID", names(info), value = TRUE))
  info <- info[setdiff(names(info), id.cols)]

  return(info)
}


#' Extend metabolite sets in M by adding neighbour-metabolites that
#' are connected (neighbours) according to given PPI. Only
#' metabolite-metabolite neighbours are considered.
#'
#' @export
extend_metabolite_sets <- function(M, ppi, add = TRUE, postfix = "(extended)",
                                   maxcost = 0.33) {
  ## get metabolite-metabolite edges from GRAPHITE
  ## ppi <- playdata::GRAPHITE_PPI
  ppi[, 1] <- ifelse(grepl("CHEBI", ppi[, 1]), ppi[, 1], paste0("SYMBOL:", ppi[, 1]))
  ppi[, 2] <- ifelse(grepl("CHEBI", ppi[, 2]), ppi[, 2], paste0("SYMBOL:", ppi[, 2]))
  sel <- which(grepl("CHEBI", ppi[, 1]) & grepl("CHEBI", ppi[, 2]) & ppi[, 3] <= maxcost)
  gr <- igraph::graph_from_edgelist(as.matrix(ppi[sel, 1:2]), directed = FALSE)
  MMI <- as.matrix(gr)
  table(colnames(M) %in% rownames(MMI))
  table(rownames(MMI) %in% colnames(M))

  ## build extended sparse matrix
  idx <- Matrix::which(M != 0, arr.ind = TRUE)
  x <- M[idx]
  features <- unique(c(colnames(M), rownames(MMI)))
  dims <- c(nrow(M), length(features))
  dimnames <- list(rownames(M), features)
  M2 <- Matrix::sparseMatrix(idx[, 1], idx[, 2],
    x = x,
    dims = dims, dimnames = dimnames
  )

  ## MMI neigborhood matrix for metabolites. match with extended matrix M2
  mmi0 <- rbind("na" = 0, cbind("na" = 0, MMI))
  ii <- match(colnames(M2), rownames(mmi0))
  ii[is.na(ii)] <- 1
  B <- mmi0[ii, ii]
  diag(B) <- 1
  colnames(B) <- rownames(B) <- colnames(M2)

  ## propagate neighbor
  extM <- M2 %*% B
  extM <- 1 * (extM != 0)
  extM

  ## row merge extended gene sets with original
  rownames(extM) <- paste(rownames(extM), postfix)
  if (add == TRUE) {
    extM <- Matrix::t(merge_sparse_matrix(Matrix::t(M), Matrix::t(extM)))
  }
  extM <- extM[order(rownames(extM)), ]
  return(extM)
}


#' Extend metabolite sets in M by adding neighbouring-metabolites that
#' are connected (neighbours) with members in metabolite set according
#' to given PPI. Any metabolite-metabolite or protein-metabolite
#' neighbours are considered.
#'
#' @export
extend_metabolite_sets2 <- function(M, ppi, add = TRUE, postfix = "(extended)", maxcost = 0.33) {
  ## get metabolite-metabolite edges from GRAPHITE
  ppi[, 1] <- ifelse(grepl("CHEBI", ppi[, 1]), ppi[, 1], paste0("SYMBOL:", ppi[, 1]))
  ppi[, 2] <- ifelse(grepl("CHEBI", ppi[, 2]), ppi[, 2], paste0("SYMBOL:", ppi[, 2]))
  sel <- which((grepl("CHEBI", ppi[, 1]) | grepl("CHEBI", ppi[, 2])) & ppi[, 3] <= maxcost)
  gr <- igraph::graph_from_edgelist(as.matrix(ppi[sel, 1:2]), directed = FALSE)
  PMI <- as.matrix(gr)

  ## build extended sparse matrix
  idx <- Matrix::which(M != 0, arr.ind = TRUE)
  x <- M[idx]
  features <- unique(c(colnames(M), rownames(PMI)))
  dims <- c(nrow(M), length(features))
  dimnames <- list(rownames(M), features)
  M2 <- Matrix::sparseMatrix(idx[, 1], idx[, 2],
    x = x,
    dims = dims, dimnames = dimnames
  )

  ## PMI neigborhood matrix. match with extended matrix M2
  mmi0 <- rbind("na" = 0, cbind("na" = 0, PMI))
  ii <- match(colnames(M2), rownames(mmi0))
  ii[is.na(ii)] <- 1 ## no matches will map to row/col 1 (all zero)
  B <- mmi0[ii, ii] ## B is now aligned to M2
  diag(B) <- 1
  colnames(B) <- rownames(B) <- colnames(M2)

  ## propagate all nodes to neighbors
  extM <- M2 %*% B
  extM <- 1 * (extM != 0)
  dim(extM)

  ## do not allow gene-gene neigbors. Add extended metabolites to
  ## original matrix M2
  ii <- grepl("SYMBOL", colnames(extM))
  extM[, ii] <- 0
  extM <- 1 * ((M2 + extM) != 0)

  ## row merge extended gene sets with original
  rownames(extM) <- paste(rownames(extM), postfix)
  if (add == TRUE) {
    extM <- Matrix::t(merge_sparse_matrix(Matrix::t(M2), Matrix::t(extM)))
  }
  extM <- extM[order(rownames(extM)), ]

  return(extM)
}


##======================================================================
##==================== LIPIDOMICS ======================================
##======================================================================


## detectLipidID <- function(id, min.match=0.5) {
##   if(mean(grepl("^HMDB",id,ignore.case=TRUE)) > min.match) return("HMDB")
##   if(mean(grepl("^CHEBI",id,ignore.case=TRUE)) > min.match) return("CHEBI")
##   if(mean(grepl("^RM",id,ignore.case=TRUE)) > min.match) return("RefMet")
##   if(mean(grepl("^LM",id,ignore.case=TRUE)) > min.match) return("LIPID_MAPS")
##   if(mean(grepl("^C[0-9]+",id,ignore.case=TRUE)) > min.match) return("KEGG")    
##   if(mean(grepl("^[A-Z].*[0-9]+[:][0-9]+", id)) > min.match) return("LION")
##   ## try with 
##   MX <- apply(playdata::METABOLITE_ID, 2, function(s) sub("[A-Za-z]+","",s))
##   idx <- gsub("^[A-Za-z]+[:]*","",id)
##   nmatch <- apply(MX[,-1], 2, function(x) mean(idx %in% x))
##   if(max(nmatch,na.rm=TRUE) > min.match) {
##     max.id <- names(which.max(nmatch))
##     message("matched METABOLITE_ID column: ", max.id)
##     return(max.id)
##   }  
##   return("unknown")
## }

#' Two-pass annotation of lipids using RefMet and rgoslin. This is a
#' two-pass approach where missing annotations are retried using
#' synonyms IDs of the original ID. This can improve annotation
#' coverage.
#'
#' 
mx.annotateLipids <- function(name, db=c("rgoslin","refmet","ramp"),
                              harmonize=TRUE, two.pass=TRUE) {

  ## remove multi-byte chars
  name <- iconv2utf8(name)
  name[is.na(name)] <- "NA"
  name <- as.character(name)

  ## strip away postfix from unique names. REALLY??
  orig.name <- name
  name <- sub("[.][1-9]$","",name)
  names(name) <- orig.name
    
  ## First we try annotating with their original IDs.
  aa <- mx.annotateLipids.000(name, db=db, harmonize=harmonize, add_id = TRUE)
  id.cols <- c("CHEBI_ID","HMDB_ID","LIPIDMAPS_ID","KEGG_ID","REFMET_ID")

  ## If there are still unmapped IDs, we try with their synonyms.
  ii <- which( aa$Standardized.name == '-')
  if(two.pass && length(ii)>0) {
    idx <- aa$Input.name[ii]
    mxmap <- aa[ii, id.cols]
    S <- cbind(idx, unlist(mxmap))
    colnames(S) <- c("inputID","synonymID")
    S <- S[S[,"synonymID"]!='-',]
    S <- S[which(!S[,"synonymID"] %in% S[,"inputID"]),,drop=FALSE ]
    if(nrow(S)>0) {
      synonyms <- S[,"synonymID"]
      db2 <- setdiff(db, "rgoslin")
      bb <- mx.annotateLipids.000(synonyms, db=db2, harmonize=harmonize, add_id=FALSE)

      bb <- cbind( inputID=S[,"inputID"], bb )
      bb <- bb[order(rowSums(bb == '-')),]
      bb <- bb[match(idx, bb$inputID),]
      bb <- bb[,which(colnames(bb) %in% colnames(aa))]

      kk <- match(colnames(bb),colnames(aa))      
      jj <- which(aa[ii,kk]=='-' & bb!='-', arr.ind=TRUE)

      if(nrow(jj) > 0) {
        message("filled ", nrow(jj), " entries with synonym IDs")
        aa[ii,kk][jj] <- bb[jj]
      }
    }
  }
  
  return(aa)
}

#' Annotate lipids using RefMet and rgoslin. One-pass method. Only
#' annotate with given ID.
#' 
mx.annotateLipids.000 <- function(name, db = c("rgoslin","refmet","ramp"),
                                  add_id=TRUE, harmonize=TRUE, verbose=1) {
  
  COLS <- c("Input.name","Standardized.name","Formula","Exact.mass",
    "Super.class","Main.class","Sub.class","Source")
  df <- as.data.frame(matrix('-',nrow=length(name),ncol=length(COLS)))
  colnames(df) <- COLS
  df$Input.name <- name
  
  for(d in db) {   

    missing <- (df$Formula %in% c(NA,"","-") & !name %in% c(NA,"","-"))  
      
    ## 1: We use first RefMet because it handles non-standard names
    ## better.
    refmet.online <- mx.ping_refmet()
    refmet.online
    if(d == "refmet" && refmet.online && any(missing)) {
        
      id2 <- name[which(missing)]
      if(all(grepl("^[0-9]+",name))) id2 <- paste0("CHEBI:",name)  
      aa <- RefMet::refmet_map_df(id2)  ## request on API server
      aa$Source <- ifelse(aa$Formula=='-', '-', "RefMet")
      aa$Input.name <- id2
      sel.cols <- c("Input.name","Standardized.name","Formula","Exact.mass",
                    "Super.class","Main.class","Sub.class","Source")
      aa <- aa[,sel.cols]
      colnames(aa) <- COLS
      jj <- which(!aa$Formula %in% c(NA,"","-"))
      if(verbose>0) message("matching ",length(jj), " missing IDs with RefMet...")
      if(length(jj)) {
          ii <- which(missing)[jj]
          df[ii,] <- aa[jj,]
      }
    }
      
    ## 2: Any missing annotation will be attempted with rgoslin. Only
    ## works if ID is common lipid name.
    if(d == "rgoslin" && any(missing)) {

      id.missing <- name[which(missing)]
      aa <- suppressMessages(suppressWarnings(rgoslin::parseLipidNames(id.missing)))
      nparse <- sum(!is.na(aa$Normalized.Name))

      if(nparse>0) {
        sel.cols <- c("Original.Name","Normalized.Name","Sum.Formula","Mass",
                      "Lipid.Maps.Category","Lipid.Maps.Main.Class",
                      "Functional.Class.Abbr","Source")
        for(k in setdiff(sel.cols,colnames(aa))) aa[[k]] <- "-"
        aa$Functional.Class.Abbr <- gsub("\\]|\\[","",aa$Functional.Class.Abbr)
        aa$Source <- "rgoslin"
        aa <- aa[,sel.cols]
        colnames(aa) <- COLS
        jj <- which(!aa$Formula %in% c(NA,"","-"))
        if(verbose>0) message("matching ",length(jj), " missing IDs with rgoslin...")

        if(length(jj)) {
            ii <- which(missing)[jj]
            df[ii,] <- aa[jj,]
        }

      }
    }
      
    ## 3: Attempt with RAMP
    if(d == "ramp" && any(missing)) {
      id.missing <- name[which(missing)]
      aa <- ramp.annotate_metabolites(id.missing)
      if(!is.null(aa) && nrow(aa)) {
        aa$Source <- "RaMP"
        sel.cols <- c("Input.name","Standardized.name","Formula","Exact.mass",
                      "Super.class","Main.class","Sub.class","Source")
        for(k in setdiff(sel.cols,colnames(aa))) aa[[k]] <- "-"
        aa <- aa[,sel.cols]
        colnames(aa) <- colnames(df)    
        jj <- which(!aa$Formula %in% c(NA,"","-"))
        if(verbose>0) message("matching ",length(jj), " missing IDs with RaMP...")
        if(length(jj)>0) {
            ii <- which(missing)[jj]
            df[ii,] <- aa[jj,]
        }
      }
    }      
  }

  missing <- (df$Formula %in% c(NA,"","-"))
  if(verbose>0) message("annotated ",sum(!missing),"/",nrow(df)," features")

  ## Harmonize class names. Class names are quite different between
  ## MetRef and rgoslin annotations.
  if(harmonize) {
    df$Sub.class <- mx.harmonizeSubclassNames(df$Sub.class)
    df$Main.class <- mx.harmonizeMainclassNames(df$Main.class)
    df$Super.class <- mx.harmonizeSuperclassNames(df$Super.class)
  }

  ## Retrieve cross-reference mapping to other IDs
  if(add_id) {
    xref <- mx.get_metabolite_mapping(
        name, method=c("refmet","playdata","annhub"),
        verbose=verbose-1)
    xref$input_ID <- NULL
    df <- cbind(df, xref)
  }
  
  return(df)
}

#' Translates rgoslin lipid classes to RefMet lipid classes
#' 
mx.harmonizeSuperclassNames <- function(name) {
  ## substitute abbreviated superclass
  superclass <- c(
    "FA"="Fatty Acyls",
    "GL"="Glycerolipids",
    "GP"="Glycerophospholipids",
    "SP"="Sphingolipids",
    "ST"="Sterol Lipids",
    "PR"="PrenolLipids",
    "SL"="Saccharolipids",
    "PK"="Polyketides")
  name <- ifelse( name %in% names(superclass),
    superclass[name], name )
  return(name)
}

#' Translates rgoslin lipid classes to RefMet lipid classes
#' 
mx.harmonizeMainclassNames <- function(name) {
  mainclass <- c(
    "TG"="Triradylglycerols",
    "DG"="Diradylglycerols",
    "LPC"="Glycerophosphocholines",
    "Cer"="Ceramides",
    "SE 27:1"="Sterol esters",
    "CE"="Cholesteryl esters",
    "FA"="Fatty acids",
    "LPE"="Glycerophosphoethanolamines",
    "LPI"="Glycerophosphoinositols",    
    "SM"="Phosphosphingolipids",
    "CAR"="Fatty esters",
    "PC"="Glycerophosphocholines",
    "PE"="Glycerophosphoethanolamines",
    "PI"="Glycerophosphoinositols",
    "PG"="Glycerophosphoglycerols",
    "PA"="Glycerophosphates",
    "PS"="Glycerophosphoserines"
  )
  name <- ifelse( name %in% names(mainclass),
    mainclass[name], name )
  return(name)
}

#' Translates rgoslin lipid classes to RefMet lipid classes
#' 
mx.harmonizeSubclassNames <- function(name) {
  subclass <- c(
    "Glycerophosphates"="PA",
    "Glycerophosphocholines"="PC",
    "Glycerophosphoethanolamines"="PE",
    "Glycerophosphoglycerols"="PG",
    "Glycerophosphoserines"="PS",
    "Chol. esters"="CE",
    "Cholesteryl esters"="CE",
    "Diradylglycerols"="DG",
    "Triacylglycerols"="TG",
    "DAG"="DG", #??
    "TAG"="TG",
    "Acyl carnitines"="CAR"
  )
  name <- ifelse( name %in% names(subclass),
    subclass[name], name )
  return(name)
}

mx.prefix_id <- function(id, uppercase=TRUE, return.prefix=FALSE) {

  ## clean up, remove prefix
  # id <- iconv2ascii(id)
  id <- iconv2utf8(id)
  id <- sub("[A-Za-z]+:","",id, useBytes = TRUE)

  prefix <- rep("symbol",length(id))  ## default
  prefix[grep("^[0-9]+$",id)] <- "chebi"
  prefix[grep("^HMDB[0-9]{7}",id)] <- "hmdb"
  prefix[grep("^LM[A-Z]{2}[0-9]{8}",id)] <- "lipidmaps"
  prefix[grep("^C[0-9]{5}",id)] <- "kegg"
  prefix[grep("^RM[0-9]{7}",id)] <- "refmet"    
  #prefix[grep("^ENS[A-Z]+[0-9]+$",id)] <- "ensembl"
  prefix[grep("[}{]",id)] <- "generic"
  prefix[grep("[A-Za-z1-9]+ [0-9]+:[0-9]+[;]*[0-9]*$",id)] <- "shorthand"
  prefix[grep("[A-Za-z1-9]+[ \\(].*[0-9]+:[0-9]+.*",id)] <- "shorthand"    

  table(prefix)  
  if(uppercase) prefix <- toupper(prefix)
  if(return.prefix) return(prefix)
  idx <- paste0(prefix,":",id)  
  idx
}

ramp.annotate_metabolites <- function(id) {

  idx <- iconv2utf8(id)
  idx <- mx.prefix_id(idx, uppercase=FALSE)
  idx <- sub("^symbol:","gene_symbol:",idx,ignore.case=TRUE)
  idx <- sub("^lipidmaps:","LIPIDMAPS:",idx,ignore.case=TRUE)  
  idx <- sub("^gene:","gene_symbol:",idx,ignore.case=TRUE)  
  
  pfx <- mx.prefix_id(id, uppercase=FALSE, return.prefix=TRUE)  
  table(pfx)
  
  ## check connection
  res <- try(suppressMessages(suppressWarnings(
    RaMP::getChemicalProperties("chebi:12345")
  )))
  if("try-error" %in% class(res)) {
    message("WARNING:: failed to connect github server")
    return(NULL)
  }
  
  ## RaMP does not like quotes inside IDs
  idx <- gsub("['`\"]","",idx)
  idx <- iconv2utf8(idx)
  
  ## chem properties
  suppressMessages(chemprop <- RaMP::getChemicalProperties(idx))
  chemdata <- chemprop$chem_props
  colnames(chemdata)
  sel.chem <- c("chem_source_id","common_name","mol_formula","monoisotop_mass")
  ii <- match(idx, chemdata$chem_source_id)
  annot <- data.frame(input_id=id, chemdata[ii,sel.chem])
  colnames(annot) <- c("Input.name","Matched.name","Standardized.name","Formula","Exact.mass")
  rownames(annot) <- NULL
  head(annot)
  
  ## add lipid class info
  chem <- suppressMessages(RaMP::getChemClass(mets = idx, inferIdMapping = TRUE))
  metclass <- chem$met_classes
  dim(metclass)
  if(nrow(metclass)) {
    table(metclass$class_level_name)
    CFclasses <- c("ClassyFire_super_class","ClassyFire_class","ClassyFire_sub_class")
    LMclasses <- c("LipidMaps_category","LipidMaps_main_class","LipidMaps_sub_class")
    nCF <- sum(metclass$class_level_name %in% CFclasses)
    nLM <- sum(metclass$class_level_name %in% LMclasses)
    ## take 'best' classification. we cannot use both.
    if(nLM > nCF) {
      classes <- LMclasses
    } else {
      classes <- CFclasses
    }
    classdata <- list()
    a=classes[1]
    for(a in classes) {
      ii <- which(metclass$class_level_name == a)
      aa <- metclass[ii,]
      aa <- aa[ match(idx, aa$sourceId), "class_name"]
      classdata[[a]] <- aa
    }
    classdata <- data.frame(classdata)
    colnames(classdata) <- c("Super.class","Main.class","Sub.class")

    ## If a name is found, fill if still missing
    common_name <- metclass$common_names[match(idx, metclass$sourceId)]
    ii <- which(is.na(annot$Matched.name) & !common_name %in% c(NA,'','-'))
    if(length(ii)) annot$Matched.name[ii] <- common_name[ii]

    annot <- cbind(annot, classdata)
  }
  
  rownames(annot) <- NULL
  return(annot)
}


#' This function maps rownames of matrix M to corresponding symbol as
#' given in annotation dataframe 'annot'. This is for example useful
#' if M is a sparse geneset matrix with non-standard identifiers as
#' rownames and we want to convert to specific symbols as defined by
#' the annotation dataframe 'annot'. 
#'
#' @no-export
map2symbol <- function(M, annot, target.symbol=NULL) {
  if(is.null(target.symbol)) target.symbol <- annot$symbol
  id.cols <- c("feature","symbol","gene_name",grep("_ID$",colnames(annot),value=TRUE))
  id.cols
  idmat <- as.matrix(annot[,id.cols])
  symbol.gmt <- apply(idmat, 1, function(m) unique(m), simplify=FALSE)
  symbol.gmt <- lapply(symbol.gmt, function(m) setdiff(m,c("-","",NA)))
  symbol.gmt <- lapply(symbol.gmt, function(m) sub("^[A-Za-z]+:","",m))
  symbol.gmt <- lapply(symbol.gmt, function(m) unique(m))
  names(symbol.gmt) <- target.symbol
  
  symbol.map <- lapply(names(symbol.gmt), function(i) cbind(symbol.gmt[[i]],i) )
  symbol.map <- symbol.map[sapply(symbol.map,ncol)==2]
  symbol.map <- do.call(rbind, symbol.map)
  
  id <- rownames(M)
  idx <- sub("^[A-Za-z]+:","",id)
  names(idx) <- id
  jj <- match(idx, symbol.map[,1])
  ii <- which(!is.na(jj))
  if(length(ii)==0) {
    message("WARNING: no synonyms match input probes")
    return(NULL)
  }
  M1 <- M[ii,]
  matched.symbol <- symbol.map[jj[ii],2]
  rownames(M1) <- matched.symbol
  M1 <- M1[, Matrix::colSums(M1!=0)>0, drop=FALSE]
  return(M1)
}


#' This adds some 'gene sets' based on the species class of the lipids
#' according to the annotation table.
#'
#' @export
mx.create_metabolite_sets <- function(annot, gmin=0, metmin=5,
                                   as_matrix=FALSE) {

  sel <- which(!annot$symbol %in% c(NA, '-', ''))
  annot <- annot[sel,,drop=FALSE]
  if(nrow(annot) == 0) {
    message("WARNING. empty annotation")
    return(list())
  }
  
  gmt <- list()
  
  orig.symbol <- annot$symbol
  ANNOT_SYMBOL <- mx.prefix_id(annot$symbol)
  names(orig.symbol) <- ANNOT_SYMBOL

  ## This uses the genesets/pathways from playdata::MSETxMETABOLITE
  ## these were compiled using the graphite R package.
  if(TRUE) {
    M <- Matrix::t(playdata::MSETxMETABOLITE)
    M1 <- map2symbol(M, annot, target.symbol=ANNOT_SYMBOL) 
    if(!is.null(M1) && nrow(M1) && ncol(M1) ) {
      gmt1 <- mat2gmt(M1)
      names(gmt1) <- sub("^METABOLITE:","METABOLITE_PATHWAY:",names(gmt1))
      message("[create_metabolite_sets] adding ",length(gmt1)," sets from MSETxMETABOLITE")
      gmt <- c(gmt, gmt1)
    }
  }
    
  ## Create class sets from metabolite class annotation.
  if (any(grepl("_class",colnames(annot)))) {
    classes <- grep("_class",colnames(annot),value=TRUE)
    classes
    gmt1 <- list()
    for(k in classes) {
      cv <- annot[,k]
      cv[cv=='-'] <- NA
      if(!all(is.na(cv))) {
        gmt2 <- tapply(ANNOT_SYMBOL, cv, c)
        names(gmt2) <- paste0(names(gmt2), " (",k,")")
        gmt1 <- c(gmt1, gmt2)
      }
    }
    if(length(gmt1)>0) {
      message("[create_metabolite_sets] adding ",length(gmt1),
        " sets from metabolite classes")
      names(gmt1) <- paste0("METABOLITE_CLASS:", names(gmt1))
      gmt <- c(gmt, gmt1)
    }
  }

  ## Create class/pathway/ontology sets using RaMP. Some sets are
  ## probably duplicated from MSETxMETABOLITE as graphite also
  ## includes reactome/wiki metabolic pathways.
  gmt2 <- ramp.get_metabolite_sets(
    id = ANNOT_SYMBOL,
    db = c("pathway","onto","class"),
    gmin = gmin,
    metmin = metmin,
    as.orig = TRUE
  ) 
  if(length(gmt2)) {
    message("[create_metabolite_sets] adding ",length(gmt2)," sets from RaMP database")
    gmt <- c(gmt, gmt2)
  }

  ## order on size, take out duplicated pathways
  gmt <- gmt[order(-sapply(gmt,length))]
  gmt <- gmt[!duplicated(names(gmt))]
  pw.id <- gsub(".*\\[|\\]","",names(gmt))
  pw.id[!grepl("PATHWAY",names(gmt))] <- NA
  sel <- which( is.na(pw.id) | !duplicated(pw.id))
  gmt <- gmt[sel]
  
  ## Filter on size
  G <- gmt2mat(gmt, bg = ANNOT_SYMBOL)
  table(sub(":.*","",rownames(G)))
  i0 <- grepl("SYMBOL|GENE|NAME|^px|^gx",rownames(G))
  i1 <- grepl("CHEBI|LIPID|HMDB|KEGG|REFMET|SHORT|^mx",rownames(G))
  n0 <- Matrix::colSums(G[i0,,drop=FALSE])
  n1 <- Matrix::colSums(G[i1,,drop=FALSE])  
  sel <- which(n0 >= gmin & n1 >= metmin)
  G <- G[,sel,drop=FALSE]

  # restore original symbols
  rownames(G) <- orig.symbol[rownames(G)]
  ##rownames(G) <- as.character(rownames(G))

  if(!as_matrix) {
    gmt <- mat2gmt(G)
    return(gmt)
  }
  
  # normalize columns (required for some methods downstream)log2foldchange
  G <- normalize_cols(G) ## ??? 
  return(G)
}

ramp.get_metabolite_sets <- function(id, db=c("pathway","onto","class"),
                                     gmin=0, metmin=5, as.orig=TRUE) {

  ## This adds automatically a prefix because RaMP needs identifiers
  ## with prefix.
  idx <- mx.prefix_id(id, uppercase=FALSE)  
  idx <- sub("^symbol:","gene_symbol:",idx,ignore.case=TRUE)
  idx <- sub("^lipidmaps","LIPIDMAPS",idx,ignore.case=TRUE)  
  idx <- sub("^gene:","gene_symbol:",idx,ignore.case=TRUE)  
  names(id) <- idx

  ## check if there are no genes but threshold is set
  ngenes <- sum(grepl("ensembl|gene|symbol|uniprot",idx))  
  if(ngenes==0 && gmin>0) {
    message("WARNING: gmin>0 but no genes in id list. Setting gmin=0")
    gmin <- 0
  }

  ## Query RaMP for patways
  gmt.pathway <- list()
  if("pathway" %in% db) {
    suppressMessages(pw <- RaMP::getPathwayFromAnalyte(
      analytes = idx,
      findSynonym = FALSE,
      namesOrIds = "ids",
      includeRaMPids = FALSE,
      includeSmpdb = FALSE,
      minPathwaySize = (metmin+gmin),
      maxPathwaySize = 400,
    ))
    if(!is.null(pw) && nrow(pw)>0) {
      pathway_id <- paste0(pw$pathwayName, " [", pw$pathwayId, "]")
      gmt.pathway <- tapply(pw$inputId, pathway_id, c) 
      names(gmt.pathway) <- paste0("METABOLITE_PATHWAY:",names(gmt.pathway))
    }
  }

  ## Query RaMP for Ontology
  gmt.onto <- list()
  if("onto" %in% db) {
    suppressMessages(onto <- RaMP::getOntoFromMeta(mets = idx))
    if(!is.null(onto) && nrow(onto)) {
      onto_id <- paste0(onto$Ontology, " [", onto$HMDBOntologyType, "]")
      gmt.onto <- tapply(onto$sourceId, onto_id, function(s) unique(s)) 
      names(gmt.onto) <- paste0("METABOLITE_ONTOLOGY:",names(gmt.onto))
    }
  }
  
  ## Query RaMP chemical/lipid class
  gmt.class <- list()
  if("class" %in% db) {
    suppressMessages(
      chem <- RaMP::getChemClass(
        mets = idx,
        inferIdMapping = TRUE
      )
    )
    metclass <- chem$met_classes
    dim(metclass)
    if(!is.null(metclass) && nrow(metclass)) {
      levels <- unique(metclass$class_level_name)
      levels
      gmt.class <- list()
      for(lev in levels) {
        ii <- which( metclass$class_level_name == lev )
        levclass <- metclass[ii,]
        levgmt <- tapply( levclass$sourceId, levclass$class_name, function(s) unique(s))
        lev <- sub("ClassyFire_","CF ",lev)
        lev <- sub("LipidMaps_","LM ",lev)
        names(levgmt) <- paste0(names(levgmt), " (",lev,")")
        gmt.class <- c(gmt.class, levgmt)
      }
      names(gmt.class) <- paste0("METABOLITE_CHEMCLASS:",names(gmt.class))
    }
  }
  
  ## merge all gmt
  gmt <- c(gmt.pathway, gmt.onto, gmt.class)
  ##hist(sapply(gmt, length), breaks=100)
  
  ## filter on size
  filter <- (gmin > 0 || metmin > 0)
  if(filter) {
    n0 <- sapply( gmt, function(s) sum(grepl("ensembl|gene_symbol",s)))
    n1 <- sapply( gmt, function(s) sum(grepl("chebi|hmdb|kegg|LIPIDMAPS|refmet",s)))  
    gmt <- gmt[which(n0 >= gmin & n1 >= metmin)]
  }

  ## convert to original id
  if(as.orig) {
    gmt <- lapply(gmt, function(s) as.vector(id[s]))
  }
  
  return(gmt)
}
