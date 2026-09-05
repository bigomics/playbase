##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

## ================================================================================
## ====================== GO FUNCTIONS ============================================
## ================================================================================

#' Get GO gene sets for organism directly from
#' AnnotationHub/OrgDB. Restrict to genes as background.
#'
#' @export
getOrganismGO <- function(organism, features=NULL, minsize=3, batch_size=2000,
                          db = c("annothub","gprofiler"), include_iea=TRUE ) {

  gmt1=gmt2=NULL
  if(is.null(db) || "annothub" %in% db) {
    gmt1 <- getOrganismGO.ANNOTHUB(organism, use.ah=NULL, orgdb=NULL)
    if(is.null(db) && !is.null(gmt1)) db <- "annothub"
    message(paste("Got",length(gmt1),"GO terms from AnnotHub"))
  }

  if(is.null(features) && !is.null(gmt1)) {
    features <- unique(unlist(gmt1))
  }
  
  if(!is.null(features) && "gprofiler" %in% db) {    
    gmt2 <- getOrganismGO.GPROFILER(organism, features,
      batch_size=batch_size, include_iea=include_iea)
    message(paste("Got",length(gmt2),"GO terms from Gprofiler"))
  }

  ## order by largest.
  gmt <- c(gmt1, gmt2)
  if(length(gmt)==0) {
    message("WARNING: empty gene sets")
    return(NULL)
  }
  
  ## take out duplicated GO terms
  gmt.id <- gsub(".*\\(GO_|\\)$","",names(gmt))
  gmt.names <- names(gmt)
  names(gmt.names) <- gmt.id
  ndup <- sum(duplicated(gmt.id))
  message(paste("merging",ndup,"duplicated GO terms"))

  ## colllapse duplicates by set union
  gmt <- tapply(gmt, gmt.id, function(g) unique(unlist(g)))
  names(gmt) <- gmt.names[names(gmt)]
  
  ## sort on largest
  gmt <- gmt[order(-sapply(gmt,length))]
  
  ## filter on minimum size
  gmt <- gmt[sapply(gmt,length) >= minsize]

  return(gmt)
}
 
getOrganismGO.ANNOTHUB <- function(organism, features = NULL, use.ah = NULL, orgdb = NULL) {
  organism <- normalizeOrganism(organism) 

  ## Load the annotation resource.
  if (is.null(orgdb)) {
    orgdb <- getOrgDb(organism, use.ah = use.ah)
  }
  if (is.null(orgdb)) {
    message(paste0("[getOrganismGO.ANNOTHUB] WARNING : unsupported organism '", organism,"'"))
    return(NULL)
  }

  go.gmt <- list()
  ont_classes <- c("BP", "CC", "MF")
  if (!"GOALL" %in% AnnotationDbi::keytypes(orgdb)) {
    message("WARNING:: missing GO annotation in database!\n")
    return(NULL)
  }
  
  ## create GO annotets
  message(paste0("Creating GO annotation for '",organism,"'using AnnotationHub..."))
  ont_classes <- c("BP", "CC", "MF")
  k <- "BP"
  for (k in ont_classes) {
    cc <- AnnotationDbi::columns(orgdb)
    gene.column <- intersect(c("SYMBOL", "GENENAME", "MGI", "ALIAS"), cc)
    gene.column <- head(gene.column, 1)
    if (length(gene.column) > 0) {
      suppressMessages(suppressWarnings(
        go_id <- AnnotationDbi::mapIds(
          orgdb,
          keys = k, keytype = "ONTOLOGY",
          column = "GO", multiVals = "list"
        )[[1]]
      ))
      go_id <- unique(go_id)
      suppressMessages(suppressWarnings(
        sets <- AnnotationDbi::mapIds(
          orgdb,
          keys = go_id, keytype = "GOALL",
          column = gene.column, multiVals = "list"
        )
      ))
      
      ## get GO title
      sets <- sets[which(names(sets) %in% AnnotationDbi::keys(GO.db::GOTERM))]
      sets <- lapply(sets, function(s) unique(s))
      go <- sapply(GO.db::GOTERM[names(sets)], AnnotationDbi::Term)
      new_names <- paste0("GO_", k, ":", go, " (", sub("GO:", "GO_", names(sets)), ")")
      names(sets) <- new_names
      
      ## append to list
      go.gmt <- c(go.gmt, sets)
    }
  }
  go.gmt
}

getOrganismGO.GPROFILER <- function(organism, features, batch_size=2000,
                                    include_iea=TRUE, verbose=1) {

  id <- .map_gprofiler_id(organism) 
  id  
  if(is.null(id)) {
    message("[getOrganismGO.GPROFILER] WARNING: organism not found")
    return(NULL)
  }
  message(paste0("Getting GO gene sets using Gprofiler for id '", id,"'"))
  res <- NULL
  if(length(features) <= batch_size) {
    gost.out <- try(gprofiler2::gost(
      query = features,
      organism = id,
      evcodes = TRUE,   ## get gene sets
      exclude_iea = !include_iea,
      significant = FALSE,
      user_threshold = 1,
      correction_method = "fdr", 
      multi_query = FALSE
    ), silent=TRUE)
    if (!inherits(gost.out, "try-error")) {
      res <- gost.out$result
    }
  } else {

    n_genes <- length(features)
    n_batches <- ceiling(n_genes / batch_size)

    if (verbose > 0) cat("Processing ", n_genes, " features in ",
      n_batches, " batches")

    batch_results <- list()
    b=1
    for (b in 1:n_batches) {
      if (verbose > 0) cat('.')
      start_idx <- (b - 1) * batch_size + 1
      end_idx <- min(b * batch_size, n_genes)
      genes_batch <- features[start_idx:end_idx]
      
      batch_out <- try(gprofiler2::gost(
        query = genes_batch,
        organism = id,
        evcodes = TRUE,   ## get gene sets
        exclude_iea = !include_iea,
        significant = FALSE,
        user_threshold = 1,
        correction_method = "fdr", 
        multi_query = FALSE
      ), silent=TRUE)
      res <- batch_out$result
      
      ## Only keep successful results
      if (!"try-error" %in% class(batch_out) &&
            inherits(res, "data.frame") && nrow(res) > 0) {
        batch_results[[length(batch_results) + 1]] <- res
      }
    }
    if (verbose > 0) cat('\n')
    
    ## Combine all successful batch results
    if (length(batch_results) > 0) {
      res <- do.call(rbind, batch_results)
    } else {
      res <- try(stop("All batches failed"), silent = TRUE)
    }
  }

  ## create gmt list
  gmt <- strsplit(res$intersection, split=",")

  ## give names
  ss <- sub(":","_",res$source)
  nn <- res$term_name
  gg <- sub(":","_",res$term_id)
  gmt.name <- paste0(ss,":",nn," (",gg,")")
  head(gmt.name)
  names(gmt) <- gmt.name

  ## sum up batches
  sum(duplicated(gmt.name))
  if(  sum(duplicated(gmt.name)) ) {
    gmt <- tapply( 1:length(gmt), gmt.name, function(i) unlist(gmt[i]))
  }
  
  gmt <- gmt[order(-sapply(gmt,length))]
  return(gmt)
}


