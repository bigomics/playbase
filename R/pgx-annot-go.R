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
getOrganismGO <- function(organism, features=NULL, minsize=3L, batch_size=2000L,
                          db = c("annothub","gprofiler","uniprot"), include_iea=TRUE,
                          go.only = TRUE, annot = NULL) {

  gmt1=gmt2=gmt3=NULL
  if(is.null(db) || "annothub" %in% db) {
    gmt1 <- getOrganismGO.ANNOTHUB(organism, use.ah=NULL, orgdb=NULL)
    if(is.null(db) && !is.null(gmt1)) db <- "annothub"
    message(paste("Got",length(gmt1),"GO terms from AnnotHub"))
  }

  if(is.null(features) && !is.null(gmt1)) {
    features <- unique(unlist(gmt1))
  }

  ## fallback if annotHub didn't work. Skip if annothub has done it
  ## because gprofiler is quite slow.
  if(!is.null(features) && "gprofiler" %in% db && is.null(gmt1)) {    
    gmt2 <- getOrganismGO.GPROFILER(
      organism, features,
      batch_size = batch_size,
      include_iea = include_iea
    )
    if(go.only) {
      ii <- grep("\\(GO_[0-9]+\\)",names(gmt2))
      gmt2 <- gmt2[ii]
    }
    message(paste("Got",length(gmt2),"GO terms from Gprofiler"))
    table(sub(":.*","",names(gmt2)))
  }

  ## last resort: UniProt. Does not need features but downloads the
  ## whole proteome, so only if the others gave nothing.
  if("uniprot" %in% db && is.null(gmt1) && is.null(gmt2)) {
    gmt3 <- getOrganismGO.UNIPROT(organism, reviewed = FALSE)
    message(paste("Got",length(gmt3),"GO terms from UniProt"))
  }

  ## merge
  gmt <- c(gmt1, gmt2, gmt3)
  if(length(gmt)==0) {
    message("WARNING: empty gene sets")
    return(NULL)
  }
  
  ## check duplicated GO terms
  num.dup <- sum(duplicated(names(gmt)))
  num.dup
  if(num.dup>0) {
    gmt <- go.merge_duplicates(gmt)
  }

  ## convert all id to species symbol
  if(!is.null(annot) && "symbol" %in% colnames(annot)) {
    gmt <- gmt.map2symbol(gmt, annot=annot, target="symbol") 
    if(length(gmt)==0) {
      message("ERROR: empty gene sets after mapping to symbol")
      return(NULL)
    }
  }
  
  ## sort on largest
  gmt <- gmt[order(-sapply(gmt,length))]
  
  ## filter on minimum size
  gmt <- gmt[sapply(gmt,length) >= minsize]
  
  return(gmt)
}

getOrganismGO.ANNOTHUB <- function(organism, use.ah = NULL, orgdb = NULL) {
  organism <- normalizeOrganism(organism) 

  ## Load the annotation resource.
  if (is.null(orgdb)) {
    orgdb <- getOrgDb(organism, use.ah = use.ah)
  }
  if (is.null(orgdb)) {
    message(paste0("[getOrganismGO.ANNOTHUB] WARNING : unsupported organism '", organism,"'"))
    return(NULL)
  }

  gmt <- list()
  ont_classes <- c("BP", "CC", "MF")
  if (!"GOALL" %in% AnnotationDbi::keytypes(orgdb)) {
    message("WARNING:: missing GO annotation in database!\n")
    return(NULL)
  }
  
  ## create GO annotets
  message(paste0("[getOrganismGO.ANNOTHUB] Creating GO for '",organism,
    "' using AnnotationHub..."))
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
      gmt <- c(gmt, sets)
    }
  }

  dbg("[getOrganismGO.ANNOTHUB] done!")
  gmt
}

getOrganismGO.GPROFILER <- function(organism, features, batch_size=2000,
                                    include_iea=TRUE, verbose=1) {

  id <- .map_gprofiler_id(organism) 
  id  
  if(is.null(id)) {
    message("[getOrganismGO.GPROFILER] WARNING: organism not found")
    return(NULL)
  }
  message(paste0("Getting GO using Gprofiler for id '", id,"'"))
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

  if(inherits( res, "try-error")) {
    warning("[getOrganismGO.GPROFILER] Error.")
    return(NULL)
  }
  
  ## create gmt list
  gmt <- strsplit(res$intersection, split=",")

  ## standardize names (like "GO_BP:uridine kinase activity (GO_0004849)")
  ss <- gsub(":","_",res$source)
  nn <- res$term_name
  gg <- gsub(":","_",res$term_id)
  gmt.name <- paste0(ss,":",nn," (",gg,")")
  head(gmt.name)
  names(gmt) <- gmt.name

  ## sum up duplicated batches
  sum(duplicated(gmt.name))
  if(sum(duplicated(gmt.name))) {
    gmt <- tapply( 1:length(gmt), gmt.name,
      function(i) unname(unlist(gmt[i])))
  }
  
  gmt <- gmt[order(-sapply(gmt,length))]
  return(gmt)
}



#' Get GO gene sets for organism from UniProt (via UniProt.ws). Works
#' for any organism in UniProt, also those without OrgDb or g:Profiler
#' support. UniProt only gives direct annotations, so terms are
#' propagated to their GO.db ancestors (like GOALL in AnnotationHub).
getOrganismGO.UNIPROT <- function(organism, reviewed = FALSE, verbose = 1) {
  if (!requireNamespace("UniProt.ws", quietly = TRUE)) {
    message("[getOrganismGO.UNIPROT] UniProt.ws not installed; skipping")
    return(NULL)
  }
  taxid <- getTaxId.UNIPROT(organism)
  if (is.null(taxid)) {
    message("[getOrganismGO.UNIPROT] WARNING: unknown organism '", organism, "'")
    return(NULL)
  }
  if (verbose) message("[getOrganismGO.UNIPROT] Getting GO from UniProt for taxId ", taxid)

  query <- paste0("organism_id:", taxid)
  if (reviewed) query <- c(query, "reviewed:true")
  res <- try(UniProt.ws::queryUniProt(
    query = query,
    fields = c("accession", "gene_primary", "go_p", "go_c", "go_f"),
    collapse = " AND ",
    n = Inf,
    pageSize = 500L
  ), silent = TRUE)
  if (inherits(res, "try-error") || is.null(res) || nrow(res) == 0) {
    message("[getOrganismGO.UNIPROT] WARNING: no results from UniProt")
    return(NULL)
  }

  ## use gene symbol, fallback to accession
  gene <- res[[2]]
  gene <- ifelse(is.na(gene) | gene == "", res[[1]], gene)

  ancestors <- list(
    BP = as.list(GO.db::GOBPANCESTOR),
    CC = as.list(GO.db::GOCCANCESTOR),
    MF = as.list(GO.db::GOMFANCESTOR)
  )
  gmt <- list()
  for (i in 1:3) {
    k <- names(ancestors)[i]
    go <- regmatches(res[[i + 2]], gregexpr("GO:[0-9]{7}", res[[i + 2]]))
    gg <- rep(gene, sapply(go, length))
    go <- unlist(go)
    ## propagate to ancestors
    anc <- lapply(ancestors[[k]][go], setdiff, "all")
    gg <- c(gg, rep(gg, sapply(anc, length)))
    go <- c(go, unlist(anc))
    sets <- lapply(split(gg, go), unique)
    sets <- sets[names(sets) %in% AnnotationDbi::keys(GO.db::GOTERM)]
    term <- sapply(GO.db::GOTERM[names(sets)], AnnotationDbi::Term)
    names(sets) <- paste0("GO_", k, ":", term, " (", sub("GO:", "GO_", names(sets)), ")")
    gmt <- c(gmt, sets)
  }

  if (verbose) message("[getOrganismGO.UNIPROT] got ", length(gmt), " GO terms")
  gmt
}
