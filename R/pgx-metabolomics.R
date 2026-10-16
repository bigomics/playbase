##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##


#' Check if RefMet server is alive
#'
mx.ping_refmet <- function() {
  out <- try(RefMet::refmet_metadata("Tyrosine"))
  if ("try-error" %in% class(out)) message("[mx.ping_refmet] WARNING! RefMet server is down")
  !("try-error" %in% class(out))
}

#' Convert IDs to CHEBI using base R functions (DEPRECATED)
#'
#' @export
convert_probe_to_chebi <- function(probes, probe_type = NULL) {
  message("WARNING: deprecated. please use mx.convert_probe()")
  annot <- .pb_libx("metabolite_id")
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

#' Convert IDs to CHEBI using the lib-x metabolite annotation table. Can
#' handle different target columns.
#'
#' @export
mx.convert_probe <- function(probes, probe_type = NULL, target_id = "ID") {
  # for id that are "", set it to na
  probes[probes == ""] <- NA
  annot <- .pb_libx("metabolite_id")
  ## probes <- sub("[ _.-].*","",probes)  ## strip any postfix??? Not good for lipid names
  if (is.null(probe_type)) {
    ## map the package probe type to a lib-x metabolite annotation column
    ptype <- playbase.annotation::detect_probetype(
      organism = NULL, probes = probes, datatype = "metabolomics"
    )
    cols <- colnames(annot)
    probe_type <- if (ptype %in% "NAME") "RefMet" else cols[match(toupper(ptype), toupper(cols))]
  }
  if (is.null(probe_type) || is.na(probe_type)) {
    message("WARNING: could not determine probe_type. Check your names.")
    return(rep(NA, length(probes)))
  }
  if (probe_type == "CHEBI") probe_type <- "ChEBI"
  if (!is.null(probe_type) && probe_type %in% c("ChEBI", "ID")) {
    # keep only numbers in ids, as chebi ids are numeric
    probes <- gsub("[^0-9]", "", probes)
  }
  # check that probetype is valid
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
  if (!probe_type %in% colnames(annot)) {
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
  na.probes <- is.na(probes) | probes %in% c("", "-")
  ids[na.probes] <- NA
  ids[ids == ""] <- NA
  return(ids)
}


#' Extend metabolite sets in M by adding neighbour-metabolites that
#' are connected (neighbours) according to given PPI. Only
#' metabolite-metabolite neighbours are considered.
#'
#' @export
extend_metabolite_sets <- function(M, ppi, add = TRUE, postfix = "(extended)",
                                   maxcost = 0.33) {
  ## get metabolite-metabolite edges from GRAPHITE
  ## ppi <- .pb_libx("graphite_ppi")
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
