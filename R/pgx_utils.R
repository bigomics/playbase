
#' Convert probe names to symbol names
#'
#' More explanation is needed for this function
#'
#' @param probes unknown.
#' @param type unknown.
#' @param org unknown.
#' @param keep.na unknown.
#'
#' @return unknown.
#'
#' @examples
#' x <- 1
.probe2symbol <- function(probes, type = NULL, org = "human", keep.na = FALSE) {
  ## strip postfix for ensemble codes
  if (mean(grepl("^ENS", probes)) > 0.5) {
    probes <- gsub("[.].*", "", probes)
  }

  if (is.null(type)) {
    hs.list <- list(
      "human.ensembl" = unlist(as.list(org.Hs.eg.db::org.Hs.egENSEMBL)),
      "human.ensemblTRANS" = unlist(as.list(org.Hs.eg.db::org.Hs.egENSEMBLTRANS)),
      "human.refseq" = unlist(as.list(org.Hs.eg.db::org.Hs.egREFSEQ)),
      "human.accnum" = unlist(as.list(org.Hs.eg.db::org.Hs.egACCNUM)),
      "human.uniprot" = unlist(as.list(org.Hs.eg.db::org.Hs.egUNIPROT)),
      "human.symbol" = unlist(as.list(org.Hs.eg.db::org.Hs.egSYMBOL))
    )

    mm.list <- list(
      "mouse.ensembl" = unlist(as.list(org.Mm.eg.db::org.Mm.egENSEMBL)),
      "mouse.ensemblTRANS" = unlist(as.list(org.Mm.eg.db::org.Mm.egENSEMBLTRANS)),
      "mouse.refseq" = unlist(as.list(org.Mm.eg.db::org.Mm.egREFSEQ)),
      "mouse.accnum" = unlist(as.list(org.Mm.eg.db::org.Mm.egACCNUM)),
      "mouse.uniprot" = unlist(as.list(org.Mm.eg.db::org.Mm.egUNIPROT)),
      "mouse.symbol" = unlist(as.list(org.Mm.eg.db::org.Mm.egSYMBOL))
    )

    id.list <- c(hs.list, mm.list)
    mx <- sapply(id.list, function(id) mean(probes %in% id))
    max.mx <- max(mx, na.rm = TRUE)
    mx0 <- names(mx)[which.max(mx)]
    org <- sub("[.].*", "", mx0)
    type <- sub(".*[.]", "", mx0)
    if (max.mx < 0.5 && max.mx > 0) {
      warning("Low mapping ratio: r= ", max.mx)
    }
    if (max.mx == 0) {
      warning("Zero mapping ratio: r= ")
      type <- NULL
    }
  }
  if (is.null(type)) {
    warning("invalid type: ", type, "\n")
    return(NULL)
  }
  if (!type %in% c(
    "ensembl", "ensemblTRANS", "unigene", "refseq",
    "accnum", "uniprot", "symbol"
  )) {
    warning("invalid type: ", type, "\n")
    return(NULL)
  }

  if (type == "symbol") {
    if (any(grep(" /// ", probes))) {
      symbol0 <- strsplit(probes, split = " /// ")
    } else if (any(grep("[;,]", probes))) {
      symbol0 <- strsplit(probes, split = "[;,\\|]")
    } else {
      symbol0 <- probes
    }
  } else {
    if (org == "human") {
      symbol0 <- AnnotationDbi::mapIds(
        org.Hs.eg.db::org.Hs.eg.db,
        probes, "SYMBOL", toupper(type)
      )
    }
    if (org == "mouse") {
      symbol0 <- AnnotationDbi::mapIds(
        org.Mm.eg.db::org.Mm.eg.db,
        probes, "SYMBOL", toupper(type)
      )
    }
  }

  ## Unrecognize probes
  nna <- which(is.na(names(symbol0)))
  if (length(nna)) names(symbol0)[nna] <- probes[nna]

  ## What to do with unmapped/missing symbols????
  symbol <- sapply(symbol0, "[", 1) ## takes first symbol only!!!
  isnull <- which(sapply(symbol, is.null))
  symbol[isnull] <- NA
  if (keep.na) {
    sel.na <- which(is.na(symbol))
    symbol[sel.na] <- probes[sel.na]
  }
  symbol <- unlist(symbol)
  names(symbol) <- NULL

  return(symbol)
}
