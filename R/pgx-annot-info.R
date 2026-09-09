##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##


## ================================================================================
## ====================== GENE INFO FUNCTIONS======================================
## ================================================================================


#' @export
pgx.getFeatureInfo <- function(pgx, feature) {
  if (is.null(feature) || length(feature) == 0) {
    return(NULL)
  }
  if (is.na(feature) || feature == "") {
    return(NULL)
  }

  annot <- as.list(pgx$genes[feature, ])
  names(annot) <- sub("uniprot", "protein", names(annot))
  annot$gene_name <- NULL
  annot$chr <- NULL
  annot$pos <- NULL
  annot$tx_len <- NULL

  datatype <- pgx$datatype
  if (pgx$datatype == "multi-omics") {
    ## dtype <- c("mx"="metabolomics","px"="proteomics","gx"="transcriptomics")
    datatype <- ifelse(grepl("^mx:|metabolomics", feature), "metabolomics", datatype)
  }

  if (!"protein" %in% names(annot) && datatype != "metabolomics") {
    annot[["protein"]] <- gene2uniprot(annot$symbol, pgx$organism)
  }

  if (annot$human_ortholog %in% names(playdata::GENE_SUMMARY)) {
    annot.summary <- playdata::GENE_SUMMARY[annot$human_ortholog]
    annot.summary <- gsub("Publication Note.*|##.*", "", annot.summary)
    annot[["summary"]] <- annot.summary
  }

  if (datatype == "metabolomics") {
    annot <- getMetaboliteInfo(
      organism = pgx$organism,
      id = annot$symbol,
      info = annot
    )
  } else {
    annot <- info.add_hyperlinks(
      annot, feature, datatype,
      nm.symbol = "symbol",
      nm.prot = "protein",
      nm.ortholog = "human_ortholog",
      as.link = TRUE, add.summary = FALSE
    )
  }

  if (!is.null(annot)) {
    annot <- annot[!duplicated(names(annot))]
    ## reorder
    nn1 <- intersect(
      c(
        "gene_symbol", "organism", "name", "map_location",
        "uniprot", "databases", "summary", names(annot)
      ),
      names(annot)
    )
    nn2 <- setdiff(names(annot), nn1)
    annot <- annot[c(nn1, nn2)]
    names(annot) <- sub("gene_symbol", "symbol", names(annot))
    names(annot) <- sub("uniprot", "protein", names(annot))
    names(annot) <- sub("map_location", "genome location", names(annot))
  } else {
    annot <- list()
    annot$summary <- "(no info available)"
  }

  return(annot)
}

#' Retrieve gene information from different databases (AnnotHub,
#' orgdb) and wrap hyperlinks.
#'
#' @export
getOrgGeneInfo <- function(organism, gene, feature, ortholog, datatype,
                           as.link = TRUE) {
  if (is.null(gene) || length(gene) == 0) {
    return(NULL)
  }
  if (is.na(gene) || gene == "") {
    return(NULL)
  }

  orgdb <- getOrgDb(organism, use.ah = NULL)
  if (is.null(orgdb)) {
    message("[getOrgGeneInfo] WARNING could not get orgdb!")
    return(NULL)
  }
  cols <- c("SYMBOL", "UNIPROT", "GENENAME", "MAP", "OMIM", "PATH", "GO")
  cols <- intersect(cols, AnnotationDbi::keytypes(orgdb))

  if (!"SYMBOL" %in% cols) {
    keytype <- detect_probetype(organism, gene)
  } else {
    keytype <- "SYMBOL"
  }

  ## return if gene is not known
  if (!gene %in% AnnotationDbi::keys(orgdb, keytype)) {
    info <- list()
    info[["feature"]] <- feature
    info[["symbol"]] <- gene
    info[["organism"]] <- organism
    info[["summary"]] <- "(no info available)"
    return(info)
  }

  ## get info from different environments
  info <- lapply(cols, function(k) {
    tryCatch(
      {
        AnnotationDbi::select(
          orgdb,
          keys = gene,
          keytype = keytype,
          columns = k
        )[[k]]
      },
      error = function(w) {
        NULL
      }
    )
  })

  annot <- AnnotationDbi::select(
    orgdb,
    keys = gene,
    keytype = keytype,
    columns = cols
  )

  names(info) <- cols

  ## take out duplicates
  info <- lapply(info, unique)

  info[["ORGANISM"]] <- organism
  info <- info.add_hyperlinks(info, feature, datatype, as.link = TRUE)

  ## rename
  tags <- c(
    "ORGANISM", "SYMBOL", "UNIPROT", "GENENAME", "MAP", "OMIM", "PATH",
    "GO", "SUMMARY", "DATABASES"
  )
  tags <- intersect(tags, names(info))
  info <- info[tags]
  new.names <- c(
    "ORGANISM" = "organism", "SYMBOL" = "gene_symbol", "UNIPROT" = "uniprot",
    "GENENAME" = "name", "MAP" = "map_location",
    "OMIM" = "OMIM", "PATH" = "pathway", "GO" = "GO",
    "SUMMARY" = "summary", "DATABASES" = "databases"
  )
  names(info) <- new.names[tags]
  return(info)
}

info.add_hyperlinks <- function(info, feature, datatype,
                                nm.symbol = "SYMBOL", nm.prot = "UNIPROT",
                                nm.ortholog = "ORTHOLOG",
                                as.link = TRUE, add.summary = TRUE) {
  ## nm.symbol='SYMBOL';nm.prot='UNIPROT';nm.ortholog='ORTHOLOG';as.link=TRUE;add.summary=TRUE
  ## nm.symbol='symbol';nm.prot='protein';nm.ortholog='human_ortholog';as.link=TRUE;add.summary=TRUE

  symbol <- info[[nm.symbol]]
  ortholog <- info[[nm.ortholog]]
  uniprot <- info[[nm.prot]]

  if (length(uniprot) == 0) {
    this.uniprot <- NULL
  } else if (length(uniprot) == 1 && uniprot[1] == "" || is.na(uniprot[1])) {
    this.uniprot <- NULL
  } else {
    uniprot <- sort(unique(unlist(strsplit(uniprot, split = ";"))))
    this.uniprot <- uniprot[which(sapply(uniprot, function(p) grepl(p, feature)))]
    if (length(this.uniprot) == 0) this.uniprot <- rev(uniprot)[1]
  }

  if (as.link && length(symbol)) {
    gene.link <- "<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=GENE' target='_blank'>GENE</a>"
    gene.link <- sapply(symbol, function(s) gsub("GENE", s, gene.link))
    info[[nm.symbol]] <- paste(gene.link, collapse = ", ")
  }

  if (as.link && length(ortholog)) {
    gene.link <- "<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=GENE' target='_blank'>GENE</a>"
    ortho.link <- sapply(ortholog, function(s) gsub("GENE", s, gene.link))
    info[[nm.ortholog]] <- paste(ortho.link, collapse = ", ")
  }

  if (as.link && length(uniprot)) {
    prot.link <- "<a href='https://www.uniprot.org/uniprotkb/UNIPROT' target='_blank'>UNIPROT</a>"
    prot.link <- sapply(uniprot, function(s) gsub("UNIPROT", s, prot.link))
    info[[nm.prot]] <- paste(prot.link, collapse = ", ")
  }

  ## create link to external databases: OMIM, GeneCards, Uniprot
  if (as.link) {
    genecards.link <- "<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=GENE' target='_blank'>GeneCards</a>"
    uniprot.link <- "<a href='https://www.uniprot.org/uniprotkb/UNIPROT' target='_blank'>UniProtKB</a>"
    if (length(symbol)) genecards.link <- sub("GENE", symbol[1], genecards.link)
    if (length(this.uniprot)) uniprot.link <- sub("UNIPROT", this.uniprot, uniprot.link)
    db.links <- paste(c(genecards.link, uniprot.link), collapse = ", ")
    if (db.links == "") db.links <- NULL
    info[["DATABASES"]] <- db.links
  }

  if (length(this.uniprot) && grepl("proteomics", datatype, ignore.case = TRUE)) {
    ## create link to PhosphoSitePlus
    phosphositeplus.link <- "<a href='https://www.phosphosite.org/simpleSearchSubmitAction.action?searchStr=GENE' target='_blank'>PhosphoSitePlus</a>"
    phosphositeplus.link <- "<a href='https://www.phosphosite.org/uniprotAccAction?id=UNIPROT' target='_blank'>PhosphoSitePlus</a>"
    ## phosphositeplus.link <- sub("GENE", symbol[1], phosphositeplus.link)
    phosphositeplus.link <- sub("UNIPROT", this.uniprot, phosphositeplus.link)
    info[["DATABASES"]] <- paste(c(info[["DATABASES"]], phosphositeplus.link), collapse = ", ")

    ## ## create links to PhosphoELM for proten and gene: db of S/T/Y phosphorylation sites
    ## phosphoELM.link <- "<a href='http://phospho.elm.eu.org/byAccession/UNIPROT' target='_blank'>PhosphoELM</a>"
    ## phosphoELM.link <- sub("UNIPROT", uniprot, phosphoELM.link)
    ## info[["DATABASES"]] <- paste(c(info[["DATABASES"]], phosphoELM.link), collapse = ", ")
  }

  ## create link to OMIM
  if (as.link && length(info[["OMIM"]])) {
    omim.link <- "<a href='https://www.omim.org/entry/OMIM' target='_blank'>OMIM</a>"
    info[["OMIM"]] <- sapply(info[["OMIM"]], function(x) gsub("OMIM", x, omim.link))
  }

  ## create link to KEGG
  if (as.link && length(info[["PATH"]])) {
    kegg.link <- "<a href='https://www.genome.jp/kegg-bin/show_pathway?map=hsaKEGGID&show_description=show' target='_blank'>KEGGNAME (KEGGID)</a>"
    for (i in 1:length(info[["PATH"]])) {
      kegg.id <- info[["PATH"]][[i]]
      kegg.id <- setdiff(kegg.id, NA)
      if (length(kegg.id) > 0) {
        kegg.name <- AnnotationDbi::mget(kegg.id, envir = KEGG.db::KEGGPATHID2NAME, ifnotfound = NA)[[1]]
        if (!is.na(kegg.name) && as.link) {
          info[["PATH"]][[i]] <- gsub("KEGGNAME", kegg.name, gsub("KEGGID", kegg.id, kegg.link))
        } else {
          info[["PATH"]][[i]] <- kegg.name
        }
      }
    }
  }

  ## create link to GO
  if (length(info[["GO"]]) && !is.na(info[["GO"]][1])) {
    ## sometimes GO.db is broken...
    suppressWarnings(try.out <- try(AnnotationDbi::Term(AnnotationDbi::mget("GO:0000001",
      envir = GO.db::GOTERM, ifnotfound = NA
    )[[1]])))
    go.ok <- (!"try-error" %in% class(try.out))
    if (go.ok) {
      amigo.link <- "<a href='http://amigo.geneontology.org/amigo/term/GOID' target='_blank'>GOTERM (GOID)</a>"
      i <- 1
      for (i in 1:length(info[["GO"]])) {
        go_id <- info[["GO"]][i]
        term_id <- AnnotationDbi::mget(go_id, envir = GO.db::GOTERM, ifnotfound = NA)[[1]]
        if (class(term_id) == "GOTerms") {
          go_term <- AnnotationDbi::Term(term_id)
          if (as.link) {
            info[["GO"]][i] <- gsub("GOTERM", go_term, gsub("GOID", go_id, amigo.link))
          } else {
            info[["GO"]][i] <- go_term
          }
        } else {
          info[["GO"]][i] <- go_id
        }
      }
    } else {
      info[["GO"]] <- NULL
    }
  }

  if (add.summary) {
    ## pull summary
    info[["SUMMARY"]] <- "(no info available)"
    ortholog <- info[[nm.ortholog]]
    if (!is.null(ortholog) && ortholog %in% names(playdata::GENE_SUMMARY)) {
      info[["SUMMARY"]] <- playdata::GENE_SUMMARY[ortholog]
      info[["SUMMARY"]] <- gsub("Publication Note.*|##.*", "", info[["SUMMARY"]])
    }
  }

  return(info)
}

