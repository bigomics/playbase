##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##


##---------------------------------------------------------------------
##----------------------- CONTENT CREATION ----------------------------
##---------------------------------------------------------------------


table_to_content <- function(df) {
  if(is.null(df)) return(NULL)
  rownames(df) <- iconv2ascii(rownames(df))
  for(i in seq_len(ncol(df))) {
    if(is.character(df[[i]])) df[[i]] <- iconv2ascii(df[[i]])
  }
  paste(as.character(knitr::kable(df,format="markdown")),collapse="\n")
}

list_to_content <- function(a, newline=FALSE) {
  if(is.null(a)) return("")
  aa <- sapply(a, function(s) paste(unlist(s), collapse='; '))
  pp <- paste0("- **",names(a),"**: ")
  if(newline) pp <- paste0(pp, '\n\n')
  cc <- paste(paste0(pp,aa), collapse='\n')
  paste(cc,'\n')
}

collate_as_sections <- function(a, level=2) {
  if(is.null(a)) return("")
  hdr <- paste(rep("#",level),collapse="")
  for(i in 1:length(a)) {
    a[[i]] <- paste(hdr,names(a)[i],"\n\n",a[[i]],'\n')
  }
  paste(a, collapse="\n\n")
}

#'
#'
#' @export
ai.create_report <- function(pgx, ntop=20, sections=NULL, collate=FALSE) {
  
  contrasts <- pgx.getContrasts(pgx)
  samples <- rownames(pgx$samples)
  ct <- contrasts[1]  ## FOR NOW!!!!

  all.sections <- c("description", "dataset_info","compute_settings",
    "differential_expression", "geneset_enrichment",
    "drug_similarity","hub_genes","wgcna_report")
  if(is.null(sections)) {
    sections <- all.sections
  } else {
    sections <- intersect(all.sections, sections)
  }

  description <- NULL
  if("description" %in% sections) {
    description <- list(
      name = pgx$name,
      title = pgx$description,
      description = pgx$description
    )
  }

  dataset_info <- NULL
  if("dataset_info" %in% sections) {
    dataset_info <- list(
      name = pgx$name,
      organism = pgx$organism,
      datatype = pgx$datatype,
      creation_date = pgx$date,
      num_samples = ncol(pgx$X),
      num_features = nrow(pgx$X),
      num_genesets = nrow(pgx$gsetX),
      num_comparisons = length(contrasts),
      samples = samples,
      comparisons = contrasts
      # features = c(head(rownames(pgx$counts),40),ifelse(nrow(pgx$counts)>40,"...",""))
    )
  }

  compute_settings <- NULL
  if("compute_settings" %in% sections) {
    compute_settings <- list(
      gene_tests = colnames(pgx$gx.meta$meta[[1]]$p),
      geneset_tests = colnames(pgx$gset.meta$meta[[1]]$p),
      pgx_slots = sort(names(pgx))
    )
  }
  
  ##-------------------------------------------------------------------
  differential_expression <- NULL
  if("differential_expression" %in% sections) {
    F <- pgx.getMetaMatrix(pgx)$fc
    F <- rename_by2(F, pgx$genes, "symbol")
    ii <- match(rownames(F),pgx$genes$symbol)
    rownames(F) <- paste0(pgx$genes$gene_title[ii]," (",rownames(F),")")
    F.up <- head( F[order(-rowMeans(F)),,drop=FALSE], 2*ntop )
    F.dn <- head( F[order(rowMeans(F)),,drop=FALSE], 2*ntop )
    
    differential_expression <- list(
      "Up-regulated genes (top hits). The top most most positively differentially expressed genes are:" = F.up,
      "Down-regulated genes (top hits). The top most negatively differentially expressed genes are:" = F.dn
    )
    differential_expression <- lapply(differential_expression, table_to_content)
  }
  
  ##-------------------------------------------------------------------
  geneset_enrichment <- NULL
  if("geneset_enrichment" %in% sections) {
    G <- pgx.getMetaMatrix(pgx, level = "geneset")$fc
    G <- G[order(-rowMeans(G)),,drop=FALSE]  ## NEED RETHINK: rowmeans is not good...
    revtail <- function(A,n) {
      if(nrow(A)==0) return(A[0,,drop=FALSE])
      head(A[nrow(A):1,,drop=FALSE],n)
    }
    BP.up <- head( G[grep("GOBP|GO_BP",rownames(G)),,drop=FALSE], ntop)
    BP.dn <- revtail( G[grep("GOBP|GO_BP",rownames(G)),,drop=FALSE], ntop)
    MF.up <- head( G[grep("GOMF|GO_MF",rownames(G)),,drop=FALSE], ntop)
    MF.dn <- revtail( G[grep("GOMF|GO_MF",rownames(G)),,drop=FALSE], ntop)
    CC.up <- head( G[grep("GOCC|GO_CC",rownames(G)),,drop=FALSE], ntop)
    CC.dn <- revtail( G[grep("GOCC|GO_CC",rownames(G)),,drop=FALSE], ntop)
    PW.up <- head( G[grep("PATHWAY",rownames(G)),,drop=FALSE], ntop)
    PW.dn <- revtail( G[grep("PATHWAY",rownames(G)),,drop=FALSE], ntop)
    
    geneset_enrichment <- list(
      "Top most positively enriched GO biological process (BP) gene sets:" = BP.up,
      "Top most negatively enriched GO biological process (BP) gene sets:" = BP.dn,
      "Top most positively enriched GO molecular function (MF) gene sets:" = MF.up,
      "Top most negatively enriched GO molecular function (MF) gene sets:" = MF.dn,
      "Top most positively enriched GO cellular component (CC) gene sets:" = CC.up,
      "Top most negatively enriched GO cellulare component (CC) gene sets:" = CC.dn,
      "Top most positively enriched pathways:" = PW.up,
      "Top most negatively enriched pathways:" = PW.dn
    )
    geneset_enrichment <- lapply(geneset_enrichment, table_to_content)
  }
  
  ##-------------------------------------------------------------------
  drug_similarity <- NULL
  if("drug_similarity" %in% sections && !is.null(pgx$drugs)) {
    ## NEED RETHINK: Only reports first contrast at the moment
    D <- pgx.getTopDrugs(pgx, ct, n=ntop, na.rm=TRUE)    
    drug_similarity <- list(
      "Drug Mechanism of Action. Drug Connectivity Map (CMap) analysis of selected comparison. Similarity of the mechanism of action (MOA) is based on correlation enrichment with drug perturbation profiles of LINCS L1000 database. The top most similar (i.e. positively correlated) drugs are:" =
        table_to_content(pgx.getTopDrugs(pgx, ct, n=ntop, dir=+1, na.rm=TRUE)), 
      "The top most inhibitory (i.e. negative correlated) drugs are:" = 
        table_to_content(pgx.getTopDrugs(pgx, ct, n=ntop, dir=-1, na.rm=TRUE))
    )
  }
  
  ##-------------------------------------------------------------------
  pcsf_report <- NULL
  if(FALSE && "pcsf_report" %in% sections) {
    pcsf_report <- list("Identification of hub genes. Hub genes can identify important regulators. The hub score is computed using a page rank network centrality score. The most central genes are:" = table_to_content(pgx.getPCSFcentrality(pgx, ct, plot = FALSE, n = ntop))
    )
  }
  
  ##-------------------------------------------------------------------
  wgcna_report <- NULL
  if("wgcna_report" %in% sections && !is.null(pgx$wgcna)) {

    out <- wgcna.describeModules(
      pgx$wgcna,
      modules = NULL,
      multi = FALSE, 
      ntop = 40,
      annot = pgx$genes, 
      experiment = pgx$description,
      verbose = FALSE,
      model = NULL
    ) 

    names(out)
    wgcna_report <- list_to_content(out$answers)
  }
    
  ##-----------------------------------------------------------
  report <- list(
    DATASET_INFO = list_to_content(dataset_info)
  )

  content <- list(
    description = list_to_content(description),
    dataset_info = list_to_content(dataset_info),
    compute_settings = list_to_content(compute_settings),
    differential_expression = list_to_content(differential_expression),
    geneset_enrichment = list_to_content(geneset_enrichment),
    drug_similarity = list_to_content(drug_similarity),
    pcsf_report = list_to_content(pcsf_report),
    wgcna_report = wgcna_report
  )
  
  if(collate) {
    content <- collate_as_sections(content,level=2) 
  }

  return(content)
}

#' Export markdown to PDF 
#' 
#' @export
markdownToPDF <- function(text, file) {

  text <- gsub(intToUtf8("8209"),"-",text)
  ##text <- gsub("---\n","",text)
  
  ## header/frontmatter
  hdr <- paste0("---\n")
  hdr <- paste0(hdr, "format:\n")
  hdr <- paste0(hdr, "  pdf:\n")
  #hdr <- paste0(hdr, "    documentclass: report\n")  
  hdr <- paste0(hdr, "    pdf-engine: lualatex\n")
  hdr <- paste0(hdr, "    documentclass: article\n")
  hdr <- paste0(hdr, "    papersize: a4\n")
  hdr <- paste0(hdr, "    geometry:\n")
  hdr <- paste0(hdr, "      - left=25mm\n")
  hdr <- paste0(hdr, "      - right=20mm\n")
  hdr <- paste0(hdr, "      - top=25mm\n")
  hdr <- paste0(hdr, "      - bottom=25mm\n")
  hdr <- paste0(hdr, "    mainfont: Lato\n")
  hdr <- paste0(hdr,"---\n") 

  text <- sub("^#", paste0(hdr,"\n#"), text)
  
  curwd <- getwd()
  tmpdir <- tempdir()
  setwd(tmpdir)
  md.file <- file.path(tmpdir,"report.md")
  write(text, file=md.file)
  quarto::quarto_render(md.file, output_format="pdf", quiet=TRUE)
  pdf.file <- file.path(tmpdir,"report.pdf")
  setwd(curwd)
  file.copy(pdf.file, file, overwrite=TRUE)
  return(file)
}


