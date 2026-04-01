##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##


##---------------------------------------------------------------------
##----------------------- CONTENT CREATION ----------------------------
##---------------------------------------------------------------------


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


#' Format and render contents as section.
#'
#' @export
rpt.compile_sections <- function(contents, hlevel=2, shift=TRUE) {
  
  div.description <- contents$description
  div.bullets <- contents$bullets

  ## render report with slight grey background
  div.report <- contents$report
  div.report <- gsub("\n---","\n",div.report)  ## remove hline
  div.report <- paste("\n::: {style='background-color: #eee;'}",
    div.report, ":::\n",sep="\n\n")
  
  div.methods <- contents$methods
  div.settings <- list_to_list(contents$settings)

  div.references <- contents$references
  div.references <- list_to_list(
    as.list(div.references),
    type = 'ol',
    add.name = FALSE
  )

  div.figures <- as.list(contents$figures)
  div.figures <- md.list_to_figs(div.figures, labels=NULL)  

  div.tables <- list()
  for(i in 1:length(contents$tables)) {
    tbl <- table_to_content(
      contents$tables[[i]]
      #caption = names(contents$tables)[i]
    )
    tbl <- paste("\n::: {style='font-size: 8pt;'}",
      tbl, ":::\n",sep="\n\n")    
    div.tables[[i]] <- tbl
  }
  names(div.tables) <- names(contents$tables)
  div.tables <- list_to_content(div.tables, newline=TRUE) 
  
  divs <- list(
    "Description" = div.description,
    "Highlights" = div.bullets,
    "AI Report" = div.report,
    "Methods" = div.methods,  
    "Settings" = div.settings,
    "References" = div.references,    
    "Figures" = div.figures,
    "Tables" = div.tables
  )
  names(divs) <- toupper(names(divs))
  names(divs) <- paste0("[",names(divs),"]{.underline}")  ## sometimes problemati

  for(i in 1:length(divs)) {
    divname <- names(divs)[i]
    divs[[i]] <- paste(divs[[i]], "\n\n---\n\n")
    if(hlevel>1 && grepl("references|figures", divname, ignore.case=TRUE) ) {
      #divs[[i]] <- paste0(divs[[i]], "\\newpage\n")
      divs[[i]] <- paste0(divs[[i]], "\n\n{{< pagebreak >}} &nbsp; \n\n")
    }
  }

  collate_as_sections(divs, csep="", hlevel=hlevel, shift=shift)  
}




#'
#'
#' @export
rpt.compile_drugconnectivity_report <- function(obj, which.db = 1, report = NULL, 
                                                hlevel = 2, shift = TRUE,
                                                model = NULL, pgx=NULL,
                                                title = "Drug Connectivity Map Analysis") {

  rpt <- NULL
  if(!is.null(obj)) {  
    if(all(c("drugs","genes","counts") %in% names(obj))) {
      ## is pgx object
      drugs <- obj$drugs
    } else {
      ## is drugs slot object
      drugs <- obj
    }
    ## take selected database
    if(is.numeric(which.db)) {
      which.db <- names(drugs)[which.db]
    }
    rpt <- drugs[[which.db]]$report
  } else if(!is.null(report))  {
    ## check if we have valid report results    
    rpt <- report
  } else {
    stop("invalid parameters: missing obj or report")
  }
  
  has.report <- !is.null(rpt)
  has.report
  if(!has.report) {
    if(!is.null(model)) {
      message("Warning: missing report results. Computing...")
      rpt <- ai.create_report_drug_connectivity(
        pgx, model, model2 = NULL, db = which.db,
        user.prompt = NULL)
    } else {
      stop("Error: missing report results. recompute pgx$drugs or specify LLM model.")
    }
  }

  has.infographic <- has.report && "infographic" %in% names(rpt)    
  has.infographic
  if(!has.infographic) {
    message("Warning: missing infographic")
    ##return(NULL)
  }
  
  ##------- description -------------
  div.description <- 
    "The Drug Connectivity Map (Drug CMap) is a resource and tool in the field of pharmacogenomics and drug discovery. It was developed by the Broad Institute (Lamb et al, 2006), and it provides a large-scale collection of gene expression profiles in response to more than 5000 compounds, for a total of more than a million gene expression profiles. The primary goal of the Drug CMap is to help researchers identify connections between drugs, genes, and diseases, facilitating the discovery of new therapeutic targets and repurposing existing drugs for different indications."
  
  div.methods <-
    "Log2 fold-changes (log2FC) were extracted for all computed contrasts. Common genes between the log2FC matrix and L1000 database are identified. If less than 20 genes are shared, the analysis is not conducted due to lack of statistical power. Drug meta sets are defined from available drugs by assigning each drug into principal drug classes. Only drug classes containing at least 10 distinct drug reports are retained for analyses. For expression data and L1000 drug matrix the rank of genes is calculated for each sample and drug, respectively, using the colRanks function from the matrixStats R package, with 'average' as method to treat ties. Rank pairwise Pearson's correlation (R) is then computed between the two rank matrices. Small, random gaussian noise is added to the correlation matrix. A sparse, binary (absence/presence) drug model matrix is created and rank correlation between the rank correlation matrix R and the drug model matrix is computed to assess potential association with drug profile and drug-response. A p-value and FDR are also computed for each drug. Drug set enrichment is also performed for each contrast using fgsea on the pre-defined drug meta sets and the rank correlation matrix. GSEA normalized enrichment score (NES) along with p and q-values for each drug are computed. For both rank correlation and GSEA analysis, the drugs are ranked (by correlation or NES) to identify the top 1000 matching drugs for each contrast. Further analyses and visualizations are conducted using the rank correlation values (R) and GSEA NES scores."
  
  div.refs <- c("Subramanian, A., Tamayo, P., Mootha, V. K., .... & Mesirov, J. P. (2005). Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles. PNAS, 102(43), 15545-15550. https://www.pnas.org/doi/10.1073/pnas.0506580102",
    "Korotkevich, G., Sukhov, V., Budin, N., .... & Sergushichev, A. (2021). Fast gene set enrichment analysis BioRxiv. https://www.biorxiv.org/content/10.1101/060012v3")
  
  ##------- figures -------------
  
  rx=2  ## resolution parameter
  fig1 <- tempfile(fileext='.png')
  png(fig1, width=700*rx, height=300*rx, pointsize=14*rx)
  pgx.plotDrugActivationMap(
    pgx,
    dmethod = which.db,
    contrast = NULL,
    nterms = 40,
    nfc = 20,
    normalize = FALSE,
    drugs = NULL,
    colorbar = FALSE,
    rotate = FALSE,
    plotlib = "base"
  )
  dev.off()

  fig2 = NULL
  fig2 <- tempfile(fileext='.png')
  if(!is.null(rpt$infographic)) {
    img <- rpt$infographic
    if(class(img) == 'array' && length(dim(img))==3 ) {
      png::writePNG(img, target=fig2)
    } else if(is.character(rpt$infographic[1])) {
      fig2 <- rpt$infographic
      if(!file.exists(fig2)) fig2 <- NULL
    }
  } else {
    png(fig2,w=800,h=450)
    par(mar=c(0,0,0,0))
    plot.new()
    text(0.5,0.5,"missing infographic")
    dev.off()
  }

  figs <- list(
    "The Activation Matrix visualizes the activation of drug activation enrichment across the conditions. The size of the circles correspond to their relative activation, and are colored according to their upregulation (red) or downregulation (blue) in the contrast profile." = fig1,
    "Infographic (AI generated)" = fig2
  )
  figs <- figs[!sapply(figs,is.null)]
  figs <- lapply(figs, function(f) paste0('/',f))
  
  ##------- tables -------------

  df1 <- pgx.getMOAmatrix(pgx, db=which.db, type="drugClass")
  df2 <- pgx.getMOAmatrix(pgx, db=which.db, type="targetGene")
  df1 <- round(df1, digits=3)
  df2 <- round(df2, digits=3)  
  
  order1 <- order(-rowMeans(df1))
  df1.up <- head(df1[order1,,drop=FALSE],10)
  df1.dn <- head(df1[rev(order1),,drop=FALSE],10)  

  order2 <- order(-rowMeans(df2))
  df2.up <- head(df2[order2,,drop=FALSE],15)
  df2.dn <- head(df2[rev(order2),,drop=FALSE],15) 
  
  div.tables <- list(
    "Top similar MOA drug class" = df1.up,
    "Top opposite MOA drug class" = df1.dn,    
    "Top similar MOA target gene" = df2.up,
    "Top opposite MOA target gene" = df2.dn    
  )
  
  ##------- create sections -------------
  contents <- list(
    description = div.description,
    bullets = rpt$bullets,
    report = rpt$report,
    methods = div.methods,  
    ##settings = wgcna$settings,
    references = div.refs,
    figures = figs,
    tables = div.tables
  )
  txt <- rpt.compile_sections(contents, hlevel=hlevel, shift=shift)
  
  ## add title
  title <- paste(title,": ",which.db)
  txt <- paste0("# ",title,"\n\n", txt)

  return(txt)
}
