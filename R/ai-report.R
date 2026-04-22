##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##


##---------------------------------------------------------------------
##------------------------ REPORT CREATION ----------------------------
##---------------------------------------------------------------------


#' @export
ai.create_report <- function(pgx, ...) {
  message("DEPRECATED: ai.create_report()")
  rpt.create_full_report(pgx, ...) 
}

#' This is a full report of all analyses of the experiment and can
#' serve as data source for LLM questions.
#'
#' @export
rpt.create_full_report <- function(pgx, ntop=20, llm=NULL,
                                   sections=NULL, collate=FALSE) {

  if(0) {
    ntop=20;sections=NULL;collate=FALSE
  }
  
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
    differential_expression <- lapply(differential_expression, round, digits=3)
    differential_expression <- lapply(differential_expression, table_to_content)    
  }
  
  ##-------------------------------------------------------------------
  geneset_enrichment <- NULL
  if("geneset_enrichment" %in% sections) {
    G <- pgx.getMetaMatrix(pgx, level = "geneset")$fc
    G <- G[order(-rowMeans(G**2)),,drop=FALSE] 
    BP <- head( G[grep("GOBP|GO_BP",rownames(G)),,drop=FALSE], ntop)
    MF <- head( G[grep("GOMF|GO_MF",rownames(G)),,drop=FALSE], ntop)
    CC <- head( G[grep("GOCC|GO_CC",rownames(G)),,drop=FALSE], ntop)
    PW <- head( G[grep("^PATHWAY",rownames(G)),,drop=FALSE], ntop)
    #DR <- head( G[grep("^DRUG",rownames(G)),,drop=FALSE], ntop)
    #DS <- head( G[grep("^DISEASE",rownames(G)),,drop=FALSE], ntop)        
    
    geneset_enrichment <- list(
      "Top enriched GO biological process (BP) gene sets:" = BP,
      "Top enriched GO molecular function (MF) gene sets:" = MF,
      "Top enriched GO cellular component (CC) gene sets:" = CC,
      "Top enriched pathways:" = PW
    )
    geneset_enrichment <- lapply(geneset_enrichment, round, digits=3)
    geneset_enrichment <- lapply(geneset_enrichment, table_to_content)
  }

  ##-------------------------------------------------------------------
  drug_similarity <- NULL
  if("drug_similarity" %in% sections && !is.null(pgx$drugs)) {

    ## NEED RETHINK: Only reports first contrast at the moment
    D <- pgx.getTopDrugs(pgx, ct, n=ntop, na.rm=TRUE)    
    drug_similarity <- list(
      "Drug Mechanism of Action. Drug Connectivity Map (CMap) analysis of selected comparison. Similarity of the mechanism of action (MOA) is based on correlation enrichment with drug perturbation profiles of LINCS L1000 database. The top most similar (i.e. positively correlated) drugs are:\n\n" =
        pgx.getTopDrugs(pgx, ct, n=ntop, dir=+1, na.rm=TRUE), 
      "The top most inhibitory (i.e. negative correlated) drugs are:\n\n" = 
        pgx.getTopDrugs(pgx, ct, n=ntop, dir=-1, na.rm=TRUE)
    )
    #drug_similarity <- lapply(drug_similarity, round, digits=3)    
    drug_similarity <- lapply(drug_similarity, table_to_content)    


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
      ntop = 50,
      annot = pgx$genes, 
      experiment = pgx$description,
      verbose = FALSE,
      model = NULL
    ) 

    names(out)
    wgcna_report <- list_to_content(out$answers, newline=TRUE)
  }
    
  ##-----------------------------------------------------------
  report <- list(
    DATASET_INFO = list_to_content(dataset_info)
  )

  content <- list(
    description = list_to_content(description),
    dataset_info = list_to_content(dataset_info),
    compute_settings = list_to_content(compute_settings),
    differential_expression = list_to_content(differential_expression, newline=TRUE),
    geneset_enrichment = list_to_content(geneset_enrichment, newline=TRUE),
    drug_similarity = list_to_content(drug_similarity, newline=TRUE),
    pcsf_report = list_to_content(pcsf_report),
    wgcna_report = wgcna_report
  )

  ## clean-up
  #content <- lapply(content, function(s) gsub("[-=_]+{4}","",s))
  content <- lapply(content, function(s) gsub("[=]+{3}","",s))
  
  if(collate) {
    content <- collate_as_sections(
      content, csep='\\newpage')
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
rpt.compile_wgcna_report <- function(obj, report = NULL, 
                                     hlevel=2, shift=TRUE,
                                     title="WGCNA Analysis Report") {

  if(all(c("wgcna","genes","counts") %in% names(obj))) {
    ## is pgx object
    wgcna <- obj$wgcna
  } else {
    ## is wgcna slot object
    wgcna <- obj
  }
  
  rpt <- wgcna$report
  if(!is.null(report))  rpt <- report
  
  has.report <- "report" %in% names(wgcna)
  has.infographic <- has.report && "infographic" %in% names(wgcna$report)    
  has.report
  has.infographic
  
  if(!has.report) {
    message("Error: missing report results. Please run wgcna.init()")
    return(NULL)
  }
  if(!has.infographic) {
    message("Warning: missing infographic")
    ##return(NULL)
  }
  
  ##------- description -------------
  div.description <- 
    "**Weighted Gene Co-expression Network Analysis (WGCNA) is a systems biology method for finding clusters (modules) of highly correlated genes in high-dimensional data, such as RNA-seq or microarray samples. It identifies modules based on similar expression patterns, relates them to external sample traits, and finds key hub genes, enabling the transition from simple gene lists to functional biological insights.**"
  
  div.methods <-
    "WGCNA computation was done in Omics Playground (BigOmics Analytics, Switzerland) using the `WGCNAplus` R package. The optimal soft threshold β was determined by optimizing the IQR of the dendrogram heights. The adjacency matrix was transformed into a topological overlapping matrix (TOM) and used for hierarchical clustering of the features. Subsequently, the dynamic tree‐cutting algorithm was employed to identify gene modules, each comprising a minimum of genes. Similar modules were merged using a shear height threshold. Finally, modules were distinguished by different colours and represented by their module eigengene (ME). Module-trait relationship was determined based on the Pearson correlation between ME and clinical features. Gene significance (GS) values were computed as the product of the Module Membership (MM) (correlation between expression and module eigengene) and maximum Trait Significance (TS) (correlation between expression and clinical characteristics) values."
  
  div.refs <- c(
    "Langfelder, P., Horvath, S. WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics 9, 559 (2008). https://doi.org/10.1186/1471-2105-9-559",
    "Langfelder P, Luo R, Oldham MC, Horvath S (2011) Is My Network Module Preserved and Reproducible? PLoS Comput Biol 7(1): e1001057. https://doi.org/10.1371/journal.pcbi.1001057",
    "WGCNAplus R package: `https://github.com/bigomics/WGCNAplus`"
  )
  
  is.mox <- all(c("layers","lasagna") %in% names(wgcna))
  is.mox
  
  modTraits <- wgcna$modTraits
  if(is.mox) {
    mm <- lapply(wgcna$layers, function(w) w$modTraits)
    modTraits <- do.call(rbind, mm)
  }
  dim(modTraits)
  
  ##------- figures -------------
  rx=2  ## resolution parameter
  
  fig1 <- tempfile(fileext='.png')
  if(is.mox) {
    nr <- ceiling(length(wgcna$layers)/2)
    png(fig1, w=800*rx, h=nr*300*rx, pointsize=14*rx)    
    wgcna.plotMultiDendroAndColors(wgcna$layers, show.traits=1,
      main="Cluster dendrogram and module colors")
    dev.off()
  } else {
    png(fig1, w=800*rx, h=400*rx, pointsize=14*rx)
    wgcna.plotDendroAndColors(wgcna, show.traits=1,
      main="Cluster dendrogram and module colors")
    dev.off()
  }
  
  fig2 <- tempfile(fileext='.png')
  tr <- (ncol(modTraits) < nrow(modTraits))
  png(fig2, w=800*rx, h=500*rx, pointsize=14*rx)
  if(is.mox) {
    wgcna.plotModuleTraitHeatmap(
      wgcna$layers, multi=TRUE, cluster=TRUE, text=FALSE,
      transpose=tr, main="Module-trait correlation")
  } else {
    wgcna.plotModuleTraitHeatmap(
      wgcna, multi=FALSE,
      cluster=TRUE, text=FALSE,
      transpose=tr, main="Module-trait correlation")
  }
  dev.off()

  fig3 = NULL
  fig3 <- tempfile(fileext='.png')
  if(!is.null(rpt$infographic)) {
    img <- rpt$infographic
    if(class(img) == 'array' && length(dim(img))==3 ) {
      png::writePNG(img, target=fig3)
    } else if(is.character(rpt$infographic[1])) {
      fig3 <- rpt$infographic
      if(!file.exists(fig3)) fig3 <- NULL
    }
  } else {
    png(fig3,w=800,h=450)
    plot.new()
    text(0.5,0.5,"missing infographic")
    dev.off()
  }

  fig4 <- NULL
  if(!is.null(rpt$diagram)) {
    fig4svg <- tempfile(fileext='.svg')
    fig4 <- tempfile(fileext='.png')
    dd <- DiagrammeR::grViz(rpt$diagram)
    svg <- try(DiagrammeRsvg::export_svg(dd))
    if(!"try-error" %in% class(svg)) {
      write(svg, file=fig4svg)
      rsvg::rsvg_png(fig4svg, file=fig4, width=2400)
    } else {
      fig4 <- NULL
    }
  }

  if(is.null(fig4)) {
    fig4 <- tempfile(fileext='.png')
    png(fig4,w=800,h=450)
    plot.new()
    text(0.5,0.5,"missing diagram")
    dev.off()
  }
  
  fig5 = NULL
  if(FALSE && is.mox && !is.null(wgcna$graph)) {
    fig5 <- tempfile(fileext='.png')
    png(fig5, w=1200, h=700, pointsize=14)
    plotMultiPartiteGraph2(
      wgcna$graph, ntop=10, min.rho=0.1,
      edge.cex=1.6, cex.label=1.4,
      normalize.edges = TRUE,
      edge.type = "inter",
      value.name = ''
    )
    dev.off()
  }
  
  figs <- list(
    "Multi-layer LASAGNA network." = fig5,    
    "Module diagram (AI annotated)." = fig4,        
    "Graphical abstract (AI generated)." = fig3,
    "Dendrogram and module colors." = fig1,  
    "Module-trait correlation heatmap" = fig2
  )
  figs <- figs[!sapply(figs,is.null)]
  figs <- lapply(figs, function(f) paste0('/',f))
  labels=NULL
  
  ##------- tables -------------
  df1 <- data.frame(size = sapply(wgcna$me.genes,length))
  df1 <- df1[order(-df1$size),,drop=FALSE]

  topmodules <- names(wgcna$report$summaries)
  df2 <- wgcna.calculateSignificanceScore(wgcna, collapse=FALSE,
    sort.by="score", annot.cols = c("feature","symbol")) 
  df2 <- df2[topmodules]
  df2 <- lapply(df2, head, 10)
  names(df2) <- paste0("Top significance scores (",names(df2),")")
  for(i in 1:length(df2)) rownames(df2[[i]]) <- NULL

  div.tables <- c(
    list("WGCNA module sizes" = df1),
    df2
  )
  
  ##------- create sections -------------
  contents <- list(
    description = div.description,
    bullets = rpt$bullets,
    report = rpt$report,
    methods = div.methods,  
    settings = wgcna$settings,
    references = div.refs,
    figures = figs,
    tables = div.tables
  )
  txt <- rpt.compile_sections(contents, hlevel=hlevel, shift=shift)
    
  ## add title
  txt <- paste0("# ",title,"\n\n", txt)

  return(txt)
}

  
#'
#'
#' @export
rpt.compile_drugconnectivity_report <- function(obj, which.db = 1, report = NULL, 
                                                hlevel = 2, shift = TRUE,
                                                model = NULL, pgx=NULL,
                                                title = "Drug Connectivity Analysis") {

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

  rx=2  ## resolution parameter
  fig3 <- tempfile(fileext='.png')
  png(fig3, width=700*rx, height=300*rx, pointsize=14*rx)
  pgx.plotMOAactivationMap(pgx, type="drugClass", plotlib='base', rotate=TRUE)
  dev.off()

  fig4 <- tempfile(fileext='.png')
  png(fig4, width=700*rx, height=300*rx, pointsize=14*rx)
  pgx.plotMOAactivationMap(pgx, type="targetGene", plotlib='base', rotate=TRUE)
  dev.off()
  
  figs <- list(
    "The Activation Matrix visualizes the drug activation enrichment across the conditions. The size of the circles correspond to their relative enrichment, and are colored according to their positive (red) or negative (blue) correlation with the contrast profile." = fig1,
    "The Activation Matrix visualizes the drugClass MOA enrichment across the conditions. The size of the circles correspond to their relative enrichment, and are colored according to their positive (red) or negative (blue) correlation with the contrast profile." = fig3,
    "The Activation Matrix visualizes the MOA targetGene enrichment across the conditions. The size of the circles correspond to their relative enrichment, and are colored according to their positive (red) or negative (blue) correlation with the contrast profile." = fig4,
    "Infographic (AI generated)" = fig2
  )
  figs <- figs[!sapply(figs,is.null)]
  figs <- lapply(figs, function(f) paste0('/',f))
  
  ##------- tables -------------

  df1 <- pgx.getMOAmatrix(pgx, db=which.db, type="drugClass")
  df2 <- pgx.getMOAmatrix(pgx, db=which.db, type="targetGene")
  df1 <- round(df1, digits=3)
  df2 <- round(df2, digits=3)  
  
  order1 <- order(-rowMeans(df1**2))
  df1.top <- head(df1[order1,,drop=FALSE],15)

  order2 <- order(-rowMeans(df2**2))
  df2.top <- head(df2[order2,,drop=FALSE],20)
  
  div.tables <- list(
    "Top enriched MOA drug class" = df1.top,
    "Top enriched MOA target" = df2.top
  )

  settings <- list(
    database = which.db
  )
  
  ##------- create sections -------------
  contents <- list(
    description = div.description,
    bullets = rpt$bullets,
    report = rpt$report,
    methods = div.methods,  
    references = div.refs,
    settings = settings,
    figures = figs,
    tables = div.tables
  )
  txt <- rpt.compile_sections(contents, hlevel=hlevel, shift=shift)

  if(!is.null(title)) {
    title <- paste0(title," (",which.db,")") 
    txt <- paste0("# ",title,"\n\n", txt)
  }
  return(txt)
}

#' Summarize drug connectivity results
#'
#' @export
rpt.summarize_drug_connectivity <- function(pgx, ct, model, db=1, ntop=10) {
  
  if(is.numeric(db)) db <- names(pgx$drugs)[db]

  toplist <- list(
    "Top most similar (i.e. positively correlated) drugs are" =
      table_to_content(pgx.getTopDrugs(pgx, ct, n=ntop, dir=+1, db=db, na.rm=TRUE)), 
    "Top most inhibitory (i.e. negative correlated) drugs are:" = 
      table_to_content(pgx.getTopDrugs(pgx, ct, n=ntop, dir=-1, db=db, na.rm=TRUE)),
    "Top most positively enriched MOA classes are" =
      table_to_content(pgx.getTopMOA(pgx, ct, n=ntop, dir=+1, db=db, level=1)), 
    "Top most negatively (opposite) enriched MOA classes are" =
      table_to_content(pgx.getTopMOA(pgx, ct, n=ntop, dir=-1, db=db, level=1)), 
    "Top most positively enriched MOA drug target genes are" =
      table_to_content(pgx.getTopMOA(pgx, ct, n=ntop, dir=+1, db=db, level=2)), 
    "Top most negatively (opposite) enriched MOA drug target genes are" =
      table_to_content(pgx.getTopMOA(pgx, ct, n=ntop, dir=-1, db=db, level=2)), 
    "Top most positively enriched gene sets are" =
      paste(pgx.getTopGS(pgx, ct, n=50, dir=+1), collapse=';'), 
    "Top most negatively enriched gene sets are" =
      paste(pgx.getTopGS(pgx, ct, n=50, dir=-1), collapse=';') 
  )

  results=NULL
  if(grepl("sensitivity",db)) {
    results <- paste("Drug Synergy Analysis using Connectivity Map (CMap) analysis. Synergy of the mechanism of action (MOA) is based on correlation enrichment with computed drug sensitivity profiles of ",db," database. Positive correlation indicate possible synergy with the given drug. Negative correlation indicate possible antagonism with given drug. :\n\n",
    list_to_content(toplist, newline=TRUE), sep=""
    )
  } else {
    results <- paste("Drug Mechanism of Action. Drug Connectivity Map (CMap) analysis of selected comparison. Similarity of the mechanism of action (MOA) is based on correlation enrichment with drug perturbation profiles of ",db," database:\n\n",
      list_to_content(toplist, newline=TRUE), sep=""
    )
  }
  
  resp <- ai.ask( paste("Give a summary of the following results from a drug connectivity MOA analysis.
Give an integrated interpretation and a pharmacological narrative. Validate inferred drug MOA with the given (measured) enriched up/down gene sets. Do not include tables, be concise, write in prose. \n\nAnalysis results: ",results), model = model)

  resp2 <- ai.ask(paste("Give a short 3 bullet point summary of the following text. Focus on similarity with drugs. Use very short sentences, no titles. :",resp),  model = model)
  
  out <- list(
    bullet = resp2,
    summary = resp
  )
}
