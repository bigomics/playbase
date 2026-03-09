##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##


#' Summarize drug connectivity results
#'
#' @export
ai.summarize_drug_connectivity <- function(pgx, ct, model, db=1, ntop=10) {
  
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
  
  resp <- ai.ask( paste("Give a summary of the following results from a drug connectivity MOA analysis. Give an integrated interpretation and a pharmacological narrative. Validate inferred drug MOA with the given (measured) enriched up/down gene sets. Do not include tables, be concise, write in prose. \n\nAnalysis results: ",results), model = model)

  resp2 <- ai.ask(paste("Give a short 3 bullet point summary of the following text. Focus on similarity with drugs. Use very short sentences, no titles. :",resp),  model = model)
  
  out <- list(
    bullet = resp2,
    summary = resp
  )
}


rpt.compile_wgcna_section <- function(pgx,
                                      llm_model,
                                      img_model = "gemini-3.1-flash-image-preview") {

  wgcna <- pgx$wgcna

  if(!is.null(wgcna$report)) {
    message("using embedded WGCNA report")
    wgcna.rpt <- wgcna$report
  } else {
    message("creating WGCNA report...")
    wgcna.rpt <- wgcna.create_report(
      wgcna, ai_model=llm_model, graph=NULL, annot=pgx$genes, multi=FALSE,
      ntop = 100, topratio = 0.85, psig = 0.05, do.diagram = FALSE,
      userprompt='', format="markdown", verbose=1, progress=NULL)
  }
  
  if(!is.null(wgcna.rpt$infographic)) {
    ## If infographic is stored in report use this
    message("found embedded infographic")
  } else {
    message("creating new infographic...")
    tmpf <- tempfile(fileext='.png')
    tmpf <- wgcna.create_infographic(
      wgcna.rpt$report,
      diagram = wgcna.rpt$diagram,
      prompt = NULL,
      model = img_model,
      add.fallback = FALSE,
      filename = tmpf)
    wgcna.rpt$infographic <- png::readPNG(tmpf)
  }

  
  ##------- description -------------
  div.description <- 
    "**Weighted Gene Co-expression Network Analysis (WGCNA) is a systems biology method for finding clusters (modules) of highly correlated genes in high-dimensional data, such as RNA-seq or microarray samples. It identifies modules based on similar expression patterns, relates them to external sample traits, and finds key hub genes, enabling the transition from simple gene l
ists to functional biological insights.**"
  
  div.methods <-
    "WGCNA computation was done in Omics Playground (BigOmics Analytics, Switzerland) using the `WGCNAplus` R package. The optimal soft threshold β was determined by optimizing the IQR of the dendrogram heights. The adjacency matrix was transformed into a topological overlapping matrix (TOM) and used for hierarchical clustering of the features. Subsequently, the dynamic tree‐cutting algorithm was employed to identify gene modules, each comprising a minimum of genes. Similar modules were merged using a shear height threshold. Finally, modules were distinguished by different colours and represented by their module eigengene (ME). Module-trait relationship was determined based on the Pearson correlation between ME and clinical features. Gene significance was assessed by computing the Module Membership (MM) (correlation coefficient between gene significance and module significance) as well as Gene Significance (GS) values (correlation coefficient between gene expression levels and clinical characteristics)."
  
  div.refs <- c(
    "Langfelder, P., Horvath, S. WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics 9, 559 (2008). https://doi.org/10.1186/1471-2105-9-559",
    "Langfelder P, Luo R, Oldham MC, Horvath S (2011) Is My Network Module Preserved and Reproducible? PLoS Comput Biol 7(1): e1001057. https://doi.org/10.1371/journal.pcbi.1001057",
  "WGCNAplus R package: `https://github.com/bigomics/WGCNAplus`"
)
  
  ##------- figures -------------
  rx=2
  fig1 <- tempfile(fileext='.png')
  png(fig1, w=800*rx, h=400*rx, pointsize=14*rx)
  wgcna.plotDendroAndColors(wgcna, show.traits=1,
    main="Cluster dendrogram and module colors")
  dev.off()
  
  fig2 <- tempfile(fileext='.png')
  tr <- (ncol(wgcna$modTraits) < nrow(wgcna$modTraits))
  png(fig2, w=800*rx, h=400*rx, pointsize=14*rx)
  wgcna.plotModuleTraitHeatmap(wgcna, cluster=TRUE, text=FALSE,
    transpose=tr, main="Module-trait correlation")
  dev.off()
  
  fig3 <- tempfile(fileext='.png')
  img <- wgcna.rpt$infographic
  png::writePNG(img, target=fig3)    
  
  figs <- list(
    "Infographic (AI generated)." = fig3,
    "Dendrogram and module colors." = fig1,  
    "Module-trait correlation heatmap" = fig2  
  )
  figs <- lapply(figs, function(f) paste0('/',f))
  labels=NULL
  
  ##------- tables -------------
  df <- data.frame(size = sapply(wgcna$me.genes,length))
  
  ts <- wgcna$stats$traitSignificance
  mm <- wgcna$stats$moduleMembership
  mod <- wgcna$net$labels[rownames(mm)]
  max.ts <- apply(abs(ts),1,max,na.rm=TRUE)
  mod.mm <- mm[cbind(1:nrow(mm),match(mod,colnames(mm)))]
  gsig <- (max.ts * mod.mm)
  df2 <- data.frame(module=mod, MM=mod.mm, max.TS=max.ts, importance=gsig)
  df2 <- df2[order(-df2$importance),]
  df2[,2:ncol(df2)] <- round(df2[,2:ncol(df2)], digits=3)
  
  div.tables <- list(
    "WGCNA module sizes" = df,
    "Feature importance" = head(df2,40)
  )
  
  ##------- create sections -------------
  contents <- list(
    description = div.description,
    bullets = wgcna.rpt$bullets,
    report = wgcna.rpt$report,
    methods = div.methods,  
    settings = wgcna$settings,
    references = div.refs,
    figures = figs,
    tables = div.tables
  )
  md.render_section(contents)
}
