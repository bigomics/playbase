##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


#' @title Compute drug enrichment from expression data
#'
#' @param X Numeric gene expression matrix
#' @param drugs Character vector of drug names
#' @param nmin Minimum number of genes per drug set
#' @param methods Methods for enrichment (rank, gsea)
#'
#' @return List with enrichment results for each drug
#'
#' @description Computes enrichment of drug gene sets in expression data using
#' rank correlation or GSEA.
#'
#' @details This function takes a gene expression matrix and a set of drugs. It extracts
#' gene sets for each drug from MSigDB. Enrichment is computed by:
#'
#' 1) Rank correlation of drug set ranks with experiment ranks.
#' 2) GSEA using the drug sets as gene sets.
#'
#' Drugs with fewer than nmin genes are filtered out. The output is a list
#' containing the enrichment results for each drug.
#'
#' @export
pgx.computeDrugEnrichment <- function(pgx, X = NULL, xdrugs = NULL,
                                      drug_info = NULL,
                                      methods = c("GSEA", "cor"), ## deprecated
                                      nmin = 15, nprune = 250, contrast = NULL) {
  ## 'pgx'   : can be ngs object or fold-change matrix
  ## X       : drugs profiles (may have multiple for one drug)
  ## xdrugs  : drug associated with profile

  if (!all(c("gx.meta", "genes") %in% names(pgx))) {
    stop("[pgx.computeDrugEnrichment] FATAL. not a pgx object")
  }

  if (is.null(X)) {
    X <- playdata::L1000_ACTIVITYS_N20D1011    
  }
  if(is.null(xdrugs)) {
    xdrugs <- gsub("[_@].*$", "", colnames(X))
  }

  FC <- pgx.getMetaMatrix(pgx)$fc
  FC <- FC[!duplicated(rownames(FC)), , drop = FALSE]
  FC <- collapse_by_humansymbol(FC, annot = pgx$genes)

  if (is.null(contrast)) {
    contrast <- colnames(FC)
  }
  contrast <- intersect(contrast, colnames(FC))
  FC <- FC[, contrast, drop = FALSE]

  ## test for overlap
  gg <- intersect(rownames(X), rownames(FC))
  if (length(gg) < 20) {
    message("WARNING::: pgx.computeDrugEnrichment : not enough common genes!!")
    return(NULL)
  }

  colnames(X) <- paste0(xdrugs,"_",colnames(X))
  enr <- metaLINCS::computeConnectivityEnrichment(
    FC, names = NULL, mDrugEnrich = X, 
    nmin = nmin, nprune = nprune) 
  
  ## Now compute the MoA enrichment (if not done)
  moa <- NULL
  if(!is.null(drug_info)) {
    moa <- metaLINCS::computeMoaEnrichment(
      enr, annot = drug_info)
  }

  results <- list()
  results$cor <- list(NULL)
  results$GSEA <- enr[c("X","Q","P","size")]
  results$stats <- enr$stats
  results$drug <- enr$drug
  results$moa <- moa
  
  message("[pgx.computeDrugEnrichment] done!")

  return(results)
}

#' Update older versions of pgx$drugs with MOA, report and infographic
#'
#' @export
pgx.update_drugs_results <- function(pgx, model, img_model) {

  if(is.null(pgx$drugs)) {
    return(pgx)
  }
  
  if(is.null(pgx$drugs[[1]]$moa)) {
    dbg("[pgx.update_drugs_results] updating MOA enrichment...")
    for(db in names(pgx$drugs)) {
      res <- pgx$drugs[[db]]
      pgx$drugs[[db]][['moa']] <- metaLINCS::computeMoaEnrichment(res) 
    }
  }
  
  if(is.null(pgx$drugs[[1]]$report) && !is.null(model)) {
    dbg("[pgx.update_drugs_results] updating reports...")
    for(db in names(pgx$drugs)) {
      pgx$drugs[[db]]$report <- cmap.create_report(
        pgx, model=model, model2=NULL, db=db, user.prompt = NULL)
    }
  }

  ## check infographic
  rpt <- pgx$drugs[[1]]$report
  if(is.null(rpt$infographic) && !is.null(img_model)) {
    dbg("[pgx.update_drugs_results] updating infographics...")    
    for(db in names(pgx$drugs)) {    
      rpt <- pgx$drugs[[db]]$report
      outfile <- ai.cmap_create_infographic(
        report = rpt$report,
        model = img_model,
        filename = "/tmp/drug-infographic.png",
        aspectratio = "16:9",
        add.fallback = TRUE
      )
      if(grepl("jpg$",outfile,ignore.case=TRUE)) {
        img <- jpeg::readJPEG(outfile)
      }
      if(grepl("png$",outfile,ignore.case=TRUE)) {
        img <- png::readPNG(outfile)
      }
      pgx$drugs[[db]]$report$infographic <- img
    }
  }

  return(pgx)
}


#' @export
pgx.getDrugConnectivityTable <- function(pgx, contrast, db=1,
                                         drugs = NULL,
                                         filter = FALSE)
{
  if(is.null(drugs)) {
    drugs <- pgx$drugs
  }
  dr <- drugs[[db]]
  nes <- dr$X[, contrast]
  pv <- dr$P[, contrast]
  qv <- dr$Q[, contrast]
  drug <- rownames(dr$X)
  if (is.null(ncol(dr$stats))) {
    stats <- dr$stats
  } else {
    stats <- dr$stats[, contrast]
  }
  annot <- dr$annot
  
  ## ## !!!SHOULD MAYBE BE DONE IN PREPROCESSING???
  ## if (is.null(annot)) {
  ##   warning("[getActiveDSEA] WARNING:: missing drug annotation in PGX file!")
  ##   annot <- read.csv(file.path(FILESX, "cmap/L1000_repurposing_drugs.txt"),
  ##     sep = "\t", comment.char = "#"
  ##   )
  ##   rownames(annot) <- annot$pert_iname
  ## }

  nes[is.na(nes)] <- 0
  qv[is.na(qv)] <- 1
  pv[is.na(pv)] <- 1
  
  ## ## compile results matrix
  jj <- match(toupper(drug), toupper(rownames(annot)))
  annot <- annot[jj, c("moa", "target")]
  dt <- data.frame(drug = drug, NES = nes, pval = pv, padj = qv, annot)
  dt <- dt[order(-dt$NES), ]
    
  if (filter) {
    sel <- which(dt$moa != "" | dt$target != "")
    dt <- dt[sel, , drop = FALSE]
  }

  dsea <- list(table = dt, clust = dr$clust, stats = stats)
  return(dsea)
}

#' Plot drug connectivity as enrichment plots
#'
#' @export
pgx.plotDrugConnectivity <- function(pgx, contrast,
                                     db = "L1000/activity",
                                     drugs = NULL,
                                     nplots = 25) {
  if (!"drugs" %in% names(pgx)) {
    stop("pgx does not have drug enrichment results")
  }

  if (!db %in% names(pgx$drugs)) {
    stop("pgx$drugs does not have database db = ", db)
  }

  res <- pgx$drugs[[db]]  
  metaLINCS::plotDrugConnectivity(
    res, contr=contrast, drugs=drugs, nplots=nplots)
}

#' Plot MOA enrichment as barplot
#'
#' @export
pgx.plotDrugMOA <- function(pgx, contrast,
                            db = "L1000/activity",
                            type=c("drugClass","targetGene")[1],
                            ntop = 20) {
  if (!"drugs" %in% names(pgx)) {
    stop("pgx does not have drug enrichment results")
  }

  if (!db %in% names(pgx$drugs)) {
    stop("pgx$drugs does not have database db = ", db)
  }

  if (!"moa" %in% names(pgx$drugs[[db]])) {
    stop("pgx$drugs missing MOA tables")
    moa <- metaLINCS::computeMoaEnrichment(pgx$drugs[[db]])     
  } else {
    moa <- pgx$drugs[[db]]$moa
  }
  metaLINCS::plotMOA(moa, contr=contrast, ntop=ntop, type=type)
}

#' Create annotation for drug combinations
#'
#' @title Create annotation for drug combinations
#'
#' @param combo Character vector of drug combinations (separated by '+')
#' @param annot0 Data frame with annotation for individual drugs
#'
#' @return Data frame with annotation for the drug combinations
#'
#' @description Creates annotation for drug combinations by combining annotation
#'   from individual drugs.
#'
#' @details This function takes a vector of drug combinations and an annotation
#'   data frame for individual drugs. It splits the combinations into individual
#'   drugs, looks up their annotation, and combines it into annotation for the combos.
#'
#'   For example, if combo is 'DrugA + DrugB', it will look up mode of action (moa)
#'   and targets for DrugA and DrugB, and combine them into moa and targets for
#'   the combination. NA values are removed.
#'
#' @export
pgx.createComboDrugAnnot <- function(combo, annot0) {
  drugs <- strsplit(combo, split = "[+]")
  drug1 <- sapply(drugs, "[", 1)
  drug2 <- sapply(drugs, "[", 2)
  j1 <- match(drug1, rownames(annot0))
  j2 <- match(drug2, rownames(annot0))
  cmoa <- paste(annot0[j1, "moa"], "+", annot0[j2, "moa"])
  ctarget <- paste(annot0[j1, "target"], "+", annot0[j2, "target"])
  cmoa <- gsub("NA [+] NA", "", cmoa)
  ctarget <- gsub("NA [+] NA", "", ctarget)
  annot1 <- data.frame(drug = combo, moa = cmoa, target = ctarget)
  rownames(annot1) <- combo
  return(annot1)
}


#' Returns MOA enrichment matrix from pgx object
#'
#' @export
pgx.getMOAmatrix <- function(pgx, db=1, type=c("drugClass","targetGene")[1]) {
  moa <- pgx$drugs[[db]]$moa
  if(is.null(moa)) {
    message("[pgx.getMOAmatrix] WARNING: pgx object has no MOA results")
    moa <- metaLINCS::computeMoaEnrichment(pgx$drugs[[db]]) 
  }
  mm <- lapply(moa, function(m) m[[type]])
  pp <- lapply(mm, function(m) m$pathway)
  pp <- Reduce(intersect, pp)
  mm <- lapply(mm, function(m) m[match(pp,m$pathway),])
  X <- sapply(mm, function(m) m$NES)
  rownames(X) <- pp
  X
}

#' Returns MOA enrichment matrix from pgx object
#'
#' @export
pgx.getDrugMOATable <- function(pgx, contrast, db=1, drugs=NULL,
                            type=c("drugClass","targetGene")[1]) {
  if(is.null(drugs)) {
    if(is.null(pgx$drugs)) {
      message("[pgx.getDrugMOATable] ERROR: pgx object has not drugs results")
      return(NULL)
    }
    drugs <- pgx$drugs
  }
  if(is.null(drugs[[db]]$moa)) {
    message("[pgx.getDrugMOATable] ERROR: pgx$drugs object has not MOA results")
    return(NULL)
  }

  moa <- drugs[[db]]$moa[[contrast]][[type]]
  moa <- moa[order(-abs(moa$NES)), ]
  return(moa)
}


#' @export
pgx.plotMOAactivationMap <- function(pgx,
                                     dmethod = 1,
                                     type=c("drugClass","targetGene")[1],
                                     contrast = NULL,
                                     nterms = 40,
                                     nfc = 20,
                                     normalize = FALSE,
                                     colorbar = FALSE,
                                     rotate = FALSE,
                                     plotlib = "base") {
  if (is.null(pgx$drugs)) {
    return(NULL)
  }

  ## get MOA matrix
  M <- pgx.getMOAmatrix(pgx, db=dmethod, type=type)
  res <- pgx$drugs[[dmethod]]
  res$X <- M
  res$Q <- M*0  

  plt <- NULL  
  if(plotlib == "plotly") {
    plt <- metaLINCS.plotlyActivationMap(
      res, contrast=contrast, nterms=nterms, nfc=nfc,
      normalize = normalize, drugs = pgx$drugs,
      colorbar = colorbar, rotate = rotate)
  } else {    
    metaLINCS::plotActivationMap(
      res, nterms=nterms, nfc=nfc, rot = !rotate)    
  }
  return(plt)
}


#' @export
pgx.plotDrugActivationMap <- function(pgx,
                                      dmethod = 1,
                                      contrast = NULL,
                                      nterms = 40,
                                      nfc = 20,
                                      normalize = FALSE,
                                      drugs = NULL,
                                      colorbar = FALSE,
                                      rotate = FALSE,
                                      plotlib = "plotly") {
  if (is.null(pgx$drugs)) {
    return(NULL)
  }
  res <- pgx$drugs[[dmethod]]

  if(plotlib == "plotly") {
    metaLINCS.plotlyActivationMap(
      res, contrast=contrast, nterms=nterms, nfc=nfc,
      normalize = normalize, drugs = drugs,
      colorbar = colorbar, rotate = rotate)
  } else {
    metaLINCS::plotActivationMap(
      res, nterms=nterms, nfc=nfc, rot = rotate)    
  }
}


## This should perhaps go to the metaLINCS package...
metaLINCS.plotlyActivationMap <- function(res, contrast = NULL,
                                          nterms=40, nfc=20,
                                          normalize = FALSE,
                                          drugs = NULL,
                                          colorbar = FALSE,
                                          rotate = FALSE) {   

  nes <- res$X
  qv <- res$Q  
  ctx <- colnames(nes)
  ctx <- sort(ctx[!grepl("^IA:", ctx)])
  nes <- nes[, ctx, drop = FALSE]
  qv <- qv[, ctx, drop = FALSE]
  
  score <- nes * (1 - qv)**2
  score[is.na(score)] <- 0

  ## sort by score
  if(is.null(contrast)) contrast <- colnames(score)
  ss <- rowMeans(score[, contrast,drop=FALSE]**2, na.rm=TRUE)
  score <- score[order(-ss), , drop = FALSE] 
    
  ## filter with table selection/search
  if (is.null(drugs)) {
    ii <- intersect(drugs, rownames(score))
    if(length(ii)) {
      score <- score[ii, , drop = FALSE]
    }
  }
  if (nrow(score) <= 1) {
    return(NULL)
  }
  
  score <- head(score, nterms) ## max number of terms
  score <- score[, head(order(-colSums(score**2)), nfc), drop = FALSE] ## max contr
  score <- score + 1e-3 * matrix(rnorm(length(score)), nrow(score), ncol(score))
  
  if (normalize) {
    score <- t(t(score) / (1e-8 + sqrt(colMeans(score**2))))
  }
  score <- sign(score) * abs(score)**3 ## fudging
  score <- score / (1e-8 + max(abs(score), na.rm = TRUE))
  
  if (NCOL(score) > 1) {
    d1 <- dist(score) ## euclidean
    d2 <- dist(t(score))
    d1[is.na(d1)] <- mean(d1,na.rm=TRUE)
    d2[is.na(d2)] <- mean(d2,na.rm=TRUE)
    ii <- hclust(d1)$order
    jj <- hclust(d2)$order
    score <- score[ii, jj, drop = FALSE]
  } else {
    score <- score[order(-score[, 1]), , drop = FALSE]
  }
  score <- score[nrow(score):1, , drop = FALSE]
  
  if(rotate) {
    score <- t(score)
  }

  colnames(score) <- substring(colnames(score), 1, 30)
  rownames(score) <- substring(rownames(score), 1, 50)
  color.scale <- list(
    list(0, "#3a5fcd"),
    list(0.25, "#ebeffa"),
    list(0.5, "white"),
    list(0.66, "#faeeee"),
    list(0.83, "#ebbbbb"),
    list(1, "#cd5555")
  )
  x_axis <- colnames(score)
  y_axis <- rownames(score)
  fig <- plotly::plot_ly(
    x = x_axis, y = y_axis,
    z = score, type = "heatmap",
    colorscale = color.scale,
    showscale = colorbar,
    zmin = -1, zmax = 1
  )
  return(fig)
}


## ======================================================================
## ======================= AI REPORT ====================================
## ======================================================================


#' Summarize drug connectivity results
#'
cmap.summarize_results <- function(pgx, ct, model, drugs=NULL, db=1,
                                   ntop=10, ntop2=50) {

  if(is.null(drugs)) drugs <- pgx$drugs
  if(is.numeric(db)) db <- names(drugs)[db]

  toplist <- list(
    "Top most positively enriched MOA classes are" =
      table_to_content(pgx.getTopMOA(pgx, ct, n=ntop, dir=+1, db=db, level=1)), 
    "Top most negatively (opposite) enriched MOA classes are" =
      table_to_content(pgx.getTopMOA(pgx, ct, n=ntop, dir=-1, db=db, level=1)), 
    "Top most positively enriched MOA drug target genes are" =
      table_to_content(pgx.getTopMOA(pgx, ct, n=ntop, dir=+1, db=db, level=2)), 
    "Top most negatively (opposite) enriched MOA drug target genes are" =
      table_to_content(pgx.getTopMOA(pgx, ct, n=ntop, dir=-1, db=db, level=2)), 

    "Top most similar (i.e. positively correlated) drugs are" =
      table_to_content(pgx.getTopDrugs(pgx, ct, n=ntop2, dir=+1, db=db, na.rm=TRUE)), 
    "Top most inhibitory (i.e. negative correlated) drugs are:" = 
      table_to_content(pgx.getTopDrugs(pgx, ct, n=ntop2, dir=-1, db=db, na.rm=TRUE)),
    "Top most positively enriched gene sets are" =
      paste(pgx.getTopGS(pgx, ct, n=ntop2, dir=+1), collapse=';'), 
    "Top most negatively enriched gene sets are" =
      paste(pgx.getTopGS(pgx, ct, n=ntop2, dir=-1), collapse=';') 
  )

  results=NULL
  if(grepl("sensitivity",db)) {
    results <- paste("Drug Synergy Analysis using Connectivity Map (CMap) analysis. Synergy of the mechanism of action (MOA) is based on correlation enrichment with computed drug sensitivity profiles of ",db," database. Positive correlation indicate possible synergy with the given drug. Negative correlation indicate possible antagonism with given drug. Indicate clearly that this is not a regular drug similarity analysis but a drug synergy/resistance analysis using CMap.\n\n**Top tables**:",
    list_to_content(toplist, newline=TRUE), sep=""
    )
  } else {
    results <- paste("Drug Mechanism of Action. Drug Connectivity Map (CMap) analysis of selected comparison. Similarity of the mechanism of action (MOA) is based on correlation enrichment with drug perturbation profiles of ",db," database:",
      "\n\n**Top tables**:", list_to_content(toplist, newline=TRUE), sep=""
    )
  }

  prompt <- paste("**Instructions**: Give an overall summary of the following results from a drug connectivity MOA analysis. Create an integrated pharmacological narrative focussing on inferred mechanisms-of-action class and possible drug targets. Validate inferred drug MOA with the given (measured) enriched up/down gene sets. Do not describe gene sets on its own, only in connection with drug pharmacological MOA. Do not include tables, be concise, write in prose is preferred. Format in markdown with title and sections.")
  prompt <- paste(prompt, "\n\n**Analysis results**: ",results)
  summary <- ai.ask(prompt, model = model)

  bullets <- ai.ask(paste("Give a short 2-3 bullet point summary of the following text. Focus on similarity with drugs MOA. Use very short sentences, no titles. :",summary),  model = model)
  
  out <- list(
    prompt = prompt,
    bullets = bullets,
    summary = summary
  )
}

#'
#'
#' @export
cmap.create_report <- function(pgx, model, model2=NULL, db=1,
                               user.prompt=NULL, force=FALSE) {

  if(is.null(pgx$drugs)) return(NULL)  
  if(is.null(model2)) model2=model

  db.name <- db
  if(is.numeric(db)) db.name <- names(pgx$drugs)[db]
  db.name

  rpt <- pgx$drugs[[db.name]]$report
  if(!is.null(rpt) && !force) {
    message("[cmap.create_report] object already has report")
    return(rpt)
  }
  
  ## check if moa results are there
  has.moa <- !is.null(pgx$drugs[[db]]$moa)
  if(!has.moa) {
    res <- pgx$drugs[[db]]
    pgx$drugs[[db]][['moa']] <- metaLINCS::computeMoaEnrichment(res) 
  }
  
  ## Create summaries for all contrasts
  summaries <- list()
  ct=colnames(pgx$contrasts)[1]
  for(ct in colnames(pgx$contrasts)) {
    message("summarizing drug connectivity results for ",ct,"...")
    summaries[[ct]] <- cmap.summarize_results(pgx, ct, model,
      drugs=NULL, db=db, ntop=10, ntop2=50)   
  }
  names(summaries)
  
  ## Create single summary from all summaries
  S <- lapply(summaries, function(x) x$summary)
  S <- collate_as_sections(S)

  ## Ask LLM
  prompt <- paste("**Instructions**: Below are drug connectivity analysis results for different comparisons. Summarize this into a single report. Minimize use of tables, writing in prose is preferred. Use scientific style. Format in markdown with title and sections.")
  if(!is.null(user.prompt)) prompt <- paste(prompt, user.prompt)
  prompt <- paste0(prompt, "\nExperiment description: ", pgx$description)
  prompt <- paste0(prompt, "\nDatabase: ", db.name)
  prompt <- paste0(prompt, "\n\n**CMap analysis**: ", S)
  rpt <- ai.ask(prompt, model = model2)
  
  ## create bullets
  prompt2 <- paste0("**Instructions**: Extract 3 short bullet points for the following drug connectivity report. Return only the list items.")
  prompt2 <- paste0(prompt2, "\n\n**Report**: ", rpt)
  bullets <- ai.ask(prompt2, model = model)

  out <- list(
    prompt = prompt,
    bullets = bullets,
    summaries = summaries,
    report = rpt,
    database = db.name,
    llm_model = model,
    llm_model2 = model2
  )
  return(out)
}


#' Create an infographic from a drug connectivity report.
#'
#' @export
ai.cmap_create_infographic <- function(report, model, filename,
                                       aspectratio = c("4:3","16:9","3:4")[2],
                                       add.fallback = TRUE)
{
  prompt <- paste("**Instructions**: Create an infographic for the following pharmacological mechanism-of-action analysis report.\n\n**summary**:", report)
  if(is.null(model) || model=="") {
    model <- ai.get_image_models()
  }
  if(add.fallback) {
    ## add fallback models
    model <- unique(c(model,playbase::ai.get_image_models()))  
  }
  out <- ai.create_image(
    prompt,
    model = model,
    format = "file",
    filename = filename,
    aspect_ratio = aspectratio,
    size = 1024
  )
  if(is.null(out) || !file.exists(out)) return(NULL)
  return(out)
}

