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
                            ntop = ntop) {
  if (!"drugs" %in% names(pgx)) {
    stop("pgx does not have drug enrichment results")
  }

  if (!db %in% names(pgx$drugs)) {
    stop("pgx$drugs does not have database db = ", db)
  }

  if (!"moa" %in% names(pgx$drugs[[db]])) {
    stop("pgx$drugs missing MOA tables")
  }

  moa <- pgx$drugs[[db]]$moa
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
  if(is.null(pgx$drugs[[1]]$moa)) {
    message("[pgx.getMOAmatrix] ERROR: pgx object has not MOA results")
    return(NULL)
  }
  moa <- pgx$drugs[[db]][['moa']] 
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
pgx.plotDrugActivationMap <- function(pgx,
                                      dmethod = 1,
                                      contrast = NULL,
                                      nterms = 40,
                                      nfc = 20,
                                      normalize = FALSE,
                                      drugs = NULL,
                                      colorbar = FALSE,
                                      rotate = FALSE) {
  if (is.null(pgx$drugs)) {
    return(NULL)
  }
  res <- pgx$drugs[[dmethod]]

  metaLINCS.plotlyActivationMap(
    res, nterms=nterms, nfc=nfc,
    normalize = normalize, drugs = drugs,
    colorbar = colorbar, rotate = rotate)


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
ai.summarize_drug_connectivity <- function(pgx, ct, model, drugs=NULL, db=1,
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
ai.create_report_drug_connectivity <- function(pgx, model, model2=NULL, db=1,
                                               user.prompt=NULL) {

  if(is.null(model2)) model2=model

  db.name <- db
  if(is.numeric(db)) db.name <- names(pgx$drugs)[db]
  db.name
  
  ## Create summaries for all contrasts
  summaries <- list()
  for(ct in colnames(pgx$contrasts)) {
    summaries[[ct]] <- ai.summarize_drug_connectivity(pgx, ct, model,
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
  rpt

  ## create bullets
  B <- lapply(summaries, function(x) x$bullets)
  B <- collate_as_sections(B)
  prompt2 <- paste0("**Instructions**: Below are bullet points for drug connectivity analysis results for different comparisons. Extract 3 short bullet points that are consistent across all. Keep short sentences. Give only the list items.")
  prompt2 <- paste0(prompt2, "\n\n**Bullet points**: ", B)
  bullets <- ai.ask(prompt2, model = model)

  out <- list(
    prompt = prompt,
    bullets = bullets,
    report = rpt
  )
  return(out)
}


#' Create an infographic from a drug connectivity report.
#'
#' @export
ai.cmap_create_infographic <- function(report, model, filename,
                                       aspectRatio = c("4:3","16:9","3:4")[2],
                                       add.fallback = TRUE)
{
  prompt <- paste("**Instructions**: Create an infographic for the following pharmacological mechanism-of-action analysis report.\n\n**summary**:", report)
  out <- playbase::ai.create_image_gemini(
    prompt,
    model = model,
    api_key = Sys.getenv("GEMINI_API_KEY"),
    format = "file",
    filename = filename,
    aspectRatio = aspectRatio,
    imageSize = "1K"
  )
  if(is.null(out) || !file.exists(out)) return(NULL)
  return(out)
}



## ======================================================================
## ==================== DEPRECATED ======================================
## ======================================================================

pgx.computeDrugEnrichment.BAK <- function(pgx, X, xdrugs, drug_info = NULL,
                                          methods = c("GSEA", "cor"),
                                          nmin = 15, nprune = 250, contrast = NULL) {
  ## 'pgx'   : can be ngs object or fold-change matrix
  ## X       : drugs profiles (may have multiple for one drug)
  ## xdrugs  : drug associated with profile

  if (!all(c("gx.meta", "genes") %in% names(pgx))) {
    stop("[pgx.computeDrugEnrichment] FATAL. not a pgx object")
  }

  if (is.null(X)) {
    X <- playdata::L1000_ACTIVITYS_N20D1011
    dim(X)
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

  ## create drug meta sets
  meta.gmt <- tapply(colnames(X), xdrugs, list)
  meta.gmt <- meta.gmt[which(sapply(meta.gmt, length) >= nmin)]
  if (length(meta.gmt) == 0) {
    message("WARNING::: pgx.computeDrugEnrichment : no valid genesets!!")
    return(NULL)
  }

  ## gene set enrichment by rank correlation. If X is sparse we can
  ## use the faster corSparse().
  message("Calculating first level rank correlation ...")
  if (any(class(X) == "dgCMatrix")) {
    fx <- apply(FC[gg, , drop = FALSE], 2, rank)
    R1 <- qlcMatrix::corSparse(X[gg, ], fx)
  } else {
    rnk1 <- t(matrixStats::colRanks(X[gg, , drop = FALSE], ties.method = "average", na.last = "keep"))
    rnk2 <- t(matrixStats::colRanks(FC[gg, , drop = FALSE], ties.method = "average", na.last = "keep"))
    system.time(R1 <- stats::cor(rnk1, rnk2, use = "pairwise"))
  }
  R1 <- as.matrix(R1)
  R1[is.nan(R1)] <- 0
  R1[is.infinite(R1)] <- 0
  R1 <- R1 + 1e-8 * matrix(stats::rnorm(length(R1)), nrow(R1), ncol(R1))
  colnames(R1) <- colnames(FC)
  rownames(R1) <- colnames(X)

  ## experiment to drug
  results <- list()
  if ("cor" %in% methods) {
    message("Calculating drug enrichment using rank correlation ...")
    D <- Matrix::sparse.model.matrix(~ 0 + xdrugs)
    colnames(D) <- sub("^xdrugs", "", colnames(D))
    rownames(D) <- colnames(X) ## not necessary..
    rho2 <- qlcMatrix::corSparse(D, R1)
    rownames(rho2) <- colnames(D)
    colnames(rho2) <- colnames(R1)
    rho2 <- rho2[order(-rowMeans(rho2**2, na.rm = TRUE)), , drop = FALSE]
    .cor.pvalue <- function(x, n) 2 * stats::pnorm(-abs(x / ((1 - x**2) / (n - 2))**0.5))
    P <- apply(rho2, 2, .cor.pvalue, n = nrow(D))
    Q <- apply(P, 2, stats::p.adjust, method = "fdr")
    results[["cor"]] <- list(X = rho2, Q = Q, P = P)
  }

  if ("GSEA" %in% methods) {
    bpparam <- BiocParallel::MulticoreParam(1)
    message("Calculating drug enrichment using GSEA ...")
    res0 <- list()
    i <- 1
    for (i in 1:ncol(R1)) {
      suppressWarnings(res0[[i]] <- fgsea::fgseaSimple(meta.gmt, stats = R1[, i], nperm = 10000, BPPARAM = bpparam))
    }
    names(res0) <- colnames(R1)

    mNES <- sapply(res0, function(x) x$NES)
    mQ <- sapply(res0, function(x) x$padj)
    mP <- sapply(res0, function(x) x$pval)
    if (length(res0) == 1) {
      mNES <- cbind(mNES)
      mP <- cbind(mP)
      mQ <- cbind(mQ)
    }

    pw <- res0[[1]]$pathway
    rownames(mNES) <- rownames(mQ) <- rownames(mP) <- pw
    colnames(mNES) <- colnames(mQ) <- colnames(mP) <- colnames(FC)
    msize <- res0[[1]]$size
    results[["GSEA"]] <- list(X = mNES, Q = mQ, P = mP, size = msize)
  }

  ## this takes only the top matching drugs for each comparison to
  ## reduce the size of the matrices
  if (nprune > 0) {
    message("[pgx.computeDrugEnrichment] pruning : nprune = ", nprune)
    k <- 1
    for (k in 1:length(results)) {
      res <- results[[k]]
      ## reduce solution set with top-N of each comparison??

      rx <- apply(abs(res$X), 2, rank)
      rownames(rx) <- rownames(res$X)
      mtop <- names(utils::head(sort(rowMeans(rx, na.rm = TRUE), decreasing = TRUE), nprune))
      top.idx <- unique(unlist(meta.gmt[mtop]))

      results[[k]]$X <- res$X[mtop, , drop = FALSE]
      results[[k]]$P <- res$P[mtop, , drop = FALSE]
      results[[k]]$Q <- res$Q[mtop, , drop = FALSE]
      results[[k]]$size <- res$size[mtop]
    }
  }

  ## reduce large stats object
  sel.drugs <- unique(unlist(sapply(results, function(res) rownames(res$X))))
  sel <- which(xdrugs %in% sel.drugs)
  results$stats <- R1[sel, , drop = FALSE]

  message("[pgx.computeDrugEnrichment] done!")

  return(results)
}

pgx.plotDrugConnectivity.BAK <- function(pgx, contrast,
                                         db = "L1000/activity",
                                         moatype = c("class", "gene", "drug")[1],
                                         ntop = 10) {
  if (!"drugs" %in% names(pgx)) {
    stop("pgx does not have drug enrichment results")
  }

  if (!db %in% names(pgx$drugs)) {
    stop("pgx$drugs does not have database db = ", db)
  }

  ## common getData-esque function for drug connectivity plots / tables
  getActiveDSEA <- function(pgx, contrast, db) {
    dr <- pgx$drugs[[db]]
    nes <- round(dr$X[, contrast], 4)
    pv <- round(dr$P[, contrast], 4)
    qv <- round(dr$Q[, contrast], 4)
    drug <- rownames(dr$X)
    stats <- dr$stats[, contrast]
    annot <- dr$annot
    nes[is.na(nes)] <- 0
    qv[is.na(qv)] <- 1
    pv[is.na(pv)] <- 1

    ## !!!SHOULD MAYBE BE DONE IN PREPROCESSING???
    if (is.null(annot)) {
      stop("[getActiveDSEA] WARNING:: missing drug annotation in PGX file!")
    }

    ## compile results matrix
    jj <- match(toupper(drug), toupper(rownames(annot)))
    annot <- annot[jj, c("moa", "target")]
    dt <- data.frame(drug = drug, NES = nes, pval = pv, padj = qv, annot)
    dt <- dt[order(-dt$NES), ]

    filter_empty <- FALSE
    if (filter_empty) {
      sel <- which(dt$moa != "" | dt$target != "")
      dt <- dt[sel, , drop = FALSE]
    }
    dsea <- list(table = dt, clust = dr$clust, stats = stats)
    return(dsea)
  }

  getMOA.target <- function(dsea) {
    ## meta-GSEA on molecular targets
    dt <- dsea$table
    targets.list <- lapply(
      enc2utf8(as.character(dt$target)),
      function(s) trimws(strsplit(s, split = "[\\|;,]")[[1]])
    )
    names(targets.list) <- rownames(dt)
    targets <- setdiff(unique(unlist(targets.list)), c(NA, "", " "))
    gmt <- lapply(targets, function(g) {
      names(which(sapply(targets.list, function(t) (g %in% t))))
    })
    names(gmt) <- targets

    rnk <- dt$NES
    names(rnk) <- rownames(dt)
    suppressWarnings(
      moa.target <- fgsea::fgsea(gmt, rnk, nperm = 20000)
    )
    moa.target <- moa.target[order(-abs(moa.target$NES)), ]
    return(moa.target)
  }

  getMOA.class <- function(dsea) {
    ## meta-GSEA on MOA terms
    dt <- dsea$table
    moa.list <- lapply(
      enc2utf8(as.character(dt$moa)),
      function(s) trimws(strsplit(s, split = "[\\|;,]")[[1]])
    )
    names(moa.list) <- rownames(dt)
    moa <- setdiff(unlist(moa.list), c("", NA, " "))
    gmt <- lapply(moa, function(g) names(which(sapply(moa.list, function(t) (g %in% t)))))
    names(gmt) <- moa
    rnk <- dt$NES
    names(rnk) <- rownames(dt)
    suppressWarnings(
      moa.class <- fgsea::fgsea(gmt, rnk)
    )
    moa.class <- moa.class[order(-abs(moa.class$NES)), ]
    return(moa.class)
  }

  plotTopEnriched <- function(res, ntop) {
    res$score <- res$NES

    qweight <- FALSE
    if (qweight) {
      res$score <- res$NES * (1 - res$padj) * (1 - 1 / res$size**1)
      yaxistitle <- "score (qNES)"
    }
    jj <- unique(c(head(order(-res$score), ntop), tail(order(-res$score), ntop)))
    moa.top <- res$score[jj]
    names(moa.top) <- res$pathway[jj]

    df <- data.frame(
      x = factor(names(moa.top), levels = names(moa.top)),
      y = as.numeric(moa.top)
    )

    barplot(
      height = df$y,
      names.arg = df$x,
      ylab = "score (NES)",
      xlab = "",
      las = 3,
      ylim = c(-1.1, 1.1) * max(abs(as.numeric(moa.top)))
    )
  }

  plotTopDrugs <- function(db, ntop = 15) {
    dr <- pgx$drugs[[db]]
    sel <- 1:nrow(dr$annot)
    dd <- rownames(dr$X)
    ##  sel <- grepl("[a-z]{4}", dd) & !is.na(dr$annot[dd, "moa"])
    sel <- grepl("[a-z]{4}", dd)
    dx <- sort(dr$X[sel, 1], decreasing = TRUE)
    dx.top <- c(head(dx, ntop), tail(dx, ntop))
    barplot(dx.top, las = 3, horiz = FALSE, ylab = "similarity (NES)")
  }

  ## this should move to pgx.computeDrugEnrichment...
  dsea <- getActiveDSEA(pgx, contrast, db)
  if (moatype == "gene") {
    res <- getMOA.target(dsea)
    plotTopEnriched(res, ntop = ntop)
  } else if (moatype == "class") {
    res <- getMOA.class(dsea)
    plotTopEnriched(res, ntop = ntop)
  } else if (moatype == "drug") {
    plotTopDrugs(db, ntop = ntop)
  }
}
