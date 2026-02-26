##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##


ai.create_report_drug_connectivity <- function(pgx, model, db=1) {
  
}



#' Summarize drug connectivity results
#'
#' @export
ai.summarize_drug_connectivity <- function(pgx, ct, model, drugs=NULL, db=1, ntop=10) {

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
      table_to_content(pgx.getTopDrugs(pgx, ct, n=ntop, dir=+1, db=db, na.rm=TRUE)), 
    "Top most inhibitory (i.e. negative correlated) drugs are:" = 
      table_to_content(pgx.getTopDrugs(pgx, ct, n=ntop, dir=-1, db=db, na.rm=TRUE)),
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

  prompt <- paste("**Instructions**: Give a summary of the following results from a drug connectivity MOA analysis. Create an integrated interpretation and a pharmacological narrative. Validate inferred drug MOA with the given (measured) enriched up/down gene sets. Do not describe gene sets on its own, only in connection with drug pharmacological MOA. Do not include tables, be concise, write in prose. \n\n**Analysis results**: ",results)
  resp1 <- ai.ask(prompt, model = model)

  resp2 <- ai.ask(paste("Give a short 2-3 bullet point summary of the following text. Focus on similarity with drugs MOA. Use very short sentences, no titles. :",resp1),  model = model)
  
  out <- list(
    prompt = prompt,
    bullets = resp2,
    summary = resp1
  )
}
