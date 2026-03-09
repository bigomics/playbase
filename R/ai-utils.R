##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

table_to_content <- function(df, caption=NULL, label=NULL) {
  if(is.null(df) || nrow(df)==0) return("<empty>")
  rownames(df) <- iconv2ascii(rownames(df))
  for(i in seq_len(ncol(df))) {
    if(is.character(df[[i]])) df[[i]] <- iconv2ascii(df[[i]])
  }
  tbl <- paste(as.character(knitr::kable(df,format="markdown")),collapse="\n")
  if(!is.null(caption)) tbl <- paste0(tbl,"\n\n: ",caption)
  if(!is.null(label)) tbl <- paste0(tbl," {#tbl-",label,"}")
  tbl <- paste0(tbl,"\n")
  tbl
}

list_to_content <- function(a, newline=FALSE) {
  if(is.null(a) || length(a)==0) return("<empty>")
  aa <- sapply(a, function(s) paste(unlist(s), collapse='; '))
  pp <- paste0("- **",names(a),"**: ")
  if(newline) pp <- paste0(pp, '\n\n')
  cc <- paste(paste0(pp,aa), collapse='\n')
  paste(cc,'\n')
}

collate_as_sections <- function(a, level=2, csep=NULL, hdr=NULL) {
  if(is.null(a) || length(a)==0) return("<empty>")
  hh <- paste(rep("#",level),collapse="")  
  for(i in 1:length(a)) {
    tt <- names(a)[i]
    if(!is.null(hdr)) tt <- paste0(hdr,": ",tt)
    a[[i]] <- paste0(hh," ",tt,"\n\n",a[[i]])
  }
  if(!is.null(csep)) a <- paste0(a,"\n",csep,"\n")
  a <- paste0(a, collapse="\n")
  return(a)
}

list_to_table <- function(s, name='value') {
  df <- data.frame(unlist(s))
  colnames(df)[1] <- name
  table_to_content(df) 
}

list_to_list <- function(s, type='ul', add.name=TRUE) {
  if(add.name && !is.null(names(s))) s <- paste(names(s),"=",s)
  if(type == 'ul') {
    s <- paste(paste0('- ',s),collapse='\n')
  } else {
    s <- paste(paste0('1. ',s),collapse='\n')
  }
  s
}

md.list_to_figs <- function(figs, labels=NULL) {
  captions <- names(figs)
  paths <- as.character(figs)
  if(is.null(labels)) labels <- 1:length(figs)
  aa <- paste0("![",captions,"](",paths,")")
  bb <- ""
  if(!is.null(labels)) bb <- paste0("#fig-",labels)
  #bb <- paste(bb, "fig-pos='h!'")
  aa <- paste0(aa,"{",bb,"}")
  aa <- paste0(paste(aa, collapse='\n\n'),"\n")
  aa
}

#' Export markdown to PDF 
#' 
#' @export
markdownToPDF <- function(text, file, tmpdir=NULL, engine='latex') {

  text <- gsub(intToUtf8("8209"),"-",text)
  ##text <- gsub("---\n","",text)

  ## header/frontmatter
  hdr <- paste0("---\n")
#  hdr <- paste0(hdr, "title: TITLE\n")
#  hdr <- paste0(hdr, "title-block-banner: true\n")  
  hdr <- paste0(hdr, "format:\n")
  if(engine=='latex') {
    hdr <- paste0(hdr, "  pdf:\n")
    hdr <- paste0(hdr, "    pdf-engine: lualatex\n")
    #hdr <- paste0(hdr, "    documentclass: report\n")  
    hdr <- paste0(hdr, "    documentclass: article\n")
    hdr <- paste0(hdr, "    papersize: a4\n")
    hdr <- paste0(hdr, "    geometry:\n")
    hdr <- paste0(hdr, "      - left=25mm\n")
    hdr <- paste0(hdr, "      - right=20mm\n")
    hdr <- paste0(hdr, "      - top=25mm\n")
    hdr <- paste0(hdr, "      - bottom=25mm\n")
    hdr <- paste0(hdr, "    mainfont: Lato\n")
    hdr <- paste0(hdr, "    fig-pos: 'h!'\n")    
  }
  if(engine=='typst') {
    hdr <- paste0(hdr, "  typst:\n")
    hdr <- paste0(hdr, "    papersize: a4\n")
    hdr <- paste0(hdr, "    margin:\n")
    hdr <- paste0(hdr, "      left: 10mm\n")
    hdr <- paste0(hdr, "      right: 0mm\n")
    #hdr <- paste0(hdr, "    mainfont: Lato\n")
    hdr <- paste0(hdr, "    mainfont: helvetica\n")
    hdr <- paste0(hdr, "    fontsize: 10pt\n")    
  }
  hdr <- paste0(hdr, "---\n\n")
  text <- paste0(hdr, text)

  curwd <- getwd()
  if(is.null(tmpdir)) tmpdir <- tempdir()
  setwd(tmpdir)
  md.file <- file.path(tmpdir,"report.md")
  write(text, file=md.file)
  quarto::quarto_render(md.file, output_format="pdf", quiet=TRUE)
  pdf.file <- file.path(tmpdir,"report.pdf")
  setwd(curwd)
  file.copy(pdf.file, file, overwrite=TRUE)
  return(file)
}
