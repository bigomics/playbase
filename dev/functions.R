##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

scan_packages <- function(path='R') {
  
  ## ---------------------------------------------------------------------
  ## Automatically scan all used packages and install
  ## ---------------------------------------------------------------------
  ## We use renv to detect dependencies. Renv is looking for library and
  ## require statements in the r/R source files.
  renv.out <- renv::dependencies(path = path, root = getwd(), errors = "ignored")
  pkg.used <- sort(unique(renv.out$Package))
  pkg.used <- setdiff(pkg.used, "playbase")

  ## Define remote locations or versions
  github_url <- function(repo) paste0("github::",repo)
  github_url <- function(repo) {
    if(grepl("@",repo)) {
      branch <- sub(".*@","",repo)
      repo <- sub("@.*","",repo)    
      paste0("url::https://github.com/",repo,"/archive/refs/heads/",branch,".zip")
    } else {
      paste0("url::https://github.com/",repo,"/archive/HEAD.zip")
    }
  }

  add_github <- function(repo) {
    pkg.name <- gsub(".*[/]|@.*","",repo)
    remotes.url[pkg.name] <<- github_url(repo)
  }
  remotes.url <- c(
    "KEGG.db" = "url::https://bioconductor.org/packages/3.11/data/annotation/src/contrib/KEGG.db_3.2.4.tar.gz",
    "org.Pf.plasmo.db" = "url::https://bioconductor.org/packages/3.14/data/annotation/src/contrib/org.Pf.plasmo.db_3.14.0.tar.gz",
    "Azimuth" = "url::https://github.com/satijalab/azimuth/archive/HEAD.zip"
  )
  
  ## commented out entries are now in standard CRAN/cBio repo
  add_github("bigomics/PCSF")
  add_github("bigomics/playdata")
  add_github("bigomics/playbase")
  add_github("bigomics/bigdash")
  add_github("bigomics/bigLoaders")
  ##add_github("bigomics/fgsea")
  add_github("bigomics/wizardR")
  ##add_github("bigomics/biomaRt")
  add_github("GfellerLab/EPIC")
  add_github("broadinstitute/infercnv")
  add_github("GfellerLab/SuperCell")
  add_github("linxihui/NNLM")
  add_github("Coolgenome/iTALK")
  add_github("wt2015-github/FastGGM")
  add_github("satijalab/azimuth")
  #add_github("JohnCoene/waiter")
  add_github("JohnCoene/firebase@omics")
  add_github("JohnCoene/bsutils")
  #add_github("ropensci/iheatmapr")
  #add_github("rstudio/bslib@v0.6.1")
  #add_github("rstudio/htmltools")
  add_github("Bioconductor/BiocFileCache")
  #add_github("cysouw/qlcMatrix")
  #add_github("cole-trapnell-lab/leidenbase")
  add_github('cole-trapnell-lab/monocle3')
  add_github('bartongroup/Proteus')
  add_github('cran/riverplot')
  add_github('Ironholds/rgeolocate')

  pkg.remotes <- remotes.url[names(remotes.url) %in% pkg.used]
  pkg.imports <- setdiff(pkg.used, names(pkg.remotes))

  pkg.installed <- installed.packages()[,'Package']
  pkg.missing <- setdiff( c(pkg.imports,names(pkg.remotes)), pkg.installed)
  imports.missing <- setdiff(pkg.imports, pkg.installed)
  remotes.missing <- pkg.remotes[!(names(pkg.remotes) %in% pkg.installed)]
  
  list(
    pkg.used = pkg.used,
    pkg.installed = pkg.installed,
    pkg.missing = pkg.missing,
    imports = pkg.imports,
    remotes = pkg.remotes,    
    imports.missing = imports.missing,
    remotes.missing = remotes.missing
  )

}

install_dependencies <- function(use.remotes=FALSE) {
  
  require <- function(pkg) (pkg %in% installed.packages()[,'Package'])
  remove.pkg <- function(p) if(require(p)) try(remove.packages(p))

  if(!require("remotes")) install.packages('remotes')
  if(!require("BiocManager")) {
    remotes::install_version('BiocManager', version='1.30.23')
    if(!BiocManager::version()=="3.18") BiocManager::install(version='3.18')
  }

  if(use.remotes) {
    # install dependencies using remotes
    remotes::install_deps('.', dependencies = c("Imports","Remotes"))
  } else {
    ## Here we handle missing dependencies ourselves. Better control of
    ## skipping packages that are already installed.
    pkg <- scan_packages('R')
    if( length(pkg$imports.missing) || length(pkg$remotes.missing) ) {
      for(p in pkg$imports.missing) {
        if(!require(p)) BiocManager::install(p, ask=FALSE, dependencies=TRUE)
      }
      for(p in pkg$remotes.missing) {
        if(!require(p)) remotes::install_url(p, ask=FALSE, dependencies=TRUE)
      }
    } else {
      message("All dependencies installed. Nothing to install!")
    }
  }

}

